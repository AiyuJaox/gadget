#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*! \file density.c 
 *  \brief SPH density computation and smoothing length determination
 *
 *  This file contains the "first SPH loop", where the SPH densities and
 *  some auxiliary quantities are computed.  If the number of neighbours
 *  obtained falls outside the target range, the correct smoothing
 *  length is determined iteratively, if needed.
 */


#ifdef PERIODIC
static double boxSize, boxHalf;

#ifdef LONG_X
static double boxSize_X, boxHalf_X;
#else
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#endif
#ifdef LONG_Y
static double boxSize_Y, boxHalf_Y;
#else
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#endif
#ifdef LONG_Z
static double boxSize_Z, boxHalf_Z;
#else
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#endif
#endif


/*! This function computes the local density for each active SPH particle,
 *  the number of neighbours in the current smoothing radius, and the
 *  divergence and curl of the velocity field.  The pressure is updated as
 *  well.  If a particle with its smoothing region is fully inside the
 *  local domain, it is not exported to the other processors. The function
 *  also detects particles that have a number of neighbours outside the
 *  allowed tolerance range. For these particles, the smoothing length is
 *  adjusted accordingly, and the density computation is executed again.
 *  Note that the smoothing length is not allowed to fall below the lower
 *  bound set by MinGasHsml.
 */
void density(void)
{
  long long ntot, ntotleft;
  int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
  int i, j, n, ndone, npleft, maxfill, source, iter = 0;
  int level, ngrp, sendTask, recvTask, place, nexport;
  double dt_entr, tstart, tend, tstart_ngb = 0, tend_ngb = 0;
  double sumt, sumcomm, timengb, sumtimengb;
  double timecomp = 0, timeimbalance = 0, timecommsumm = 0, sumimbalance;
  MPI_Status status;

#ifdef NGB_MULTI_SEARCH
  int idx, idy, p, p1, p2;
  int len;
  //double ro_part;
#ifdef NGB_MULTI_SEARCH_GROUP  
  double d;
  int ii, main_flag, pt;
#endif
#endif

#ifdef PERIODIC
  boxSize = All.BoxSize;
  boxHalf = 0.5 * All.BoxSize;
#ifdef LONG_X
  boxHalf_X = boxHalf * LONG_X;
  boxSize_X = boxSize * LONG_X;
#endif
#ifdef LONG_Y
  boxHalf_Y = boxHalf * LONG_Y;
  boxSize_Y = boxSize * LONG_Y;
#endif
#ifdef LONG_Z
  boxHalf_Z = boxHalf * LONG_Z;
  boxSize_Z = boxSize * LONG_Z;
#endif
#endif


  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);

  for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {
      SphP[n].Left = SphP[n].Right = 0;

      if(P[n].Ti_endstep == All.Ti_Current)
	NumSphUpdate++;
    }

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);

  /* we will repeat the whole thing for those particles where we didn't
   * find enough neighbours
   */
  do
    {
      i = 0;			/* beginn with this index */
      ntotleft = ntot;		/* particles left for all tasks together */

      while(ntotleft > 0)
	{
	  for(j = 0; j < NTask; j++)
	    nsend_local[j] = 0;

	  /* do local particles and prepare export list */
	  tstart = second();
#ifdef NGB_MULTI_SEARCH
          memset(NgblistFlag, 0, All.MaxPart * (long) sizeof(int));
#endif
	  for(nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeDensity - NTask; i++)
	    if(P[i].Ti_endstep == All.Ti_Current)
	      {
#ifndef NGB_MULTI_SEARCH
		ndone++;
#endif

#ifdef NGB_MULTI_SEARCH
                if (NgblistFlag[i]) {
                  continue;
                } else {
                  ndone++;

                  NgblistCount = 0;
                  density_evaluate(i, 0);
                  NgblistFlag[i] = 1;
                  
                  len = NgblistCount;
                  
                  //printf("Ro: %d, %f, %f\n", len, SphP[i].Hsml, ((double)len / SphP[i].Hsml));
                  // Density of neighbors
                  /*if ((double)len / SphP[i].Hsml < 0.015) {
                    continue;
                  }*/
       
                  idy = 0;
                  for(idx = 0; idx < len; idx++) {
                    p = Ngblist[idx];
                    if (NgblistFlag[p] || P[p].Ti_endstep != All.Ti_Current) {
	              continue;
                    }
                    NgblistMulti[idy] = p;
                    idy++;  
                  }
                  len = idy;
                  if (len < 2) continue;
                  if (len & 1) len--;

#ifdef NGB_MULTI_SEARCH_GROUP   
                  // Fandorin Get two particle
                  main_flag = 0; 
                  for(idx = 0; idx < len; idx+=2) {
                    p1 = NgblistMulti[idx];
                    p2 = NgblistMulti[idx + 1];

                    d = 0;
                    for(ii = 0; ii < 3; ii++) {
                      d += (P[p1].Pos[ii] - P[p2].Pos[ii]) * (P[p1].Pos[ii] - P[p2].Pos[ii]);
                    }
                    d = sqrt(d);
                    // Check distance between of the centers
                    if (d > 0.5 * (SphP[p1].Hsml + SphP[p2].Hsml)) {
                     // printf("Continue!\n");
                      continue;
                    } 
                    main_flag = 1;
                    break;
                  }

                  if (main_flag == 0 && len != 0) {
                    printf("OPA!!!\n");
                  }
                  if (main_flag) {
                    density_evaluate_group(p1, p2);
                    /*idy = 0; 
                    for(idx = 0; idx < len; idx++) {
                      pt = NgblistMulti[idx];
                      if (pt == p1 || pt == p2) continue;
                      for(ii = 0; ii < 3; ii++) {
                        d += (P[pt].Pos[ii] - P[p_m].Pos[ii]) * (P[pt].Pos[ii] - P[p_m].Pos[ii]);
                      } 
                      if () 
                      NgblistMultiGroup[idy] = pt;
                      idy++;
                    }*/
                  }
#endif
                  for(idx = 0; idx < len; idx+=2) {
                    
                    p1 = NgblistMulti[idx];           
                    p2 = NgblistMulti[idx + 1];
                    
                    /*if (Father[p1] != Father[p2]) {
                      continue;
                    }*/
                    if (density_evaluate_multi(p1, p2, 0) == 0) {
                      density_evaluate(p1, 0);
                      density_evaluate(p2, 0);
                    }
                    NgblistFlag[p1] = 1;
                    NgblistFlag[p2] = 1;
                    ndone += 2;
                    
                  }
	          
                  continue;
                }
#else
                density_evaluate(i, 0);
#endif
                for(j = 0; j < NTask; j++)
                  Exportflag[j] = 0;

		for(j = 0; j < NTask; j++)
		  {
		    if(Exportflag[j])
		      {
                  	DensDataIn[nexport].Pos[0] = P[i].Pos[0];
			DensDataIn[nexport].Pos[1] = P[i].Pos[1];
			DensDataIn[nexport].Pos[2] = P[i].Pos[2];
			DensDataIn[nexport].Vel[0] = SphP[i].VelPred[0];
			DensDataIn[nexport].Vel[1] = SphP[i].VelPred[1];
			DensDataIn[nexport].Vel[2] = SphP[i].VelPred[2];
			DensDataIn[nexport].Hsml = SphP[i].Hsml;
			DensDataIn[nexport].Index = i;
			DensDataIn[nexport].Task = j;
			nexport++;
			nsend_local[j]++;
		      }
		  }
	      }
	  tend = second();
	  timecomp += timediff(tstart, tend);
          TimeMultiSearch += timediff(tstart, tend);

          printf("Multi search time: %.4f, %.4f\n", timediff(tstart, tend), TimeMultiSearch);

	  qsort(DensDataIn, nexport, sizeof(struct densdata_in), dens_compare_key);

	  for(j = 1, noffset[0] = 0; j < NTask; j++)
	    noffset[j] = noffset[j - 1] + nsend_local[j - 1];

	  tstart = second();

	  MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

	  tend = second();
	  timeimbalance += timediff(tstart, tend);


	  /* now do the particles that need to be exported */

	  for(level = 1; level < (1 << PTask); level++)
	    {
	      tstart = second();
	      for(j = 0; j < NTask; j++)
		nbuffer[j] = 0;
	      for(ngrp = level; ngrp < (1 << PTask); ngrp++)
		{
		  maxfill = 0;
		  for(j = 0; j < NTask; j++)
		    {
		      if((j ^ ngrp) < NTask)
			if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
			  maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		    }
		  if(maxfill >= All.BunchSizeDensity)
		    break;

		  sendTask = ThisTask;
		  recvTask = ThisTask ^ ngrp;

		  if(recvTask < NTask)
		    {
		      if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
			{
			  /* get the particles */
			  MPI_Sendrecv(&DensDataIn[noffset[recvTask]],
				       nsend_local[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				       recvTask, TAG_DENS_A,
				       &DensDataGet[nbuffer[ThisTask]],
				       nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_in),
				       MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
			}
		    }

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
		}
	      tend = second();
	      timecommsumm += timediff(tstart, tend);


	      tstart = second();
	      for(j = 0; j < nbuffer[ThisTask]; j++) {
		density_evaluate(j, 1);
#ifdef NGB_MULTI_SEARCH
                // TODO
#endif
              }
	      tend = second();
	      timecomp += timediff(tstart, tend);

	      /* do a block to explicitly measure imbalance */
	      tstart = second();
	      MPI_Barrier(MPI_COMM_WORLD);
	      tend = second();
	      timeimbalance += timediff(tstart, tend);

	      /* get the result */
	      tstart = second();
	      for(j = 0; j < NTask; j++)
		nbuffer[j] = 0;
	      for(ngrp = level; ngrp < (1 << PTask); ngrp++)
		{
		  maxfill = 0;
		  for(j = 0; j < NTask; j++)
		    {
		      if((j ^ ngrp) < NTask)
			if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
			  maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		    }
		  if(maxfill >= All.BunchSizeDensity)
		    break;

		  sendTask = ThisTask;
		  recvTask = ThisTask ^ ngrp;

		  if(recvTask < NTask)
		    {
		      if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
			{
			  /* send the results */
			  MPI_Sendrecv(&DensDataResult[nbuffer[ThisTask]],
				       nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_out),
				       MPI_BYTE, recvTask, TAG_DENS_B,
				       &DensDataPartialResult[noffset[recvTask]],
				       nsend_local[recvTask] * sizeof(struct densdata_out),
				       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);

			  /* add the result to the particles */
			  for(j = 0; j < nsend_local[recvTask]; j++)
			    {
			      source = j + noffset[recvTask];
			      place = DensDataIn[source].Index;

			      SphP[place].NumNgb += DensDataPartialResult[source].Ngb;
			      SphP[place].Density += DensDataPartialResult[source].Rho;
			      SphP[place].DivVel += DensDataPartialResult[source].Div;

			      SphP[place].DhsmlDensityFactor += DensDataPartialResult[source].DhsmlDensity;

			      SphP[place].Rot[0] += DensDataPartialResult[source].Rot[0];
			      SphP[place].Rot[1] += DensDataPartialResult[source].Rot[1];
			      SphP[place].Rot[2] += DensDataPartialResult[source].Rot[2];
			    }
			}
		    }

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
		}
	      tend = second();
	      timecommsumm += timediff(tstart, tend);

	      level = ngrp - 1;
	    }

	  MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
	  for(j = 0; j < NTask; j++)
	    ntotleft -= ndonelist[j];
	}



      /* do final operations on results */
      tstart = second();
      for(i = 0, npleft = 0; i < N_gas; i++)
	{
	  if(P[i].Ti_endstep == All.Ti_Current)
	    {
	      {
		SphP[i].DhsmlDensityFactor =
		  1 / (1 + SphP[i].Hsml * SphP[i].DhsmlDensityFactor / (NUMDIMS * SphP[i].Density));

		SphP[i].CurlVel = sqrt(SphP[i].Rot[0] * SphP[i].Rot[0] +
				       SphP[i].Rot[1] * SphP[i].Rot[1] +
				       SphP[i].Rot[2] * SphP[i].Rot[2]) / SphP[i].Density;

		SphP[i].DivVel /= SphP[i].Density;

		dt_entr = (All.Ti_Current - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;

		SphP[i].Pressure =
		  (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, GAMMA);
	      }


	      /* now check whether we had enough neighbours */

	      if(SphP[i].NumNgb < (All.DesNumNgb - All.MaxNumNgbDeviation) ||
		 (SphP[i].NumNgb > (All.DesNumNgb + All.MaxNumNgbDeviation)
		  && SphP[i].Hsml > (1.01 * All.MinGasHsml)))
		{
		  /* need to redo this particle */
		  npleft++;

		  if(SphP[i].Left > 0 && SphP[i].Right > 0)
		    if((SphP[i].Right - SphP[i].Left) < 1.0e-3 * SphP[i].Left)
		      {
			/* this one should be ok */
			npleft--;
			P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive */
			continue;
		      }

		  if(SphP[i].NumNgb < (All.DesNumNgb - All.MaxNumNgbDeviation))
		    SphP[i].Left = dmax(SphP[i].Hsml, SphP[i].Left);
		  else
		    {
		      if(SphP[i].Right != 0)
			{
			  if(SphP[i].Hsml < SphP[i].Right)
			    SphP[i].Right = SphP[i].Hsml;
			}
		      else
			SphP[i].Right = SphP[i].Hsml;
		    }

		  if(iter >= MAXITER - 10)
		    {
		      printf
			("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
			 i, ThisTask, (int) P[i].ID, SphP[i].Hsml, SphP[i].Left, SphP[i].Right,
			 (float) SphP[i].NumNgb, SphP[i].Right - SphP[i].Left, P[i].Pos[0], P[i].Pos[1],
			 P[i].Pos[2]);
		      fflush(stdout);
		    }

		  if(SphP[i].Right > 0 && SphP[i].Left > 0)
		    SphP[i].Hsml = pow(0.5 * (pow(SphP[i].Left, 3) + pow(SphP[i].Right, 3)), 1.0 / 3);
		  else
		    {
		      if(SphP[i].Right == 0 && SphP[i].Left == 0)
			endrun(8188);	/* can't occur */

		      if(SphP[i].Right == 0 && SphP[i].Left > 0)
			{
			  if(P[i].Type == 0 && fabs(SphP[i].NumNgb - All.DesNumNgb) < 0.5 * All.DesNumNgb)
			    {
			      SphP[i].Hsml *=
				1 - (SphP[i].NumNgb -
				     All.DesNumNgb) / (NUMDIMS * SphP[i].NumNgb) * SphP[i].DhsmlDensityFactor;
			    }
			  else
			    SphP[i].Hsml *= 1.26;
			}

		      if(SphP[i].Right > 0 && SphP[i].Left == 0)
			{
			  if(P[i].Type == 0 && fabs(SphP[i].NumNgb - All.DesNumNgb) < 0.5 * All.DesNumNgb)
			    {
			      SphP[i].Hsml *=
				1 - (SphP[i].NumNgb -
				     All.DesNumNgb) / (NUMDIMS * SphP[i].NumNgb) * SphP[i].DhsmlDensityFactor;
			    }
			  else
			    SphP[i].Hsml /= 1.26;
			}
		    }

		  if(SphP[i].Hsml < All.MinGasHsml)
		    SphP[i].Hsml = All.MinGasHsml;
		}
	      else
		P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive */
	    }
	}
      tend = second();
      timecomp += timediff(tstart, tend);


      numlist = malloc(NTask * sizeof(int) * NTask);
      MPI_Allgather(&npleft, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
      for(i = 0, ntot = 0; i < NTask; i++)
	ntot += numlist[i];
      free(numlist);

      if(ntot > 0)
	{
	  if(iter == 0)
	    tstart_ngb = second();

	  iter++;

	  if(iter > 0 && ThisTask == 0)
	    {
	      printf("ngb iteration %d: need to repeat for %d%09d particles.\n", iter,
		     (int) (ntot / 1000000000), (int) (ntot % 1000000000));
	      fflush(stdout);
	    }

	  if(iter > MAXITER)
	    {
	      printf("failed to converge in neighbour iteration in density()\n");
	      fflush(stdout);
	      endrun(1155);
	    }
	}
      else
	tend_ngb = second();
    }
  while(ntot > 0);


  /* mark as active again */
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep < 0)
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;

  free(ndonelist);
  free(nsend);
  free(nsend_local);
  free(nbuffer);
  free(noffset);


  /* collect some timing information */
  if(iter > 0)
    timengb = timediff(tstart_ngb, tend_ngb);
  else
    timengb = 0;

  MPI_Reduce(&timengb, &sumtimengb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecomp, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecommsumm, &sumcomm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      All.CPU_HydCompWalk += sumt / NTask;
      All.CPU_HydCommSumm += sumcomm / NTask;
      All.CPU_HydImbalance += sumimbalance / NTask;
      All.CPU_EnsureNgb += sumtimengb / NTask;
    }
}

#ifdef NGB_MULTI_SEARCH
FLOAT ngb_center(FLOAT pos_m[3], FLOAT pos_m1[3], FLOAT pos[3], FLOAT pos1[3], FLOAT h, FLOAT h1, FLOAT *h2) {
  FLOAT d, d1, d2, d3;
  int i;
  FLOAT p3[3], p4[3], p5[3], p6[3];

  d = 0;
  for(i = 0; i < 3; i++) {
    d += (pos[i] - pos1[i]) * (pos[i] - pos1[i]);
  }
  d = sqrt(d);

/*
  if (d == 0) {
    printf("NGB center warning: distance is zero, maybe stop?\n");
    return 0;
  }
*/

  d1 = d2 = d3 = 0;
  for(i = 0; i < 3; i++) {
    p3[i] = ((h + d) / d) * (pos[i] - pos1[i]) + pos1[i];
    p4[i] = ((h1 + d) / d) * (pos1[i] - pos[i]) + pos[i];
  
    p5[i] = (h / d) * (pos1[i] - pos[i]) + pos[i];
    p6[i] = (h1 / d) * (pos[i] - pos1[i]) + pos1[i];

    pos_m[i] = (p3[i] + p4[i]) / 2;
    pos_m1[i] = (p5[i] + p6[i]) / 2;

    d1 += (pos[i] - p4[i]) * (pos[i] - p4[i]);
    d2 += (pos1[i] - p3[i]) * (pos1[i] - p3[i]);
    d3 += (pos_m1[i] - p5[i]) * (pos_m1[i] - p5[i]);
  }
  *h2 = d3;//sqrt(d3);
 
  /*for(i = 0; i < 3; i++) {
    pos_m[i] = (p3[i] + p4[i]) / 2;
  }*/
 
  /*d1 = 0;
  for(i = 0; i < 3; i++) {
    d1 += (pos[i] - p4[i]) * (pos[i] - p4[i]); 
  }

  d2 = 0;
  for(i = 0; i < 3; i++) {
    d2 += (pos1[i] - p3[i]) * (pos1[i] - p3[i]);
  }
*/
  if (d1 < h * h) {
    for(i = 0; i < 3; i++) {
      pos_m[i] = pos[i];
      pos_m1[i] = pos1[i];
    }
    d = h;
    *h2 = h1;
  } else if (d2 < h1 * h1) {
    for(i = 0; i < 3; i++) {
      pos_m[i] = pos1[i];
      pos_m1[i] = pos[i];
    }
    d = h1;
    *h2 = h;
  } else {
    d = 0;
    for(i = 0; i < 3; i++) {
      d += (p4[i] - pos_m[i]) * (p4[i] - pos_m[i]);
    }
    d = sqrt(d);
  }

  return d;
}

#ifdef NGB_MULTI_SEARCH_GROUP
int density_evaluate_multi_group(int p1, int p2, int mode) {
  int i, j, n, startnode, startnode1, numngb, numngb_inbox;
  double d, h, h1, h2, h12, h_m, fac, hinv, hinv1, hinv3, hinv4, hinv13, hinv14;
  double rho, divv, wk, dwk, rho1, divv1;
  double dx, dy, dz, dx1, dy1, dz1, r, r2, r12, r22, u, mass_j;
  double dvx, dvy, dvz, rotv[3], rotv1[3];
  double weighted_numngb, dhsmlrho, weighted_numngb1, dhsmlrho1;
  FLOAT *pos, *pos1, *vel, *vel1;
  FLOAT pos_m[3], pos_m1[3];
  FLOAT h3;

  int target = p1;
  int target1 = p2;

  if(mode == 0) {
    pos = P[target].Pos;
    vel = SphP[target].VelPred;
    h = SphP[target].Hsml;

    pos1 = P[target1].Pos;
    vel1 = SphP[target1].VelPred;
    h1 = SphP[target1].Hsml;
  } else {
    pos = DensDataGet[target].Pos;
    vel = DensDataGet[target].Vel;
    h = DensDataGet[target].Hsml;

    pos1 = DensDataGet[target1].Pos;
    vel1 = DensDataGet[target1].Vel;
    h1 = DensDataGet[target1].Hsml;
  }



}
#endif

int density_evaluate_multi(int p1, int p2, int mode) {
  int i, j, n, startnode, startnode1, numngb, numngb_inbox;
  //int m, flag, diff1, diff2, numngb_inbox1, numngb_inbox2;
  double d, h, h1, h2, h12, h_m, fac, hinv, hinv1, hinv3, hinv4, hinv13, hinv14;
  double rho, divv, wk, dwk, rho1, divv1;
  double dx, dy, dz, dx1, dy1, dz1, r, r2, r12, r22, u, mass_j;
  double dvx, dvy, dvz, rotv[3], rotv1[3];
  double weighted_numngb, dhsmlrho, weighted_numngb1, dhsmlrho1;
  FLOAT *pos, *pos1, *vel, *vel1;
  FLOAT pos_m[3], pos_m1[3];
  FLOAT h3;

  int *multi_res, *p1_res, *p2_res;
/*
  multi_res = (int*)malloc(MAX_NGB * sizeof(int));
  p1_res = (int*)malloc(MAX_NGB * sizeof(int));
  p2_res = (int*)malloc(MAX_NGB * sizeof(int));
*/
  int target = p1;
  int target1 = p2;

  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = SphP[target].VelPred;
      h = SphP[target].Hsml;

      pos1 = P[target1].Pos;
      vel1 = SphP[target1].VelPred;
      h1 = SphP[target1].Hsml;
    }
  else
    {
      pos = DensDataGet[target].Pos;
      vel = DensDataGet[target].Vel;
      h = DensDataGet[target].Hsml;

      pos1 = DensDataGet[target1].Pos;
      vel1 = DensDataGet[target1].Vel;
      h1 = DensDataGet[target1].Hsml;
    }

  d = 0;
  for(i = 0; i < 3; i++) {
    d += (pos[i] - pos1[i]) * (pos[i] - pos1[i]);
  }
  d = sqrt(d); 
  
  if (d > 0.25 * (h + h1)) {
    return 0;
  }

  //MultiAlgCount++; 
  //printf("MultiAlgCount: %d\n", MultiAlgCount);

  h3 = 0;

  h_m = ngb_center(&pos_m[0], &pos_m1[0], &pos[0], &pos1[0], h, h1, &h3);

// target
  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
  hinv4 = hinv3 * hinv;

  rho = divv = rotv[0] = rotv[1] = rotv[2] = 0;
  weighted_numngb = 0;
  dhsmlrho = 0;

  // target1
  h12 = h1 * h1;
  hinv1 = 1.0 / h1;
#ifndef  TWODIMS
  hinv13 = hinv1 * hinv1 * hinv1;
#else
  hinv13 = hinv1 * hinv1 / boxSize_Z;
#endif
  hinv14 = hinv13 * hinv1;

  rho1 = divv1 = rotv1[0] = rotv1[1] = rotv1[2] = 0;
  weighted_numngb1 = 0;
  dhsmlrho1 = 0;

  startnode = All.MaxPart;
  numngb = 0;

// target
  do
    {

/*
      numngb_inbox1 = ngb_treefind_variable(&pos[0], h, &startnode, target);

      for(n = 0; n < numngb_inbox1; n++) {
        p1_res[n] = Ngblist[n];
      }

      startnode = All.MaxPart;
      numngb_inbox2 = ngb_treefind_variable(&pos1[0], h1, &startnode, target);
      for(n = 0; n < numngb_inbox2; n++) {
        p2_res[n] = Ngblist[n];
      }
*/
      numngb_inbox = ngb_treefind_variable(&pos_m[0], h_m, &startnode, 0);
 /*
      for(n = 0; n < numngb_inbox; n++) {
        multi_res[n] = Ngblist[n];
      }

      diff1 = 0;
      for(n = 0; n < numngb_inbox1; n++) {
        flag = 0;
        for(m = 0; m < numngb_inbox; m++) {
          if (p1_res[n] == multi_res[m]) {
            flag = 1;
            break;
          }
        }
        if (flag == 0) {
          diff1++;
        }
      }

      diff2 = 0;
      for(n = 0; n < numngb_inbox2; n++) {
        flag = 0;
        for(m = 0; m < numngb_inbox; m++) {
          if (p2_res[n] == multi_res[m]) {
            flag = 1;
            break;
          }
        }
        if (flag == 0) {
          diff2++;
        }
      }
      
      printf("Ngb (%d, %d, %d), ", numngb_inbox1, numngb_inbox2, numngb_inbox);
      printf("diff1 = %d, diff2 = %d\n", diff1, diff2);
      if (diff1 != 0 || diff2 != 0) {
        printf("NGB center warning: different of neighbors, maybe stop?\n");
      }
*/
      int flg = 0;

      for(n = 0; n < numngb_inbox; n++)
        {
          j = Ngblist[n];
          //flg = 0;

          // Intersection
          dx = pos_m1[0] - P[j].Pos[0];
          dy = pos_m1[1] - P[j].Pos[1];
          dz = pos_m1[2] - P[j].Pos[2];
          
          r22 = dx * dx + dy * dy + dz * dz;

          dx = pos[0] - P[j].Pos[0];
          dy = pos[1] - P[j].Pos[1];
          dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC
          if(dx > boxHalf_X)
            dx -= boxSize_X;
          if(dx < -boxHalf_X)
            dx += boxSize_X;
          if(dy > boxHalf_Y)
            dy -= boxSize_Y;
          if(dy < -boxHalf_Y)
            dy += boxSize_Y;
          if(dz > boxHalf_Z)
            dz -= boxSize_Z;
          if(dz < -boxHalf_Z)
            dz += boxSize_Z;
#endif

          r2 = dx * dx + dy * dy + dz * dz;

          dx1 = pos1[0] - P[j].Pos[0];
          dy1 = pos1[1] - P[j].Pos[1];
          dz1 = pos1[2] - P[j].Pos[2];

#ifdef PERIODIC
          if(dx1 > boxHalf_X)
            dx1 -= boxSize_X;
          if(dx1 < -boxHalf_X)
            dx1 += boxSize_X;
          if(dy1 > boxHalf_Y)
            dy1 -= boxSize_Y;
          if(dy1 < -boxHalf_Y)
            dy1 += boxSize_Y;
          if(dz1 > boxHalf_Z)
            dz1 -= boxSize_Z;
          if(dz1 < -boxHalf_Z)
            dz1 += boxSize_Z;
#endif
          r12 = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;

          if (r22 <= h3 || (r2 < h2 && r12 < h12)) {
            flg = 1;
          } else {
            if (r2 < h2) {
              flg = 2;
            } else if (r12 < h12) {
              flg = 3;
            }  
          }
         
          if (flg == 0) continue;
                    
              if (flg != 3) {
                r = sqrt(r2);
                u = r * hinv;

                if(u < 0.5)
                  {
                    wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                    dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
                  }
                else
                  {
                    wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
                    dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
                  }

                mass_j = P[j].Mass;
                rho += mass_j * wk;
                weighted_numngb += NORM_COEFF * wk / hinv3;
                dhsmlrho += -mass_j * (NUMDIMS * hinv * wk + u * dwk);

                if(r > 0)
                  {
                    fac = mass_j * dwk / r;
                    dvx = vel[0] - SphP[j].VelPred[0];
                    dvy = vel[1] - SphP[j].VelPred[1];
                    dvz = vel[2] - SphP[j].VelPred[2];
                    divv -= fac * (dx * dvx + dy * dvy + dz * dvz);
                    rotv[0] += fac * (dz * dvy - dy * dvz);
                    rotv[1] += fac * (dx * dvz - dz * dvx);
                    rotv[2] += fac * (dy * dvx - dx * dvy);
                  }
              }
 
              if (flg != 2) {
                r = sqrt(r12);
                u = r * hinv1;

                if(u < 0.5)
                  {
                    wk = hinv13 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                    dwk = hinv14 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
                  }
                else
                  {
                    wk = hinv13 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
                    dwk = hinv14 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
                  }

                mass_j = P[j].Mass;
                rho1 += mass_j * wk;
                weighted_numngb1 += NORM_COEFF * wk / hinv13;
                dhsmlrho1 += -mass_j * (NUMDIMS * hinv1 * wk + u * dwk);
                if(r > 0)
                  {
                    fac = mass_j * dwk / r;
                    dvx = vel1[0] - SphP[j].VelPred[0];
                    dvy = vel1[1] - SphP[j].VelPred[1];
                    dvz = vel1[2] - SphP[j].VelPred[2];
                    divv1 -= fac * (dx1 * dvx + dy1 * dvy + dz1 * dvz);
                    rotv1[0] += fac * (dz1 * dvy - dy1 * dvz);
                    rotv1[1] += fac * (dx1 * dvz - dz1 * dvx);
                    rotv1[2] += fac * (dy1 * dvx - dx1 * dvy);
                  }
              }
           flg = 0;
        }


    }
  while(startnode >= 0);

  if(mode == 0)
    {
      SphP[target].NumNgb = weighted_numngb;
      SphP[target].Density = rho;
      SphP[target].DivVel = divv;
      SphP[target].DhsmlDensityFactor = dhsmlrho;
      SphP[target].Rot[0] = rotv[0];
      SphP[target].Rot[1] = rotv[1];
      SphP[target].Rot[2] = rotv[2];
    }
  else
    {
      DensDataResult[target].Rho = rho;
      DensDataResult[target].Div = divv;
      DensDataResult[target].Ngb = weighted_numngb;
      DensDataResult[target].DhsmlDensity = dhsmlrho;
      DensDataResult[target].Rot[0] = rotv[0];
      DensDataResult[target].Rot[1] = rotv[1];
      DensDataResult[target].Rot[2] = rotv[2];
    }


// target1
/*  numngb = 0;
  
      startnode1 = -1; 
      for(n = 0; n < numngb_inbox; n++)
        {
          j = Ngblist[n];
 
          dx = pos1[0] - P[j].Pos[0];
          dy = pos1[1] - P[j].Pos[1];
          dz = pos1[2] - P[j].Pos[2];

#ifdef PERIODIC 
          if(dx > boxHalf_X)
            dx -= boxSize_X;
          if(dx < -boxHalf_X)
            dx += boxSize_X;
          if(dy > boxHalf_Y)
            dy -= boxSize_Y;
          if(dy < -boxHalf_Y)
            dy += boxSize_Y;
          if(dz > boxHalf_Z)
            dz -= boxSize_Z;
          if(dz < -boxHalf_Z)
            dz += boxSize_Z;
#endif


          r2 = dx * dx + dy * dy + dz * dz;

          if(r2 < h12)
            {
              numngb++;

              r = sqrt(r2);
              u = r * hinv1;

              if(u < 0.5)
                {
                  wk = hinv13 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                  dwk = hinv14 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
                }
              else
                {
                  wk = hinv13 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
                  dwk = hinv14 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
                }

              mass_j = P[j].Mass;
              rho1 += mass_j * wk;
              weighted_numngb1 += NORM_COEFF * wk / hinv13;
              dhsmlrho1 += -mass_j * (NUMDIMS * hinv1 * wk + u * dwk);
              if(r > 0)
                {
                  fac = mass_j * dwk / r;
                  dvx = vel1[0] - SphP[j].VelPred[0];
                  dvy = vel1[1] - SphP[j].VelPred[1];
                  dvz = vel1[2] - SphP[j].VelPred[2];
                  divv1 -= fac * (dx * dvx + dy * dvy + dz * dvz);
                  rotv1[0] += fac * (dz * dvy - dy * dvz);
                  rotv1[1] += fac * (dx * dvz - dz * dvx);
                  rotv1[2] += fac * (dy * dvx - dx * dvy);
                }
            }
        }
*/  
  if(mode == 0)
    {
      SphP[target1].NumNgb = weighted_numngb1;
      SphP[target1].Density = rho1;
      SphP[target1].DivVel = divv1;
      SphP[target1].DhsmlDensityFactor = dhsmlrho1;
      SphP[target1].Rot[0] = rotv1[0];
      SphP[target1].Rot[1] = rotv1[1];
      SphP[target1].Rot[2] = rotv1[2];
    }
  else
    {
      DensDataResult[target1].Rho = rho1;
      DensDataResult[target1].Div = divv1;
      DensDataResult[target1].Ngb = weighted_numngb1;
      DensDataResult[target1].DhsmlDensity = dhsmlrho1;
      DensDataResult[target1].Rot[0] = rotv1[0];
      DensDataResult[target1].Rot[1] = rotv1[1];
      DensDataResult[target1].Rot[2] = rotv1[2];
    }
/*
  free(multi_res);
  free(p1_res); 
  free(p2_res);
*/
  return 1;
}
#endif

/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
void density_evaluate(int target, int mode)
{
  int j, n, startnode, numngb, numngb_inbox;
  double h, h2, fac, hinv, hinv3, hinv4;
  double rho, divv, wk, dwk;
  double dx, dy, dz, r, r2, u, mass_j;
  double dvx, dvy, dvz, rotv[3];
  double weighted_numngb, dhsmlrho;
  FLOAT *pos, *vel;

  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = SphP[target].VelPred;
      h = SphP[target].Hsml;
    }
  else
    {
      pos = DensDataGet[target].Pos;
      vel = DensDataGet[target].Vel;
      h = DensDataGet[target].Hsml;
    }

  //BaseAlgCount++;
  //printf("BaseAlgCount: %d\n", BaseAlgCount);

  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
  hinv4 = hinv3 * hinv;

  rho = divv = rotv[0] = rotv[1] = rotv[2] = 0;
  weighted_numngb = 0;
  dhsmlrho = 0;

  startnode = All.MaxPart;
  numngb = 0;
  
  do
    {
     
#ifdef BOTTOM_UP_WALK
      ngb_search_startnode_variable(&pos[0], h, &startnode, target);
#endif      
      numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode, target);
#ifdef NGB_MULTI_SEARCH
      NgblistCount = numngb_inbox;
#endif
      for(n = 0; n < numngb_inbox; n++)
	{
	  j = Ngblist[n];

	  dx = pos[0] - P[j].Pos[0];
	  dy = pos[1] - P[j].Pos[1];
	  dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	  if(dx > boxHalf_X)
	    dx -= boxSize_X;
	  if(dx < -boxHalf_X)
	    dx += boxSize_X;
	  if(dy > boxHalf_Y)
	    dy -= boxSize_Y;
	  if(dy < -boxHalf_Y)
	    dy += boxSize_Y;
	  if(dz > boxHalf_Z)
	    dz -= boxSize_Z;
	  if(dz < -boxHalf_Z)
	    dz += boxSize_Z;
#endif
	  r2 = dx * dx + dy * dy + dz * dz;

	  if(r2 < h2)
	    {
	      numngb++;

	      r = sqrt(r2);

	      u = r * hinv;

	      if(u < 0.5)
		{
		  wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
		  dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
		}
	      else
		{
		  wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
		  dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		}

	      mass_j = P[j].Mass;

	      rho += mass_j * wk;

	      weighted_numngb += NORM_COEFF * wk / hinv3;

	      dhsmlrho += -mass_j * (NUMDIMS * hinv * wk + u * dwk);

	      if(r > 0)
		{
		  fac = mass_j * dwk / r;

		  dvx = vel[0] - SphP[j].VelPred[0];
		  dvy = vel[1] - SphP[j].VelPred[1];
		  dvz = vel[2] - SphP[j].VelPred[2];

		  divv -= fac * (dx * dvx + dy * dvy + dz * dvz);

		  rotv[0] += fac * (dz * dvy - dy * dvz);
		  rotv[1] += fac * (dx * dvz - dz * dvx);
		  rotv[2] += fac * (dy * dvx - dx * dvy);
		}
	    }
	}
    }
  while(startnode >= 0);

  if(mode == 0)
    {
      SphP[target].NumNgb = weighted_numngb;
      SphP[target].Density = rho;
      SphP[target].DivVel = divv;
      SphP[target].DhsmlDensityFactor = dhsmlrho;
      SphP[target].Rot[0] = rotv[0];
      SphP[target].Rot[1] = rotv[1];
      SphP[target].Rot[2] = rotv[2];
    }
  else
    {
      DensDataResult[target].Rho = rho;
      DensDataResult[target].Div = divv;
      DensDataResult[target].Ngb = weighted_numngb;
      DensDataResult[target].DhsmlDensity = dhsmlrho;
      DensDataResult[target].Rot[0] = rotv[0];
      DensDataResult[target].Rot[1] = rotv[1];
      DensDataResult[target].Rot[2] = rotv[2];
    }
}




/*! This routine is a comparison kernel used in a sort routine to group
 *  particles that are exported to the same processor.
 */
int dens_compare_key(const void *a, const void *b)
{
  if(((struct densdata_in *) a)->Task < (((struct densdata_in *) b)->Task))
    return -1;

  if(((struct densdata_in *) a)->Task > (((struct densdata_in *) b)->Task))
    return +1;

  return 0;
}

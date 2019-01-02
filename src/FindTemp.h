/* ===============================================================================================
   FindTemp.h  (source code generated 2017-05-23)
   Copyright (c) Jiahui Guan.
   Copyright (c) Patrice Koehl.
   =============================================================================================== */

#ifndef _FINDTEMP_
#define _FINDTEMP_

#include <math.h>
#include <cstdlib>
#include "mt.h"

using namespace std;

/* ===============================================================================================
   Prototypes
   =============================================================================================== */

void FindTemp(double *DistMat, int N, int M, int P, int Neigen, double *Ensemble, double *Q,
                double *Degree, double *WorkD, int *WorkI, double *Centers, double Cutoff,
                int iter, int Nrepeat, bool *Membership, int *ClusterID, int Nmin, int *Nc, double *Temp,
                int *Ntemp, int nthreads, MersenneTwister **mt);

double DefineScale(int N, double *DistMat);

double FindTmax(double *DistMat, int N, int M, int P, int Neigen, double *Ensemble, double *Q,
	double *Degree, double *WorkD, int *WorkI, double *Centers, double Cutoff,
	int iter, int Nrepeat, bool *Membership, int *ClusterID, double Tinit,
	int nthreads, MersenneTwister **mt);

double FindT0(int *Nmin, double *DistMat, int N, int M, int P, int Neigen, double *Ensemble, double *Q,
	double *Degree, double *WorkD, int *WorkI, double *Centers, double Cutoff,
	int iter, int Nrepeat, bool *Membership, int *ClusterID, double Tinit,
	int nthreads, MersenneTwister **mt);

double BracketTemp(double T1, double T2, int Ntarget, double *DistMat, int N, int M, int P, int Neigen,
        double *Ensemble, double *Q,double *Degree, double *WorkD, int *WorkI, double *Centers, double Cutoff,
        int iter, int Nrepeat, bool *Membership, int *nc, int *ClusterID, int nthreads, MersenneTwister **mt);


#endif

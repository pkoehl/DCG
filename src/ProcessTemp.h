/* ===============================================================================================

ProcessTemp.h -part of DCG 
Copyright (c) Jiahui Guan.
Copyright (c) Patrice Koehl.

================================================================================================== */

#ifndef _PROCESSTEMP_
#define _PROCESSTEMP_

#include <iostream>
#include <cmath>
#include "mt.h"

using namespace std;

/* ===============================================================================================
   Prototypes for all procedures related to ProcessTemp
================================================================================================== */

void ProcessTemp(double *DistMat, int N, double Temp, int M, int P, int Neigen, double *Ensemble, 
	double *Q, double *Degree, double *WorkD, int *WorkI, double *Centers, double Cutoff, 
	int iter, int Nrepeats, bool *Membership , int *nc, int *ClusterID, 
	int nthread, MersenneTwister **mt, bool iprint); 

void NormalizeEnsemble(double *Ensemble, double *B, int N);
void NormalizeEnsemble2(double *Ensemble, double *B, int N);

#endif

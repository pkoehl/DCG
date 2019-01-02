/* ===============================================================================================
 RandomWalk.h part of DCG  (source code generated 2017-05-21)
 Copyright (c) Jiahui Guan.
 Copyright (c) Patrice Koehl.

================================================================================================== */

#ifndef _RandomWalk_
#define _RandomWalk_

#include <iostream>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include "mt.h"

/* ===============================================================================================
   Prototypes for the different functions needed to generate a RandomWalk
================================================================================================== */

void RandomWalk(double *DistMat, double *Q, double *Degree, int N, int M, int P, double T, double *Ensemble, int nthread, MersenneTwister **mt);

void TransMatrix(int N, double Temp, double *DistMat, double *Q, double *Degree);

int sampleGenerator(int N, double *P, MersenneTwister *mt);

void SymmetricEnsemble (double *Ensemble, int N);

void* worker_thread(void* data);

#endif

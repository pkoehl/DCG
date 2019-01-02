/* ===============================================================================================
 SpectralClustering.h part of DCG  (source code generated 2017-05-21)
 Copyright (c) Jiahui Guan.
 Copyright (c) Patrice Koehl.

================================================================================================== */

#ifndef _SPECTRALCLUSTERING_
#define _SPECTRALCLUSTERING_

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include "mt.h"

/* ===============================================================================================
 We will use some Fortran routines.... we need to pay attention of transfers between the two
 languages

 Naming convention between C and Fortran

 The name of the Fortran subroutine (foo) does not contain an underscore (_)
 In the C program, we must add ONE underscore to the Fortran name:
 the C program refers to the subroutine as:          foo_(...)
 while the Fortran code writes it as:
               foo(...)
 This is independent of the compiler pair (at least for gcc/f77 and
 Intel icc/ifort)
 ================================================================================================= */

#define F77Name(x) x##_   /* Need to add one underscore to Fortran program */

/* ===============================================================================================
   Prototypes for the different functions needed to generate a RandomWalk
================================================================================================== */

extern "C" void F77Name(eigendrv)(int *N, int *Neigen, double *mat,
         double *eigen, double *eigen0, double *eigenvect, double *WorkD, int *WorkI);

int find_Ncluster(int Neigen, double *eigenval, double Cutoff);

void SpectralClustering(int N, int Neigen, double *Laplacian, double *eigenvect, double *eigenval,
                double *WorkD, int *WorkI, double *Centers, double Cutoff, int iter, int Nrepeat,
                int *nc, int *ClusterID, MersenneTwister *mt, bool iprint);

extern "C" void F77Name(kmnsdrv)(int *M, int *N, double *A, int *K, double *C, double *W, int *IW, int *iter, double *score);

void KmeansSpectral(int N, int Ncoord, double *Data, int Ncluster, int iter, int Nrepeat, 
double *WorkD, int *WorkI, double *Centers, int *ClusterID, MersenneTwister *mt);

void init_cluster_center(int N, int Ncoord, double *Data, double *WorkD, int Ncluster, double *Centers, MersenneTwister *mt);

void random_cluster_center(int N, int Ncoord, double *Data, int *WorkI, int Ncluster, double *Centers, MersenneTwister *mt);
double Dist2Center(int N, int Ncoord, double *Data, int i, int Ncluster, double *Centers, int ic);

void Spectral_2_Clusters(int N, double *Data, int *ClusterID);
double GetMedian(double *Vect, int N);
double GetMean(double *Vect, int N);

#endif

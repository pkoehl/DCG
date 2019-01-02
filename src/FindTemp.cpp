/* ===============================================================================================
   FindTemp, part of DCG  (source code generated 2016-07-11)
   Copyright (c) Jiahui Guan.
   Copyright (c) Patrice Koehl.

   Contains a set of routines for finding the array of temperatures to be used in DCG:

	- FindTemp:	the main routine
	- DefineScale:  defines the "scale" of the data
	- FindTmax:	finds the highest temperature to be used (corresponds to only one cluster)
	- FindT0:	finds the lowest temperature to be used (corresponds to Nc clusters, where
			Nc is given as input to DCG)
	- BracketTemp	finds a temperature using a bracketing procedure

 =============================================================================================== */

#include "FindTemp.h"
#include "ProcessTemp.h"

using namespace std;

/* ===============================================================================================
   FindTemp(DistMat, N, Nwalks, P, Neigen, Ensemble, Q, Degree, WorkD, WorkI,
	Centers, Cutoff, iter, Nrepeat, Membership, ClusterID, Nmin, Nc, Temp, NTemp, nthreads, mt)

   What it does:
	Find the different temperatures to be used by DCG procedure

   Input:
	DistMat : distance matrix
	N       : number of data points (also number of rows and columns of DistMat)
	Nwalks  : parameter; number of steps in a random walk
	P       : parameter; number of random walks per point
	Ensemble: Work array, size N*N (will store Ensemble matrix)
	Q       : Work array, size N*N(will store transition probability matrix)
	Degree	: Work array, size N (will store degree of each point)
	WorkD	: Work array (double), for eigen computation and kmeans
	WorkI	: Work array (int), for eigen computation and kmeans
	Centers : Work array (double), for kmeans
	iter	: # of iterations for Kmeans algorithm
	Nrepeat : # of kmeans repeats
	Membership: Work array (bool, size N*N)
	ClusterID:  Work array (int, size N)
	Nmin	: minimal number of clusters at the lowest temperature
	nthreads: number of threads for parallel computing of random walk
	mt	: random number generators
    Output:
  	Ntemp	: number of temperatures
	Temp    : vector of temperatures
	Nc	: vector of # of clusters for each temperature

   =============================================================================================== */

void FindTemp(double *DistMat, int N, int M, int P, int Neigen, double *Ensemble, double *Q,
		double *Degree, double *WorkD, int *WorkI, double *Centers, double Cutoff, 
		int iter, int Nrepeat, bool *Membership, int *ClusterID, int Nmin, int *Nc, double *Temp,
		int *Ntemp, int nthreads, MersenneTwister **mt)
{

/* 	==========================================================================================
   	Find "scale" of Distance matrix
	========================================================================================== */

	double scale = DefineScale(N, DistMat);

/* 	=========================================================================================
   	Find Largest temperature: ProcessTemp should only find one cluster at that temperature
   	========================================================================================= */

	double Tmax = FindTmax(DistMat, N, M, P, Neigen, Ensemble, Q, Degree, WorkD, WorkI,
	Centers, Cutoff, iter, Nrepeat, Membership, ClusterID, scale, nthreads, mt);

/* 	=========================================================================================
   	Find smallest temperature: ProcessTemp should find at least Ncluster at that temperature
   	========================================================================================= */

	double T0 = FindT0(&Nmin, DistMat, N, M, P, Neigen, Ensemble, Q, Degree, WorkD, WorkI,
	Centers, Cutoff, iter, Nrepeat, Membership, ClusterID, scale, nthreads, mt);

/* 	=========================================================================================
   	Find temperatures between T0 and Tmax: find one temperature per cluster number, using a 
	bracketing procedure
   	========================================================================================= */

	*Ntemp = Nmin;
	Temp[0] = T0;
	Nc[0]   = Nmin;
	Temp[*Ntemp-1] = Tmax;
	Nc[*Ntemp-1] = 1;

	double T1 = T0;
	double T2 = Tmax;
	double Tk;
	int Ntarget;
	int ncluster;

	for(int k = 1; k < *Ntemp -1; k++) {
		Ntarget = k+1;
		Tk = BracketTemp(T1, T2, Ntarget, DistMat, N, M, P, Neigen, Ensemble,
			Q, Degree, WorkD, WorkI, Centers, Cutoff, iter, Nrepeat, Membership, 
			&ncluster,ClusterID, nthreads, mt);
		Temp[Nmin-k-1] = Tk;
		Nc[Nmin-k-1] = ncluster;
		T2 = Tk;
	}
}

/* ===============================================================================================
   DefineScale.cpp

    Identify "scale" of data, as the average of the min distance between points and their neighbours

    Input:
		N:	 number of data points
		DistMat: the distance matrix for the datapoints considered
    ========================================================================================== */

double DefineScale(int N, double *DistMat)
{
	double d,temp_min, sum; 

	sum = 0.0;
	for (int i = 0; i < N ; i ++)
	{
		temp_min=DistMat[i*N+0];
		for (int j =1; j < N ; j ++)
		{
			d = DistMat[i*N+j];
			if(d !=0 && d <temp_min) {
				temp_min= d;     
			}
		}
		sum = sum + temp_min;
	}
    
	sum = sum / N;

	return sum; 

}
/* ===============================================================================================
   FindTmax.cpp

   Finds the maximum temperature for DCG: at this temperature, the Random Walk procedure should
   only identify one cluster

   =============================================================================================== */

double FindTmax(double *DistMat, int N, int M, int P, int Neigen, double *Ensemble, double *Q,
		double *Degree, double *WorkD, int *WorkI, double *Centers, double Cutoff, 
		int iter, int Nrepeat, bool *Membership, int *ClusterID, double Tinit,
		int nthreads, MersenneTwister **mt)
{
	int ncluster = -1;
	double Tmax = Tinit/2.0;
	bool iprint=false;
	while (ncluster != 1) {
		Tmax = Tmax*2.0;
		ProcessTemp(DistMat, N, Tmax, M, P, Neigen, Ensemble, Q, Degree, WorkD, WorkI,
		 Centers, Cutoff, iter, Nrepeat, Membership, &ncluster, ClusterID, nthreads, mt, iprint);
	}

	return Tmax;
}

/* ===============================================================================================
   FindT0.cpp

   Finds the minimum temperature for DCG: at this temperature, the Random Walk procedure should
   identify at least Nmin clusters

   =============================================================================================== */

double FindT0(int *Nmin, double *DistMat, int N, int M, int P, int Neigen, double *Ensemble, double *Q,
		double *Degree, double *WorkD, int *WorkI, double *Centers, double Cutoff, 
		int iter, int Nrepeat, bool *Membership, int *ClusterID, double Tinit,
		int nthreads, MersenneTwister **mt)
{
	int ncluster = -1;
	double T0 = 4.0*Tinit/3.0;
	bool iprint=false;
	while (ncluster < *Nmin) {
		T0 = T0*3.0/4.0;
		ProcessTemp(DistMat, N, T0, M, P, Neigen, Ensemble, Q, Degree, WorkD, WorkI,
		 Centers, Cutoff, iter, Nrepeat, Membership, &ncluster, ClusterID, nthreads, mt, iprint);
	}

	*Nmin = ncluster;
	return T0;
}

/* ===============================================================================================
   BracketTemp
   Find the threshold temperatures using bracketing 
   =============================================================================================== */

double BracketTemp(double T1, double T2, int Ntarget, double *DistMat, int N, int M, int P, int Neigen, 
	double *Ensemble, double *Q,double *Degree, double *WorkD, int *WorkI, double *Centers, double Cutoff,
        int iter, int Nrepeat, bool *Membership, int *nc, int *ClusterID, int nthreads, MersenneTwister **mt)
{
	double tol = 0.005;
	double scale = (T2 - T1)/T2;
	double T_try;
	int ncluster;
	bool iprint = false;
	
	while (scale > tol) {

		T_try = (T2+T1)/2.0;

		ProcessTemp(DistMat, N, T_try, M, P, Neigen, Ensemble, Q, Degree, WorkD, WorkI,
		Centers, Cutoff, iter, Nrepeat, Membership, &ncluster, ClusterID, nthreads, mt, iprint);

		if(ncluster == Ntarget) {
			*nc = ncluster;
			return T_try;
		}
		else if (ncluster < Ntarget) {
			T2 = T_try;
		}
		else {
			T1 = T_try;
		}
		scale = (T2-T1)/T2;

	}

	*nc = ncluster;
	return T1;

}

/* ===============================================================================================

 ProcessTemp , part of DCG  (source code generated 2017-05-20)
 Copyright (c) Jiahui Guan.
 Copyright (c) Patrice Koehl.

================================================================================================== */

#include "ProcessTemp.h"
#include "RandomWalk.h"
#include "SpectralClustering.h"


using namespace std;

/* ===============================================================================================
 Function ProcessTemp

 What it does:
 This function reads in a distance matrix, perform multiple random walks on the corresponding
 network, derive an Ensemble matrix, finds the top eigenvalues of that matrix, use that information
 to derive a number of clusters for that matrix, and generate membership table

 Input:
	DistMat:	the distance matrix
	N      :	Number of row/col for DistMat (i.e. DistMat is a matrix of size NxN)
	Temp   :	The temperature considered
	M      :	Number of steps for a given random walk
	P      :	Number of repeated random walks for a given data point
	Neigen:		Number of eigenvalues of Ensemble matrix to compute
	Ensemble:	Storage for the Ensemble matrix
	Q       :	Storage for Transition matrix
	Deg	:	Storage for Degree of Transition matrix
	WorkD	:	Storage space (double)
	WorkI	:	Storage space (int)
	Centers:	Storage space (doubles)
	Cutoff  : 	parameter to decide on the number of clusters for that temperature
	iter:		# of iterations for Kmeans
	Nrepeat:	# or repeats for Kmeans
	nthread:	Number of threads to use for parallel computing
	mt:		the different random number generators
	iprint:		flag: if true, print information about eigenvalue
 Output
	Membership:	D array that defines if two nodes are in the same cluster or not,
			for that given temperature
	ClusterID:	Cluster # for each data point

================================================================================================== */

void ProcessTemp(double *DistMat, int N, double Temp, int M, int P,
		int Neigen, double *Ensemble, double *Q, double *Degree, double *WorkD,
		int *WorkI, double *Centers, double Cutoff, int iter, int Nrepeat, bool *Membership,
		int *nc, int *ClusterID, int nthread, MersenneTwister **mt, bool iprint)
{

/*  	==========================================================================================
    	Generate RandomWalks
    	========================================================================================== */

	RandomWalk(DistMat, Q, Degree, N, M, P, Temp, Ensemble, nthread, mt);
    
/*  	==========================================================================================
    	Laplacian normalization
    	========================================================================================== */

	NormalizeEnsemble(Ensemble,Degree,N); //Note that Degree is just a work array now

/*  	==========================================================================================
	Perform Spectral Clustering on the Ensemble matrix (that is now a Laplacian)

    	Note: we use Degree to store eigenvalues, and Q for eigenvectors, to save space
    	========================================================================================== */

	SpectralClustering(N, Neigen, Ensemble, Q, Degree, WorkD, WorkI, Centers, Cutoff,
		iter, Nrepeat, nc, ClusterID, mt[0], iprint);

/* 	==========================================================================================
   	Build membership table based on cluster assignment
   	========================================================================================== */

	for(int i = 0; i < N*N; i++) {
		Membership[i] = 0;
	}
	for(int i = 0; i < N; i++)
	{
		for(int j = i; j < N; j++)
		{
			if(ClusterID[i]==ClusterID[j])
			{
				Membership[N*i+j] = 1;
				Membership[N*j+i] = 1;
			}
		}
	}

}

/* ===============================================================================================
 Function NormalizeEnsemble(Ensemble, B, N)

 What it does:
	Normalize the Ensemble matrix to produce a positive semi-definite matrix:
	I - B^(-1/2) * E * B^{-1/2} where B is the diagonal matrix such that B(i,i)
	is the sum of the elements on row i of E.

 Input:
	Ensemble: ensemble matrix
	B       : work array
	N       : size of Ensemble (i.e. Ensemble is a matrix of size NxN
 Output:
	Ensemble: normalized ensemble matrix

================================================================================================== */
void NormalizeEnsemble(double *Ensemble, double *B, int N)
{

/*  	==========================================================================================
    	Find array B (we only store the diagonal)
    	========================================================================================== */

	double temp;
	for(int i=0; i< N; i++)
	{
		temp=0;
		for(int j=0; j<N; j++)
		{
			temp = temp + Ensemble[N*i+j];
		}
		B[i]=1/sqrt(temp);
	}
    
/*  	==========================================================================================
    	Normalize: E = I - B^(-1/2) * E * B^(-1/2)
    	========================================================================================== */

	int idx;
	for(int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			idx = N*i+j;
			Ensemble[idx]=-B[i]*B[j]*Ensemble[idx];
		}
	}

	for(int i=0; i<N; i++)
	{
		idx = N*i+i;
		Ensemble[idx]=Ensemble[idx] + 1;
	}

}


/* ===============================================================================================
   Function NormalizeEnsemble2(Ensemble, B, N)

   What it does:
	Normalize the Ensemble matrix to produce a positive semi-definite matrix:
	B^(-1/2) * E * B^{-1/2} where B is the diagonal matrix such that B(i,i)
	is the sum of the elements on row i of E.

   Input:
	Ensemble: ensemble matrix
	B       : work array
	N       : size of Ensemble (i.e. Ensemble is a matrix of size NxN
   Output:
	Ensemble: normalized ensemble matrix
    ========================================================================================== */
   
void NormalizeEnsemble2(double *Ensemble, double *B, int N)
{

/*  	==========================================================================================
    	Find array B (we only store the diagonal)
    	========================================================================================== */

	double temp;
	int idx;
	for(int i=0; i< N; i++)
	{
		temp=0;
		for(int j=0; j<N; j++)
		{
			idx = i*N+j;
			temp = temp + Ensemble[idx];
		}
		B[i]=1.0/sqrt(temp);
	}
    
/*  	==========================================================================================
    	Normalize: E = B^(-1/2) * E * B^(-1/2)
    	========================================================================================== */

	for(int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			idx = N*i+j;
		    	Ensemble[idx]=B[i]*B[j]*Ensemble[idx];
		}
	}

}

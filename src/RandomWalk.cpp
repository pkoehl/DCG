/* ===============================================================================================
 RandomWalk.cpp part of DCG  (source code generated 2017-05-21)
 Copyright (c) Jiahui Guan.
 Copyright (c) Patrice Koehl.
   ================================================================================================== */

#include "RandomWalk.h"
#include "mt.h"

using namespace std;

/* ===============================================================================================
 Prepare structure to pass data for multi-threaded Random Walk generation

 The structure is defined as jobs_data, and contains:

 SM:            similarity matrix
 N:             number of nodes
 N1:		first node considered
 N2:		last node considered
 Nsteps:        number of steps for a given random walk
 Nwalks:	number of repeats for a given node
 NeachThread:   Number of MC random walk(iteration) in each thread
 *Ensemble:     the Ensemble matrix.
 ================================================================================================= */

#define NUM_THREADS 16
pthread_t threads[NUM_THREADS];

typedef struct jobs_data {
        double *SM;
        double *Degree;
        int N;
        int N1;
        int N2;
        int Nsteps;
        int Nwalks;
        double *Ensemble;
        MersenneTwister *mt;
} jobs_data;

jobs_data jobs[NUM_THREADS];

int threadids[NUM_THREADS];


void RandomWalk(double *DistMat, double *Q, double *Degree, int N, int M, int P, double T, 
double *Ensemble, int nthreads, MersenneTwister **mt)
{
/* 	==========================================================================================
   	Compute the transition matrix
   	========================================================================================== */

	TransMatrix(N, T, DistMat, Q, Degree);

/* 	==========================================================================================
   	Perform random walks (in parallel)
   	========================================================================================== */

	int nval = N/nthreads;
	int N1,N2;

	for (int i=0; i < nthreads; i++)
	{
		N1 = i*nval;
		N2 = N1+nval;
		if(i == nthreads-1) N2 = N;

		threadids[i] = i;

		jobs[i].SM = Q;
		jobs[i].Degree = Degree;
		jobs[i].N  = N;
		jobs[i].N1  = N1;
		jobs[i].N2  = N2;
		jobs[i].Nsteps  = M;
		jobs[i].Nwalks = P;
		jobs[i].Ensemble = &Ensemble[N1*N];
		jobs[i].mt = mt[i];

		pthread_create(&threads[i], NULL, worker_thread, (void*) &threadids[i]);
        }

/*  	==========================================================================================
    	Join all the threads (to make sure they are all finished)
    	========================================================================================== */

	for (int i=0; i < nthreads; i++)
	{
		pthread_join(threads[i], NULL);
	}

/*  	==========================================================================================
    	Scale resulting Ensemble matrix
    	========================================================================================== */

	for (int i=0; i<N*N; i++)
	{
		Ensemble[i] = Ensemble[i]/M;
	}

/*  	==========================================================================================
    	Symmetrize the matrix 
    	========================================================================================== */

	SymmetricEnsemble(Ensemble,N); 	

}

/* ===============================================================================================
   Procedure to convert Distance Matrix into a Transition Probability Matrix

   Let D be the original distance matrix.
   We first compute the transition matrix W as:
	W(i,j) = exp( -D(i,j) / T) for i != j, and W(i,i) = 0
   Then we normalize W based on the degree of each node:
	Deg(i) = Sum_{i=1}^N W(i,j)
   and
	Q(i,j) = W(i,j)/Deg(i)

   To save space, the matrix Q overwrites the matrix W

   Input:
	Ndim: size of the matrix
	Temp: temperature considered
	DistMat: the original distance matrix
   Output:
	Q: the transition probability matrix
	Deg: the degree of each node

   =============================================================================================== */

void TransMatrix(int N, double Temp, double *DistMat, double *Q, double *Degree)
{

/* 	==========================================================================================
	Compute transition matrix
  	========================================================================================== */

	int idx, idx1, idx2;
	double x;

	for(int i=0; i< N; i++)
	{
		idx = N*i+i;
		Q[idx] = 0;
		for(int j=i+1; j< N; j++)
		{
			idx1 = i*N+j;
			idx2 = j*N+i;
			x = DistMat[i*N+j]/Temp;
			Q[idx1] = exp(-x);
//			Q[idx1]= exp(-x*x);
			Q[idx2]=Q[idx1];
		}
	}

/* 	==========================================================================================
	Compute degree of each vertex (node)
   	========================================================================================== */

	for(int i=0;i < N; i++)
	{
		double temp=0;
		for(int j=0; j< N; j++)
		{
			idx = N*i + j;
			temp=temp + Q[idx];
		}
		Degree[i]=temp;
	}

/* 	==========================================================================================
	Convert Q into a transition probability matrix (i.e. normalize with degree)
   	========================================================================================== */

	for (int i=0;i< N; i++)
	{
		for (int j=0; j< N; j++)
		{
			idx = N*i + j;
			Q[idx]=Q[idx]/Degree[i];
		}
	}
}


/* ===============================================================================================
 Procedure	void sampleGenerator(N, P, cdf, Idx, mt)

 What it does:
	Choose a node based on a given probability distribution function
	1. Computes the cdf (cumulative distribution function) given pdf
	2. Generate a random number uniformly between 0 and 1
	3. Find in which interval it belongs
	4. Pick the node ID accordingly

 Input:
	N: number of nodes considered
	P:  Probability distribution
 Output:
	new node

================================================================================================== */

int sampleGenerator(int N, double *P, MersenneTwister *mt)
{
/* 	==========================================================================================
   	Generate a random number between 0 and 1
   	========================================================================================== */

	double rnd = mt->genrand_res53();

/* 	==========================================================================================
   	Scan the cdf
   	========================================================================================== */

	int i;
	double cdf=P[0];
	if(cdf > rnd) {
		return 0;
	}
	for(i=1; i<N;i++)
	{
		cdf=cdf+P[i];
		if(cdf > rnd) {
			return i;
		}
	}
	i = (int) (N*rnd);
	if(i >(N-1)) i = N-1;
	return i;

}

/* ===============================================================================================
 Function void* worker_thread

 What it does:
	Multi-threaded version of the random walks on a given network.
	Each thread should execute a subset of the random walks, for a subset of
	the Nodes
================================================================================================== */
void* worker_thread(void* data)
{

/* 	==========================================================================================
	Get thread id and collect information from corresponding structure jobs
   	========================================================================================== */

	int threadid = *((int *) data);

	int Nwalks, Nsteps, N, N1, N2, nval;
	int node;

	N = jobs[threadid].N;
	N1 = jobs[threadid].N1;
	N2 = jobs[threadid].N2;
	Nwalks = jobs[threadid].Nwalks; 
	Nsteps = jobs[threadid].Nsteps; 
	nval = N2-N1;

	for(int k=0; k< N*nval; k++)
	{
		jobs[threadid].Ensemble[k]=0;
	}
  
/* 	==========================================================================================
   	Perform the random walks and accumulate the corresponding Ensemble matrix
   	========================================================================================== */

	for (int k=0; k < Nwalks; k++)
	{
		for (int i = 0; i < nval; i++)
		{
			node = i+N1;
			for (int j = 0; j < Nsteps; j++)
			{
		    		node = sampleGenerator(N, &jobs[threadid].SM[N*node], jobs[threadid].mt);
		    		jobs[threadid].Ensemble[i*N+node]=jobs[threadid].Ensemble[i*N+node]+1; 
		    }
		}
	}

	return 0;
}

/* ===============================================================================================
 Procedure	void SymmetricEnsemble (double *Ensemble, int N)

 What it does:
    Symmetrize the matrix by adding E[i][j]+E[j][i] at the entry of (i,j)

 Input:
	Ensemble: A non-symmetric matrix
	N:  Its dimension: N by N 
 Output:
	Ensemble: A symmetric matrix

================================================================================================== */

void SymmetricEnsemble (double *Ensemble, int N)
{
	double temp; 
	int idx1,idx2;

	for (int i=0; i<N; i++)
	{
		for (int j=i+1; j<N; j++)
		{
			idx1 = i*N+j;
			idx2 = j*N+i;
			temp = (Ensemble[idx1]+Ensemble[idx2])/2.0;
			Ensemble[idx1]=temp;
			Ensemble[idx2]=temp;
		}
	}
}

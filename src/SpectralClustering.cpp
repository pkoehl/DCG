/* ===============================================================================================

 SpectralClustering, part of DCG  (source code generated 2017-05-21)
 Copyright (c) Jiahui Guan.
 Copyright (c) Patrice Koehl.

================================================================================================== */

#include "SpectralClustering.h"

using namespace std;


/* ===============================================================================================
 Procedure SpectralClustering (N, NEigen, Laplacian, eigenvect, eigenval,
	WorkD, WorkI, Centers, Cutoff, iter, Nrepeat, Ncluster, ClusterID, mt, iprint)

 What is does:
 Performs spectral clustering on a Laplacian matrix

 Input:
	N: 	  	number of data points
	Neigen:		number of eigenvalues computed for the Laplacian
	Laplacian:    	Laplacian matrix
	Workd		an array of doubles, used as work space
	WorkI		an array of int, used as work space
	Centers		centers of the clusters (work array)
	Cutoff:		cutoff to decide to keep an eigenvalue
	iter		number of iterations for kmeans (usually 30)
	Nrepeat 	number of times Kmeans is repeated
	mt		pointer to a Mersenne Twister random number generator
	iprint:		flag for printing
 Output:
	Ncluster	Number of clusters
	ClusterID	for each data points, cluster it is assigned to

================================================================================================== */
void SpectralClustering(int N, int Neigen, double *Laplacian, double *eigenvect, double *eigenval,
		double *WorkD, int *WorkI, double *Centers, double Cutoff, int iter, int Nrepeat,
		int *nc, int *ClusterID, MersenneTwister *mt, bool iprint)
{
/*  	==========================================================================================
    	Compute the Neigen largest eigenvalue of the Ensemble matrix using ARPACK
    	compute eigenvalues using ARPACK
		- find first largest eigenvalue, eigen0
		- then compute the Neigen smallest eigenvalues, stored in eigenval
		- finally, compute eigenval(i) = 1 - eigenval(i)/eigen0
    	========================================================================================== */

	double eigen0;
	F77Name(eigendrv)(&N, &Neigen, Laplacian, eigenval, &eigen0, eigenvect, WorkD, WorkI);

	if(iprint) {
		cout << " " << endl;
		cout << "            Number of eigenvalues computed: " << Neigen << endl;
		cout << "            Scaled Eigenvalues            : " ;
		for(int i=0; i< Neigen; i++){
			cout << eigenval[i] << " ";
		}
		cout << endl;
		cout << "eigen0 = " << eigen0 << endl;
	}


	for(int i=0; i< Neigen; i++){
		eigenval[i]= 1-eigenval[i]/eigen0;
	}

	if(iprint) {
		cout << " " << endl;
		cout << "            Number of eigenvalues computed: " << Neigen << endl;
		cout << "            Scaled Eigenvalues            : " ;
		for(int i=0; i< Neigen; i++){
			cout << eigenval[i] << " ";
		}
		cout << endl;
	}

/*  	==========================================================================================
    	Decide on the number of clusters: "plot" scaled eigenvalues, and count how many are above
    	a threshold.
    	Note that find_Ncluster will destroy the vector eigenval
    	========================================================================================== */

	int Ncluster = find_Ncluster(Neigen, eigenval, Cutoff);
	if(iprint) {
		cout << "            Number of clusters            : " << Ncluster << endl;
	}

/* 	==========================================================================================
	Assign points to cluster using Kmeans clustering
   	========================================================================================== */

	if(Ncluster == 1) {
		for (int i = 0; i < N; i++) {
			ClusterID[i] = 1;
		}
	}
//	else if(Ncluster == 2) {
//		Spectral_2_Clusters(N, eigenvect, ClusterID);
//	}
	else {
		KmeansSpectral(N, Ncluster, eigenvect, Ncluster, iter, Nrepeat, WorkD, WorkI,
		Centers, ClusterID, mt);
	}
	(*nc) = Ncluster;

}
/* ===============================================================================================
 Procedure KmeansSpectral (N, Eigen, Ncluster)

 What is does:
 Perform kmeans clustering on a set of data points based on the eigenvectors of the ensemble matrix
 computed for those data points.

 Input:
	N: 	  	number of data points
	Ncoord:		number of features describing data points
	Data:    	Data matrix, with points as columns, and features as rows
	Ncluster: 	number of clusters
	iter		number of iterations for kmeans (usually 30)
	Nrepea 		number of times Kmeans is repeated
	Workd		an array of doubles, used as work space
	WorkI		an array of int, used as work space
	mt		pointer to a Mersenne Twister random number generator
 Output:
	ClusterID	for each data points, cluster it is assigned to
	Centers		centers of the clusters

================================================================================================== */

void KmeansSpectral(int N, int Ncoord, double *Data, int Ncluster, int iter, int Nrepeat, 
	double *WorkD, int *WorkI, double *Centers, int *ClusterID, MersenneTwister *mt)
{ 

	double score, best;

/* 	==========================================================================================
   	First clustering is based on Kmeans++ procedure of Arthur and Vassilvitskii
   	========================================================================================== */

	init_cluster_center(N, Ncoord, Data, WorkD, Ncluster, Centers, mt); // initialize cluster centers

	F77Name(kmnsdrv)(&N, &Ncoord, Data, &Ncluster, Centers, WorkD, WorkI, &iter, &score); // performs kmeans

   	best=score; 
	for(int i = 0; i < N; i++) {
		ClusterID[i] = WorkI[i];
	}

/* 	==========================================================================================
   	Now repeat Nrepeat-1 times, either with new Kmeans++ assignment, or pure random...
   	========================================================================================== */

	int type = 1; 	// 1: Kmeans++; 0: random

	for(int k = 0; k < Nrepeat-1 ; k++)
	{
		if(type==1) {
			init_cluster_center(N, Ncoord, Data, WorkD, Ncluster, Centers, mt);
		}
		else {
			random_cluster_center(N, Ncoord, Data, WorkI, Ncluster, Centers, mt);
		}
		F77Name(kmnsdrv)(&N, &Ncoord, Data, &Ncluster, Centers, WorkD, WorkI, &iter, &score);

		if(score < best) {
			best = score;
			for(int i = 0; i < N; i++) {
				ClusterID[i] = WorkI[i];
			}
		}
	}
}

/* ===============================================================================================
  Procedure init_cluster_center(N, Ncoord, Data, D, Ncluster, Centers, mt2);

 What is does:
	Initialize kmeans cluster by picking Ncluster Centers, using the procedure kmeans++
	of Arthur and Vassilvitskii:
	Arthur, D.; Vassilvitskii, S. (2007). "k-means++: the advantages of careful seeding"
	Proceedings of the eighteenth annual ACM-SIAM symposium on Discrete algorithms. 
	Society for Industrial and Applied Mathematics Philadelphia, PA, USA. pp. 1027–1035

 Input:
	N: 	  number of data points
	Ncoord:	  number of coordinates per Data point
	Data:     Data matrix, needed to compute distances; note that rows are coord, columns are points
	D:	  work array (double); size N
	Ncluster: number of clusters
	mt2:	  the Mersenne twister number generator
 Output:
	Centers: positions of the initial centers

================================================================================================== */

void init_cluster_center(int N, int Ncoord, double *Data, double *Dist, int Ncluster, double *Centers, 
		MersenneTwister *mt2)
{
/*	==========================================================================================
	Randomly pick a point to be the first center
	========================================================================================== */

	double rnd = mt2->genrand_res53();

	int ipick = (int) (rnd*N);
	for(int i = 0; i < Ncoord; i++)
	{
		Centers[i*Ncluster] = Data[ipick+i*N];
	}

	int icluster = 1;

/*	==========================================================================================
	Apply kmeans++ procedure
	========================================================================================== */

	double x, Dmin;
	double Sum;
	double d;

	while(icluster < Ncluster)
	{
/*		==================================================================================
		Compute Dist(i), the minimal distance between data point i and all existing 
		Cluster centers
		================================================================================== */

		Sum = 0;
		for(int i = 0; i < N; i++)
		{
			for(int ic = 0; ic < icluster; ic++)
			{
				d = Dist2Center(N, Ncoord, Data, i, Ncluster, Centers, ic);
				if(ic==0) {
					Dmin = d;
				}
				else {
					Dmin = min(Dmin,d);
				}
			}
			Dist[i]=Dmin;
			Sum = Sum + Dist[i];
		}

/*		==================================================================================
		Convert Dist(i) into a cumulative distribution function
		Pick random point based on this cumulative distribution function
		================================================================================== */

		for(int i = 0; i < N; i++)
		{
			Dist[i] = Dist[i]/Sum;
		}

		ipick = 0;
		rnd = mt2->genrand_res53();
		for(int i = 1; i < N; i++)
		{
			Dist[i] = Dist[i] + Dist[i-1];
			if(Dist[i] > rnd) {
				ipick = i-1;
				break;
			}
		}

/*		==================================================================================
		Define cluster centers
		================================================================================== */

		for(int i = 0; i < Ncoord; i++)
		{
			Centers[icluster + Ncluster*i] = Data[ipick+i*N];
		}

		icluster = icluster + 1;
	}
}

/* ===============================================================================================
  Procedure random_cluster_center(N, Ncoord, Data, D, Ncluster, Centers, mt2);

 What is does:
	Initialize kmeans cluster by picking Ncluster Centers at random

 Input:
	N: 	  number of data points
	Ncoord:	  number of coordinates per Data point
	Data:     Data matrix, needed to compute distances; note that rows are coord, columns are points
	Ncluster: number of clusters
	mt2:	  the Mersenne twister number generator
 Output:
	Centers: positions of the initial centers

================================================================================================== */

void random_cluster_center(int N, int Ncoord, double *Data, int *Ipick, int Ncluster, double *Centers, 
		MersenneTwister *mt2)
{
/*	==========================================================================================
	Randomly pick a point to be the first center
	========================================================================================== */

	double rnd;
	int j, ipck;

	for (int i = 0; i < Ncluster; i++) {
		j = 1;
		while (j > 0) {
			j = 0;
			rnd = mt2->genrand_res53();
			ipck = (int) (rnd*N);
			for(int k = 0; k < i; k++) {
				if(Ipick[k] == ipck) j = 1;
			}
		}
		Ipick[i] = ipck;
	}

	for(int k = 0; k < Ncluster; k++) {
		for(int i = 0; i < Ncoord; i++)
		{
			Centers[k+i*Ncluster] = Data[Ipick[k]+i*N];
		}
	}
}

/* ===============================================================================================
   Procedure Dist2Center(N, Ncoord, Data, idx, Ncluster, Centers, jdx)

 What is does:
	Compute the distance between one point, and one of the cluster centers

 Input:
	N	: number of datapoints
	Ncoord	: number of coordinates
	Data	: Data matrix
	idx	: data point considered
	Ncluster: number of clusters
	Centers : coordinates of the centers
	jdx	: cluster center considered
 Output:
	distance (squared) between point idx and center jdx

================================================================================================== */

double Dist2Center(int N, int Ncoord, double *Data, int idx, int Ncluster, double *Centers, int jdx)
{
	double dist = 0;
	double x;
	for (int i = 0; i < Ncoord; i++)
	{
		x = Data[idx + i*N] - Centers[jdx+i*Ncluster];
		dist = dist + x*x;
	}
	return dist;
}

/* ===============================================================================================
 Procedure Spectral_2_Clusters (N, Data, ClusterID)

 What is does:
 Separate data points into two clusters, based on largest eigenvector of the 
 normalized ensemble matrix computed for those data points.

 Input:
	N: 	  	number of data points
	Data:    	eigenvectors considered
 Output:
	ClusterID	cluster assignments for the points

================================================================================================== */

void Spectral_2_Clusters(int N, double *Data, int *ClusterID)
{ 

/*	==========================================================================================
	Find mean of coord 1
	========================================================================================== */

	double mean = GetMean(Data,N);

/*	==========================================================================================
	Assign each point to a cluster based on its position with respect to mean
	========================================================================================== */

	for (int i = 0; i < N; i++)
	{
		if(Data[i] < mean) {
			ClusterID[i] = 1;
		}
		else {
			ClusterID[i] = 2;
		}
	}

}

/* ===============================================================================================
 Procedure GetMedian(Vect, N)

 What is does:
 Get the median value from an array of unsorted values

 Input:
	Vect:	array considered
	N:	number of points in the array
 Output:
	median

================================================================================================== */

double GetMedian(double Array[], int Size) 
{

/*	==========================================================================================
	Allocate a vector of the same size and sort it.
	========================================================================================== */

	double *sorted = new double[Size];
	for(int i = 0; i < Size; i++)
	{
		sorted[i] = Array[i];
	}
	sort(sorted,sorted+Size);

/*	==========================================================================================
	Get median from sorted array
	========================================================================================== */

	double Median = 0.0;
	if ((Size % 2) == 0) {
		Median = (sorted[Size/2] + sorted[(Size/2) - 1])/2.0;
	} 
	else {
		Median = sorted[Size/2];
	}

	delete [] sorted;

	return Median;
}

/* ===============================================================================================
 Procedure GetMean(Vect, N)

 What is does:
 Get the mean value from an array of unsorted values

 Input:
	Vect:	array considered
	N:	number of points in the array
 Output:
	mean

================================================================================================== */

double GetMean(double Array[], int Size) 
{

/*	==========================================================================================
	Get mean from array
	========================================================================================== */

	double Mean = 0.0;
	for(int i = 0; i < Size; i++) {
		Mean = Mean + Array[i];
	}
	Mean = Mean/Size;

	return Mean;
}
/* ===============================================================================================
 Procedure find_Ncluster(Neigen, eigenval, Cutoff)

 What it does:
	Find number of clusters for the given temperature, based on spectrum of eigenvalues

 Input:
	Neigen:	number of eigenvalues
	eigenval: eigenvalues, ranked from the largest to the smallest
 Output:
	Ncluster: estimated number of cluster

   =============================================================================================== */

int find_Ncluster(int Neigen, double *eigenval, double Cutoff)
{

	double small_threshold= 1.60922e-17; 
    
/* 	===============================================================================================
   	To improve discrimination, replace eigenval with log(eigenval)
   	=============================================================================================== */

	for(int k=0; k< Neigen; k++)
	{
		if(eigenval[k] < small_threshold) {
			eigenval[k] = 35;
		}
    		else {
    			eigenval[k] = abs(log(eigenval[k]));
    		}
	}

/* 	===============================================================================================
   	Find first index for which eigenvalue is above threshold
   	=============================================================================================== */

	int i=0;
	while(eigenval[i]<=Cutoff) i++;

	return i; 
}


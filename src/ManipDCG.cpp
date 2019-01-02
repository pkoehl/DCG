/* ===============================================================================================

   ManipDCG, part of DCG  (source code generated 2017-05-20)
 Copyright (c) Jiahui Guan.
 Copyright (c) Patrice Koehl.
================================================================================================== */

#include "ManipDCG.h"

/* ===============================================================================================
Procedure:    UpdateDCG(N, DCG, Status, Member, Temp)

What it does:
updates current DCG matrix: for any pair of points that are still "alife" (with Status true),
if they do not belong to the same cluster at the current temperature, then set their distance
to the temperature Temp, and set their status to false

Input:
	N	Number of data points
	DCG	Current DCG matrix
	Status	Status of each pair of points
	Member	membership at current temperature
	Temp	previous temperature
Output:
	DCG	New DCG matrix
	Status	Updated Status of each pair of points
================================================================================================== */

void UpdateDCG(int N, double *DCG, bool *Status, bool *Member, double Temp) 
{

	int ipos;
	for(int i=0; i<N; i++)
	{
		for(int j=i+1; j< N; j++)
		{
			ipos = i*(N-1) - (i*(i+1))/2 + j - 1;
			if(Status[ipos]==true) {
				if(Member[N*i+j]==false) {
					DCG[ipos] = Temp;
					Status[ipos] = false;
				}
			}
		}
	}
}

/* ===============================================================================================
Procedure:    FinalDCG(N, DCG, Status, Temp)

What it does:
Finalize DCG matrix: for any pair of points that are still "alife" (with Status true),
set their distance to the temperature Temp;

Input:
	N	Number of data points
	DCG	Current DCG matrix
	Status	Status of each pair of points
	Temp	previous temperature
Output:
	DCG	New DCG matrix
	Status	Updated Status of each pair of points
================================================================================================== */

void FinalDCG(int N, double *DCG, bool *Status, double Temp) 
{

	int ipos;
	for(int i=0; i<N; i++)
	{
		for(int j=i+1; j< N; j++)
		{
			ipos = i*(N-1) - (i*(i+1))/2 + j - 1;
			if(Status[ipos]==true) {
				DCG[ipos] = Temp;
				Status[ipos] = false;
			}
		}
	}
}

/* ===============================================================================================
Procedure:    CopyDCG(N, DCG, Mat)

What it does:
Copies DCG matrix in compressed format into full DCG matrix

Input:
	N	Number of data points
	DCG	Compressed DCG matrix
Output:
	Mat	Uncompressed DCG matrix
================================================================================================== */

void CopyDCG(int N, double *DCG, double *Mat) 
{

	int ipos, idx, jdx;
	for(int i=0; i<N; i++)
	{
		for(int j=i+1; j< N; j++)
		{
			ipos = i*(N-1) - (i*(i+1))/2 + j - 1;
			idx = N*i+j;
			jdx = N*j+i;
			Mat[idx] = DCG[ipos];
			Mat[jdx] = DCG[ipos];
		}
		Mat[i*N+i] = 0.0;
	}
}

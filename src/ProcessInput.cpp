/* ===============================================================================================

   ProcessInput contains:

	- count_data: 		counts number of subject /features contain in the input file
	- read_data:  		reads the data
	- HammingDist: 		computes Hamming distance for binary data
	- EuclideanDist: 	computes Euclidean distance for continuous data

   Copyright (c) Jiahui Guan.
   Copyright (c) Patrice Koehl.

   ================================================================================================== */

#include "ProcessInput.h"

/* ===============================================================================================
   Procedure to counting the number of rows and columns in the data matrix
   =============================================================================================== */

bool count_data(string filename, int *row, int *col)
{
	string line;

	ifstream fileIN;
	fileIN.open(filename.c_str());

/* 	==========================================================================================
   	First check that file exists!
   	========================================================================================== */

	if(fileIN.fail())
	{
		cerr <<"The DATA file you are trying to access cannot be found \n";
		cout <<"Try again\n";
		return false;
	}

/* 	==========================================================================================
   	Go over each line in the file
   	========================================================================================== */

	double x;
	int nval=0;
	*row = 0;

	while(getline(fileIN, line)) // until reach the end of file
	{
		if (line.substr(0,2) != "//")
		{
    			(*row)++;
			if(nval==0)
			{
    				istringstream stream(line);
    				while(stream >>x)
				{
    					nval++;
    				}
    			}
		}
	}
	*col = nval;

/* 	==========================================================================================
   	Close the file and return
   	========================================================================================== */

	fileIN.close();

	return true;

}

/* ===============================================================================================
   Procedure to read in data matrix 
   =============================================================================================== */

void read_data(string filename, double **DATA)
{
	string line;

	ifstream fileIN;
	fileIN.open(filename.c_str());

/* 	==========================================================================================
	No need to check if file exists: was done in count_data
   	Go over each line in the file
   	========================================================================================== */

	double x;
	int row=0;
	int col;

	while(getline(fileIN, line))
	{
		if (line.substr(0,2) != "//")
		{
    			istringstream stream(line);
    			col=0;
    			while(stream >>x){
				DATA[row][col]=x;
    				col++;
    			}
    			row++;
    		}
	}

/* 	==========================================================================================
   	Close the file and return
   	========================================================================================== */

	fileIN.close();

}


/* ===============================================================================================
   Procedure to compute the Hamming distance for a binary matrix:

   DATA is the input data matrix, with "row" rows (the subjects), and "col" columns 
   (the binary features). The corresponding Hamming distance is computed on the columns 
   of DATA. For two indices (subjects)  i and j, the Hamming distance D(i,j) is computed
   by considering the two columns C(i) and C(j) of DATA, and counting the number
   of positions in which they differ

   As matrices are stored row-wise in C++, we loop first on the rows, then on the columns

   Note that the Distance Matrix DistMat is expected to have been initialized with 0 in the calling
   program
   =============================================================================================== */

void HammingDist(int row, int col, double **DATA, double *DistMat)
{
/* 	==========================================================================================
   	Compute Hamming distance between two columns C(i) and C(j) with i < j
   	========================================================================================== */

	for (int k = 0; k < row; k++)
	{
		for (int i=0; i< col; i++)
		{
			for (int j=i; j <col;j++)
			{
				if(DATA[k][i]!= DATA[k][j])
				{
					DistMat[row*i+j] = DistMat[row*i+j]+1;
				}
			} 
		}
	}

/* 	==========================================================================================
	Make DistMat symmetric
   	========================================================================================== */

	for (int i=1; i< col; i++)
	{
		for (int j=0; j <i;j++)
		{
			DistMat[row*i+j]=DistMat[row*j+i];
		}
	}
}

/* ================================================================================================
   Calculate the Euclidean distance between each point in dataset, store results in distance matrix
   ================================================================================================ */

void EuclideanDist(int row, int col, double **DATA, double *DistMat)
{
	int i, j, k;
	double x, sum, d;

	for (i = 0; i < row; i++)
	{
		for (j = i+1; j < row; j++)
		{
			sum = 0;
			for(k = 0; k < col; k++)
			{
				x = (DATA[i][k] - DATA[j][k]);  
				sum = sum + x*x;
			}
			d = sqrt(sum);
			DistMat[row*i+j] = d;
			DistMat[row*j+i] = d;
		}
		DistMat[row*i+i] = 0;
	}

} 

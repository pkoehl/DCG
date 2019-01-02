/* ===============================================================================================

DCG  (source code generated 2017-05-20)
Copyright (c) Jiahui Guan.
Copyright (c) Patrice Koehl.

>>> SOURCE LICENSE >>>

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

>>> END OF LICENSE >>>

================================================================================================== */

#include "DCG.h"
#include "ManipDCG.h"
#include "ProcessInput.h"
#include "ProcessTemp.h"
#include "FindTemp.h"
#include "mt.h"

using namespace std;

/* ===============================================================================================
 Main 
 Input: A data matrix (continuous), a binary matrix or a distance matrix (symmetric)
 Output: a DCG distance matrix (symmetric)
 ================================================================================================= */

int main(int argc, char *argv[])
{

/* 	==========================================================================================
   	Show usage if needed
   	========================================================================================== */

	int narg = argc;

	if( argc < 2 )
	{
		usage(argv);
		return -1;
	}

	string input = argv[1];
	if( input == "-h" || input == "-help" )
	{
		usage(argv);
		return -1;
	}

/* 	==========================================================================================
   	Read in all inputs (some values may have been preset)
   	========================================================================================== */

	string fileIN;
	string fileOUT;
	string fileLOG;
	string fileDCG;

	int data_type=2;
	int nthreads=1;

	int Nmin=10; 
	int M=1000;
    	int P=10; 

	double Cutoff=0.2;	// Cutoff for defining meaningful eigenvalue
	
	if (!parse_args(argc, argv, &fileIN, &fileOUT, &Nmin, &M, &P, &Cutoff, &nthreads)) return 1; 


	int Neigen = Nmin;

	fileLOG = fileOUT;
	fileLOG.append(".log");
	fileDCG = fileOUT;
	fileDCG.append(".res");

	ofstream logfile(fileLOG.c_str());

/* 	==========================================================================================
	Detect automatically data type based on file extension:
	 if ext is .crd, data_type = 2
	 if ext is .dist, data_type = o
	 if ext is .bin, data_type = 1
   	========================================================================================== */

	if(fileIN.substr(fileIN.find_last_of(".") + 1) == "crd") {
		data_type = 2;
	}
	else if(fileIN.substr(fileIN.find_last_of(".") + 1) == "bin") {
		data_type = 1;
	}
	else if(fileIN.substr(fileIN.find_last_of(".") + 1) == "dist") {
		data_type = 0;
	}
	else {
		cout << " Could not identify data type!" << endl;
		return -1;
	}

/* 	==========================================================================================
   	Read data matrix in two passes: first get size (number of rows and number of columns),
	then read into matrix DATA

        If data matrix is actually a distance matrix, row should then be equal to col
	If data matrix contains subjects / features, then row is the number of subjects, and
	col is the number of features for each subjects
   	========================================================================================== */

	int row = 0;
	int col = 0;
	if(!count_data(fileIN, &row, &col)) return 1;

	Neigen = min(Neigen, row/2);

	string *types = new string[3];
	types[0] = "Distance matrix";
	types[1] = "Binary data";
	types[2] = "Continuous data";

	cout << " " << endl;
	cout << "Data file: " << fileIN << endl;
	cout << "Data type: " << types[data_type] << endl; 
	cout << "Data matrix size: nrow = " << row << " ncol = " << col << endl;
	cout << " " << endl;

	logfile << " " << endl;
	logfile << "Input file:" << endl;
	logfile << "===========" << endl;
	logfile << "Data file            : " << fileIN << endl;
	logfile << "Data type            : " << types[data_type] << endl; 
	logfile << "Number of data points: " << row << endl;
	logfile << "Number of features   : " << col << endl;
	logfile << "Initial distance     : Euclidean " << endl;
	logfile << " " << endl;
	logfile << "Parameters:" << endl;
	logfile << "===========" << endl;
	logfile << "Maximum number of clusters    : " << Neigen << endl;
	logfile << "Number of random walks / point: " << P << endl;
	logfile << "Number of steps / random walk : " << M << endl;
	logfile << "Cutoff for selecting clusters : " << Cutoff << endl;
	logfile << " " << endl;

	double **DATA = new double*[row];
	for(int i = 0; i < row; i++)
	{
		DATA[i] = new double[col];
	}
	read_data(fileIN, DATA);

/* 	==========================================================================================
   	Generate distance matrix: if data_type = 1, nothing to do but copy array, otherwise compute
	Hamming distance for binary data or Euclidean distance for continuous data
   	========================================================================================== */

	double *DistMat = new double[row*row];

	if (data_type==0)
	{
		for (int i=0; i<row; i++)
		{
			memcpy(&DistMat[i*row],DATA[i],col*sizeof(double));
		}
	}
	else if (data_type==1)
	{
		HammingDist(row,col,DATA,DistMat);
	}
	else if (data_type==2)
	{
		EuclideanDist(row,col,DATA,DistMat);
	}

	for (int i=0; i<row; i++)
	{
		delete[] DATA[i];
	}
	delete[] DATA;

/* 	==========================================================================================
	Create all arrays needed for the runs
	For simplicity, all matrices are stored as one-dimensional arrays
   	========================================================================================== */

	int N = row;				// Size of distance matrix (# of subjects)
	double *Ensemble = new double[N*N];	// Ensemble matrix generated by the random walks
	double *Q = new double[N*N];		// work array
	double *Degree = new double[N];		// Work array
	bool *Membership = new bool[N*N];	// Membership table
	int *WorkI = new int[2*N+4*Neigen];
	double *WorkD = new double[2*N*Neigen+4*N+4*Neigen*(Neigen+6)];
	double *Centers = new double[Neigen*Neigen];
	int *ClusterID = new int[N];		// Cluster ID for each data point

	for(int i = 0; i < N*N; i++) 		// Initialize matrices to 0
	{
		Ensemble[i] = 0;
		Q[i] = 0;
		Membership[i] = 0;
	}
	int nc;					// Will contain number of clusters found
	
/* 	==========================================================================================
	Initialize clock
   	========================================================================================== */

	clock_t start_s = clock();

/* 	==========================================================================================
   	Generate all random number generators for all threads considered
   	========================================================================================== */

	MersenneTwister **mt = new MersenneTwister*[nthreads];
	for(int i = 0; i < nthreads; i++) 
	{
		mt[i]= new MersenneTwister();
		unsigned long seed = (unsigned long) time(NULL);
		mt[i]-> init_genrand(seed);
	}

/* 	==========================================================================================
   	Find all temperatures for DCG run
   	========================================================================================== */

	int iter = 50;
	int Nrepeat = 50;

	double *TempArray = new double[Neigen+1];
	int *NcArray = new int[Neigen+1];
	int NTemp;

	FindTemp(DistMat, N, M, P, Neigen, Ensemble, Q, Degree, WorkD, WorkI,
        Centers, Cutoff, iter, Nrepeat, Membership, ClusterID, Nmin, NcArray, TempArray, 
	&NTemp, nthreads, mt);

	cout << " " <<endl;
	cout << "NTemp = " << NTemp << endl;
	cout << " " <<endl;

/* 	==========================================================================================
   	Process temperatures:
		- if two temperatures have the same number of clusters, only keep the highest one
		- define actual temperatures as middle of intervals of temperatures found by FindTemp
   	========================================================================================== */

	double alpha = 100./(TempArray[NTemp-1]-TempArray[0]);
	double beta = - alpha*TempArray[0];

	int NTemp2=0;
	double *Temps = new double[3*NTemp];

	for (int i = 0; i < NTemp -1 ; i++) 
	{
		Temps[NTemp2]=TempArray[i];
		NTemp2++;
		Temps[NTemp2] = (TempArray[i]+TempArray[i+1])/2.0;
		NTemp2++;
	}
	Temps[NTemp2] = TempArray[NTemp-1];
	NTemp = NTemp2 + 1;

	cout << "new NTemp = " << NTemp << endl;
	for (int i = 0; i < NTemp; i++)
	{
		cout << "i = " << i << " Temp[i] = " << Temps[i] << endl;
	}
	cout << " " <<endl;

	delete [] TempArray;
	delete [] NcArray;

/* 	==========================================================================================
   	Process data at the chosen temperatures
   	========================================================================================== */

	double Temp;
	bool iprint=false;

	cout << "NTemp = " << NTemp << " N = " << N << endl;

	bool *Status = new bool[(N*(N-1))/2];
	double *DCG  = new double[(N*(N-1))/2];
	int *Nclust = new int[NTemp];

	for(int i = 0; i < (N*(N-1))/2; i++)
	{
		DCG[i] = alpha*Temps[NTemp-1]+beta;
		Status[i] = true;
	}

	cout << endl;
	cout << "Number of temperatures identified : " << NTemp << endl;
	cout << "Temperatures, and number of clusters associated with them: " << endl;
	logfile << "DCG procedure: " << endl;
	logfile << "============== " << endl;
	logfile << "Number of temperatures identified : " << NTemp << endl;
	logfile << "Temperatures, and number of clusters associated with them: " << endl;
        logfile << "            Idx    Temperature   Ncluster" << endl;

	double TempUp;
	Nclust[NTemp-1] = 1;
	for(int k = NTemp-2; k > -1; k--) 
	{
		Temp = Temps[k];
		TempUp = alpha*Temps[k+1]+beta;
		ProcessTemp(DistMat, N, Temp, M, P, Neigen, Ensemble, Q, Degree, WorkD, WorkI, Centers,
		Cutoff, iter, Nrepeat, Membership, &nc, ClusterID, nthreads, mt, iprint);
		Nclust[k] = nc;
		UpdateDCG(N, DCG, Status, Membership, TempUp);
		cout << "idx: " << k << " Temp: " << Temps[k] << " Number of clusters: " << Nclust[k] << endl;
		logfile << "TEMP" << setw(10) << k << " " << setw(10) << setprecision(3) << alpha*Temps[k]+beta << " " << setw(10) << Nclust[k] << endl;
	}
	logfile << " " << endl;
	FinalDCG(N, DCG, Status, alpha*Temps[0]+beta);

/* 	==========================================================================================
	Free up some space
   	========================================================================================== */

	delete [] Q;
	delete [] Degree;
	delete [] DistMat;
	delete [] Status;
	delete [] Membership;
	delete [] WorkI;
	delete [] WorkD;
	delete [] Centers;
	delete [] ClusterID;

/* 	==========================================================================================
	Write resulting DCG matrix into a ASCII file
   	========================================================================================== */

	logfile << "Final DCG matrix:" << endl;
	logfile << "=================" << endl;
	logfile << "ASCII file containing DCG matrix: " << fileDCG << endl;
	logfile << " " << endl;

	CopyDCG(N, DCG, Ensemble);
	WriteDCG(fileDCG, N, Ensemble);

/* 	==========================================================================================
	Stop clock to get computing time
   	========================================================================================== */

	clock_t stop_s = clock();

	cout << " " << endl;
	cout << "Total running time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " seconds" << endl;
	cout << " " << endl;

	logfile << "CPU info:" << endl;
	logfile << "=========" << endl;
	logfile << "Number of threads used: " << nthreads << endl;
	logfile << "Total running time    : " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " seconds" << endl;
	logfile << " " << endl;

	logfile.close();

/* 	==========================================================================================
	Terminate program
   	========================================================================================== */

	return 0;
}

/* ===============================================================================================
   Usage
 ================================================================================================= */

static void usage(char** argv)
{
    cout << "\n\n" <<endl;
    cout << "     " << "================================================================================================"<<endl;
    cout << "     " << "================================================================================================"<<endl;
    cout << "     " << "=                                                                                              ="<<endl;
    cout << "     " << "=                                         DCG                                                  ="<<endl;
    cout << "     " << "=                                                                                              ="<<endl;
    cout << "     " << "=    This program reads in a data matrix, either a distance matrix or a data matrix,           ="<<endl;
    cout << "     " << "=    and automatically clusters the data points by generating a random walk to build           ="<<endl;
    cout << "     " << "=    an Ensemble matrix, followed by Spectral Clustering                                       ="<<endl;
    cout << "     " << "=                                                                                              ="<<endl;
    cout << "     " << "=    Usage is:                                                                                 ="<<endl;
    cout << "     " << "=                                                                                              ="<<endl;
    cout << "     " << "=    DCG.exe -i FILE.in  -o FILE.res -n Ncluster -m M -p P -c Cutoff -a Temp  -t nthread       ="<<endl;
    cout << "     " << "=                                                                                              ="<<endl;
    cout << "     " << "=    where:                                                                                    ="<<endl;
    cout << "     " << "=                 -i FILE.in    --> Input file, contains the data matrix                       ="<<endl;
    cout << "     " << "=                 -o FILE.res   --> Output file for DCG distance matrix                        ="<<endl;
    cout << "     " << "=                 -n NCluster   --> Maximum number of clusters expected (default 20)           ="<<endl;
    cout << "     " << "=                 -p P          --> Number of random walks / point (default 10)                ="<<endl;
    cout << "     " << "=                 -m M          --> Number of steps / random walk (default 1000)               ="<<endl;
    cout << "     " << "=                 -c Cutoff     --> Cutoff value in eigenvalues selection (default 0.2)        ="<<endl;
    cout << "     " << "=                 -t nthread    --> Number of threads for multithreading (default 1)           ="<<endl;
    cout << "     " << "================================================================================================"<<endl;
    cout << "     " << "================================================================================================"<<endl;
    cout << "\n\n" <<endl;
}

/* ===============================================================================================
   Parse Argument from command line:

    -i FILE.in  --> input file for data matrix
    -o FILE.res --> Output file for DCG distance matrix
    -n Neigen   --> Number of eigenvalues considered  
    -m M        --> Number of steps in random walks (usually set to be 1000)
    -p P        --> Number of random walk / point (usually set to be 10) 
    -c Cutoff   --> Cutoff value in eigenvalues' selection(default 0.2) 
    -t nthread  --> Number of threads for multithreading
   =============================================================================================== */

bool parse_args(int argc, char **argv, string *fileIN,  string *fileOUT, int *Neigen, int *M,int *P, 
double *Cutoff, int *nthread)
{
//
// Make sure we have at least two parameters....
//
    string param;
    if (argc == 1) 
    {
        return false;
    } 
    else 
    {
        for (int i = 1; i < argc - 1; i = i + 2) 
        {
            param = argv[i];

            if (param == "-i") {
                *fileIN = argv[i+1];
            }

            else if (param == "-n") {
                *Neigen = atoi(argv[i + 1]);
            }
            else if (param == "-o") {
                *fileOUT = argv[i + 1];
            }
            else if (param == "-m") {
                *M = atoi(argv[i + 1]);
            }
            else if (param == "-p") {
                *P = atoi(argv[i + 1]);
            }
            else if (param == "-c") {
                *Cutoff = atof(argv[i + 1]);
            }
            else if (param == "-t") {
                *nthread = atoi(argv[i + 1]);
            }
        }
    }
    return true;
}

/* ===============================================================================================
   WriteDCG(fileOUT, N Ensemble)

What it does: writes DCG matrix into an ascii file
================================================================================================== */

void WriteDCG(string fileOUT, int N, double *Ensemble)
{
	ofstream output;
	output.open(fileOUT.c_str());

	int idx;
	for(int i = 0; i < N; i++)
	{
		for(int j=0; j < N; j++)
		{
			idx = N*i+j;
			output << Ensemble[idx] << " ";
		}
		output << endl;
	}

	output.close();
}
	

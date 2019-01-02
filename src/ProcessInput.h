/* ===============================================================================================

ProcessInput.h  (source code generated 2017-05-20)
Copyright (c) Jiahui Guan.
Copyright (c) Patrice Koehl.

================================================================================================== */

#ifndef _ProcessInput_
#define _ProcessInput_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <iomanip>
#include <cmath>

using namespace std;

/* ===============================================================================================
 Define prototypes for all functions defined in ProcessInput
 ================================================================================================= */

bool count_data(string fileIN, int *row, int *col);
void read_data(string fileIN, double **DATA);

void HammingDist(int row, int col, double **DATA, double *DistMat);
void EuclideanDist(int row, int col, double **DATA, double *DistMat);

#endif

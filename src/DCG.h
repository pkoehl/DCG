/* ===============================================================================================

DCG.h  (source code generated 2016-07-11)
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

#ifndef _DCG_
#define _DCG_


/* ===============================================================================================
 Set all includes for the program
 ================================================================================================= */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <vector>
#include <ctime> 


using namespace std;

/* ===============================================================================================
 Define prototypes for all functions defined in main, and 
 ================================================================================================= */

static void usage(char** argv);

bool parse_args(int argc, char **argv, string *fileIN,  string *fileOUT, int *Neigen,
int *M,int *P, double *Cutoff, int *nthreads);

void WriteDCG(string fileOUT, int N, double *Ensemble);

#endif

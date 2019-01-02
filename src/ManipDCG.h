/* ===============================================================================================

ManipDCG.h  (source code generated 2017-05-20)
Copyright (c) Jiahui Guan.
Copyright (c) Patrice Koehl.

================================================================================================== */

#ifndef _ManipDCG_
#define _ManipDCG_

/* ===============================================================================================
 Define prototypes for all functions defined in ManipDCG
 ================================================================================================= */

void UpdateDCG(int N, double *DCG, bool *Status, bool *Member, double Temp);
void FinalDCG(int N, double *DCG, bool *Status, double Temp);
void CopyDCG(int N, double *DCG, double *Mat);

#endif

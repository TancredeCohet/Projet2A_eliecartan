/*
 *  mexfiles ClearSparseMatrix
 *
 *  April 2008
 *
 *  Contributors:   Jean-Baptiste TAVERNIER - INRIA - CORIDA Project Team
 *                  Nicolas SAUGNIER - INRIA - ALICE Project Team
 *
 *  MATLAB Usage : ClearSparseMatrix()
 *
 *  See ClearSparseMatrix.m for more information.
 *
 */

#include "mex.h"
#include "libBL.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    (void) nlhs; /* Unused parameters */
    (void) plhs;
    (void) nrhs;
    (void) prhs;
    
    ClearSparseMatrix();
}

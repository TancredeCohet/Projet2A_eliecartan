/*
 *  mexfiles DeleteSparseMatrix
 *
 *  April 2008
 *
 *  Contributors:   Jean-Baptiste TAVERNIER - INRIA - CORIDA Project Team
 *                  Nicolas SAUGNIER - INRIA - ALICE Project Team
 *
 *  MATLAB Usage : DeleteSparseMatrix(MatrixID)
 *
 *  See DeleteSparseMatrix.m for more information.
 *
 */

#include "mex.h"
#include "libBL.h"
#include "common_tests.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    (void) nlhs; /* Unused parameters */
    (void) plhs;
    
    int MatrixID;
    
    // nbr of arguments
    TestNbrArgs(nrhs, 1, 1);
    
    // Type check for matrix ID
    TestValidMatrixID(prhs[0], "The first argument");
   
    
    MatrixID=(int)mxGetScalar(prhs[0]);
    
    DeleteSparseMatrix(MatrixID);
    
    
}

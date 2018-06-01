/*
 *  mexfiles NnzSparseMatrix
 *
 *  April 2008
 *
 *  Contributors:   Jean-Baptiste TAVERNIER - INRIA - CORIDA Project Team
 *                  Nicolas SAUGNIER - INRIA - ALICE Project Team
 *
 *  MATLAB Usage : nnz = NnzSparseMatrix(MatrixID)
 *      nnz numeric scalar.
 *
 *  See NnzSparseMatrix.m for more information.
 *
 */

#include "mex.h"
#include "libBL.h"
#include "common_tests.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    int MatrixID, nnz;
    
    // nbr of arguments
    TestNbrArgs(nrhs, 1, 1);
    
    // Type check for matrix ID
    TestValidMatrixID(prhs[0], "The first argument");
    
    
    MatrixID=(int)mxGetScalar(prhs[0]);
    
    nnz=NnzSparseMatrix(MatrixID);
    
    if (nlhs==1){
        plhs[0]=mxCreateDoubleScalar(nnz);
    }
    
}

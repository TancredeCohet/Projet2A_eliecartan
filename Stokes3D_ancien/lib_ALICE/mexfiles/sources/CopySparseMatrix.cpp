/*
 *  mexfile CopySparseMatrix
 *
 *  INRIA Fevrier 2012
 *
 *
 *  MATLAB Usage : A_c=CopySparseMatrix(A)
 *      A is the matrix ID (numeric scalar).
 *      A_c is the matrix ID of the copy
 *
 */

#include "mex.h"
#include "libBL.h"
#include "common_tests.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    int MatrixID, newID;
    double i,j;
    
    // nbr of arguments
    TestNbrArgs(nrhs, 1, 1);
    
    // Type check for matrix ID
    TestValidMatrixID(prhs[0], "The first argument");
    
    MatrixID=(int)mxGetScalar(prhs[0]);
    
    newID=CopySparseMatrix(MatrixID);

    if (nlhs==1){
        plhs[0]=mxCreateDoubleScalar(newID);
    }
}


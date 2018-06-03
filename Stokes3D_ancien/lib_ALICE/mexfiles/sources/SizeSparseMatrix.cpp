/*
 *  mexfile SizeSparseMatrix
 *
 *  INRIA Janvier 2012
 *
 *
 *  MATLAB Usage : [m,n] = SizeSparseMatrix(MatrixID)
 *      I,J numeric scalars.
 *
 */

#include "mex.h"
#include "libBL.h"
#include "common_tests.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    int MatrixID;
    double val;
    
    // nbr of arguments
    TestNbrArgs(nrhs, 1, 1);
    
    // Type check for matrix ID
    TestValidMatrixID(prhs[0], "The first argument");
    
    MatrixID=(int)mxGetScalar(prhs[0]);
    
    gx_size_T m = SparseMatrixM(MatrixID);
    gx_size_T n = SparseMatrixN(MatrixID);
    if (nlhs == 2){
        plhs[0] = mxCreateDoubleScalar(static_cast<double>(m));
        plhs[1] = mxCreateDoubleScalar(static_cast<double>(n));
    }
}

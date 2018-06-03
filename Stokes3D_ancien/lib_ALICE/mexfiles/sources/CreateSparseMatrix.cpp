/*
 *  mexfiles CreateSparseMatrix
 *
 *  April 2008
 *
 *  Contributors:   Jean-Baptiste TAVERNIER - INRIA - CORIDA Project Team
 *                  Nicolas SAUGNIER - INRIA - ALICE Project Team
 *
 *  MATLAB Usage : A=CreateSparseMatrix(I,J)
 *      A is the matrix ID (numeric scalar).
 *      I,J numeric scalars.
 *
 *  See CreateSparseMatrix.m for more information.
 *
 */

#include "mex.h"
#include "libBL.h"
#include "common_tests.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    int MatrixID;
    double i,j;
    
    // nbr of arguments
    TestNbrArgs(nrhs, 2, 2);
    
    // Type check for I and J
    TestNumber(prhs[0], "The first argument");
    TestNumber(prhs[1], "The second argument");
    
    
    i=mxGetScalar(prhs[0]);
    j=mxGetScalar(prhs[1]);
    TestVectStrictPositive(&i, 1, "the first argument");
    TestVectStrictPositive(&j, 1, "the second argument");
    
    
    MatrixID=NewSparseMatrix((int)i,(int)j);
    
    if (nlhs==1){
        plhs[0]=mxCreateDoubleScalar(MatrixID);
    }
}


/*
 *  mexfiles coefficient
 *
 *  April 2008
 *
 *  Contributors:   Jean-Baptiste TAVERNIER - INRIA - CORIDA Project Team
 *                  Nicolas SAUGNIER - INRIA - ALICE Project Team
 *
 *  MATLAB Usage : a = coefficient(MatrixID,I,J)
 *      I,J numeric scalars.
 *      a numeric scalar
 *
 *  See coefficient.m for more information.
 *
 */

#include "mex.h"
#include "libBL.h"
#include "common_tests.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    int MatrixID;
    double i,j;
    double val;
    
    // nbr of arguments
    TestNbrArgs(nrhs, 3, 3);
    
    // Type check for matrix ID
    TestValidMatrixID(prhs[0], "The first argument");
    
    // Type check for I, J
    TestNumber(prhs[1], "The second argument");
    TestNumber(prhs[2], "The third argument");
    
    
    MatrixID=(int)mxGetScalar(prhs[0]);
    
    i=mxGetScalar(prhs[1])-1;
    j=mxGetScalar(prhs[2])-1;
    TestVectRange( &i, 1, 0, SparseMatrixM(MatrixID)-1, "the second argument");
    TestVectRange( &j, 1, 0, SparseMatrixN(MatrixID)-1, "the third argument");
    
    val=coefficient(MatrixID,(int)i,(int)j);
    
    if (nlhs == 1){
        plhs[0]=mxCreateDoubleScalar(val);
    }

}

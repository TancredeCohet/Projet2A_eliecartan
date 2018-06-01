/*
 *  mexfiles SparseMatrixToMatlab
 *
 *  April 2008
 *
 *  Contributors:   Jean-Baptiste TAVERNIER - INRIA - CORIDA Project Team
 *                  Nicolas SAUGNIER - INRIA - ALICE Project Team
 *
 *  MATLAB Usage : A=SparseMatrixToMatlab(MatrixID)
 *      A numeric sparse matrix.
 *
 *  See SparseMatrixToMatlab.m for more information.
 *
 */

#include "mex.h"
#include "libBL.h"
#include "common_tests.h"


#if defined(NAN_EQUALS_ZERO)
#define IsNonZero(d) ((d)!=0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d)!=0.0)
#endif


void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    int MatrixID;
    
    double *sr;
    
    mwIndex *irs, *jcs;
    
    // nbr of arguments
    TestNbrArgs(nrhs, 1, 1);
    
    // Type check for matrix ID
    TestValidMatrixID(prhs[0], "The first argument");
    
    
    MatrixID=(int)mxGetScalar(prhs[0]);
    
    if (nlhs==1){
        plhs[0] = mxCreateSparse(
        SparseMatrixM(MatrixID),
        SparseMatrixN(MatrixID),
        NnzSparseMatrix(MatrixID),
        mxREAL
        );
        
        sr  = mxGetPr(plhs[0]);
        
        irs = mxGetIr(plhs[0]);
        
        jcs = mxGetJc(plhs[0]);
        
        SparseMatrixToMatlab(MatrixID, irs, jcs, sr);
    }
    
    
    
}








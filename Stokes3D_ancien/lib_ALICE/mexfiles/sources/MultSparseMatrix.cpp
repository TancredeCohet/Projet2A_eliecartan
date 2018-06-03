/*
 *  mexfiles MultSparseMatrix
 *
 *  April 2008
 *
 *  Contributors:   Nicolas SAUGNIER - INRIA - ALICE Project Team
 *
 *  MATLAB Usage : C_ID=MultSparseMatrix( A_ID, B_ID, A_Transpose, B_Transpose )
 *      A_ID, B_ID and C_ID sparse matrix indices.
 *      A_Transpose and B_Transpose facultative scalars or logicals. 
 *      
 *
 *  See MultSparseMatrix.m for more information.
 *
 */

#include "mex.h"
#include "libBL.h"
#include "common_tests.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    int A_ID, B_ID, C_ID;
    bool A_trans = false;
    bool B_trans = false;
    
    int C_m, C_n, nA, mB;
    
    // nbr of arguments
    TestNbrArgs(nrhs, 2, 4);
    
    // Type check for matrix IDs
    TestValidMatrixID(prhs[0], "The first argument");
    TestValidMatrixID(prhs[1], "The second argument");
    
    // Type check and treatement for Transpose
    if (nrhs >= 3){
        if ( !mxIsLogicalScalar(prhs[2]) ){
            TestNumericRealNotEmpty(prhs[2], "The third argument");
            A_trans = ( mxGetPr(prhs[2])[0] != 0.0 );
        } else {
            A_trans = mxIsLogicalScalarTrue(prhs[2]);
        }
        
        if (nrhs == 4){
            if ( !mxIsLogicalScalar(prhs[3]) ){
                TestNumericRealNotEmpty(prhs[3], "The fourth argument");
                B_trans = ( mxGetPr(prhs[3])[0] != 0.0 );
            } else {
                B_trans = mxIsLogicalScalarTrue(prhs[3]);
            }
        }
    }
    
    
    A_ID=(int)mxGetScalar(prhs[0]);
    B_ID=(int)mxGetScalar(prhs[1]);
    
    if(!A_trans){
        nA   = SparseMatrixN(A_ID);
        C_m = SparseMatrixM(A_ID); 
    } else {
        nA   = SparseMatrixM(A_ID);
        C_m = SparseMatrixN(A_ID); 
    }
    
    if(!B_trans){
        mB   = SparseMatrixM(B_ID);
        C_n = SparseMatrixN(B_ID); 
    } else {
        mB   = SparseMatrixN(B_ID);
        C_n = SparseMatrixM(B_ID); 
    }
    
    if (nA != mB){
        mexErrMsgIdAndTxt("OpenNL:MultSparseMatrix:MatrixSizeMismatch",
            "Matrix dimensions mismatch. Taking transpose into account, first matrix' N = %d, second matrix' M = %d.",
            nA,mB);
    }
    
    // we should test the possibility to create a new matrix...
    C_ID=NewSparseMatrix(C_m,C_n);
    
    
    if (!A_trans){
        if (!B_trans){
            MultABSparseMatrix(A_ID,B_ID,C_ID);
        } else {
            MultABtSparseMatrix(A_ID,B_ID,C_ID);
        }
    } else {
        if (!B_trans){
            MultAtBSparseMatrix(A_ID,B_ID,C_ID);
        } else {
            MultAtBtSparseMatrix(A_ID,B_ID,C_ID);
        }
    }
    
    
    if (nlhs==1){
        plhs[0]=mxCreateDoubleScalar(C_ID);
    }
    
}

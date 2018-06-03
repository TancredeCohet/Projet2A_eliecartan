/*
 *  mexfiles RemoveRowsColsSparseMatrix
 *
 *  April 2008
 *
 *  Contributors:   Nicolas SAUGNIER - INRIA - ALICE Project Team
 *
 *  MATLAB Usage :  RemoveRowsColsSparseMatrix(MatrixID, I)
 *                  B2=RemoveRowsColsSparseMatrix(MatrixID, I, B)
 *          I, B, B2 numeric vectors.
 *
 *  See RemoveRowsColsSparseMatrix.m for more information.
 *
 */
#include "stdafx.h"
#include "mex.h"
#include "libBL.h"
#include "common_tests.h"
#include <algorithm>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    int MatrixID;
    double* indexVect; 
    
    gx_size_T nbIndex;
    gx_index_T* indexInt;

    
    // nbr of arguments
    TestNbrArgs(nrhs, 2, 3);
    
    // Type check for matrix ID
    TestValidMatrixID(prhs[0], "The first argument");
    
    // Type check for I
    TestVector(prhs[1], "The second argument");
    if (mxIsSparse(prhs[1])){
         mexErrMsgIdAndTxt("OpenNL:AddMatElem:SparseUnsupported",
         "The second argument must be a full matrix.");
    }
    
    // TypeCheck for B
    if (nrhs == 3){
    	TestVector(prhs[2], "The third argument");
        if (mxIsSparse(prhs[2]) && mxGetN(prhs[2])!=1 ){
            mexErrMsgIdAndTxt("OpenNL:AddMatElem:SparseUnsupported",
            "The third argument must be either a full vector or a sparse "
            "column vector.");
        }
    }
    
    
    MatrixID = (int)mxGetScalar(prhs[0]);
    
    nbIndex = mxGetNumberOfElements(prhs[1]);
    
    indexInt = new gx_index_T[nbIndex];
    
    indexVect = mxGetPr(prhs[1]);
    TestVectRange(indexVect,nbIndex,1,Geex::gx_min(SparseMatrixM(MatrixID),SparseMatrixN(MatrixID)),"the second argument");
    
    
    for(unsigned int i = 0; i < nbIndex; i++){
        indexInt[i] = (gx_index_T) indexVect[i] - 1 ; // first matlab indice is 1...
    }
    
    std::sort(indexInt, indexInt+nbIndex);
    RemoveColsSparseMatrix(MatrixID, nbIndex, indexInt);
    RemoveRowsSparseMatrix(MatrixID, nbIndex, indexInt);
    for(unsigned int i = 0; i< nbIndex; i++){
        AddSparseMatrix(MatrixID, indexInt[i],indexInt[i],1);
    }
    
    
    
    if ( nrhs == 3 && nlhs == 1){
        if ( mxIsSparse(prhs[2]) ){
            mwSize m = mxGetM(prhs[2]);
            mwSize n = mxGetN(prhs[2]);
            mwSize nzmax = mxGetNzmax(prhs[2]);
            
            plhs[0] = mxCreateSparse(m, n, nzmax, mxREAL);
            
            double* inputValue = mxGetPr(prhs[2]);
            mwIndex*   inputIr = mxGetIr(prhs[2]);
            mwIndex*   inputJc = mxGetJc(prhs[2]);
            
            /* Jc has 2 elements : Jc[0] must be 0, and
             * Jc[1] = nnz; and nnz <= nzmax..
             */
            mwSize nnz = inputJc[1];
            
            double* outputValue = mxGetPr(plhs[0]);
            mwIndex*   outputIr = mxGetIr(plhs[0]);
            mwIndex*   outputJc = mxGetJc(plhs[0]);
        
            // prhs[2] may not be sorted... :(
            
            outputJc[0] = 0;
                        
            bool copyThisValue;
            size_t i_out = 0;
            for( size_t i_in = 0; i_in<nnz; i_in++ ){
                copyThisValue = true;
                for ( size_t k = 0; k<nbIndex; k++){
                   if ( indexInt[k] >= inputIr[i_in]  ){
                       if (inputIr[i_in] == indexInt[k]){
                           copyThisValue = false;
                       }
                       break;
                   }
                }
                if ( copyThisValue ){
                    outputValue[i_out] = inputValue[i_in];
                    outputIr[i_out] = inputIr[i_in];
                    i_out++;
                }
            }
            
            outputJc[1] = i_out ; // i_out == nnz for the output
            
            
        } else {
            plhs[0] = mxCreateDoubleMatrix( mxGetM(prhs[2]), mxGetN(prhs[2]), mxREAL);

            double* inputValue = mxGetPr(prhs[2]);
            double* outputValue = mxGetPr(plhs[0]);

            size_t lastSeenIndex = 0;
            size_t nbInputValue = mxGetNumberOfElements(prhs[2]);

            for (size_t j = 0; j < nbInputValue ; j++){
                if ( indexInt[lastSeenIndex] != j ){
                    outputValue[j] = inputValue[j];
                } else {
                    outputValue[j] = 0;
                    lastSeenIndex++;
                }
            }
        }
    }
}



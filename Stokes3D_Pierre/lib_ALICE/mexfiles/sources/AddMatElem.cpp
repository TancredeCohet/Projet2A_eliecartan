/*
 *  mexfiles AddMatElem
 *
 *  April 2008
 *
 *  Contributors:   Jean-Baptiste TAVERNIER - INRIA - CORIDA Project Team
 *                  Nicolas SAUGNIER - INRIA - ALICE Project Team
 *
 *  MATLAB Usage : AddMatElem(MatrixID,I,J,VAL)
 *      I,J numeric full vectors, VAL numeric matrix.
 *
 *  See AddMatElem.m for more information.
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
    int long_i, long_j;
    int long_valM, long_valN;
    int k=0;
    
    double *vec_i, *vec_j;
    
    double *vec_val;
    
    // nbr of arguments
    TestNbrArgs(nrhs, 4, 4);
    
    // Type check for matrix ID
    TestValidMatrixID(prhs[0], "The first argument");
    
    // Type check for I and J
    TestVector(prhs[1], "The second argument");
    TestVector(prhs[2], "The third argument");
    if (mxIsSparse(prhs[1])){
        mexErrMsgIdAndTxt("OpenNL:AddMatElem:SparseUnsupported",
         "The second argument must be a full vector.");
    }
    if (mxIsSparse(prhs[2])){
        mexErrMsgIdAndTxt("OpenNL:AddMatElem:SparseUnsupported",
         "The third argument must be a full vector.");
    }
    
    // Type check for Val
    TestNumericRealNotEmpty(prhs[3], "The fourth argument");
    
    
    
    MatrixID=(int)mxGetScalar(prhs[0]);
   
    vec_i=mxGetPr(prhs[1]);
    long_i=mxGetNumberOfElements(prhs[1]);
    
    vec_j=mxGetPr(prhs[2]);
    long_j=mxGetNumberOfElements(prhs[2]);
    
    vec_val=mxGetPr(prhs[3]);
    long_valM=mxGetM(prhs[3]);
    long_valN=mxGetN(prhs[3]);
    
    // Consistance check
    if ((long_i!=long_valM) || (long_j!=long_valN)){
       mexErrMsgIdAndTxt("OpenNL:AddMatElem:InconsistantInput","Dim(VAL) = %d * %d\n"
       "Dim(I) is %d, should be equal to first Dim of VAL.\n"
       "Dim(J) is %d, should be equal to second Dim of VAL.",
       long_valM,long_valN,long_i,long_j);
    }
    TestVectRange(vec_i,long_i,1,SparseMatrixM(MatrixID),"the second argument");
    TestVectRange(vec_j,long_j,1,SparseMatrixN(MatrixID),"the third argument");
    
    
    if ( mxIsSparse(prhs[3]) ){
        
        mwIndex* ir = mxGetIr(prhs[3]);
        mwIndex* jc = mxGetJc(prhs[3]);
        
        int i;
        int j = 0;
        // jc[long_valN] == nnz...
        for (unsigned int l=0; l< jc[long_valN] ; l++){
            i = ir[l];
            
            while( (size_t)l >= jc[j+1] ){
                j++;
            }
            //mexPrintf("J'ajoute %f à (%d,%d)\n",vec_val[l], (int)vec_i[i]-1, (int)vec_j[j]-1);
            if (vec_val[k] != 0.0){
                AddSparseMatrix(MatrixID,(int)vec_i[i]-1,(int)vec_j[j]-1,vec_val[l]);
            }
        }
       
    } else {
        // matrix is full
        for (int j=0; j<long_j; j++){
            for (int i=0; i<long_i; i++){
                if (vec_val[k] != 0.0){
                    AddSparseMatrix(MatrixID,(int)vec_i[i]-1,(int)vec_j[j]-1,vec_val[k]);
                }
                k++;
            }
        }
    }
    
    
}


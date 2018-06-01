/*
 *  mexfiles AddSparseMatrix
 *
 *  April 2008
 *
 *  Contributors:   Jean-Baptiste TAVERNIER - INRIA - CORIDA Project Team
 *                  Nicolas SAUGNIER - INRIA - ALICE Project Team
 *
 *  MATLAB Usage : AddSparseMatrix(MatrixID,I,J,VAL)
 *      I,J,VAL numeric scalars.
 *
 *  See AddSparseMatrix.m for more information.
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
    TestNbrArgs(nrhs, 4, 4);
    
    // Type check for matrix ID
    TestValidMatrixID(prhs[0], "The first argument");
    
    // Type check for I, J, VAL
    TestVector(prhs[1], "The second argument");
    TestVector(prhs[2], "The third argument");
    TestVector(prhs[3], "The fourth argument");
   
    if (mxGetNumberOfElements(prhs[1]) != mxGetNumberOfElements(prhs[2]) ||
        (mxGetNumberOfElements(prhs[2]) != mxGetNumberOfElements(prhs[3]) ) ||
        (mxGetNumberOfElements(prhs[1]) != mxGetNumberOfElements(prhs[3]) ) ) {
            mexErrMsgIdAndTxt("OpenNL:AddSparseMatrix:SparseUnsupported",
            "the three last arguments must have the same length");
    }     
    if (mxIsSparse(prhs[1]) || mxIsSparse(prhs[2]) || mxIsSparse(prhs[3]) ) {
            mexErrMsgIdAndTxt("OpenNL:AddSparseMatrix:SparseUnsupported",
            "arguments must be full vectors");
    }
    MatrixID=(int)mxGetScalar(prhs[0]);
    
    // double pointers on arguments
    double * i_p = mxGetPr(prhs[1]);
    double * j_p = mxGetPr(prhs[2]);
    gx_size_T taille = mxGetNumberOfElements(prhs[1]) ;
  
    // local copies for shifting
    gx_index_T * i = new gx_index_T[taille]; 
    gx_index_T * j = new gx_index_T[taille]; 


    // c tables begin to zero
    for (gx_size_T k = 0; k  < taille ; ++k) {
        i[k] = (gx_index_T) i_p[k] -1;
        j[k] = (gx_index_T) j_p[k] -1;
    }

    TestVectRange( i, taille, 0, SparseMatrixM(MatrixID)-1, "the second argument");
    TestVectRange( j, taille, 0, SparseMatrixN(MatrixID)-1, "the third argument");
    
    double * Valeur = mxGetPr(prhs[3]);
    
    for (gx_size_T k = 0; k  < taille ; ++k) {
       AddSparseMatrix(MatrixID,i[k],j[k],Valeur[k]);
    }
    delete [] i;
    delete [] j;
}


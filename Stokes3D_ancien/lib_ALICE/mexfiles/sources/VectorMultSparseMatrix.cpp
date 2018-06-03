/*
 *  mexfiles VectorMultSparseMatrix
 *
 *  April 2008
 *
 *  Contributors:   Nicolas SAUGNIER - INRIA - ALICE Project Team
 *
 *  MATLAB Usage : Y=VectorMultSparseMatrix(MatrixID,X,Transpose)
 *      X,Y numeric vectors. If full : row or column. If sparse : column.
 *      Transpose facultative scalar or logical. 
 *      
 *
 *  See VectorMultSparseMatrix.m for more information.
 *
 */

#include "mex.h"
#include "libBL.h"
#include "common_tests.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    (void) nlhs; // unused parameter
    
    int MatrixID;
    double *x,*y;
    bool trans = false;
    
    
    // nbr of arguments
    TestNbrArgs(nrhs, 2, 3);
    
    // Type check for matrix ID
    TestValidMatrixID(prhs[0], "The first argument");
    
    // Type check for X
    TestVector(prhs[1], "The second argument");
    
    // Type check and treatement for Transpose
    if (nrhs == 3){
        if ( !mxIsLogicalScalar(prhs[2]) ){
            TestNumericRealNotEmpty(prhs[2], "The third argument");
            trans = ( mxGetPr(prhs[2])[0] != 0.0 );
        } else {
            trans = mxIsLogicalScalarTrue(prhs[2]);
        }
    }
    
    
    MatrixID=(int)mxGetScalar(prhs[0]);
    x=mxGetPr(prhs[1]);
    
    if (mxIsSparse(prhs[1])){
        mexErrMsgIdAndTxt("OpenNL:AddMatElem:NYI",
         "Mult matrix * sparse vector : Not yet implemented.");
    } else {
        //if (nlhs==1){
            plhs[0] =  mxCreateDoubleMatrix(
                mxGetM(prhs[1])==1?1:SparseMatrixM(MatrixID),
                mxGetN(prhs[1])==1?1:SparseMatrixN(MatrixID),
                mxREAL);
            y = mxGetPr(plhs[0]);
        //}

        if (trans){
            if ( mxGetNumberOfElements(prhs[1]) != SparseMatrixM(MatrixID) ){
                mexErrMsgIdAndTxt("OpenNL:VectorMultSparseMatrix:InconsistantInput",
                "Vector length do not match first matrix dimension. Len = %d, should be %d.",
                mxGetNumberOfElements(prhs[1]),SparseMatrixM(MatrixID));
            }
            VectorMultTransposeSparseMatrix(MatrixID, x, y);
        } else {
            if ( mxGetNumberOfElements(prhs[1]) != SparseMatrixN(MatrixID) ){
                mexErrMsgIdAndTxt("OpenNL:VectorMultSparseMatrix:InconsistantInput",
                "Vector length do not match second matrix dimension. Len = %d, should be %d.",
                mxGetNumberOfElements(prhs[1]),SparseMatrixN(MatrixID));
            }
            VectorMultSparseMatrix(MatrixID, x, y);
        }

    }
    
}



/*
 *  mexfiles SolveSparseMatrix
 *
 *  April 2008
 *
 *  Contributors:   Nicolas SAUGNIER - INRIA - ALICE Project Team
 *
 *  MATLAB Usage : X=SolveSparseMatrix(MatrixID, SolverName ,B)
 *      Let A be the matrix which ID is MatrixID, SolveSparseMatrix solves
 *      the system AX=B with the solver named SolverName.
 *      SolverName char vector. See SolveSparseMatrix.m for valid values.
 *      B numeric full vector or numeric sparse column vector.
 *
 *
 *  NOT YET IMPLEMENTED !!!
 *
 *  See SolveSparseMatrix.m for more information.
 *
 */

#include "mex.h"
#include "libBL.h"
#include "common_tests.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    (void) plhs;
    
    int MatrixID;
    
    // nbr of arguments
    TestNbrArgs(nrhs, 3, 3);
    
    // Type check for matrix ID
    TestValidMatrixID(prhs[0], "The first argument");
   
    // TODO test Second argument : vector of char
    // TODO est second arg value (one of the available solvers)
    
    // TODO test third argument : vector of numeric
    // TODO test size of the third argument compared to the matrix
    // TODO test full or sparse
    
    MatrixID = (int)mxGetScalar(prhs[0]);
   
    
    // TODO treat B according to his format (full/sparse)
    // TODO create the result vector to be filled by the solver
    
    // TODO call the solver... 
    
    if (nlhs==1){
        // TODO affect result to plhs[0]
    }
}


/*
 *  mexfiles SetBoundaryConditions
 *
 *  April 2008
 *
 *  Contributors:   Nicolas SAUGNIER - INRIA - ALICE Project Team
 *                  Marc Fuentes - INRIA - SED
 *  MATLAB Usage :  SetBoundaryConditions(MatrixID, I)
 *                  B2=SetBoundaryConditions(MatrixID, I, B)
 *                  B2=SetBoundaryConditions(MatrixID, I, B, V)
 *          I, B, B2, V numeric vectors.
 *
 *  See SetBoundaryConditions.m for more information.
 *
 */
// attention a mettre en premier pour eviter les redefinitiions de max
#include "coupled_sort.h"

#include "mex.h"
#include "libBL.h"
#include "timing.h"
#include "common_tests.h"
#include <algorithm>
#include <functional>
#include <vector>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    int MatrixID;
    double* indexVect; 
    
    gx_size_T nbIndex;
    gx_index_T* indexInt;

    // nbr of arguments
    TestNbrArgs(nrhs, 2, 4);
    
    // Type check for matrix ID
    TestValidMatrixID(prhs[0], "The first argument");
   
    MatrixID = (int)mxGetScalar(prhs[0]);
    
    gx_size_T M = SparseMatrixM(MatrixID);
    gx_size_T N = SparseMatrixN(MatrixID);

    // Type check for I
    TestVector(prhs[1], "The second argument");
    if (mxIsSparse(prhs[1])){
         mexErrMsgIdAndTxt("OpenNL:AddMatElem:SparseUnsupported",
         "The second argument must be a full matrix.");
    }
    
    // TypeCheck for B
    if (nrhs >= 3){
    	TestVector(prhs[2], "The third argument");
        if (mxGetNumberOfElements(prhs[2])!=M) {
            mexErrMsgIdAndTxt("OpenNL:AddMatElem:SparseUnsupported",
            "The third argument must have the same number of elements"
            " that the number of rows of the matrix.");
	}
	if (mxIsSparse(prhs[2]) && mxGetN(prhs[2])!=1   ){
            mexErrMsgIdAndTxt("OpenNL:AddMatElem:SparseUnsupported",
            "The third argument must be either a full vector or a sparse "
            "column vector.");
        }
        // TypeCheck for V
        if (nrhs == 4){
    	   TestVector(prhs[3], "The fourth argument");
           if (M!=N) {
            mexErrMsgTxt("OpenNL:SetBoundaryConditions FIXME : non homogenous second member not supported yet in the rectangular case M!=N");
	   }
	   if (mxIsSparse(prhs[3]) || mxGetNumberOfElements(prhs[3])!=mxGetNumberOfElements(prhs[1]) ){
	    mexErrMsgTxt("OpenNL:SetBoundaryConditions: The forth argument must be full vector and have same length that second argument");
           }
       
        }
    }
    
    
    nbIndex = mxGetNumberOfElements(prhs[1]);
    
    indexInt = new gx_index_T[nbIndex];
    
    indexVect = mxGetPr(prhs[1]);
    TestVectRange(indexVect,nbIndex,1,Geex::gx_min(M,N),"indices of the second argument are not valid");
    
    
    for(unsigned int i = 0; i < nbIndex; i++){
        indexInt[i] = (gx_index_T) indexVect[i] - 1 ; // first matlab indice is 1...
    }
   
    if (nrhs <= 2) {
        std::sort(indexInt, indexInt+nbIndex);
        RemoveColsSparseMatrix(MatrixID, nbIndex, indexInt);
        RemoveRowsSparseMatrix(MatrixID, nbIndex, indexInt);
        for(unsigned int i = 0; i< nbIndex; i++){
            AddSparseMatrix(MatrixID, indexInt[i],indexInt[i],1);
        }
    }
    // third arguments
    else if ( nrhs == 3  && nlhs == 1){
        double *outputValue;
	
	mwSize m_f = mxGetM(prhs[2]);
        mwSize n_f = mxGetN(prhs[2]);
        // rhs size 
	size_t nbInputValue = mxGetNumberOfElements(prhs[2]);
	
	if ( mxIsSparse(prhs[2])){
            mwSize nzmax = mxGetNzmax(prhs[2]);
            plhs[0] = mxCreateSparse(m_f, n_f, nzmax, mxREAL);
        }
	
	else {
            plhs[0] = mxCreateDoubleMatrix( m_f, n_f, mxREAL);
        }
	
	outputValue = mxGetPr(plhs[0]);
	
	std::sort(indexInt, indexInt+nbIndex);
   
        // suppress cols and rows of with indices in indexInt
        RemoveColsSparseMatrix(MatrixID, nbIndex, indexInt);
        RemoveRowsSparseMatrix(MatrixID, nbIndex, indexInt);
        
	// put ones on the diagonal corresponding to the suppress indices
	for(unsigned int i = 0; i< nbIndex; i++){
            AddSparseMatrix(MatrixID, indexInt[i],indexInt[i],1);
        } 
        	
	
	// third argument is sparse
	if ( mxIsSparse(prhs[2]) ){
           
           double* inputValue = mxGetPr(prhs[2]);
           mwIndex*   inputIr = mxGetIr(prhs[2]);
           mwIndex*   inputJc = mxGetJc(prhs[2]);
	   mwSize         nnz = inputJc[1];
           
           // for the case where
	   // nrhs ==3 and the third argument is sparse
           size_t i_out ;
           mwIndex*   outputIr ; //Ai
           mwIndex*   outputJc ; //Ap

           /* Jc has 2 elements : Jc[0] must be 0, and
            * Jc[1] = nnz; and nnz <= nzmax..
            */
           outputIr = mxGetIr(plhs[0]); 
           outputJc = mxGetJc(plhs[0]); 
           outputJc[0] = 0;
           i_out = 0;
	   // loop on the all nnz entries of B
	   for( size_t i_in = 0; i_in<nnz; i_in++ ){
               // copy only the indices which are not in IndexInt
	       if ( !std::binary_search(indexInt, indexInt + nbIndex, inputIr[i_in])) {
	              outputValue[i_out] = inputValue[i_in];
                      outputIr[i_out] = inputIr[i_in];
                      i_out++;
	           }
           }
	   outputJc[1] = i_out ; // i_out == nnz for the output
        } 
	// third  argument is full
	else {

           double* inputValue = mxGetPr(prhs[2]);

	   // copy nnz values and put to zero
	   // the values from indices(indexInt)
	   size_t lastSeenIndex = 0;
           for (size_t j = 0; j < nbInputValue ; j++){
              if ( indexInt[lastSeenIndex] != j ){
                 outputValue[j] = inputValue[j];
              } 
	      else {
                 outputValue[j] = 0;
                 lastSeenIndex++;
              }
           }
	   
	}
        
    }	
    // non homogenous case 
    // subtract the vector res = A(:,ind)*g
    else if (nrhs == 4 && nlhs == 1) {
        
	double duree;
	mwSize m_f = mxGetM(prhs[2]);
        mwSize n_f = mxGetN(prhs[2]);
        size_t nbInputValue = mxGetNumberOfElements(prhs[2]);
        
	// the left hand side
	plhs[0] = mxCreateDoubleMatrix( m_f, n_f, mxREAL);
	double *outputValue = mxGetPr(plhs[0]);
        std::fill(outputValue, outputValue + nbInputValue, 0);
        
        
        // sort indexInt and stock the permutation
	std::vector<size_t> perm;
	sort_with_perm(indexInt, indexInt+nbIndex, perm);
        double *gValue = mxGetPr(prhs[3]);
	size_t nbgValue =  mxGetNumberOfElements(prhs[3]); 
        
	// local copy to avoid to modify directly gValue
	double *gCopy = new double [nbgValue];
        std::copy(gValue, gValue + nbgValue, gCopy);
             
        // apply the permutation to gValue
        apply_perm( gCopy, gCopy +nbgValue, perm);
            
        // copy values of g at the right positions in gFull
        double * gFull = new double [N];
        std::fill( gFull, gFull+N, 0);
        for( size_t j = 0 ; j < nbgValue ; j++) { 
               gFull[indexInt[j]] = gCopy[j]; 
        }
 
        // res = A(:,index) * g(index)
        double * res = new double [M];
        std::fill( res,res+M, 0);
        VectorMultSparseMatrix( MatrixID, gFull , res) ;
    
        delete [] gFull;
        
        // suppress cols and rows of with indices in indexInt
        RemoveColsSparseMatrix(MatrixID, nbIndex, indexInt);
        RemoveRowsSparseMatrix(MatrixID, nbIndex, indexInt);
        
        // put ones on the diagonal corresponding to the suppress indices
        for(unsigned int i = 0; i< nbIndex; i++){
            AddSparseMatrix(MatrixID, indexInt[i],indexInt[i],1);
        } 
        // third argument is sparse
        if ( mxIsSparse(prhs[2]) ){
           double* inputValue = mxGetPr(prhs[2]);
           mwIndex*   inputIr = mxGetIr(prhs[2]);
           mwIndex*   inputJc = mxGetJc(prhs[2]);
           mwSize         nnz = inputJc[1];
           // straight copy nnz values from Input into Output
           // copy the right hand side
           for( size_t i_in = 0; i_in<nnz; i_in++ ) {
              outputValue[inputIr[i_in]] = inputValue[i_in];
           }
        } 
        else {
           // third argument is full
           double* inputValue = mxGetPr(prhs[2]);
           std::copy( inputValue, inputValue+nbInputValue, outputValue );
        }	   
        // put to g the values corresponding to indices
        // substract res to OutputValue 
        std::transform( outputValue, outputValue + nbInputValue, res,
                        outputValue, std::minus<double>() );
        
	// and now put the g's values at the indices coordinates!
	for( size_t j = 0 ; j < nbIndex ;  ++j)
                   outputValue[indexInt[j]] = gCopy[j];
        
	delete [] gCopy;
        delete [] res;
    }
    else {
       mexErrMsgTxt("Wrong number of arguments or lhs void");
    }
    delete [] indexInt ;
}



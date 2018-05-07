#include "mex.h"

/* interface MATLAB pour COLLOCATION :

    [Ap,Ai,Ax] = CCS(M) 

    INPUT : 
     - M sparse matrix to "explore" 
    OUTPUT : 
     - Ap column index vector
     - Ai row index vector 
     - Ax data vector
  */


   
void mult(mwIndex * row , mwIndex * col , double * val , double * vect ,double * result,int N) {
	int i;
	int j; 
	int k,indval;
        double sum;
#pragma omp parallel for schedule(dynamic) private (indval,sum,j,k) 
	for (i = 0 ; i< N;i++) {
		indval = row[i];
		sum = 0.0;
		
    		for (j = row[i]; j< row[i+1];j++) {
			k=col[j];
	 		sum += val[indval] * vect[k];
			indval++;
    		}
		result[i]  = sum ;
	}
}






void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
   
   int o = 0;
 
   mwIndex n = 0;

  
   mwIndex *Ap,
          *Ai;
 
   double * Ax;

  double      * y;

/*
  if (nrhs != 2 ) {
      mexErrMsgIdAndTxt("ccs:rhs","\n\r bad number of arguments 2 expected");
  }
  
  if (!mxIsSparse(prhs[0])) {
      mexErrMsgIdAndTxt("ccs:rhs","\n\r arg 1 must be a sparse matrix");
  }  
 */
  Ax = mxGetPr(prhs[0]);
  Ai = mxGetIr(prhs[0]);
  Ap = mxGetJc(prhs[0]);
  
  n   = mxGetN(prhs[0]);
 
  double * X = mxGetPr(prhs[1]);


  /* first return value : double matrix of size NbPxNbP  */

  plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
  y  = mxGetPr(plhs[0]);


  mult(Ap, Ai, Ax, X, y ,n );

}


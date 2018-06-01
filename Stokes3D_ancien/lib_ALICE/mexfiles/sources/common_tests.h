/*
 *  Common tests used in the OpenNL mexfiles.
 *
 *  April 2008
 *
 *  Contributors:   Nicolas SAUGNIER - INRIA - ALICE Project Team
 *
 *
 */
#include "mex.h"

// define OPENNL_NO_TESTS when compiling if you don't want all the tests

#ifdef OPENNL_NO_TESTS

#define TestNumericRealNotEmpty(array,object)              
#define TestVector(array,object)                           
#define TestNumber(array,object)                          
#define TestValidMatrixID(array,object)                   
#define TestNbrArgs(nbArgs,nbMin,nbMax)                    
#define TestVectStrictPositive(vect,nbElements,object)     
#define TestVectRange(vect,nbElements,nMin,nMax,object)    

#else

#define TestNumericRealNotEmpty(array,object)              DTestNumericRealNotEmpty(array,object)
#define TestVector(array,object)                           DTestVector(array,object)
#define TestNumber(array,object)                           DTestNumber(array,object)
#define TestValidMatrixID(array,object)                    DTestValidMatrixID(array,object)
#define TestNbrArgs(nbArgs,nbMin,nbMax)                    DTestNbrArgs(nbArgs,nbMin,nbMax)
#define TestVectStrictPositive(vect,nbElements,object)     DTestVectStrictPositive(vect,nbElements,object)
#define TestVectRange(vect,nbElements,nMin,nMax,object)    DTestVectRange(vect,nbElements,nMin,nMax,object)

#endif



/*
    In all functions, 'object' is a string used to represent the tested object in the messages.
 */

// Tests if the array is numeric, real and not empty
void DTestNumericRealNotEmpty(const mxArray* array, const char* object){
    if ( mxIsEmpty(array) || mxIsComplex(array) || !mxIsNumeric(array) ){
        mexErrMsgIdAndTxt("OpenNL:Test:NotNumericRealNotEmpty","%s must be non-empty, numeric and real.", object);
    }    
}

// Tests if the array is a vector of scalar
void DTestVector(const mxArray* array, const char* object){
    DTestNumericRealNotEmpty(array, object);
    if ( mxGetM(array)!=1 && mxGetN(array)!=1 ){
        mexErrMsgIdAndTxt("OpenNL:Test:NotVector","%s must be a vector (at least one of its dimentions should be 1).", object);
    }
}

// Tests if the array is a scalar
void DTestNumber(const mxArray* array, const char* object){
    DTestNumericRealNotEmpty(array, object);
    if ( mxGetNumberOfElements(array) > 1 ){
        mexWarnMsgIdAndTxt("OpenNL:Test:NotNumber","%s has more than one elements. The first one will be used.", object);
    }
}

// Tests if the array is valide matrix ID
void DTestValidMatrixID(const mxArray* array, const char* object){
    DTestNumber(array, object);
    if ( !SparseMatrixExists((int)mxGetScalar(array)) ){
        mexErrMsgIdAndTxt("OpenNL:Test:NotValidMatrixID","The specified sparse matrix does not exist.");
    }
}

// Tests number of args. Must be between min and max, inclusive.
void DTestNbrArgs(int nbArgs, int nbMin, int nbMax){
    if (nbArgs<nbMin){
        mexErrMsgIdAndTxt("OpenNL:Test:NotEnoughArgs","Not enough input arguments");
    }
    if (nbArgs>nbMax){
        mexErrMsgIdAndTxt("OpenNL:Test:TooManyArgs","Too many input arguments.");
    }
}

// Tests if all elements are strict positive
void DTestVectStrictPositive( double* vect, int nbElements, const char* object){
    for (int i=0; i<nbElements; i++){
        if( vect[i] <= 0 ){
           mexErrMsgIdAndTxt("OpenNL:Test:VectorNotStrictPositiv", "Element %d of %s is out of range : must be strict positive.", i+1, object);
        }
    }
}

// Tests if all elements are between min and max, inclusive.
template <class T>
void DTestVectRange( T* vect, int nbElements, int nMin, int nMax , const char* object){
    for (int i=0; i<nbElements; i++){
        if( static_cast<int>(vect[i]) < nMin || static_cast<int>(vect[i]) > nMax ){
           mexErrMsgIdAndTxt("OpenNL:Test:VectorNotInRange", "Element %d of %s is out of range.", i+1, object);
        }
    }
}


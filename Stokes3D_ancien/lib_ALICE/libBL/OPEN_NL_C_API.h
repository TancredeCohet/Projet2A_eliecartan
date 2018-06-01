#ifndef OPEN_NL_C_API_H
#define OPEN_NL_C_API_H

#include <stdlib.h>

#include "sparse_matrix.h"
#include "geex_defs.h"
#include "oNL_types.h"


#ifdef __cplusplus
extern "C" {  // only need to export C interface if
              // used by C++ source code
#endif

// Create a new Sparse Matrix
MatrixID_T NewSparseMatrix(gx_size_T m,gx_size_T n);

// Copy a Sparse Matrix
MatrixID_T CopySparseMatrix(MatrixID_T MatrixID);

// Delete a Sparse Matrix
void DeleteSparseMatrix(MatrixID_T matrixID);


// Delete all SparseMatrix
void ClearSparseMatrix(void);


// Add a coefficient to a matrix
void AddSparseMatrix(MatrixID_T matrixID,gx_index_T i,gx_index_T j,double value);


// Get matrix NNZ
gx_size_T NnzSparseMatrix(MatrixID_T matrixID);


// Get matrix' N
gx_size_T SparseMatrixN(MatrixID_T matrixID);


// Get matrix' M
gx_size_T SparseMatrixM(MatrixID_T MatrixID);


// Get a matrix' coefficient
double coefficient(MatrixID_T matrixID, gx_index_T i, gx_index_T j);


// Convert a Sparse Matrix to a Matlab matrix
void SparseMatrixToMatlab(MatrixID_T matrixID, size_t* irs, size_t* jcs, double* sr);


// put 0 on the row (resp. column) of the matrix at the given index. 
// indexes tab must be sorted !
void RemoveRowsSparseMatrix(MatrixID_T MatrixID, gx_size_T nbi, gx_index_T* indexes);
void RemoveColsSparseMatrix(MatrixID_T MatrixID, gx_size_T nbi, gx_index_T* indexes);


// Check the existance of a matrix
bool SparseMatrixExists(MatrixID_T MatrixID); 


// Matrix-vector  multiplication : A.x = y
void VectorMultSparseMatrix(MatrixID_T MatrixID, const double* x, double* y ); 

// Matrix-vector  multiplication : At.x = y
void VectorMultTransposeSparseMatrix(MatrixID_T MatrixID, const double* x, double* y ); 


// Matrix-Matrix  multiplication : C = A.B
void MultABSparseMatrix(MatrixID_T A_ID, MatrixID_T B_ID, MatrixID_T C_ID ); 

// Matrix-Matrix  multiplication : C = At.B
void MultAtBSparseMatrix(MatrixID_T A_ID, MatrixID_T B_ID, MatrixID_T C_ID ); 

// Matrix-Matrix  multiplication : C = A.Bt
void MultABtSparseMatrix(MatrixID_T A_ID, MatrixID_T B_ID, MatrixID_T C_ID ); 

// Matrix-Matrix  multiplication : C = At.Bt
void MultAtBtSparseMatrix(MatrixID_T A_ID, MatrixID_T B_ID, MatrixID_T C_ID ); 


#ifdef __cplusplus
}
#endif


#endif

%Mexfile VectorMultSparseMatrix, OpenNL, 2008.
%
%VectorMultSparseMatrix(MatrixID,X,Transpose) returns the result of the
%   multiplication of the X vector by the given matrix. If Transpose is
%   true, the transpose of the matrix is multiplied by X.
%
%   X can be either column or row vector. The result will have the same
%   orientation. If X is sparse, then it must be a column vector.
%   Transpose is optionnal, and can be either a logical scalar or a real
%   scalar. Non-zero values are considered true. Default is false.
%
%Exemple:
%   A=CreateSparseMatrix(2,3);
%   AddSparseMatrix(A, 2, 2, 1)
%   AddSparseMatrix(A, 1, 1, 2)
%   AddSparseMatrix(A, 2, 1, 3)
%   AddSparseMatrix(A, 1, 3, 4)
%   B=SparseMatrixToMatlab(A); full(B)
%
%   ans =
%       2     0     4
%       3     1     0
%
%
%   VectorMultSparseMatrix(A, [0.1 10], 1 )
%
%   ans =
%       30.2000   10.0000    0.4000
%
%   VectorMultSparseMatrix(A, [1 ; 10; 100] )
%
%   ans =
%       402
%        13
%


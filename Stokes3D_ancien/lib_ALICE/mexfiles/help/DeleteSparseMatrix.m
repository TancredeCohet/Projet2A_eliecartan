%Mexfile DeleteSparseMatrix, OpenNL, 2008.
%
%DeleteSparseMatrix(MatrixID) deletes a sparse matrix.
%   
%   The allocated memory is freed and the slot it used is available again.
%
%Exemple:
%   A=CreateSparseMatrix(2,2)
%
%   A =
%       0
%
%   B=CreateSparseMatrix(5,5)
%
%   B =
%       1
%
%   DeleteSparseMatrix(A);
%   C=CreateSparseMatrix(3,3)
%
%   C =
%       0
%
%   D=CreateSparseMatrix(4,4)
%
%   D =
%       2
%

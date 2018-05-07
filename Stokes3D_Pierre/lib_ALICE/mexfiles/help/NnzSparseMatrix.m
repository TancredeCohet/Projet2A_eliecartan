%Mexfile NnzSparseMatrix, OpenNL, 2008.
%
%NnzSparseMatrix(MatrixID) is the Number of Non Zero coefficients in a matrix.
%
%
%Exemple:
%   A=CreateSparseMatrix(3,3);
%   AddSparseMatrix(A,2,1);
%   NnzSparseMatrix(A)
%   
%   ans =
%       1
%
%   AddSparseMatrix(A,3,2);
%   AddSparseMatrix(A,1,1);
%   NnzSparseMatrix(A)
%   
%   ans =
%       3
%


%Mexfile coefficient, OpenNL, 2008.
%
%coefficient(MatrixID,I,J) is the (I,J) coefficient of a matrix.
%   I and J should be scalar.
%
%Exemple:
%   A=CreateSparseMatrix(2,2);
%   AddSparseMatrix(A, 1, 2, 3.14);
%   coefficient(A,1,1)
%   
%   ans =
%           0
%
%   coefficient(A,1,2)
%
%   ans =
%      3.1400
%

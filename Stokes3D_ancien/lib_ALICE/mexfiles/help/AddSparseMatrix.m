%Mexfile AddSparseMatrix, OpenNL, 2008.
%
%AddSparseMatrix(MatrixID,I,J,VAL) adds the scalar VAL at position (I,J) to
%   a matrix.
%
%   I, J and VAL should be scalar.
%
%Exemple:
%   A=CreateSparseMatrix(4, 4);
%   AddSparseMatrix(A, 2, 3, 3.14);
%   RES=SparseMatrixToMatlab(A);
%   full(RES)
%       
%   RES = 
%     0       0       0       0
%     0       0  3.1400       0
%     0       0       0       0
%     0       0       0       0
%

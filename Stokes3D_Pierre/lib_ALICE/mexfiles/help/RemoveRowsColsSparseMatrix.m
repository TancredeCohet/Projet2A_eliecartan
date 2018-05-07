%Mexfile RemoveRowsColsSparseMatrix, OpenNL, 2008.
%
%RemoveRowsColsSparseMatrix(MatrixID,I) removes all the non-zero
%   coefficients of the rows and columns given in I, and sets the 
%   corresponding elements on the diagonal to 1.
%
%   I must be a numeric full vector.
%
%B2=RemoveRowsColsSparseMatrix(MatrixID,I,B) does the same to the matrix,
%   but returns a copy of the B vector where all the values at the indexes
%   given in I are zeroed. 
%
%   B must be a numeric vector at least as long as the smallest dimension of
%   the matrix. If B is sparse, it must be a column vector.
%
%Exemple:
%   A=CreateSparseMatrix(5, 5);
%   AddSparseMatrix(A, 1, 1, 1.1);
%   AddSparseMatrix(A, 1, 5, 4.567);
%   AddSparseMatrix(A, 2, 3, 3.14);
%   AddSparseMatrix(A, 4, 1, 1.5);
%   AddSparseMatrix(A, 4, 4, 2.7);
%   AddSparseMatrix(A, 5, 4, 8);
%   M1=SparseMatrixToMatlab(A);
%   full(M1)
%       
%   M1 = 
%1.1000       0       0       0  4.5670
%     0       0  3.1400       0       0
%     0       0       0       0       0
%1.5000       0       0  2.7000       0
%     0       0       0  8.0000       0
%   
%   B=[1 2 3 4 5];
%
%   B2=RemoveRowsColsSparseMatrix(A,[1 5],B)
%
%   B2 =
%       0 2 3 4 0
%
%   M2=SparseMatrixToMatlab(A);
%   full(M2)
%       
%   M2 = 
%1.0000       0       0       0       0
%     0       0  3.1400       0       0
%     0       0       0       0       0
%     0       0       0  2.7000       0
%     0       0       0       0  1.0000
%   







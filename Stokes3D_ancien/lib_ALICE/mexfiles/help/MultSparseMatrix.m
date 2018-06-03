%Mexfile MultSparseMatrix, OpenNL, 2008.
%
%VectorMultSparseMatrix(A_ID, B_ID, A_Transpose, B_Transpose) computes the 
%   product of the A and B matrix. Using A_Transpose and B_Transpose, you 
%   can get the products AtB, ABt and AtBt.
%
%   A_Transpose and B_Transpose are optionnal, and can be either logical 
%   or real scalars. Non-zero real values are considered true. Default is 
%   false.
%
%Exemple:
%   A=CreateSparseMatrix(2,3);
%   AddSparseMatrix(A, 2, 2, 1)
%   AddSparseMatrix(A, 1, 1, 2)
%   AddSparseMatrix(A, 2, 1, 3)
%   AddSparseMatrix(A, 1, 3, 4)
%   a=SparseMatrixToMatlab(A); full(a)
%
%   ans =
%       2     0     4
%       3     1     0
%
%   B=CreateSparseMatrix(2,3);
%   AddSparseMatrix(B, 2, 3, 0.1)
%   AddSparseMatrix(B, 1, 1, 0.3)
%   AddSparseMatrix(B, 2, 1, 0.5)
%   AddSparseMatrix(B, 1, 2, 0.7)
%   b=SparseMatrixToMatlab(B); full(b)
%   
%   ans =
%
%    0.3000    0.7000         0
%    0.5000         0    0.1000
%    
%   C1=MultSparseMatrix(A, B, true, false);
%
%   c1=SparseMatrixToMatlab(C1);full(c1)
%
%   ans =
%
%    2.1000    1.4000    0.3000
%    0.5000         0    0.1000
%    1.2000    2.8000         0
%
%   C2=MultSparseMatrix(A, B, 0, 1);
%
%   c2=SparseMatrixToMatlab(C2);full(c2)
%
%   ans =
%
%    0.6000    1.4000
%    1.6000    1.5000
%


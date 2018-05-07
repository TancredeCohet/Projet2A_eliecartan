%Mexfile AddMatElem, OpenNL, 2008.
%
%AddMatElem(MatrixID,I,J,VAL) adds the elements of the VAL to a matrix,
%   at the indexes given by I and J.
%   
%   I and J can be either column or row full vectors.
%   The first dimension of VAL must match the dimension of I, and its
%   second dimension must match the dimension of J.
%   VAL must be a numeric matrix (either sparse or full).
%
%Exemple:
%   A=CreateSparseMatrix(4, 4);
%   AddMatElem(A, [1 4], [2 3], [6 7;8 9]);
%   RES=SparseMatrixToMatlab(A);
%   full(RES)
%       
%   RES = 
%     0     6     7     0
%     0     0     0     0
%     0     0     0     0
%     0     8     9     0
%

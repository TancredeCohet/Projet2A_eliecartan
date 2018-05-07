function []=printSparseMatrix(MatrixID)
[m,n]=SizeSparseMatrix(MatrixID);
fprintf('[');
for i=1:m
  for j=1:n
   coeff=coefficient(MatrixID,i,j);
   fprintf('%f ',coeff);
  end
  fprintf('\n');
end  
fprintf(']');

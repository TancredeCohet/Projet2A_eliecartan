% Test_AddMatElem.m, OpenNL, 2008.
%
% Unit test for AddMatElem mexfile.
%

% $Revision : 1.0.0 $


Test_A=CreateSparseMatrix(4, 4) ;
AddMatElem(Test_A, [1 4], [2 3], [6 7;8 9]);
Test_RES=SparseMatrixToMatlab(Test_A);

Test_COMP = zeros(4,4);
Test_COMP(1,2)=6;
Test_COMP(1,3)=7;
Test_COMP(4,2)=8;
Test_COMP(4,3)=9;

Test_Z = not( Test_RES == Test_COMP );

if( sum(sum(Test_Z)) )
   disp('TEST AddMatElem 1 : ERREUR : le resultat attendu est :')
   Test_RES
end


ClearSparseMatrix;


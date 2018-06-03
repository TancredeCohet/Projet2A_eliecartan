function [] = libBL_build( compileOptions )
%***********************************************************************
% This function aims at building the binaries required by the sushi3d script
% it builds the libBL library distributed by the ALICE INRIA team 
% and compile all the required mexfiles including those who handle VTK functions

	if (~exist('compileOptions','var'))
		compileOptions = ' ';
	end

	if (ismac)
		lib_file ='libBL.dylib';
	else
	    	lib_file ='libBL.so';
	end

	% On cherche la librairie dans le repertoire courant
	if (~exist(lib_file,'file'))
		fprintf('No library for sparse matrices in this directory, copy from ../libBL\n');
		% On cherche la librairie dans le repertoire ALICE
		if ~exist(strcat('../libBL/',lib_file),'file')
			fprintf('No BL in ../libBL/, we need to compile it\n');
			% construction de la bibliotheque libBL
			cd '../libBL';
			bl_build;
			if (ismac)
				command = sprintf('install_name_tool -id @rpath/%s %s',lib_file,lib_file);
				system(command);
			end
			cd '../mexfiles';
		end
		% on copie la bibliotheque
		command=strcat('cp ../libBL/',lib_file,' sources');
		system(command);
	end
	compile_mex;
	% on les copie dans bin
	command='mv sources/*.mex* bin';
	system(command);
    
end

function [] = compile_mex()

cd 'sources';
fprintf('Compiling mexfiles \n');

pathlibBL=fullfile(pwd, '../../libBL');
IpathlibBL=['-I' pathlibBL];
LpathlibBL=['-L' pathlibBL];
LDFLAGS=['LDFLAGS=\$LDFLAGS -Xlinker -rpath -Xlinker ', fullfile(pwd)];

mex(LDFLAGS,'-largeArrayDims','-I.',IpathlibBL,LpathlibBL,'-lbBL','AddMatElem.cpp');
mex(LDFLAGS,'-largeArrayDims','-I.',IpathlibBL,LpathlibBL,'-lBL','AddSparseMatrix.cpp');
mex(LDFLAGS,'-largeArrayDims','-I.',IpathlibBL,LpathlibBL,'-lBL','ClearSparseMatrix.cpp');
mex(LDFLAGS,'-largeArrayDims','-I.',IpathlibBL,LpathlibBL,'-lBL','coefficient.cpp');
mex(LDFLAGS,'-largeArrayDims','-I.',IpathlibBL,LpathlibBL,'-lBL','CopySparseMatrix.cpp');
mex(LDFLAGS,'-largeArrayDims','-I.',IpathlibBL,LpathlibBL,'-lBL','CreateSparseMatrix.cpp');
mex(LDFLAGS,'-largeArrayDims','-I.',IpathlibBL,LpathlibBL,'-lBL','DeleteSparseMatrix.cpp');
mex(LDFLAGS,'-largeArrayDims','-I.',IpathlibBL,LpathlibBL,'-lBL','MultSparseMatrix.cpp');
mex(LDFLAGS,'-largeArrayDims','-I.',IpathlibBL,LpathlibBL,'-lBL','NnzSparseMatrix.cpp');
mex(LDFLAGS,'-largeArrayDims','-I.',IpathlibBL,LpathlibBL,'-lBL','RemoveRowsColsSparseMatrix.cpp');
mex(LDFLAGS,'-largeArrayDims','-I.',IpathlibBL,LpathlibBL,'-lBL','SetBoundaryConditions.cpp');
mex(LDFLAGS,'-largeArrayDims','-I.',IpathlibBL,LpathlibBL,'-lBL','SizeSparseMatrix.cpp');
mex(LDFLAGS,'-largeArrayDims','-I.',IpathlibBL,LpathlibBL,'-lBL','SolveSparseMatrix.cpp');
mex(LDFLAGS,'-largeArrayDims','-I.',IpathlibBL,LpathlibBL,'-lBL','SparseMatrixToMatlab.cpp');
mex(LDFLAGS,'-largeArrayDims','-I.',IpathlibBL,LpathlibBL,'-lBL','VectorMultSparseMatrix.cpp');

cd '../';

end

% % contains a list of the usual paraview include path
% paraviewCommonIncPaths = { '/usr/include/paraview',...
% 			'/',...
% 			'/home',...
% 			'/usr'};
% 
% %find which one corresponds to a valid paraview inc path
% 
% paraviewIncPath = getPath(paraviewCommonIncPaths,'/vtkPolyData.h','include/paraview');
% 
% % contains a list of the usual paraview lib path
% paraviewCommonLibPaths = { '/usr/lib/paraview',...
% 			'/',...
% 			'/home',...
% 			'/usr'};
% %find which one corresponds to a valid paraview lib path
% 
% paraviewLibPath = getPath(paraviewCommonLibPaths,'libvtkCommon.so','lib/paraview');




% fprintf('Compiling mexfiles \n');
% setenv('LIBPATH',pwd)
% 
% %setenv('PARALIBPATH',paraviewLibPath)
% %setenv('PARAINCPATH',paraviewIncPath)
% 
% mex LDFLAGS="\$LDFLAGS -Xlinker -rpath -Xlinker \$LIBPATH" -largeArrayDims  AddMatElem.cpp -v -I./ -I../libBL  -L. -lBL; 
% mex LDFLAGS="\$LDFLAGS -Xlinker -rpath -Xlinker \$LIBPATH" -largeArrayDims  AddSparseMatrix.cpp  -I./ -I../libBL -L.  -lBL; 
% mex LDFLAGS="\$LDFLAGS -Xlinker -rpath -Xlinker \$LIBPATH" -largeArrayDims  ClearSparseMatrix.cpp  -I./ -I../libBL -L.  -lBL; 
% mex LDFLAGS="\$LDFLAGS -Xlinker -rpath -Xlinker \$LIBPATH" -largeArrayDims  coefficient.cpp  -I./ -I../libBL  -L. -lBL; 
% mex LDFLAGS="\$LDFLAGS -Xlinker -rpath -Xlinker \$LIBPATH" -largeArrayDims  CopySparseMatrix.cpp  -I./ -I../libBL  -L. -lBL; 
% mex LDFLAGS="\$LDFLAGS -Xlinker -rpath -Xlinker \$LIBPATH" -largeArrayDims  CreateSparseMatrix.cpp  -I./ -I../libBL  -L. -lBL; 
% mex LDFLAGS="\$LDFLAGS -Xlinker -rpath -Xlinker \$LIBPATH" -largeArrayDims  DeleteSparseMatrix.cpp  -I./ -I../libBL  -L. -lBL; 
% mex LDFLAGS="\$LDFLAGS -Xlinker -rpath -Xlinker \$LIBPATH" -largeArrayDims  MultSparseMatrix.cpp  -I./ -I../libBL  -L. -lBL; 
% mex LDFLAGS="\$LDFLAGS -Xlinker -rpath -Xlinker \$LIBPATH" -largeArrayDims  NnzSparseMatrix.cpp  -I./ -I../libBL  -L. -lBL; 
% mex LDFLAGS="\$LDFLAGS -Xlinker -rpath -Xlinker \$LIBPATH" -largeArrayDims  RemoveRowsColsSparseMatrix.cpp -I./ -I../libBL  -L. -lBL;
% mex LDFLAGS="\$LDFLAGS -Xlinker -rpath -Xlinker \$LIBPATH" -largeArrayDims  SetBoundaryConditions.cpp  -I./ -I../libBL  -L. -lBL; 
% mex LDFLAGS="\$LDFLAGS -Xlinker -rpath -Xlinker \$LIBPATH" -largeArrayDims  SizeSparseMatrix.cpp  -I./ -I../libBL  -L. -lBL; 
% mex LDFLAGS="\$LDFLAGS -Xlinker -rpath -Xlinker \$LIBPATH" -largeArrayDims  SolveSparseMatrix.cpp  -I./ -I../libBL  -L. -lBL; 
% mex LDFLAGS="\$LDFLAGS -Xlinker -rpath -Xlinker \$LIBPATH" -largeArrayDims  SparseMatrixToMatlab.cpp  -I./ -I../libBL  -L. -lBL; 
% mex LDFLAGS="\$LDFLAGS -Xlinker -rpath -Xlinker \$LIBPATH" -largeArrayDims  VectorMultSparseMatrix.cpp  -I./ -I../libBL  -L. -lBL; 
% %mex LDFLAGS="\$LDFLAGS -Xlinker -rpath -Xlinker \$PARALIBPATH" -largeArrayDims  IsInside_VTK.cpp -v -I. -I"\$PARAINCPATH" -L"\$PARALIBPATH" -lvtkFiltering -lvtkCommon -lvtkGraphics; 
% %mex LDFLAGS="\$LDFLAGS -Xlinker -rpath -Xlinker \$PARALIBPATH" -largeArrayDims  FindTetra_VTK.cpp -v -I. -I"\$PARAINCPATH" -L"\$PARALIBPATH" -lvtkFiltering -lvtkCommon -lvtkGraphics; 
% mex -v -O -largeArrayDims CFLAGS="\$CFLAGS -fopenmp" CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" CCS.c;


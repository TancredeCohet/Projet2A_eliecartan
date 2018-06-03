function [] = bl_build( compileOptions )

if (~exist('compileOptions','var'))
    compileOptions = ' ';
end

if (ismac)
    lib_file ='libBL.dylib';
elseif (isunix)
    lib_file ='libBL.so';
end

if (~exist(lib_file,'file'))
    fprintf('compiling BL library\n')
    set_makefile;
    make;
end

function [] = make()
% calling make
command='make';
system(command);

function [] = set_makefile()
% select the suited Makefile.inc according to the arch
if (ispc)
   error('system not supported');
end
if (ismac)  
    inc_file = 'Makefile_MacOS64.inc';
else % c'est forcement un linux
    inc_file = 'Makefile_Linux.inc';
end
system('rm -f Makefile.inc');
command = ['ln -s',' ',inc_file,' ', 'Makefile.inc'];
system(command);

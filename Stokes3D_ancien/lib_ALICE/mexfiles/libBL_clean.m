function []=libBL_clean
fprintf('cleaning mexfiles\n');
command = 'rm -f bin/*';
system(command);
% clean library
if (ismac)
    lib_file ='libBL.dylib';
end
if (isunix)
    lib_file ='libBL.so';
end
fprintf('cleaning library\n');
command = sprintf('rm  -f %s',lib_file);
system(command);
cd '../libBL';
bl_clean;
cd '../mexfiles';

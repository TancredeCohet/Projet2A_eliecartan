function []=bl_clean
command = 'make clean_lib';
system(command);
command = 'make clean';
system(command);
command = 'rm -f Makefile.inc';
system(command);

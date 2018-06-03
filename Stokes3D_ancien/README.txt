Building ALICE/OPEN_NL library (libBL) and compiling mex-files for the 
efficient management of sparse matrices.
========================================================================

1) Prerequisites:

    - g++ > 4.2.0 

2) Start matlab 

   $ matlab 

3) Build MEX-Files (libBL + mex-functions)
   >> cd lib_ALICE/mexfiles
   >> libBL_clean
   >> libBL_build

if any problem occurs involving the libstdc++.so.6.0.10 just go in your 
matlab installation directory and remove 

matlab/sys/os/glnxa64/libstdc++.so.6.0.10

 


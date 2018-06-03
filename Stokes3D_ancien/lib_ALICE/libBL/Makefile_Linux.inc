# C config 
CXX=g++
OPENMP=-fopenmp
CXXFLAGS= -O3 -fPIC  $(OPENMP) -DNDEBUG
LDNAME=LD_LIBRARY_PATH
LIBEXT=so
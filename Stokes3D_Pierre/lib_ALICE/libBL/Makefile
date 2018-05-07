include Makefile.inc

#objects
objects =OPEN_NL_C_API.o  sparse_matrix.o  

#rules for library compiling
libBL.$(LIBEXT) : $(objects)
	$(CXX) $(OPENMP) -fPIC -shared -o $@ $(objects)	

.cpp.o :
	$(CXX) $(CXXFLAGS) -I../mexfiles/ -c $<

OPEN_NL_C_API.o : OPEN_NL_C_API.cpp 

sparse_matrix.o : sparse_matrix.cpp 

test_OpenNLCAPI.o : test_OpenNLCAPI.cpp 

# test target
check : libBL.$(LIBEXT) test_OpenNLCAPI.o
	$(CXX) -o test_BL test_OpenNLCAPI.o -L. -lBL

# clean targets
clean_all :: clean clean_lib clean_tests 

clean_lib:
	@echo cleaning lib 
	@rm -f libBL.$(LIBEXT)
clean:
	@echo cleaning object files 
	@rm -f *.o 
clean_tests:
	@rm -f test_BL 

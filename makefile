################################################################################
#                                                                              #
#   makefile						         	       #
#   SParse Quadrature Routines (SPQR)				               # 
#                                                                              #
#   created 9/11/2016                                                          #
#                                                                              #
################################################################################
#
VERSION=1.0
SHELL=/bin/sh

# set your C++ compiler
CPP= g++
# set your MEX compiler
MEX= mex

# set C++ compile flags
CPPFLAGS= -O3 -fPIC -Wall -W -pedantic
# set Eigen include path
EIGENPATH= /usr/include/eigen3/
#
#

OBJ = MXsparseQuadrature.cpp SparseQuadrature.o TensorProductQuadrature.o\
		TDindexSet.o

all:			$(OBJ)
				$(MEX) $(OBJ) -I$(EIGENPATH)
				

# tell make how to create a .o from a .cpp

%.o:%.cpp
	$(CPP)  $(CPPFLAGS) -I$(EIGENPATH) -c $<


.PHONY: clean

clean:
	rm -f *.o

# DO NOT DELETE THIS LINE -- make depends on it.

# NOTE ONly the following must be prexist
#
#./tetgen/*
#
# you can download updated versions from internet or just 
# use the current ones.
# if library libtet.a is not in the above folder the current makefile
# will automatically create it! 
#
# 
%.o: %.f90
	$(F90) $(F90_CFLAGS) $^ -o $@
%.o: %.cxx
	$(CPPC) $(CPPC_CFLAGS) -I$(INCS) $^ -o $@  

# COMPILER_PATH = /opt/gcc/bin/
COMPILER_PATH =
F90 = $(COMPILER_PATH)gfortran
# F90 = $(COMPILER_PATH)ifort
#F90_CFLAGS = -c -O3
F90_CFLAGS = -c -Wall
CPPC = $(COMPILER_PATH)g++
# CPPC = $(COMPILER_PATH)icpc
CPPC_CFLAGS = -c -DTETLIBRARY
LIBS =  ./tetgen/
INCS =  $(LIBS)
FCLFLAGS = 
LFLAGS = -ltet  -lstdc++
# LFLAGS = -ltet
TETGEN_LIB = $(LIBS)libtet.a

curved.run: $(TETGEN_LIB) lag_basis.o tetgen_wrapper.o tetmesher.o tet_props.o curved_tet.o
	$(F90) $(FCLFLAGS) -I$(INCS) -L$(LIBS) lag_basis.o tetgen_wrapper.o tetmesher.o tet_props.o curved_tet.f90 $(LFLAGS)

clean:
	rm -f *.o *.mod *.out *~ dumped*

$(TETGEN_LIB):
	(cd ./tetgen/; $(MAKE) tetlib )



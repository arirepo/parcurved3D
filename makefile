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
%.o: %.f
	$(F77) -c $^ -o $@

%.o: %.f90
	$(F90) $(F90_CFLAGS) $^ -o $@
%.o: %.cxx
	$(CPPC) $(CPPC_CFLAGS) -I$(INCS) -I$(INCS2) $^ -o $@  

# COMPILER_PATH = /opt/gcc/bin/
COMPILER_PATH = 
# COMPILER_PATH = /home/aghasemi/Desktop/openmpi-1.10.0/install_dir/bin/
# F90 = $(COMPILER_PATH)gfortran
# F90 = $(COMPILER_PATH)mpifort
F90 = $(COMPILER_PATH)mpif90
# F90 = $(COMPILER_PATH)ifort
#F90_CFLAGS = -c -O3
# F90_CFLAGS = -c -check all -g -CB
F90_CFLAGS = -c -W1
F77 = $(F90) 
# CPPC = $(COMPILER_PATH)g++
# CPPC = $(COMPILER_PATH)mpicxx
CPPC = $(COMPILER_PATH)mpic++
# CPPC = $(COMPILER_PATH)icpc
CPPC_CFLAGS = -c -DTETLIBRARY
LIBS =  ./tetgen/
INCS =  $(LIBS)
FCLFLAGS = 
# LFLAGS = $(LIBS2)lib* -ltet  -lstdc++ $(OMP_FLAG) -lmetis
LFLAGS = -locasdirty -ltet  -lstdc++ $(OMP_FLAG) -lmetis -lblas -llapack
# LFLAGS = -ltet
TETGEN_LIB = $(LIBS)libtet.a

# LIBS2 =  /home/aghasemi/Desktop/fortran_opencascade_wrapper/opencascade-6.9.0/install/lin64/gcc/lib/
#LIBS2 =  /home/aghasemi/Desktop/fortran_opencascade_wrapper/opencascade-6.9.1/tmp_o_same_gcc/
LIBS2 =  /home/aghasemi/Desktop/fortran_opencascade_wrapper/opencascade-6.9.1/tmp_o_rivermont/

# INCS2 =  /home/aghasemi/Desktop/fortran_opencascade_wrapper/opencascade-6.9.0/install/inc/
INCS2 =  /home/aghasemi/Desktop/fortran_opencascade_wrapper/opencascade-6.9.1/install_rivermont/inc/

# LIBS3 =  /home/aghasemi/metis-5.1.0/install_SUSE_13.1_64/lib/
# LIBS3 =  /home/aghasemi/metis-5.1.0/install_cerberus/lib/
# LIBS3 =  /home/aghasemi/metis-5.1.0/install_ibm_blue/lib/
LIBS3 =  /home/aghasemi/metis-5.1.0/install_intel_cerberus/lib/

# OMP_FLAG = -fopenmp
OMP_FLAG =

curved.run: $(TETGEN_LIB) tet_props.o var_array.o timing.o ocas_hooks.o op_cascade.o lag_basis.o tetgen_wrapper.o prism_mesher.o tetmesher.o mpi_comm_mod.o curved_tet.o
	$(F90) $(FCLFLAGS) -I$(INCS) -I$(INCS2) -L$(LIBS) -L$(LIBS2) -L$(LIBS3) tet_props.o var_array.o timing.o ocas_hooks.o op_cascade.o lag_basis.o tetgen_wrapper.o prism_mesher.o tetmesher.o mpi_comm_mod.o curved_tet.f90 $(LFLAGS)

vl.run: $(TETGEN_LIB) tet_props.o renka_trimesh_lib.o renka_trimesh.o gen_basis.o var_array.o ocas_hooks.o op_cascade.o lag_basis.o tetgen_wrapper.o prism_mesher.o tetmesher.o mpi_comm_mod.o
	$(F90) $(FCLFLAGS) -I$(INCS) -I$(INCS2) -L$(LIBS) -L$(LIBS2) -L$(LIBS3) tet_props.o renka_trimesh_lib.o renka_trimesh.o gen_basis.o var_array.o ocas_hooks.o op_cascade.o lag_basis.o tetgen_wrapper.o prism_mesher.o tetmesher.o mpi_comm_mod.o curved_prism.f90 $(LFLAGS)

clean:
	rm -f *.o *.mod *.smod *.out *~ dumped* *.tec $(LIBS)libtet.a *.txt tmp.m opencascade_faces.m

$(TETGEN_LIB):
	(cd ./tetgen/; $(MAKE) tetlib )



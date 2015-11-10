#!/bin/tcsh -f
#$ -N sobol
#$ -pe papertape_openmpi 200
#$ -S /bin/tcsh
##$ -hard -l virtual_free=3.5G,cpu_rem=3.50000
cd /home/aghasemi/Dropbox/delit/cerberus/curved_3d_elems
#$ -cwd
#$ -o stdio
#$ -e stderr
echo "JOB = $JOB_ID at `date`" >> Submitted
setenv MPICH_PROCESS_GROUP no

# setenv LD_LIBRARY_PATH "/home/aghasemi/Dropbox/delit/cerberus/curved_3d_elems/hack_lib/glibc2:$LD_LIBRARY_PATH"
# setenv LD_LIBRARY_PATH /home/aghasemi/Desktop/fortran_opencascade_wrapper/opencascade-6.9.0/install_cerberus/lin64/gcc/lib:${LD_LIBRARY_PATH}
# setenv LD_LIBRARY_PATH /usr/lib:/usr/lib64:${LD_LIBRARY_PATH}
# setenv PATH "/usr/bin:$PATH"
# setenv LD_LIBRARY_PATH /usr/local/lib64:${LD_LIBRARY_PATH}
# setenv LD_LIBRARY_PATH /usr/local/gcc/lib:/usr/local/gcc/lib64:${LD_LIBRARY_PATH}
setenv LD_LIBRARY_PATH /usr/local/gcc-4.3.1/lib:/usr/local/gcc-4.3.1/lib64:${LD_LIBRARY_PATH}

# open MPI
#setenv PATH "/home/aghasemi/Desktop/openmpi-1.10.0/install_cerberus/bin:$PATH"
#setenv LD_LIBRARY_PATH "/home/aghasemi/Desktop/openmpi-1.10.0/install_cerberus/lib:$LD_LIBRARY_PATH"

# gcc compiler
#setenv PATH "/home/aghasemi/gcc-trunk/install_cerberus/bin/:$PATH"

#setenv LD_LIBRARY_PATH "/lib64:/usr/lib64:$LD_LIBRARY_PATH"
#setenv LD_LIBRARY_PATH "/home/aghasemi/Dropbox/delit/cerberus/curved_3d_elems/hack_lib/glibc1:$LD_LIBRARY_PATH"


#setenv LD_LIBRARY_PATH "/home/aghasemi/gcc-trunk/install_cerberus/lib:/home/aghasemi/gcc-trunk/install_cerberus/lib64:$LD_LIBRARY_PATH"
#setenv LD_LIBRARY_PATH "/home/aghasemi/Dropbox/delit/cerberus/curved_3d_elems/hack_lib:$LD_LIBRARY_PATH"

# if `uname -i` != "i386" then
#    source /usr/local/intel/composerxe/bin/compilervars.csh intel64
# endif

echo "LD_LIBRARY_PATH:" >> Submitted
echo $LD_LIBRARY_PATH >> Submitted
echo "The following is the location of runtime on machine:"
which ifort
which mpirun

pwd
ldd "/home/aghasemi/Dropbox/delit/cerberus/curved_3d_elems/a.out"
#mpirun --mca pls_rsh_agent rsh -np 4  ./a.out
mpirun -np 200  ./a.out

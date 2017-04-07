############################################################################
# Makefile for GCOE server
############################################################################

#Define compiler
#  F_COMP=ifort
   F_COMP=mpif90
   
#Set compile options (Remove needless options after developing phase for )
   #-fast      :Fast option
   #-openmp    :Activate OpenMP indications
   #-traceback :Display line numer when the program was terminated abnormally
   #-warn all  :Display all error message during compiling
   OPT_COMPILES= -fast -warn all

#Link option 1
   #= -openmp :Activate OpenMP indications
   OPT_LINK1= 
   
#Link option 2
   #= -pg :Activate profiler
   OPT_LINK2=
   
#Define execution file name 
   EXE= go.out
   
#Define source files
   F_SRC=  modules.f90 \
           start_mpi.f90 \
           main_mpi.f90 \
           initialize.f90 \
           metabolic.f90 \
           output.f90 \
           physics.f90 \
           population_regu.f90 \
           spatial_calc.f90 \
           etc.f90
   
$(EXE): $(F_SRC)
	$(F_COMP) -c $(OPT_COMPILES) $(?)
	$(F_COMP) $(OPT_LINK1) -o $(EXE)  *.o $(OPT_LINK2)
   
# Definition for make clean
clean:
	rm -f $(EXE) *.o *.L gmon.out *.mod core

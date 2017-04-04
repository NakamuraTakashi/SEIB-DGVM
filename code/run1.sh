#
# Runscript for SEIB-DGVM:
#
# ------------------ COMPILE & RUN ----------------------------
#
 gfortran -fbounds-check mod_grid.f90 start_point.f90 -I/usr/include -L/usr/lib -lnetcdff -Wl,--stack,838860800 -fopenmp -O2 -o seib.exe
#gfortran -fbounds-check mod_grid.f90 start_point.f90 -I/usr/include -L/usr/lib -lnetcdff -O2 -fopenmp -o seib.exe
#gfortran -fbounds-check mod_grid.f90 start_point.F90 -I/usr/include -L/usr/lib -lnetcdff -O2 -o seib.exe
#
 export OMP_NUM_THREADS=12
#
 ./seib.exe
#

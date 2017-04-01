#
# Runscript for SEIB-DGVM:
#
# ------------------ COMPILE & RUN ----------------------------
#
 gfortran -fbounds-check mod_grid.f90 start_point.f90 -I/usr/include -L/usr/lib -lnetcdff -Wl,--stack,8388608 -O2 -o seib.exe
#
 ./seib.exe
#

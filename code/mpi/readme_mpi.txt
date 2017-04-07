!!!!! Important Note !!!!! 
The MPI related codes in this folder works with version 2.71.
Modifications will be required to work with version 2.80.


********************************************************************
Notes for SEIB-DGVM MPI-version (17/Nov/2011)

Hisashi SATO (JAMSTEC)
http://seib-dgvm.com/hsato/
********************************************************************

(1) This package enables the SEIB-DGVM to compute multigrid cells by using MPI parallel computation. This package also contains a R code for visualizing result of the multi-grid simulation. Ask your system administrator wether your environment satisfies this condition. You also have to prepare climate data. I do not distribute it online, because of the large size! If you want to get a copy of climatic data,  contact me, then I may I can help.

(2) Following files consist this package
  readme_mpi.txt       This document
  Makefile             Makefile for compiling program code
  go.bat               Specification for MPI jobs
  go.sh                Execusion script
  start_mpi.f90        The alternative of start.f90 in the one-grid simulation version
  main_mpi.f90         The alternative of main.f90 in the one-stand simulation version
  visualization.R      R code for mapping simulation results
  landmask_0.5deg.txt  Landmask for simulations @ 0.5 degree grid mesh
  landmask_1.0deg.txt  Landmask for simulations @ 1.0 degree grid mesh
  landmask_2.0deg.txt  Landmask for simulations @ 2.0 degree grid mesh
  index.html           A HTML file that shows visualized simulation results
  index_lai.html       A HTML file that shows visualized simulation results
 
(3) Following files must be modified to meet machine environment and simulation design. Refer comment lines in these files sequentially.
  Makefile
  go.bat
  go.sh
  start_mpi.f90
  main_mpi.f90
  parameter.f90

(4) For executing code, type following command sequantially. The first command compiles the code, and the secound command conducts simulation and visualization.
  makefile
  sh go.sh&

(5) If simulation year is longer than climatic data, this code use climate data repeatedly. To change this feature, modify folllowing line, which will be found in the main_mpi.f90.
  if (year_climate==YearMaxClimate+1) year_climate=1

For example, if you want to repeat climate data of last year when simulation year exceeds data lenght, change as follows:
  if (year_climate==YearMaxClimate+1) year_climate=YearMaxClimate

(6) Simulation assumes that atmospheric CO2 concentration is fixed at 368.0 ppm. To employ the historical record or IPCC A1B scenario, activate one of the following comment line, which will be found in the main_mpi.f90.
!co2atm = aco2_1750to2100(150+year_climate) !Starts from 1900
!co2atm = aco2_1750to2100(250+year_climate) !Starts from 2000
!co2atm = aco2_1750to2100(150)              !Fixed at 1900
!co2atm = aco2_1750to2100(350)              !Fixed at 2100

(7) Following parameters in parameter.txt will NOT be used for MPI simulations. Corresponding parameters are defined in start_mpi.f90.
  Fn_climate
  Fn_location
  Fn_spnin
  Fn_spnout

(8) In the default code, distribution maps are generated using average results during last 10 years of the simulation. To change this period for average, alter the value in parameter YearForMean in start_mpi.f90.

(9) Followings are the default simulation conditions
Simualtion grid. If you change this setting, you also have to change line 10 of visualization.R
-> 2.0deg

Simualtion area.  If you change this setting, you also have to change line 10 of visualization.R
-> Global

Simualtion year
-> 1901~2000

Atmospheric CO2 concentration
-> Fixed at 368ppm

Location where code and parameter files exist
-> '/home/hsato/work/'

Location where climatic data of 0.5deg exists
-> '/data/climate_2.0deg/

Location where result and restart files will be outputted
-> '/home/hsato/work/result/'

Location where summary results and visualized images will be outputted. If you change this setting, you also have to change line 10 of visualization.R
-> '/home/hsato/work/result_visualized/'

Number of processors to employ for MPI computation
-> 8 per each node * 5 nodes = total 40 processors
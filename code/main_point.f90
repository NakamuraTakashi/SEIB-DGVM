!******************************************************************************
!* Spatially Explicit Individual-Based Dynamic Global Vegetation Model        *
!* SEIB-DGVM ver. 3.01                                                        *
!*                                                                            *
!*   All rights are reserved by Dr. Hisashi SATO (JAMSTEC)                    *
!*                                                                            *
!*   Modified by                                                              *
!*        Drs. Tomo'omi KUMAHGAI & Taro NAKAI                                 *
!*        @Nagoya University                                                  *
!*                                                                            *
!*   The file SFLXALL_SRC_VER_2.7.1.f90 was the main program of the NOAH-LSM, *
!*   which is developped by US researchers. I obtained the code from the      *
!*   following URL, and modified for incoorprating into the SEIB-DGVM.        *
!*   ftp://ftp.emc.ncep.noaa.gov/mmb/gcp/ldas/noahlsm/                        *
!*                                                                            *
!*   Latest information and code using policy can be obtained at              *
!*   http://seib-dgvm.com/                                                    *
!******************************************************************************

!*************************************************************************************************
! MAIN SIMULATION LOOP
!*************************************************************************************************
SUBROUTINE main_loop ( &
           LAT, LON, GlobalZone, YearMaxClimate, YearMaxCO2, &
           tmp_air, tmp_air_range, prec, rad_short, rad_long, &
           wind, rh, tmp_soil, &
           aco2_annual, ALT, Albedo_soil0, W_fi, W_wilt, W_sat, W_mat, SoilClass)
   
!_____________ Global Variables
   USE data_structure
   USE time_counter
   USE vegi_status_current1
   USE vegi_status_current2
   USE grid_status_current1
   USE grid_status_current2
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none
   
!_____________ Set Local Parameters
   !Interval for intensive computation of radiation distribution (day)
   integer,parameter::Days_LightComp_Int = 14
   
!_____________ NOAH-LSM Integration (Phase 1: Set parameters for NOAH-LSM)
   INTEGER     ,parameter::NSOLD       = NumSoil !Max number of soil layer
   REAL        ,parameter::DT          = 1800.0  !Scound in a Time step (<=3600)
   INTEGER     ,parameter::TimeStepMax = 48      !Number of Time step in a day, 60*60*24/int(DT)
   
!Model Configuration:
   INTEGER,parameter:: NSOIL    = NumSoil   !Number of soil layers (2-20)
   REAL   ,parameter:: Z        = 6.0000    !Height (above ground) of the forcing wind vector(m)
   REAL   ,parameter,Dimension(1:NSOLD):: & !Thickness of each soil layer (m)
   SLDPTH = (/0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, &
              0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100 /) 
   
   INTEGER,parameter:: SLOPETYP = 1   !Slope type index 1-9 (default = 1)
!  INTEGER,parameter:: SLOPETYP = 9   !Activate to disable bottom runoff
   
   INTEGER,parameter:: ICE      = 0   !Sea ice flag (keep as integer 0 to designate non-sea)
   
!Other parameters
!MEMO by H.Sato: Thoese parametes are NOT used in simulations
   !Minimum value for Monthly SHDFAC (green vegetation fraction)
   REAL,parameter:: SHDMIN = 0
   !Veg. type specific threshold-snow-depth (in water equivalent m) that implies 100% snow cover
   REAL,parameter,Dimension(1:NumSoil):: SNUPX = (/  &
   0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.040, &
   0.040, 0.040, 0.040, 0.040, 0.040, 0.010, 0.013, 0.020, 0.013, 0.020 /)
   !A parameter for snow-distribution-shape
   REAL,parameter::SALP        = 4.0
   
!_____________ Set Augments
!Coordination
   real   ,intent(IN):: LAT !latitude  (degree)
   real   ,intent(IN):: LON !longitude (degree)
   
!ID number of global zone
   integer,intent(IN):: GlobalZone 
   
!Length of Inputted Climate and CO2 Data (yr)
   integer,intent(IN):: YearMaxClimate
   integer,intent(IN):: YearMaxCO2
   
!Climatic data
   real,dimension(Day_in_Year,YearMaxClimate),intent(IN)::&
   tmp_air       , & !Surface air temperature (Celcius)
   tmp_air_range , & !Daily range of tmp_air (Celcius)
   prec          , & !Precipitation (mm day-1)
   rad_short     , & !Shortwave radiation @ midday (W m-2)
   rad_long      , & !Daily mean of longwave radiation (W m-2)
   wind          , & !Wind velocity (m s-1)
   rh                !Relative humidity (%)
   
   real,dimension(Day_in_Year, YearMaxClimate, NumSoil),intent(IN):: &
   tmp_soil      !Soil temperature for each layers (Celcius)
   
!Atomospheric CO2 time-series @ ppm
!(1850~2000 from historical record, +2001~2100 from RCP8.5 scenario)
   real,dimension(YearMaxCO2),intent(IN)::aco2_annual !Atomospheric co2 concentration (ppm)
   
!Location data
   real:: &
   ALT          ,& !Altitude  (m above MSL)
   Albedo_soil0 ,& !Soil albedo
   W_fi         ,& !Filed capacity   (m3/m3, 0.0 -> 1.0)
   W_wilt       ,& !Wilting point    (m3/m3, 0.0 -> 1.0)
   W_sat        ,& !Saturate point   (m3/m3, 0.0 -> 1.0)
   W_mat           !Matrix potential 
   
   Integer,intent(IN):: SoilClass !Zobler(1986)'s slope type index 1-9
   
!_____________ Set Local variables

!Daily mean meteological variables
   real :: &
   tmp_air_Today       , & !Surface air temperature (Celcius)
   tmp_air_range_Today , & !Daily range of tmp_air  (Celcius)
   prec_Today          , & !Precipitation           (mm day-1)
   rad_short_Today     , & !Shortwave radiation @ midday (W m-2)
   rad_long_Today      , & !Daily mean of longwave radiation (W m-2)
   wind_Today          , & !Wind velocity     (m s-1)
   rh_Today                !Relative humidity (%)
   
   real,dimension(NumSoil):: &
   tmp_soil_Today=0.0   !Soil temperature for each layers (Celcius)
   
!Daily meteological variables at 30min time-resolution
   real,dimension(1:TimeStepMax):: &
   AirTmp_30min   , & !Air temperature               (C)
   Prec_30min     , & !Precipitation amount          (mm)
   RadShort_30min , & !Downward short wave radiation (W/m^2)
   RadLong_30min  , & !Downward long  wave radiation (W/m^2)
   rh_30min           !Relative humidity             (%)
   
   !Relative intensity of the short-wave radiation for each time-step of the day (0.0~1.0)
   !daily maximum value (midday value) == 1.0
   real,dimension(1:365, 1:TimeStepMax)::RadShortWeight
   
!Other variables
   integer year_climate !Current climate year
   real    cloud        !cloudness
   integer i, j, p, no  !Loop counters
   integer n, m         !For General Usage
   real    x, y         !For General Usage

!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
!For standard output
   integer time1
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add

   real,dimension(20)::out_timeseries
   real UnitConv
   
!_____________ NOAH-LSM Integration (Phase 2: Define augments)
!Parameter variables whose initial values will be given by the parameter file
   !Physical parameters:
   real::  TBOT      =0.0 !Annual constant bottom boundary soil temperature (K)
   real::  DF1_given =0.0 !THERMAL CONDUCTIVITY for top of the soil-layer (W m-1 k-1)
   
   !Vegation state
   integer:: VEGTYP=0  !Vegetation type index 1-13 (@after clear cut)
   real   :: SNOALB=0. !Max albedo over very deep snow (@after clear cut)
   
   !Initial state variables:
   real:: T1          =0. !Initial skin temperature (K)
   real:: STC(NSOLD)  =0. !SOIL TEMP (K)
   real:: SMC(NSOLD)  =0. !TOTAL SOIL MOISTURE CONTENT (VOLUMETRIC FRACTION)
   real:: SH2O(NSOLD) =0. !UNFROZEN SOIL MOISTURE CONTENT (VOLUMETRIC FRACTION)
   real:: CMC         =0. !Canopy water content (m)
   real:: SNOWH       =0. !Actual snow depth (m)
   real:: SNEQV       =0. !Water equiv snow depth (m)
   
   !Meteorological status inputed from forcing data
   real:: SFCSPD     = 0. !average wind vector speed at 6m (m/s)
   real:: RHumid     = 0. !relative humidity at 3 m        (%)
   real:: SOLDN      = 0. !Incoming short wave radiation   (W/m2)
   real:: LWDN       = 0. !Downward long wave radiation    (W/m2)
   real:: PRCP       = 0. !Precipitation (same as precip.rate in mm/sec)
   real:: SFCTMP     = 0. !Air temperature at 1st mdl lvl abv skin (K)
   real:: SFCPRS     = 0. !Pressure at 1st mdl lvl abv skin (PASCALS)
   
   !Loop counter related
   integer:: TimeStep =0  !Time step counter
   
   integer:: IMONTH   =0  !Month (1-12)
   integer:: IDAY     =0  !Day in the month (1-31)
   
   integer:: INDI=0     !For summing up time counter
   integer:: IIDAY = 0  !For summing up simulation days
   
   !Other augment
   integer:: NROOT = 0  !NO. OF SOIL LAYERS IN ROOT ZONE (1.LE.NROOT.LE.NSOIL) 
   
   real:: FLX1    = 0.0
   real:: FLX2    = 0.0
   real:: FLX3    = 0.0
   real:: EC      = 0.0
   real:: EDIR    = 0.0
   real:: ETT     = 0.0
   real:: DEW     = 0.0
   real:: FFROZP  = 0.0
   real:: RES     = 0.0
   real:: RUNOFF2 = 0.0
   real:: RUNOFF3 = 0.0
   real:: AET     = 0.0
   real:: ALB     = 0.0
   real:: ALBEDO  = 0.0
   REAL   DQSDT        
   real:: DQSDT2  = 0.0
   real:: ESAT    = 0.0
   real:: ETA     = 0.0
   real:: ETP     = 0.0
   real:: FUP     = 0.0
   real:: SHEAT   = 0.0
   real:: Q2      = 0.0
   real:: Q2SAT   = 0.0
   real:: RUNOFF1 = 0.0
   real:: SHDFAC  = 0.0
   real:: SNOMLT  = 0.0
   real:: SOILW   = 0.0
   real:: SOILM   = 0.0
   real:: S       = 0.0
   real:: T1V     = 0.0
   real:: T2V     = 0.0
   real:: TH2     = 0.0
   real:: TH2V    = 0.0
   real:: SNCOVR  = 0.0
   real:: XLAI    = 0.0
   real:: SOLNET  = 0.0
   real:: LVH2O   = 0.0
   
   !MEMO by H.Sato: Values in the following two variable does not change the results
   REAL   :: CH   = 1.E-4 !SFC EXCHANGE COEFFICIENT FOR HEAT/MOISTURE
   REAL   :: CM   = 1.E-4 !SFC EXCHANGE COEFFICIENT FOR MOMENTUM
   
   !MEMO by H.Sato: This subtoutine does not concern these variable except declare
   real:: BETA      = 0.0
   real:: DRIP      = 0.0
   real:: ET(NSOLD) = 0.0
   real:: ESNOW     = 0.0
   real:: RSMIN     = 0.0
   real:: RC        = 0.0
   real:: PC        = 0.0
   real:: RCS       = 0.0
   real:: RCT       = 0.0
   real:: RCQ       = 0.0
   real:: RCSOIL    = 0.0
   real:: SMCWLT    = 0.0
   real:: SMCMAX    = 0.0
   real:: SMCREF    = 0.0
   real:: SMCDRY    = 0.0
   real:: PTU       = 0.0
   
   !Daily sumupのために佐藤が追加した変数
   real:: ETP_sum     = 0.0 !潜在蒸発散量 (W M-2をmm/dayに変換した)
   real:: EC_sum      = 0.0 !CANOPY WATER EVAPORATION  (同上)
   real:: EDIR_sum    = 0.0 !DIRECT SOIL EVAPORATION   (同上)
   real:: ETT_sum     = 0.0 !TOTAL PLANT TRANSPIRATION (同上)
   real:: ESNOW_sum   = 0.0 !SUBLIMATION FROM SNOWPACK (同上)
   real:: SNOMLT_sum  = 0.0 !SNOW MELT (同上)
   real:: RUNOFF1_sum = 0.0
   real:: RUNOFF2_sum = 0.0
   real:: SHEAT_sum   = 0.0
   real:: S_sum       = 0.0
   real:: FUP_sum     = 0.0
   real:: ETA_sum     = 0.0
   
!Exchange soil properties when NOAH-module is activated
   INTEGER SOILTYP  !Soil type index 1-9 (10番に独自の値を入力した)
! SOIL TYPES   ZOBLER (1986)	  COSBY ET AL (1984) (quartz cont.(1))
!  1	    COARSE	      LOAMY SAND	 (0.82)
!  2	    MEDIUM	      SILTY CLAY LOAM	 (0.10)
!  3	    FINE	      LIGHT CLAY	 (0.25)
!  4	    COARSE-MEDIUM     SANDY LOAM	 (0.60)
!  5	    COARSE-FINE       SANDY CLAY	 (0.52)
!  6	    MEDIUM-FINE       CLAY LOAM 	 (0.35)
!  7	    COARSE-MED-FINE   SANDY CLAY LOAM	 (0.60)
!  8	    ORGANIC	      LOAM		 (0.40)
!  9	    GLACIAL LAND ICE  LOAMY SAND	 (NA using 0.82)
! 10        Observed value @ Sppaskaya-pad, Siberia -Inserted by H.Sato-
   
   If (Flag_land_physics) then
      SOILTYP = SoilClass
      Select Case (SOILTYP)
         Case(1); W_fi = 0.421; W_sat = 0.421; W_wilt = 0.029
         Case(2); W_fi = 0.464; W_sat = 0.464; W_wilt = 0.119
         Case(3); W_fi = 0.468; W_sat = 0.468; W_wilt = 0.139
         Case(4); W_fi = 0.434; W_sat = 0.434; W_wilt = 0.047
         Case(5); W_fi = 0.406; W_sat = 0.406; W_wilt = 0.100
         Case(6); W_fi = 0.465; W_sat = 0.465; W_wilt = 0.103
         Case(7); W_fi = 0.404; W_sat = 0.404; W_wilt = 0.069
         Case(8); W_fi = 0.439; W_sat = 0.439; W_wilt = 0.066
         Case(9); W_fi = 0.421; W_sat = 0.421; W_wilt = 0.029
         Case default; W_fi = 0.406; W_sat = 0.406; W_wilt = 0.100 !Same value for SOILTYP=5
      End Select
   End If
   
!_____________ NOAH-LSM Integration (Phase 3: Initialize state variables)
!If (observed_data is available) then
!   !Observed values in Sppaskaya-pad, eastern-Siberia
!   TBOT = 270.65 !Annual constant bottom boundary soil temperature (K)
!                 !MEMO by H.Sato: 270.65K=-2.5 C taken from Vasiliev @ Fedorov (2003)
!   T1    = 263.6909 !Initial Skin Temperature (K)
!else
   !Adjust to mean annual air-temperature
   TBOT = sum(tmp_air(:,:)) / real(YearMaxClimate) / real(Day_in_Year) + ABS_ZERO
   T1   = sum(tmp_air(:,:)) / real(YearMaxClimate) / real(Day_in_Year) + ABS_ZERO
!endif
   
   STC   = (/264.45, 269.35, 269.65, 270.05, 270.35, &
             270.75, 271.15, 271.55, 271.95, 272.25, &
             272.55, 272.75, 272.85, 272.95, 273.05, &
             273.05, 273.05, 273.05, 273.05, 273.05 /)
   
   SMC   = (/0.356, 0.278, 0.266, 0.259, 0.236, &
             0.211, 0.193, 0.185, 0.184, 0.188, &
             0.194, 0.202, 0.210, 0.220, 0.230, &
             0.239, 0.247, 0.253, 0.259, 0.269 /)
   
   SH2O  = (/0.146, 0.149, 0.153, 0.153, 0.152, &
             0.149, 0.149, 0.151, 0.157, 0.165, &
             0.175, 0.186, 0.200, 0.213, 0.224, &
             0.234, 0.242, 0.246, 0.245, 0.242 /)
   
   CMC   = 0.0000000E+00 !Initial canopy water content (m)
!  CMC   = 3.9353027E-04 !Initial canopy water content (m)
   
   SNOWH = 0.6301540     !Initial actual snow depth (m)
!  SNOWH = 1.0600531E-03 !Initial actual snow depth (m)
   
   SNEQV = 5.6783587E-02 !Initial water equiv snow depth (m)
!  SNEQV = 2.0956997E-04 !Initial water equiv snow depth (m)
   
   SOILM = 0.0
   do i = 1, NSOIL
      SOILM = SOILM + SLDPTH(i)*SMC(i) !TOTAL SOIL COLUMN WATER CONTENT in M
   end do
   
   LVH2O =  2.4501000E+6 !Vaporization heat for water at 20C (J/Kg H2O）
   
!_____________ Intialize variables
!Initialize variables
   Call init_value (W_fi, tmp_air(:,1), tmp_soil(:,1,:), prec(:,1))
   Call radiation_seasonal_change (LAT)
   
!Read spinup files
   Spinup_year  = 0
   if (Flag_spinup_read) then
      open (File_no(1), file=Fn_spnin, &
      access="sequential", form="unformatted", status="old", action="read")
      
      Call spinup_in (File_no(1), NSOLD, CMC, SNOWH, SNEQV, T1, STC, SMC, SH2O)
      close (File_no(1))
      
      Call direct_radiation (LAT)
      Call floor_radiation  ()
      Call crown_coverage   ()
      
      SOILM = 0.0
      do i = 1, NSOIL
         SOILM = SOILM + SLDPTH(i)*SMC(i) !TOTAL SOIL COLUMN WATER CONTENT in M
      end do
   endif
   
!_____________ Open output files
   open ( File_no( 1), file = 'output/log.txt'          )
   open ( File_no( 2), file = 'output/output.txt'       )
   open ( File_no( 3), file = 'output/climate.txt'      )
   open ( File_no( 4), file = 'output/air.txt'          )
   open ( File_no( 5), file = 'output/radiation.txt'    )
   open ( File_no( 6), file = 'output/water.txt'        )
   open ( File_no( 7), file = 'output/grass.txt'        )
   open ( File_no( 8), file = 'output/lai.txt'          )
   open ( File_no( 9), file = 'output/cflux.txt'        )
   open ( File_no(10), file = 'output/netradiation.txt' )
   open ( File_no(11), file = 'output/wflux.txt'        )
   open ( File_no(12), file = 'output/ld_vertical.txt'  )
   open ( File_no(13), file = 'output/annual.txt'       )
   open ( File_no(14), file = 'output/forest.txt'       )
   open ( File_no(15), file = 'output/biomass.txt'      )
   
!_____________ NOAH-LSM Integration (Phase 5: Compute daily chages of radiation)
IF (Flag_land_physics) THEN
   Call radiation_daily_change (LAT, TimeStepMax, RadShortWeight)
END IF
   
!_____________ Initialize various counters for the daily simulation loop
   doy          = 0
   year         = 1
   year_climate = 1
   
   if (Flag_spinup_read) then
      counter_begin = 1 + Spinup_year*Day_in_Year
      counter_end   = (Simulation_year+Spinup_year) * Day_in_Year
   else
      counter_begin = 1
      counter_end   = Simulation_year * Day_in_Year
   endif
   
!*****************************************************************************
!This loop corresponds to a simulation day, which calls subroutines sequentially.
!For readability of the code, I tried to avoid to call another subroutine 
!from the subroutine that was called from this loop.
DO counter = counter_begin, counter_end
   
!_____________ Daily update of field statuses
!Time counters1
   doy = doy + 1
   if (doy == Day_in_Year + 1) then
      doy          = 1
      year         = year + 1
      year_climate = year_climate + 1
      if (year_climate==YearMaxClimate+1) year_climate=1
   endif
   
   IIDAY  = IIDAY + 1          !Total Simulation day counter (Initial value = 0)
   IMONTH = Month (DOY)        !Month (1~12)
   IDAY   = Day_of_Month (DOY) !Day in the month (1~31)
   
!Time counters2 (wild fire and phenology related)
   dfl_fire = dfl_fire + 1 !Day from the last fire
   
   Do p=1, PFT_no
      dfl_leaf_onset(p) = dfl_leaf_onset(p) + 1 !Day from the last leaf onset
      dfl_leaf_shed (p) = dfl_leaf_shed (p) + 1 !Day from the last leaf shedding
   End do
   
!_____________ Prepare Climatic data for this cycle
!Daily mean meteological properties
   tmp_air_Today       = tmp_air       (doy, year_climate)
   tmp_air_range_Today = tmp_air_range (doy, year_climate)
   prec_Today          = prec          (doy, year_climate)
   rad_short_Today     = rad_short     (doy, year_climate)
   rad_long_Today      = rad_long      (doy, year_climate)
   wind_Today          = wind          (doy, year_climate)
   rh_Today            = rh            (doy, year_climate)
   
!_____________ NOAH-LSM Integration (Phase 6: Prepare climate data at 30min time-resoluition)
IF (Flag_land_physics) THEN
   call climate_convert ( &
        year_climate, YearMaxClimate, TimeStepMax, &
        RadShortWeight, &
        tmp_air, tmp_air_range, prec, rad_short, rad_long, rh, &
        AirTmp_30min, Prec_30min, RadShort_30min, RadLong_30min, rh_30min )
ELSE
   do i=1, NumSoil
    tmp_soil_Today(i)  = tmp_soil (doy, year_climate, i)
   end do
ENDIF
   
!_____________ Daily update of metabolic status
!stat_water(1:PFT_no), a vegetation growth limitter due to soil water shortage
   Do p = 1, PFT_no
      x  = 0.0
      no = 0
      
      do i=1, max(1,RootDepth(p))
      if (tmp_soil_Today(i)<=0.0) cycle
         no = no + 1
         x  = x + (pool_w(i)/Depth - W_wilt) / max(W_fi-W_wilt, 0.001)
      enddo
      x = x / max(1,no) ; x = max(min(x, 1.0), 0.0)
      stat_water(p) = x !Default
      
   End do
   
!Growth suppresss regulator for each tree
!(Suppressed trees do not conduct tree growth)
   if (doy==Day_in_Year) then
      flag_suppress(:) = 0
      do no=1, Max_no
         if ( tree_exist(no) ) then
         if ( mort_regu1(no)<10.0 .and. age(no)>3 ) then
            flag_suppress(no)=1
         endif
         endif
      enddo
   end if
   
!_____________ Daily update of field records (running recorders, annual means, etc..)
!Daily update of Climate statuses
!Calculate climatic statistics, their running record, and their running mean
   write(*,*) 'stat_climate' !!!<<<<<<<<<<<<TN:add
   Call stat_climate (prec_Today, tmp_air_Today, tmp_soil_Today)
   
!Daily update of Carbon pools and Carbon Fluxes statuses
!Reset variables, update their record in array variables, and calculate their running mean
   write(*,*) 'stat_carbon' !!!<<<<<<<<<<<<TN:add
   Call stat_carbon ()
   
!Daily update of vegetation statuses
!Reset variables, update their record in array variables, and calculate their running mean
   write(*,*) 'stat_vegetation' !!!<<<<<<<<<<<<TN:add
   Call stat_vegetation ()
   
   !mortality related
   if (doy==1) then
      mort_regu1(:) = 0.0 !NPP annual (g / individual)
      mort_regu2(:) = 0.0 !average leaf area of last year (m2/day) {update on: growth_wood}
      mort_regu4(:) = 0.0 !stem diameter increament in last year (m year-1)
   endif
   
!_____________ DAILY PSYSICAL PROCESSES
   !Time series of atmospheric CO2
   !co2atm = aco2_annual(min(year_climate,YearMaxCO2)) !Progress
    co2atm = aco2_annual(150)                          !Fixed at 2000
   !co2atm = aco2_annual(1)                            !Fixed at 1850
   
   !Calculate variables in relation to atmospheric physics (e.g. Air pressure, Vapor density)
   write(*,*) 'air' !!!<<<<<<<<<<<<TN:add
   Call air (tmp_air_Today, rh_Today, ALT)
   
   write(*,*) 'radiation' !!!<<<<<<<<<<<<TN:add
   Call radiation (LAT, rad_short_Today, cloud)
   par_RunningRecord(2:Day_in_Year) = par_RunningRecord(1:Day_in_Year-1)
   par_RunningRecord(1)             = par
   
   write(*,*) 'diffused_radiation' !!!<<<<<<<<<<<<TN:add
   Call diffused_radiation ()
   
   if (mod(doy, Days_LightComp_Int)==1) then
      !Compute relative intensity of direct radiation for each crown disk of each tree.
      write(*,*) 'direct_radiation' !!!<<<<<<<<<<<<TN:add
      Call direct_radiation (LAT)
      
      !Compute relative intensity of radiation for each cell of tree establishment.
      write(*,*) 'floor_radiation' !!!<<<<<<<<<<<<TN:add
      Call floor_radiation ()
      
      write(*,*) 'crown_coverag' !!!<<<<<<<<<<<<TN:add
      Call crown_coverage ()
      
   endif
   
   !Wild fire subroutines
   if     (GlobalZone==1) then
     !African continent
!     Call fire_regime2 (wind_Today) !!!>>>>>>>>>>>>>>>>>>>TN: rm
   else
     !Default
!     Call fire_regime (W_fi) !!!>>>>>>>>>>>>>>>>>>>TN: rm
   endif
   
!*****************************************************************************
!_____________ NOAH-LSM integration (Phase 7)
IF (Flag_land_physics) THEN
   write(*,*) 'NOAH-LSM integration' !!!<<<<<<<<<<<<TN:add
   x      = sum(lai_RunningRecord(1,:))        !LAI for all PFTs
   SHDFAC = min(1.00,max(0.01, LOG(x*4.5+1) )) !SHDFAC: Green vegetation fraction
   XLAI   = x / SHDFAC                         !XLAI: LAI of green vegetation area (m2/m2)
   
!   !Alternative method (Under testing)
!   !In the caluculation of 'Vegetatino Fraction', which is employed by NOAH-LSM, 
!   !'Vegetation' is determined for grass cells, of which relative PAR intensity is
!   !less than 95.12% (This corresponds to LAI>0.1).
!   !
!   !SHDFAC: Green vegetation fraction
!   !XLAI  : LAI of green vegetation area (m2/m2), LAI = SHDFAC * XLAI
!   n = 0
!   do i=1, Dived
!   do j=1, Dived
!       if (par_floor_rel(i,j)<=0.951229) n=n+1 !LAI=0.1以上に相当
!   enddo
!   enddo
!   if (n<=Dived) then
!      !in case of sparse foliages
!      SHDFAC = 0.0
!      XLAI   = 0.0
!   else
!      !in case of thick foliages
!      SHDFAC = real(n) / (real(Dived)*real(Dived))
!      XLAI   = lai / SHDFAC
!   endif
   
!Update land surface properties, which relate vegetation structures
!This section contains all procedures those provide variables from SEIB-DGVM to NOAH-LSM
   
   !Vegetation type index 1-13
      select case (biome)
      case ( 0); VEGTYP = 13! 13: GLACIAL (THE SAME PARAMETERS AS FOR TYPE 11)
      case ( 1); VEGTYP = 11! 11: BARE SOIL
      case ( 2); VEGTYP = 10! 10: DWARF TREES AND SHRUBS WITH GROUNDCOVER (TUNDRA)
      case ( 3); VEGTYP =  1!  1: BROADLEAF-EVERGREEN TREES  (TROPICAL FOREST)
      case ( 4); VEGTYP =  1!  1: BROADLEAF-EVERGREEN TREES  (TROPICAL FOREST)
      case ( 5); VEGTYP =  4!  4: NEEDLELEAF-EVERGREEN TREES
      case ( 6); VEGTYP =  1!  1: BROADLEAF-EVERGREEN TREES  (TROPICAL FOREST)
      case ( 7); VEGTYP =  2!  2: BROADLEAF-DECIDUOUS TREES
      case ( 8); VEGTYP =  4!  4: NEEDLELEAF-EVERGREEN TREES
      case ( 9); VEGTYP =  5!  5: NEEDLELEAF-DECIDUOUS TREES (LARCH)
      case (10); VEGTYP = 10! 10: DWARF TREES AND SHRUBS WITH GROUNDCOVER (TUNDRA)
      case (11); VEGTYP =  7!  7: GROUNDCOVER ONLY (PERENNIAL)
      case (12); VEGTYP = 11! 11: BARE SOIL
      case  default; VEGTYP = 11! 11: BARE SOIL
      end select
   
   !★NOAHのパラメーターをカラマツ林で固定したい場合には、以下の行をActivate
!     VEGTYP = 5       !@forest (5 -> NEEDLELEAF-DECIDUOUS TREES (LARCH))
   
   if (VEGTYP==13 .or. VEGTYP==11 .or. VEGTYP==10 .or. VEGTYP==7) then
      !Max albedo over very deep snow
      SNOALB = 0.69  !@after clear cut
      !Snow free Albedo for larch grassland (Hollinger et al. 2009)
      ALB    = 0.19  !@after clear cut
      !Root depth (number of soil layers)
      NROOT  = 1     !@after clear cut
   else
      !Max albedo over very deep snow
      SNOALB = 0.37    !@forest
      !Snow free Albedo for larch forests (Hollinger et al. 2009)
      ALB    = 0.14    !@forest
      !Root depth (number of soil layers)
      NROOT = 5        !@forest
   endif
   
   !For regulating heat conductance of the top-soil-layer by litter amount
      !x: Aboveground litter content (KgDM/m2)
      x =     pool_litter_leaf + pool_litter_ag + 0.5*pool_litter_trunk !(in gDM/stand)
!      x = x / Max_loc / Max_loc / 1000.0                                !(convert to Kg DM/m2)  !!!>>>>>>>>>>>>>>TN: rm
      x = x / real(GRID%Area) / 1000.0                                !(convert to Kg DM/m2)  !!!<<<<<<<<<<<<<<<TN: add
      
      DF1_given = 0.17 - 0.05 * (x/10.0)
      DF1_given = max( 0.01 , DF1_given )
   
!Each TimeStep Loop
DO TimeStep= 1, TimeStepMax
   !Total simulation time-step counter (Initial value = 0)
      INDI=INDI+1
   
   !Canopy Conductance for this time step
   !canopy_cond: Canopy Conductance (mol H2O m-2 s-1)
   !RC         : Canopy Resistence (s m-1)
   !Conver with the equation "40 mol m-2 s-1 = 1 m s-1", which is desctirbed in 
   !Jones (1992), Plant and Microclimate, 2nd edition, Equation 3.23, Cambridge press
   !   i = int( TimeStepMax * 0.5 ) !Time Step @ Mid day
   !   x = RadShort_30min(TimeStep) / max(0.1, RadShort_30min(i) ) !relative sun light intensity
   !   y = -0.14 * ( AirTmp_30min(TimeStep) ) + 43.927 !conversion coefficient 
   !   RC = y / max(canopy_cond * x, 0.01)
   
!Prepare physical environment
   !Give climate variables for each time step
   SFCTMP =   AirTmp_30min   (TimeStep) + ABS_ZERO !Air temperature (K)
   PRCP   =   Prec_30min     (TimeStep) / DT  !Precipitation amount (mm/sec)
   SOLDN  =   RadShort_30min (TimeStep)       !Shortwave radiation, Downward (W/m^2)
   LWDN   =   RadLong_30min  (TimeStep)       !Longwave radiation, Downward  (W/m^2)
   RHumid =   rh_30min       (TimeStep)       !Relative humidity    (%)
   SFCSPD =   wind_Today                      !Average wind speed   (m s-1)
   SFCPRS =   ap*100.0                        !Pressure at 1st mdl lvl abv skin (Pascals)
   
   !CALCULATE A SATURATION MIX RATIO (Q2SAT)
   !NEED Q2 (FROM REL.HUMID.) USE SUBROUTINE QDATAP
   !Q2 = specific humidity (g/kg)
   CALL QDATAP (SFCTMP, SFCPRS, RHumid, Q2, Q2SAT, ESAT, LVH2O)
   IF (Q2 < 0.1E-5 ) Q2 = 0.1E-5
   IF (Q2 >=  Q2SAT) Q2 = Q2SAT*0.99
   
   !CALCULATE SLOPE OF SAT SPECIFIC HUMIDITY CURVE FOR PENMAN: DQSDT2
   DQSDT2 = DQSDT (SFCTMP, SFCPRS, LVH2O)
   
   !CALC VIRTUAL TEMPS AND POTENTIAL TEMPS AT GRND (SUB 1) AND AT
   !THE 1ST MDL LVL ABV THE GRND (SUB 2). EXPON IS CP DIVD BY R.
   TH2 = SFCTMP + ( 0.0098 * Z ) !AIR POTENTIAL TEMPERATURE (K) AT HEIGHT ZLVL ABOVE GROUND
   T2V = SFCTMP * (1.0 + 0.61 * Q2 )
   
   T1V  =    T1 * (1.0 + 0.61 * Q2 )
   TH2V =   TH2 * (1.0 + 0.61 * Q2 )
   
   ! Determine precipitation type
   !When snowing FFROZP=1.0 ; Otherwise FFROZP=0.0
   FFROZP = 0.0
   if (PRCP > 0.0)         then !PRCP = precipitation 
   if (SFCTMP <= ABS_ZERO) then !Air temperature (K) <= 273.15
      FFROZP = 1.0
   endif
   endif
   
   !SNOW COVER, ALBEDO OVER LAND
   ! This equation is taken from P29 of Oleson, K. W., et al. (2004)
   ! Technical description of the Community Land Model (CLM),
   ! Tech. Note NCARTN‐461+STR, 174 pp., Natl. Cent. for Atmos. Res., Boulder, Colo.
   IF (SNEQV == 0.0) THEN
      SNCOVR = 0.0
      ALBEDO = ALB
   ELSE
      !10.0:a factor constant , 0.01:momentam roughness length (m)
      SNCOVR = SNOWH / (SNOWH+10.0*0.01)
      ALBEDO = ALB + SNCOVR*(SNOALB-ALB)
   ENDIF
   
   !SNOW COVER, ALBEDO OVER SEA-ICE
   IF (ICE .EQ. 1) THEN
      SNCOVR = 1.0
      ALBEDO = 0.65
   ENDIF
   
   SOLNET = SOLDN*(1.0-ALBEDO)
   
!CALL LAND-SURFACE PHYSICS
   CALL SFLX ( &
   FFROZP,ICE,DT,Z,NSOIL,SLDPTH, &
   LWDN,SOLDN,SOLNET,SFCPRS,PRCP,SFCTMP,Q2,SFCSPD, &
   TH2,Q2SAT,DQSDT2, &
   VEGTYP,SOILTYP,SLOPETYP,SHDFAC,SHDMIN,PTU,ALB,SNOALB,TBOT, &
   CMC,T1,STC,SMC,SH2O,SNOWH,SNEQV,ALBEDO,CH,CM, &
   ETA,SHEAT, &
   EC,EDIR,ET,ETT,ESNOW,DRIP,DEW, &
   BETA,ETP,S, &
   FLX1,FLX2,FLX3, &
   SNOMLT,SNCOVR,SNUPX,SALP, &
   RUNOFF1,RUNOFF2,RUNOFF3, &
   RC,PC,RSMIN,XLAI,RCS,RCT,RCQ,RCSOIL, &
   SOILW,SOILM,SMCWLT,SMCDRY,SMCREF,SMCMAX,NROOT,LVH2O, &
   DF1_given) 
   
   !AET: ACTUAL EVAPOTRANSPIRATIVE ENERGY
   !ETA: ACTUAL LATENT HEAT FLUX
   AET = ETA
   
   !CALCULATE UPWARD LONGWAVE RAD USING UPDATED SKIN TEMPERATURE
   FUP = 5.67E-8 * T1 * T1 * T1 * T1
   
   ! CALCULATE RESIDUAL OF ALL SURFACE ENERGY BALANCE EQN TERMS.
   S      = -S
   x      = SOLDN*(1.0-ALBEDO) + LWDN
   RES    = x - SHEAT - S - AET - FUP - FLX1 - FLX2 - FLX3
   
!出力関係の処理
   !1日のシミュレーション分割数
   x = real(TimeStepMax)
   
   !Sumup daily values in mm for water balance
   !Following valiables are initialized at the begginig of each day
   ETP_sum     = ETP_sum    + ETP   *DT/LVH2O  !Potential Evapotranspiration (mm)
   EC_sum      = EC_sum     + EC    *DT/LVH2O  !CANOPY WATER EVAPORATION (mm)
   EDIR_sum    = EDIR_sum   + EDIR  *DT/LVH2O  !DIRECT SOIL EVAPORATION (mm)
   ETT_sum     = ETT_sum    + ETT   *DT/LVH2O  !TOTAL PLANT TRANSPIRATION (mm)
   ESNOW_sum   = ESNOW_sum  + ESNOW *DT/LVH2O  !SUBLIMATION FROM SNOWPACKM (mm)
   SNOMLT_sum  = SNOMLT_sum + SNOMLT*DT/LVH2O  !SNOW MELT (mm)
   RUNOFF1_sum = RUNOFF1_sum + RUNOFF1 * DT * 1000. !Runoff1 (mm)
   RUNOFF2_sum = RUNOFF2_sum + RUNOFF2 * DT * 1000. !Runoff2 (mm)
   
   SHEAT_sum  = SHEAT_sum  + SHEAT  /x !Sensible heat flux (W M-2)
   S_sum      = S_sum      + S      /x !Grnd SFC flux (W M-2)
   FUP_sum    = FUP_sum    + FUP    /x !Upward grnd LW radiation (W M-2)
   ETA_sum    = ETA_sum    + ETA    /x !Actual Latent heat flux (W M-2)
   
!END OF DAILY TIME LOOP
END DO !TimeStep loop
   
!This paragraph contains all procedures for providing variables from NOAH-LSM to SEIB-DGVM
   do i=1, NumSoil
      pool_w        (i) = SMC(i)*Depth     
      tmp_soil_Today(i) = STC(i) - ABS_ZERO
   end do
   
   x = SNEQV * 1000. - pool_snow !changes of snow pack
   pool_snow = SNEQV * 1000. !water equivalent snow depth (mm)
   
   flux_ro1_RunningRecord(2:Day_in_Year) = flux_ro1_RunningRecord(1:Day_in_Year-1)
   flux_ro2_RunningRecord(2:Day_in_Year) = flux_ro2_RunningRecord(1:Day_in_Year-1)
   flux_ic_RunningRecord (2:Day_in_Year) = flux_ic_RunningRecord (1:Day_in_Year-1)
   flux_ev_RunningRecord (2:Day_in_Year) = flux_ev_RunningRecord (1:Day_in_Year-1)
   flux_tr_RunningRecord (2:Day_in_Year) = flux_tr_RunningRecord (1:Day_in_Year-1)
   flux_sl_RunningRecord (2:Day_in_Year) = flux_sl_RunningRecord (1:Day_in_Year-1)
   flux_tw_RunningRecord (2:Day_in_Year) = flux_tw_RunningRecord (1:Day_in_Year-1)
   flux_sn_RunningRecord (2:Day_in_Year) = flux_sn_RunningRecord (1:Day_in_Year-1)
   ev_pot_RunningRecord  (2:Day_in_Year) = ev_pot_RunningRecord  (1:Day_in_Year-1)
   pool_w1_RunningRecord (2:Day_in_Year) = pool_w1_RunningRecord (1:Day_in_Year-1)
   
   flux_ro1_RunningRecord(1) = RUNOFF1_sum
   flux_ro2_RunningRecord(1) = RUNOFF2_sum
   flux_ic_RunningRecord (1) = EC_sum
   flux_ev_RunningRecord (1) = EDIR_sum
   flux_tr_RunningRecord (1) = ETT_sum
   flux_sl_RunningRecord (1) = ESNOW_sum
   flux_tw_RunningRecord (1) = SNOMLT_sum
   flux_sn_RunningRecord (1) = max(0.0, x + ESNOW_sum + SNOMLT_sum)
   ev_pot_RunningRecord  (1) = ETP_sum
   pool_w1_RunningRecord (1) = sum(pool_w(1:5)) / 5.0
   
ENDIF !Endo oof the NOAH-LSM process

!*****************************************************************************
!_____________ DAILY BIOLOGICAL PROCESSES
   !Calculate photosynthesis rate controlling variables
   Call photosynthesis_condition (tmp_air_Today)
   
if (Flag_photosynthesis_type) then
   !Calculate photosynthesis rate controlling variables for the Faquhar's equation
   Call photosynthesis_condition_farquhar  (tmp_air_Today, DT)
   
   !Photosynthesis process Using the the Faquhar's equation
   do TimeStep= 1, TimeStepMax
   Call photosynthesis_timestep (TimeStep, TimeStepMax, RadShort_30min)
   enddo
   
else
   !Photosynthesis process
   Call photosynthesis ()
end if
   
   !Calculate optimal leaf area index for grass layer
   write(*,*) 'lai_optimum' !!!<<<<<<<<<<<<TN:add
   Call lai_optimum (tmp_air_Today, sum(tmp_soil_Today(1:5))/5.0, TimeStepMax, RadShort_30min  )
   
   !Maintenance respiration
   write(*,*) 'maintenance_resp' !!!<<<<<<<<<<<<TN:add
   Call maintenance_resp (tmp_air_Today, sum(tmp_soil_Today(1:5))/5.0 )
   
   !Turnover of plant organs
   write(*,*) 'turnover' !!!<<<<<<<<<<<<TN:add
   Call turnover ()
   
   !Penology controller, and biological processes before and after penology change
   !(such as gradual release of stock biomass after onset of foliage phase).
   write(*,*) 'leaf_season' !!!<<<<<<<<<<<<TN:add
   Call leaf_season (LAT, tmp_soil_Today )
   
   !Daily growth procedure for woody PFTs (foliation, fine root growth, and reproduction)
   write(*,*) 'growth_wood' !!!<<<<<<<<<<<<TN:add
   Call growth_wood ()
   
   !Growth and reproduction procedure for grass PFTs.
   write(*,*) 'growth_grass' !!!<<<<<<<<<<<<TN:add
   Call growth_grass ()
   
   !Calculate free space around each individual tree. 
   !This information will be employed by following subroutine growth_trunc, 
   !where stem and crown expansion occurs.
   write(*,*) 'spatial_limitation' !!!<<<<<<<<<<<<TN:add
   Call spatial_limitation ()
   
   !Expands stem diameter, which cause stem height and crown diameter increases
   write(*,*) 'growth_trunk' !!!<<<<<<<<<<<<TN:add
   Call growth_trunk ()
   
   !Heterotrophic respiration (Litter and soil organic matter decomposition)
   write(*,*) 'decomposition' !!!<<<<<<<<<<<<TN:add
   Call decomposition (W_fi)
   
!_____________ ANNUALY BIOLOGICAL PROCESSES
   
IF ( doy==1 ) then
   !Listing all PFT that exist in this virtual forest
   Call pft_present ()
ENDIF
   
IF ( doy==Day_in_Year ) then
   
   !Determine which tree die
   write(*,*) 'mortality' !!!<<<<<<<<<<<<TN:add
   Call mortality ()
   
   !Adjust crown depth by purging crown disks from the bottom of the crown layer
   write(*,*) 'crown_adjust' !!!<<<<<<<<<<<<TN:add
   Call crown_adjust ()
   
   !Cause horizontal movement of crown of each tree to the open direction
   write(*,*) 'crown_shake' !!!<<<<<<<<<<<<TN:add
   Call crown_shake ()
   
   !Establishment of woody PFTs
   x = sum(prec_RunningRecord(1:Day_in_Year))                  !annual precipitation (mm/year)
   y = sum(tmp_air_RunningRecord(1:Day_in_Year)) / Day_in_Year !annual mean temperature (Celcius)
   if ( x > max(100.0, 20.0*y) ) then                          !from Koppen's criteria (1936)
      
      !Determine ground mesh point, where newly saplings can establish.
      !This subroutine should be called before subroutine 'establish'.
      write(*,*) 'ground_vacant' !!!<<<<<<<<<<<<TN:add
      Call ground_vacant ()
      
      !Recruitment of trees
      write(*,*) 'establish' !!!<<<<<<<<<<<<TN:add
      Call establish (GlobalZone)
      
   end if
   
   !Compute fraction of tree crown coverage.
   write(*,*) 'crown_coverage' !!!<<<<<<<<<<<<TN:add
   Call crown_coverage ()
   
   !Determine dominant grass PFT in the next year
   write(*,*) 'grass_priority' !!!<<<<<<<<<<<<TN:add
   Call grass_priority ()
   
   !Determine current biome type of this grid cell
   write(*,*) 'biome_determine' !!!<<<<<<<<<<<<TN:add
   Call biome_determine ()
   
END IF
   
!_____________ ANNUALY Physics PROCESSES
!Comment out to disable automatic adjustment of soil bottom temperature
IF ( doy==Day_in_Year ) then
   x = sum (tmp_soil_RunningRecord(1:Day_in_Year, NumSoil)) / real(Day_in_Year)
   x = x + ABS_ZERO
   TBOT = 0.9*TBOT + 0.1*x
END IF
   
!_____________ DAILY UPDATE NET-RADIATION & SOIL-WATER STATUS
   
   !Compute Albedo
   write(*,*) 'albedo_calc' !!!<<<<<<<<<<<<TN:add
   Call albedo_calc (Albedo_soil0)
   
   !Compute net radiation
   write(*,*) 'net_radiation' !!!<<<<<<<<<<<<TN:add
   Call net_radiation (tmp_air_Today, rad_short_Today, rad_long_Today, cloud)
   
   !Compute water cycle
   if (.not. Flag_land_physics) then
   write(*,*) 'waterbudget' !!!<<<<<<<<<<<<TN:add
   Call waterbudget &
   (W_fi, W_wilt, prec_Today, wind_Today, tmp_air_Today, tmp_soil_Today(:) )
   endif
   
!_____________ DAILY INCREMENT of LITTER AND FUEL POOLS
   ! Increment of litters pools and fuel pools are basically conducted here
   ! by considering litter fluxes, which are computed in descendent subroutines.
   ! This way of convergent-computation is for handling 
   ! a variety of litter pools and fuel pools in the model.
   
   flux_c_lit_RR(1) =  flux_litter_trunk + flux_litter_leaf + flux_litter_root + &
                       flux_litter_ag + flux_litter_bg
   
   pool_litter_leaf  = max(0.0, pool_litter_leaf  + flux_litter_leaf )
   pool_litter_trunk = max(0.0, pool_litter_trunk + flux_litter_trunk)
   pool_litter_root  = max(0.0, pool_litter_root  + flux_litter_root )
   pool_litter_ag    = max(0.0, pool_litter_ag    + flux_litter_ag   )
   pool_litter_bg    = max(0.0, pool_litter_bg    + flux_litter_bg   )
   
   pool_fuel_standT  = max(0.0, pool_fuel_standT  + flux_litter_leaf )
   pool_fuel_standG  = max(0.0, pool_fuel_standG  + flux_litter_ag   )
   
   !Reset sumup-variables for daily ourputs
   ETP_sum     = 0.0
   EC_sum      = 0.0 
   EDIR_sum    = 0.0 
   ETT_sum     = 0.0 
   ESNOW_sum   = 0.0 
   SNOMLT_sum  = 0.0 
   RUNOFF1_sum = 0.0
   RUNOFF2_sum = 0.0
   SHEAT_sum   = 0.0
   S_sum       = 0.0
   FUP_sum     = 0.0
   ETA_sum     = 0.0
   
!_____________ MAKE OUTPUT FILES (Write simulation results in output files)
IF (Flag_output_write) then
   write(*,*) 'Write simulation results' !!!<<<<<<<<<<<<TN:add
   
   !Daily output (everyday)
   Call output_for_viewer &
   (File_no(2), LAT, LON, ALT, Albedo_soil0, W_fi, W_wilt, W_sat, &
   tmp_air_Today, &
   sum(tmp_soil_Today(1:5))/5.0, sum(tmp_soil_Today(6:10))/5.0, sum(tmp_soil_Today(11:20))/10.0, &
   prec_Today, rh_Today, cloud, rad_short_Today, rad_long_Today, wind_Today)
   
   !Daily output (on each specifyed year)
   Call output_climate &
   (File_no(3), cloud, prec_Today, rh_Today, wind_Today, tmp_air_Today, &
   sum(tmp_soil_Today(1:5))/5.0, rad_short_Today, rad_long_Today)
   Call output_air          (File_no(4))
   Call output_radiation    (File_no(5), rad_short_Today)
   Call output_water        (File_no(6), W_fi, W_wilt, W_sat)
   Call output_grass        (File_no(7))
   Call output_lai          (File_no(8))
   Call output_cflux        (File_no(9))
   
   !When the NOAH-LSM treats land-physics, 
   !variables concerining hydrology and net-radiation in the SEIB-DGVM are meaningless
   if ( .not. Flag_land_physics ) then
      Call output_netradiation (File_no(10))
      Call output_wflux        (File_no(11), prec_Today)
   endif
   
   !Monthly output
   if ( Day_of_Month(doy) == Day_in_month(Month(doy)) ) then
      Call output_ld_vertical (File_no(12))
   endif
   
   !Annual output (At the end of the year)
   if (doy == Day_in_Year) then
!  if (doy == 213) then !@ 1 Aug
      Call output_annual (File_no(13))
      Call output_forest (File_no(14))
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
     Call system_clock(time1)
     write(*,*) time1
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   endif
   
   !Annual output (At the middle of gwoeing season in Northern Hemisphere)
   if (doy == 195) then
      Call output_biomass(File_no(15))
   endif
   
ENDIF
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
 write(*,*) 'Year: ', year, 'Day: ', doy
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add

END DO !The end of daily loop
!*****************************************************************************

!_____________ After Main-simulation-loop procedures
   close (File_no( 1))   !output_log.txt
   close (File_no( 2))   !output.txt
   close (File_no( 3))   !output_climate.txt
   close (File_no( 4))   !output_air.txt
   close (File_no( 5))   !output_radiation.txt
   close (File_no( 6))   !output_water.txt
   close (File_no( 7))   !output_grass.txt
   close (File_no( 8))   !output_lai.txt
   close (File_no( 9))   !output_cflux.txt
   close (File_no(10))   !output_netradiation.txt
   close (File_no(11))   !output_wflux.txt
   close (File_no(12))   !output_ld_vertical.txt
   close (File_no(13))   !output_annual.txt
   close (File_no(14))   !output_forest.txt
   close (File_no(15))   !output_biomass.txt
   
IF ( Flag_spinup_write ) then
   !Write spinup data for restarting simualtion
!  open(File_no(1), file = Fn_spnout)
   open(File_no(1), file = Fn_spnout, &
        access="sequential", form="unformatted", status="replace", action="write")
   Call spinup_out (File_no(1), NSOLD, CMC, SNOWH, SNEQV, T1, STC, SMC, SH2O)
   close (File_no(1))
END IF

END SUBROUTINE main_loop



!**************************************************************************************************
! Compute seasonal changes in radiation related variables
!**************************************************************************************************
SUBROUTINE radiation_daily_change (LAT, TimeStepMax, RadShortWeight)
   USE data_structure
   USE grid_status_current2
   implicit none
   
!_________________ Define augments
   !LATITUDE..(N > 0.00 (+); S < 0.00 (-))
   real,intent(IN)::  LAT
   
   !Number of Time step in a day, 60*60*24/int(DT)
   Integer,intent(IN):: TimeStepMax
   
   !Relative intensity of the short-wave radiation for each time-step of the day (0.0~1.0)
   !daily maximum value (midday value) == 1.0
   real,dimension(1:Day_in_Year, 1:TimeStepMax)::RadShortWeight
   
!_________________ Define Local Variables
   integer TimeStep !time counter
   integer doy      !time counter
   real    x, y, a1, a2 !for general usage
   
!_________________ RadShortWeight(day,TimeStep)
   !Initialize
   RadShortWeight (:,:) = 0.0
   
Do doy  =    1,  365
   if (dlen(doy)<1.0) cycle
   
   !Sun break time @ TimeStepMax
   a1 = 0.5*( 24.0-dlen(doy) ) * real(TimeStepMax)/24.0
   
   !Sunset time @ TimeStepMax
   a2 = 0.5*( 24.0+dlen(doy) ) * real(TimeStepMax)/24.0
   
   !RadShortWeight: Weighting factor for daily changes in shortwave radiation（range：0.0～1.0）
   Do TimeStep = 1, TimeStepMax
      x = real(TimeStep)
      If (x < a1) cycle !From midnight until sun break
      If (x > a2) cycle !From sunset to midnight
      
      y = 1.0 - (a2-x)/(a2-a1)  !Linerly change from 0.0 (sun break) to 1.0 (sunset)
      y = sin(y*PI) * sin(y*PI) !Daily change in shortwave radiation
      RadShortWeight(doy,TimeStep) = y
   End Do
   
End Do

RETURN
END SUBROUTINE radiation_daily_change



!*****************************************************************************
SUBROUTINE climate_convert ( &
           year_climate, YearMaxClimate, TimeStepMax, &
           RadShortWeight, &
           tmp_air, tmp_air_range, prec, rad_short, rad_long, rh, &
           AirTmp_30min, Prec_30min, RadShort_30min, RadLong_30min, rh_30min )
   USE data_structure
   USE time_counter
   USE grid_status_current2
   
   implicit none
   
!_________________ Define Augments
   !Time counters
   Integer,intent(IN):: year_climate     !Climate year
   Integer,intent(IN):: YearMaxClimate !Climate year
   Integer,intent(IN):: TimeStepMax      !Number of Time step in a day, 60*60*24/int(DT)
   
   !Relative intensity of the short-wave radiation for each time-step of the day (0.0~1.0)
   !daily maximum value (midday value) == 1.0
   real,dimension(1:365, 1:TimeStepMax),intent(IN)::RadShortWeight
   
   !Input daily climate data
   real,dimension(365,YearMaxClimate),intent(IN):: &
    tmp_air      , & !Surface air temperature (C)
    tmp_air_range, & !Daily range of tmp_air  (C)
    prec         , & !Precipitation (mm day-1)
    rad_short    , & !Midday Shortwave rad. (W/m2)
    rad_long     , & !Daily mean longwave rad.(W/m2)
    rh               !Relative humidity (kg kg-1)
   
   !Output climate data (30min interval)
   real,dimension(1:TimeStepMax),intent(OUT):: &
    AirTmp_30min   , & !Air temperature (C)
    Prec_30min     , & !Precipitation amount (mm)
    RadShort_30min , & !Short-wave-radiation, Downward (W/m^2)
    RadLong_30min  , & !Long-wave-radiation, Downward (W/m^2)
    rh_30min        !Relative humidity (%)
   
!_________________ Define local variables
   !Variables for estimating daily changes in climatic conditions
   integer yearP, yearN               !Years for refering (P:Previous, N:Next)
   real    time1, time2, time3        !Duration of time
   real    temp1, temp2, temp3, temp4 !Temperature at various timing]
   real vp_          !vapour pressure (hPa)
   real vp_sat_      !saturated vapour pressure (hPa)
   
   !For general usage
   real    x, y
   integer TimeStep 
   
!_________________ Define 'previous' and 'next' years for interpolating air temperature
   yearP = year_climate - 1 !Previous Year
   yearN = year_climate + 1 !Next Year
   if (year_climate==1)              yearP = YearMaxClimate
   if (year_climate==YearMaxClimate) yearN = 1
   
!_________________ Initialization
   AirTmp_30min   (:) = 0.0
   Prec_30min     (:) = 0.0
   RadShort_30min (:) = 0.0
   RadLong_30min  (:) = 0.0
   rh_30min    (:) = 0.0
   
!_________________ (1) Air Temperature
   !Durations of time (0.0~real(TimeStepMax))
   time1 = 0.5 * (24.0-dlen(doy)) * real(TimeStepMax)/24.0 !From 0AM To sunbreak
   time2 = (dlen(doy) * 2.0/3.0)  * real(TimeStepMax)/24.0 !From sunbreak To max air-tmp of the day
   time3 = real(TimeStepMax) - time1 - time2               !From max air-tmp of the day To 12PM
   
   !Air-temperatures for reference
   ! temp1: Maximum air-temperature of previous day
   ! temp2: Minimum air-temperature of today
   ! temp3: Maximum air-temperature of today
   ! temp4: Minimum air-temperature of tomorrow
   temp2 = tmp_air (doy,year_climate) - 0.5 * tmp_air_range (doy,year_climate) 
   temp3 = tmp_air (doy,year_climate) + 0.5 * tmp_air_range (doy,year_climate) 
   if (doy==1) then
      temp1 = tmp_air (365,   yearP       ) + 0.5 * tmp_air_range (365,          yearP)
      temp4 = tmp_air (doy+1, year_climate) - 0.5 * tmp_air_range (doy+1, year_climate)
   elseif (doy==365) then
      temp1 = tmp_air (doy-1, year_climate) + 0.5 * tmp_air_range (doy-1, year_climate)
      temp4 = tmp_air (1,     yearN        ) - 0.5 * tmp_air_range (1,     yearN        )
   else
      temp1 = tmp_air (doy-1, year_climate) + 0.5 * tmp_air_range (doy-1, year_climate)
      temp4 = tmp_air (doy+1, year_climate) - 0.5 * tmp_air_range (doy+1, year_climate)
   endif
   
   !Interpolate daily changes in Air-temperature
   Do TimeStep = 1, TimeStepMax
      x = real(TimeStep)
      If     ( x < time1 ) then
      !From 0AM To sunbreak
         y = (time3+x) / (time3+time1)
         y = sin( PI * 0.5 * y )
         AirTmp_30min (TimeStep) = temp1 + (temp2-temp1)*y
      Elseif ( x < (time1+time2) ) then
      !From sunbreak To max air-tmp of the day
         y = (x-time1) / time2
         y = sin( PI * 0.5 * y )
         AirTmp_30min (TimeStep) = temp2 + (temp3-temp2)*y
      Else
      !From max air-tmp of the day To 12PM
         y = (x-time1-time2) / (time3+time1)
         y = sin( PI * 0.5 * y )
         AirTmp_30min (TimeStep) = temp3 + (temp4-temp3)*y
      Endif
   End Do
   
!_________________ (2) Precipitation (mm/daystep)
   x = prec(doy, year_climate) / 2.0 / 2.0
   Prec_30min( 1) = x
   Prec_30min( 2) = x
   Prec_30min(25) = x
   Prec_30min(26) = x
  
!Prec_30min(:) = Prec_30min(:) * 2.0

!_________________ (3) Shortwave Radiation
   Do TimeStep = 1, TimeStepMax
      RadShort_30min(TimeStep) = &
      rad_short(doy,year_climate) * RadShortWeight(doy,TimeStep)
   End Do
   
!_________________ (4) Longwave Radiation
   RadLong_30min(:) = rad_long(doy,year_climate)
   
!_________________ (5) Relative Humidity (%)
   Do TimeStep = 1, TimeStepMax
      
      !vp_sat_: saturated vapour pressure for each time step (hPa)
      vp_sat_ = Sat_vp (AirTmp_30min(TimeStep))
      
      if (vp_sat_<0.01) then
         RH_30min (TimeStep) = 0.0
      else
         RH_30min (TimeStep) = 100.*(vp/vp_sat_)
         RH_30min (TimeStep) = max(min( RH_30min(TimeStep), 100. ), 0.0 )
      endif
      
   End Do
   
RETURN
END SUBROUTINE climate_convert



! ★ここから下のサブルーチンは、コード整理を行っていない
!*****************************************************************************
FUNCTION DQS (T,LVH2O) 
   IMPLICIT NONE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C  PURPOSE:  TO CALCULATE VALUES OF VAPOR PRESSURE (E)
!C            AND P * DQS/DT (P TIMES CHG IN SAT MXG RATIO WITH RESPECT
!C            TO THE CHG IN TEMP) IN SUBSTITUTION TO THE LOOK-UP TABLES.
!C
!C            SUBSTITUTES LOOK-UP TABLES ASSOCIATED WITH THE DATA
!C            BLOCK  /CHMXR/ .
!C
!C            FORMULAS AND CONSTANTS FROM ROGERS AND YAU, 1989.
!C
!C                         ADDED BY PABLO J. GRUNMANN, 6/30/97.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      REAL DQS, T, LVH2O
      REAL DESDT, LW, ES, CPV, RV, CW, EPS, ESO, TO
      
!K    PARAMETER  (CP = 1005.)
!K    PARAMETER  (CV = 718.)
!K    PARAMETER (CVV = 1410.)
      PARAMETER (CPV = 1870.)
      PARAMETER  (RV = 461.5)
      PARAMETER ( CW = 4187.)
      PARAMETER (EPS =  0.622)
      PARAMETER (ESO =  611.2)
      PARAMETER ( TO =  273.15)
! 
!     ABOUT THE PARAMETERS:
!      
!     EPS ---------- WATER - DRY AIR MOLECULAR MASS RATIO, EPSILON
!      
!   VALUES FOR SPECIFIC HEAT CAPACITY AND INDIVIDUAL GAS CONSTANTS 
!   IN [JOULES/(KG*KELVIN)] UNITS.
!
!     DRY AIR: 
!             CP, CV
!     WATER VAPOR:
!                 CVV = 1410. 
!                 CPV = 1870.
!                 RV  =  461.5
!     LIQUID WATER:
!                  CW = 4187.
!
!     ESO = ES(T=273.15 K) = SAT. VAPOR PRESSURE (IN PASCAL) AT T=TO
!      TO = 273.15
!      
!     SAT. MIXING  RATIO: QS ~= EPS*ES/P
!     CLAUSIUS-CLAPEYRON: DES/DT = L*ES/(RV*T^2)
!     @QS/@T =  (EPS/P)*DES/DT
!    
   LW = LVH2O - ( CW - CPV ) * ( T - TO )
   ES = ESO*EXP (LW*(1/TO - 1/T)/RV)  
   DESDT = LW*ES/(RV*T*T)
!
!      FOR INSERTION IN DQSDT FUNCTION: 
!      DQSDT = DQS/P , WHERE DQS = EPS*DESDT  
!
   DQS = EPS*DESDT
   
RETURN
END

! ______________________________________________________________________
FUNCTION DQSDT ( SFCTMP, SFCPRS, LVH2O )
   IMPLICIT NONE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C    PURPOSE:  TO RETRIEVE THE APPROPRIATE VALUE OF DQSDT (THE CHANGE
!C    =======   OF THE SATURATION MIXING RATIO WITH RESPECT TO THE 
!C              CHANGE IN TEMPERATURE) FROM: 
!C
!C               FORMULAS INTRODUCED IN NEW FUNCTION DQS 
!C                                  (MODIFIED BY PABLO GRUNMANN, 7/9/97).
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL SFCTMP, SFCPRS, LVH2O, DQS, DQSDT

      IF ((SFCTMP >= 173.0) .AND. (SFCTMP  <=  373.0)) THEN

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!       IF THE INPUT SFC AIR TEMP IS BTWN 173 K AND 373 K, USE
!       FUNCTION DQS TO DETERMINE THE SLOPE OF SAT.MIX RATIO FUNCTION
!                                 -ADAPTED TO USE NEW DQS, 7/9/97.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        DQSDT = DQS (SFCTMP,LVH2O) / SFCPRS

      ELSE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!       OTHERWISE, SET DQSDT EQUAL TO ZERO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        DQSDT = 0.0

      END IF

RETURN
END

! ______________________________________________________________________
FUNCTION E(T,LVH2O) 
   IMPLICIT NONE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C  PURPOSE:  TO CALCULATE VALUES OF SAT. VAPOR PRESSURE (E)
!C            SUBSTITUTES LOOK-UP TABLES ASSOCIATED WITH THE DATA
!C            BLOCK  /VAPPRS/ .
!C            FORMULAS AND CONSTANTS FROM ROGERS AND YAU, 1989.
!C
!C                         ADDED BY PABLO J. GRUNMANN, 7/9/97.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
   REAL E, T, LVH2O
   REAL LW, CPV, RV, CW, ESO, TO
   
   PARAMETER (CPV = 1870.)
   PARAMETER  (RV = 461.5)
   PARAMETER  (CW = 4187.)
   PARAMETER (ESO = 611.2)
   PARAMETER  (TO = 273.15)      
! 
!     ABOUT THE PARAMETERS:
!      
!     EPS ---------- WATER - DRY AIR MOLECULAR MASS RATIO, EPSILON
!      
!   VALUES FOR SPECIFIC HEAT CAPACITY AND INDIVIDUAL GAS CONSTANTS 
!   IN [JOULES/(KG*KELVIN)] UNITS.
!
!     DRY AIR: 
!             CP, CV
!     WATER VAPOR:
!                 CVV = 1410. 
!                 CPV = 1870.
!                 RV  =  461.5
!     LIQUID WATER:
!                  CW = 4187.
!
!     ESO = ES(TO) = SAT. VAPOR PRESSURE (IN PASCAL) AT T=TO. [TO = 273.15]
!_______________________________________________________________________

!CLAUSIUS-CLAPEYRON: DES/DT = L*ES/(RV*T^2)
    LW = LVH2O - ( CW - CPV ) * ( T - TO )
    E = ESO * EXP (LW*(1/TO - 1/T)/RV)  
    
RETURN
END 

! ______________________________________________________________________
SUBROUTINE QDATAP (T,P,RH,QD,QS,ES,LVH2O) 
IMPLICIT NONE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C  PURPOSE:  OBTAIN SPECIFIC HUMIDITY (q) FROM RELATIVE HUMIDITY 
!C            AND GIVEN PRESSURE AND TEMPERATURE.
!C            
!C
!C            FORMULAS AND CONSTANTS FROM ROGERS AND YAU, 1989: 'A 
!C            SHORT COURSE IN CLOUD PHYSICS', PERGAMON PRESS, 3rd ED.
!C
!C                                   Pablo J. Grunmann, 3/6/98.
!C                Updated to eliminate subroutine SVP, 6/24/98.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!----------------------------------------
! In:
!        T    Temperature (K)
!        P    Pressure (Pa)
!        RH   Relative humidity (%)
!-----------------------------------------
! Out:
!        QD   Specific humidity (Kg/Kg)
!        QS   Saturation Specific humidity (Kg/Kg)
!        ES   Saturation vapor pressure for water (Pa)
!----------------------------------------
   REAL T, P, RH, RHF
   REAL QD, QS, ES, EP, EPS, E, LVH2O
   PARAMETER (eps=0.622 )

!     ABOUT THE PARAMETER:
!      
!     eps ---------- (Water)/(dry air) molecular mass ratio, epsilon
!
!    function E(T) = Sat. vapor pressure (in Pascal) at 
!                    temperature T (uses Clausius-Clapeyron).
          Es = E(T,LVH2O)
!  CONVERT REL. HUMIDITY (%) TO THE FRACTIONAL VALUE
          RHF = RH/100.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C      CALCULATE SATURATION MIXING RATIO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!      QS = 0.622 * ES /P   was substituted by a more precise
! formula:                              -PABLO J. GRUNMANN, 05/28/98.
        QS = 0.622 * ES /(P - (1.-0.622)*ES)

!
!  CONVERSION FROM REL. HUMIDITY:
!     (Rogers, pg. 17)
!
        EP = (P*Es*RHF)/(P - Es*(1. - RHF))
        QD = eps*EP/(P - (1. - eps)*EP)
!     
RETURN
END

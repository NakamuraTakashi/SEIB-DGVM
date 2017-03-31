!**************************************************************************
!* Spatially Explicit Individual-Based Dynamic Global Vegetation Model    *
!* SEIB-DGVM ver. 2.81                                                    *
!*                                                                        *
!*   All rights are reserved by Dr. Hisashi SATO (JAMSTEC)                *
!*                                                                        *
!*   Latest information and code using policy can be obtained at          *
!*   http://seib-dgvm.com/                                                *
!**************************************************************************

!*************************************************************************************************
! MAIN SIMULATION LOOP
!*************************************************************************************************
SUBROUTINE main_loop ( LAT, LON, GlobalZone, YearMaxClimate, YearMaxCO2, &
           tmp_air, prec, cloud, wind, humid, tmp_air_range, tmp_soil, &
           aco2_annual, ALT, Albedo_soil0, W_fi, W_wilt, W_sat, W_mat )
   
!_____________ Global Variables
   USE data_structure
   USE time_counter
   USE vegi_status_current1
   USE vegi_status_current2
   USE grid_status_current1
   USE grid_status_current2
   implicit none
   
!_____________ Set Local Parameters
   !Interval for intensive computation of radiation distribution (day)
   integer,parameter::Days_LightComp_Int = 14
   
   !Default CO2 concentration in the atmosphere @ ppm
   real,parameter::CO2atm_default = 368.0
   
!_____________ Set Augments
!Coordination
   real   ,intent(IN):: LAT       !latitude  (degree)
   real   ,intent(IN):: LON       !longitude (degree)
   
!ID number of global zone
   integer,intent(IN):: GlobalZone 
   
!Length of Inputted Climate Data (yr)
   integer,intent(IN):: YearMaxClimate
   
!Length of Inputted CO2 Data (yr)
   integer,intent(IN):: YearMaxCO2
   
!Climatic data
   real,dimension(Day_in_Year, YearMaxClimate),intent(IN) :: &
    tmp_air       , & !Surface air temperature (Celcius)
    prec          , & !Precipitation (mm day-1)
    cloud         , & !Cloudness (fraction)
    wind          , & !Wind velocity (m s-1)
    humid         , & !Specific humidity (kg kg-1)
    tmp_air_range     !Daily range of tmp_air (Celcius)
   
   real,dimension(Day_in_Year, YearMaxClimate, NumSoil),intent(IN):: &
    tmp_soil      !Soil temperature for each layers (Celcius)
   
!Atomospheric CO2 time-series @ ppm
   real,dimension(YearMaxCO2),intent(IN)::aco2_annual !Atomospheric co2 concentration (ppm)
   
!Location data
   real,intent(IN):: &
    ALT          ,& !Altitude  (m above MSL)
    Albedo_soil0 ,& !Soil albedo
    W_fi         ,& !Filed capacity   (m3/m3, 0.0 -> 1.0)
    W_wilt       ,& !Wilting point    (m3/m3, 0.0 -> 1.0)
    W_sat        ,& !Saturate point   (m3/m3, 0.0 -> 1.0)
    W_mat           !Matrix potential 
   
!_____________ Set Local variables

!Daily mean meteological variables
   real :: &
    tmp_air_Today       , & !Surface air temperature (Celcius)
    prec_Today          , & !Precipitation           (mm day-1)
    cloud_Today         , & !Cloudness (fraction)
    wind_Today          , & !Wind velocity     (m s-1)
    humid_Today         , & !Wind velocity     (m s-1)
    tmp_air_range_Today     !Daily range of tmp_air  (Celcius)
   
   real,dimension(NumSoil):: &
    tmp_soil_Today=0.0   !Soil temperature for each layers (Celcius)
   
!Other variables
   integer year_climate !Current climate year
   integer i, p, no     !Loop counters
   real    x, y         !For General Usage
   
!_____________ Intialize variables
!Initialize variables
   Call init_value (W_fi, tmp_air(:,1), tmp_soil(:,1,:), prec(:,1))
   Call radiation_seasonal_change (LAT)
   
!Read spinup files
   Spinup_year  = 0
   if (Flag_spinup_read) then
      open (File_no(1), file=Fn_spnin)
      Call spinup_in (File_no(1))
      close (File_no(1))
   endif
   
!_____________ Open output files
   open ( File_no( 1), file = 'output_log.txt'          )
   open ( File_no( 2), file = 'output.txt'              )
   open ( File_no( 3), file = 'output_climate.txt'      )
   open ( File_no( 4), file = 'output_air.txt'          )
   open ( File_no( 5), file = 'output_radiation.txt'    )
   open ( File_no( 6), file = 'output_water.txt'        )
   open ( File_no( 7), file = 'output_grass.txt'        )
   open ( File_no( 8), file = 'output_lai.txt'          )
   open ( File_no( 9), file = 'output_cflux.txt'        )
   open ( File_no(10), file = 'output_netradiation.txt' )
   open ( File_no(11), file = 'output_wflux.txt'        )
   open ( File_no(12), file = 'output_ld_vertical.txt'  )
   open ( File_no(13), file = 'output_annual.txt'       )
   open ( File_no(14), file = 'output_forest.txt'       )
   open ( File_no(15), file = 'output_biomass.txt'      )
   
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
   if (Logging==1) write (File_no(1),*) 'Daily update of field statuses'
   
!Time counters1
   doy = doy + 1
   if (doy == Day_in_Year + 1) then
      doy   = 1
      year  = year + 1
      
      year_climate = year_climate + 1
      if (year_climate==YearMaxClimate+1) year_climate=1
   endif
   
!Time counters2 (wild fire and phenology related)
   dfl_fire = dfl_fire + 1 !Day from the last fire
   
   Do p=1, PFT_no
      dfl_leaf_onset(p) = dfl_leaf_onset(p) + 1 !Day from the last leaf onset
      dfl_leaf_shed (p) = dfl_leaf_shed (p) + 1 !Day from the last leaf shedding
   End do
   
!_____________ Prepare Climatic data for this cycle
!Daily mean meteological properties
   tmp_air_Today       = tmp_air       (doy, year_climate)
   prec_Today          = prec          (doy, year_climate)
   cloud_Today         = cloud         (doy, year_climate)
   wind_Today          = wind          (doy, year_climate)
   humid_Today         = humid         (doy, year_climate)
   tmp_air_range_Today = tmp_air_range (doy, year_climate)
   
   do i=1, NumSoil
    tmp_soil_Today(i)  = tmp_soil (doy, year_climate, i)
   end do
   
!_____________ Daily update of metabolic status
!stat_water(1:PFT_no), a vegetation growth limitter due to soil water shortage
   Do p = 1, PFT_no
      x  = 0.0
      y  = 0.0
      no = 0  
      
      do i=1, max(1,RootDepth(p))
      if (tmp_soil_Today(i)<=0.0) cycle
         no = no + 1
         x  =     x + (pool_w(i)/Depth - W_wilt) / max(W_fi-W_wilt, 0.001)    !Method1
         y  = max(y,  (pool_w(i)/Depth - W_wilt) / max(W_fi-W_wilt, 0.001) )  !Method2
      enddo
      x = x / max(1,no)
      
      x = max(min(x, 1.0), 0.0)
      y = max(min(y, 1.0), 0.0)
      
      if (Life_type(p)==2) then
         stat_water(p) = y !for larch
      else
         stat_water(p) = x !for other PFTs
      endif
      
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
   if (Logging==1) write (File_no(1),*) 'Daily update of field records'
   
!Daily update of Climate statuses
!Calculate climatic statistics, their running record, and their running mean
   Call stat_climate (prec_Today, tmp_air_Today, tmp_soil_Today)
   
!Daily update of Carbon pools and Carbon Fluxes statuses
!Reset variables, update their record in array variables, and calculate their running mean
   Call stat_carbon ()
   
!Daily update of vegetation statuses
!Reset variables, update their record in array variables, and calculate their running mean
   Call stat_vegetation ()
   
   !mortality related
   if (doy==1) then
      mort_regu1(:) = 0.0 !NPP annual (g / individual)
      mort_regu2(:) = 0.0 !average leaf area of last year (m2/day) {update on: growth_wood}
      mort_regu4(:) = 0.0 !stem diameter increament in last year (m year-1)
   endif
   
!_____________ DAILY PSYSICAL PROCESSES
   if (Logging==1) write (File_no(1),*) 'DAILY PSYSIOLOGICAL PROCESSES'
   
   !Specify CO2 concentration in the atmosphere
    co2atm = CO2atm_default                                 !Fixed @ Default value
   !co2atm = aco2_annual(max(year_climate+51, YearMaxCO2))  !Same as year_climate
   
   !Calculate variables in relation to atmospheric physics (e.g. Air pressure, Vapor density)
   Call air (tmp_air_Today, humid_Today, ALT)
   
   Call radiation (LAT, cloud)
   par_RunningRecord(2:Day_in_Year) = par_RunningRecord(1:Day_in_Year-1)
   par_RunningRecord(1)             = par
   
   Call diffused_radiation ()
   
   if (mod(doy, Days_LightComp_Int)==1) then
      !Compute relative intensity of direct radiation for each crown disk of each tree.
      Call direct_radiation ()
      
      !Compute relative intensity of radiation for each cell of tree establishment.
      Call floor_radiation ()
      
      Call crown_coverage ()
      
   endif
   
   !Calculate coefficients ir_tree and ir_grass, which describe fractions of incoming radiation
   !absorbed by crown layer and grass layer, respectively. 
   !This subroutine also calculates a variable par_grass, 
   !intensity of photosynthesis active radiation on the top of grass layer.
   Call ir_index ()
   
   !Wild fire subroutines
   if     (GlobalZone==1) then
     !African continent
     Call fire_regime2 (wind)
   else
     !Default
     Call fire_regime (W_fi)
   endif
   
!*****************************************************************************
   
!_____________ DAILY BIOLOGICAL PROCESSES
   if (Logging==1) write (File_no(1),*) 'DAILY BIOLOGICAL PROCESSES'
   
   !Calculate photosynthesis rate controlling variables
   Call photosynthesis_condition (tmp_air_Today)
   
   !Photosynthesis process
   Call photosynthesis ()
   
   !Calculate optimal leaf area index for grass layer
   Call lai_optimum (tmp_air_Today, sum(tmp_soil_Today(1:5))/5.0 )
   
   !Maintenance respiration
   Call maintenance_resp (tmp_air_Today, sum(tmp_soil_Today(1:5))/5.0 )
   
   !Turnover of plant organs
   Call turnover ()
   
   !Penology controller, and biological processes before and after penology change
   !(such as gradual release of stock biomass after onset of foliage phase).
   Call leaf_season (LAT, tmp_soil_Today )
   
   !Daily growth procedure for woody PFTs (foliation, fine root growth, and reproduction)
   Call growth_wood ()
   
   !Growth and reproduction procedure for grass PFTs.
   Call growth_grass ()
   
   !Calculate free space around each individual tree. 
   !This information will be employed by following subroutine growth_trunc, 
   !where stem and crown expansion occurs.
   Call spatial_limitation ()
   
   !Expands stem diameter, which cause stem height and crown diameter increases
   Call growth_trunk ()
   
   !Heterotrophic respiration (Litter and soil organic matter decomposition)
   Call decomposition (W_fi)
   
!_____________ ANNUALY BIOLOGICAL PROCESSES
   if (Logging==1) write (File_no(1),*) 'ANNUALY BIOLOGICAL PROCESSES'
   
IF ( doy==1 ) then
   !Listing all PFT that exist in this virtual forest
   Call pft_present ()
ENDIF
   
IF ( doy==Day_in_Year ) then
   
   !Determine which tree die
   Call mortality ()
   
   !Adjust crown depth by purging crown disks from the bottom of the crown layer
   Call crown_adjust ()
   
   !Cause horizontal movement of crown of each tree to the open direction
   Call crown_shake ()
   
   !Establishment of woody PFTs
   x = sum(prec_RunningRecord(1:Day_in_Year))                  !annual precipitation (mm/year)
   y = sum(tmp_air_RunningRecord(1:Day_in_Year)) / Day_in_Year !annual mean temperature (Celcius)
   if ( x > max(100.0, 20.0*y) ) then                          !from Koppen's criteria (1936)
      
      !Determine ground mesh point, where newly saplings can establish.
      !This subroutine should be called before subroutine 'establish'.
      Call ground_vacant ()
      
      !Recruitment of trees
      Call establish (GlobalZone)
      
   end if
   
   !Compute fraction of tree crown coverage.
   Call crown_coverage ()
   
   !Determine dominant grass PFT in the next year
   Call grass_priority ()
   
   !Determine current biome type of this grid cell
   Call biome_determine ()
   
END IF
   
!_____________ DAILY UPDATE NET-RADIATION & SOIL-WATER STATUS
   if (Logging==1) write (File_no(1),*) 'DAILY UPDATE NET-RADIATION & SOIL-WATER STATUS'
   
   !Compute Albedo
   Call albedo_calc (Albedo_soil0)
   
   !Compute net radiation
   Call net_radiation (tmp_air_Today, cloud_Today)
   
   !Compute water cycle
   Call waterbudget &
   (W_fi, W_wilt, prec_Today, wind_Today, tmp_air_Today, tmp_soil_Today(:) )
   
!_____________ DAILY INCREMENT of LITTER AND FUEL POOLS
   ! Increment of litters pools and fuel pools are basically conducted here
   ! by considering litter fluxes, which are computed in descendent subroutines.
   ! This way of convergent-computation is for handling 
   ! a variety of litter pools and fuel pools in the model.
   
   pool_litter_leaf  = max(0.0, pool_litter_leaf  + flux_litter_leaf )
   pool_litter_trunk = max(0.0, pool_litter_trunk + flux_litter_trunk)
   pool_litter_root  = max(0.0, pool_litter_root  + flux_litter_root )
   pool_litter_ag    = max(0.0, pool_litter_ag    + flux_litter_ag   )
   pool_litter_bg    = max(0.0, pool_litter_bg    + flux_litter_bg   )
   
   pool_fuel_standT  = max(0.0, pool_fuel_standT  + flux_litter_leaf )
   pool_fuel_standG  = max(0.0, pool_fuel_standG  + flux_litter_ag   )
   
!_____________ MAKE OUTPUT FILES (Write simulation results in output files)
IF (Flag_output_write) then
   
   !Daily output (everyday)
   Call output_for_viewer &
   (File_no(2), LAT, LON, ALT, Albedo_soil0, W_fi, W_wilt, W_sat, W_mat, &
   cloud, prec_Today, humid_Today, wind_Today, tmp_air_Today, sum(tmp_soil_Today(1:5))/5.0 )
   
   !Daily output (on each specifyed year)
   Call output_climate &
   (File_no(3), &
   cloud, prec_Today, humid_Today, wind_Today, tmp_air_Today, sum(tmp_soil_Today(1:5))/5.0)
   Call output_air          (File_no(4))
   Call output_radiation    (File_no(5))
   Call output_water        (File_no(6), W_fi, W_wilt, W_sat)
   Call output_grass        (File_no(7))
   Call output_lai          (File_no(8))
   Call output_cflux        (File_no(9))
   Call output_netradiation (File_no(10))
   Call output_wflux        (File_no(11), prec_Today)
   
   !Monthly output
   if ( Day_of_Month(doy) == Day_in_month(Month(doy)) ) then
      Call output_ld_vertical (File_no(12))
   endif
   
   !Annual output (At the end of the year)
   if (doy == Day_in_Year) then
!  if (doy == 213) then !@ 1 Aug
      Call output_annual (File_no(13))
      Call output_forest (File_no(14))
   endif
   
   !Annual output (At the middle of gwoeing season in Northern Hemisphere)
   if (doy == 195) then
      Call output_biomass(File_no(15))
   endif
   
ENDIF

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
   open(File_no(1), file = Fn_spnout)
   Call spinup_out (File_no(1))
   close (File_no(1))
END IF

END SUBROUTINE main_loop

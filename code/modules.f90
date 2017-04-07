!**************************************************************************************************
! Define common parameters and common functions
!**************************************************************************************************
MODULE data_structure

!_____________ Patameters that must be defined at the beginning, Part 1
integer,parameter::PFT_no  = 4   !Number of plant functional types!!!>>>>>>>>>>>>TN:changed 16 -> 4

integer,parameter::Max_no  = 10000  !maximum individual number in a plot

integer,parameter::Max_hgt = 600  !Maximum number of vertical layer (in STEP).
! This parameter is referred when calculating sunlight distribution in virtual forest.
! Specify larger value of maximum tree height in the forest.

!integer,parameter::Dived   = 30   !Resolution of establishment site for woody PFTs.  !!!>>>>>>>>>>>>TN:rm not used for mangrove
! The specified value divides one side of virtual forest.
!For example, when this value is 60 and one side length of virtual forest is 30.0m,
! then establishment site exist for 0.5m interval
! (for this case, 60*60 = 3600 establishment site exist).

!integer,parameter::DivedG = 30 !Resolution of Grass cell.   !!!>>>>>>>>>>>>TN:rm not used for mangrove

integer,parameter::NumSoil = 30 !Number of Soil layer

!Number of light intensity class for computing photosynthesis rate and transpiration rate (Only used for the original way for phtosynthesis)
integer,parameter::MaxParClass = 10

!Number of light intensity class for computing photosynthesis rate and transpiration rate (Only used for Farquhar's method for phtosynthesis)
integer,parameter::MaxLightClass = 30 

integer,parameter::Logging = 0 !Make log file or not; 0->off, 1->on

!_____________ Fixed patameters (no needs to modify)
real,parameter   ::PI          = 3.141592 !pai
real,parameter   ::DtoR        = PI/180.0 !angle conversion, from degree to radian [dTr]
real,parameter   ::RtoD        = 180.0/PI !angle conversion, from radian to degree [rTd]
real,parameter   ::ZAT         = 273.15   !zero degree centigrade in absolute temperature
real,parameter   ::GasConst    = 8.31451  !gas constant (J K-1 mol-1)
integer,parameter::Day_in_Year = 365      !number of day containing a year (day)

!number of day containing each month (day)
integer,dimension(1:12),parameter::&
Day_in_month = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

integer,dimension(1:Day_in_Year),parameter::& !for 'doy' -> 'day_of_month'
Day_of_Month = (/&
1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,&
1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,         &
1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,&
1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,   &
1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,&
1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,   &
1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,&
1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,&
1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,   &
1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,&
1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,   &
1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31/)

integer,dimension(1:Day_in_Year),parameter::& !for 'doy' -> 'month'
Month = (/&
 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,          &
 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,    &
 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,    &
 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, &
 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, &
 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,    &
10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10, &
11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,    &
12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12 /)

!_____________ Patameters that whill be intialized with paramete.txt
!Control
integer,save:: Simulation_year   !simulation year (yr)

logical,save:: Flag_spinup_read         !flag: restart from spinup files or not
logical,save:: Flag_spinup_write        !flag: make spinup files or not
logical,save:: Flag_output_write        !flag: make output files or not
logical,save:: Flag_randomization       !flag: True  -> employ different random seed for each run

!integer,save          ::Max_loc         !side lenght of forest stand (m)  !!!>>>>>>>>>>>>TN:rm
real                  ::Depth           !depth of each soil layer
real                  ::STEP            !vertical resolution of tree canopy (m)
real                  ::C_in_drymass    !proportion of carbon in biomass (W/W)

integer,dimension (25)::File_no         !device number for I/O

character*128::Fn_climate  !input  file name for climate data
character*128::Fn_CO2      !input  file name for CO2 data (annualy)
character*128::Fn_location !input  file name for soil data
character*128::Fn_spnin    !input  file name for spinup data
character*128::Fn_spnout   !output file name for spinup data

!PFT type
integer,dimension(PFT_no)::Life_type     !0  : Other tree species
                                         !1  : tropical rain trees (TrBE)
                                         !2  : larch (BoNS)
                                         !3  : C3 grass
                                         !4  : C4 grass

integer,dimension(PFT_no)::Phenology_type !0  : evergreen
                                          !1-3: temperature-controlling deciduous
                                          !4  : water-controlling deciduous

!Respiration
real,dimension(PFT_no)::RM           !maintenance respiration rate: foliage (gDM/gN day-1 @15C)
real,dimension(PFT_no)::RG_f         !Specific growth respiration rate for foliage (g dm /g dm)
real,dimension(PFT_no)::RG_f_suck    !Fraction of leaf mass suck when defoliation  (g dm /g dm)
real,dimension(PFT_no)::RG_s         !Specific growth respiration rate for stem    (g dm /g dm)
real,dimension(PFT_no)::RG_r         !Specific growth respiration rate for root    (g dm /g dm)
real                  ::RG_stock_in  !from sucrose to polysaccharide       (g dm /g dm)
real                  ::RG_stock_out !from polysaccharide to sucrose       (g dm /g dm)

!Turnover
real,dimension(PFT_no)::TO_f !Fixed turn over time for Foliage (yr-1)
real,dimension(PFT_no)::TO_s !Fixed turn over time for Sapwood (yr-1)
real,dimension(PFT_no)::TO_r !Fixed turn over time for Root (yr-1)

!Metabolic Characteristic
real   ,dimension(PFT_no)::ALM1 !Allometry index of LA vs dbh_sapwood (m2/m)
real   ,dimension(PFT_no)::ALM2 !Allometry index of crown_area
real   ,dimension(PFT_no)::ALM3 !Allometry index of mass_trunk (g dm/ m3)
real   ,dimension(PFT_no)::ALM4 !Allometry index of crown adjustment (proportion)
real   ,dimension(PFT_no)::ALM5 !Allometry index, diameter proportion of sapwood
real   ,dimension(PFT_no)::ALM6 !Allometry index, crown depth fraction to tree height (fraction)
real   ,dimension(PFT_no)::DBH_limit !limitation on trunk diameter (m)
real   ,dimension(PFT_no)::HGT_max   !Maximum tree height (m)
real   ,dimension(PFT_no)::HGT_s     !Initial value of relative growth rate of height to dbh (m/m)
real   ,dimension(PFT_no)::CD_max    !Maximum crown diameter (m)
real   ,dimension(PFT_no)::LA_max    !Maximum leaf area per canopy surface (m2/ m2)
real   ,dimension(PFT_no)::FR_ratio  !ratio of mass_leaf vs mass_root (g dm/g dm)
integer,dimension(PFT_no)::RootDepth !(*Depth mm)

real   ,dimension(PFT_no)::PN_f    !proportion of nitrogen:leaf      (gN / gDM)
real   ,dimension(PFT_no)::PN_s    !proportion of nitrogen:sapwood   (gN / gDM)
real   ,dimension(PFT_no)::PN_r    !proportion of nitrogen:fine root (gN / gDM)
real   ,dimension(PFT_no)::SLA     !specific leaf area          (one sided m2 / g DM)

!Photosynthesis
real,dimension(PFT_no)::Pmax    !potential maximum photosynthesisrate (micro mol CO2 m-2 s-1)
real,dimension(PFT_no)::EK0     !light attenuation coefficient to vertical direction (no dimension)
real,dimension(PFT_no)::Lue0    !control light dependence coefficient (mol CO2 mol photon-1)
real,dimension(PFT_no)::Topt0   !optimum temperature (deg C)
real,dimension(PFT_no)::Tmin    !minimum temperature (deg C)
real,dimension(PFT_no)::Tmax    !maximum temperature (deg C)
real,dimension(PFT_no)::GS_b1   !parameters of stomatal conductance (mol H2O m-2 s-1)
real,dimension(PFT_no)::GS_b2   !parameters of stomatal conductance (no dimension)
real,dimension(PFT_no)::GS_b3   !parameters of stomatal conductance (hPa)
real,dimension(PFT_no)::KM      !effect of photosynthesis on intercellular CO2 concentration (ppmv)
real,dimension(PFT_no)::CO2cmp0 !CO2 compensation point (ppmv)

!Dynamic characters
real   ,dimension(PFT_no)::M1        !Mortality; asymptotic maximum mortality rate
real   ,dimension(PFT_no)::M2        !Mortality; a parameter of back ground mortality
real   ,dimension(PFT_no)::M4        !Mortality; a parameter of back ground mortality
real   ,dimension(PFT_no)::M5        !Mortality; a parameter of back ground mortality
real   ,dimension(PFT_no)::Msal1     !Mortality; a parameter of salinity effect
real   ,dimension(PFT_no)::Msal2     !Mortality; a parameter of salinity effect
real   ,dimension(PFT_no)::TC_min         !Minimum coldest month temperature for persisting

real   ,dimension(PFT_no)::P_establish !Establishment probability in safe_site (year-1 m-2)
real   ,dimension(PFT_no)::TC_max   !Maximum coldest month temperature for establishment (cecius)
real   ,dimension(PFT_no)::GDD_max  !Maximum degree-day sum (5 Celcius base) for establishment
real   ,dimension(PFT_no)::GDD_min  !Minimum degree-day sum (5 Celcius base) for establishment
real   ,dimension(PFT_no)::PAR_min  !Minimum midday PAR (micro mol photon m-2 s-1) to establish
integer,dimension(PFT_no)::DM_max   !Maximum drought month for establishment (month, 1-12)
real   ,dimension(PFT_no)::AGE_max  !Maximum age above which trees die (year)

integer                  ::Est_scenario !Scenario for establishment (0 or 1 or 2 or 3)
   !0 -> no trees can establish
   !1 -> only specified woody PFTs can establish
   !2 -> every potentially eatablishable PFT have same chance of establishment
   !3 -> in proportion of existing biomass after specific year
   !4 -> in proportion of existing biomass after specific year,
   !     while little portion of establishment was randomly selected
   !     from all potentialy establishable PFTs
   
integer                  ::Est_year_change !for Scenario 3 and 4:
                                           !lenght of year until establishment pattern will change
logical,dimension(PFT_no)::Est_pft_OnOff   !for Scenario 1: establishment switch
real                     ::Est_frac_random !for Scenario 4: fraction of random establishment

!Fire regime
real                  ::Fuel_min            !Minimum fuel load for fire spread (g C m-2)
real,dimension(PFT_no)::Moisture_extinction !Litter moisture of extinction (proportion)
real,dimension(PFT_no)::M3            !Fire resistance parameter for woody PFTs (no dimension)

!Decomposition
real TO_litter         !turnover rate of litter   (yr-1)
real TO_fast           !turnover rate of fast decomposition pool of soil organic matter (yr-1)
real TO_slow           !turnover rate of slow decomposition pool of soil organic matter (yr-1)
real F_air             !fraction of decomposed litter that directly goes into air
real F_inter           !fraction of soil incorporated matter that goes into pool_som_inter
real,dimension(PFT_no)::F_resp !fraction of decomposed litter that directly goes into air
real,dimension(PFT_no)::F_fast !fraction of soil incorporated matter that goes into pool_som_inter

!For employing different random seed for each run (by Shigeki Ikeda @ Kyodo Univ.)
integer :: size_count,clock
integer :: seedsize
integer,allocatable :: seed(:)

!***********************  set NameList  *************************
    namelist /Control/        &
        Simulation_year, Flag_spinup_read, Flag_spinup_write, Flag_output_write, &
        Flag_randomization, &
        Max_loc, Depth, STEP, &
        C_in_drymass, File_no, &
        Fn_climate, Fn_CO2, Fn_location, Fn_spnin, Fn_spnout
    
    namelist /PFT_type/       &
        Life_type, Phenology_type
    
    namelist /Respiration/    &
        RM, RG_f, RG_f_suck, RG_s, RG_r, RG_stock_in, RG_stock_out
    
    namelist /Turnover_n/       &
        TO_f,  TO_s,  TO_r
    
    namelist /Metabolic/      &
        ALM1, ALM2, ALM3, ALM4, ALM5, ALM6,&
        CD_max, LA_max, FR_ratio, RootDepth, DBH_limit, HGT_max, HGT_s, SLA, PN_f, PN_s, PN_r
    
    namelist /Assimilation/ &
        Pmax, EK0, Lue0, Topt0, Tmin, Tmax, &
        GS_b1, GS_b2, GS_b3, KM, CO2cmp0
    
    namelist /Dynamics/       &
        M1, M2, M4, M5, Msal1, Msal2, TC_min, &
        P_establish, TC_max, GDD_max, GDD_min, PAR_min, DM_max, AGE_max, &
        Est_scenario, Est_year_change, Est_pft_OnOff, Est_frac_random
    
    namelist /Disturbance/    &
        Fuel_min, Moisture_extinction, M3
    
    namelist /Soil_resp/  &
        TO_litter, TO_fast, TO_slow, F_air, F_inter, F_resp, F_fast
    
!_____________ Common functions
CONTAINS

!Maximum tree_height for a given diameter (step) 
INTEGER FUNCTION h(dbh, pft)
   real    dbh !Diamter at Breath Height (m)
   integer pft !PFT
   real    hgt !tree height (m)
   
   !FORSKA way of modeling
   !h = int(   (HGT_max(pft) * (1 - exp( (-HGT_s(pft)*dbh)/HGT_max(pft) )))   /STEP)
   
   hgt = 1.0 / ( 1.0/(HGT_s(pft)*dbh) + 1.0/HGT_max(pft) )
   h = int( (hgt-1.3) / STEP )
   h = max(2, h)
END FUNCTION h

!Stem weight (g DM)
REAL FUNCTION stem_weight(dbh,height,pft)
   real    dbh     !diameter at breast height (m)
   integer height  !height of tree above 1.3 m (STEP)
   integer pft     !pft number
   real    x       !for temporal usage
   
   !default
   x = ALM3(pft) * (1.3+height*STEP) * ((dbh**2) * PI / 4.0)
   
   if ( Life_type(pft)==1 ) then
   !for tropical evergreen trees (by Huth & Ditzer 2000)
!      x = ALM3(pft) * (1.3+height*STEP) * ((dbh**2)*PI/4.0) * (1.0/0.7) * (0.5941-0.0108*log(x))
!      x = ALM3(pft) * (1.3+height*STEP) * ((dbh**2)*PI/4.0) * (1.0/0.7) * (0.5941-0.0108*log(x))
!      x = ALM3(pft) * (1.3+height*STEP) * ((dbh**2)*PI/4.0) * (1.0/0.7) * (0.5941-0.0108*log(x))
!      x = ALM3(pft) * (1.3+height*STEP) * ((dbh**2)*PI/4.0) * (1.0/0.7) * (0.5941-0.0108*log(x))
      
      !dbh‚¾‚¯‚ÌŠÖ”‚Æ‚µ‚ÄŠ²d‚ªŒˆ‚Ü‚é‚Æ‚µ‚½
      x = ALM3(pft)* (1.3+h(dbh,pft)*STEP)* ((dbh**2)*PI/4.0)* (1.0/0.7)* 0.5
      x = ALM3(pft)* (1.3+h(dbh,pft)*STEP)* ((dbh**2)*PI/4.0)* (1.0/0.7)* (0.5941-0.0108*log(x))
      x = ALM3(pft)* (1.3+h(dbh,pft)*STEP)* ((dbh**2)*PI/4.0)* (1.0/0.7)* (0.5941-0.0108*log(x))
      x = ALM3(pft)* (1.3+h(dbh,pft)*STEP)* ((dbh**2)*PI/4.0)* (1.0/0.7)* (0.5941-0.0108*log(x))
      x = ALM3(pft)* (1.3+h(dbh,pft)*STEP)* ((dbh**2)*PI/4.0)* (1.0/0.7)* (0.5941-0.0108*log(x))
      
   elseif ( Life_type(pft)==5 ) then
   !For Tropical Broadleaf trees in Africa
      ![from aDGVM]
      !x = 27523 * 1000.0 *  (dbh**2.55) 
      
      !from Type II model of Chave et al. (2005)
      x = 207.0 * ((ALM3(pft)/1000000.0)**1.036) * ((dbh*100)**2.179) 
      
      !convert AGB to tunc biomass, which includes coarse root, Deans at al. (1996)
      x = x * 1.33
      
   elseif ( Life_type(pft)==2 .or. Life_type(pft)==6 ) then
   !for boreal deciduous needle leaf trees
      x = ALM3(pft) * (1.3+height*STEP) * ((dbh**2) * PI / 6.0)
      
      !Sato et al. (2010) for larch
      !x = 0.19*( (100*dbh)**1.81 ) + 0.0428*( (100*dbh)**1.79 ) + 0.171*( (100*dbh)**1.67 )
      !x = 1000 * x
      !y = -11300 + 4670000 * 3.14*((0.5*dbh)**2)
      !y = y * 1.5 !(convert: AGB-> tunc biomass)
      !x = max(x,y)
   end if
   stem_weight = x
   
END FUNCTION

!Random value generator (output:1.0-0.0)
REAL FUNCTION randf()
   !integer I_seed
   !I_seed = I_seed*48828125
   !if(I_seed.le.0) then
   !   I_seed = (I_seed+2147483647)+1
   !else
   !   I_seed = I_seed
   !end if
   !randf = float(I_seed)/2147483647
   call random_number(randf)
   return
END FUNCTION randf

END MODULE



!**************************************************************************************************
! Define global variables
!**************************************************************************************************
!Simulation loop counter
MODULE time_counter
   integer counter       !counter number of present day
   integer counter_begin !counter number for start of the simulation
   integer counter_end   !counter number for termination of the simulation
   integer doy           !day of the year (1-Day_in_Year)
   integer year          !loop counter for year (1-Simulation_year)
   integer Spinup_year   !Year length of spin up that has been conducted.
                         !This value is used when simulation starts from spin up file.
END MODULE time_counter

!Vegetation status (current)
MODULE vegi_status_current1
   USE data_structure
   
   !Variables for whole vegetation
   integer biome          !biome type (Definition ->subroutine biome_type)
   real    lai            !Total LAI @ update everyday             (m2/m2)
   
   !Variables for each PFT
   logical,dimension(PFT_no)::phenology      !Phenological status (T:foliation, F:dormance)
   logical,dimension(PFT_no)::pft_exist      !flag: PFT exist or not
   integer,dimension(PFT_no)::dfl_leaf_onset !day from last leaf onset
   integer,dimension(PFT_no)::dfl_leaf_shed  !counter number of last leaf shedding
   real   ,dimension(PFT_no)::stat_water     !state of water satisfactory (0.0-1.0)
   
   !Variables for for each individual tree
   logical,dimension(Max_no)::tree_exist     !tree exist flag
   integer,dimension(Max_no)::pft            !Plant Functional Type
   integer,dimension(Max_no)::age            !tree age (year)
   integer,dimension(Max_no)::height         !tree height from 1.3m above ground (STEP)
   integer,dimension(Max_no)::bole           !bole length from 1.3m above ground (STEP)
   integer,dimension(Max_no)::height_limit   !maximum tree height that proximate trees permit(step)
   integer,dimension(Max_no)::flag_suppress
   
   real   ,dimension(Max_no)::dbh_heartwood  !heartwood diameter at 1.3m (m)
   real   ,dimension(Max_no)::dbh_sapwood    !sapwood   diameter at 1.3m (m)
   real   ,dimension(Max_no)::crown_diameter !crown diameter      (m)
   real   ,dimension(Max_no)::crown_area     !crown area          (m2)
   real   ,dimension(Max_no)::bole_x         !x_location of bole  (m)
   real   ,dimension(Max_no)::bole_y         !y_location of bole  (m)
   real   ,dimension(Max_no)::crown_x        !x_location of crown (m)
   real   ,dimension(Max_no)::crown_y        !y_location of crown (m)
   real   ,dimension(Max_no)::radius_limit   !maximum canopy radius that proximate trees permit (m)
   real   ,dimension(Max_no)::la             !leaf area           (m2)
   real   ,dimension(Max_no)::mass_leaf      !biomass of foliage  (gDM/tree)
   real   ,dimension(Max_no)::mass_trunk     !biomass of trunk    (gDM/tree)
   real   ,dimension(Max_no)::mass_root      !biomass of fine root(gDM/tree)
   real   ,dimension(Max_no)::mass_stock     !biomass for stock   (gDM/tree)
   real   ,dimension(Max_no)::mass_available !biomass available   (gDM/tree)
   real   ,dimension(Max_no)::mort_regu1               !sum of NPP within the previous year (gDM/tree)
   real   ,dimension(Max_no)::mort_regu2               !average leaf area of last year (m2)
   real   ,dimension(Max_no)::mort_regu4               !stem diameter increament in last year (m year-1)
   real   ,dimension(Max_no)   ::npp_crowntop    !annual NPP at the top CrownDisk (gDM/year/step)
   real   ,dimension(Max_no,10)::npp_crownbottom !annual NPP of btm CrownDisks    (gDM/year/step)
   
   !Variables for each floor cell
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:rm
!   real,dimension(DivedG,DivedG):: gmass_leaf      !c3/c4 grass foliage   biomass (gDM/cell)
!   real,dimension(DivedG,DivedG):: gmass_root      !c3/c4 grass root      biomass (gDM/cell)
!   real,dimension(DivedG,DivedG):: gmass_available !c3/c4 grass available biomass (gDM/cell)
!   real,dimension(DivedG,DivedG):: gmass_stock     !c3/c4 grass stock     biomass (gDM/cell)
!   real,dimension(DivedG,DivedG):: lai_grass       !c3/c4 leaf area index of grass (m2/m2)
!   
!   real,dimension(20,DivedG,DivedG)::lai_opt_grass_RunningRecord
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:rm
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   real,allocatable,dimension(:,:):: gmass_leaf      !c3/c4 grass foliage   biomass (gDM/cell)
   real,allocatable,dimension(:,:):: gmass_root      !c3/c4 grass root      biomass (gDM/cell)
   real,allocatable,dimension(:,:):: gmass_available !c3/c4 grass available biomass (gDM/cell)
   real,allocatable,dimension(:,:):: gmass_stock     !c3/c4 grass stock     biomass (gDM/cell)
   real,allocatable,dimension(:,:):: lai_grass       !c3/c4 leaf area index of grass (m2/m2)
   
   real,allocatable,dimension(:,:,:)::lai_opt_grass_RunningRecord
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   
   !Running Records of plant properties
   logical,dimension(Day_in_Year,PFT_no)::phenology_RunningRecord  !RR of Phenology Flag
   real   ,dimension(Day_in_Year,PFT_no)::npp_RunningRecord        !RR of NPP (g DM/ stand/ day)
   real   ,dimension(Day_in_Year,PFT_no)::gpp_RunningRecord        !RR of GPP (g DM/ stand/ day)
   real   ,dimension(Day_in_Year,PFT_no)::stat_water_RunningRecord !RR of Stat water
   real   ,dimension(Day_in_Year,PFT_no)::lai_RunningRecord        !RR of LAI (m2/m2)
   
   real,dimension(20,PFT_no)::mass_sum_RunningRecord !20yrs RR of biomass for each PFT (gDM/m2)
   
END MODULE vegi_status_current1



MODULE vegi_status_current2
!following variables are not necessary to wirte in spinup files
   USE data_structure
   
   !Constants, which are initialized at the beggining of each simulation run
   integer C3g_no       !PFT number of C3 grass
   integer C4g_no       !PFT number of C4 grass
   
   !Metabolic status
   real                    canopy_cond  !canopy conductance (mol H2O m-2 s-1)
   
   real,dimension(PFT_no,MaxLightClass)::assim_rate !Carbon assimiration rate (gDM m-2 TimeStep-1)
   real,dimension(PFT_no,MaxLightClass)::gtc_rate   !Vapor conductances on leaf surface (mol H2O m-2 s-1)
   
real,dimension(PFT_no)::eK       !light attenuation coefficient to sun direction, no dimension
real,dimension(Max_no)::lue      !light use efficiency (mol CO2 mol photon-1)
real,dimension(Max_no)::co2cmp   !CO2 compensation point (ppmv)
real,dimension(Max_no)::psat     !light-saturated photosynthesis rate (micro mol CO2 m-2 s-1)
real,dimension(Max_no)::npp_crownbottom_daily
                                 !daily stat_leaf on the bottom of CrownDisk

real,dimension(10)::lue_grass    !light use efficiency for each par intensity class (mol CO2 mol photon-1)
real,dimension(10)::co2cmp_grass !CO2 compensation point for each par intensity class (ppmv)
real,dimension(10)::psat_grass   !light-saturated photosynthesis rate for each par intensity class
                                 ! (micro mol CO2 m-2 s-1)   
   
   !Carbon flux
   real,dimension(PFT_no)::gpp           !GPP for each pft (g dm / stand / day  )     (reset@everday)
   real,dimension(PFT_no)::npp           !NPP for each pft (g dm / stand / day  )     (reset@everday)
   real,dimension(Max_no)::gpp_daily_ind !GPP for each individual (g dm / tree / day) (reset@everday)
   
   real,dimension(Max_no)::resp_trunk    !respiration of woody trunk       (gDM/tree/day)
   real,dimension(Max_no)::resp_leaf     !respiration of woody leaf        (gDM/tree/day)
   real,dimension(Max_no)::resp_root     !respiration of woody fine root   (gDM/tree/day)
   real                    resp_grass_ag !respiration of grass aboveground (gDM/stand/day)
   real                    resp_grass_bg !respiration of grass underground (gDM/stand/day)
   
   !Litter flux
   real flux_litter_trunk !litter flux of woody trunk       (gDM/stand/day)
   real flux_litter_leaf  !litter flux of woody leaf        (gDM/stand/day)
   real flux_litter_root  !litter flux of woody fine root   (gDM/stand/day)
   real flux_litter_ag    !litter flux of grass aboveground (gDM/stand/day)
   real flux_litter_bg    !litter flux of grass underground (gDM/stand/day)
   
   !Patch status
!   logical,dimension(Dived,Dived)::patch_vacant !safe site for establishment of woody PFTs !!!>>>>>>>>>>>>TN:rm
   logical,allocatable,dimension(:,:)::patch_vacant !safe site for establishment of woody PFTs !!!<<<<<<<<<<<<TN:add
   
END MODULE vegi_status_current2



!Grid status1 (current status)
MODULE grid_status_current1
   USE data_structure
   
   !Varibles for whole vegetation
   integer dfl_fire    !day from last fire (day)
   integer fire_number !ccumlative fire intensity
   
   real tmp_coldest_20yr_ave !20yr running average of tmp_coldest_month (Celcius)
   real tmp_hottest_20yr_ave !20yr running average of tmp_hottest_month (Celcius)
   real gdd_20yr_ave         !20yr running average of growing-degree-day-sum
   
   real pool_litter_trunk    !litter pool in a stand (trunk of woody PFTs)         (gDM / stand)
   real pool_litter_leaf     !litter pool in a stand (leaves of woody PFTs)        (gDM / stand)
   real pool_litter_root     !litter pool in a stand (fine root of woody PFTs)     (gDM / stand)
   real pool_litter_ag       !litter pool in a stand (aboveground of grass)        (gDM / stand)
   real pool_litter_bg       !litter pool in a stand (belowground of grass)        (gDM / stand)
   
   real pool_som_int         !pool of s.o.m. with intermediated decomposition rate (gDM / stand)
   real pool_som_slow        !pool of s.o.m. with slow decomposition rate          (gDM / stand)
   
   real pool_fuel_standT !standing dead biomass of tree  leaves (for Africa extension)(gDM / stand)
   real pool_fuel_standG !standing dead biomass of grass leaves (for Africa extension)(gDM / stand)
   
   real,dimension(1:NumSoil):: pool_w    !soil water content of each soil layer (mm)
   real                        pool_snow !water equivalent snow depth (mm)
   
   !Running Records for 20yrs
   real,dimension(20):: gdd_20yr_RunningRecord     !Growth Degree day
   real,dimension(20):: tmp_coldest_RunningRecord  !Coldest month temperature (Celcius)
   real,dimension(20):: tmp_hottest_RunningRecord  !Hottest month temperature (Celcius)
   real,dimension(20):: tmp_ave_GrassGrowth_RR     !average temperature of growth phase of grass (Celcius)
   
   !Running Records of Carbon Pool
   real,dimension(Day_in_Year):: pool_c_RR !1yr RR of total carbon pool in the virtual forest (Kg C / m2)
   
   !Running Records of Carbon Flux
   real,dimension(Day_in_Year):: flux_c_uptake_RR !carbon uptake   due to photosyhthesis            (gDM/day/stand)
   real,dimension(Day_in_Year):: flux_c_mnt_RR    !carbon emission due to maintenance   respiration (gDM/day/stand)
   real,dimension(Day_in_Year):: flux_c_gro_RR    !carbon emission due to growth        respiration (gDM/day/stand)
   real,dimension(Day_in_Year):: flux_c_htr_RR    !carbon emission due to heterotrophic respiration (gDM/day/stand)
   real,dimension(Day_in_Year):: flux_c_fir_RR    !carbon emission due to fire incident             (gDM/day/stand)
   
   !Running Records of Water Flux
   real,dimension(Day_in_Year):: flux_ro_RunningRecord !1yr RR runoff water      (mm)
   real,dimension(Day_in_Year):: flux_ic_RunningRecord !1yr RR interception water(mm)
   real,dimension(Day_in_Year):: flux_ev_RunningRecord !1yr RR evaporated water  (mm)
   real,dimension(Day_in_Year):: flux_tr_RunningRecord !1yr RR tranpirated water (mm)
   
   !Running Records of Environmental Factors
   real,dimension(Day_in_Year        )::tmp_air_RunningRecord   !1yr RR (Celcius)
   real,dimension(Day_in_Year,NumSoil)::tmp_soil_RunningRecord  !1yr RR (Celcius), soil_tmp RR@
   real,dimension(Day_in_Year        )::prec_RunningRecord      !1yr RR (mm/day)
   real,dimension(Day_in_Year        )::ev_pot_RunningRecord    !1yr RR of potential avapotranspiration (mm/day)
   real,dimension(Day_in_Year)::par_RunningRecord       !1yrs RR
   real,dimension(Day_in_Year)::pool_w1_RunningRecord   !1yr RR
   
END MODULE grid_status_current1

!Grid status2
!These variables are recalculated at least once in a year without referring old values,
!or they are calcualted each time when the the program runs,
!and hence they are not to be wirtten in spinup files
MODULE grid_status_current2
   USE data_structure
   
   !Physical status of air
   real co2atm      !ambient (canopy) CO2 concentration (ppmv)
   real ap          !air pressure (hPa)
   real vp          !vapour pressure (hPa)
   real vpd         !vapour pressure deficit (hPa)
   real dnsa        !density of air (kg m-3)
   real vp_sat      !saturation vapour pressure (hPa)
   real slope_vps   !slope of saturation vapour pressure related to tempertaure (hPa deg C-1)
   
   !Physical status of radiation
   real sl_hgt           (Day_in_Year)  !solar hight at midday (degree)
   real dlen             (Day_in_Year)  !day length (hour)
   real rad_stratosphere (Day_in_Year)  ! shortwave radiation at the atmosphere-top (W/m2) 
   
   real rad         !shortwave radiation at mid-day over canopy (W/m2)
   real par         !photosynthetically active radiation of mid-day (micro mol photon m-2 s-1)
   real par_direct  !direct PAR on canopy top of mid-day (micro mol photon m-2 s-1)
   real par_diffuse !difussed PAR on canopy top of mid-day (micro mol photon m-2 s-1)
   real albedo_mean !albedo, averaged over the virtual forest
   real albedo_soil !soil surface albedo
   real albedo_leaf !leaf albedo
   real radnet_long !longwave net radiation, upward (W m-2)
   real radnet_veg  !canopy net radiation (W m-2, day-time mean, short+long waves)
   real radnet_soil !soil surface net radiation (W m-2, whole day mean, short+long waves)
   real ir_tree     !radiation-interruption-coefficient by tree corwn (0.0~1.0) lower=more absorp.
   real ir_grass    !radiation-interruption-coefficient by grass leaf
   
   !Growing-degree-day-sum (reset on Jan 1)
   real gdd0      !0 Celcius base
   real gdd5      !5 Celcius base
   
   !Relative PAR intensity (0.0-1.0)
   real,dimension(Max_no,Max_hgt)::par_direct_rel  !direct PAR for each CrownDisk
   real,dimension(Max_hgt)       ::par_diffuse_rel !diffused PAR intensity for each forest layer
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:rm
!   real,dimension(Dived ,Dived ) ::par_floor_rel   !total PAR on each estaclishment cell of forest floor
!   real,dimension(DivedG,DivedG) ::par_grass_rel   !total PAR on each grass         cell of forest floor
!   
!   !1yr sum of PAR intensity for establishment cells
!   real,dimension(Dived, Dived)::sum_par_floor
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:rm
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   real,allocatable,dimension(:,:) ::par_floor_rel   !total PAR on each estaclishment cell of forest floor
   real,allocatable,dimension(:,:) ::par_grass_rel   !total PAR on each grass         cell of forest floor
   
   !1yr sum of PAR intensity for establishment cells
   real,allocatable,dimension(:,:)::sum_par_floor
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   
   !Fraction of crown coverage
   real frac_crown_coverage
   
END MODULE grid_status_current2

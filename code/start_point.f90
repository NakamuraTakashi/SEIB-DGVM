!In case of the Intel Fortran, append the compile option '-assume byterecl' as fillows.
!　　　ifort start_point.f90 -o go.out -assume byterecl

!*************************************************************************************************
! Start up procedures for selected sites
! (read common parameters and climate data)
!*************************************************************************************************
   Include 'modules.f90'        
   Include 'main_point.f90'     
   Include 'initialize.f90'     
   Include 'metabolic.f90'      
   Include 'output.f90'         
   Include 'physics.f90'        
   Include 'population_regu.f90'
   Include 'spatial_calc.f90'   
   Include 'etc.f90'            
   Include 'SFLXALL_SRC_VER_2.7.1.f90'
   
PROGRAM start_point
   USE data_structure
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE vegi_status_current1
   USE vegi_status_current2
   USE grid_status_current2
   USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none
   
!_____________ Set Variables
   
!Set year number of climate CO2 data
   integer YearMaxClimate !year length of the climate data
   integer YearMaxCO2     !year length of the CO2 data
   
!Climate data
   real,allocatable,dimension(:,:)  ::&
   tmp_air       ,& ! 1. Surface air temperature (Celcius)
   tmp_air_range ,& ! 2. Daily range of tmp_air (Celcius)
   prec          ,& ! 3. Precipitation (mm day-1)
   rad_short     ,& ! 4. Shortwave radiation, downward, @ midday (W m-2)
   rad_long      ,& ! 5. Daily mean of longwave radiation, downward (W m-2)
   wind          ,& ! 6. Wind velocity (m s-1)
   rh            ,& ! 7. Relative humidity (%)
   tmp_soil1     ,& ! 8. Soil temperature   0- 10cm depth (Celcius)
   tmp_soil2     ,& ! 9. Soil temperature  10-200cm depth (Celcius)
   tmp_soil3        !10. Soil temperature 200-300cm depth (Celcius)
   
   real,allocatable,dimension(:,:,:)::tmp_soil   !Soil temperature for each layers (Celcius)
   
!Atomospheric CO2 time-series @ ppm
   real,allocatable,dimension(:)::aco2_annual
   
!Location data
   integer Mask         !Land ocean mask (1:land, 0:ocean)
   real    ALT          !altitude (m above MSL)
   real    Albedo_soil0 !albedo, default
   real    W_fi         !filed capacity   (m3/m3, 0.0 -> 1.0)
   real    W_wilt       !wilting point    (m3/m3, 0.0 -> 1.0)
   real    W_sat        !saturate point   (m3/m3, 0.0 -> 1.0)
   real    W_mat        !matrix potential (m, -0.0001 -> -3.0)
   
   integer SoilClass    !Zobler(1986)'s slope type index 1-9 (default = 1)
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
   !  0	    Ocean / Water body
   
!Others
   real    LAT, LON       !latitude and logitude for simulate
   integer GlobalZone     !ID number of global zone
   integer point          !
   integer i, j, count    !for general usage
   real    x              !for general usage
   character(len=256):: chr_tmp !for counting number of data item
   
!_____________ Set location
!Set location
!  LAT    north:+, south:- (decimalized)
!  LON    east:+, west:-   (decimalized)
!   !Sppaskaya-pad
!    LAT    = 62.017  
!    LON    = 129.717

!   !Pasoh
!    LAT    = 2.97  
!    LON    = 102.3 
!!!! Fukido@ Ishigaki Island, Japan <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN changed
    LAT    =  24.48  
    LON    = 124.23 

   
!Set cordinate variables
   point = (90-int(LAT)-1)*360 + int(LON+180) + 1
   
!GlobalZone: Set Location category
   if (LON>=-20 .and. 60>=LON .and. 23.0>=LAT ) then
      !African continent
      GlobalZone = 1
   elseif (LON>=100 .and. 170>=LON .and. 50.0<=LAT) then
      !Eastern Siberia
      GlobalZone = 2
   else
      !Default
      GlobalZone = 0
   endif
   
!_____________ Read Parameters
!Read Parameter files
   write(*,*) 'Read Parameter files' !!!<<<<<<<<<<<<TN:add
   open (1, file='parameter_mangrove.txt', action='READ', status='OLD')
      read ( unit=1, nml=Control)       
      read ( unit=1, nml=PFT_type)      
      read ( unit=1, nml=Respiration)   
      read ( unit=1, nml=Turnover_n)    
      read ( unit=1, nml=Metabolic)     
      read ( unit=1, nml=Assimilation)  
      read ( unit=1, nml=Dynamics)      
      read ( unit=1, nml=Disturbance)   
      read ( unit=1, nml=Soil_resp)     
   close (1)
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:rm
!!Read Location data
!   open ( File_no(1), file=Fn_location, status='OLD')
!   do i=1, point
!      read(File_no(1),*) Mask, ALT, Albedo_soil0, W_sat, W_fi, W_mat, W_wilt, SoilClass
!   end do
!   close( File_no(1) )
!   
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:rm
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:add
!Grid data
   write(*,*) 'Read ROMS netcdf files' !!!<<<<<<<<<<<<TN:add
   Call read_ROMS_files('D:/ROMS/Data/Fukido/fukido_grd_v7.nc'           &
   &                   ,'D:/ROMS/output/Fukido/test/ocean_his2.nc'    &
   &                   , 4, 253  )
   
    Mask    = 1         !Land ocean mask (1:land, 0:ocean)
    ALT     = 0.0       !altitude (m above MSL)
    Albedo_soil0 = 0.1  !albedo, default                      !!!値てきとう
    W_sat   = 0.4565    !saturate point   (m3/m3, 0.0 -> 1.0) !!!値てきとう
    W_fi    = 0.1       !filed capacity   (m3/m3, 0.0 -> 1.0) !!!値てきとう
    W_mat   =-0.3       !matrix potential (m, -0.0001 -> -3.0)!!!値てきとう
    W_wilt  = 0.16642   !wilting point    (m3/m3, 0.0 -> 1.0) !!!値てきとう
    SoilClass = 1       !Zobler(1986)'s slope type index 1-9 (default = 1) !!!値てきとう
    write(*,*) GRID%N_x,GRID%N_y
    write(*,*) GRID%Max_x,GRID%Max_y
    
   ! Allocate variables for each floor cell
   allocate (gmass_leaf     (GRID%N_x, GRID%N_y)     )
   allocate (gmass_root     (GRID%N_x, GRID%N_y)     )
   allocate (gmass_available(GRID%N_x, GRID%N_y)     )
   allocate (gmass_stock    (GRID%N_x, GRID%N_y)     )
   allocate (lai_grass      (GRID%N_x, GRID%N_y)     )
   allocate (lai_opt_grass_RunningRecord(20, GRID%N_x, GRID%N_y) )
   
   allocate (patch_vacant   (GRID%N_x, GRID%N_y)     )

   allocate (par_floor_rel  (GRID%N_x, GRID%N_y)     )
   allocate (par_grass_rel  (GRID%N_x, GRID%N_y)     )
   allocate (par_floor_dir_rel  (GRID%N_x, GRID%N_y)     )
   allocate (par_floor_dif_rel  (GRID%N_x, GRID%N_y)     )
   allocate (par_grass_dir_rel  (GRID%N_x, GRID%N_y)     )
   allocate (par_grass_dif_rel  (GRID%N_x, GRID%N_y)     )
   allocate (sum_par_floor  (GRID%N_x, GRID%N_y)     )
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:add
   !data processing
   if (W_fi   > W_sat ) W_fi   = W_sat
   if (W_wilt > W_sat ) W_wilt = W_sat
   
!TMP: Soil properties by Hazeltine & Prentice (1996), which is also employed by Sato et al (2010)
!   W_sat  = 0.443
!   W_fi   = 0.399
!   W_wilt = 0.211
   
!TMP: Soil properties of SOILTYP=6 of NOAH-LSM
!   W_sat  = 0.465
!   W_fi   = 0.465
!   W_wilt = 0.103
   
!_____________ Prepare Climate Data
!Count line numbers in climate-data-file
   count = 0
   open (File_no(1), file=Fn_climate, status='OLD')
   do i=1,Day_in_Year * 10000
      count = count + 1
      read(File_no(1),'(a)',end=100) chr_tmp
   enddo
   100 continue
   close (File_no(1))
   YearMaxClimate = int(count/Day_in_Year)
   
!Count item number in each line of the climate-data-file
   chr_tmp = trim(chr_tmp)
   do count=1, 10
      j = SCAN (chr_tmp,',')
      if (j==0) exit
      chr_tmp = chr_tmp((j+1):)
   enddo
   
!Set sizes of allocatable climate data table
   allocate (tmp_air       (Day_in_Year, YearMaxClimate)     )
   allocate (tmp_air_range (Day_in_Year, YearMaxClimate)     )
   allocate (prec          (Day_in_Year, YearMaxClimate)     )
   allocate (rad_short     (Day_in_Year, YearMaxClimate)     )
   allocate (rad_long      (Day_in_Year, YearMaxClimate)     )
   allocate (wind          (Day_in_Year, YearMaxClimate)     )
   allocate (rh            (Day_in_Year, YearMaxClimate)     )
   allocate (tmp_soil1     (Day_in_Year, YearMaxClimate)     )
   allocate (tmp_soil2     (Day_in_Year, YearMaxClimate)     )
   allocate (tmp_soil3     (Day_in_Year, YearMaxClimate)     )
   
   allocate (tmp_soil      (Day_in_Year, YearMaxClimate, NumSoil) )
   
!Read Climatic data
   write(*,*) 'Read Climatic data' !!!<<<<<<<<<<<<TN:add
   open (File_no(1), file=Fn_climate, status='OLD')
   do i=1, YearMaxClimate
   do j=1, Day_in_Year
      if (count==8) then
         !In case soil temperatures data is NOT available
         read (File_no(1),*) tmp_air(j,i), tmp_air_range(j,i), &
                             prec(j,i), rad_short(j,i), rad_long(j,i), wind(j,i), rh(j,i)
         tmp_soil1(j,i) = tmp_air(j,i)
         tmp_soil2(j,i) = tmp_air(j,i)
         tmp_soil3(j,i) = tmp_air(j,i)
      else
         !In case soil temperatures data is available
         read (File_no(1),*) tmp_air(j,i), tmp_air_range(j,i), &
                             prec(j,i), rad_short(j,i), rad_long(j,i), wind(j,i), rh(j,i), &
                             tmp_soil1(j,i), tmp_soil2(j,i), tmp_soil3(j,i)
      endif
      
      !Give adhock values for soil temperature
      Call tmp_soil_interpolate (tmp_soil1(j,i), tmp_soil2(j,i), tmp_soil3(j,i), tmp_soil(j,i,:))
   enddo
   enddo
   close (File_no(1))
   
   !______________ Read time-series of atmospheric CO2
!Scan climate data length
   YearMaxCO2 = 0
   open (File_no(1), file=Fn_CO2, status='OLD')
   do i=1, 10000
      YearMaxCO2 = YearMaxCO2 + 1
      read(File_no(1),*,end=200)
   enddo
   200 continue
   close (File_no(1))
   YearMaxCO2 = YearMaxCO2 - 1
   
!Set sizes of allocatable CO2 data table
   allocate ( aco2_annual(YearMaxCO2) )
   
!Read CO2 data
   Open (File_no(1), file=Fn_CO2, status='OLD')
   do i = 1, YearMaxCO2
      read(File_no(1), *) aco2_annual(i)
   end do
   Close (File_no(1))
   
!___________ For employing different random seed for each run (by Shigeki Ikeda @ Kyoto Univ.)
   IF (Flag_randomization) then
      call random_seed(size=seedsize)
      allocate(seed(seedsize))
      
      do size_count=1,seedsize
      call system_clock(count=clock)
      seed(size_count)=clock
      end do
      
      call random_seed(put=seed)
   EndIf
   
!_____________ Display Simulation Conditions
   !Print location properties
   write (*,*)
   write (*,*) '*********  Coodinate Configurations  *********'
   write (*,*) 'Latitude, Longtitude   :', LAT, LON
   write (*,*) 'Point nomal            :', point
   write (*,*)
   
   !Print soil properties
   write (*,*) '*********  Location properties  *********'
   write (*,*) 'Altitude    :', ALT         
   write (*,*) 'Albedo_soil0:', Albedo_soil0
   write (*,*) 'W_sat       :', W_sat       
   write (*,*) 'W_fi        :', W_fi        
   write (*,*) 'W_wilt      :', W_wilt      
   write (*,*)                              
   
   !Print climate properties
   x = real(Day_in_Year * YearMaxClimate)
   
   write (*,*) '*********  Wether statistics (annual mean of the all years)  *********'
   write (*,*) '2m air temperature  (Cecius)  :', sum(tmp_air  (:,:)) / x
   write (*,*) 'precipitation       (mm/year) :', sum(prec     (:,:)) / real(YearMaxClimate)
   write (*,*) 'Relative humidity   (%)       :', sum(rh       (:,:)) / x
   write (*,*) 'wind                (m/s)     :', sum(wind     (:,:)) / x
   write (*,*) 
   
!_____________ Simulation
   !Call simulation loop
   IF (Mask==0 .or. Albedo_soil0<=0.0 .or. W_sat<=0.0 .or. W_fi<=0.0 .or. W_wilt<=0.0 ) then
      write(*,*) 'Error: Invalid location properties'
   ELSE
      write (*,*) '*********  Now simulating  *********'
      Call main_loop ( &
      LAT, LON, GlobalZone, YearMaxClimate, YearMaxCO2, &
      tmp_air(:,:), tmp_air_range(:,:), prec(:,:), rad_short(:,:), rad_long(:,:), &
      wind(:,:), rh(:,:), tmp_soil(:,:,:), &
      aco2_annual(:), ALT, Albedo_soil0, W_fi, W_wilt, W_sat, W_mat, SoilClass)
      write (*,*) '*********  Done  *********'
   END IF
   
   STOP
   
END PROGRAM start_point

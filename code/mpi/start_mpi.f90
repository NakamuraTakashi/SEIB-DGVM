!*************************************************************************************************
! Start up procedure for multiple grid cells
! using parallel computation with MPI
!*************************************************************************************************
PROGRAM omp_kickoff
   USE mpi
   USE data_structure
   implicit none
   
!_____________ Set Parameters
!Directory of spin up files
   character(len=*),parameter::Loc_spnin_files='/data/hsato/siberia/spnin1850/'
!  character(len=*),parameter::Loc_spnin_files='/data/hsato/siberia/spnin2006/'
 
!Directory of aCO2 data
   character(len=*),parameter::Loc_CO2_files='co2_1850_2100_rcp85.dat' !RCP8.5
!  character(len=*),parameter::Loc_CO2_files='co2_1850_2100_rcp26.dat' !RCP2.6
 
!Directory of climate data for future projection
 character(len=*),parameter::Loc_climate_data2='/data/climate_MirocAR5Base_daily_0.5deg_V3/RCP85/'
!character(len=*),parameter::Loc_climate_data2='/data/climate_MirocAR5Base_daily_0.5deg_V3/RCP26/'
 
!Directory of climate data for histrical reconstruction
 character(len=*),parameter::Loc_climate_data1='/data/climate_MirocAR5Base_daily_0.5deg_V3/historical/'
   
!Set year number of climate & CO2 data
   integer,parameter::YearMaxClimate1 = 156  !historical(1850-2005年)
   integer,parameter::YearMaxClimate2 =  95  !RCPs (2006-2100年)
   integer,parameter::YearMaxCO2      = 251  !year length of the CO2 data (1850~2100)
   
!Land mask file
   character(len=*),parameter::Fn_landmask = &
                               '/data/landmask_0.5deg.txt'
   
!Directory for writing an output file of each simulation grid
   character(len=*),parameter::Loc_result_files = &
                               '/data/hsato/siberia4/result/'
   
!Directory for writing analysis data for whole simulation grid
   character(len=*),parameter::Loc_analysis_files= &
                               '/home/hsato/work_distribute/result_visualize/'
   
!Location of soil property data
!A same data-set is employed at every simulation grids
!This data-set was obtained from GSWP2 (http://www.iges.org/gswp2/)
   character(len=*),parameter::Fn_landprop = 'land_prop.txt'
   
!Maximum grid number for longitude and latitude
  integer,parameter::LatMax = 360 !Grid number for latitude  (@ 0.5deg Grid System)
  integer,parameter::LonMax = 720 !Grid number for longitude (@ 0.5deg Grid System)
   
!Grid length for simulation
  !(East Siberia @ 0.5deg grid mesh)
  integer,parameter::LatNoStart =  20 !(N10)
  integer,parameter::LatNoEnd   = 100 !(N40)
  integer,parameter::LonNoStart = 521 
  integer,parameter::LonNoEnd   = 720 
  
!_____________ Set Variables
   !MPI control variables
   integer myid, numprocs, ierr
   
   !Coodinate variables for parallel computation
   integer gridNo, gridNoMax        !Sequential number for simulation grid cell
   integer latNo, lonNo             !Sequential number for north-south and west-east
   real    lat,lon                  !Latitude and Longitude
   integer,dimension(LatMax*LonMax)::latNo_save, lonNo_save
   
   !For preparing Land-Ocean mask
   logical,dimension(LonMax,LatMax)::Landmask     !landmask; true->land, false->water
   integer,dimension(1:LonMax)     ::readerLonMax !For reading land mask
   integer,dimension(360*180)::Mask, ALT, Albedo_soil0, W_sat, W_fi, W_mat, W_wilt, SoilClass
   integer point
   
   !Atomospheric CO2 time-series @ ppm
   real,dimension(251)::aco2_1850to2100
   
   !Counter
   integer i
   
!_____________ Read Parameters
!Read Parameter files
   open (1, file='parameter.txt', action='READ', status='OLD')
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
   
!______________ Prepare landmask (extension for wide area simulation)
!Read Landmask
   landmask(:,:) = .false.
   open (1, file=Fn_landmask, status='OLD')
   do latNo=1, LatMax
      read(1,*) readerLonMax(1:LonMax)
      do lonNo=1, LonMax
         if (readerLonMax(lonNo)==1) landmask(lonNo,latNo)=.true.
      end do
   end do
   close(1)
   
   Mask        (:)=0
   ALT         (:)=0.0
   Albedo_soil0(:)=0.0
   W_sat       (:)=0.0
   W_fi        (:)=0.0
   W_mat       (:)=0.0
   W_wilt      (:)=0.0
   SoilClass   (:)=0
   open (1, file='land_prop.txt', status='OLD')
   do i=1, 180*360
      read(1,*) Mask(i), ALT(i), Albedo_soil0(i), &
                W_sat(i), W_fi(i), W_mat(i), W_wilt(i), SoilClass(i)
   end do
   close(1)
   
   !Landmask correction
   do latNo=1, LatMax
   do lonNo=1, LonMax
      if (.not. landmask(lonNo,latNo)) cycle
      
      lat   =    90.0 - (real(latNo)-0.5) * 0.5
      lon   = - 180.0 + (real(lonNo)-0.5) * 0.5
      point = (90-int(lat)-1)*360 + int(lon+180) + 1 !grid point number @1.0 deg system
      
      if (Mask     (point)==0) landmask(lonNo,latNo)=.false.
      if (SoilClass(point)==0) landmask(lonNo,latNo)=.false.
      if (SoilClass(point)==9) landmask(lonNo,latNo)=.false.
   end do
   end do
   
!______________ Read time-series of atmospheric CO2 concentration
   Open (1, file=trim(Loc_CO2_files), status='OLD')
   do i = 1, 251 !1850~2100
      read(1, *) aco2_1850to2100(i)
   end do
   Close (1)
   
!______________ Provide reference IDs for each grid cell
   gridNo = 0
   do latNo = LatNoStart, LatNoEnd
   do lonNo = LonNoStart, LonNoEnd
      if (.not. landmask(lonNo,latNo)) cycle
      
      gridNo             = gridNo + 1
      gridNoMax          = gridNo
      latNo_save(gridNo) = latNo
      lonNo_save(gridNo) = lonNo
   end do
   end do
   if (gridNoMax==0) Stop
   
!______________ Initialization of MPI
   Call MPI_INIT( ierr )
   Call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
   Call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
   
!______________ Conduct parallel computation with MPI
   if (myid==1) write(*,*) 'Start parallel computation....'
   
   DO gridNo = 1, gridNoMax
   if ( myid == mod(gridNo-1, numprocs) ) then
      Call start(myid, latNo_save(gridNo), lonNo_save(gridNo), aco2_1850to2100, &
                 YearMaxClimate1, YearMaxClimate2, YearMaxCO2, &
                 Loc_climate_data1, Loc_climate_data2, Loc_result_files, Loc_spnin_files)
   endif
   END DO
   
!______________ Termination procudure of MPI compuation
   Write(*,*) 'Terminate processor:', myid, ', Total processor:', numprocs
   Call MPI_Barrier(MPI_COMM_WORLD, ierr) !Wait until end of all simulation
   Call MPI_FINALIZE(ierr)                !Finalize procudure for MPI
   
!_____________ Conversion of output files for each simulation grids into maps
!グリッド別のファイル出力を、ここで地図に加工する
   if (myid==0) then
      write(*,*) 'Converting result files....'
      !Analysis output values of Last 10 yers
      Call after_sim1(LatMax, LonMax, LatNoStart, LatNoEnd, LonNoStart, LonNoEnd, &
                      Loc_result_files, Loc_analysis_files, landmask)
   endif
   
STOP
END PROGRAM omp_kickoff



!*************************************************************************************************
! Start up procedure for each simulation grid
!*************************************************************************************************
Subroutine start (myid, latNo, lonNo, aco2_1850to2100, &
                  YearMaxClimate1, YearMaxClimate2, YearMaxCO2, &
                  Loc_climate_data1, Loc_climate_data2, &
                  Loc_result_files, Loc_spnin_files)
   
   USE data_structure
   implicit none
   
!_____________ Set Augment
   integer            ,intent(IN):: myid, latNo, lonNo
   integer            ,intent(IN):: YearMaxClimate1, YearMaxClimate2, YearMaxCO2
   character(len=*)   ,intent(IN):: Loc_climate_data1, Loc_climate_data2
   character(len=*)   ,intent(IN):: Loc_result_files, Loc_spnin_files
   real,dimension(YearMaxCO2),intent(IN):: aco2_1850to2100 !Atomospheric co2 concentration (ppm)
   
!_____________ Set Variables
!Climate data
   real,allocatable,dimension(:,:)  ::&
    tmp_air       ,& !1. Surface air temperature (Celcius)
    tmp_air_range ,& !2. Daily range of tmp_air (Celcius)
    prec          ,& !3. Precipitation (mm day-1)
    rad_short     ,& !4. Shortwave radiation, downward @ midday (W m-2)
    rad_long      ,& !5. Daily mean of longwave radiation, downward (W m-2)
    wind          ,& !6. Wind velocity (m s-1)
    rh               !7. Relative humidity (%)
   
   real,allocatable,dimension(:,:,:)::&
    tmp_soil      !Soil temperature for each layers (Celcius)
   
!Location data
   integer::Mask         !Land ocean mask at 1.0 deg (1:land, 0:ocean)
   real   ::ALT          !altitude (m above MSL)
   real   ::Albedo_soil0 !albedo, default
   real   ::W_fi         !filed capacity   (m3/m3, 0.0 -> 1.0)
   real   ::W_wilt       !wilting point    (m3/m3, 0.0 -> 1.0)
   real   ::W_sat        !saturate point   (m3/m3, 0.0 -> 1.0)
   real   ::W_mat        !matrix potential (m, -0.0001 -> -3.0)
   integer::SoilClass    !
   
!Others
   real    LAT, LON              !latitude and logitude for simulate
   integer GlobalZone            !ID number of global zone
   integer year, doy, dat, point !Counters
   integer i, j                  !for general usage
   
!_____________ Set Variables, for wide area computations
!For reading data
   real,dimension(Day_in_Year, YearMaxClimate1, 1:8)::dataREAD1  !1850~2005
   real,dimension(Day_in_Year, YearMaxClimate2, 1:8)::dataREAD2  !2006~2100
   
!For I/O
   character(len= 3) nam_lat
   character(len= 3) nam_lon
   character(len= 2) nam_dat
   integer i1, i2, i3
   integer file_no_grid1, file_no_grid2, file_no_grid3
   
!______________________Set varieties of reference number for this grid cell
!Device number for I/O
   file_no_grid1 = myid + 100 !For output files1, and for general usage
   file_no_grid2 = myid + 200 !For reading spinup files
   file_no_grid3 = myid + 300 !For writing spinup files
   
!Reference number for location
   ! LAT: north +, south - (decimalized)
   ! LON: east  +, west  - (decimalized)
   LAT   =    90.0 - (real(latNo)-0.5) * 0.5
   LON   = - 180.0 + (real(lonNo)-0.5) * 0.5
   point = (90-int(LAT)-1)*360 + int(LON+180) + 1 !grid point number @1.0 deg system
   
!Characters for designating latitude and longitude
   i1 = int(  latNo               /100 ) ; nam_lat(1:1) = char(i1+48)
   i2 = int( (latNo -i1*100)      /10  ) ; nam_lat(2:2) = char(i2+48)
   i3 =    (  latNo -i1*100-i2*10      ) ; nam_lat(3:3) = char(i3+48)
   
   i1 = int(  lonNo               /100 ) ; nam_lon(1:1) = char(i1+48)
   i2 = int( (lonNo -i1*100)      /10  ) ; nam_lon(2:2) = char(i2+48)
   i3 =    (  lonNo -i1*100-i2*10      ) ; nam_lon(3:3) = char(i3+48)
   
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
   
!______________ Read Location data
   open (file_no_grid1, file='land_prop.txt', status='OLD')
   do i=1, point
      read(file_no_grid1,*) Mask, ALT, Albedo_soil0, W_sat, W_fi, W_mat, W_wilt, SoilClass
   end do
   close(file_no_grid1)
   if (W_fi   > W_sat ) W_fi   = W_sat
   if (W_wilt > W_sat ) W_wilt = W_sat
   
!_____________ Prepare Climate Data
!Set sizes of allocatable climate data table
   allocate (tmp_air       (Day_in_Year, YearMaxClimate1+YearMaxClimate2)          ) !1
   allocate (tmp_air_range (Day_in_Year, YearMaxClimate1+YearMaxClimate2)          ) !2
   allocate (prec          (Day_in_Year, YearMaxClimate1+YearMaxClimate2)          ) !3
   allocate (rad_short     (Day_in_Year, YearMaxClimate1+YearMaxClimate2)          ) !4
   allocate (rad_long      (Day_in_Year, YearMaxClimate1+YearMaxClimate2)          ) !5
   allocate (wind          (Day_in_Year, YearMaxClimate1+YearMaxClimate2)          ) !6
   allocate (rh            (Day_in_Year, YearMaxClimate1+YearMaxClimate2)          ) !7
   allocate (tmp_soil      (Day_in_Year, YearMaxClimate1+YearMaxClimate2, NumSoil) ) ! 
   
   !Read climate data, 1850~2005
   i = 7 * 365 * YearMaxClimate1 * 4 !file size for each variables = 4byte
   Open (file_no_grid1, file=Loc_climate_data1//nam_lat//'/'//nam_lon//'.dat', &
         access='direct', recl=i, status='old')
      read (file_no_grid1,rec=1) &
      (((dataREAD1(doy,year,dat), dat=1,7), doy=1,365), year=1, YearMaxClimate1)
   Close(file_no_grid1)
   
   !Read climate data, 2006~2100
   i = 7 * 365 * YearMaxClimate2  * 4 !データサイズ * 4byte
   Open (file_no_grid1, file=Loc_climate_data2//nam_lat//'/'//nam_lon//'.dat', &
         access='direct', recl=i, status='old')
      read (file_no_grid1,rec=1) &
      (((dataREAD2(doy,year,dat), dat=1,7), doy=1,365), year=1, YearMaxClimate2)
   Close(file_no_grid1)
   
!_____________ Climate data Conversion
   do year=1,YearMaxClimate1
   do doy=1,Day_in_Year
      tmp_air      (doy,year)   = dataREAD1(doy,year,1)
      tmp_soil     (doy,year,:) = dataREAD1(doy,year,1) !substitute air-temp.
      tmp_air_range(doy,year)   = dataREAD1(doy,year,2)
      prec         (doy,year)   = dataREAD1(doy,year,3)
      rad_short    (doy,year)   = dataREAD1(doy,year,4)
      rad_long     (doy,year)   = dataREAD1(doy,year,5)
      wind         (doy,year)   = dataREAD1(doy,year,6)
      rh           (doy,year)   = dataREAD1(doy,year,7)
   enddo
   enddo
   
   do year = YearMaxClimate1+1, YearMaxClimate1+YearMaxClimate2
   do doy=1,Day_in_Year
      tmp_air      (doy, year)   = dataREAD2(doy,year-YearMaxClimate1,1)
      tmp_soil     (doy, year,:) = dataREAD2(doy,year-YearMaxClimate1,1) !substitute air-temp.
      tmp_air_range(doy, year)   = dataREAD2(doy,year-YearMaxClimate1,2)
      prec         (doy, year)   = dataREAD2(doy,year-YearMaxClimate1,3)
      rad_short    (doy, year)   = dataREAD2(doy,year-YearMaxClimate1,4)
      rad_long     (doy, year)   = dataREAD2(doy,year-YearMaxClimate1,5)
      wind         (doy, year)   = dataREAD2(doy,year-YearMaxClimate1,6)
      rh           (doy, year)   = dataREAD2(doy,year-YearMaxClimate1,7)
   enddo
   enddo
   
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
   
!_____________ !Call simulation loop for each grid
!Open I/O files
   !for writing output files
   if (Flag_output_write) then
      open (file_no_grid1, &
      file = Loc_result_files//nam_lat//'_'//nam_lon//'_out.txt')
   endif
   
   !for reading spinup files
   if (Flag_spinup_read) then
     open ( file_no_grid2, &
     file = Loc_spnin_files//nam_lat//'_'//nam_lon//'_spnin.dat', &
     access="sequential", form="unformatted", status="old", action="read")
   endif
   
   !for writing spinup files
   if (Flag_spinup_write) then
      open ( file_no_grid3, &
      file = Loc_result_files//nam_lat//'_'//nam_lon//'_spnout.dat', &
      access="sequential", form="unformatted", action="write")
   endif
   
!Call main program
   Call main_loop ( &
   LAT, LON, GlobalZone, YearMaxClimate1+YearMaxClimate2, YearMaxCO2, &
   file_no_grid1, file_no_grid2, file_no_grid3, &
   tmp_air(:,:), tmp_air_range(:,:), prec(:,:), rad_short(:,:), rad_long(:,:), &
   wind(:,:), rh(:,:), tmp_soil(:,:,:), &
   aco2_1850to2100(:), ALT, Albedo_soil0, W_fi, W_wilt, W_sat, W_mat, SoilClass)
   
!Close I/O files
   if (Flag_output_write) close ( file_no_grid1 ) !取りまとめ地図用の出力
   if (Flag_spinup_read ) close ( file_no_grid2 ) !for reading spinup files
   if (Flag_spinup_write) close ( file_no_grid3 ) !for writing spinup files
   
END Subroutine start



!*************************************************************************************************
! Analysis output values of Last 10 yers
!*************************************************************************************************
Subroutine after_sim1(LatMax, LonMax, LatNoStart, LatNoEnd, LonNoStart, LonNoEnd, &
                      Loc_result_files, Loc_analysis_files, landmask)
   
   implicit none
   
!_____________ Set Augment
   integer         ,intent(IN):: LatMax, LonMax, LatNoStart, LatNoEnd, LonNoStart, LonNoEnd
   character(len=*),intent(IN):: Loc_result_files, Loc_analysis_files
   
   logical,dimension(LonMax, LatMax),intent(IN):: landmask !Land Ocean mask
   
!_____________ Set parameters
   !Period for writing monthly outputs @ year
   !Monthly data is written for last YearForMean years of the simulation
   integer,parameter::YearForMean = 10 !This value has to be same as one in the main_mpi.f90
   
!_____________ Set variables
!Array for inputing result files
   real,dimension(1:LatMax, 1:LonMax):: Larch_area !Larch Area (fraction)
   
   integer,dimension(                                               2)::data1_read
   integer,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd, 12,  2)::data1
   
   real   ,dimension(                                              28)::data2_read
   real   ,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd, 12, 28)::data2
   
!Variables for outputting geographic distribution
   real,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd):: &
      out_water        , & !Water content @ top soil layer, annual average
      out_precipitation, & !Annual precip.
      out_gpp          , & !GPP (kg C/ m2/ year)
      out_npp          , & !NPP (kg C/ m2/ year)
      out_nep          , & !NEP (kg C/ m2/ year)
      out_csoil        , & !
      out_hr           , & !HR  (kg C/ m2/ year)
      out_lai_amean    , & !Annual mean of LAI (All   PFTs)
      out_lai_amean_t  , & !Annual mean of LAI (Woody PFTs)
      out_lai_amean_g  , & !Annual mean of LAI (Grass PFTs)
      out_lai_max      , & !Annula maximum of LAI (All   PFTs)
      out_lai_max_t    , & !Annula maximum of LAI (Woody PFTs)
      out_lai_max_g    , & !Annula maximum of LAI (Grass PFTs)
      out_ald_max      , & !Annula Maximum of Active Layer Depth (m)
      out_water_JJA        !Available soil water at top 5 layers during JJA(mm)
   
   !For summary statistics
   real sum_wbiomass  ! 1, Carbon in Woody biomass [kg C / m2]
   real sum_gbiomass  ! 2, Carbon in Grass biomass [kg C / m2]
   real sum_litter    ! 3, Carbon in litter        [kg C / m2]
   real sum_som_int   ! 4, Carbon in som_int       [kg C / m2]
   real sum_som_slow  ! 5, Carbon in som_slow      [kg C / m2]
   real sum_gpp       ! 9, GPP                     [kg C / m2 / month]
   real sum_npp       !10, NPP                     [kg C / m2 / month]
   real sum_nep       !11, NEP                     [kg C / m2 / month]
   real sum_hr        !12, Heterotrophic resp.     [kg C / m2 / month]
   real sum_runoff    !15, runoff                  [mm/month]
   real sum_intercept !16, interception            [mm/month]
   real sum_evapor    !17, evaporation             [mm/month]
   real sum_transpi   !18, transpiration           [mm/month]
   
   !etc
   character(len=7) fname            !For file name
   integer i, i1, i2, i3, i4, i5, i6    !For processing file name
   integer lat, lon, month, year     !Loop counter
   real    x, y, sum_weight
   real    a1, a2, a3, a4, a5, a6
   
!_______________ Read result files and make average value
write(*,*) "Reading result files and making output variables"
   !Initialize
   data1(:,:,:,:) =   0
   data2(:,:,:,:) = 0.0
   
   !Loop for each simulation grid
   DO lat = LatNoStart, LatNoEnd
   DO lon = LonNoStart, LonNoEnd
   IF (landmask(lon,lat) .eqv. .false.) cycle
      
      !Make letter string for designating a file
      i1 = int(  lat              /100 ) ; fname(1:1)=char(i1+48)
      i2 = int( (lat-i1*100)      / 10 ) ; fname(2:2)=char(i2+48)
      i3 =    (  lat-i1*100-i2*10      ) ; fname(3:3)=char(i3+48)
                                           fname(4:4)='_'
      i4 = int(  lon              /100 ) ; fname(5:5)=char(i4+48)
      i5 = int( (lon-i4*100)      / 10 ) ; fname(6:6)=char(i5+48)
      i6 =    (  lon-i4*100-i5*10      ) ; fname(7:7)=char(i6+48)
      
      !Read output files and make average values
      open (1, file=Loc_result_files//fname//'_out.txt', status='OLD')
         
         !Sumup all values for YearForMean year (except for Biome code)
         do year  =1, YearForMean
            do month =1, 12
               read(1,*) data1_read(:), data2_read(:)
               data1(lat,lon,month,1) =                          data1_read(1) !Biome code
               data1(lat,lon,month,2) = data1(lat,lon,month,2) + data1_read(2) !Drought days
               data2(lat,lon,month,:) = data2(lat,lon,month,:) + data2_read(:) !Other Variables
            enddo
         enddo
         
         !Make average values by dividing YearForMean (except for Biome code)
         do month =1, 12
            data1(lat,lon,month,2) = data1(lat,lon,month,2) / real(YearForMean)
            data2(lat,lon,month,:) = data2(lat,lon,month,:) / real(YearForMean)
         enddo
         
      close (1)
      
   END DO
   END DO
   
!_______________ Prepare output variables
write(*,*) "Prepare output variables"
   !initialize
   out_water         (:,:) = 0.0
   out_precipitation (:,:) = 0.0
   out_csoil         (:,:) = 0.0
   
   out_gpp           (:,:) = 0.0
   out_npp           (:,:) = 0.0
   out_nep           (:,:) = 0.0
   out_hr            (:,:) = 0.0
   
   out_lai_amean     (:,:) = 0.0
   out_lai_amean_t   (:,:) = 0.0
   out_lai_amean_g   (:,:) = 0.0
   
   out_lai_max       (:,:) = 0.0
   out_lai_max_t     (:,:) = 0.0
   out_lai_max_g     (:,:) = 0.0
   
   out_ald_max       (:,:) = 0.0
   out_water_JJA     (:,:) = 0.0
   
   Do lat   = LatNoStart, LatNoEnd
   Do lon   = LonNoStart, LonNoEnd
      
      Do month = 1, 12
      !Soil water contents
      out_water(lat,lon)= out_water(lat,lon)+ data2(lat,lon,month, 6) / 12.0
      
      !Annual GPP
      out_gpp  (lat,lon)= out_gpp  (lat,lon)+ data2(lat,lon,month, 9)
      
      !Annual NPP
      out_npp  (lat,lon)= out_npp  (lat,lon)+ data2(lat,lon,month,10)
      
      !Annual NEP
      out_nep  (lat,lon)= out_nep  (lat,lon)+ data2(lat,lon,month,11)
      
      !Annual Heterotrophic resp.
      out_hr   (lat,lon)= out_hr   (lat,lon)+ data2(lat,lon,month,25)
      
      !Annual means of LAI
      out_lai_amean  (lat,lon) = out_lai_amean  (lat,lon) + data2(lat,lon,month,13) / 12.0 &
                                                          + data2(lat,lon,month,14) / 12.0
      out_lai_amean_t(lat,lon) = out_lai_amean_t(lat,lon) + data2(lat,lon,month,13) / 12.0
      out_lai_amean_g(lat,lon) = out_lai_amean_g(lat,lon) + data2(lat,lon,month,14) / 12.0
      
      !Annual maximums of LAI
      out_lai_max  (lat,lon) = max( out_lai_max  (lat,lon), data2(lat,lon,month,13) &
                                                          + data2(lat,lon,month,14) )
      out_lai_max_t(lat,lon) = max( out_lai_max_t(lat,lon), data2(lat,lon,month,13) )
      out_lai_max_g(lat,lon) = max( out_lai_max_g(lat,lon), data2(lat,lon,month,14) )
      
      !Annual Precipitation
      out_precipitation(lat,lon) = out_precipitation(lat,lon)+ sum(data2(lat,lon,month,15:18))
      
      !Hydorology Related
      out_ald_max   (lat,lon) = max ( out_ald_max(lat,lon), data2(lat,lon,month,27) )
      
      End Do
      
      !Hydorology Related
      Do month = 6, 8
      out_water_JJA (lat,lon) = out_water_JJA (lat,lon) + data2(lat,lon,month,28) / 3.0
      End Do
      
   End Do
   End Do
   
!_______________ Write output files 1
write(*,*) "Writing result maps 1"
   Open ( 1, file=Loc_analysis_files//'out_biome.txt'        )
   Open ( 2, file=Loc_analysis_files//'out_wbiomass.txt'     )
   Open ( 3, file=Loc_analysis_files//'out_water1.txt'       )
   Open ( 4, file=Loc_analysis_files//'out_gpp.txt'          )
   Open ( 5, file=Loc_analysis_files//'out_npp.txt'          )
   Open ( 6, file=Loc_analysis_files//'out_nep.txt'          )
   Open ( 7, file=Loc_analysis_files//'out_hr.txt'           )
   
   Open ( 8, file=Loc_analysis_files//'out_precipitation.txt')
   Open ( 9, file=Loc_analysis_files//'out_fire.txt'         )
   
   Open (10, file=Loc_analysis_files//'out_lai_max.txt'      )
   Open (11, file=Loc_analysis_files//'out_lai_max_t.txt'    )
   Open (12, file=Loc_analysis_files//'out_lai_max_g.txt'    )
   
   Open (20, file=Loc_analysis_files//'out_lai_amean.txt'    )
   Open (21, file=Loc_analysis_files//'out_lai_amean_t.txt'  )
   Open (22, file=Loc_analysis_files//'out_lai_amean_g.txt'  )
   
   Open (23, file=Loc_analysis_files//'out_ald_max.txt'      )
   Open (24, file=Loc_analysis_files//'out_water_JJA.txt'    )
   
   Open (25, file=Loc_analysis_files//'out_csoil.txt'        )
   
   Do lat = LatNoStart, LatNoEnd
      Do lon = LonNoStart, LonNoEnd
      write( 1,'(  i2,a)', advance='no') data1     (lat,lon,12, 1),  ',' !Biome
      write( 2,'(f7.2,a)', advance='no') data2     (lat,lon,12, 1),  ',' !Woody biomass
      write( 3,'(f9.1,a)', advance='no') out_water (lat,lon)      ,  ',' !Water content @ top soil layer, annual average
      write( 4,'(f8.3,a)', advance='no') out_gpp          (lat,lon), ',' !GPP (kg C/ m2/ year)
      write( 5,'(f8.3,a)', advance='no') out_npp          (lat,lon), ',' !NPP (kg C/ m2/ year)
      write( 6,'(f8.3,a)', advance='no') out_nep          (lat,lon), ',' !NEP (kg C/ m2/ year)
      write( 7,'(f8.3,a)', advance='no') out_hr           (lat,lon), ',' !HR  (kg C/ m2/ year)
      
      write( 8,'(f9.1,a)', advance='no') out_precipitation(lat,lon      ), ',' !Annual precip.
      write( 9,'(f5.3,a)', advance='no') data2            (lat,lon,12,21), ',' !Fire frequency
      
      !Annual maximum LAI
      write(10,'(f4.1,a)', advance='no') out_lai_max      (lat,lon), ',' !All PFTs
      write(11,'(f4.1,a)', advance='no') out_lai_max_t    (lat,lon), ',' !Woody PFTs
      write(12,'(f4.1,a)', advance='no') out_lai_max_g    (lat,lon), ',' !Grass PFTs
      
      !Annual average LAI
      write(20,'(f4.1,a)', advance='no') out_lai_amean    (lat,lon), ',' !All PFTs
      write(21,'(f4.1,a)', advance='no') out_lai_amean_t  (lat,lon), ',' !Woody PFTs
      write(22,'(f4.1,a)', advance='no') out_lai_amean_g  (lat,lon), ',' !Grass PFTs
      
      !Annual maximum of activeaverage LAI
      write(23,'(f5.3,a)', advance='no') out_ald_max      (lat,lon), ',' !Active Layer Maximum
      
      !Annual maximum of activeaverage LAI
      write(24,'(f5.1,a)', advance='no') out_water_JJA    (lat,lon), ',' !JJA available water [mm]
      
      !Soil Carbon
      write(25,'(f9.2,a)', advance='no') sum(data2(lat,lon,12,3:5)),  ',' !Soil Carbon [kg C / m2]
      
      End Do
      
      !Insert feed code
      write( 1,*); write( 2,*); write( 3,*); write( 4,*); write( 5,*); write( 6,*); write( 7,*)
      write( 8,*); write( 9,*)
      write(10,*); write(11,*); write(12,*)
      write(20,*); write(21,*); write(22,*)
      write(23,*); write(24,*); write(25,*)
      
   End Do
   Close ( 1); Close ( 2); Close ( 3); Close ( 4); Close ( 5); Close ( 6); Close ( 7)
   Close ( 8); Close ( 9)
   Close (10); Close (11); Close (12)
   Close (20); Close (21); Close (22)
   Close (23); Close (24); Close (25)
   
!!_______________ Output of latitudinal averages within a rectangular area
!   Open (1, file=Loc_analysis_files//'SEIB_LatGradient_2005_NoFreezingNoBottomleak.txt')
!  !Open (1, file=Loc_analysis_files//'SEIB_LatGradient_2005_noPermadrost.txt')
!  !Open (1, file=Loc_analysis_files//'SEIB_LatGradient_2100.txt') !
!      Do lat = 36, 100    
!         i=0 ; a1=0.0 ; a2=0.0 ; a3=0.0 ; a4=0.0 ; a5=0.0 ; a6=0.0
!         Do lon = 521,690 
!            if (landmask(lon,lat)==.false.) cycle
!            i  =  i + 1
!            
!            a1 = a1 + out_gpp       (lat,lon)       !GPP                [kgC/m2/yr]
!            a2 = a2 + out_npp       (lat,lon)       !NPP                [kgC/m2/yr]
!           if (out_ald_max(lat,lon)<2.0) then
!            a3 = a3 + 1.0                           !永久凍土の存在比率 [frac]
!           endif
!            a4 = a4 + out_water_JJA (lat,lon)/500.0 !0~50cm含水率@JJA   [frac]
!            a5 = a5 + sum(data2(lat,lon,12, 1:2))   !Biomass            [kgC/m2]
!            a6 = a6 + data2(lat,lon, 7,13)          !July LAI of trees  [m2/m2]
!         Enddo
!         
!         If (i>=1) then
!            a1=a1/real(i);a2=a2/real(i);a3=a3/real(i);a4=a4/real(i);a5=a5/real(i);a6=a6/real(i)
!         Else
!            a1=0.0;a2=0.0;a3=0.0;a4=0.0;a5=0.0;a6=0.0
!         Endif
!         
!         write(1, '(f8.3,a, f8.3,a, f5.3,a, f5.3,a, f5.2,a, f4.2)') &
!                 a1,',',a2,',',a3,',',a4,',',a5,',',a6
!      Enddo
!   Close (1)
!   
!!_______________ Write output files 2 (Monthly LAI of woody PFTs)
!write(*,*) "Writing result maps 2 (Monthly LAI of woody PFTs)"
!   Open ( 1, file=Loc_analysis_files//'out_lai_wood_01.txt')
!   Open ( 2, file=Loc_analysis_files//'out_lai_wood_02.txt')
!   Open ( 3, file=Loc_analysis_files//'out_lai_wood_03.txt')
!   Open ( 4, file=Loc_analysis_files//'out_lai_wood_04.txt')
!   Open ( 5, file=Loc_analysis_files//'out_lai_wood_05.txt')
!   Open ( 6, file=Loc_analysis_files//'out_lai_wood_06.txt')
!   Open ( 7, file=Loc_analysis_files//'out_lai_wood_07.txt')
!   Open ( 8, file=Loc_analysis_files//'out_lai_wood_08.txt')
!   Open ( 9, file=Loc_analysis_files//'out_lai_wood_09.txt')
!   Open (10, file=Loc_analysis_files//'out_lai_wood_10.txt')
!   Open (11, file=Loc_analysis_files//'out_lai_wood_11.txt')
!   Open (12, file=Loc_analysis_files//'out_lai_wood_12.txt')
!   
!   Do lat = LatNoStart, LatNoEnd
!      Do lon = LonNoStart, LonNoEnd
!      write( 1,'(f4.1,a)', advance='no') data2(lat,lon, 1,13), ','
!      write( 2,'(f4.1,a)', advance='no') data2(lat,lon, 2,13), ','
!      write( 3,'(f4.1,a)', advance='no') data2(lat,lon, 3,13), ','
!      write( 4,'(f4.1,a)', advance='no') data2(lat,lon, 4,13), ','
!      write( 5,'(f4.1,a)', advance='no') data2(lat,lon, 5,13), ','
!      write( 6,'(f4.1,a)', advance='no') data2(lat,lon, 6,13), ','
!      write( 7,'(f4.1,a)', advance='no') data2(lat,lon, 7,13), ','
!      write( 8,'(f4.1,a)', advance='no') data2(lat,lon, 8,13), ','
!      write( 9,'(f4.1,a)', advance='no') data2(lat,lon, 9,13), ','
!      write(10,'(f4.1,a)', advance='no') data2(lat,lon,10,13), ','
!      write(11,'(f4.1,a)', advance='no') data2(lat,lon,11,13), ','
!      write(12,'(f4.1,a)', advance='no') data2(lat,lon,12,13), ','
!      End Do
!      
!      write( 1, *); write( 2, *); write( 3, *); write( 4, *); write( 5, *); write( 6, *)
!      write( 7, *); write( 8, *); write( 9, *); write(10, *); write(11, *); write(12, *)
!   End Do
!   Close ( 1); Close ( 2); Close ( 3); Close ( 4); Close ( 5); Close ( 6)
!   Close ( 7); Close ( 8); Close ( 9); Close (10); Close (11); Close (12)
!   
!!_______________ Write output files 3 (Monthly LAI of all PFTs)
!write(*,*) "Writing result maps 2 (Monthly LAI of all PFTs)"
!   Open ( 1, file=Loc_analysis_files//'out_lai_01.txt')
!   Open ( 2, file=Loc_analysis_files//'out_lai_02.txt')
!   Open ( 3, file=Loc_analysis_files//'out_lai_03.txt')
!   Open ( 4, file=Loc_analysis_files//'out_lai_04.txt')
!   Open ( 5, file=Loc_analysis_files//'out_lai_05.txt')
!   Open ( 6, file=Loc_analysis_files//'out_lai_06.txt')
!   Open ( 7, file=Loc_analysis_files//'out_lai_07.txt')
!   Open ( 8, file=Loc_analysis_files//'out_lai_08.txt')
!   Open ( 9, file=Loc_analysis_files//'out_lai_09.txt')
!   Open (10, file=Loc_analysis_files//'out_lai_10.txt')
!   Open (11, file=Loc_analysis_files//'out_lai_11.txt')
!   Open (12, file=Loc_analysis_files//'out_lai_12.txt')
!   
!   Do lat = LatNoStart, LatNoEnd
!      Do lon = LonNoStart, LonNoEnd
!      write( 1,'(f4.1,a)', advance='no') sum(data2(lat,lon, 1,13:14)), ','
!      write( 2,'(f4.1,a)', advance='no') sum(data2(lat,lon, 2,13:14)), ','
!      write( 3,'(f4.1,a)', advance='no') sum(data2(lat,lon, 3,13:14)), ','
!      write( 4,'(f4.1,a)', advance='no') sum(data2(lat,lon, 4,13:14)), ','
!      write( 5,'(f4.1,a)', advance='no') sum(data2(lat,lon, 5,13:14)), ','
!      write( 6,'(f4.1,a)', advance='no') sum(data2(lat,lon, 6,13:14)), ','
!      write( 7,'(f4.1,a)', advance='no') sum(data2(lat,lon, 7,13:14)), ','
!      write( 8,'(f4.1,a)', advance='no') sum(data2(lat,lon, 8,13:14)), ','
!      write( 9,'(f4.1,a)', advance='no') sum(data2(lat,lon, 9,13:14)), ','
!      write(10,'(f4.1,a)', advance='no') sum(data2(lat,lon,10,13:14)), ','
!      write(11,'(f4.1,a)', advance='no') sum(data2(lat,lon,11,13:14)), ','
!      write(12,'(f4.1,a)', advance='no') sum(data2(lat,lon,12,13:14)), ','
!      End Do
!      
!      !Insert feeding code
!      write( 1, *); write( 2, *); write( 3, *); write( 4, *); write( 5, *); write( 6, *)
!      write( 7, *); write( 8, *); write( 9, *); write(10, *); write(11, *); write(12, *)
!   End Do
!   Close ( 1); Close ( 2); Close ( 3); Close ( 4); Close ( 5); Close ( 6)
!   Close ( 7); Close ( 8); Close ( 9); Close (10); Close (11); Close (12)
!   
!!_______________ Write output files 7 (Summup for whole simulation grids)
!write(*,*) "Phase 9"
!   !Initialization
!   sum_wbiomass  = 0.0
!   sum_gbiomass  = 0.0
!   sum_litter    = 0.0
!   sum_som_int   = 0.0
!   sum_som_slow  = 0.0
!   sum_gpp       = 0.0
!   sum_npp       = 0.0
!   sum_nep       = 0.0
!   sum_hr        = 0.0
!   sum_runoff    = 0.0
!   sum_intercept = 0.0
!   sum_evapor    = 0.0
!   sum_transpi   = 0.0
!   
!   sum_weight = 0
!   
!   !Sum up & Unit conversion: Carbon pool [kgC/m2] -> [PgC], Carbon flux [kgC/m2/month] -> [PgC/yr]!
!   Do lat = LatNoStart, LatNoEnd
!      !x: Unit conversion factor
!      !   Area of this grid cell (m2) * (10^-12)
!      !   entire circumference of the earth is 40000km
!      x = (4.0/real(LonMax)) * (2.0/real(LatMax)) *100.0
!      
!      !y: Adjuestment values for latitude
!      y = cos( 3.141592 * (real(lat)/real(LatMax)-0.5) )
!      
!      Do lon = LonNoStart, LonNoEnd
!         sum_wbiomass  = sum_wbiomass  + x * y * data2(lat,lon,12,1) ! 1, Carbon in Woody biomass
!         sum_gbiomass  = sum_gbiomass  + x * y * data2(lat,lon,12,2) ! 2, Carbon in Grass biomass
!         sum_litter    = sum_litter    + x * y * data2(lat,lon,12,3) ! 3, Carbon in litter
!         sum_som_int   = sum_som_int   + x * y * data2(lat,lon,12,4) ! 4, Carbon in som_int
!         sum_som_slow  = sum_som_slow  + x * y * data2(lat,lon,12,5) ! 5, Carbon in som_slow
!         sum_gpp       = sum_gpp       + x * y * sum(data2(lat,lon,1:12, 9)) ! 9, GPP
!         sum_npp       = sum_npp       + x * y * sum(data2(lat,lon,1:12,10)) !10, NPP
!         sum_nep       = sum_nep       + x * y * sum(data2(lat,lon,1:12,11)) !11, NEP
!         sum_hr        = sum_hr        + x * y * sum(data2(lat,lon,1:12,25)) !12, HR
!         sum_runoff    = sum_runoff    +     y * sum(data2(lat,lon,1:12,15)) !15, runoff
!         sum_intercept = sum_intercept +     y * sum(data2(lat,lon,1:12,16)) !16, interception
!         sum_evapor    = sum_evapor    +     y * sum(data2(lat,lon,1:12,17)) !17, evaporation
!         sum_transpi   = sum_transpi   +     y * sum(data2(lat,lon,1:12,18)) !18, transpiration
!         
!         sum_weight = sum_weight + y
!      End Do
!   End Do
!   
!   !Water flux (grids average, not grid sumup)
!   sum_runoff    = sum_runoff    / sum_weight
!   sum_intercept = sum_intercept / sum_weight
!   sum_evapor    = sum_evapor    / sum_weight
!   sum_transpi   = sum_transpi   / sum_weight
!   
!   !Write variables into a file
!   Open (1, file=Loc_analysis_files//'out_analysis0.txt')
!   write(1,*) 'Carbon in Woody biomass [PgC]   ', sum_wbiomass
!   write(1,*) 'Carbon in Grass biomass [PgC]   ', sum_gbiomass
!   write(1,*) 'Carbon in Soil          [PgC]   ', sum_litter + sum_som_int + sum_som_slow  
!   write(1,*) 'GPP                     [PgC/yr]', sum_gpp
!   write(1,*) 'NPP                     [PgC/yr]', sum_npp
!   write(1,*) 'NEP                     [PgC/yr]', sum_nep
!   write(1,*) 'Heterotrophic resp.     [PgC/yr]', sum_hr 
!   write(1,*) 'Runoff                  [mm/yr] ', sum_runoff    
!   write(1,*) 'Intercepted             [mm/yr] ', sum_intercept 
!   write(1,*) 'Evaporation             [mm/yr] ', sum_evapor    
!   write(1,*) 'Transpiration           [mm/yr] ', sum_transpi   
!   Close (1)
   
!********************************************************************
! Format of result file of each simulation grid
! * Written in text format
! * Monthly data: line number = YearForMean * 12
! * Each line contains following values (tab delimited)
! 
! @data1(Lat, Lon, Mon, DataNo): Integer variables
! 1: Biome no (classfication)
! 2: Number of no water stress days [day/year]
!
! @data2(Lat, Lon, Mon, DataNo): Real variables
! 1: Carbon in Woody biomass [kg C / m2]
! 2: Carbon in Grass biomass [kg C / m2]
! 3: Carbon in litter        [kg C / m2]
! 4: Carbon in som_int       [kg C / m2]
! 5: Carbon in som_slow      [kg C / m2]
! 6: Water in top soil layer [mm]
! 7: Water in soil layer 2   [mm]
! 8: Water in snow           [mm]
! 9: GPP                     [kg C / m2 / month]
!10: NPP                     [kg C / m2 / month]
!11: NEP                     [kg C / m2 / month]
!12: mean_woody_mass         [kg C / tree]
!13: LAI of woody PFTs       [m2/m2]
!14: LAI of grass PFTs       [m2/m2]
!15: runoff                  [mm/month]
!16: interception            [mm/month]
!17: evaporation             [mm/month]
!18: transpiration           [mm/month]
!19: tree_density            [N / ha]
!20: monthly stomatal conductance [mol H2O m-2 s-1]
!21: Fire frequency
!22: Woody    NPP [kg C / m2 / month]
!23: C3 Grass NPP [kg C / m2 / month]
!24: C4 Grass NPP [kg C / m2 / month]
!25: HR           [kg C / m2 / month]
!26: Fraction of crown coverage  [fraction]
!27: Active Layer Depth [m]
!28: Available water on the top 5 soil layers [mm]
!********************************************************************

END Subroutine after_sim1


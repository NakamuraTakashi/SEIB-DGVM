!*************************************************************************************************
! Start up procedure for multiple grid cells
! using parallel computation with MPI
!*************************************************************************************************
PROGRAM omp_kickoff
   use mpi
   implicit none
   
!_____________ Set Parameters
!Directory for writing an output file of each simulation grid
   character*128,parameter::Loc_result_files  = '/home/hsato/work/result/'
   
!Directory for writing analysis data for whole simulation grid
   character*128,parameter::Loc_analysis_files= '/home/hsato/work/result_visualized/'
   
!Specify Land mask file
   character*128,parameter::Fn_landmask = '/data/landmask_2.0deg.txt'
  !character*128,parameter::Fn_landmask = '/data/landmask_1.0deg.txt'
  !character*128,parameter::Fn_landmask = '/data/landmask_0.5deg.txt'
   
!Specify soil property data
   !A same data-set is employed at every simulation grids
   !This data-set was obtained from GSWP2 (http://www.iges.org/gswp2/)
   character*128,parameter::Fn_landprop       = 'land_prop.txt'
   
!●
!Directory of climate data
   character*128,parameter::Loc_climate_data  = '/data/climate_2.0deg/'
  !character*128,parameter::Loc_climate_data  = '/data/climate_0.5deg/'
  !character*128,parameter::Loc_climate_data  = '/data/climate_future_0.5deg/'
   
!●
!Set year number of climate data
   integer,parameter::YearMaxClimate    =  105 !20c data (1901~2005, CRU-NCEP/NCAR mix)
  !integer,parameter::YearMaxClimate    =  100 !21c data (2001~2100, CRU-MIROC-NCEP/NCAR mix)
   
!Maximum grid number for longitude and latitude
   integer,parameter::LatMax        =   90 !Grid number for latitude  (@2.0deg)
   integer,parameter::LonMax        =  180 !Grid number for longitude (@2.0deg)
  !integer,parameter::LatMax        =  180 !Grid number for latitude  (@1.0deg)
  !integer,parameter::LonMax        =  360 !Grid number for longitude (@1.0deg)
  !integer,parameter::LatMax        =  360 !Grid number for latitude  (@0.5deg)
  !integer,parameter::LonMax        =  720 !Grid number for longitude (@0.5deg)
   
!Grid length for simulation
  !(Whole globe @ every Grid mesh)
  integer,parameter::LatNoStart    =      1 !
  integer,parameter::LatNoEnd      = LatMax !
  integer,parameter::LonNoStart    =      1 !
  integer,parameter::LonNoEnd      = LonMax !
  !
  !(African continent @ 2.0deg grid mesh)
  !integer,parameter::LatNoStart    =     27 !
  !integer,parameter::LatNoEnd      =     63 !
  !integer,parameter::LonNoStart    =     80 !
  !integer,parameter::LonNoEnd      =    116 !
  
  !(African continent @ 1.0deg grid mesh)
  !integer,parameter::LatNoStart    =     54 !
  !integer,parameter::LatNoEnd      =    126 !
  !integer,parameter::LonNoStart    =    160 !
  !integer,parameter::LonNoEnd      =    232 !
  !
  !(African continent @ 0.5deg grid mesh)
  !integer,parameter::LatNoStart    =    108 !
  !integer,parameter::LatNoEnd      =    252 !
  !integer,parameter::LonNoStart    =    320 !
  !integer,parameter::LonNoEnd      =    464 !
  
!Period for writing monthly outputs @ year
  !Monthly data is written for last YearForMean years of the simulation
  integer,parameter::YearForMean = 10
  
!_____________ Set Variables
   !MPI control variables
   integer myid, numprocs, ierr
   
   !Coodinate variables for parallel computation
   integer gridNo, gridNoMax        !Sequential number for simulation grid cell
   integer latNo, lonNo             !Sequential number for north-south and west-east
   real    lat,lon                  !Latitude and Longitude
   integer point                    !Sequential number for 1.0deg grids to refer soil properties
   integer,dimension(LatMax*LonMax)::latNo_save, lonNo_save
   
   !Land Ocean mask
   logical,dimension(LonMax,LatMax)::Landmask     !landmask; true->land, false->water
   integer,dimension(1:LonMax)     ::readerLonMax !For reading land mask
   
   !Location data
   integer,dimension(360*180)::Mask         !Land ocean mask at 1.0 deg (1:land, 0:ocean)
   real   ,dimension(360*180)::ALT          !altitude (m above MSL)
   real   ,dimension(360*180)::Albedo_soil0 !albedo, default
   real   ,dimension(360*180)::W_fi         !filed capacity   (m3/m3, 0.0 -> 1.0)
   real   ,dimension(360*180)::W_wilt       !wilting point    (m3/m3, 0.0 -> 1.0)
   real   ,dimension(360*180)::W_sat        !saturate point   (m3/m3, 0.0 -> 1.0)
   real   ,dimension(360*180)::W_mat        !matrix potential (m, -0.0001 -> -3.0)
   
   !Atomospheric CO2 time-series @ ppm
   !(1850~2000 from historical record, +2001~2100 from a1b scenario)
   real,dimension(351)::aco2_1750to2100
   
   !Counter
   integer i
   
!______________ Prepare landmask (extension for wide area simulation)
!Read soil property data 
   open (1, file=Fn_landprop, status='OLD')
   do i=1, 360*180
      read(1,*) Mask(i), ALT(i), Albedo_soil0(i), W_sat(i), W_fi(i), W_mat(i), W_wilt(i)
      if (W_fi  (i) > W_sat(i) ) W_fi  (i) = W_sat(i)
      if (W_wilt(i) > W_sat(i) ) W_wilt(i) = W_sat(i)
   end do
   close(1)
   
!Read Landmask, which is made from climate data
   landmask(:,:) = .false.
   open (1, file=Fn_landmask, status='OLD')
   do latNo=1, LatMax
      read(1,*) readerLonMax(1:LonMax)
      do lonNo=1, LonMax
         if (readerLonMax(lonNo)==1) landmask(lonNo,latNo)=.true.
      end do
   end do
   close(1)
   
!Adjust landmask using soil property data
   do latNo = LatNoStart, LatNoEnd
   do lonNo = LonNoStart, LonNoEnd
      !Refer landmask, which is made from limatic data
      if (.not. landmask(lonNo,latNo)) cycle
      
      !Refer soil property, which is made from GSWP2 data
      lat    =    90.0 - (real(latNo)-0.5) * (180.0/real(LatMax))
      lon    = - 180.0 + (real(lonNo)-0.5) * (360.0/real(LonMax)) 
      point  = (90-int(lat)-1)*360 + int(lon+180) + 1 !grid point number @1.0 deg system
      
      !Turn off landmasks for grid cells whose land properties are inappropriate
      if (Mask        (point) ==    0) landmask(lonNo,latNo) = .false.
      if (Albedo_soil0(point) <=  0.0) landmask(lonNo,latNo) = .false.
      if (W_sat       (point) <=  0.0) landmask(lonNo,latNo) = .false.
      if (W_fi        (point) <=  0.0) landmask(lonNo,latNo) = .false.
      if (W_wilt      (point) <=  0.0) landmask(lonNo,latNo) = .false.
   end do
   end do
   
!______________ Read atmospheric CO2 time-series during years for 1750~2100
   Open (1, file='co2_1750-2100_a1b.txt', status='OLD')
   do i = 1, 351
      read(1, *) aco2_1750to2100(i)
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
   
!______________ Initialization procedure for MPI computation
   Call MPI_INIT( ierr )
   Call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
   Call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
   
!______________ Conduct parallel computation with MPI
   if (myid==1) write(*,*) 'Start parallel computation....'
   
   DO gridNo = 1, gridNoMax
   if ( myid == mod(gridNo-1, numprocs) ) then
      Call start(gridNo, latNo_save(gridNo), lonNo_save(gridNo), aco2_1750to2100, &
                 LatMax, LonMax, YearMaxClimate, Loc_climate_data, Loc_result_files, YearForMean)
   endif
   END DO
   
!______________ Termination procudure for MPI compuation
   Write(*,*) 'Terminate processor:', myid, ', Total processor:', numprocs
   Call MPI_Barrier(MPI_COMM_WORLD, ierr) !Wait until end of all simulation
   Call MPI_FINALIZE(ierr)                !Finalize procudure for MPI
   
!______________ Conversion of output files of simulation grids
   if (myid==0) then
      write(*,*) 'Converting result files....'
      !Analysis output values of Last 10 yers
      Call after_sim1(LatMax, LonMax, LatNoStart, LatNoEnd, LonNoStart, LonNoEnd, &
                      Loc_result_files, Loc_analysis_files, landmask)
      !Analysis output values of time series
      Call after_sim2(LatMax, LonMax, LatNoStart, LatNoEnd, LonNoStart, LonNoEnd, &
                      Loc_result_files, Loc_analysis_files, landmask)
   endif
   
STOP
END PROGRAM omp_kickoff



!*************************************************************************************************
! Start up procedure for each simulation grid
!*************************************************************************************************
Subroutine start (gridNo, latNo, lonNo, aco2_1750to2100, &
                  LatMax, LonMax, YearMaxClimate, Loc_climate_data, Loc_result_files, YearForMean)
   
   USE data_structure
   implicit none
   
!_____________ Set Augment
   integer,intent(IN):: gridNo, latNo, lonNo, LatMax, LonMax, YearMaxClimate, YearForMean
   character*128,intent(IN):: Loc_climate_data, Loc_result_files
   real,dimension(351),intent(IN)::aco2_1750to2100 !Atomospheric co2 concentration (ppm)
   
!_____________ Set Variables
!Climate data
   real,dimension(Day_in_Year,YearMaxClimate)::tmp_air   !2m air temperature(Celcius)
   real,dimension(Day_in_Year,YearMaxClimate)::tmp_soil1 !soil temperature   0- 10cm depth(Celcius)
   real,dimension(Day_in_Year,YearMaxClimate)::tmp_soil2 !soil temperature  10-200cm depth(Celcius)
   real,dimension(Day_in_Year,YearMaxClimate)::tmp_soil3 !soil temperature 200-300cm depth(Celcius)
   real,dimension(Day_in_Year,YearMaxClimate)::cloud     !total cloudiness (fraction)
   real,dimension(Day_in_Year,YearMaxClimate)::prec      !precipitation (mm day-1)
   real,dimension(Day_in_Year,YearMaxClimate)::humid     !specific humidity (kg kg-1)
   real,dimension(Day_in_Year,YearMaxClimate)::wind      !wind velocity (m s-1)
   real,dimension(Day_in_Year,YearMaxClimate)::tmp_air_dr !daily range of air temperature (Celcius)
   
!Location data
   integer::Mask         !Land ocean mask at 1.0 deg (1:land, 0:ocean)
   real   ::ALT          !altitude (m above MSL)
   real   ::Albedo_soil0 !albedo, default
   real   ::W_fi         !filed capacity   (m3/m3, 0.0 -> 1.0)
   real   ::W_wilt       !wilting point    (m3/m3, 0.0 -> 1.0)
   real   ::W_sat        !saturate point   (m3/m3, 0.0 -> 1.0)
   real   ::W_mat        !matrix potential (m, -0.0001 -> -3.0)
   
!Others
   real    LAT, LON                   !latitude and logitude for simulate
   integer year, doy, data_num, point !Counters
   integer i                          !for general usage
   
!_____________ Set Variables, for wide area computations
!For reading data
   real,dimension(Day_in_Year, YearMaxClimate, 1:9)::dataREAD  !1グリッドの気象データ読出用
   
!For I/O
   character(len= 3) nam_lat
   character(len= 3) nam_lon
   integer i1, i2, i3
   integer file_no_grid1, file_no_grid2, file_no_grid3, file_no_grid4
   
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
   
!______________________Set varieties of reference number for this grid cell
!Device number for I/O
   file_no_grid1 = gridNo + 100 !For output files1, and for general usage
   file_no_grid2 = gridNo + 300 !For reading spinup files
   file_no_grid3 = gridNo + 500 !For writing spinup files
   file_no_grid4 = gridNo + 700 !For output files2
   
!Reference number for location
   ! LAT: north +, south - (decimalized)
   ! LON: east  +, west  - (decimalized)
   LAT   =    90.0 - (real(latNo)-0.5) * (180.0/real(LatMax))
   LON   = - 180.0 + (real(lonNo)-0.5) * (360.0/real(LonMax)) 
   point = (90-int(LAT)-1)*360 + int(LON+180) + 1 !grid point number @1.0 deg system
   
!Characters for designating latitude and longitude
   i1 = int(  latNo               /100 ) ; nam_lat(1:1) = char(i1+48)
   i2 = int( (latNo -i1*100)      /10  ) ; nam_lat(2:2) = char(i2+48)
   i3 =    (  latNo -i1*100-i2*10      ) ; nam_lat(3:3) = char(i3+48)
   i1 = int(  lonNo               /100 ) ; nam_lon(1:1) = char(i1+48)
   i2 = int( (lonNo -i1*100)      /10  ) ; nam_lon(2:2) = char(i2+48)
   i3 =    (  lonNo -i1*100-i2*10      ) ; nam_lon(3:3) = char(i3+48)
   
!______________ Read Location data
   open (file_no_grid1, file=Fn_location, status='OLD')
   do i=1, point
      read(file_no_grid1,*) Mask, ALT, Albedo_soil0, W_sat, W_fi, W_mat, W_wilt
   end do
   close(file_no_grid1)
   if (W_fi   > W_sat ) W_fi   = W_sat
   if (W_wilt > W_sat ) W_wilt = W_sat
   
!_____________ Prepare Climate Data
   !Read binary data
   i = 9 * 365 * YearMaxClimate
   Open (file_no_grid1, file=trim(Loc_climate_data)//nam_lat//'/'//nam_lon//'.dat',&
   access='direct',recl=i)
      read (file_no_grid1,rec=1) &
      (((dataREAD(doy,year,data_num), data_num=1,9), doy=1,365), year=1,YearMaxClimate)
   Close(file_no_grid1)
   
   !Name change of climatic variables
   do year=1,YearMaxClimate
   do doy=1,Day_in_Year
      wind      (doy,year) = dataREAD(doy,year,1)
      tmp_air   (doy,year) = dataREAD(doy,year,2)
      tmp_soil1 (doy,year) = dataREAD(doy,year,3)
      tmp_soil2 (doy,year) = dataREAD(doy,year,4)
      tmp_soil3 (doy,year) = dataREAD(doy,year,5)
      prec      (doy,year) = dataREAD(doy,year,6)
      humid     (doy,year) = dataREAD(doy,year,7)
      cloud     (doy,year) = dataREAD(doy,year,8)
      tmp_air_dr(doy,year) = dataREAD(doy,year,9)
   enddo
   enddo
   
!_____________ !Call simulation loop for each grid
   !Open I/O files
   !for writing output files
   if (Flag_output_write) then
      open (file_no_grid1, &
      file = trim(Loc_result_files)//nam_lat//'_'//nam_lon//'_out.txt')
      
      open (file_no_grid4, &
      file = trim(Loc_result_files)//nam_lat//'_'//nam_lon//'_out2.txt')
   endif
  
   !for reading spinup files
   if (Flag_spinup_read) then
     open ( file_no_grid2, &
!●
!    file = trim(Loc_result_files)//nam_lat//'_'//nam_lon//'_spnin.txt', &
     file = trim(Loc_result_files)//'spinin_2000/'//nam_lat//'_'//nam_lon//'_spnin.txt', &
     status='OLD')
   endif
   
   !for writing spinup files
   if (Flag_spinup_write) then
      open ( file_no_grid3, &
      file = trim(Loc_result_files)//nam_lat//'_'//nam_lon//'_spnout.txt')
   endif
   
   Call main_loop ( &
   LAT, LON, YearMaxClimate, &
   tmp_air(:,:), tmp_soil1(:,:), tmp_soil2(:,:), tmp_soil3(:,:), &
   cloud(:,:), prec(:,:), humid(:,:), wind(:,:), aco2_1750to2100, &
   ALT, Albedo_soil0, W_fi, W_wilt, W_sat, W_mat, &
   file_no_grid1, file_no_grid2, file_no_grid3, file_no_grid4, YearForMean)
   
   !Close I/O files
   close ( file_no_grid1 )
   close ( file_no_grid2 )
   close ( file_no_grid3 )
   close ( file_no_grid4 )
   
END Subroutine start



!*************************************************************************************************
! Output file maker1 for multi-grid computations
! (Called at the end of each simulation month during last YearForMean years of simulation)
!*************************************************************************************************
SUBROUTINE output_global (Fn, W_fi)

!_____________ Set variables
!Namespace
   USE data_structure
   USE time_counter
   USE vegi_status_current1
   USE vegi_status_current2
   USE grid_status_current1
   USE grid_status_current2
   implicit none
   
!Arguments
   integer,intent(IN)::Fn   !File I/O number
   real   ,intent(IN)::W_fi !filed capacity   (m3/m3, 0.0 -> 1.0)
   
!Local variables
   !for output variables
   real    woody_carbon, grass_carbon, som1, som2, som3   !(kg C / m2)
   real    gpp_out, npp_out, nep_out, npp_tree, npp_C3g, npp_C4g !(kg C / m2 / month)
   real    mean_woody_mass      !(g C / tree)
   real    tree_density         !(N / ha)
   real    lai_wood, lai5, lai6 !(m2/m2)
   integer dry_days             !(day/year)
   
   !for general usage
   integer  no, i, j, k, p
   real     x, y
   
!_____________ Prepare output variables
   !common procedure
   x = 0.0 !sum of tree biomass (g dm)
   i = 0   !sum of tree number        
   do no=1, Max_no
   if (tree_exist(no)) then
      x = x + mass_leaf(no)+mass_trunk(no)+mass_root(no)+mass_stock(no)+mass_available(no)
      i = i + 1
   endif
   enddo
   
   !tree_density (N / ha)
   tree_density = 10000.0 * real(i) / Max_loc / Max_loc
   
   !mean_woody_mass (kg C / tree)
   mean_woody_mass = x * C_in_drymass / max(1.0, real(i)) / 1000.0
   
   !woody_carbon (kg C / m2)
   woody_carbon = x * C_in_drymass / Max_loc / Max_loc / 1000.0
   
   !grass_carbon (kg C / m2)
   grass_carbon = sum(gmass_leaf(:,:)) + sum(gmass_root(:,:)) + sum(gmass_stock(:,:)) + sum(gmass_available(:,:))
   grass_carbon = grass_carbon * C_in_drymass / 1000.0 / Max_loc / Max_loc
   
   !soil_carbon (kg C / m2)
   som1 = (pool_litter_trunk + pool_litter_leaf + pool_litter_root + pool_litter_ag + pool_litter_bg) &
                        * C_in_drymass / Max_loc / Max_loc / 1000.0
   som2 = pool_som_int  * C_in_drymass / Max_loc / Max_loc / 1000.0
   som3 = pool_som_slow * C_in_drymass / Max_loc / Max_loc / 1000.0
   
   !GPP & NPP (kg C / m2 / month)
   i = Day_in_month(Month(doy)) !number of day of the current month (day)
   gpp_out   = Sum(gpp_RunningRecord(1:i, :)) * C_in_drymass / Max_loc / Max_loc / 1000.0
   npp_out   = Sum(npp_RunningRecord(1:i, :)) * C_in_drymass / Max_loc / Max_loc / 1000.0
   
   npp_tree  = 0.0
   do p=1, PFT_no
      if (p==C3g_no .or. p==C4g_no) cycle
      npp_tree = npp_tree + Sum( npp_RunningRecord(1:i, p) )
   enddo
   npp_tree = npp_tree                             * C_in_drymass / Max_loc / Max_loc / 1000.0
   
   npp_C3g   = Sum(npp_RunningRecord(1:i, C3g_no)) * C_in_drymass / Max_loc / Max_loc / 1000.0
   npp_C4g   = Sum(npp_RunningRecord(1:i, C4g_no)) * C_in_drymass / Max_loc / Max_loc / 1000.0
   
   !NEP (kg C / m2 / month)
   i  = Day_in_month(Month(doy))
   nep_out = sum(flux_c_uptake_RR(1:i)) &
           - sum(flux_c_mnt_RR(1:i)) - sum(flux_c_gro_RR(1:i)) &
           - sum(flux_c_htr_RR(1:i)) - sum(flux_c_fir_RR(1:i))  
   nep_out = nep_out * C_in_drymass / Max_loc / Max_loc / 1000.0
   
   
   !Tree LAI (m2/m2)
   lai_wood = sum(la) / Max_loc / Max_loc
   
   !LAI of PFT5 and PFT6
   lai5 = 0.0
   lai6 = 0.0
   do no=1, Max_no
   if (tree_exist(no)) then
      p = pft(no)
      if (p==5) lai5=lai5+la(no)
      if (p==6) lai6=lai6+la(no)
   endif
   enddo
   lai5 = lai5 / Max_loc / Max_loc
   lai6 = lai6 / Max_loc / Max_loc
   
   !dry_days: Number of dry days based on Priestley-Taylor model
   dry_days=0
   k=0
   do i=1, 12                 !i: month
      x = 0.0
      y = 0.0
      do j=1, Day_in_month(i)            !j: day form the begginig of this month
         k = k + 1                       !k: doy
         x = x + Prec_RunningRecord  (k) !x: precipitation                 (mm/month)
         y = y + ev_pot_RunningRecord(k) !x: potential ecapotranspirationn (mm/month)
      enddo
      if (x<y) dry_days = dry_days + Day_in_month(i)
   enddo
   
   !dry_days: Number of dry days based on Priestley-Taylor model
   dry_days=0
   do i=1, Day_in_Year
      if (i>30) then
         x = sum(Prec_RunningRecord   (i-30:i)) !monthly precipitation         (mm/month)
         y = sum(ev_pot_RunningRecord (i-30:i)) !potential ecapotranspirationn (mm/month)
      else
         x = sum(Prec_RunningRecord   (1:i)) + sum(Prec_RunningRecord   (Day_in_Year-30+i:Day_in_Year))
         y = sum(ev_pot_RunningRecord (1:i)) + sum(ev_pot_RunningRecord (Day_in_Year-30+i:Day_in_Year))
      endif
      
      if (x<y) dry_days = dry_days + 1
   enddo
   
!_____________ Write output-data
write (Fn,'( 2(1x,i3), 7(1x,f7.2), 1(1x,f9.1), 3(1x,f8.3), 1(1x,f7.1), 2(1x,f4.1), 4(1x,f6.1), 1x, f6.0, 1x, f6.2, 1x,f5.3, 3(1x,f8.3), 2(1x,f4.1), 1x, f5.3 )') &
   biome               , & ! 1: Biome no (classfication)
   Day_in_Year-dry_days, & ! 2: Number of no water stress days [day/year]
   
   woody_carbon     , & ! 1: Carbon in Woody biomass [kg C / m2]
   grass_carbon     , & ! 2: Carbon in Grass biomass [kg C / m2]
   som1             , & ! 3: Carbon in litter        [kg C / m2]
   som2             , & ! 4: Carbon in som_int       [kg C / m2]
   som3             , & ! 5: Carbon in som_slow      [kg C / m2]
   
   sum(pool_w(1: 5)) / ( 5*W_fi*Depth), & ! 6: Water in soil layer 1 [fraction]
   sum(pool_w(6:15)) / (10*W_fi*Depth), & ! 7: Water in soil layer 2 [fraction]
   
   pool_snow        , & ! 8: Water in snow         [mm]
   gpp_out          , & ! 9: GPP [kg C / m2 / month]
   npp_out          , & !10: NPP [kg C / m2 / month]
   nep_out          , & !11: NEP [kg C / m2 / month]
   
   mean_woody_mass  , & !12: mean_woody_mass   [kg C / tree]
   
   lai_wood                             , & !13: LAI of woody PFTs [m2/m2]
   sum(lai_grass(:,:)) /DivedG /DivedG  , & !14: LAI of grass PFTs [m2/m2]
   
   sum(flux_ro_RunningRecord(1:Day_of_Month(doy))), & !15: runoff        [mm/month]
   sum(flux_ic_RunningRecord(1:Day_of_Month(doy))), & !16: interception  [mm/month]
   sum(flux_ev_RunningRecord(1:Day_of_Month(doy))), & !17: evaporation   [mm/month]
   sum(flux_tr_RunningRecord(1:Day_of_Month(doy))), & !18: transpiration [mm/month]
   
   tree_density                                   , & !19: tree_density         [N / ha]
   canopy_cond                                    , & !20: stomatal conductance [mol H2O m-2 s-1]
   real(fire_number)/(real(counter)/real(Day_in_Year)), & !21: Fire frequency   [n/year]
   
   npp_tree                                       , & !22: Woody    NPP   [kg C / m2 / month]
   npp_C3g                                        , & !23: C3 grass NPP   [kg C / m2 / month]
   npp_C4g                                        , & !24: C4 grass NPP   [kg C / m2 / month]
   lai5                                           , & !25: LAI of PFT5    [m2/m2]
   lai6                                           , & !26: LAI of PFT6    [m2/m2]
   frac_crown_coverage                                !27: Crown coverage [fraction]
   
END SUBROUTINE output_global



!*************************************************************************************************
! Output file maker1 for multi-grid computations
! (Called at the end of each simulation year)
!*************************************************************************************************
SUBROUTINE output_global2 (Fn)

!_____________ Set variables
!Namespace
   USE data_structure
   USE time_counter
   USE vegi_status_current1
   USE vegi_status_current2
   USE grid_status_current1
   USE grid_status_current2
   implicit none
   
!Arguments
   integer,intent(IN)::Fn   !File I/O number
   
!Local variables
   real,dimension(12)::out_timeseries
   
!_____________ Prepare output variables
   !Biomass
   out_timeseries(1) = sum(mass_leaf(:)) + sum(mass_trunk(:)) + sum(mass_root(:)) + &
                       sum(mass_stock(:)) + sum(mass_available(:)) + &
                       sum(gmass_leaf(:,:)) + sum(gmass_root(:,:)) + &
                       sum(gmass_available(:,:)) + sum(gmass_stock(:,:))
   !Litter woody
   out_timeseries(2) = pool_litter_trunk + pool_litter_leaf + pool_litter_root
   
   !Litter grass
   out_timeseries(3) = pool_litter_ag + pool_litter_bg
   
   !SOM fast
   out_timeseries(4) = pool_som_int
   
   !SOM slow
   out_timeseries(5) = pool_som_slow
   
   !Carbon fluxes
   out_timeseries( 6) = sum(flux_c_uptake_RR(:)) !uptake
   out_timeseries( 7) = sum(flux_c_mnt_RR(:)   ) !due to maintenance   respiration
   out_timeseries( 8) = sum(flux_c_gro_RR(:)   ) !due to growth        respiration
   out_timeseries( 9) = sum(flux_c_htr_RR(:)   ) !due to heterotrophic respiration
   out_timeseries(10) = sum(flux_c_fir_RR(:)   ) !due to fire incident 
   
   !LAI max
   out_timeseries(11) = lai_max_wood   !Wood  LAI @ annual maximum (m2/m2)
   out_timeseries(12) = lai_max_grass  !Grass LAI @ annual maximum (m2/m2)
   
   !Unit conversion (gDM/stand -> kgC/m2)
   out_timeseries(1:10) = out_timeseries(1:10) * C_in_drymass / Max_loc / Max_loc / 1000
   
!_____________ Write output-data
   write(Fn, '( 5(f8.3,a), 5(f9.4,a), f4.1, a, f4.1 )' ) &
   out_timeseries( 1), ',', &
   out_timeseries( 2), ',', &
   out_timeseries( 3), ',', &
   out_timeseries( 4), ',', &
   out_timeseries( 5), ',', &
   out_timeseries( 6), ',', &
   out_timeseries( 7), ',', &
   out_timeseries( 8), ',', &
   out_timeseries( 9), ',', &
   out_timeseries(10), ',', &
   out_timeseries(11), ',', &
   out_timeseries(12)
   
END SUBROUTINE output_global2



!*************************************************************************************************
! Analysis output values of Last 10 yers
!*************************************************************************************************
Subroutine after_sim1(LatMax, LonMax, LatNoStart, LatNoEnd, LonNoStart, LonNoEnd, &
                      Loc_result_files, Loc_analysis_files, landmask)
   
   implicit none
   
!_____________ Set Augment
   integer      ,intent(IN):: LatMax, LonMax, LatNoStart, LatNoEnd, LonNoStart, LonNoEnd
   character*128,intent(IN):: Loc_result_files, Loc_analysis_files
   
   logical,dimension(LonMax, LatMax),intent(IN):: landmask !Land Ocean mask
   
!_____________ Set parameters
!○
   !潜在植生の分布図（解像度に応じて異なるデータセットを選択すること）
   character*128,parameter::Fn_potveg = 'potential_veg_hd.asc'
   
!_____________ Set variables
!Year lenght for monthly outputs in simulation files
   integer YearForMean
   
!Array for inputing result files
   integer,dimension(                                               2)::data1_read
   integer,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd, 12,  2)::data1
   
   real   ,dimension(                                              27)::data2_read
   real   ,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd, 12, 27)::data2
   
!Variables for output
   integer,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd)::out_biome_consistant
   
   real,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd)::out_water
   real,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd)::out_precipitation
   real,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd)::out_gpp
   real,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd)::out_npp
   
   !Annual mean of LAI
   real,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd)::out_lai_amean    !All PFT
   real,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd)::out_lai_amean_t  !Woody PFTs
   real,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd)::out_lai_amean_g  !Grass PFTs
   real,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd)::out_lai_amean_5  !PFT5
   real,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd)::out_lai_amean_6  !PFT6
   
   !Annula maximum of LAI
   real,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd)::out_lai_max      !All PFT
   real,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd)::out_lai_max_t    !Woody PFTs
   real,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd)::out_lai_max_g    !Grass PFTs
   real,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd)::out_lai_max_5    !PFT5
   real,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd)::out_lai_max_6    !PFT6
   
!For summary statistics
   real sum_weight
   real sum_wbiomass  ! 1, Carbon in Woody biomass [kg C / m2]
   real sum_gbiomass  ! 2, Carbon in Grass biomass [kg C / m2]
   real sum_litter    ! 3, Carbon in litter        [kg C / m2]
   real sum_som_int   ! 4, Carbon in som_int       [kg C / m2]
   real sum_som_slow  ! 5, Carbon in som_slow      [kg C / m2]
   real sum_gpp       ! 9, GPP                     [kg C / m2 / month]
   real sum_npp       !10, NPP                     [kg C / m2 / month]
   real sum_nep       !11, NEP                     [kg C / m2 / month]
   real sum_runoff    !15, runoff                  [mm/month]
   real sum_intercept !16, interception            [mm/month]
   real sum_evapor    !17, evaporation             [mm/month]
   real sum_transpi   !18, transpiration           [mm/month]
   
!etc
   integer,dimension(100)::biome_counter
   character(len=7) fname            !For file name
   integer i1, i2, i3, i4, i5, i6    !For processing file name
   integer lat, lon, month, year, i  !Loop counter
   real    x, y                      !For general usage
   
!!_____________ Optional解析のための変数群
!   
!   real FR_ratio_c3, FR_ratio_c4
!   real RG_f_c3, RG_f_c4
!   real a1, a2, a3, a4, a5, a6
!   
!   logical,dimension(LonMax,LatMax)::Veg_pot !Potential vegetation
   
!_______________ Survey lengh of YearForMean
   !Loop for each simulation grid
   DO lat = LatNoStart, LatNoEnd
   DO lon = LonNoStart, LonNoEnd
   IF (landmask(lon,lat)==.false.) cycle
      
      !Make letter string for designating a file
      i1 = int(  lat              /100 ) ; fname(1:1)=char(i1+48)
      i2 = int( (lat-i1*100)      / 10 ) ; fname(2:2)=char(i2+48)
      i3 =    (  lat-i1*100-i2*10      ) ; fname(3:3)=char(i3+48)
                                           fname(4:4)='_'
      i4 = int(  lon              /100 ) ; fname(5:5)=char(i4+48)
      i5 = int( (lon-i4*100)      / 10 ) ; fname(6:6)=char(i5+48)
      i6 =    (  lon-i4*100-i5*10      ) ; fname(7:7)=char(i6+48)
      
      !Read output files and make average values
      open (1, file=trim(Loc_result_files)//fname//'_out.txt' )
         do YearForMean =1, 10000000
         do month =1, 12
            read(1, *, end=100) 
         enddo
         enddo
      close (1)
      
   END DO
   END DO
   
100 close (1)
   YearForMean = YearForMean - 1
   
!_______________ Read result files and make average value
write(*,*) "Reading result files and making output variables"
   !Initialize
   data1(:,:,:,:) =   0
   data2(:,:,:,:) = 0.0
   out_biome_consistant(:,:) = 0
   
   !Loop for each simulation grid
   DO lat = LatNoStart, LatNoEnd
   DO lon = LonNoStart, LonNoEnd
   IF (landmask(lon,lat)==.false.) cycle
      
      !Make letter string for designating a file
      i1 = int(  lat              /100 ) ; fname(1:1)=char(i1+48)
      i2 = int( (lat-i1*100)      / 10 ) ; fname(2:2)=char(i2+48)
      i3 =    (  lat-i1*100-i2*10      ) ; fname(3:3)=char(i3+48)
                                           fname(4:4)='_'
      i4 = int(  lon              /100 ) ; fname(5:5)=char(i4+48)
      i5 = int( (lon-i4*100)      / 10 ) ; fname(6:6)=char(i5+48)
      i6 =    (  lon-i4*100-i5*10      ) ; fname(7:7)=char(i6+48)
      
      !Read output files and make average values
      open (1, file=trim(Loc_result_files)//fname//'_out.txt' )
         
         !Sumup all values for YearForMean year (except for Biome code)
         biome_counter(:) = 0
         do year  =1, YearForMean
            do month =1, 12
               read(1,*) data1_read(:), data2_read(:)
               data1(lat,lon,month,1) =                          data1_read(1) !Biome code
               data1(lat,lon,month,2) = data1(lat,lon,month,2) + data1_read(2) !Drought days
               data2(lat,lon,month,:) = data2(lat,lon,month,:) + data2_read(:) !Other Variables
            enddo
            
            i = data1_read(1)
            biome_counter(i) = biome_counter(i) + 1
         enddo
         
         do i=1, 100
            if ( biome_counter(i) > (YearForMean/2) ) then
               out_biome_consistant(lat,lon) = i
            endif
         end do
         
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
   
   out_gpp           (:,:) = 0.0
   out_npp           (:,:) = 0.0
   
   out_lai_amean     (:,:) = 0.0
   out_lai_amean_t   (:,:) = 0.0
   out_lai_amean_g   (:,:) = 0.0
   out_lai_amean_5   (:,:) = 0.0
   out_lai_amean_6   (:,:) = 0.0
   
   out_lai_max       (:,:) = 0.0
   out_lai_max_t     (:,:) = 0.0
   out_lai_max_g     (:,:) = 0.0
   out_lai_max_5     (:,:) = 0.0
   out_lai_max_6     (:,:) = 0.0
   
   Do lat   = LatNoStart, LatNoEnd
   Do lon   = LonNoStart, LonNoEnd
      
      Do month = 1, 12
      !Soil water contents
      out_water(lat,lon)= out_water(lat,lon)+ data2(lat,lon,month, 6) / 12.0
      
      !Annual GPP
      out_gpp  (lat,lon)= out_gpp  (lat,lon)+ data2(lat,lon,month, 9)
      
      !Annual NPP
      out_npp  (lat,lon)= out_npp  (lat,lon)+ data2(lat,lon,month,10)
      
      !Annual means of LAI
      out_lai_amean  (lat,lon) = out_lai_amean  (lat,lon) + data2(lat,lon,month,13) / 12.0 &
                                                          + data2(lat,lon,month,14) / 12.0
      out_lai_amean_t(lat,lon) = out_lai_amean_t(lat,lon) + data2(lat,lon,month,13) / 12.0
      out_lai_amean_g(lat,lon) = out_lai_amean_g(lat,lon) + data2(lat,lon,month,14) / 12.0
      out_lai_amean_5(lat,lon) = out_lai_amean_5(lat,lon) + data2(lat,lon,month,25) / 12.0
      out_lai_amean_6(lat,lon) = out_lai_amean_6(lat,lon) + data2(lat,lon,month,26) / 12.0
      
      !Annual maximums of LAI
      out_lai_max  (lat,lon) = max( out_lai_max  (lat,lon), data2(lat,lon,month,13) &
                                                          + data2(lat,lon,month,14) )
      out_lai_max_t(lat,lon) = max( out_lai_max_t(lat,lon), data2(lat,lon,month,13) )
      out_lai_max_g(lat,lon) = max( out_lai_max_g(lat,lon), data2(lat,lon,month,14) )
      out_lai_max_5(lat,lon) = max( out_lai_max_5(lat,lon), data2(lat,lon,month,25) )
      out_lai_max_6(lat,lon) = max( out_lai_max_6(lat,lon), data2(lat,lon,month,26) )
      
      !Annual Precipitation
      out_precipitation(lat,lon) = out_precipitation(lat,lon)+ sum(data2(lat,lon,month,15:18))
      End Do
      
   End Do
   End Do
   
!_______________ Write output files 1
write(*,*) "Writing result maps 1"
   Open ( 1, file=trim(Loc_analysis_files)//'out_biome.txt'        )
   Open ( 2, file=trim(Loc_analysis_files)//'out_wbiomass.txt'     )
   Open ( 3, file=trim(Loc_analysis_files)//'out_water1.txt'       )
   Open ( 4, file=trim(Loc_analysis_files)//'out_gpp.txt'          )
   Open ( 5, file=trim(Loc_analysis_files)//'out_npp.txt'          )
   Open ( 6, file=trim(Loc_analysis_files)//'out_precipitation.txt')
   Open ( 7, file=trim(Loc_analysis_files)//'out_fire.txt'         )
   
   Open (10, file=trim(Loc_analysis_files)//'out_lai_max.txt'      )
   Open (11, file=trim(Loc_analysis_files)//'out_lai_max_t.txt'    )
   Open (12, file=trim(Loc_analysis_files)//'out_lai_max_g.txt'    )
   Open (13, file=trim(Loc_analysis_files)//'out_lai_max_5.txt'    )
   Open (14, file=trim(Loc_analysis_files)//'out_lai_max_6.txt'    )
   
   Open (20, file=trim(Loc_analysis_files)//'out_lai_amean.txt'    )
   Open (21, file=trim(Loc_analysis_files)//'out_lai_amean_t.txt'  )
   Open (22, file=trim(Loc_analysis_files)//'out_lai_amean_g.txt'  )
   Open (23, file=trim(Loc_analysis_files)//'out_lai_amean_5.txt'  )
   Open (24, file=trim(Loc_analysis_files)//'out_lai_amean_6.txt'  )
   
   Open (30, file=trim(Loc_analysis_files)//'out_biome_consist.txt'  )
   
   Do lat = LatNoStart, LatNoEnd
      Do lon = LonNoStart, LonNoEnd
      write( 1,'(  i2,a)', advance='no') data1     (lat,lon,12, 1),  ',' !Biome
      write( 2,'(f7.2,a)', advance='no') data2     (lat,lon,12, 1),  ',' !Woody biomass
      write( 3,'(f9.1,a)', advance='no') out_water (lat,lon)      ,  ',' !Water content @ top soil layer, annual average
      write( 4,'(f8.3,a)', advance='no') out_gpp          (lat,lon), ',' !GPP (kg C/ m2/ year)
      write( 5,'(f8.3,a)', advance='no') out_npp          (lat,lon), ',' !NPP (kg C/ m2/ year)
      write( 6,'(f9.1,a)', advance='no') out_precipitation(lat,lon), ',' !Annual precipitation
      write( 7,'(f5.3,a)', advance='no') data2            (lat,lon,12,21), ',' !Fire frequency
      
      !Annual maximum LAI
      write(10,'(f4.1,a)', advance='no') out_lai_max      (lat,lon), ',' !All PFTs
      write(11,'(f4.1,a)', advance='no') out_lai_max_t    (lat,lon), ',' !Woody PFTs
      write(12,'(f4.1,a)', advance='no') out_lai_max_g    (lat,lon), ',' !Grass PFTs
      write(13,'(f4.1,a)', advance='no') out_lai_max_5    (lat,lon), ',' !PFT5
      write(14,'(f4.1,a)', advance='no') out_lai_max_6    (lat,lon), ',' !PFT6
      
      !Annual average LAI
      write(20,'(f4.1,a)', advance='no') out_lai_amean    (lat,lon), ',' !All PFTs
      write(21,'(f4.1,a)', advance='no') out_lai_amean_t  (lat,lon), ',' !Woody PFTs
      write(22,'(f4.1,a)', advance='no') out_lai_amean_g  (lat,lon), ',' !Grass PFTs
      write(23,'(f4.1,a)', advance='no') out_lai_amean_5  (lat,lon), ',' !PFT5
      write(24,'(f4.1,a)', advance='no') out_lai_amean_6  (lat,lon), ',' !PFT6
      
      write(30,'(  i2,a)', advance='no') out_biome_consistant(lat,lon), ',' !Biome_consistent
      
      End Do
      
      !Insert feed code
      write( 1,*); write( 2,*); write( 3,*); write( 4,*); write( 5,*); write( 6,*); write( 7,*)
      write(10,*); write(11,*); write(12,*); write(13,*); write(14,*)
      write(20,*); write(21,*); write(22,*); write(23,*); write(24,*)
      write(30,*)
      
   End Do
   Close ( 1); Close ( 2); Close ( 3); Close ( 4); Close ( 5); Close ( 6); Close ( 7)
   Close (10); Close (11); Close (12); Close (13); Close (14)
   Close (20); Close (21); Close (22); Close (23); Close (24)
   Close (30)
   
!_______________ Write Summary for whole simulation grids
write(*,*) "Writing a Summry"
   !Initialization
   sum_wbiomass  = 0.0
   sum_gbiomass  = 0.0
   sum_litter    = 0.0
   sum_som_int   = 0.0
   sum_som_slow  = 0.0
   sum_gpp       = 0.0
   sum_npp       = 0.0
   sum_nep       = 0.0
   sum_runoff    = 0.0
   sum_intercept = 0.0
   sum_evapor    = 0.0
   sum_transpi   = 0.0
   
   sum_weight = 0
   
   !Sum up & Unit conversion: Carbon pool [kgC/m2] -> [PgC], Carbon flux [kgC/m2/month] -> [PgC/yr]
   Do lat = LatNoStart, LatNoEnd
      !x: Unit conversion factor
      !   Area of this grid cell (m2) * (10^-12)
      !   entire circumference of the earth is 40000km
      x = (4.0/real(LonMax)) * (2.0/real(LatMax)) *100.0
      
      !y: Adjuestment values for latitude
      y = cos( 3.141592 * (real(lat)/real(LatMax)-0.5) )
      
      Do lon = LonNoStart, LonNoEnd
         sum_wbiomass  = sum_wbiomass  + x * y * data2(lat,lon,12,1) ! 1, Carbon in Woody biomass
         sum_gbiomass  = sum_gbiomass  + x * y * data2(lat,lon,12,2) ! 2, Carbon in Grass biomass
         sum_litter    = sum_litter    + x * y * data2(lat,lon,12,3) ! 3, Carbon in litter
         sum_som_int   = sum_som_int   + x * y * data2(lat,lon,12,4) ! 4, Carbon in som_int
         sum_som_slow  = sum_som_slow  + x * y * data2(lat,lon,12,5) ! 5, Carbon in som_slow
         sum_gpp       = sum_gpp       + x * y * sum(data2(lat,lon,1:12, 9)) ! 9, GPP
         sum_npp       = sum_npp       + x * y * sum(data2(lat,lon,1:12,10)) !10, NPP
         sum_nep       = sum_nep       + x * y * sum(data2(lat,lon,1:12,11)) !11, NEP
         sum_runoff    = sum_runoff    +     y * sum(data2(lat,lon,1:12,15)) !15, runoff
         sum_intercept = sum_intercept +     y * sum(data2(lat,lon,1:12,16)) !16, interception
         sum_evapor    = sum_evapor    +     y * sum(data2(lat,lon,1:12,17)) !17, evaporation
         sum_transpi   = sum_transpi   +     y * sum(data2(lat,lon,1:12,18)) !18, transpiration
         
         sum_weight = sum_weight + y
      End Do
   End Do
   
   !Water flux (grids average, not grid sumup)
   sum_runoff    = sum_runoff    / sum_weight
   sum_intercept = sum_intercept / sum_weight
   sum_evapor    = sum_evapor    / sum_weight
   sum_transpi   = sum_transpi   / sum_weight
   
   !Write variables into a file
   Open (1, file=trim(Loc_analysis_files)//'out_summary.txt')
   write(1,*) 'Carbon in Woody biomass [PgC]   ', sum_wbiomass
   write(1,*) 'Carbon in Grass biomass [PgC]   ', sum_gbiomass
   write(1,*) 'Carbon in Soil          [PgC]   ', sum_litter + sum_som_int + sum_som_slow  
   write(1,*) 'GPP                     [PgC/yr]', sum_gpp
   write(1,*) 'NPP                     [PgC/yr]', sum_npp
   write(1,*) 'NEP                     [PgC/yr]', sum_nep       
   write(1,*) 'Runoff                  [mm/yr] ', sum_runoff    
   write(1,*) 'Intercepted             [mm/yr] ', sum_intercept 
   write(1,*) 'Evaporation             [mm/yr] ', sum_evapor    
   write(1,*) 'Transpiration           [mm/yr] ', sum_transpi   
   Close (1)
   
!_______________ Write Monthly LAI of woody PFTs
write(*,*) "Writing result maps 2 (Monthly LAI of woody PFTs)"
   Open ( 1, file=trim(Loc_analysis_files)//'out_lai_month_01.txt')
   Open ( 2, file=trim(Loc_analysis_files)//'out_lai_month_02.txt')
   Open ( 3, file=trim(Loc_analysis_files)//'out_lai_month_03.txt')
   Open ( 4, file=trim(Loc_analysis_files)//'out_lai_month_04.txt')
   Open ( 5, file=trim(Loc_analysis_files)//'out_lai_month_05.txt')
   Open ( 6, file=trim(Loc_analysis_files)//'out_lai_month_06.txt')
   Open ( 7, file=trim(Loc_analysis_files)//'out_lai_month_07.txt')
   Open ( 8, file=trim(Loc_analysis_files)//'out_lai_month_08.txt')
   Open ( 9, file=trim(Loc_analysis_files)//'out_lai_month_09.txt')
   Open (10, file=trim(Loc_analysis_files)//'out_lai_month_10.txt')
   Open (11, file=trim(Loc_analysis_files)//'out_lai_month_11.txt')
   Open (12, file=trim(Loc_analysis_files)//'out_lai_month_12.txt')
   
   Do lat = LatNoStart, LatNoEnd
      Do lon = LonNoStart, LonNoEnd
      write( 1,'(f4.1,a)', advance='no') sum( data2(lat, lon,  1, 13:14) ), ','
      write( 2,'(f4.1,a)', advance='no') sum( data2(lat, lon,  2, 13:14) ), ','
      write( 3,'(f4.1,a)', advance='no') sum( data2(lat, lon,  3, 13:14) ), ','
      write( 4,'(f4.1,a)', advance='no') sum( data2(lat, lon,  4, 13:14) ), ','
      write( 5,'(f4.1,a)', advance='no') sum( data2(lat, lon,  5, 13:14) ), ','
      write( 6,'(f4.1,a)', advance='no') sum( data2(lat, lon,  6, 13:14) ), ','
      write( 7,'(f4.1,a)', advance='no') sum( data2(lat, lon,  7, 13:14) ), ','
      write( 8,'(f4.1,a)', advance='no') sum( data2(lat, lon,  8, 13:14) ), ','
      write( 9,'(f4.1,a)', advance='no') sum( data2(lat, lon,  9, 13:14) ), ','
      write(10,'(f4.1,a)', advance='no') sum( data2(lat, lon, 10, 13:14) ), ','
      write(11,'(f4.1,a)', advance='no') sum( data2(lat, lon, 11, 13:14) ), ','
      write(12,'(f4.1,a)', advance='no') sum( data2(lat, lon, 12, 13:14) ), ','
      End Do
      
      write( 1, *); write( 2, *); write( 3, *); write( 4, *); write( 5, *); write( 6, *)
      write( 7, *); write( 8, *); write( 9, *); write(10, *); write(11, *); write(12, *)
   End Do
   Close ( 1); Close ( 2); Close ( 3); Close ( 4); Close ( 5); Close ( 6)
   Close ( 7); Close ( 8); Close ( 9); Close (10); Close (11); Close (12)
   
!!■■■■■■■■■■■　ここから下はOptionalです　■■■■■■■■■■■
!
!!_______________ Potential biome地図(from ISLSCP2)を読み出す
!write(*,*) "Reading Potential Vegetation Map"
!   
!   Veg_pot(:,:) = 0 !initialize
!   open (1, file=trim(Fn_potveg), status='OLD') !解像度に応じてデータを変えること
!   do lat=1, LatMax
!      read(1,*) Veg_pot(1:LonMax, lat)
!   end do
!   close(1)
!!1 Tropical Evergreen Forest/Woodland 
!!2 Tropical Deciduous Forest/Woodland 
!!3 Temperate Broadleaf Evergreen Forest/Woodland 
!!4 Temperate Needleleaf Evergreen Forest/Woodland 
!!5 Temperate Deciduous Forest/Woodland 
!!6 Boreal Evergreen Forest/Woodland 
!!7 Boreal Deciduous Forest/Woodland 
!!8 Mixed Forest 
!!9 Savanna 
!!10 Grassland/Steppe 
!!11 Dense Shrubland 
!!12 Open Shrubland 
!!13 Tundra  
!!14 Desert 
!!15 Polar desert/Rock/Ice 
!
!!_______________ Write output files 1.1
!write(*,*) "Writing result maps 1.1"
!   Open (1, file=trim(Loc_analysis_files)//'out_crowncover.txt')
!   Do lat = LatNoStart, LatNoEnd
!      Do lon = LonNoStart, LonNoEnd
!         write( 1,'(f7.3,a)', advance='no') data2(lat,lon,12,27),  ','
!      End Do
!      write( 1,*)
!   End Do
!   Close ( 1)
!   
!!_______________ Write output files 1.2 for sensitivity test
!write(*,*) "Writing result maps 1.2"
!   Open (1, file=trim(Loc_analysis_files)//'out_sensible_test.txt')
!   Do lat = LatNoStart, LatNoEnd
!   Do lon = LonNoStart, LonNoEnd
!      
!      !Make variables for output
!      a1 =      data2(lat,lon, 12,   1)     !Woody biomass [kg C / m2]
!      a2 =      data2(lat,lon, 12,   2)     !Grass biomass [kg C / m2]
!      a3 = a1 + a2                          !Total biomass [kg C / m2]
!      a4 = sum( data2(lat,lon, 12,   3:5) ) !Soil carbon   [kg C / m2]
!      a5 = sum( data2(lat,lon, 1:12, 10)  ) !NPP           [kg C / m2 / year]
!      a6 = sum( data2(lat,lon, 1:12, 15)  ) !Runoff        [mm/year]
!      
!      !Make variables for output
!      write(1,'(5(f7.3,a), f8.2 )') a1,',',a2,',',a3,',',a4,',',a5,',',a6
!   End Do
!   End Do
!   Close (1)
!   
!
!!_______________ Write output files 3 (Monthly LAI of all PFTs)
!write(*,*) "Writing result maps 2 (Monthly LAI of all PFTs)"
!   Open ( 1, file=trim(Loc_analysis_files)//'out_lai_01.txt')
!   Open ( 2, file=trim(Loc_analysis_files)//'out_lai_02.txt')
!   Open ( 3, file=trim(Loc_analysis_files)//'out_lai_03.txt')
!   Open ( 4, file=trim(Loc_analysis_files)//'out_lai_04.txt')
!   Open ( 5, file=trim(Loc_analysis_files)//'out_lai_05.txt')
!   Open ( 6, file=trim(Loc_analysis_files)//'out_lai_06.txt')
!   Open ( 7, file=trim(Loc_analysis_files)//'out_lai_07.txt')
!   Open ( 8, file=trim(Loc_analysis_files)//'out_lai_08.txt')
!   Open ( 9, file=trim(Loc_analysis_files)//'out_lai_09.txt')
!   Open (10, file=trim(Loc_analysis_files)//'out_lai_10.txt')
!   Open (11, file=trim(Loc_analysis_files)//'out_lai_11.txt')
!   Open (12, file=trim(Loc_analysis_files)//'out_lai_12.txt')
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
!!_______________ Write output files 4 (Values for adjusting grass NPP)
!write(*,*) "writing out_analysis1.txt"
!   Open (1, file=trim(Loc_analysis_files)//'out_analysis1.txt')
!   Do lat = LatNoStart, LatNoEnd
!   Do lon = LonNoStart, LonNoEnd
!      if (landmask(lon,lat)==.false.) cycle
!      write (1,'( 2(i4,a), f8.1,a, i4,a, 2(f8.1,a), i3,a, f7.1 )') &
!      lat                          ,',', & !
!      lon                          ,',', & !
!      out_precipitation(lat,lon)   ,',', & !precipitation (mm/year)
!      data1(lat,lon,12,2)          ,',', & !Number of no water stress days (day/year)
!      out_gpp(lat,lon)* 2.0* 1000.0,',', & !Annual GPP (g DM/ m2/ year)
!      out_npp(lat,lon)* 2.0* 1000.0,',', & !Annual NPP (g DM/ m2/ year)
!      Veg_pot(lon,lat)             ,',', & !Potential vegetation code number
!      sum(data2(lat, lon, 12, 3:5))        !Soil carbon (Kg C/ m2)
!   End Do
!   End Do
!   Close (1)
!   
!!_______________ Write output files 5 (Annual above ground production & precipitation)
!write(*,*) "Phase 7"
!   !Write Annual-precipitation and Annual above ground production
!   !for every grids whose annual precipiation is within 200~1500mm.
!   !ちなみに、Higgins et al. (2000)によると、アフリカのサバナ地帯において
!   !草本の地上部生産量(KgDM/ha/yr) = 3.37 × 降水量(mm/yr)
!   
!   !Parameters
!   FR_ratio_c3 = 0.33
!   FR_ratio_c4 = 0.33
!   RG_f_c3     = 1.5
!   RG_f_c4     = 1.5
!   
!   Open (1, file=trim(Loc_analysis_files)//'out_analysis2.txt')
!   Do lat = LatNoStart, LatNoEnd
!   Do lon = LonNoStart, LonNoEnd
!      if ( landmask(lon,lat) == .false.       ) cycle
!      if ( out_precipitation(lat,lon) > 1500.0) cycle
!      if ( out_precipitation(lat,lon) <  200.0) cycle
!      
!      a1 = sum(data2(lat,lon,:,23))*2.0*100.0*100.0 !C3 Grass NPP (Kg DM/ha/yr)
!      a2 = sum(data2(lat,lon,:,24))*2.0*100.0*100.0 !C4 Grass NPP (Kg DM/ha/yr)
!      if (a1<=0.0 .and. a2<=0.0) cycle
!      
!      !Annual aboveground yield
!      x = a1* (1.0-FR_ratio_c3)/RG_f_c3 + a2* (1.0-FR_ratio_c4)/RG_f_c4
!      
!      write (1,'( f8.1,a,f8.1 )') out_precipitation(lat,lon), ',', x
!   End Do
!   End Do
!   Close (1)
!   
!!_______________ Write output files 6 (Precipitation & Fraction of crown coverage)
!write(*,*) "Phase 8"
!   Open (1, file=trim(Loc_analysis_files)//'out_analysis3.txt')
!   Do lat = LatNoStart, LatNoEnd
!   Do lon = LonNoStart, LonNoEnd
!      if (landmask(lon,lat)==.false.) cycle
!      
!      x=0.0
!      do month=1,12
!         x = max(x, data2(lat,lon,month,27))
!      end do
!      
!      write (1,'( f8.1,a,f7.4 )') out_precipitation(lat,lon), ',', & !Precipitation [mm/year]
!                                  x !Fraction of crown coverage, Annual maximum {fraction}
!   End Do
!   End Do
!   Close (1)
!   

!********************************************************************
! (MEMO) Format of result file of each simulation grid
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
!25: LAI of PFT5  [m2/m2]
!26: LAI of PFT6  [m2/m2]
!27: Fraction of crown coverage  [fraction]
!********************************************************************

END Subroutine after_sim1





!*************************************************************************************************
! Analysis output values of time series
!*************************************************************************************************
Subroutine after_sim2(LatMax, LonMax, LatNoStart, LatNoEnd, LonNoStart, LonNoEnd, &
                      Loc_result_files, Loc_analysis_files, landmask)
   implicit none
   
!_____________ Set Augment
   integer,intent(IN):: LatMax, LonMax, LatNoStart, LatNoEnd, LonNoStart, LonNoEnd
   character*128,intent(IN):: Loc_result_files, Loc_analysis_files
   logical,dimension(LonMax, LatMax),intent(IN):: landmask !Land Ocean mask
   
!_____________ Set variables
   !input and output
   real,dimension(:,:,:,:),allocatable:: data_read
   real,dimension(1:12)               :: data_out
   
   !etc
   character(len=7) fname          !For file names
   integer Simulation_year         !length of simulation year
   integer i1, i2, i3, i4, i5, i6  !For file names
   integer lat, lon, dat, year     !Loop counter
   integer count                   !landgrid counter
   
!_______________ Survey lengh of simulation_year
   DO lat = LatNoStart, LatNoEnd
   DO lon = LonNoStart, LonNoEnd
   IF (landmask(lon,lat)==.false.) cycle
      
      !Prepare letter string of file name
      i1 = int(  lat              /100 ) ; fname(1:1)=char(i1+48)
      i2 = int( (lat-i1*100)      / 10 ) ; fname(2:2)=char(i2+48)
      i3 =    (  lat-i1*100-i2*10      ) ; fname(3:3)=char(i3+48)
                                           fname(4:4)='_'
      i4 = int(  lon              /100 ) ; fname(5:5)=char(i4+48)
      i5 = int( (lon-i4*100)      / 10 ) ; fname(6:6)=char(i5+48)
      i6 =    (  lon-i4*100-i5*10      ) ; fname(7:7)=char(i6+48)
      
      !Read output files
      open (1, file=trim(Loc_result_files)//fname//'_out2.txt' )
      do Simulation_year=1, 100000
          read(1, *, end=100)
      enddo
      
   END DO
   END DO
   
100 close (1)
   Simulation_year = Simulation_year - 1
   allocate ( data_read(1:12, LatNoStart:LatNoEnd, LonNoStart:LonNoEnd, Simulation_year) )
   
!_______________ Read output files
write(*,*) "Reading time series data"
   count = 0
   data_read(:,:,:,:) = 0.0
   DO lat = LatNoStart, LatNoEnd
   DO lon = LonNoStart, LonNoEnd
   IF (landmask(lon,lat)==.false.) cycle
      !increament of landgrid counter
      count=count+1
      
      !Prepare letter string of file name
      i1 = int(  lat              /100 ) ; fname(1:1)=char(i1+48)
      i2 = int( (lat-i1*100)      / 10 ) ; fname(2:2)=char(i2+48)
      i3 =    (  lat-i1*100-i2*10      ) ; fname(3:3)=char(i3+48)
                                           fname(4:4)='_'
      i4 = int(  lon              /100 ) ; fname(5:5)=char(i4+48)
      i5 = int( (lon-i4*100)      / 10 ) ; fname(6:6)=char(i5+48)
      i6 =    (  lon-i4*100-i5*10      ) ; fname(7:7)=char(i6+48)
      
      !Read output files
      open (1, file=trim(Loc_result_files)//fname//'_out2.txt' )
         do year=1, Simulation_year
            read(1,*) data_read(1:12, lat, lon, year)
         enddo
      close (1)
      
   END DO
   END DO
   
!_______________ Calculate and write averages of output variables
write(*,*) "Phase 2"
   Open (1, file=trim(Loc_analysis_files)//'out_timeseries.txt')
      Do year= 1, Simulation_year
         
         !initialize output variables
         data_out(:) = 0.0
         
         !calculate output variables for this year
         Do dat = 1, 12
            do lat   = LatNoStart, LatNoEnd
            do lon   = LonNoStart, LonNoEnd
                  data_out(dat) = data_out(dat) + data_read(dat,lat,lon,year)
            end do
            end do
         data_out(dat) = data_out(dat) / real(count)
         End Do
         
         !writing output variables
         ! 1: Biomass [kgC/m2]
         ! 2: Litter woody [kgC/m2]
         ! 3: Litter grass [kgC/m2]
         ! 4: SOM fast [kgC/m2]
         ! 5: SOM slow [kgC/m2]
         ! 6: Carbon fluxes (uptake) [kgC/m2/year]
         ! 7: Carbon fluxes (due to maintenance   respiration) [kgC/m2/year]
         ! 8: Carbon fluxes (due to growth        respiration) [kgC/m2/year]
         ! 9: Carbon fluxes (due to heterotrophic respiration) [kgC/m2/year]
         !10: Carbon fluxes (due to fire incident) [kgC/m2/year]
         !11: LAI max (wood)
         !12: LAI max (grass)
         write(1,'(12f7.3)') data_out(:)
         
      End Do
   Close (1)
   
   
END Subroutine after_sim2

!*************************************************************************************************
! Climate statistics recorder
!*************************************************************************************************
SUBROUTINE stat_climate (prec, tmp_air, tmp_soil)

!_____________ Set variables
!Namespace
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
   
!Augments
   real,intent(IN)::prec              !precipitation (mm day-1)
   real,intent(IN)::tmp_air           !2m air temperature (Celcius)
   real,intent(IN)::tmp_soil(NumSoil) !soil temperature for each soil layer (Celcius)
   
!Local variables
   real    tmp_air_month_ave      !monthly average of air temperature (Celcius)
   real    x
   integer count, i, j, p
   
!_____________ Main part
!air temperature
   tmp_air_RunningRecord(2:Day_in_Year) = tmp_air_RunningRecord(1:Day_in_Year-1)
   tmp_air_RunningRecord(1)             = tmp_air
   
!soil temperature
   do i=1, NumSoil
      do j=Day_in_Year, 2, -1
         tmp_soil_RunningRecord(j,i) = tmp_soil_RunningRecord(j-1,i)
      enddo
      tmp_soil_RunningRecord(1,i)              = tmp_soil(i)
   enddo
   
!precipitation
   prec_RunningRecord(2:Day_in_Year) = prec_RunningRecord(1:Day_in_Year-1)
   prec_RunningRecord(1)             = prec
   
!coldest and hottest month of the year
   !update coldest and hottest month of the year
   if (Day_of_Month(doy) == Day_in_month(Month(doy))) then	!on the end of month
       tmp_air_month_ave = &
       sum( tmp_air_RunningRecord(1:Day_in_month(Month(doy))) ) / real( Day_in_month(Month(doy)) )
       
       if (Month(doy) == 1) then	!initialize on 31 January
           tmp_coldest_RunningRecord(2:20) = tmp_coldest_RunningRecord(1:19)
           tmp_coldest_RunningRecord(1)    = tmp_air_month_ave
           tmp_hottest_RunningRecord(2:20) = tmp_hottest_RunningRecord(1:19)
           tmp_hottest_RunningRecord(1)    = tmp_air_month_ave
       else
           tmp_coldest_RunningRecord(1) = min(tmp_coldest_RunningRecord(1), tmp_air_month_ave)
           tmp_hottest_RunningRecord(1) = max(tmp_hottest_RunningRecord(1), tmp_air_month_ave)
       end if
   end if
   
	!update 20 year running mean
	if (doy == Day_in_Year) then	!on the end of year
		if (Flag_spinup_read) then
			tmp_coldest_20yr_ave = sum( tmp_coldest_RunningRecord (1:20) ) / 20.0
			tmp_hottest_20yr_ave = sum( tmp_hottest_RunningRecord (1:20) ) / 20.0
		else
			tmp_coldest_20yr_ave = sum( tmp_coldest_RunningRecord (1:min(20,year)) ) /min(20,year)
			tmp_hottest_20yr_ave = sum( tmp_hottest_RunningRecord (1:min(20,year)) ) /min(20,year)
		endif
	end if
   
!annual sum of PAR intensity on establishment cells and grass cells
   if (doy==1) then
      sum_par_floor(:,:) = 0.0
   endif
   
!   do i=1, Dived !!!>>>>>>>>>>>>>>>TN: rm
!   do j=1, Dived !!!>>>>>>>>>>>>>>>TN: rm
   do i=1, GRID%N_x !!!<<<<<<<<<<<<<<<TN: add
   do j=1, GRID%N_y !!!<<<<<<<<<<<<<<<TN: add
      sum_par_floor(i,j) = sum_par_floor(i,j) + par_floor_rel(i,j) * par
   enddo
   enddo
   
!Growth degree day
   if (doy==1) then
      gdd_20yr_RunningRecord(2:20) = gdd_20yr_RunningRecord(1:19)
      gdd_20yr_RunningRecord(1)    = gdd5
      
      if (Flag_spinup_read) then
         gdd_20yr_ave = sum( gdd_20yr_RunningRecord(1:20) ) / 20.0
      else
         gdd_20yr_ave = sum( gdd_20yr_RunningRecord(1:min(20, year)) ) / real(min(20,year))
      endif
      
      gdd0         = 0.0
      gdd5         = 0.0
   end if
   gdd0 = gdd0 + max(0.0, tmp_air    )
   gdd5 = gdd5 + max(0.0, tmp_air-5.0)
   
!Temperature recorder for growth phase of grass
   if (doy==1) then
      if (pft_exist(C3g_no)) then
         p=C3g_no
      else
         p=C4g_no
      endif
     
      x = 0.0                               !reset
      count = 0                             !reset
      do i=1, Day_in_Year
         if (phenology_RunningRecord(i,p)) then
            x     = x + tmp_air_RunningRecord(i)
            count = count + 1
         endif
      end do
      
      tmp_ave_GrassGrowth_RR(2:20) = tmp_ave_GrassGrowth_RR(1:19)
      tmp_ave_GrassGrowth_RR(1)    = x / max(1, count)
      
   endif
   
!stat_water
   do p = 1, PFT_no
      do i = Day_in_Year, 2, -1
         stat_water_RunningRecord(i,p) = stat_water_RunningRecord(i-1,p)
      end do
      stat_water_RunningRecord(1,p) = stat_water(p)
   end do
   
END SUBROUTINE stat_climate



!*************************************************************************************************
! Carbon statistics recorder
!*************************************************************************************************
SUBROUTINE stat_carbon ()

!_____________ Set variables
   !Namespace
   USE data_structure
   USE vegi_status_current1
   USE vegi_status_current2
   USE grid_status_current1
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none
   
   !Local variables
   integer no, p, i
   real x, y, z
   
!_____________ Update of Ecosystem Total Carbon Pool running record
   x = pool_litter_trunk + pool_litter_leaf + pool_litter_root + &
       pool_litter_ag + pool_litter_bg + pool_som_int + pool_som_slow
   
   y = 0.0
   do no=1, Max_no
   if (tree_exist(no)) then
      y = y + mass_leaf(no) + mass_trunk(no) + mass_root(no) + mass_stock(no) + mass_available(no)
   endif
   end do
   
   z = sum(gmass_leaf(:,:)) + sum(gmass_root(:,:)) + &
       sum(gmass_available(:,:)) + sum(gmass_stock(:,:))
   
   pool_c_RR (2:Day_in_Year) = pool_c_RR (1:Day_in_Year-1) 
   
!   pool_c_RR (1) = (x+y+z) * C_in_drymass / Max_loc / Max_loc /1000.0  !!!>>>>>>>>>>>>TN:rm
   pool_c_RR (1) = (x+y+z) * C_in_drymass / real(GRID%Area) /1000.0 !!!<<<<<<<<<<<<TN:add
   !                        Unit conversion: (g DM / forest) -> (Kg C / m2)
   
!_____________ Carbon Fluxes
   !Update of NPP and GPP running record
   do p = 1, PFT_no
   do i = Day_in_Year, 2, -1
      npp_RunningRecord(i,p) = npp_RunningRecord(i-1,p)
      gpp_RunningRecord(i,p) = gpp_RunningRecord(i-1,p)
   end do
   end do
   
   do p = 1, PFT_no
      npp_RunningRecord(1,p) = npp(p)
      gpp_RunningRecord(1,p) = gpp(p)
   end do
   
   !Initializing NPP and GPP daily-sumup-variables
   gpp_daily_ind (:) = 0.0
   gpp           (:) = 0.0
   npp           (:) = 0.0
   
   !Carbon flux for each ecological compartment (reset at the beggining of each day)
   flux_c_uptake_RR(2:Day_in_Year) = flux_c_uptake_RR(1:Day_in_Year-1) ; flux_c_uptake_RR(1) = 0.0
   flux_c_mnt_RR   (2:Day_in_Year) = flux_c_mnt_RR   (1:Day_in_Year-1) ; flux_c_mnt_RR   (1) = 0.0
   flux_c_gro_RR   (2:Day_in_Year) = flux_c_gro_RR   (1:Day_in_Year-1) ; flux_c_gro_RR   (1) = 0.0
   flux_c_htr_RR   (2:Day_in_Year) = flux_c_htr_RR   (1:Day_in_Year-1) ; flux_c_htr_RR   (1) = 0.0
   flux_c_fir_RR   (2:Day_in_Year) = flux_c_fir_RR   (1:Day_in_Year-1) ; flux_c_fir_RR   (1) = 0.0
   flux_c_lit_RR   (2:Day_in_Year) = flux_c_lit_RR   (1:Day_in_Year-1) ; flux_c_lit_RR   (1) = 0.0
   flux_c_som_RR   (2:Day_in_Year) = flux_c_som_RR   (1:Day_in_Year-1) ; flux_c_som_RR   (1) = 0.0
   
   !Carbon flux for each organ (reset at the beggining of each day)
   resp_trunk(:) = 0.0
   resp_leaf (:) = 0.0
   resp_root (:) = 0.0
   resp_grass_ag = 0.0
   resp_grass_bg = 0.0
   
   !Litter flux (reset at the beggining of each day)
   flux_litter_trunk = 0.0
   flux_litter_leaf  = 0.0
   flux_litter_root  = 0.0
   flux_litter_ag    = 0.0
   flux_litter_bg    = 0.0
   
END SUBROUTINE stat_carbon



!*************************************************************************************************
! Vegetation statistics recorder
!*************************************************************************************************
SUBROUTINE stat_vegetation ()

!_____________ Set variables
!Namespace
   USE time_counter
   USE data_structure
   USE vegi_status_current1
   USE vegi_status_current2
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none
   
!Local variables
   integer                 no, p, i
   real,dimension(PFT_no)::sumup
   real,dimension(PFT_no)::lai_each
   
!_____________ Update the total biomass pool recorder, mass_sum_RunningRecord
If (doy==150) then
   
   sumup(:) = 0.0
   
   !woody PFTs
   Do no = 1, Max_no
   if (tree_exist(no)) then
      p = pft(no)
      sumup(p) = sumup(p) + mass_leaf(no) + mass_trunk(no) + mass_root(no) &
               + mass_available(no) + mass_stock(no)
   end if
   End do
   
   !herbaceous PFTs
   if (pft_exist(C3g_no)) then; p=C3g_no
   else                       ; p=C4g_no
   endif
   
   sumup(p) = &
   sum(gmass_leaf(:,:)) +sum(gmass_root(:,:)) +sum(gmass_available(:,:)) +sum(gmass_stock(:,:))
   
   !unit conversion ( gDM/forest -> gDM/m2 )
!   sumup(:) = sumup(:) / Max_loc / Max_loc !!!>>>>>>>>>>>>>>>TN: rm
   sumup(:) = sumup(:) / real(GRID%Area) !!!<<<<<<<<<<<<<<<TN: add
   
   !update running record
   do p=1, PFT_no
      do i=20,2,-1
      mass_sum_RunningRecord (i,p) = mass_sum_RunningRecord (i-1,p) 
      end do
      mass_sum_RunningRecord (1,p) = sumup(p)
   end do
   
Endif
   
!_____________ Update LAI records
   !Calculate LAI for each woody PFT
   lai_each(:) = 0.0
   do no=1, Max_no
   if (.not. tree_exist(no)) cycle
      p      = pft(no)
      lai_each(p) = lai_each(p) + la(no)
   enddo
   
   do p = 1, PFT_no
!      lai_each(p) = lai_each(p) / Max_loc / Max_loc !!!>>>>>>>>>>>>>>>TN: rm
      lai_each(p) = lai_each(p) / real(GRID%Area) !!!<<<<<<<<<<<<<<<TN: add
   end do
   
   !Calculate LAI for grass PFTs
   if (pft_exist(C3g_no)) then
      p=C3g_no
   else
      p=C4g_no
   endif
   
!   lai_each(p) = sum(lai_grass(:,:)) / DivedG / DivedG !!!>>>>>>>>>>>>>>>TN: rm
   lai_each(p) = sum(lai_grass(:,:)) / real(GRID%N_tot) !!!<<<<<<<<<<<<<<<TN: add
   
   !Update of LAI running record
   do p = 1, PFT_no
   do i = Day_in_Year, 2, -1
      lai_RunningRecord(i,p) = lai_RunningRecord(i-1,p)
   end do
   end do
   
   do p = 1, PFT_no
      lai_RunningRecord(1,p) = lai_each(p)
   end do
   
!_____________ Leaf phenology
   do p = 1, PFT_no
      do i = Day_in_Year, 2, -1
         phenology_RunningRecord(i,p) = phenology_RunningRecord(i-1,p)
      enddo
      phenology_RunningRecord(1,p) = phenology(p)
   enddo
   
END SUBROUTINE stat_vegetation



!**************************************************************************************************
! Search existing Plant Functiona Type in the grid
!**************************************************************************************************
SUBROUTINE pft_present ()

!_____________ Set variables
!Namespace
   USE data_structure
   USE vegi_status_current1
   USE vegi_status_current2
   implicit none
   
!Local variables
   integer	p, no	!loop counter
   
!_____________ Reset varialbes
   !p: save PFT number of available grass PFT
   if (pft_exist(C3g_no)) then
      p=C3g_no
   else
      p=C4g_no
   endif
   
   !reset and restore
   pft_exist(:)      = .false. !reset flag of PFT existence
   pft_exist(p)      = .true.  !restore available grass PFT
   
!_____________ Servey for existance of woody PFTs
   do no=1, Max_no
   if (tree_exist(no)) then
      pft_exist(pft(no)) = .true. 
   end if
   end do
   
END SUBROUTINE pft_present



!*************************************************************************************************
! Biome detemination routine (Must be called @ doy=Day_in_Year)
!
!   Criteria was taken from BIOME3 of A. Haxeltine & I.C. Prentice (1996) with some modifications.
!   Table B9 in Sato et al. (2007) is an abstract of this criteria.
!*************************************************************************************************
SUBROUTINE biome_determine ()

!_____________ Set variables
!Namespace
   USE data_structure
   USE vegi_status_current1
   USE vegi_status_current2
   USE grid_status_current1
   USE grid_status_current2
   implicit none
   
!Local variables
   integer dominant1, dominant2 !
   real    lai_max              !
   integer p, doy               !for general usage
   real    x                    !for general usage
   
!_____________ Preparare indices
!lai_max
   lai_max = 0.0
   do doy = 1, Day_in_Year
      x = sum(lai_RunningRecord(doy,:)) !LAI of this day
      lai_max = max(lai_max, x)
   end do
   
!First dominant PFT
   dominant1 = 1   !initialize
   x         = 0.0 !initialize
   do p = 1, PFT_no
      if ( sum(npp_RunningRecord(:,p)) > x ) then
         x         = sum(npp_RunningRecord (:,p))
         dominant1 = p
      endif
   enddo
   
!Secound dominant PFT
   dominant2 = 1   !initialize
   x         = 0.0 !initialize
   do p = 1, PFT_no
      if ( p .ne. dominant1               ) then
      if (sum(npp_RunningRecord(:,p)) > x ) then
         x         = sum(npp_RunningRecord (:,p))
         dominant2 = p
      endif
      endif
   enddo
   
!_____________ Determine biome type

!----------- Group 1 -----------
   !1: Polar desert
   if (gdd0 < 150) then
      biome = 1; return
   endif
   
!----------- Group 2 -----------
   !2: Arctic/Alpine-tundra
   if (gdd5 < 370) then
      biome = 2; return
   endif
   
!----------- Group 3 -----------
SELECT CASE (dominant1)
   !3: tropical evergreen forest
   CASE (1:5)
      if (lai_max>=2.5) then
            biome = 3; return
      endif
   
   !4: tropical deciduous forest
   CASE (6)
      if (lai_max>=2.5) then
         biome = 4; return
      endif
   
   !5: temperate conifer forest
   CASE (7)
      if (lai_max>=1.5) then
         biome = 5; return
      endif
   
   !6: temperate broad-leaved evergreen forest
   CASE (8)
      if (lai_max>=3.0) then
         biome = 6; return
      endif
   
   !7: temperate deciduous forest
   CASE (9)
      if (lai_max>=2.5 ) then
         biome = 7; return
      endif
   
   !8: boreal evergreen forest
   CASE (10:12)
         biome = 8; return
   
   !9: boreal deciduous forest
   CASE (13:14)
         biome = 9; return
   
END SELECT
   
!----------- Group 4 -----------
   !10: xeric woodland / scrub
   if (dominant1==1 .or. dominant1==2 .or. dominant1==3 .or. &
       dominant1==4 .or. dominant1==5 .or. dominant1==6 .or. dominant1==8) then
   if (lai_max>=1.0 ) then
      biome = 10; return
   endif
   endif
   
   if (dominant1==7 .or. dominant1==8 .or. dominant1==10 .or. &
       dominant1==11 .or. dominant1==12 .or. dominant1==13 .or. dominant1==14) then
   if (lai_max>=1.5 ) then
      biome = 10; return
   endif
   endif
   
!----------- Group 5 -----------
   !11: Grassland / steppe / Savanna
   !12: Desert
   if ( lai_max>=0.2 ) then
      biome = 11
   else
      biome = 12
   endif
   
END SUBROUTINE biome_determine



!*************************************************************************************************
! Output spinup files
!*************************************************************************************************
SUBROUTINE spinup_out (Fn, NSOLD, CMC, SNOWH, SNEQV, T1, STC, SMC, SH2O)

!_____________ Set variables
!Namespace
	USE data_structure
	USE time_counter
	USE vegi_status_current1
	USE grid_status_current1
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
	USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
	implicit none
	
!Arguments
	integer,intent(IN)::Fn !File I/O number
	
	!variables for NOAH-LSM
	INTEGER,intent(IN)::NSOLD       !NumSoil !Max number of soil layer
	REAL   ,intent(IN)::CMC         !Canopy water content (m)
	REAL   ,intent(IN)::SNOWH       !Actual snow depth (m)
	REAL   ,intent(IN)::SNEQV       !Water equiv snow depth (m)
	REAL   ,intent(IN)::T1          !Initial skin temperature (K)
	REAL   ,intent(IN)::STC(NSOLD)  !SOIL TEMP (K)
	REAL   ,intent(IN)::SMC(NSOLD)  !TOTAL SOIL MOISTURE CONTENT (VOLUMETRIC FRACTION)
	REAL   ,intent(IN)::SH2O(NSOLD) !UNFROZEN SOIL MOISTURE CONTENT (VOLUMETRIC FRACTION)
	
!Local variables
	integer no, p, i, j
	integer pft_available_max
	integer pft_available (1:PFT_no)
   
!_____________ Main part1 (variables of vegi_status_current1)
	!These two variables must be read at the beggining
	write (Fn) &
	pft_exist      (:), &
	tree_exist     (:)   
	
	!List up available PFTs
	i = 0
	do p = 1, PFT_no
		if (pft_exist(p)) then
			i = i + 1
			pft_available(i) = p
		endif
	enddo
	pft_available_max = i
	
	write (Fn) &
	Spinup_year+Simulation_year, & !Variable 'Spinup_year' for the continued simulation
	biome                          
	
	if (pft_available_max >= 1) then
		do p=1, pft_available_max
			write (Fn)                         &
			phenology      (pft_available(p)), &
			dfl_leaf_onset (pft_available(p)), &
			dfl_leaf_shed  (pft_available(p)), &
			stat_water     (pft_available(p))   
		enddo
	endif
	
	Do no=1, Max_no
	if (tree_exist(no)) then
		write (Fn) &
		pft            (no), &
		age            (no), &
		height         (no), &
		bole           (no), &
		height_limit   (no), &
		flag_suppress  (no), &
		dbh_heartwood  (no), &
		dbh_sapwood    (no), &
		crown_diameter (no), &
		crown_area     (no), &
		bole_x         (no), &
		bole_y         (no), &
		crown_x        (no), &
		crown_y        (no), &
		radius_limit   (no), &
		la             (no), &
		mass_leaf      (no), &
		mass_trunk     (no), &
		mass_root      (no), &
		mass_stock     (no), &
		mass_available (no), &
		mort_regu1     (no), &
		mort_regu2     (no), &
		mort_regu4     (no), &
!		MSR_plant      (no), &	! R-SEIB
!		mort_yearly_total (no), &	! R-SEIB
		npp_crowntop   (no), &
		npp_crownbottom(no,:) 
	endif
	End do
	
!	do i=1, DivedG !!!>>>>>>>>>>>>>>>TN: rm
	do i=1, GRID%N_x !!!<<<<<<<<<<<<<<<TN: add
		write (Fn)             &
		gmass_leaf      (i,:), &
		gmass_root      (i,:), &
		gmass_available (i,:), &
		gmass_stock     (i,:), &
		lai_grass       (i,:)   
    enddo
	
!	do i=1, DivedG    !!!>>>>>>>>>>>>>>>TN: rm
!    do j=1, DivedG !!!>>>>>>>>>>>>>>>TN: rm
	do i=1, GRID%N_x    !!!<<<<<<<<<<<<<<<TN: add
    do j=1, GRID%N_y  !!!<<<<<<<<<<<<<<<TN: add
		write (Fn) lai_opt_grass_RunningRecord(:,i,j)
	enddo
	enddo
	
	do p=1, pft_available_max
	do i=1, Day_in_Year
		write (Fn)                                     &
		phenology_RunningRecord (i, pft_available(p)), &
		npp_RunningRecord       (i, pft_available(p)), &
		gpp_RunningRecord       (i, pft_available(p)), &
		stat_water_RunningRecord(i, pft_available(p)), &
		lai_RunningRecord       (i, pft_available(p))   
	enddo
	enddo
	
	do p=1, pft_available_max
	do i=1, 20
		write (Fn) mass_sum_RunningRecord(i, pft_available(p))
	enddo
	enddo
	
!	! R-SEIB ***********************************
!	do p=1,  pft_available_max
!		write (Fn) &
!			assim_rate0(pft_available(p),:), &
!			gtc_rate0(pft_available(p),:), &
!			assim_rate1(pft_available(p),:), &
!			gtc_rate1(pft_available(p),:), &
!			assim_rate2(pft_available(p),:), &
!			gtc_rate2(pft_available(p),:), &
!			assim_rate3(pft_available(p),:), &
!			gtc_rate3(pft_available(p),:), &
!			assim_rate4(pft_available(p),:), &
!			gtc_rate4(pft_available(p),:), &
!			assim_rate5(pft_available(p),:), &
!			gtc_rate5(pft_available(p),:)
!	enddo
!	
!	do i=1, DivedG
!		write (Fn)             &
!			asim_sat     (i,:), &
!			asimlue      (i,:)
!    enddo
!	
!	do no=1, Max_no
!		if (tree_exist(no)) then
!			write (Fn) &
!				G_plant        (no), &
!				G_p_pre        (no), &
!				HF_plant       (no)
!		endif
!	enddo
!	! R-SEIB end *******************************
	
!_____________ Main part2 (variables of grid_status_current1)
	write (Fn) &
	dfl_fire            , &
	fire_number         , &
	tmp_coldest_20yr_ave, &
	tmp_hottest_20yr_ave, &
	gdd_20yr_ave        , &
	pool_litter_trunk   , &
	pool_litter_leaf    , &
	pool_litter_root    , &
	pool_litter_ag      , &
	pool_litter_bg      , &
	pool_som_int        , &
	pool_som_slow       , &
	pool_fuel_standT    , &
	pool_fuel_standG    , &
	pool_w(:)           , &
	pool_snow              
	
	write (Fn) &
	gdd_20yr_RunningRecord    (:), &
	tmp_coldest_RunningRecord (:), &
	tmp_hottest_RunningRecord (:), &
	tmp_ave_GrassGrowth_RR    (:)   
	
	write (Fn) &
	pool_c_RR            (:), &
	flux_c_uptake_RR     (:), &
	flux_c_mnt_RR        (:), &
	flux_c_gro_RR        (:), &
	flux_c_htr_RR        (:), &
	flux_c_fir_RR        (:), &
	flux_c_lit_RR        (:), &
	flux_c_som_RR        (:), &
    flux_ro1_RunningRecord (:), &
    flux_ro2_RunningRecord (:), &
    flux_ic_RunningRecord  (:), &
    flux_ev_RunningRecord  (:), &
    flux_tr_RunningRecord  (:), &
    flux_sl_RunningRecord  (:), &
    flux_tw_RunningRecord  (:), &
    flux_sn_RunningRecord  (:)   
	
	write (Fn) &
	tmp_air_RunningRecord (:)
!	tmp_0m_RunningRecord (:), &	! R-SEIB
!	tmp_soil_RunningRecord(:)   
	
	do i=1, NumSoil
		write (Fn) tmp_soil_RunningRecord(:,i)
	enddo
	
	write (Fn)                 &
	prec_RunningRecord    (:), &
	ev_pot_RunningRecord  (:), &
	par_RunningRecord     (:), &
	pool_w1_RunningRecord (:)   
!	soil_theta_RunningRecord (:)	! R-SEIB
	
!_____________ Main part3 (variables for NOAH-LSM)
	write (Fn) &
	CMC    , &
	SNOWH  , &
	SNEQV  , &
	T1     , &
	STC (:), &
	SMC (:), &
	SH2O(:)   

END SUBROUTINE spinup_out



!*************************************************************************************************
! Input spinup files
!*************************************************************************************************
SUBROUTINE spinup_in (Fn, NSOLD, CMC, SNOWH, SNEQV, T1, STC, SMC, SH2O)

!_____________ Set variables
!Namespace
	USE data_structure
	USE time_counter
	USE vegi_status_current1
	USE grid_status_current1
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
	USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
	implicit none
	
!Arguments
	integer,intent(IN)::Fn !File I/O number
	
	!variables for NOAH-LSM
	INTEGER,intent(IN) ::NSOLD       !NumSoil !Max number of soil layer
	REAL   ,intent(OUT)::CMC         !Canopy water content (m)
	REAL   ,intent(OUT)::SNOWH       !Actual snow depth (m)
	REAL   ,intent(OUT)::SNEQV       !Water equiv snow depth (m)
	REAL   ,intent(OUT)::T1          !Initial skin temperature (K)
	REAL   ,intent(OUT)::STC(NSOLD)  !SOIL TEMP (K)
	REAL   ,intent(OUT)::SMC(NSOLD)  !TOTAL SOIL MOISTURE CONTENT (VOLUMETRIC FRACTION)
	REAL   ,intent(OUT)::SH2O(NSOLD) !UNFROZEN SOIL MOISTURE CONTENT (VOLUMETRIC FRACTION)
	
!Local variables
	integer no, p, i, j
	integer pft_available_max
	integer pft_available (1:PFT_no)
	
!_____________ Main part1 (variables of vegi_status_current1)
	!These two variables must be read at the beggining
	read (Fn) &
	pft_exist      (:), &
	tree_exist     (:)   
	
	!List up available PFTs
	i = 0
	do p = 1, PFT_no
		if (pft_exist(p)) then
			i = i + 1
			pft_available(i) = p
		endif
	enddo
	pft_available_max = i
	
	read (Fn) &
	Spinup_year, & !Variable 'Spinup_year' for the continued simulation
	biome          
	
	if (pft_available_max >= 1) then
		do p=1, pft_available_max
			read (Fn)                         &
			phenology      (pft_available(p)), &
			dfl_leaf_onset (pft_available(p)), &
			dfl_leaf_shed  (pft_available(p)), &
			stat_water     (pft_available(p))   
		enddo
	endif
	
	Do no=1, Max_no
	if (tree_exist(no)) then
		read (Fn) &
		pft            (no), &
		age            (no), &
		height         (no), &
		bole           (no), &
		height_limit   (no), &
		flag_suppress  (no), &
		dbh_heartwood  (no), &
		dbh_sapwood    (no), &
		crown_diameter (no), &
		crown_area     (no), &
		bole_x         (no), &
		bole_y         (no), &
		crown_x        (no), &
		crown_y        (no), &
		radius_limit   (no), &
		la             (no), &
		mass_leaf      (no), &
		mass_trunk     (no), &
		mass_root      (no), &
		mass_stock     (no), &
		mass_available (no), &
		mort_regu1     (no), &
		mort_regu2     (no), &
		mort_regu4     (no), &
!		MSR_plant      (no), &	! R-SEIB
!		mort_yearly_total (no), &	! R-SEIB
		npp_crowntop   (no), &
		npp_crownbottom(no,:) 
	endif
	End do
	
!	do i=1, DivedG !!!>>>>>>>>>>>>>>>TN: rm
	do i=1, GRID%N_x !!!<<<<<<<<<<<<<<<TN: add
		read (Fn)             &
		gmass_leaf      (i,:), &
		gmass_root      (i,:), &
		gmass_available (i,:), &
		gmass_stock     (i,:), &
		lai_grass       (i,:)   
    enddo
	
!	do i=1, DivedG    !!!>>>>>>>>>>>>>>>TN: rm
!    do j=1, DivedG !!!>>>>>>>>>>>>>>>TN: rm
	do i=1, GRID%N_x    !!!<<<<<<<<<<<<<<<TN: add
    do j=1, GRID%N_y  !!!<<<<<<<<<<<<<<<TN: add
		read (Fn) lai_opt_grass_RunningRecord(:,i,j)
	enddo
	enddo
	
	do p=1, pft_available_max
	do i=1, Day_in_Year
		read (Fn)                                      &
		phenology_RunningRecord (i, pft_available(p)), &
		npp_RunningRecord       (i, pft_available(p)), &
		gpp_RunningRecord       (i, pft_available(p)), &
		stat_water_RunningRecord(i, pft_available(p)), &
		lai_RunningRecord       (i, pft_available(p))   
	enddo
	enddo
   
	do p=1, pft_available_max
	do i=1, 20
		read (Fn) mass_sum_RunningRecord(i, pft_available(p))
	enddo
	enddo
	
!	! R-SEIB
!	do p=1,  pft_available_max
!		write (Fn) &
!			assim_rate0(pft_available(p),:), &
!			gtc_rate0(pft_available(p),:), &
!			assim_rate1(pft_available(p),:), &
!			gtc_rate1(pft_available(p),:), &
!			assim_rate2(pft_available(p),:), &
!			gtc_rate2(pft_available(p),:), &
!			assim_rate3(pft_available(p),:), &
!			gtc_rate3(pft_available(p),:), &
!			assim_rate4(pft_available(p),:), &
!			gtc_rate4(pft_available(p),:), &
!			assim_rate5(pft_available(p),:), &
!			gtc_rate5(pft_available(p),:)
!	enddo
!	
!	do i=1, DivedG
!		write (Fn)             &
!			asim_sat     (i,:), &
!			asimlue      (i,:)
!    enddo
!	
!	do no=1, Max_no
!		if (tree_exist(no)) then
!			write (Fn) &
!				G_plant        (no), &
!				G_p_pre        (no), &
!				HF_plant       (no)
!		endif
!	enddo
!	
!_____________ Main part2 (variables of grid_status_current1)
	read (Fn)             &
	dfl_fire            , &
	fire_number         , &
	tmp_coldest_20yr_ave, &
	tmp_hottest_20yr_ave, &
	gdd_20yr_ave        , &
	pool_litter_trunk   , &
	pool_litter_leaf    , &
	pool_litter_root    , &
	pool_litter_ag      , &
	pool_litter_bg      , &
	pool_som_int        , &
	pool_som_slow       , &
	pool_fuel_standT    , &
	pool_fuel_standG    , &
	pool_w(:)           , &
	pool_snow              
	
	read (Fn)                      &
	gdd_20yr_RunningRecord    (:), &
	tmp_coldest_RunningRecord (:), &
	tmp_hottest_RunningRecord (:), &
	tmp_ave_GrassGrowth_RR    (:)   
	
	read (Fn)                 &
	pool_c_RR            (:), &
	flux_c_uptake_RR     (:), &
	flux_c_mnt_RR        (:), &
	flux_c_gro_RR        (:), &
	flux_c_htr_RR        (:), &
	flux_c_fir_RR        (:), &
	flux_c_lit_RR        (:), &
	flux_c_som_RR        (:), &
    flux_ro1_RunningRecord (:), &
    flux_ro2_RunningRecord (:), &
    flux_ic_RunningRecord  (:), &
    flux_ev_RunningRecord  (:), &
    flux_tr_RunningRecord  (:), &
    flux_sl_RunningRecord  (:), &
    flux_tw_RunningRecord  (:), &
    flux_sn_RunningRecord  (:)   
	
	read (Fn) &
	tmp_air_RunningRecord (:)
!	tmp_0m_RunningRecord (:), &	! R-SEIB
!	tmp_soil_RunningRecord(:)   
	
	do i=1, NumSoil
		read (Fn) tmp_soil_RunningRecord(:,i)
	enddo
	
	read (Fn)                  &
	prec_RunningRecord    (:), &
	ev_pot_RunningRecord  (:), &
	par_RunningRecord     (:), &
	pool_w1_RunningRecord (:)   
!	soil_theta_RunningRecord (:)	! R-SEIB
	
!_____________ Main part3 (variables for NOAH-LSM)
	read (Fn) &
	CMC    , &
	SNOWH  , &
	SNEQV  , &
	T1     , &
	STC (:), &
	SMC (:), &
	SH2O(:)   
	
END SUBROUTINE spinup_in



!***********************************************************************************************
! Climate statistics recorder
!***********************************************************************************************
SUBROUTINE tmp_soil_interpolate (tmp_soil1, tmp_soil2, tmp_soil3, tmp_soil)
   USE data_structure
   
!_____________ Set variables
!Augments
   real,intent(IN)   ::tmp_soil1         !soil temperature @ top layer    (Celcius)
   real,intent(IN)   ::tmp_soil2         !soil temperature @ 20th layer   (Celcius)
   real,intent(IN)   ::tmp_soil3         !soil temperature @ bottom layer (Celcius)
   real              ::tmp_soil(NumSoil) !soil temperatures from top of the soil surface (Celcius)
   
!Local variables
   real    x
   integer i
   
!_____________ Main part
   do i=1, min(20, NumSoil)
     x = real(i) / 20.
     tmp_soil(i) = (1.0-x)*tmp_soil1 + x*tmp_soil2  
   enddo
   
   if (NumSoil>20) then
      do i=21, min(30, NumSoil)
         x = (real(i)-20.0) / 10.0
         tmp_soil(i) = (1.0-x)*tmp_soil2 + x*tmp_soil3  
      enddo
   endif
   
END SUBROUTINE tmp_soil_interpolate

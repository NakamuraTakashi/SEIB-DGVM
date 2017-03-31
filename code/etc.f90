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
      tmp_coldest_20yr_ave = sum( tmp_coldest_RunningRecord (1:min(20,year)) ) /min(20, year)
      tmp_hottest_20yr_ave = sum( tmp_hottest_RunningRecord (1:min(20,year)) ) /min(20, year)
   end if
   
!annual sum of PAR intensity on establishment cells and grass cells
   if (doy==1) then
      sum_par_floor(:,:) = 0.0
   endif
   
   do i=1, Dived
   do j=1, Dived
      sum_par_floor(i,j) = sum_par_floor(i,j) + par_floor_rel(i,j) * par
   enddo
   enddo
   
!Growth degree day
   if (doy==1) then
      gdd_20yr_RunningRecord(2:20) = gdd_20yr_RunningRecord(1:19)
      gdd_20yr_RunningRecord(1)    = gdd5
      gdd_20yr_ave = sum( gdd_20yr_RunningRecord(1:min(20, year)) ) / min(20, year)
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
   
   pool_c_RR (1) = (x+y+z) * C_in_drymass / Max_loc / Max_loc /1000.0
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
   sumup(:) = sumup(:) / Max_loc / Max_loc
   
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
      lai_each(p) = lai_each(p) / Max_loc / Max_loc
   end do
   
   !Calculate LAI for grass PFTs
   if (pft_exist(C3g_no)) then
      p=C3g_no
   else
      p=C4g_no
   endif
   
   lai_each(p) = sum(lai_grass(:,:)) / DivedG / DivedG
   
   !Update of LAI running record
   do p = 1, PFT_no
   do i = Day_in_Year, 2, -1
      lai_RunningRecord(i,p) = lai_RunningRecord(i-1,p)
   end do
   end do
   
   do p = 1, PFT_no
      lai_RunningRecord(1,p) = lai_each(p)
   end do
   
   !Update total LAI
   lai = sum(lai_each(:))
   
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
   p = dominant1
   !10: xeric woodland / scrub
   if (p==1 .or. p==2 .or. p==3 .or. p==4 .or. p==5 .or. p==6 .or. p==8) then
   if (lai_max>=1.0 ) then
      biome = 10; return
   endif
   endif
   
   if (p==7 .or. p==8 .or. p==10 .or. p==11 .or. p==12 .or. p==13 .or. p==14) then
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
! output spinup files
!*************************************************************************************************
SUBROUTINE spinup_out (Fn)
!_____________ Set variables
!Namespace
   USE data_structure
   USE time_counter
   USE vegi_status_current1
   USE grid_status_current1
   implicit none
   
!Arguments
   integer,intent(IN)::Fn !File I/O number
   
!Local variables
   integer           n ,i, j, k
   character(len=3)  string_pft, string_year, string_dived, string_divedG
   
!_____________ preparation for strings
i = int(  Day_in_Year       /100 ) ; string_year  (1:1)= char(i+48)
j = int( (Day_in_Year-i*100)/ 10 ) ; string_year  (2:2)= char(j+48)
k =    (  Day_in_Year-i*100-j*10 ) ; string_year  (3:3)= char(k+48)

i = int(  PFT_no            /100 ) ; string_pft   (1:1)= char(i+48)
j = int( (PFT_no     -i*100)/ 10 ) ; string_pft   (2:2)= char(j+48)
k = int(  PFT_no     -i*100-j*10 ) ; string_pft   (3:3)= char(k+48)

i = int(  Dived             /100 ) ; string_dived (1:1)= char(i+48)
j = int( (Dived      -i*100)/ 10 ) ; string_dived (2:2)= char(j+48)
k = int(  Dived      -i*100-j*10 ) ; string_dived (3:3)= char(k+48)

i = int(  DivedG            /100 ) ; string_divedG(1:1)= char(i+48)
j = int( (DivedG     -i*100)/ 10 ) ; string_divedG(2:2)= char(j+48)
k = int(  DivedG     -i*100-j*10 ) ; string_divedG(3:3)= char(k+48)

!_____________ Main part1 (variables of vegi_status_current1)
write (Fn,*) Spinup_year+Simulation_year  !Variable 'Spinup_year' for the continued simulation
write (Fn,*) tree_exist(:)                !This variable must be read at the beggining

write (Fn,'(i5)')       biome
write (Fn,'(f10.5,1x)') lai

write (Fn,*)                             phenology      (:)
write (Fn,*)                             pft_exist      (:)
write (Fn,'('//string_pft//'(i7  ,1x))') dfl_leaf_onset (:)
write (Fn,'('//string_pft//'(i7  ,1x))') dfl_leaf_shed  (:)
write (Fn,'('//string_pft//'(f7.3,1x))') stat_water     (:)

Do n=1, Max_no
if (tree_exist(n)) then
   write (Fn,'(6(i5,1x) )'  ) pft(n), age(n), height(n), bole(n), height_limit(n), flag_suppress(n)
   write (Fn,'(4(f10.5,1x))') dbh_heartwood(n), dbh_sapwood(n), crown_diameter(n), crown_area(n)
   write (Fn,'(5(f10.5,1x))') bole_x(n), bole_y(n), crown_x(n), crown_y(n), radius_limit(n)
   write (Fn,'(6(f12.2,1x))') la(n), mass_leaf(n), mass_trunk(n), mass_root(n), mass_stock(n), mass_available(n)
   write (Fn,'(f12.3,1x, f10.5,1x, f10.8)') mort_regu1(n), mort_regu2(n), mort_regu4(n)
   write (Fn,'(11(f10.5,1x))') npp_crowntop(n), npp_crownbottom(n,:)
endif
End do

Do i=1, DivedG
   write (Fn,'('//string_divedG//'(f12.3,1x))') gmass_leaf      (i,:)
   write (Fn,'('//string_divedG//'(f12.3,1x))') gmass_root      (i,:)
   write (Fn,'('//string_divedG//'(f12.3,1x))') gmass_available (i,:)
   write (Fn,'('//string_divedG//'(f12.3,1x))') gmass_stock     (i,:)
   write (Fn,'('//string_divedG//'(f12.3,1x))') lai_grass       (i,:)
End do

Do i=1, DivedG
Do j=1, DivedG
   write (Fn,'(20(f5.2,1x))')  lai_opt_grass_RunningRecord(:,i,j)
End do
End do

Do i=1, Day_in_Year
   write (Fn, *)                             phenology_RunningRecord (i,:)
   write (Fn,'('//string_pft//'(f10.1,1x))') npp_RunningRecord       (i,:)
   write (Fn,'('//string_pft//'(f10.1,1x))') gpp_RunningRecord       (i,:)
   write (Fn,'('//string_pft//'(f7.3,1x))' ) stat_water_RunningRecord(i,:)
   write (Fn,'('//string_pft//'(f7.3,1x))' ) lai_RunningRecord       (i,:)
End do

Do i=1, 20
   write (Fn,'('//string_pft//'(f10.3,1x))') mass_sum_RunningRecord(i,:)
End do

!_____________ Main part2 (variables of vegi_status_current2)
write (Fn,'(2(i7,1x))'    ) dfl_fire, fire_number
write (Fn,'(3(f10.5,1x))' ) tmp_coldest_20yr_ave, tmp_hottest_20yr_ave, gdd_20yr_ave
write (Fn,'(3(f17.1,1x))' ) pool_litter_trunk, pool_litter_leaf, pool_litter_root
write (Fn,'(4(f17.1,1x))' ) pool_litter_ag, pool_litter_bg, pool_som_int, pool_som_slow
write (Fn,'(2(f17.1,1x))' ) pool_fuel_standT, pool_fuel_standG
write (Fn,'(31(f10.2,1x))') pool_w(:), pool_snow

write (Fn,'( 20(f7.2,1x))') gdd_20yr_RunningRecord    (:)
write (Fn,'( 20(f6.2,1x))') tmp_coldest_RunningRecord (:)
write (Fn,'( 20(f6.2,1x))') tmp_hottest_RunningRecord (:)
write (Fn,'( 20(f6.2,1x))') tmp_ave_GrassGrowth_RR    (:)

write (Fn,'('//string_year//'(f15.5,1x))') pool_c_RR (:)

write (Fn,'('//string_year//'(f15.1,1x))') flux_c_uptake_RR (:)
write (Fn,'('//string_year//'(f15.1,1x))') flux_c_mnt_RR    (:)
write (Fn,'('//string_year//'(f15.1,1x))') flux_c_gro_RR    (:)
write (Fn,'('//string_year//'(f15.1,1x))') flux_c_htr_RR    (:)
write (Fn,'('//string_year//'(f15.1,1x))') flux_c_fir_RR    (:)

write (Fn,'('//string_year//'(f10.5,1x))') flux_ro_RunningRecord(:)
write (Fn,'('//string_year//'(f10.5,1x))') flux_ic_RunningRecord(:)
write (Fn,'('//string_year//'(f10.5,1x))') flux_ev_RunningRecord(:)
write (Fn,'('//string_year//'(f10.5,1x))') flux_tr_RunningRecord(:)

write (Fn,'('//string_year//'(f6.2,1x))') tmp_air_RunningRecord (:)
   Do i=1, NumSoil
   write (Fn,'('//string_year//'(f6.2,1x))') tmp_soil_RunningRecord(:,i)
   End do
write (Fn,'('//string_year//'(f6.2,1x))') prec_RunningRecord    (:)
write (Fn,'('//string_year//'(f6.2,1x))') ev_pot_RunningRecord  (:)
write (Fn,'('//string_year//'(f6.1,1x))') par_RunningRecord     (:)
write (Fn,'('//string_year//'(f6.1,1x))') pool_w1_RunningRecord (:)

END SUBROUTINE spinup_out



!*************************************************************************************************
! output spinup files
!*************************************************************************************************
SUBROUTINE spinup_in (Fn)

!_____________ Set variables
!Namespace
   USE data_structure
   USE time_counter
   USE vegi_status_current1
   USE grid_status_current1
   implicit none
   
!Arguments
   integer,intent(IN)::Fn !File I/O number
   
!Local variables
   integer n ,i, j
   
!_____________ Main part1 (variables of vegi_status_current1)
read (Fn,*) Spinup_year
read (Fn,*) tree_exist(:)                !This variable must be read at the beggining

read (Fn,*) biome
read (Fn,*) lai

read (Fn,*) phenology      (:)
read (Fn,*) pft_exist      (:)
read (Fn,*) dfl_leaf_onset (:)
read (Fn,*) dfl_leaf_shed  (:)
read (Fn,*) stat_water     (:)

Do n=1, Max_no
if (tree_exist(n)) then
   read (Fn,*) pft(n), age(n), height(n), bole(n), height_limit(n), flag_suppress(n)
   read (Fn,*) dbh_heartwood(n), dbh_sapwood(n), crown_diameter(n), crown_area(n)
   read (Fn,*) bole_x(n), bole_y(n), crown_x(n), crown_y(n), radius_limit(n)
   read (Fn,*) la(n), mass_leaf(n), mass_trunk(n), mass_root(n), mass_stock(n), mass_available(n)
   read (Fn,*) mort_regu1(n), mort_regu2(n), mort_regu4(n)
   read (Fn,*) npp_crowntop(n), npp_crownbottom(n,:)
endif
End do

Do i=1, DivedG
   read (Fn,*) gmass_leaf      (i,:)
   read (Fn,*) gmass_root      (i,:)
   read (Fn,*) gmass_available (i,:)
   read (Fn,*) gmass_stock     (i,:)
   read (Fn,*) lai_grass       (i,:)
End do

Do i=1, DivedG
Do j=1, DivedG
   read (Fn,*) lai_opt_grass_RunningRecord(:,i,j)
End do
End do

Do i=1, Day_in_Year
   read (Fn,*) phenology_RunningRecord (i,:)
   read (Fn,*) npp_RunningRecord       (i,:)
   read (Fn,*) gpp_RunningRecord       (i,:)
   read (Fn,*) stat_water_RunningRecord(i,:)
   read (Fn,*) lai_RunningRecord       (i,:)
End do

Do i=1, 20
   read (Fn,*) mass_sum_RunningRecord(i,:)
End do

!_____________ Main part2 (variables of vegi_status_current2)
read (Fn,*) dfl_fire, fire_number
read (Fn,*) tmp_coldest_20yr_ave, tmp_hottest_20yr_ave, gdd_20yr_ave

read (Fn,*) pool_litter_trunk, pool_litter_leaf, pool_litter_root
read (Fn,*) pool_litter_ag, pool_litter_bg, pool_som_int, pool_som_slow
read (Fn,*) pool_fuel_standT, pool_fuel_standG
read (Fn,*) pool_w(:), pool_snow

read (Fn,*) gdd_20yr_RunningRecord    (:)
read (Fn,*) tmp_coldest_RunningRecord (:)
read (Fn,*) tmp_hottest_RunningRecord (:)
read (Fn,*) tmp_ave_GrassGrowth_RR    (:)

read (Fn,*) pool_c_RR (:)

read (Fn,*) flux_c_uptake_RR (:)
read (Fn,*) flux_c_mnt_RR    (:)
read (Fn,*) flux_c_gro_RR    (:)
read (Fn,*) flux_c_htr_RR    (:)
read (Fn,*) flux_c_fir_RR    (:)

read (Fn,*) flux_ro_RunningRecord(:)
read (Fn,*) flux_ic_RunningRecord(:)
read (Fn,*) flux_ev_RunningRecord(:)
read (Fn,*) flux_tr_RunningRecord(:)

read (Fn,*) tmp_air_RunningRecord (:)
   Do i=1, NumSoil
   read (Fn,*) tmp_soil_RunningRecord(:,i)
   End do
read (Fn,*) prec_RunningRecord    (:)
read (Fn,*) ev_pot_RunningRecord  (:)
read (Fn,*) par_RunningRecord     (:)
read (Fn,*) pool_w1_RunningRecord (:)

END SUBROUTINE spinup_in



!***********************************************************************************************
! Climate statistics recorder
!***********************************************************************************************
SUBROUTINE tmp_soil_interpolate (tmp_soil1, tmp_soil2, tmp_soil3, tmp_soil)
   USE data_structure
   implicit none
   
!_____________ Set variables
!Augments
   real,intent(IN) ::tmp_soil1         !soil temperature @ top layer    (Celcius)
   real,intent(IN) ::tmp_soil2         !soil temperature @ 20th layer   (Celcius)
   real,intent(IN) ::tmp_soil3         !soil temperature @ bottom layer (Celcius)
   real,intent(OUT)::tmp_soil(NumSoil) !soil temperatures from top of the soil surface (Celcius)
   
!Local variables
   real    x
   integer i
   
!_____________ Main part
   do i=1, 20
     x = real(i) / 20.0
     tmp_soil(i) = (1.0-x)*tmp_soil1 + x*tmp_soil2  
   enddo
    do i=21, NumSoil
      x = (real(i)-20.0) / 10.0
      tmp_soil(i) = (1.0-x)*tmp_soil2 + x*tmp_soil3  
    enddo
   
END SUBROUTINE tmp_soil_interpolate

!*************************************************************************************************
! Output for Viewer
!*************************************************************************************************
SUBROUTINE output_for_viewer &
   (Fn, LAT, LON, ALT, Albedo_soil0, W_fi, W_wilt, W_sat, &
   tmp_air, tmp_soil1, tmp_soil2, tmp_soil3, prec, rh, cloud, rad_short, rad_long, wind)
   
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
   
!Arguments
   !File I/O number
   integer,intent(IN)::Fn
   
   !Location data (single precision)
   real,intent(IN)::LAT          !latitude  (degree)
   real,intent(IN)::LON          !longitude  (degree)
   real,intent(IN)::ALT          !altitude  (m above MSL)
   real,intent(IN)::Albedo_soil0 !soil albedo
   real,intent(IN)::W_fi         !filed capacity  (m3/m3, 0.0 -> 1.0)
   real,intent(IN)::W_wilt       !wilting point   (m3/m3, 0.0 -> 1.0)
   real,intent(IN)::W_sat        !saturate point  (m3/m3, 0.0 -> 1.0)
   
   !Climatic data (single precision)
   real,intent(IN)::tmp_air    !2m air temperature                (Celcius) 
   real,intent(IN)::tmp_soil1  !Soil temperatures for layers 1~5   (Celcius) 
   real,intent(IN)::tmp_soil2  !Soil temperatures for layers 6~10  (Celcius) 
   real,intent(IN)::tmp_soil3  !Soil temperatures for layers 11~20 (Celcius) 
   real,intent(IN)::prec       !Precipitation                     (mm day-1)
   real,intent(IN)::rh         !Relative Humidity                 (%)
   real,intent(IN)::cloud      !Total cloudiness                  (fraction)
   real,intent(IN)::rad_short  !Shortwave radiation @ midday      (W m-2)
   real,intent(IN)::rad_long   !Daily mean of longwave radiation  (W m-2)
   real,intent(IN)::wind       !Wind velocity                     (m s-1)   
   
!Local variables
   !character(len=3)                  PFT_no_string
   real   ,dimension(PFT_no)       ::c_leaf, c_trunk, c_root, c_stock, c_available
   integer,dimension(PFT_no,0:41)::size_dist
   
   integer i, p, no, type_grass, count
   real    a1, a2, a3, a4, a5, a6, b1, b2, b3, c1, c2, x
   real    tree_density, Dist_step
   
!_____________ Define local parameter, set local variable
   !D.B.H. size step for size class output
   Dist_step = 2.5
   
   !type_grass: type of grass available 
   if (pft_exist(C3g_no)) then
      type_grass = 3
   else
      type_grass = 4
   endif
   
!!_____________ Write File header (Write Only Onece at the Begginign of the simulation)
If ( counter==1 .or.  (Flag_spinup_read .and. counter==1+Spinup_year*Day_in_Year) ) then
   
   !write out setting
   write (Fn, '( 3(f7.2,a), i5, a,i5, a,i4,  a,i4 )' )  &
!   LAT,',',LON,',',ALT,',',Simulation_year,',',PFT_no,',',Max_loc,',',NumSoil  !!!>>>>>>>>>>>>TN:rm
   LAT,',',LON,',',ALT,',',Simulation_year,',',PFT_no,',',GRID%N_x,',',NumSoil  !!!<<<<<<<<<<<<TN:add y方向の情報がない。要検討
   
   write (Fn, '( f7.5,a,f7.5,a,f7.5,a,f7.5 )' ) &
   Albedo_soil0,',',W_fi,',',W_wilt,',',W_sat
   
   do p=1, PFT_no; write (Fn,'(i1,a$)'  ) Life_type(p),     ','; end do
   write (Fn,*)
   
   do p=1, PFT_no; write (Fn,'(i1,a$)'  ) Phenology_type(p),','; end do
   write (Fn,*)
   
   do p=1, PFT_no; write (Fn,'(f6.4,a$)') PN_s(p),          ','; end do
   write (Fn,*)
   
   do p=1, PFT_no; write (Fn,'(f6.4,a$)') PN_r(p),          ','; end do
   write (Fn,*)
   
Endif

!_____________ Write Annual output
If ( doy==1 ) then
   !Flag for existence of PFT
   write (Fn,*)
   do p=1, PFT_no
      i = 0
      if (pft_exist(p)) i=1
      write (Fn,'(i1,a$)' ) i,','
   end do
   write (Fn,*)
   write (Fn,*)
   
Endif

!_____________ Write Daily output
!*** 1st line (Simulation day counter)
   write (Fn,'(i7)') &
   counter

!*** 2nd line (Climatic data)
   write (Fn,'( 4(f5.1,a), f4.2,a, f6.2,a, f5.1,a, f6.2 ,a, f6.1 ,a, f6.1 )') &
   tmp_air,',',tmp_soil1,',',tmp_soil2,',',tmp_soil3,',',cloud,',',prec,',',rh,',',wind,',',rad_short,',',rad_long
   
!*** 3rd line (Radiation properties on vegetation surface)
   write (Fn,'( f4.2,a,f4.2,a,f4.2 )') &
   albedo_mean,',',albedo_soil,',',albedo_leaf !,',',ir_tree,',',ir_grass

!*** 4th line (Physical status of radiation)
   write (Fn,'( f6.1,a, f6.1,a, f6.1,a, f6.1,a, f7.1,a, f6.1 )') &
!   par_direct,',',par_diffuse,',',par*sum(par_grass_rel(:,:))/DivedG/DivedG,',',radnet_soil,',',radnet_veg, ',', radlong_up   !!!>>>>>>>>>>>>TN:rm
   par_direct,',',par_diffuse,',',par*sum(par_grass_rel(:,:))/real(GRID%N_tot),',',radnet_soil,',',radnet_veg, ',', radlong_up !!!<<<<<<<<<<<<TN:add

!*** 5th line (Hydrogical varibales)
   write (Fn,'( f7.2, 3(a,f6.1), a,f8.1 )') &
   canopy_cond,',',sum(pool_w(1:5)),',',sum(pool_w(6:10)),',',sum(pool_w(11:20)),',',pool_snow
   
!*** 6th~ lines (PFT specific variables)
   Do p=1, PFT_no
   if (pft_exist(p)) then
      write (Fn,'( f5.2)') lai_RunningRecord(1,p)
   endif
   End do
   
!_____________ Write Monthly output
IF ( Day_of_Month(doy) == Day_in_month(Month(doy)) ) then
   
!*** 1st line
   !prepare tree density (N/ha)
   i = 0
   do no=1, Max_no
      if ( tree_exist(no) ) i = i + 1
   end do
!   tree_density = real(i) * ( (100.0/Max_loc)**2 ) !!!>>>>>>>>>>>>TN:rm
   tree_density = real(i) * ( 10000.0/real(GRID%Area) ) !!!<<<<<<<<<<<<TN:add 単位要チェック
   
   write (Fn, '(i5, a, i2, a, i2, a, f8.1, a, i2, a, f6.1)')  &
         year,',',Month(doy),',',biome,',',tree_density,',',type_grass,',',co2atm
   
!*** 2nd line (Carbon strage & flux)
   !Sumup carbon properties
   i  = Day_in_month(Month(doy))
!   x  = C_in_drymass/Max_loc/Max_loc/100.0 !(gDM/Stand) -> (MgC/ha) !!!>>>>>>>>>>>>>>>>>TN:rm
   x  = C_in_drymass/real(GRID%Area)/100.0 !(gDM/Stand) -> (MgC/ha) !!!<<<<<<<<<<<<<<<<<TN: add
   
   !Carbon fluxes (MgC/ha/month)
   a1 = sum(flux_c_uptake_RR (1:i))                           * x !Carbon uptake
   a2 = ( sum(flux_c_mnt_RR(1:i)) + sum(flux_c_gro_RR(1:i)) ) * x !Plant respiration
   a3 = sum(flux_c_htr_RR(1:i))                               * x !Heterotphic respiration
   a4 = sum(flux_c_lit_RR(1:i))                               * x !Biomass to Litter flux
   a5 = sum(flux_c_som_RR(1:i))                               * x !Litter to SOM flux
   a6 = sum(flux_c_fir_RR(1:i))                               * x !Fire emission
   
   !Litter and SOM carbon (Mg C/ha)
   b1 = (pool_litter_trunk+ pool_litter_leaf+ pool_litter_root+ pool_litter_ag+ pool_litter_bg) * x
                          !Litter  (Mg C/ha)
   b2 = pool_som_int  * x !SOM int (Mg C/ha)
   b3 = pool_som_slow * x !SOM slow(Mg C/ha)
   
   !Woody carbon (Mg C/ha)
   c1 = 0.0
   Do no=1, Max_no
   If (tree_exist(no)) then
      c1 = c1 + mass_leaf(no)+mass_trunk(no)+mass_root(no)+mass_stock(no)+mass_available(no)
   End if
   End do
!   c1 = c1 * C_in_drymass / Max_loc / Max_loc / 100.0 !!!>>>>>>>>>>>>TN:rm
   c1 = c1 * C_in_drymass / real(GRID%Area) / 100.0 !!!<<<<<<<<<<<<TN:add 単位要チェック
   
   !Grass carbon (Mg C/ha)
   c2 = sum(gmass_leaf(:,:)) + sum(gmass_root(:,:)) + sum(gmass_available(:,:)) + sum(gmass_stock(:,:))
!   c2 = c2 * C_in_drymass / 100.0 / (Max_loc**2) !!!>>>>>>>>>>>>TN:rm
   c2 = c2 * C_in_drymass / 100.0 / real(GRID%Area) !!!<<<<<<<<<<<<TN:add 単位要チェック
   
   write (Fn, '( 6(f6.1,a), 5(f8.1,a), f6.1 )') &
         a1,',',a2,',',a3,',',a4,',',a5,',',a6,',',b1,',',b2,',',b3,',',c1,',',c2
   
!*** 3rd line (Water flux)
   write (Fn, '( 8(f8.2,a), f8.2 )') &
      sum(  prec_RunningRecord  (1:Day_of_Month(doy))),',', &
      sum(flux_ro1_RunningRecord(1:Day_of_Month(doy))),',', &
      sum(flux_ro2_RunningRecord(1:Day_of_Month(doy))),',', &
      sum(flux_ic_RunningRecord (1:Day_of_Month(doy))),',', &
      sum(flux_ev_RunningRecord (1:Day_of_Month(doy))),',', &
      sum(flux_tr_RunningRecord (1:Day_of_Month(doy))),',', &
      sum(flux_sl_RunningRecord (1:Day_of_Month(doy))),',', &
      sum(flux_tw_RunningRecord (1:Day_of_Month(doy))),',', &
      sum(flux_sn_RunningRecord (1:Day_of_Month(doy)))       
   
!*** 4th line ~ (for each PFT)
   !Sumup woody biomass (Mg C / ha)
   !  Reset variables
   c_leaf     (:) = 0.0
   c_trunk    (:) = 0.0
   c_root     (:) = 0.0
   c_stock    (:) = 0.0
   c_available(:) = 0.0
   
   !  Woody PFTs
   Do no=1, Max_no
   If (tree_exist(no)) then
      p = pft(no)
      c_leaf     (p) = c_leaf     (p) + mass_leaf(no)
      c_trunk    (p) = c_trunk    (p) + mass_trunk(no)
      c_root     (p) = c_root     (p) + mass_root(no)
      c_stock    (p) = c_stock    (p) + mass_stock(no)
      c_available(p) = c_available(p) + mass_available(no)
   End if
   End do
   
!   x = C_in_drymass / Max_loc / Max_loc / 100.0 !!!>>>>>>>>>>>>TN:rm
   x = C_in_drymass / real(GRID%Area) / 100.0 !!!<<<<<<<<<<<<TN:add 単位要チェック
   c_leaf     (:) = c_leaf     (:) * x
   c_trunk    (:) = c_trunk    (:) * x
   c_root     (:) = c_root     (:) * x
   c_stock    (:) = c_stock    (:) * x
   c_available(:) = c_available(:) * x
   
   !Sunup grass biomass (Mg C / ha)
!   x = C_in_drymass / 100.0 / (Max_loc**2) !!!>>>>>>>>>>>>TN:rm
   x = C_in_drymass / 100.0 / real(GRID%Area) !!!<<<<<<<<<<<<TN:add 単位要チェック
   if (pft_exist(C3g_no)) then
      c_leaf     (C3g_no) = sum(gmass_leaf(:,:))      * x
      c_trunk    (C3g_no) = 0.0                
      c_root     (C3g_no) = sum(gmass_root(:,:))      * x
      c_stock    (C3g_no) = sum(gmass_stock(:,:))     * x
      c_available(C3g_no) = sum(gmass_available(:,:)) * x
      
      c_leaf     (C4g_no) = 0.0
      c_trunk    (C4g_no) = 0.0
      c_root     (C4g_no) = 0.0
      c_stock    (C4g_no) = 0.0
      c_available(C4g_no) = 0.0
   else
      c_leaf     (C3g_no) = 0.0
      c_trunk    (C3g_no) = 0.0
      c_root     (C3g_no) = 0.0
      c_stock    (C3g_no) = 0.0
      c_available(C3g_no) = 0.0
      
      c_leaf     (C4g_no) = sum(gmass_leaf(:,:))      * x
      c_trunk    (C4g_no) = 0.0                
      c_root     (C4g_no) = sum(gmass_root(:,:))      * x
      c_stock    (C4g_no) = sum(gmass_stock(:,:))     * x
      c_available(C4g_no) = sum(gmass_available(:,:)) * x
   endif
   
   i = Day_in_month(Month(doy)) !number of day of the current month (day)
   DO p=1, PFT_no
   if (pft_exist(p)) then
      write (Fn,'( 7(f9.4,a) )') &
!      sum(gpp_RunningRecord (1:i,p)) * C_in_drymass / Max_loc / Max_loc / 100.0,',', & !!!>>>>>>>>>>>>TN:rm
!      sum(npp_RunningRecord (1:i,p)) * C_in_drymass / Max_loc / Max_loc / 100.0,',', & !!!>>>>>>>>>>>>TN:rm
      sum(gpp_RunningRecord (1:i,p)) * C_in_drymass / real(GRID%Area) / 100.0,',', & !!!<<<<<<<<<<<<TN:add 単位要チェック
      sum(npp_RunningRecord (1:i,p)) * C_in_drymass / real(GRID%Area) / 100.0,',', & !!!<<<<<<<<<<<<TN:add 単位要チェック
      c_leaf(p),',',c_trunk(p),',',c_root(p),',',c_stock(p),',',c_available(p)
   endif
   END DO
   
END IF

!_____________ Write Annual output
IF ( Doy==Day_in_Year ) then
   
   !Prepare size class distribution for each PFT
   size_dist(:,:) = 0.0 !reset
   DO no=1, Max_no
   If (tree_exist(no)) then
      
      p = pft(no)
      i = int(   (100*dbh_heartwood(no) + 100*dbh_sapwood(no)) / Dist_step   )
      
      if (i>=41) then
            size_dist(p,41) = size_dist(p,41) + 1
      elseif (i<1) then
            size_dist(p, 0) = size_dist(p, 0) + 1
      else
            size_dist(p, i) = size_dist(p, i) + 1
      end if
      
   End if
   END DO
   
   !Write size class ditribution
   DO p=1, PFT_no
   if ( (p.ne.c3g_no) .and. (p.ne.c4g_no) ) then
   if ( pft_exist(p) ) then
      do i = 0, 10
      write (Fn,'(i4,a$)') size_dist(p,i),','
      end do
      
      do i = 11, 20
      write (Fn,'(i3,a$)') size_dist(p,i),','
      end do
      
      do i = 21, 41
      write (Fn,'(i2,a$)') size_dist(p,i),','
      end do
      
      write (Fn,*)
   endif
   endif
   END DO
   
   !Write wild fire flag
   if (dfl_fire <= 365 .and. year>=2 ) then
      write (Fn,'(i1)') 1
   else
      write (Fn,'(i1)') 0
   endif
   
END IF

!_____________ Write Annual output (Forest Structure)
IF ( Doy==Day_in_Year ) then
!Count tree number
   count = 0
   Do i=1, Max_no
      if ( tree_exist(i) ) count = count +1
   End do
   
!Write data
   write (Fn,'(i4,a1)') count,','
   Do i=1, Max_no
   if (tree_exist(i)) then
      write (Fn,'( 4(f6.1,a), 2(f7.3,a), 2(f7.4,a), i2,a )') &
      bole_x(i)                           ,',',& !1: bole location x
      bole_y(i)                           ,',',& !2: bole location y
      crown_x(i)                          ,',',& !3: crown location x
      crown_y(i)                          ,',',& !4: crown location y
      real(bole(i)) * STEP + 1.3          ,',',& !5: bole height
      real( height(i) - bole(i) ) * STEP  ,',',& !6: foliage height
      dbh_heartwood(i) + dbh_sapwood(i)   ,',',& !7: dbh
      crown_diameter(i)/2.0               ,',',& !8: crown_radius
      pft(i)                              ,','   !9: pft
      
   end if
   End do
End if

END SUBROUTINE output_for_viewer



!**************************************************************************************************
! Forest status at the every 31 December
!**************************************************************************************************
SUBROUTINE output_annual (Fn)

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
   
!Arguments
   integer,intent(IN)::Fn   !File I/O number
   
!Local variables
   real    unit_conv      !Unit converter: gDM/stand -> MgC/ha
   real    mass_wood1     !total woody  biomass (g C / m2)
   real    mass_wood2     !Aboveground woody  biomass (g C / m2)
   real    mass_grass     !grass  biomass (g C / m2)
   real    nee
   real    height_mean    !mean tree height (m)
   real    dbh_mean       !basal area (g C / m2)
   real    ba             !basal area (g C / m2)
   real    pool_litter    !litter pool
   
   real    lai_max_wood   !Annual max  of woody LAI (m2 m-2)
   real    lai_max_grass  !Annual max  of grass LAI (m2 m-2)
   real    lai_mean_wood  !Annual mean of woody LAI (m2 m-2)
   real    lai_mean_grass !Annual mean of grass LAI (m2 m-2)
   
   integer ald            !Annual maximum of Active layer Depth (Step)
   integer ald_doy        !Day of Year when Active-layer-Depth reaches max
   integer tree_counter   !tree counter
   integer no,i,j         !loop counters
   
   real x, y              !For general usage
   
!_____________ Unit converter: gDM/stand -> MgC/ha
!   unit_conv = C_in_drymass * (100/Max_loc) * (100/Max_loc) / (1000*1000) !!!>>>>>>>>>>>>TN:rm
   unit_conv = C_in_drymass / real(GRID%Area) / 100.0 !!!<<<<<<<<<<<<TN:add 単位要チェック
   
!_____________ Prepare output variables
   tree_counter =   0
   mass_wood1   = 0.0
   mass_wood2   = 0.0
   height_mean  = 0.0
   dbh_mean     = 0.0
   ba           = 0.0
   do no=1, Max_no
   if (tree_exist(no)) then
      tree_counter = tree_counter + 1
      mass_wood1   = mass_wood1   + &
                     mass_leaf(no)+mass_trunk(no)+mass_root(no)+mass_stock(no)+mass_available(no)
      mass_wood2   = mass_wood2   + &
                     mass_leaf(no)+mass_trunk(no)*0.667        +mass_stock(no)+mass_available(no)
      
      height_mean  = height_mean  + height(no)*STEP + 1.3
      dbh_mean     = dbh_mean     + dbh_heartwood(no) + dbh_sapwood(no)
      ba           = ba           + 3.14 * ( 0.5*100*(dbh_heartwood(no)+dbh_sapwood(no)) )**2
   endif
   end do
   height_mean = height_mean / max(1, tree_counter) !Mean tree height (m)
   dbh_mean    = dbh_mean    / max(1, tree_counter) !Mean DBH         (m)
!   ba          = ba /Max_loc/Max_loc                !BA (cm2 m-2) !!!>>>>>>>>>>>>TN:rm
   ba          = ba / real(GRID%Area)                !BA (cm2 m-2) !!!<<<<<<<<<<<<TN:add
   
!mass_grass: Grass biomass (gDM/stand)
   mass_grass = sum(gmass_leaf(:,:)) + sum(gmass_root(:,:)) + sum(gmass_stock(:,:)) + sum(gmass_available(:,:))
   
   lai_max_wood   = 0.0
   lai_max_grass  = 0.0
   lai_mean_wood  = 0.0
   lai_mean_grass = 0.0
   do i=1, Day_in_Year
      x = lai_RunningRecord(i,C3g_no) + lai_RunningRecord(i,C4g_no) !Grass LAI
      y = sum(lai_RunningRecord(i,:)) - x                           !Tree  LAI
      lai_max_wood   = max(lai_max_wood ,  y)
      lai_max_grass  = max(lai_max_grass,  x)
      lai_mean_wood  = lai_mean_wood  + y / real(Day_in_Year)
      lai_mean_grass = lai_mean_grass + x / real(Day_in_Year)
   enddo
   
!TMP
!ald: maximum seasonal depth of active layer
!ald_doy: DOY when ald reaches
   ald = 0
   ald_doy = 1
   do i=1, Day_in_Year
      do j=1,NumSoil ; if (tmp_soil_RunningRecord(i,j)<0.0) exit ; enddo
      
      if (j-1>ald) then
         ald     = j-1
         ald_doy = 366-i
      endif
      
      if (ald==NumSoil) exit
   enddo
   
!Others
   pool_litter = pool_litter_trunk + pool_litter_leaf + pool_litter_root + pool_litter_ag + pool_litter_bg
   nee         = sum(flux_c_uptake_RR(:)) &
               - sum(flux_c_mnt_RR(:)) - sum(flux_c_gro_RR(:)) - sum(flux_c_htr_RR(:)) - sum(flux_c_fir_RR(:))
   
!_____________ Write title and data ***
   if (year==1) then
   write (Fn, '(80a)', advance='no') &
   'Year, TreeDensity, Heigh, DBH, BA, MassW, MassW_ag, MassG, Litter, SOM, GPP, NPP'
   write (Fn, '(79a)') &
   ' , NEE, LAImax_W, LAImax_G, LAImean_W, LAImean_G, ccon, RO, Biome, ALD, ALD_DOY'
   end if
   
   write (Fn,'( i4,a, 2(f7.1,a), 2(f9.3,a), 5(f6.1,a), 2(f5.1,a), f6.1,a, 4(f5.1,a), f7.3,a, f6.1,a, 2(i3,a), i3 )') &
   year                                      , ',', & ! 1 Simulation year
!   tree_counter * ( (100.0/Max_loc)**2 )     , ',', & ! 2 Tree density (ha-1) !!!>>>>>>>>>>>>TN:rm
   tree_counter * ( 10000.0/real(GRID%Area) )     , ',', & ! 2 Tree density (ha-1) !!!<<<<<<<<<<<<TN:add 単位要チェック
   height_mean                               , ',', & ! 3 Mean tree height (m)
   dbh_mean                                  , ',', & ! 4 Mean DBH   (m)
   ba                                        , ',', & ! 5 Basal area (cm2 m-2) (m2 ha-1)
   mass_wood1                     * unit_conv, ',', & ! 6 Total woody biomass (MgC ha-1)
   mass_wood2                     * unit_conv, ',', & ! 7 Above ground woody biomass (MgC ha-1)
   mass_grass                     * unit_conv, ',', & ! 8 Total grass biomass (MgC ha-1)
   pool_litter                    * unit_conv, ',', & ! 9 Litter carbon       (MgC ha-1)
   (pool_som_int + pool_som_slow) * unit_conv, ',', & !10 SOM carbon          (MgC ha-1)
   sum(gpp_RunningRecord (:,:))   * unit_conv, ',', & !11 GPP                 (MgC ha-1 yr-1)
   sum(npp_RunningRecord (:,:))   * unit_conv, ',', & !12 NPP                 (MgC ha-1 yr-1)
   nee                            * unit_conv, ',', & !13 NEE                 (MgC ha-1 yr-1)
   lai_max_wood                              , ',', & !14 Annual max  of woody LAI (m2 m-2)
   lai_max_grass                             , ',', & !15 Annual max  of grass LAI (m2 m-2)
   lai_mean_wood                             , ',', & !16 Annual mean of woody LAI (m2 m-2)
   lai_mean_grass                            , ',', & !17 Annual mean of grass LAI (m2 m-2)
   canopy_cond                               , ',', & !18 Canopy conductance
     sum(flux_ro1_RunningRecord(:))                 &
   + sum(flux_ro2_RunningRecord(:))          , ',', & !19 Runoff (mm year-1)
   biome                                     , ',', & !20 Biome type
   ald                                       , ',', & !21 Annual Max of Active Layer Depth (Step)
   ald_doy                                            !22 DOY of the ALD reached
   
END SUBROUTINE output_annual



!**************************************************************************************************
! Output forest 3D structure
!**************************************************************************************************
SUBROUTINE output_forest (Fn)

!_____________ Set variables
!Namespace
   USE data_structure
   USE time_counter
   USE vegi_status_current1
   implicit none
   
!Arguments
   integer,intent(IN)::Fn   !File I/O number
   
!Local variable
   integer count, i
   
!_____________ Main part
!Count tree number
   count = 0
   Do i=1, Max_no
      if ( tree_exist(i) ) count = count +1
   End do
   
!Write data
!   write (Fn,'(i4,a1)') count,',' !!!>>>>>>>>>>>>TN:rm
   write (Fn,'(i9,a1)') count,',' !!!<<<<<<<<<<<<TN:add
   write (*,*) 'Number of tree:', count !!!<<<<<<<<<<<<TN:add
   Do i=1, Max_no
   if (tree_exist(i)) then
!      write (Fn,'( 4(f6.1,a), 2(f7.3,a), 2(f7.4,a), i2,a )') & !!!>>>>>>>>>>>>TN:rm
      write (Fn,'( 4(f8.2,a), 2(f7.3,a), 2(f7.4,a), i2,a )') & !!!<<<<<<<<<<<<TN:add
      bole_x(i)                           ,',',& !1: bole location x
      bole_y(i)                           ,',',& !2: bole location y
      crown_x(i)                          ,',',& !3: crown location x
      crown_y(i)                          ,',',& !4: crown location y
      real(bole(i)) * STEP + 1.3          ,',',& !5: bole height
      real( height(i) - bole(i) ) * STEP  ,',',& !6: foliage height
      dbh_heartwood(i) + dbh_sapwood(i)   ,',',& !7: dbh
      crown_diameter(i)/2.0               ,',',& !8: crown_radius
      pft(i)                              ,','   !9: pft
      
   end if
   End do
   
   !For separate files
   !character :: file_name*32
   !write(file_name,*) year
   !file_name = trim(file_name) // '.txt'
   !open (Fn, file=file_name)
   !close (Fn)
   
END SUBROUTINE output_forest



!**************************************************************************************************
! Output air variables
!**************************************************************************************************
SUBROUTINE output_air (Fn)

!_____________ Set variables
!Namespace
   USE data_structure
   USE time_counter
   USE grid_status_current2
   implicit none
   
!Augments
   integer,intent(IN)::Fn   !File I/O number
   
!_____________ Main part
   !Write title
   if (year == 1 .and. doy==1) then
   write (Fn,*) &
   '  Yr  doy  ap    vps      vp      vpd    dnsa  slope_vps'
   end if
   
   !Write data
   write (Fn,'(2i4, 1f8.1, 5f8.3, f8.2)') &
   year, doy, ap, vp_sat, vp, vpd, dnsa, slope_vps
   
END SUBROUTINE output_air



!**************************************************************************************************
! Output climate variables
!**************************************************************************************************
SUBROUTINE output_climate (Fn, cloud, prec, rh, wind, tmp_air, tmp_soil, rad_short, rad_long)

!_____________ Set variables
!Namespace
   USE data_structure
   USE time_counter
   USE grid_status_current2
   implicit none
   
!Arguments
   integer,intent(IN)::Fn     !File I/O number
   
   real,intent(IN)::tmp_air   ! 2m air temperature (Celcius)
   real,intent(IN)::tmp_soil  ! average soil temperature for 50cm depth (Celcius)
   real,intent(IN)::cloud     ! total cloudness (fraction)
   real,intent(IN)::prec      ! precipitation (mm day-1)
   real,intent(IN)::rh        ! relative humidity (%)
   real,intent(IN)::wind      ! wind velocity (m s-1)
   real,intent(IN)::rad_short ! 
   real,intent(IN)::rad_long  ! 
   
!_____________ Main part
!Write title
   if (year == 1 .and. doy==1) then
   write (Fn,*) &
   '  yr  doy    air   soil  cloud    prec   rh  vp/vpsat   wind  rad_short  rad_long'
   end if
   
!Write data
   write (Fn,'(2i5, 2f7.1, f5.2, f8.2, f6.1, f8.3, 3f7.1 )') &
   year, doy, tmp_air, tmp_soil, cloud, prec, rh, vp/vp_sat, wind, rad_short, rad_long
   
END SUBROUTINE output_climate



!**************************************************************************************************
! Output air variables
!**************************************************************************************************
SUBROUTINE output_radiation (Fn, rad_short)

!_____________ Set variables
!Namespace
   USE data_structure
   USE time_counter
   USE grid_status_current2
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none
   
!Augments
   integer,intent(IN)::Fn        !File I/O number
   real   ,intent(IN)::rad_short 
   
!_____________ Main part
!Write title
   if (year == 1 .and. doy==1) then
   write (Fn,*) &
   ' Yr doy  s_hgt dlen  rad_short    par     par_grass'
   end if
   
!Write data
   write (Fn,'(2i4, 2f6.1, 1f7.1, 2f8.1)') &
!   year, doy, sl_hgt(doy), dlen(doy), rad_short, par, par*sum(par_grass_rel(:,:))/DivedG/DivedG !!!>>>>>>>>>>>>TN:rm
   year, doy, sl_hgt(doy), dlen(doy), rad_short, par, par*sum(par_grass_rel(:,:))/real(GRID%N_tot) !!!<<<<<<<<<<<<TN:add
   
END SUBROUTINE output_radiation



!**************************************************************************************************
! Output net radiation
!**************************************************************************************************
SUBROUTINE output_netradiation (Fn)

!_____________ Set variables
!Namespace
   USE data_structure
   USE time_counter
   USE grid_status_current2
   implicit none
   
!Augments
   integer,intent(IN)::Fn   !File I/O number
   
!_____________ Main part
!Write title
   1 format(E10.4e1,a)  !Output format (Sample: " 0.2846E+3"+",")
   2 format(E10.4e1  )  !Output format (Sample: " 0.2846E+3"    )
   
   if (year == 1 .and. doy==1) then
      write (Fn,'(a100)') &
     'Yr, doy, albedo_soil, albedo_leaf, albedo_mean, radnet_soil, radnet_veg, radlong_up'
   end if
   
!Write data
   write (Fn,'(i4,a)'   , advance='no') year       , ','
   write (Fn,'(i3,a)'   , advance='no') doy        , ','
   write (Fn,'(f5.2,a)' , advance='no') albedo_soil, ','
   write (Fn,'(f5.2,a)' , advance='no') albedo_leaf, ','
   write (Fn,'(f5.2,a)' , advance='no') albedo_mean, ','
   write (Fn, 1         , advance='no') radnet_soil, ','
   write (Fn, 1         , advance='no') radnet_veg , ','
   write (Fn, 1         , advance='no') radlong_up
   write (Fn,*)
   
END SUBROUTINE output_netradiation



!**************************************************************************************************
! Output water status
!**************************************************************************************************
SUBROUTINE output_water (Fn, W_fi, W_wilt, W_sat)

!_____________ Set variables
   USE data_structure
   USE time_counter
   USE grid_status_current1
   implicit none
   
!Arguments
   integer,intent(IN)::Fn   !File I/O number
   
   real,intent(IN)::W_fi         !filed capacity   (m3/m3, 0.0 -> 1.0)
   real,intent(IN)::W_wilt       !wilting point    (m3/m3, 0.0 -> 1.0)
   real,intent(IN)::W_sat        !saturate point   (m3/m3, 0.0 -> 1.0)
   
!Local valuable
   real a1, a2
   
!_____________ Main part
!Prepare values
   !a1 = pool_w1 / (W_sat * Depth)
   !a2 = pool_w2 / (W_sat * Depth)
   
   a1 = sum(pool_w(1:5)) / 5.0
   a1 = (a1/Depth-W_wilt) / (W_fi-W_wilt)
   a1 = max(0.0, min(1.0, a1))
   
   a2 = sum(pool_w(6:15)) / 10.0
   a2 = (a2/Depth-W_wilt) / (W_fi-W_wilt)
   a2 = max(0.0, min(1.0, a2))
   
!Write title
   if (year == 1 .and. doy==1) then
   write (Fn,*) &
   'Yr doy   snow   pool_w1   pool_w2  stat_w1   stat_w2'
   end if
   
!Write data
   write (Fn,'(2i4, 5f10.2, 2f10.3)') &
   year, doy, pool_snow, sum(pool_w(1:5))/5.0, sum(pool_w(6:15))/10.0, a1, a2
   
END SUBROUTINE output_water



!**************************************************************************************************
! Output carbon flux status
!**************************************************************************************************
SUBROUTINE output_cflux (Fn)

!_____________ Set variables
!Namespace
   USE data_structure
   USE time_counter
   USE grid_status_current1
   USE vegi_status_current1
   USE vegi_status_current2
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none
   
!Arguments
   integer,intent(IN)::Fn   !File I/O number
   
!Local valiables
   real unit_conv !Unit converters
   
!_____________ Unit converter2
!   unit_conv = C_in_drymass * ( (100/real(Max_loc))**2 ) / (1000*1000) !gDM/stand -> MgC/ha!!!>>>>>>>>>>>>TN:rm
   unit_conv = C_in_drymass / real(GRID%Area) / 100.0 !gDM/stand -> MgC/ha !!!<<<<<<<<<<<<TN:add 単位要チェック
   
!_____________ Main part
!Write title
!   if (year==1 .and. doy==1) then
!      write (Fn,'(154a)') 'Yr, DOY, flux_c_uptake, flux_c_mnt, flux_c_gro, flux_c_htr, flux_c_fir, GPP, NPP, r_trunk, r_leaf, r_root, r_ag, r_bg, l_trunk, l_leaf, l_root, l_ag, l_bg'
!   end if
   
!Write data, all units except time counter are (Mg C/day/ha)
   write (Fn,'( i4, a, i3, 17(a, f8.5) )')      & 
      year                                     ,',', & !Year
      doy                                      ,',', & !doy
      flux_c_uptake_RR(1)          * unit_conv ,',', & !carbon uptake  
      flux_c_mnt_RR(1)             * unit_conv ,',', & !carbon emission due to maintenance  respiration
      flux_c_gro_RR(1)             * unit_conv ,',', & !carbon emission due to growth       respiration
      flux_c_htr_RR(1)             * unit_conv ,',', & !carbon emission due to heterotropic respiration
      flux_c_fir_RR(1)             * unit_conv ,',', & !carbon emission due to fire incident
      sum(gpp_RunningRecord (1,:)) * unit_conv ,',', & !GPP
      sum(npp_RunningRecord (1,:)) * unit_conv ,',', & !NPP
      sum(resp_trunk(:))           * unit_conv ,',', & !respiration of woody trunk      
      sum(resp_leaf (:))           * unit_conv ,',', & !respiration of woody leaf       
      sum(resp_root (:))           * unit_conv ,',', & !respiration of woody fine root  
      resp_grass_ag                * unit_conv ,',', & !respiration of grass aboveground
      resp_grass_bg                * unit_conv ,',', & !respiration of grass underground
      flux_litter_trunk            * unit_conv ,',', & !litter flux of woody trunk      
      flux_litter_leaf             * unit_conv ,',', & !litter flux of woody leaf       
      flux_litter_root             * unit_conv ,',', & !litter flux of woody fine root  
      flux_litter_ag               * unit_conv ,',', & !litter flux of grass aboveground
      flux_litter_bg               * unit_conv         !litter flux of grass underground
      
END SUBROUTINE output_cflux



!**************************************************************************************************
! Output water flux status
!**************************************************************************************************
SUBROUTINE output_wflux (Fn, prec)

!_____________ Set variables
!Namespace
   USE data_structure
   USE time_counter
   USE grid_status_current1
   USE grid_status_current2
   implicit none
   
!Arguments
   integer,intent(IN)::Fn   !File I/O number
   real   ,intent(IN)::prec !precipitation (mm/day)
   
!Local variables
   real aet, lh
   
!_____________ Main part
!Preparation
   !aet: actual evapotranspiration (mm/m2/day or kg/m2/day)
   aet = flux_tr_RunningRecord(1) + flux_ev_RunningRecord(1) + flux_ic_RunningRecord(1) + flux_sl_RunningRecord(1)
   
   !lh: latent heat of water vaporization (MJ/kg H2O)
   !    2.259  (= latent heat at 100 degree Celsius in MJ/kg H2O)
   !    0.0042 (= Energy requirement for warming H2O by 1deg-Celcius in MJ/Kg H20)
   lh = 2.259 + 0.0042 * (100-tmp_air_RunningRecord(1))
   
!Write title
   if (year == 1 .and. doy==1) then
      write (Fn,*) 'Yr, doy, Precip, RunOff, Intcep, Evapor, Transp, SensiHeat, LatenHeat'
      write (Fn,*) ' -,   -, mm/day, mm/day, mm/day, mm/day, mm/day, MJ/m2/day, MJ/m2/day'
   end if
   
!Write data
   1 format(f6.1,a)  !Output format
   2 format(f6.1  )  !Output format
   
   write (Fn,'(i5,a)', advance="no") year,','
   write (Fn,'(i3,a)', advance="no") doy ,','
   write (Fn, 1, advance="no") prec                    ,','
   write (Fn, 1, advance="no") flux_ro1_RunningRecord(1)+flux_ro2_RunningRecord(1), ','
   write (Fn, 1, advance="no") flux_ic_RunningRecord(1), ','
   write (Fn, 1, advance="no") flux_ev_RunningRecord(1), ','
   write (Fn, 1, advance="no") flux_tr_RunningRecord(1), ','
   write (Fn, 1, advance="no") (radnet_veg*dlen(doy) + radnet_soil*24)*60*60/1000000 - lh*aet, ','
   write (Fn, 2, advance="no") lh*aet
   write (Fn,               *) 
   
END SUBROUTINE output_wflux


!**************************************************************************************************
! Output vertical leaf distribution
!出力するPFTに応じて、コードを変更する必要があります
!**************************************************************************************************
SUBROUTINE output_ld_vertical (Fn)

!_____________ Set variables
!Namespace
   USE data_structure
   USE time_counter
   USE vegi_status_current1
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none
   
!Arguments
   integer,intent(IN)::Fn   !File I/O number
   
!Local parameters
   integer,parameter::Vertical_Interval = 10 !number of STEP for each vertical layer for output
   
!Local variables
   integer        no, p, l1, l2
   real           x 
   real,dimension(PFT_no, Max_hgt )                       ::ld1 !sum of leaf area for each step
   real,dimension(PFT_no, int(Max_hgt/Vertical_Interval) )::ld2 !sum of leaf area for each interval
   
!_________________ Sumup leaf area for each tree
   ld1(:,:) = 0.0
   Do no=1, Max_no
   if ( .not. tree_exist(no) ) cycle
      p = pft(no)
      x = la(no) / real(height(no) - bole(no)) !leaf area for each crown disk
      do l1 = bole(no)+1, height(no)
         ld1(p,l1) = ld1(p,l1) + x
      end do
   End do
!   ld1(:,:) = ld1(:,:) / Max_loc / Max_loc !adjust unit: (m2 step-1 stand-1) --> (m2 step-1 m-2) !!!>>>>>>>>>>>>TN:rm
   ld1(:,:) = ld1(:,:) / real(GRID%Area) !adjust unit: (m2 step-1 stand-1) --> (m2 step-1 m-2) !!!<<<<<<<<<<<<TN:add
   
!_________________ Convert vertical resolution
   ld2(:,:) = 0.0
   Do p=1, PFT_no
   Do l1=1, Max_hgt
      l2 = int( (l1-1)/Vertical_Interval ) + 1
      ld2(p,l2) = ld2(p,l2) + ld1(p,l1)
   End do
   End do
   
!_________________ Write Time stamp
   write (Fn, '( 1(i5,a),1(i5,a),1(i7,a) )') year,',',doy,',',counter
   
!_________________ Write distribution
   write (Fn, *)
   Do l2=1, int(Max_hgt/Vertical_Interval)
      !write (Fn,'( f8.5,a,f8.5,a,f8.5,a,f8.5 )') &
      !      ( ld2(1,l2),',',ld2(2,l2),',',ld2(3,l2),',',ld2(4,l2) )
      
      write (Fn,'(f8.5)') ld2(1,l2)!!!<<<<<<<<<<<<TN:changed 10->1
      
   End do
   
END SUBROUTINE output_ld_vertical



!**************************************************************************************************
! Output grass
!**************************************************************************************************
SUBROUTINE output_grass (Fn)

!_____________ Set variables
!Namaspace
   USE data_structure
   USE time_counter
   USE vegi_status_current1
   USE vegi_status_current2
   USE grid_status_current1
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none
   
!Arguments
   integer,intent(IN)::Fn   !File I/O number
   
!Local variable
   integer   p, type_grass
   character phenology_symbol
   
!_________________ set phenology_symbol
   if (pft_exist(C3g_no)) then
      p = C3g_no
      type_grass = 3
   else
      p = C4g_no
      type_grass = 4
   endif
   
   phenology_symbol='-' ; if (phenology(p)) phenology_symbol='*'
   
!_________________ Write title
   if (year==1 .and. doy==1) then
   write (Fn,*) ' Yr doy pt       F     R     S    Av    Sd   lai lai_o'
   end if
   
!_________________ Write data
   write (Fn,'( i3 , i4 , i3 , a3, 5f6.0, 2f5.1 )') &
   year, doy, type_grass, phenology_symbol, &
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:rm
!   sum(gmass_leaf      (:,:))             /Max_loc/Max_loc, & !Biomass Aboveground
!   sum(gmass_root      (:,:))             /Max_loc/Max_loc, & !Biomass Root       
!   sum(gmass_stock     (:,:))             /Max_loc/Max_loc, & !Biomass Stock      
!   sum(gmass_available (:,:))             /Max_loc/Max_loc, & !Biomass Available  
!   pool_fuel_standG                       /Max_loc/Max_loc, & !Standing dead
!   sum(lai_grass    (:,:))                /DivedG /DivedG , & !Grass LAI
!   sum(lai_opt_grass_RunningRecord(1,:,:))/DivedG /DivedG     !Optimum grass LAI
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:rm
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:add
   sum(gmass_leaf      (:,:))             /real(GRID%Area ), & !Biomass Aboveground
   sum(gmass_root      (:,:))             /real(GRID%Area ), & !Biomass Root       
   sum(gmass_stock     (:,:))             /real(GRID%Area ), & !Biomass Stock      
   sum(gmass_available (:,:))             /real(GRID%Area ), & !Biomass Available  
   pool_fuel_standG                       /real(GRID%Area ), & !Standing dead
   sum(lai_grass    (:,:))                /real(GRID%N_tot), & !Grass LAI
   sum(lai_opt_grass_RunningRecord(1,:,:))/real(GRID%N_tot)    !Optimum grass LAI
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:add
   
END SUBROUTINE output_grass



!**************************************************************************************************
! Biomass composition
!**************************************************************************************************
SUBROUTINE output_biomass (Fn)

!_____________ Set variables
!Namespace
   USE data_structure
   USE vegi_status_current1
   USE vegi_status_current2
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none
   
!Arguments
   integer,intent(IN)::Fn   !File I/O number
   
!Local variables
   real,dimension(PFT_no)::biomass !biomass for each PFT (g C / m2)
   integer                 no      !for loop counter
   real                    x       !for gengeral usage
   
   integer                 i, j, k !for string procudure
   character(len=3)        string  !for string procudure
   
!_____________ Main part
!reset output variables
   biomass(:)=0.0
   
!woody biomass (g dm / stand)
   Do no=1, Max_no
   If (tree_exist(no)) then
      biomass(pft(no)) = biomass(pft(no))+&
      mass_leaf(no)+mass_trunk(no)+mass_root(no)+mass_stock(no)+mass_available(no)
   End if
   End do
   
!grass biomass (g dm / stand)
   x = sum(gmass_leaf(:,:)) + sum(gmass_root(:,:)) + sum(gmass_stock(:,:)) + sum(gmass_available(:,:))
   if (pft_exist(C3g_no)) then
      biomass(C3g_no) = x
      biomass(C4g_no) = 0.0
   else
      biomass(C3g_no) = 0.0
      biomass(C4g_no) = x
   endif
   
!adjust unit ( g dm / stand ---> Mg C / ha )
!   biomass(:) = biomass(:) * C_in_drymass / Max_loc / Max_loc / 100.0 !!!>>>>>>>>>>>>TN:rm
   biomass(:) = biomass(:) * C_in_drymass / real(GRID%Area) / 100.0 !!!<<<<<<<<<<<<TN:add 単位要チェック
   
!Preare string for output variables
   i = int(  PFT_no       /100 ) ; string(1:1)= char(i+48)
   j = int( (PFT_no-i*100)/ 10 ) ; string(2:2)= char(j+48)
   k =    (  PFT_no-i*100-j*10 ) ; string(3:3)= char(k+48)
   
!write title and data
   write (Fn,'( '//string//'(f6.2,1x) )') biomass(:)
   
END SUBROUTINE output_biomass



!**************************************************************************************************
! Biomass composition
!**************************************************************************************************
SUBROUTINE output_lai(Fn)

!_____________ Set variables
!Namespace
   USE data_structure
   USE vegi_status_current1
   USE vegi_status_current2
   implicit none
   
!Arguments
   integer,intent(IN)::Fn   !File I/O number
   
!Local variables
   integer                 i, j, k  !for string procudure
   character(len=3)        string   !for string procudure
   
!_____________ Main part
!Preare string for output variables
   i = int(  PFT_no       /100 ) ; string(1:1)= char(i+48)
   j = int( (PFT_no-i*100)/ 10 ) ; string(2:2)= char(j+48)
   k =    (  PFT_no-i*100-j*10 ) ; string(3:3)= char(k+48)
   
!Write title and data
   write (Fn,'( '//string//'(f4.1,1x) )') lai_RunningRecord(1,:)
   
END SUBROUTINE output_lai



!*************************************************************************************************
! output file maker (for global grid) (called @ the end of each month)
!*************************************************************************************************
SUBROUTINE output_global (Fn, W_fi, W_wilt, tmp_soil_Today)

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
   
!Arguments
   integer,intent(IN)::Fn                      !File I/O number
   real   ,intent(IN)::W_fi                    !Filed capacity   (m3/m3, 0.0 -> 1.0)
   real   ,intent(IN)::W_wilt                  !Wilting point    (m3/m3, 0.0 -> 1.0)
   real   ,intent(IN)::tmp_soil_Today(NumSoil) !Soil temperature for each layers (Celcius)
   
!Local variables
   !for output variables
   real    woody_carbon, grass_carbon, som1, som2, som3   !(kg C / m2)
   real    gpp_out, npp_out, nep_out, htr_out, npp_tree, npp_C3g, npp_C4g !(kg C / m2 / month)
   real    mean_woody_mass      !(g C / tree)
   real    tree_density         !(N / ha)
   real    lai_wood             !(m2/m2)
   real    ald                  !Active layer depth (m)
   real    water_available      !Available water on the top 5 soil layers [mm]
   
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
!   tree_density = 10000.0 * real(i) / Max_loc / Max_loc !!!>>>>>>>>>>>>>>>TN: rm
   tree_density = 10000.0 * real(i) / real(GRID%Area) !!!<<<<<<<<<<<<<<<TN: add
   
   !mean_woody_mass (kg C / tree)
   mean_woody_mass = x * C_in_drymass / max(1.0, real(i)) / 1000.0
   
   !woody_carbon (kg C / m2)
!   woody_carbon = x * C_in_drymass / Max_loc / Max_loc / 1000.0 !!!>>>>>>>>>>>>>>>TN: rm
   woody_carbon = x * C_in_drymass / real(GRID%Area) / 1000.0 !!!<<<<<<<<<<<<<<<TN: add
   
   !grass_carbon (kg C / m2)
   grass_carbon = sum(gmass_leaf(:,:)) + sum(gmass_root(:,:)) + sum(gmass_stock(:,:)) + sum(gmass_available(:,:))
!   grass_carbon = grass_carbon * C_in_drymass / 1000.0 / Max_loc / Max_loc !!!>>>>>>>>>>>>>>>TN: rm
   grass_carbon = grass_carbon * C_in_drymass / 1000.0 / real(GRID%Area) !!!<<<<<<<<<<<<<<<TN: add
   
   !soil_carbon (kg C / m2)
!   som1 = (pool_litter_trunk + pool_litter_leaf + pool_litter_root + pool_litter_ag + pool_litter_bg) & !!!>>>>>>>>>>>>>>>TN: rm
!                        * C_in_drymass / Max_loc / Max_loc / 1000.0 !!!>>>>>>>>>>>>>>>TN: rm
!   som2 = pool_som_int  * C_in_drymass / Max_loc / Max_loc / 1000.0 !!!>>>>>>>>>>>>>>>TN: rm
!   som3 = pool_som_slow * C_in_drymass / Max_loc / Max_loc / 1000.0 !!!>>>>>>>>>>>>>>>TN: rm
   som1 = (pool_litter_trunk + pool_litter_leaf + pool_litter_root + pool_litter_ag + pool_litter_bg) & !!!<<<<<<<<<<<<<<<TN: add
                        * C_in_drymass / real(GRID%Area) / 1000.0 !!!<<<<<<<<<<<<<<<TN: add
   som2 = pool_som_int  * C_in_drymass / real(GRID%Area) / 1000.0 !!!<<<<<<<<<<<<<<<TN: add
   som3 = pool_som_slow * C_in_drymass / real(GRID%Area) / 1000.0 !!!<<<<<<<<<<<<<<<TN: add
   
   !GPP & NPP (kg C / m2 / month)
   i = Day_in_month(Month(doy)) !number of day of the current month (day)
!   gpp_out   = Sum(gpp_RunningRecord(1:i, :)) * C_in_drymass / Max_loc / Max_loc / 1000.0 !!!>>>>>>>>>>>>>>>TN: rm
!   npp_out   = Sum(npp_RunningRecord(1:i, :)) * C_in_drymass / Max_loc / Max_loc / 1000.0 !!!>>>>>>>>>>>>>>>TN: rm
   gpp_out   = Sum(gpp_RunningRecord(1:i, :)) * C_in_drymass / real(GRID%Area) / 1000.0 !!!<<<<<<<<<<<<<<<TN: add
   npp_out   = Sum(npp_RunningRecord(1:i, :)) * C_in_drymass / real(GRID%Area) / 1000.0 !!!<<<<<<<<<<<<<<<TN: add
   
   npp_tree  = 0.0
   do p=1, PFT_no
      if (p==C3g_no .or. p==C4g_no) cycle
      npp_tree = npp_tree + Sum( npp_RunningRecord(1:i, p) )
   enddo
!   npp_tree = npp_tree                             * C_in_drymass / Max_loc / Max_loc / 1000.0 !!!>>>>>>>>>>>>>>>TN: rm
!   
!   npp_C3g   = Sum(npp_RunningRecord(1:i, C3g_no)) * C_in_drymass / Max_loc / Max_loc / 1000.0 !!!>>>>>>>>>>>>>>>TN: rm
!   npp_C4g   = Sum(npp_RunningRecord(1:i, C4g_no)) * C_in_drymass / Max_loc / Max_loc / 1000.0 !!!>>>>>>>>>>>>>>>TN: rm
   npp_tree = npp_tree                             * C_in_drymass / real(GRID%Area) / 1000.0 !!!<<<<<<<<<<<<<<<TN: add
   
   npp_C3g   = Sum(npp_RunningRecord(1:i, C3g_no)) * C_in_drymass / real(GRID%Area) / 1000.0 !!!<<<<<<<<<<<<<<<TN: add
   npp_C4g   = Sum(npp_RunningRecord(1:i, C4g_no)) * C_in_drymass / real(GRID%Area) / 1000.0 !!!<<<<<<<<<<<<<<<TN: add
   
   !NEP (kg C / m2 / month)
   i  = Day_in_month(Month(doy))
   nep_out = sum(flux_c_uptake_RR(1:i)) &
           - sum(flux_c_mnt_RR(1:i)) - sum(flux_c_gro_RR(1:i)) &
           - sum(flux_c_htr_RR(1:i)) - sum(flux_c_fir_RR(1:i))  
!   nep_out = nep_out * C_in_drymass / Max_loc / Max_loc / 1000.0 !!!>>>>>>>>>>>>>>>TN: rm
   nep_out = nep_out * C_in_drymass / real(GRID%Area) / 1000.0 !!!<<<<<<<<<<<<<<<TN: add
   
   !Hetertrophic respiration (kg C / m2 / month)
   i  = Day_in_month(Month(doy))
   htr_out = sum( flux_c_htr_RR(1:i) )
!   htr_out = htr_out * C_in_drymass / Max_loc / Max_loc / 1000.0 !!!>>>>>>>>>>>>>>>TN: rm
   htr_out = htr_out * C_in_drymass / real(GRID%Area) / 1000.0 !!!<<<<<<<<<<<<<<<TN: add
   
   !Tree LAI (m2/m2)
!   lai_wood = sum(la) / Max_loc / Max_loc !!!>>>>>>>>>>>>>>>TN: rm
   lai_wood = sum(la) / real(GRID%Area) !!!<<<<<<<<<<<<<<<TN: add
   
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
   
   !Active layer depth [m]
   do no=1, NumSoil
      if ( tmp_soil_Today(no)<0.0 ) exit
   enddo
   ald = Depth * (no-1) / 1000.
   
   !Available water on the top 5 soil layers [mm]
   water_available = 0.0
   do no=1, 5
      if ( tmp_soil_Today(no)<0.0 ) exit
      water_available = water_available + (pool_w(no) - W_wilt*Depth)
   enddo
   
!_____________ Write output-data
write (Fn,'( 2(1x,i3), 7(1x,f7.2), 1(1x,f9.1), 3(1x,f8.3), 1(1x,f7.1), &
2(1x,f4.1), 4(1x,f6.1), 1x, f6.0, 1x, f6.2, 1x,f5.3, 4(1x,f8.3), 2(1x,f5.3), 1x,f5.1 )') &
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
!   sum(lai_grass(:,:)) /DivedG /DivedG  , & !14: LAI of grass PFTs [m2/m2]!!!>>>>>>>>>>>>>>>TN: rm
   sum(lai_grass(:,:)) /real(GRID%N_tot)  , & !14: LAI of grass PFTs [m2/m2] !!!<<<<<<<<<<<<<<<TN: add
   
   sum(flux_ro1_RunningRecord(1:Day_of_Month(doy)))+ &
   sum(flux_ro2_RunningRecord(1:Day_of_Month(doy))), & !15: runoff        [mm/month]
   sum(flux_ic_RunningRecord (1:Day_of_Month(doy))), & !16: interception  [mm/month]
   sum(flux_ev_RunningRecord (1:Day_of_Month(doy))), & !17: evaporation   [mm/month]
   sum(flux_tr_RunningRecord (1:Day_of_Month(doy))), & !18: transpiration [mm/month]
   
   tree_density                                   , & !19: tree_density         [N / ha]
   canopy_cond                                    , & !20: stomatal conductance [mol H2O m-2 s-1]
   real(fire_number)/(real(counter)/real(Day_in_Year)), & !21: Fire frequency   [n/year]
   
   npp_tree                                       , & !22: Woody    NPP   [kg C / m2 / month]
   npp_C3g                                        , & !23: C3 grass NPP   [kg C / m2 / month]
   npp_C4g                                        , & !24: C4 grass NPP   [kg C / m2 / month]
   htr_out                                        , & !25: HetroT. Resp.  [kg C / m2 / month]
   frac_crown_coverage                            , & !26: Crown coverage [fraction]
   ald                                            , & !27: Active Layer Depth [m]
   water_available                                    !28: Available water on the top 5 soil layers [mm]
   
END SUBROUTINE output_global

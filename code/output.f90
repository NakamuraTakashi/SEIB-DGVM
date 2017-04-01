!*************************************************************************************************
! Output for Viewer
!*************************************************************************************************
SUBROUTINE output_for_viewer &
   (Fn, LAT, LON, ALT, Albedo_soil0, W_fi, W_wilt, W_sat, W_mat, &
   cloud, prec, humid, wind, tmp_air, tmp_soil)

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
   real,intent(IN)::W_mat        !matrix potential
   
   !Climatic data (single precision)
   real,intent(IN)::tmp_air  !2m air temperature              (Celcius) 
   real,intent(IN)::tmp_soil !average soil temperature for 50cm depth (Celcius) 
   real,intent(IN)::cloud    !total cloudiness                (fraction)
   real,intent(IN)::prec     !precipitation                   (mm day-1)
   real,intent(IN)::humid    !Specific humidity               (kg kg-1)
   real,intent(IN)::wind     !wind velocity                   (m s-1)   
   
!Local variables
   !character(len=3)                  PFT_no_string
   real   ,dimension(PFT_no)       ::tmp0, tmp1, tmp2, tmp3
   real   ,dimension(PFT_no)       ::c_leaf, c_trunk, c_root, c_stock, c_available
   integer,dimension(PFT_no,0:41)::size_dist
   
   integer i, p, no, type_grass
   real    a1, a2, a3, b1, b2, b3, c1, c2, x
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
   
!_____________ Write File header
If ( counter==1 .or.  (Flag_spinup_read .and. counter==1+Spinup_year*Day_in_Year) ) then
   
   !write out setting
   a1 =  500.0 !A variable for providing consistency with older versions of the SEIB-Viewer
   a2 = 1000.0 !A variable for providing consistency with older versions of the SEIB-Viewer
   a3 = 1500.0 !A variable for providing consistency with older versions of the SEIB-Viewer
   write (Fn, '( 3(f7.2,a), i5, a, i5, 3(a,f8.1), a,i4 )' )  &
!   LAT,',',LON,',',ALT,',',Simulation_year,',',PFT_no,',',a1,',',a2,',',a3,',',Max_loc !!!>>>>>>>>>>>>TN:rm
   LAT,',',LON,',',ALT,',',Simulation_year,',',PFT_no,',',a1,',',a2,',',a3,',',GRID%N_x !!!<<<<<<<<<<<<TN:add y方向の情報がない。要検討
   
   write (Fn, '( f7.5,a,f7.5,a,f7.5,a,f7.5,a,f7.5 )' ) &
   Albedo_soil0,',',W_fi,',',W_wilt,',',W_sat,',',W_mat
   
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
   write (Fn,'( 2(f5.1,a), f5.2,a, f6.2,a, f9.7,a, f5.2 )') &
   tmp_air,',',tmp_soil,',',cloud,',',prec,',',humid,',',wind

!*** 3rd line (Physical status of air)
   write (Fn,'( f6.1,a, 2(f8.3,a), f6.3,a, f6.3 )') &
   ap,',',vp_sat,',',vp,',',dnsa,',',slope_vps

!*** 4th (Radiation properties on vegetation surface)
   write (Fn,'( 5(f4.2,a) )') &
   albedo_mean,',',albedo_soil,',',albedo_leaf,',',ir_tree,',',ir_grass

!*** 5th line (Physical status of radiation)
   write (Fn,'( f6.1,a, f6.1,a, f6.1,a, f6.1,a, f6.1,a, f6.1 )') &
!   par_direct,',',par_diffuse,',',par*sum(par_grass_rel(:,:))/DivedG/DivedG,',',rad,',',radnet_soil,',',radnet_veg !!!>>>>>>>>>>>>TN:rm
   par_direct,',',par_diffuse,',',par*sum(par_grass_rel(:,:))/GRID%N_tot,',',rad,',',radnet_soil,',',radnet_veg !!!<<<<<<<<<<<<TN:add

!*** 6th line (Hydrogical varibales)
   write (Fn,'( f7.2, 3(a,f6.1), a,f8.1 )') &
   canopy_cond,',',sum(pool_w(1:5)),',',sum(pool_w(6:10)),',',sum(pool_w(11:20)),',',pool_snow
   
!*** 7th~ lines (PFT specific variables)
   !Write each PFT data
   tmp0(:)=0.0 ; tmp1(:)=0.0 ; tmp2(:)=0.0 ; tmp3(:)=0.0
   Do no=1, Max_no
      if ( .not. tree_exist(no)     ) cycle
      if ( .not. phenology(pft(no)) ) cycle
      if ( la(no) <=0.0             ) cycle
      
      p = pft(no)
      tmp0(p) = tmp0(p) + 1
      tmp1(p) = tmp1(p) + lue   (no)
      tmp2(p) = tmp2(p) + co2cmp(no)
      tmp3(p) = tmp3(p) + psat  (no)
   End do
   
   Do p=1, PFT_no
      tmp1(p) = tmp1(p) / Max(1.0, tmp0(p))
      tmp2(p) = tmp2(p) / Max(1.0, tmp0(p))
      tmp3(p) = tmp3(p) / Max(1.0, tmp0(p))
   End do
   
!   i = int( sum(par_grass_rel(:,:)) / DivedG / DivedG ) !!!>>>>>>>>>>>>TN:rm
   i = int( sum(par_grass_rel(:,:)) / GRID%N_tot ) !!!<<<<<<<<<<<<TN:add！
   i = max(1,i)
   if (pft_exist(C3g_no)) then
      tmp1(C3g_no) = lue_grass   (i) ; tmp1(C4g_no) = 0.0
      tmp2(C3g_no) = co2cmp_grass(i) ; tmp2(C4g_no) = 0.0
      tmp3(C3g_no) = psat_grass  (i) ; tmp3(C4g_no) = 0.0
   else
      tmp1(C4g_no) = lue_grass   (i) ; tmp1(C3g_no) = 0.0
      tmp2(C4g_no) = co2cmp_grass(i) ; tmp2(C3g_no) = 0.0
      tmp3(C4g_no) = psat_grass  (i) ; tmp3(C3g_no) = 0.0
   endif
   
   Do p=1, PFT_no
   if (pft_exist(p)) then
      write (Fn,'( f5.2, a, f6.3, a, f7.3, a, f6.1 )') &
      lai_RunningRecord(1,p),',',tmp1(p),',',tmp2(p),',',tmp3(p)
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
   
   write (Fn, '(i5, a, i2, a, i2, a, f8.1, a, i2)')  &
         year,',',Month(doy),',',biome,',',tree_density,',',type_grass
   
!*** 2nd line (Carbon strage & flux)
   !Sumup carbon properties
   i  = Day_in_month(Month(doy))
!   a1 = sum(flux_c_uptake_RR (1:i)) * C_in_drymass /Max_loc/Max_loc/ 100.0 !carbon uptake   (Mg C/month/ha) !!!>>>>>>>>>>>>TN:rm
   a1 = sum(flux_c_uptake_RR (1:i)) * C_in_drymass /real(GRID%Area)/ 100.0 !carbon uptake   (Mg C/month/ha) !!!<<<<<<<<<<<<TN:add 単位要チェック
   a2 = ( sum(flux_c_mnt_RR(1:i)) + sum(flux_c_mnt_RR(1:i)) + sum(flux_c_htr_RR(1:i)) + sum(flux_c_fir_RR(1:i)) ) &
!                                    * C_in_drymass /Max_loc/Max_loc/ 100.0 !carbon emission (Mg C/month/ha) !!!>>>>>>>>>>>>TN:rm
                                    * C_in_drymass /real(GRID%Area)/ 100.0 !carbon emission (Mg C/month/ha) !!!<<<<<<<<<<<<TN:add 単位要チェック
   
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:rm
!   b1 = (pool_litter_trunk + pool_litter_leaf + pool_litter_root + pool_litter_ag + pool_litter_bg) &
!                      * C_in_drymass / Max_loc / Max_loc / 100.0 !Litter  (Mg C/ha)
!   b2 = pool_som_int  * C_in_drymass / Max_loc / Max_loc / 100.0 !SOM int (Mg C/ha)
!   b3 = pool_som_slow * C_in_drymass / Max_loc / Max_loc / 100.0 !SOM slow(Mg C/ha)
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:rm
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:add
   b1 = (pool_litter_trunk + pool_litter_leaf + pool_litter_root + pool_litter_ag + pool_litter_bg) &
                      * C_in_drymass / real(GRID%Area) / 100.0 !Litter  (Mg C/ha)
   b2 = pool_som_int  * C_in_drymass / real(GRID%Area) / 100.0 !SOM int (Mg C/ha)
   b3 = pool_som_slow * C_in_drymass / real(GRID%Area) / 100.0 !SOM slow(Mg C/ha)
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:add
   
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
   
   write (Fn, '( 2(f6.1,a), 4(f8.1,a), f6.1 )') &
         a1,',',a2,',',b1,',',b2,',',b3,',',c1,',',c2
   
!*** 3rd line (Water flux)
   write (Fn, '( f8.2, a, f8.2, a, f8.2, a, f8.2, a, f8.2 )') &
      sum(  prec_RunningRecord(1:Day_of_Month(doy))),',', &
      sum(flux_ro_RunningRecord(1:Day_of_Month(doy))),',', &
      sum(flux_ic_RunningRecord(1:Day_of_Month(doy))),',', &
      sum(flux_ev_RunningRecord(1:Day_of_Month(doy))),',', &
      sum(flux_tr_RunningRecord(1:Day_of_Month(doy)))       
   
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
   
END IF

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
   
   integer tree_counter   !tree counter
   integer no,i           !loop counters
   
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
   
!Others
   pool_litter = pool_litter_trunk + pool_litter_leaf + pool_litter_root + pool_litter_ag + pool_litter_bg
   nee         = sum(flux_c_uptake_RR(:)) &
               - sum(flux_c_mnt_RR(:)) - sum(flux_c_gro_RR(:)) - sum(flux_c_htr_RR(:)) - sum(flux_c_fir_RR(:))
   
!_____________ Write title and data ***
   if (year==1) then
   write (Fn, '(82a)', advance='no') &
   'Year, TreeDensity, Heigh, DBH, BA, MassW, MassW_ag, MassG, Litter, SOM, GPP, NPP, '
   write (Fn, '(76a)') &
   'NEE, LAImax_W, LAImax_G, LAImean_W, LAImean_G, ccon, RO, Biome, ALD, ALD_DOY'
   end if
   
   write (Fn,'( i4,a, 2(f7.1,a), 2(f9.3,a), 5(f6.1,a), 2(f5.1,a), f6.1,a, 4(f5.1,a), f7.3,a, f6.1,a, i3 )') &
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
   sum(flux_ro_RunningRecord(:))             , ',', & !19 Runoff (mm year-1)
   biome                                              !20 Biome type
   
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
SUBROUTINE output_climate (Fn, cloud, prec, humid, wind, tmp_air, tmp_soil)

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
   real,intent(IN)::humid     ! humidity (kg kg-1)
   real,intent(IN)::wind      ! wind velocity (m s-1)
   
!_____________ Main part
!Write title
   if (year == 1 .and. doy==1) then
   write (Fn,*) &
   '  yr  doy    air   soil  cloud    prec   humid  vp/vpsat   wind'
   end if
   
!Write data
   write (Fn,'(2i5, 3f7.1, f8.2, f9.7, f8.3, f7.1)') &
   year, doy, tmp_air, tmp_soil, cloud, prec, humid, vp/vp_sat, wind
   
END SUBROUTINE output_climate



!**************************************************************************************************
! Output air variables
!**************************************************************************************************
SUBROUTINE output_radiation (Fn)

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
   integer,intent(IN)::Fn   !File I/O number
   
!_____________ Main part
!Write title
   if (year == 1 .and. doy==1) then
   write (Fn,*) &
   ' Yr doy  s_hgt dlen  rad    par     par_grass'
   end if
   
!Write data
   write (Fn,'(2i4, 2f6.1, 1f7.1, 2f8.1)') &
!   year, doy, sl_hgt(doy), dlen(doy), rad, par, par*sum(par_grass_rel(:,:))/DivedG/DivedG !!!>>>>>>>>>>>>TN:rm
   year, doy, sl_hgt(doy), dlen(doy), rad, par, par*sum(par_grass_rel(:,:))/real(GRID%N_tot) !!!<<<<<<<<<<<<TN:add
   
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
     'Yr, doy, albedo_soil, albedo_leaf, albedo_mean, radnet_soil, radnet_veg, radnet_long'
   end if
   
!Write data
   write (Fn,'(i4,a)'   , advance='no') year       , ','
   write (Fn,'(i3,a)'   , advance='no') doy        , ','
   write (Fn,'(f5.2,a)' , advance='no') albedo_soil, ','
   write (Fn,'(f5.2,a)' , advance='no') albedo_leaf, ','
   write (Fn,'(f5.2,a)' , advance='no') albedo_mean, ','
   write (Fn, 1         , advance='no') radnet_soil, ','
   write (Fn, 1         , advance='no') radnet_veg , ','
   write (Fn, 1         , advance='no') radnet_long
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
   aet = flux_tr_RunningRecord(1) + flux_ev_RunningRecord(1) + flux_ic_RunningRecord(1)
   
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
   write (Fn, 1, advance="no") flux_ro_RunningRecord(1), ','
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
      
      write (Fn,'(f8.5)') ld2(10,l2)
      
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

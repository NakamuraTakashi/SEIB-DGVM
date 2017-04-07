!**************************************************************************************************
! Pre-simulation procedure; initialization valuables
!**************************************************************************************************
SUBROUTINE init_value (W_fi, tmp_air, tmp_soil, prec)

!_____________ Set variables
!Namespace
   USE data_structure
   USE vegi_status_current1
   USE vegi_status_current2
   USE grid_status_current1
   USE grid_status_current2
   implicit none
   
!Augments
   real                               ,intent(IN)::W_fi      !field capacity  (m3/m3, 0.0 -> 1.0)
   real,dimension(Day_in_Year        ),intent(IN)::tmp_air   !2m air temperature (Celcius)
   real,dimension(Day_in_Year,NumSoil),intent(IN)::tmp_soil  !soil temperature (Celcius)
   real,dimension(Day_in_Year        ),intent(IN)::prec      !precipitation (mm day-1)
   
!Local variables
   integer i,j,p
   
!________ Initialize Parameters
   !PFT number of C3 and C4 grass
   do p=1, PFT_no
     if (Life_type(p)==3) C3g_no = p
     if (Life_type(p)==4) C4g_no = p
   end do
   
!________ Initialize vegi_status_current1
   !Variables about whole vegetation
   biome             = 12 !biome type (12->desert)
   lai            = 0.0 !Total LAI @ update every day            (m2/m2)
   
   !Variables for each PFT
   do p=1, PFT_no
      if (Phenology_type(p)==0 .or. Life_type(p)==3 .or. Life_type(p)==4) then
         phenology(p) =.true.
      else
         phenology(p) =.false.
      endif
   end do
   
   pft_exist(:)      = .false.
   pft_exist(C4g_no) = .true.
   
   dfl_leaf_onset (:) = 0
   dfl_leaf_shed  (:) = 0
   stat_water     (:) = 0.0
   
   !Variables for for each individual tree
   tree_exist(:)         = .false. 
   pft(:)                = 1
   age(:)                = 0
   height(:)             = 0
   bole(:)               = 0
   height_limit(:)       = 0
   flag_suppress(:)      = 0
   dbh_heartwood(:)      = 0.0
   dbh_sapwood(:)        = 0.0
   crown_diameter(:)     = 0.0
   crown_area(:)         = 0.0
   bole_x(:)             = 0.0
   bole_y(:)             = 0.0
   crown_x(:)            = 0.0
   crown_y(:)            = 0.0
   radius_limit(:)       = 0.0
   la(:)                 = 0.0
   mass_leaf(:)          = 0.0
   mass_trunk(:)         = 0.0
   mass_root(:)          = 0.0
   mass_stock(:)         = 0.0
   mass_available(:)     = 0.0
   mort_regu1(:)         = 0.0
   mort_regu2(:)         = 0.0
   mort_regu4(:)         = 0.0
   npp_crowntop   (:)    = 0.0
   npp_crownbottom(:,:)  = 0.0
   
   !Status of grass PFT (Biomass density approach)
   gmass_leaf     (:,:) = 0.0
   gmass_root     (:,:) = 0.0
   gmass_available(:,:) = 0.0
   gmass_stock    (:,:) = 50.0
   lai_grass      (:,:) = gmass_leaf(1,1) * SLA(C4g_no)
   lai_opt_grass_RunningRecord(:,:,:) = 0.0
   
   !Running Records of plant properties
   do p=1, PFT_no
      phenology_RunningRecord(:,p) = phenology(p)
   end do
   npp_RunningRecord          (:,:) = 1.0
   gpp_RunningRecord          (:,:) = 1.0
   stat_water_RunningRecord   (:,:) = 0.0
   
   mass_sum_RunningRecord     (:,:) = 0.0
   lai_RunningRecord          (:,:) = 0.0
   
!________ Initialize grid_status_current1
   dfl_fire    = 0  
   fire_number = 0
   
   tmp_coldest_20yr_ave = 0.0
   tmp_hottest_20yr_ave = 0.0
   
   gdd_20yr_ave = 0.0
   do i = 1, Day_in_Year
      gdd_20yr_ave = gdd_20yr_ave + max(0.0, tmp_air(i)-5.0)
   end do
   
   pool_litter_trunk= 0.2 * Max_loc * Max_loc
   pool_litter_leaf = 0.2 * Max_loc * Max_loc
   pool_litter_root = 0.2 * Max_loc * Max_loc
   pool_litter_ag   = 0.2 * Max_loc * Max_loc
   pool_litter_bg   = 0.2 * Max_loc * Max_loc
   pool_som_int     = 1.0 * Max_loc * Max_loc
   pool_som_slow    = 1.0 * Max_loc * Max_loc
   
   pool_fuel_standT = pool_litter_leaf
   pool_fuel_standG = pool_litter_ag
   
   pool_w(:) = 1.0 * W_fi * Depth !Assume soil water saturated at the beggining
   pool_snow = 0.000
   
   gdd_20yr_RunningRecord    (:) = gdd_20yr_ave
   tmp_coldest_RunningRecord (:) = 0.0
   tmp_hottest_RunningRecord (:) = 0.0
   tmp_ave_GrassGrowth_RR    (:) = 0.5 * (Topt0(C3g_no)+Topt0(C4g_no))
   
   pool_c_RR(:) = 0.0
   
   flux_c_uptake_RR(:) = 0.0
   flux_c_mnt_RR   (:) = 0.0
   flux_c_gro_RR   (:) = 0.0
   flux_c_htr_RR   (:) = 0.0
   flux_c_fir_RR   (:) = 0.0
   
   flux_ro_RunningRecord(:) = 0.0
   flux_ic_RunningRecord(:) = 0.0
   flux_ev_RunningRecord(:) = 0.0
   flux_tr_RunningRecord(:) = 0.0
   
   do i = 1, Day_in_Year
   tmp_air_RunningRecord (i) = tmp_air  (i)
   prec_RunningRecord    (i) = prec     (i)
      do j = 1, NumSoil
      tmp_soil_RunningRecord(i,j) = tmp_soil(i,j)
      enddo
   end do
   
   ev_pot_RunningRecord             (:)     = 0.0
   par_RunningRecord                (:)     = 0.0
   pool_w1_RunningRecord(:) = pool_w(1)
   
END SUBROUTINE init_value

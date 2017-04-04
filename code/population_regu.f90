!**************************************************************************************************
! Establishment of woody PFTs
!**************************************************************************************************
SUBROUTINE establish (GlobalZone)

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
   
!Local parameter
   real,parameter::Thres_stat_water_dry = 0.3 !Threshold stat_water to determine dry month
   
!Argument
   integer,intent(IN):: GlobalZone !ID number of global zone
   
!Local variables
   real   ,dimension(PFT_no)::est_possibility  !see above for definition
   logical,dimension(PFT_no)::est_capacity     !see above for definition
   real   ,dimension(PFT_no)::mass_sum         !total biomass of examining tree (DM g)
   integer                    no_possible      !number of woody PFTs that can establish here
   integer                    dry_month        !
   integer                    yr  !Simulation year from initialization (including spinup period)
   
   integer est_scenario_current !establishment scenario for current year
   integer new                  !number of new sapling
   integer i, j, k, p, no       !for loop counter
   real    x, y                 !for general usage
   
!Definitions of establishment control parameters
!   P_establish
!      Actual establishment rate when a establishment patch was assigned to a PFT
!   
!   Est_pft_OnOff:
!      establishment switch (for Scenerio 1)
!      .true.  -> allow establish
!      .false. -> do not allow establish
!   
!   Est_year_change
!      lenght of year until establishment pattern will change (for Scenerio 3, 4)
!   
!   Est_frac_random
!      Fraction of random establishment after Est_year_change (for Scenerio  4)

!Definitions of establishment regulator (PFT specific variables)
!   est_capacity
!      It indicates wheather a PFT can establish or not on this site (true or false)
!   
!   est_possibility
!      Frequency that a establishment patch is allocated to a PFT
!      Thus, sum of est_possibility for all PFTs is always 1.0


!_____________ calculate year of this simulation
   yr = int( (counter-1)/Day_in_Year ) + 1

!_____________ Determine to prceed establishment process
!If tree number is maximum, omit entire subroutine
   no = 0
   do i=1, Max_no
     if (tree_exist(i)) no=no+1
   end do
   if (no >= Max_no) return
   
!If fire occured within 10 years, omit entire subroutine
!   if ( dfl_fire <= 10*Day_in_Year ) return

!_____________ Compute variables that control establishment
!dry_month: Number of dry month based on Priestley-Taylor model
   dry_month = 0
   k = Day_in_Year + 1
   do i=1, 12                 !i: month
      x = 0.0
      y = 0.0
      do j=1, Day_in_month(i)            !j: day form the begginig of this month
         k = k - 1                       !k: doy
         x = x + prec_RunningRecord  (k) !x: precipitation                 (mm/month)
         y = y + ev_pot_RunningRecord(k) !y: potential evapotranspirationn (mm/month)
      enddo
      if (x<y) dry_month = dry_month + 1
   enddo
   
!est_capacity(1:PFT_no): establishment potential for each PFT based on climatic range
   est_capacity(:) = .true.
!   x = sum(sum_par_floor(:,:)) / Dived / Dived / Day_in_Year !!!>>>>>>>>>>>>TN:rm
   x = sum(sum_par_floor(:,:)) / real(GRID%N_tot) / Day_in_Year !!!<<<<<<<<<<<<TN:add
   do p = 1, PFT_no
      if ( p==C3g_no .or. p==C4g_no                  ) then; est_capacity(p)=.false. ;endif
      if ( tmp_coldest_20yr_ave >=   TC_max  (p)     ) then; est_capacity(p)=.false. ;endif
      if ( gdd_20yr_ave         >=   GDD_max (p)     ) then; est_capacity(p)=.false. ;endif
      if ( gdd_20yr_ave         <=   GDD_min (p)     ) then; est_capacity(p)=.false. ;endif
      if ( dry_month            >    DM_max  (p)     ) then; est_capacity(p)=.false. ;endif
      if ( Life_type(p)==2 .and. dfl_fire>50*365     ) then; est_capacity(p)=.false. ;endif
      
      if     (GlobalZone==1) then
        !African continent
        if ( Life_type(p).ne.5 ) then; est_capacity(p)=.false. ;endif
        
      elseif (GlobalZone==2) then
        !Eastern Siberia
        if (.not. (Life_type(p)==2 .or. Life_type(p)==6)) then; est_capacity(p)=.false. ;endif
        
      else
        !Default
        if ( Life_type(p) == 5 ) then; est_capacity(p)=.false. ;endif
      endif
   
   end do
   
!mass_sum(1:PFT_no), Sum of biomass of each woody PFT, Last 20 years means
   mass_sum(:) = 0.0
   do p=1, PFT_no
      if (p==C3g_no) cycle
      if (p==C4g_no) cycle
      mass_sum(p) = sum(mass_sum_RunningRecord(1:20,p)) / 20.0
     !mass_sum(p) = sum(mass_sum_RunningRecord(1:1 ,p)) /  1.0 !Previous method
   enddo
   
!_____________ Compute establishment probability for each PFT
!Initialize variables
   est_possibility(1:PFT_no)= 0.0 !establishment probability, a scenario specific value
   no_possible              = 0   !number of PFT that establish at this site

!Determine establishment pattern of current year
Select case (Est_scenario)
   Case (0)
      est_scenario_current = 0
   Case (1)
      est_scenario_current = 1
   Case (2)
      est_scenario_current = 2
   Case (3)
!      if (yr>=Est_year_change .and. dfl_fire>=Est_year_change*Day_in_Year) then
!         est_scenario_current = 3
!      else
!         est_scenario_current = 2
!      endif
      if (yr>=Est_year_change) then
         est_scenario_current = 3
      else
         est_scenario_current = 2
      endif
   Case (4)
      if (yr>=Est_year_change) then
         est_scenario_current = 4
      else
         est_scenario_current = 2
      endif
End select

!Give establishment probability according to the establishment pattern
Select case (est_scenario_current)
   Case (1) !only specified woody PFTs can establish
      do p = 1, PFT_no
      if ( est_capacity(p) .and. Est_pft_OnOff(p) ) then
         no_possible        = no_possible + 1
         est_possibility(p) = 1.0
      endif
      enddo
      !Adjust est_possibility, so that sum of value will be 1.0
      est_possibility(:) = est_possibility(:) / real( max(1, no_possible) )
      
   Case (2) !every potentially eatablishable PFT have same chance of establishment
      do p = 1, PFT_no
      if ( est_capacity(p) ) then
         no_possible        = no_possible + 1
         est_possibility(p) = 1.0
      end if
      end do
      !Adjust est_possibility, so that sum of value will be 1.0
      est_possibility(:) = est_possibility(:) / real( max(1, no_possible) )
      
   Case (3) !in proportion of existing biomass
      do p = 1, PFT_no
      if ( est_capacity(p) ) then
         est_possibility(p) = mass_sum(p) / max( 1.0, sum(mass_sum(:)) )
         if (est_possibility(p) > 0.00) no_possible=no_possible+1
      endif
      end do
      
   Case (4) !in proportion of existing biomass, while little portion of establishment
            !was randomly selected from all potentialy establishable PFTs
      
      do p = 1, PFT_no
         if (est_capacity(p)) no_possible=no_possible+1
      enddo
      
      do p = 1, PFT_no
      if (est_capacity(p)) then
         x = Est_frac_random / max(1, no_possible)                          !fraction stochastic
         y = (1.0-Est_frac_random) * (mass_sum(p) / max(1.0,sum(mass_sum))) !fraction deterministic
         
        !Option: No stochastic establishment when tree biomass is less than 100.0g (DM/m2)
        !if ( sum(mass_sum) < 100.0 ) y=0.0
         
         est_possibility(p) = x + y
      endif
      end do
      
End select

!_____________ Establishment processes
if (no_possible == 0) return

!loop for each ground mesh
new = 1
!Do i = 1, Dived !!!>>>>>>>>>>>>TN:rm
!Do j = 1, Dived !!!>>>>>>>>>>>>TN:rm
Do i = 1, GRID%N_x !!!<<<<<<<<<<<<TN:add
Do j = 1, GRID%N_y !!!<<<<<<<<<<<<TN:add
   
   !K: determin PFT of newly establish tree
   x = randf()
   k = 1
   do while ( sum(est_possibility(1:k)) < x .and. k<PFT_no)
      k = k+1
   end do
   
   !x: adjusted establishment rate of PFT k
!   x = P_establish(k) * real(Max_loc**2) / real(Dived**2) !!!>>>>>>>>>>>>TN:rm
   x = P_establish(k) * real(GRID%Area) / real(GRID%N_tot) !!!<<<<<<<<<<<<TN:add
   
   !Omit establishment
   if ( .not. patch_vacant(i,j)                       ) cycle
   if ( sum_par_floor(i,j) / Day_in_Year < PAR_min(k) ) cycle
   if ( x < randf()                                   ) cycle
   
   !determine new sapling ID-number
   do while ( tree_exist(new) )
      new = new + 1
      if ( new > Max_no ) return !!!<<<<<<<<<<<<TN:add ??? bug fix ????
   end do
!   if ( new > Max_no ) return !!!>>>>>>>>>>>>TN:rm ??? bug fix ????
   
   !initiate new sapling
   pft           (new) = k
   tree_exist    (new) = .true.
   
!   bole_x        (new) = real(Max_loc) * (real(i)-0.5) / real(Dived) !!!>>>>>>>>>>>>TN:rm
!   bole_y        (new) = real(Max_loc) * (real(j)-0.5) / real(Dived) !!!>>>>>>>>>>>>TN:rm
   bole_x        (new) = real(GRID%Max_x) * (real(i)-0.5) / real(GRID%N_x) !!!<<<<<<<<<<<<TN:add
   bole_y        (new) = real(GRID%Max_y) * (real(j)-0.5) / real(GRID%N_y) !!!<<<<<<<<<<<<TN:add
   
   crown_x       (new) = bole_x (new)
   crown_y       (new) = bole_y (new)
   
   age           (new) = 1
   dbh_sapwood   (new) = 0.01   ! <--- this value determine initial size of new sapling
   dbh_heartwood (new) = 0.00
   height        (new) = h(dbh_sapwood(new) + dbh_heartwood(new), pft(new))
   bole          (new) = 0
   crown_area    (new) = 0.06  ! <--- this value determine initial size of crown area
   crown_diameter(new) = 2 * sqrt(crown_area(new)/PI)
   la            (new) = 0.0   !alternative:  2.0 * crown_area(new) * height(new) * STEP
   mort_regu1    (new) = 0.0
   mort_regu2    (new) = la(new) / Day_in_Year
   mass_leaf     (new) = 0.0   !alternative:  la(new) / SLA(pft(new))
   mass_root     (new) = 0.0   !alternative:  mass_leaf(new) / FR_ratio(pft(new))
   mass_trunk    (new) = stem_weight(dbh_sapwood(new)+dbh_heartwood(new), height(new), pft(new))
   mass_stock    (new) = 250.0
   mass_available(new) = 250.0
   
   height_limit  (new) = height(new)
   radius_limit  (new) = sqrt(crown_area(new)/PI)
   
   npp_crowntop     (new)   = 0.00
   npp_crownbottom  (new,:) = 0.00
   
   !Adjust carbon balance
   x                 = min(mass_leaf     (new), pool_litter_leaf )
   flux_litter_leaf  = flux_litter_leaf - x
   
   x                 = min(mass_root     (new), pool_litter_root )
   flux_litter_root  = flux_litter_root  - x
   
   x                 = min(mass_stock    (new), pool_litter_trunk)
   flux_litter_trunk = flux_litter_trunk - x
   
   x                 = min(mass_available(new), pool_litter_trunk)
   flux_litter_trunk = flux_litter_trunk - x
   
   x                 = min(mass_trunk    (new), pool_litter_trunk)
   flux_litter_trunk = flux_litter_trunk - x
   
End do
End do

END SUBROUTINE establish



!**************************************************************************************************
! Change priority grass PFT
!**************************************************************************************************
SUBROUTINE grass_priority ()

!_____________ Set variables
!Namespace
   USE data_structure
   USE vegi_status_current1
   USE vegi_status_current2
   USE grid_status_current1
   USE grid_status_current2
   implicit none
   
!Local variables
   integer i, j, dom, mon
   real    x, y
   logical l, flag_swap
   logical,dimension(12):: flag_growth_month
   real   ,dimension(12):: mean_temperature_month
   
!_____________ Prepare variables
   !mean_temperature_month(1:12) : Monthly mean air temperature
   !flag_growth_month     (1:12) : growth month flag
   
   mean_temperature_month(:) = 0.0
   flag_growth_month     (:) = .false.
   
   i = 0
   Do mon= 1, 12
      x = 0.0
      y = 0.0
      Do dom  = 1, Day_in_month(mon)
         i = i + 1                                      !i: day of the year
         x = x + tmp_air_RunningRecord(Day_in_Year-i+1) !Aggregates air temperature (Celcius)
         y = y + prec_RunningRecord   (Day_in_Year-i+1) !Aggregates precipitation   (mm/day) 
      End do
      mean_temperature_month(mon) = x / real(Day_in_month(mon))
      
      if ( mean_temperature_month(mon)>5.0 .and. y>25.0 ) then
         flag_growth_month(mon) = .true.
      endif
   End do
   
   !x: partial pressure of atomospheric CO2 (Pa)
   x = co2atm * 101325.0 / 1000000.0 
   
   !y: Cross over temperature for C3 and C4 dominates (Celcius)
   if (x>13.0) then ; y = 1.0 / ( 1.0/(x-10.0) + 1/68.0 )
   else             ; y = 0.0
   endif                     
   
!_____________ Determine which PFT will dominate
   i = 0; !i:Number of C4 advantageous C4 month
   j = 0; !j:Number of grass growing month
   Do mon= 1, 12
      if (.not. flag_growth_month(mon)) cycle
      if (mean_temperature_month(mon)>y) i = i + 1
                                         j = j + 1
   End do
   
   flag_swap=.false.
   if (i>=(j-i)) then
   !under C4 dominant environment
      if (pft_exist(C3g_no)) flag_swap=.true.
      pft_exist(C3g_no) = .false.
      pft_exist(C4g_no) = .true.
   else
   !under C3 dominant environment
      if (pft_exist(C4g_no)) flag_swap=.true.
      pft_exist(C3g_no) = .true.
      pft_exist(C4g_no) = .false.
   endif
   
!_____________ Swap dominant Grass PFT if necessary
   if (flag_swap) then
      l=phenology     (C3g_no) ; phenology     (C3g_no)=phenology     (C4g_no) ; phenology     (C4g_no)=l
      i=dfl_leaf_onset(C3g_no) ; dfl_leaf_onset(C3g_no)=dfl_leaf_onset(C4g_no) ; dfl_leaf_onset(C4g_no)=i
      i=dfl_leaf_shed (C3g_no) ; dfl_leaf_shed (C3g_no)=dfl_leaf_shed (C4g_no) ; dfl_leaf_shed (C4g_no)=i
   end if
   
END SUBROUTINE grass_priority



!**************************************************************************************************
! Fire regime (for each year calculation)
!**************************************************************************************************
Subroutine fire_regime (W_fi)

!_____________ Set variables
!Namespace
   USE data_structure
   USE time_counter
   USE vegi_status_current1
   USE vegi_status_current2
   USE grid_status_current1
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none
   
!Local parameter
   !Fraction of litter moisture to soil moisture at the top layer
   real,parameter::Frac_Moisture_litter = 0.5
   
   !minimum stock mass of grass layer after fire (g DM/m2)
   real,parameter::Min_gmass_stock = 50.0
   
   !Fraction of trunk biomass that exist above ground
   real,parameter::Frac_trunk_AG = 0.7
   
   !Fraction of Trunk (and its litter) lost when fire occurs (0.0~1.0)
   real,parameter::Frac_TrunkLostAtFire = 0.5
   
!Augment
   real,intent(IN):: W_fi        !field capacity   (m3/m3, 0.0 -> 1.0)
   
!Local variables
   real moist_factor !adjusted moisture factor for fire spread
   real fire_factor  !factor for fire probability calculation
   
   real,dimension(PFT_no)::mass_sum     !total aboveground biomass for each PFT  (gDM/stand)
   real                    mass_overall !total aboveground biomass in this stand (gDM/stand)
   real                    mass_combust !carbon conbusted fue to fire (gDM/stand)
   
   integer no          !loop counter for individual trees
   integer i, j        !loop counter
   integer p           !PFT numer
   real    fire_prob   !probability of fire occurance in a year
   real    x, y        !for general usage
   
!_____________ Return when fire occured within a year
   if (dfl_fire < Day_in_Year*0.5) return !minimum fire interval is assumed to be half year
   
!_____________ Return when this site does not have satisfies minimum-fuel-load
   !x: total fuel load in this forest stand (g C / m2)
   x = Frac_trunk_AG * pool_litter_trunk + pool_litter_leaf + pool_litter_ag
   Do no=1, Max_no
      if ( .not. tree_exist(no) ) cycle
      x = x + Frac_trunk_AG * mass_trunk(no) + mass_leaf(no) + mass_available(no)
   End do
   x = x + sum(gmass_leaf(:,:)) + sum(gmass_available(:,:))
!   x = x * C_in_drymass / (Max_loc * Max_loc) !!!>>>>>>>>>>>>TN:rm
   x = x * C_in_drymass / real(GRID%Area) !!!<<<<<<<<<<<<TN:add
   if ( x < Fuel_min  ) return
   
!_____________ moist_factor calculation
   !Reset local variables
   mass_sum(:)  = 0.0
   mass_overall = 0.0
   moist_factor = 0.0
   
   !mass_total(PFT_no): total above ground biomass for each PFT
      !woody PFTs (Assumed that 70% biomass of trunk exists aboveground)
      Do no = 1, Max_no
      if (tree_exist(no)) then
         p = pft(no)
         mass_sum(p) = mass_sum(p) + &
                       mass_leaf(no) + Frac_trunk_AG*mass_trunk(no) + mass_available(no)
      end if
      End do
      
      !herbaceous PFTs
      if (pft_exist(C3g_no)) then; p=C3g_no
      else                       ; p=C4g_no
      endif
      mass_sum(p) = sum(gmass_leaf(:,:)) + sum(gmass_available(:,:))
   
   !mass_overall: total biomass in this grid
   mass_overall = sum(mass_sum) !grand total biomass in this forest stand
   mass_overall = max(0.01, mass_overall)
   
   !moist_factor: adjusted litter moisture weighting factor (0.01~1.0)
   Do p = 1, PFT_no
      moist_factor = moist_factor + Moisture_extinction(p) * (mass_sum(p) / mass_overall)
   End do
   moist_factor = min(1.0, max(0.01, moist_factor))
   
!_____________ Fire_factor (for disturbance computation)
   x = Frac_Moisture_litter * pool_w1_RunningRecord(1) / (Depth*W_fi)
   fire_factor = exp( -PI * ((x/moist_factor)**2) )
   
!_____________ Fire_prob: calculate fire probability [A(s) in eq.37]
   !Note: fire_prob and fire_factor are very close
   x         = fire_factor - 1.0
   fire_prob = fire_factor * EXP( x / (0.45*x*x*x + 2.83*x*x + 2.96*x + 1.04) )
   fire_prob = fire_prob / Day_in_Year !Convert fire probability from year-1 to day-1
   
   if ( fire_prob <= randf()    ) return 
   
!_____________ When fire occurs
   if (Logging==1) write(*,*) 'Fire occured!, ', counter
   dfl_fire     = 0
   fire_number  = fire_number + 1
   mass_combust = 0.0
   
   !Grass PFTs -> total above ground biomass lost
   mass_combust= mass_combust + sum(gmass_leaf(:,:))      ; gmass_leaf     (:,:)= 0.0
   mass_combust= mass_combust + sum(gmass_available(:,:)) ; gmass_available(:,:)= 0.0
   
!   x = Min_gmass_stock * ((Max_loc/DivedG)**2) !x: minumum grass stock biomass (gDM/cell) !!!>>>>>>>>>>>>TN:rm
   x = Min_gmass_stock * (real(GRID%Area)/real(GRID%N_tot)) !x: minumum grass stock biomass (gDM/cell) !!!<<<<<<<<<<<<TN:add
!   do i=1, DivedG !!!>>>>>>>>>>>>TN:rm
!   do j=1, DivedG !!!>>>>>>>>>>>>TN:rm
   do i=1, GRID%N_x !!!<<<<<<<<<<<<TN:add
   do j=1, GRID%N_y !!!<<<<<<<<<<<<TN:add
      y                = max(0.0, x - gmass_stock(i,j) ) !y: shortfall stock biomass (gDM/cell)
      gmass_stock(i,j) = max(y  , gmass_stock(i,j)     ) !Assure minimu stock grass mass for each grass cell
      pool_litter_ag   = max(0.0, pool_litter_ag - y   ) !This mass will be taken from aboveground grass litter
      pool_fuel_standG = max(0.0, pool_litter_ag - y   ) !This mass will be taken from aboveground grass litter
   end do
   end do
   
   !Woody PFT: A fraction of trees were killed, and their above ground biomass lost
   Do no = 1, Max_no
   if (tree_exist(no)) then
   if ( randf() <  M3(pft(no)) ) then
      !for dead trees
      mass_combust  = mass_combust + mass_leaf(no) + mass_stock(no) + mass_available(no)
      mass_leaf     (no) = 0.0
      mass_stock    (no) = 0.0
      mass_available(no) = 0.0
      
      mass_combust       = mass_combust      + mass_trunk(no) * Frac_TrunkLostAtFire
      flux_litter_trunk  = flux_litter_trunk + mass_trunk(no) * (1.0 - Frac_TrunkLostAtFire)
      mass_trunk(no)     = 0.0
      
      flux_litter_root   = flux_litter_root + mass_root(no)
      mass_root(no)      = 0.0
      
      tree_exist(no)     = .false.
      pft(no)            = 0
      age(no)            = 0
      height(no)         = 0
      bole(no)           = 0
      dbh_sapwood(no)    = 0.0
      dbh_heartwood(no)  = 0.0
      crown_area(no)     = 0.0
      crown_diameter(no) = 0.0
      la(no)             = 0.0
   else
      !for survived trees
      mass_combust  = mass_combust + mass_leaf(no)
      mass_leaf(no) = 0.0
      la(no)        = 0.0
   end if
   end if
   End do
   
   !Deciduous PFTs -> change to dormance phase
   Do i=1, PFT_no
   if ( Phenology_type(i) .ne. 0.0 ) then 
      phenology(i) = .false.
      dfl_leaf_shed(i) = 0
   end if
   End do
   
   !Lost of litter-pool
   mass_combust = mass_combust + &
   pool_litter_leaf + Frac_TrunkLostAtFire * pool_litter_trunk + pool_litter_ag
   
   pool_litter_leaf  = 0.0
   pool_litter_trunk = pool_litter_trunk * (1.0 - Frac_TrunkLostAtFire)
   pool_litter_ag    = 0.0
   
   !Standind dead litter pool -> all lost
   pool_fuel_standT  = 0.00
   pool_fuel_standG  = 0.00
   
   !Update fire flux recorder
   flux_c_fir_RR(1) = flux_c_fir_RR(1) + mass_combust
   
End subroutine fire_regime



!**************************************************************************************************
! Fire regime (for each year calculation)
!**************************************************************************************************
Subroutine fire_regime2 (wind)

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
   
!Local parameters
   !minimum stock mass of grass layer (g DM/m2)
   real,parameter::Min_gmass_stock= 50.0
   
   !Fraction of trunk burn when fire occcured
   real,parameter::Frac_trunk_burn = 0.50
   
!Augment
   real,intent(IN):: wind !wind velocity (m s-1)
   
!Local variables (NEW)
   real    Unit_conv1      !Unit converter ( stand-1 -> m-2      )
   real    Unit_conv2      !Unit converter ( grass_cell-1 -> m-2 )
   real    fuel            !fuel        load !(g DM/m2)
   real    fuel_dead       !fuel dead   load !(g DM/m2)
   real    fuel_live       !fuel living load !(g DM/m2)
   real    moisture_cont   !moisture content in fuel biomass
   real    prob_fire       !probability of fire occurance in a year
   real    prob_topkill    !Top kill probability
   real    mass_combust    !carbon conbusted fue to fire (gDM/stand)
   
!   real,dimension(DivedG,DivedG)::grass_avandance !relative avandance of grass (fraction) !!!>>>>>>>>>>>>TN:rm
!   real,dimension(DivedG,DivedG)::fire_potential  !potential fire intensity (KJ s-1 m-1)  !!!>>>>>>>>>>>>TN:rm
   real,dimension(GRID%N_x,GRID%N_y)::grass_avandance !relative avandance of grass (fraction) !!!<<<<<<<<<<<<TN:add
   real,dimension(GRID%N_x,GRID%N_y)::fire_potential  !potential fire intensity (KJ s-1 m-1)  !!!<<<<<<<<<<<<TN:add
   
   real    x, y            !for general usage
   integer no              !loop counter for individual trees
   integer i, j            !loop counter
   
!_____________ Return when fire occured within a year
   if (dfl_fire < int(Day_in_Year*0.5) ) return !minimum fire interval is assumed to be half year
   
!_____________ Calculate potential fire intensity
   !Unit converters
!   Unit_conv1 = 1.0 / Max_loc  / Max_loc    !(stand-1 ---> m-2) !!!>>>>>>>>>>>>TN:rm
!   Unit_conv2 = 1.0 / ((Max_loc/DivedG)**2) !(grass_cell-1 ---> m-2) !!!>>>>>>>>>>>>TN:rm
   Unit_conv1 = 1.0 / GRID%Area    !(stand-1 ---> m-2) !!!<<<<<<<<<<<<TN:add
   Unit_conv2 = 1.0 / (real(GRID%Area)/real(GRID%N_tot)) !(grass_cell-1 ---> m-2) !!!<<<<<<<<<<<<TN:add
   
   !Grass_avandance(Dived_G,Dived_G)
   !Relative grass avandance for each grass cell
!   x = sum(gmass_leaf(:,:))/ DivedG/ DivedG  !!!>>>>>>>>>>>>TN:rm
   x = sum(gmass_leaf(:,:))/ real(GRID%N_tot)  !!!<<<<<<<<<<<<TN:add
   x = max(x , 0.001)                       !x: Average ABG of grass (gDM/grass_cell)
!   do i=1, DivedG !!!>>>>>>>>>>>>TN:rm
!   do j=1, DivedG !!!>>>>>>>>>>>>TN:rm
   do i=1, GRID%N_x !!!<<<<<<<<<<<<TN:add
   do j=1, GRID%N_y !!!<<<<<<<<<<<<TN:add
      grass_avandance(i,j) = gmass_leaf(i,j) / x
   end do
   end do
   
   !fire_potential(Dived_G,Dived_G)
   !Potential fire intensity for each grass cell (KJ s-1 m-1)
   !Scheiter and Higgins (2008) with some simplification
!   do i=1, DivedG !!!>>>>>>>>>>>>TN:rm
!   do j=1, DivedG !!!>>>>>>>>>>>>TN:rm
   do i=1, GRID%N_x !!!<<<<<<<<<<<<TN:add
   do j=1, GRID%N_y !!!<<<<<<<<<<<<TN:add
      
      !fuel_live: living fuel biomass for each grass cell (gDM/m2)
      fuel_live = gmass_leaf(i,j) * Unit_conv2 + &
                  0.5 * pool_fuel_standG * grass_avandance(i,j) * Unit_conv1
      
      !fuel_dead: average standing dead fuel biomass (g/m2)
      fuel_dead = 0.5*pool_fuel_standG                          * grass_avandance(i,j) &
                + max(0.0, pool_litter_ag   - pool_fuel_standG) * grass_avandance(i,j) &
                + max(0.0, pool_litter_leaf - pool_fuel_standT)                         
      fuel_dead = fuel_dead * Unit_conv1
      
      !moisture_cont: moisture content in fuel mass 
      moisture_cont = (vp/vp_sat) * fuel_live / max(0.01, fuel_live+fuel_dead)
      
      !fuel: total fuel load
      fuel = fuel_live + fuel_dead
      
      !fire_potential: Potential fire intensity of the current cell
      x = fuel * ( fuel/(fuel+119.7) )
      y = atan(wind) * 301.0 / ( 2600.0*moisture_cont + 161000.0*(1.0-moisture_cont) )
      fire_potential(i,j) = 16890.0 * x * y
   end do
   end do
   
!_____________ Determine whether fire occurs
   !return when average potential intensity is not sufficient (Wilgen and Scholes, 1997)
!   if ( sum(fire_potential(:,:)) / DivedG / DivedG < 300.0) return  !!!>>>>>>>>>>>>TN:rm
   if ( sum(fire_potential(:,:)) / real(GRID%N_tot) < 300.0) return  !!!<<<<<<<<<<<<TN:add
   
   !prob_fire: probability of fire occurs in this year (in stand-1 day-1)
   prob_fire = 0.0015 !from aDGVM (when crown coverage is 0%)
   
   !from Archibald et al (2009)
   if     (frac_crown_coverage>0.65) then; prob_fire = 0.0
   elseif (frac_crown_coverage>0.40) then; prob_fire = prob_fire * 0.20
   endif
   
   if ( prob_fire < randf() ) return 
   
!_____________ When fire occurs
   dfl_fire     = 0               !Reset 'dfl' (Day From Last Fire)
   fire_number  = fire_number + 1 !Fire number from the beggining of the simulation
   mass_combust = 0.0             !Mass to be burned (gDM/stand)`
   
!  !Deciduous PFTs -> change to dormance phase
!   Do i=1, PFT_no
!   if ( Phenology_type(i) .ne. 0.0 ) then 
!      phenology(i) = .false.
!      dfl_leaf_shed(i) = 0
!   end if
!   End do
   
   !Grass PFTs -> change to dormance phase
   Do i=1, PFT_no
   if ( Life_type(i)==3 .or. Life_type(i)==4 ) then
      phenology(i) = .false.
      dfl_leaf_shed(i) = 0
   end if
   End do
   
   !Grass PFTs -> All above ground biomass will be burned
   mass_combust= mass_combust+ sum(gmass_leaf     (:,:)) ; gmass_leaf     (:,:)= 0.0
   mass_combust= mass_combust+ sum(gmass_available(:,:)) ; gmass_available(:,:)= 0.0
   
   !Grass PFTs -> give minimum stock biomass for recovery
!   x = Min_gmass_stock * ((Max_loc/DivedG)**2) !x: minumum grass stock biomass (gDM/cell) !!!>>>>>>>>>>>>TN:rm
   x = Min_gmass_stock * (real(GRID%Area)/real(GRID%N_tot)) !x: minumum grass stock biomass (gDM/cell) !!!<<<<<<<<<<<<TN:add
!   do i=1, DivedG !!!>>>>>>>>>>>>TN:rm
!   do j=1, DivedG !!!>>>>>>>>>>>>TN:rm
   do i=1, GRID%N_x !!!<<<<<<<<<<<<TN:add
   do j=1, GRID%N_y !!!<<<<<<<<<<<<TN:add
      !y: shortage of stock biomass (gDM/cell)
      y                = max(0.0, x - gmass_stock(i,j) )
      if (y==0.0) cycle
      
      !Assure minimum stock grass mass for each grass cell
      gmass_stock(i,j) = max(x  ,     gmass_stock(i,j) )
      
      !This mass will be taken from aboveground standing dead mass to maintain carbon balance
      y = min(y, pool_fuel_standG)
      pool_fuel_standG = pool_fuel_standG - y
      pool_litter_ag   = pool_litter_ag   - y
   end do
   end do
   
   !Litter pool consumption
   mass_combust = mass_combust + &
                  pool_litter_ag + pool_litter_leaf + pool_litter_trunk * Frac_trunk_burn
   
   pool_litter_ag    = 0.0
   pool_litter_leaf  = 0.0
   pool_litter_trunk = pool_litter_trunk * (1.0 - Frac_trunk_burn)
   
   pool_fuel_standT  = 0.0
   pool_fuel_standG  = 0.0
   
   !Woody PFTs: A fraction of trees will be killed, 
   !and above ground biomass of dead trees will be burned
   Do no = 1, Max_no
   if (tree_exist(no)) then
      
      !i, j: x and y location of the mentioned tree (Grass cells coodinate)
!      i = int( DivedG * (bole_x(no) / Max_loc) ) !!!>>>>>>>>>>>>TN:rm
!      j = int( DivedG * (bole_y(no) / Max_loc) ) !!!>>>>>>>>>>>>TN:rm
      i = int( GRID%N_x * (bole_x(no) / GRID%Max_x) ) !!!<<<<<<<<<<<<TN:add
      j = int( GRID%N_y * (bole_y(no) / GRID%Max_y) ) !!!<<<<<<<<<<<<TN:add
      
      !prob_topkill: Probability of tree death due to fire from aDGVM
      x = exp( 4.3 - 5.003*log(height(no)*STEP+1.3) + 0.004408*sqrt(fire_potential(i,j)*1000.0) )
      prob_topkill = x / (1.0+x)
      
      !for dead trees
      if (randf() < prob_topkill) then
         !carbon balance 
         mass_combust       = mass_combust        + mass_leaf  (no)
         
         x = mass_trunk (no) + mass_available(no) + mass_stock(no)
         mass_combust       = mass_combust        + x * Frac_trunk_burn
         flux_litter_trunk  = flux_litter_trunk   + x * (1.0-Frac_trunk_burn)
         flux_litter_root   = flux_litter_root    + mass_root(no)
         
         !initialize variables that was allocated for the dead tree
         tree_exist    (no) = .false.
         pft           (no) = 0
         age           (no) = 0
         height        (no) = 0
         bole          (no) = 0
         dbh_sapwood   (no) = 0.0
         dbh_heartwood (no) = 0.0
         crown_area    (no) = 0.0
         crown_diameter(no) = 0.0
         la            (no) = 0.0
         mort_regu1    (no) = 0.0
         mort_regu2    (no) = 0.0
         mass_leaf     (no) = 0.0
         mass_root     (no) = 0.0
         mass_trunk    (no) = 0.0
         mass_stock    (no) = 0.0
         mass_available(no) = 0.0
         
         height_limit  (no) = 0
         radius_limit  (no) = 0.0
         
         npp_crowntop   (no)   = 0.0 
         npp_crownbottom(no,:) = 0.0
         
      end if
   end if
   End do
   
   !Update fire flux recorder
   flux_c_fir_RR(1) = flux_c_fir_RR(1) + mass_combust
   
End subroutine fire_regime2



!**************************************************************************************************
! Tree death process due to variety of mortality reason except wild fire.
! Death by fire is coded on subroutines fire_regime and fire_regime2.
!**************************************************************************************************
SUBROUTINE mortality ()

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
   
!Local parameters (Gap formation parameters from Huth & Ditzer (2000) )
   ! bottom height of forest layer no 2,3,4,5 (m) 
   !(top    height of forest layer no 1,2,3,4 (m) )
   real,dimension(4),parameter::Layer_top = (/1.3, 15.0, 25.0, 36.0/)
   
   !Threshold fraction of crown layer overlapping to declare 'crowded layer' (frac)
   real,parameter::Frac_crowded         = 1.00
   
   !Threshold Tree height to have possibility to form gap (m)
   real,parameter::Thres_height_GapForm = 25.0
   
   !Fraction of gap formation when large tree die
   real,parameter::Frac_GapForm         = 0.20
   
   !Fraction of gap formation when large tree die (m)
   real,parameter::Dist_GapEntangled    = 11.3
   
   ! for tree layer no 1,2,3,4,5 (m)
   real,dimension(5),parameter::Frac_GapEntangled = (/0.0, 0.3, 0.6, 0.8, 0.4/)
   
!Local variable
   !Frags for death assigned trees
   logical,dimension(Max_no)::death_list     
   logical,dimension(Max_no)::death_list_gap 
   
   !Mortality component
   real mort_greff   !mortality as a result of low growth efficiency
   real mort_lim     !mortality as a result of bioclimic limits
   real mort_etc     !mortality as a result of other factors
   real mort_total   !mortality as a result of above all factors
   
   real sal          !top soil salinity (PSU)
   
!From embedding FORMIX3
!   !Crown area for each cohort (m2)
!   real   ,dimension((int(0.99999+Max_loc/20.0))**2, 5)::cohort_ca     
!   !crowded flag for each cohort
!   logical,dimension((int(0.99999+Max_loc/20.0))**2, 5)::cohort_crowded
!   integer,dimension(Max_no)::id_location
   !Frag wheather trees are crowded
   logical,dimension(Max_no)::tree_crowded
   !crown layer numner that each tree belongs
   integer,dimension(Max_no)::id_layer
   
!For general usage
   integer no, p, i, j                              !loop counter
   real    x, y                                     !for general usage
   real    d1, d2, d3, d4, d5, d6, d7, d8, d9, dist !
   
!Local variables for tree_crowded
integer me, you
real    dist1, dist2, r1, r2
real    frac_overlap_sum
real    cosine1, cosine2

!_____________ Determine cohort structure
!   cohort_ca      (:,:) = 0.0     !(location,layer) crown area for each cohort (m2)
!   cohort_crowded (:,:) = .false. !(location,layer) crowded flag for each cohort
!   id_location    (:)   = 0       !(tree_number)    location number that each tree belongs
   
   !determine height class for each tree
   Do no=1, Max_no 
   if ( .not. tree_exist(no) ) cycle
      !height_class_me: my height class
      if     ( height(no)*STEP+1.3<Layer_top(1) ) then ;id_layer(no) = 1
      elseif ( height(no)*STEP+1.3<Layer_top(2) ) then ;id_layer(no) = 2
      elseif ( height(no)*STEP+1.3<Layer_top(3) ) then ;id_layer(no) = 3
      elseif ( height(no)*STEP+1.3<Layer_top(4) ) then ;id_layer(no) = 4
      else                                             ;id_layer(no) = 5
      endif
   End do
   
   !determine crowded frag for each tree
   tree_crowded(:)=.false.
   Do me=1, Max_no 
   if ( .not. tree_exist(me) ) cycle
      
      !reset overlapping fraction
      frac_overlap_sum = 0.0
      
      Do you=1, Max_no 
      if ( .not. tree_exist(you)           ) cycle
      if ( me==you                         ) cycle
      if ( id_layer(me) .ne. id_layer(you) ) cycle
         
         !dist1: absolute y-axis distance between circle centers (adjusted for mirror world)
         dist1 = crown_y(me)-crown_y(you)
!         dist1 = min( abs(dist1), abs(dist1-Max_loc), abs(dist1+Max_loc) ) !!!>>>>>>>>>>>>TN:rm
         dist1 = min( abs(dist1), abs(dist1-real(GRID%Max_y)), abs(dist1+real(GRID%Max_y)) ) !!!<<<<<<<<<<<<TN:add
         
         !dist2: absolute x-axis distance between circle centers (adjusted for mirror world)
         dist2 = crown_x(me)-crown_x(you)
!         dist2 = min( abs(dist2), abs(dist2-Max_loc), abs(dist2+Max_loc) ) !!!>>>>>>>>>>>>TN:rm
         dist2 = min( abs(dist2), abs(dist2-real(GRID%Max_x)), abs(dist2+real(GRID%Max_x)) ) !!!<<<<<<<<<<<<TN:add
         
         !dist: x-y-plane distance between crown and shade center
         dist = sqrt( dist1**2 + dist2**2 )
         
         !sum up overlapping fraction
         x = dbh_heartwood(me)  + dbh_sapwood(me)  !dbh (me, unit=m)
         y = dbh_heartwood(you) + dbh_sapwood(you) !dbh (me, unit=m)
         
         d1 = max(25*x, (50-125*x)*x) !potential crown diameter (me , unit=m)
         d2 = max(25*y, (50-125*y)*y) !potential crown diameter (you, unit=m)
         
         r1 = min(d1,d2) / 2.0 !radius of smaller circle (m)
         r2 = max(d1,d2) / 2.0 !radius of larger  circle (m)
         if     (r1+r2 <= dist) then
            x = 0.000
            
         elseif (dist <= r2-r1) then
            x = PI * r1 * r1
            
         elseif (0 < r1**2 - r2**2 + dist**2 ) then
            cosine1 = (r1**2 - r2**2 + dist**2) / (2.0*r1*dist)
            cosine1 = min(1.0,max(-1.0,cosine1))
            cosine2 = (r2**2 - r1**2 + dist**2) / (2.0*r2*dist)
            cosine2 = min(1.0,max(-1.0,cosine2))
            x       = r1*r1*( acos(cosine1) - cosine1 * sqrt(1.0 - cosine1**2) ) &
                    + r2*r2*( acos(cosine2) - cosine2 * sqrt(1.0 - cosine2**2) )  
            
         else
            cosine1 = (r1**2 - r2**2 + dist**2) / (2.0*r1*dist)
            cosine1 = min(1.0,max(-1.0,cosine1))
            cosine2 = (r2**2 - r1**2 + dist**2) / (2.0*r2*dist)
            cosine2 = min(1.0,max(-1.0,cosine2))
            x       = r1*r1*acos(cosine1) + r2*r2*acos(cosine2)   &
                    - dist * r2 * sqrt(1.0 - cosine2**2)  
            
         end if
         
         !Sumup fraction overlap
         frac_overlap_sum = frac_overlap_sum + min(1.00, x / crown_area(me))
         
      End do
      
      !Determine whether this tree is in crowded layer
      if (frac_overlap_sum >= Frac_crowded) tree_crowded(me)=.true.
      
   End do
   
!   Do no=1, Max_no 
!   If ( tree_exist(no) ) then
!      
!      id_location(no) = 1+int(bole_x(no)/20) + 5*int(bole_y(no)/20) 
!      
!      x = height(no)*STEP+1.3 !x: tree height (m)
!      if     (x<Layer_top(1)) then ;id_layer(no) = 1
!      elseif (x<Layer_top(2)) then ;id_layer(no) = 2
!      elseif (x<Layer_top(3)) then ;id_layer(no) = 3
!      elseif (x<Layer_top(4)) then ;id_layer(no) = 4
!      else                         ;id_layer(no) = 5
!      endif
!      
!      cohort_ca (id_location(no),id_layer(no)) = &
!      cohort_ca (id_location(no),id_layer(no)) + &
!      ( 25*(dbh_heartwood(no)+dbh_sapwood(no)) )**2 * PI / 4.0
!      !crown_area(no)の代わりに、ポテンシャルの樹冠断面面積を使った
!      
!   End if
!   End do
!   
!   Do i=1, (int(0.99999+Max_loc/20.0))**2 !for each location
!   Do j=1, 5                              !for each layer
!      if ( cohort_ca(i,j)/400 > Frac_crowded ) cohort_crowded(i,j)=.true.
!   End do
!   End do
   
!_____________ Determine trees to be die 1
   death_list(:) = .false.
   
DO no=1, Max_no 
if ( tree_exist(no) .and. age(no)>1 ) then
   
   !reset mortality variables for each PFT
   mort_greff=0.0 ; mort_lim=0.0 ; mort_etc=0.0 ; mort_total=0.0
   
   !set pft number
    p = pft(no)
   
   !mortality 1: growth efficiency
   Select Case (Life_type(p))
      !tropical evergreen trees
      case (1) 
         
         !Under normal condition
         mort_greff = 0.0178 * exp(-242.57*mort_regu4(no))
         mort_greff = max(mort_greff, 0.0032)
         mort_greff = mort_greff * M4(p)
         
         !When crowded
!         if ( cohort_crowded(id_location(no),id_layer(no)) ) mort_greff=M5(p)
!         if (tree_crowded(no)) mort_greff=M5(p)
         
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add Mangrove feature ***tentative
      case (7,8) 
         sal = GRID%sal_ave(int(bole_x(no)*2.0)+1,int(bole_y(no)*2.0)+1)
         !Under normal condition
         x      = (mort_regu1(no)/1000.0) / max(0.01, mort_regu2(no))
            !x         : Annual NPP per leaf area (Kg dm m-2)
            !mort_regu1: NPP annual (g dm / individual)
            !mort_regu2: mean leaf area of last year (m2/day)
            
         x      = max(0.0, x)
         mort_greff = M1(p) * exp(Msal1(p)* ( Msal2(p) - sal))/ (M2(p)**x)
            !large M2 --> intensify grouth rate efficiency to mortality rate
         
         !When crowded
!         if ( cohort_crowded(id_location(no),id_layer(no)) ) mort_greff=M5(p)
!         if (tree_crowded(no)) mort_greff=M5(p)
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
         
      !other woody PFTs
      case default
         x      = (mort_regu1(no)/1000.0) / max(0.01, mort_regu2(no))
            !x         : Annual NPP per leaf area (Kg dm m-2)
            !mort_regu1: NPP annual (g dm / individual)
            !mort_regu2: mean leaf area of last year (m2/day)
            
         x      = max(0.0, x)
         mort_greff = M1(p) / (M2(p)**x)
            !large M2 --> intensify grouth rate efficiency to mortality rate
         
         !other method
         !x      = mort_regu1(no) / max(0.1, mort_regu2(no)) / 1000.0
         !x      = max(0.001, x)
         !mort_greff = 1.0 / ( 1.0 + M1(p) * ( x ** M2(p) ) )
         
   End select
   
   !mortality 2: bioclimic limits (based on 20yrs mean of coldest month temperature)
   if ( tmp_coldest_20yr_ave < TC_min(p) )                                       mort_lim=0.1
   if ( Life_type(p)==2 .and. (tmp_hottest_20yr_ave-tmp_coldest_20yr_ave)<43.0 ) mort_lim=0.1
   
   !mortality 3: mortality by other factors
   x = mass_leaf(no) + mass_trunk(no) + mass_root(no) + mass_stock(no) + mass_available(no)
   
   if ( height(no)-bole(no)<=2 .and. age(no)>3   ) mort_etc=0.1
   if ( dbh_heartwood(no)+dbh_sapwood(no) > 1.00 ) mort_etc=0.1
   if ( age(no) >= AGE_max(p)                    ) mort_etc=0.1
   
   Select case (Life_type(p))
      case (2)
         !Larch
         if ( mort_regu1(no)<0.0 .and. age(no)>3 ) mort_etc=mort_etc + mort_greff*9
      case (6)
         !Other Siberian woody PFTs
         if ( mort_regu1(no)<0.0 .and. age(no)>3 ) mort_etc=mort_etc + mort_greff*9
      case (5)
         !African trees
         if ( mort_regu1(no)<0.0 .and. age(no)>3 ) mort_etc=mort_etc + mort_greff*9
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add Mangrove feature ***tentative
      case (7,8)
         !Other woody PFTs
         if ( mort_regu1(no)<0.0 .and. age(no)>3 ) mort_etc=mort_etc
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
      case default
         !Other woody PFTs
         if ( mort_regu1(no)<0.0 .and. age(no)>3 ) mort_etc=1.0
   End Select
   
   !Sum up all mortality components, and calculate probabilty of death
   mort_total   = min(1.0, mort_greff + mort_lim + mort_etc)
   
   !Determine whether die
   if ( mort_total>randf() ) death_list(no) = .true.
   
end if
END DO
   
!_____________ Determine trees to be die 2 (Gap formation)
   death_list_gap(:) = .false.
   
DO no=1, Max_no 
   
   !When gap forms
   if ( (tree_exist(no)                            ) .and. &
        (death_list(no)                            ) .and. &
        (Life_type(pft(no)) == 1                   ) .and. &
        (height(no)*STEP+1.3 > Thres_height_GapForm) .and. &
        (Frac_GapForm > randf()                    )         ) then
      
      !determine center location of new gap
      x = crown_x(no) !randf() * Max_loc
      y = crown_y(no) !randf() * Max_loc
      
      !determine which tree will die due to the gap formation
      do i=1, Max_no
      if ( tree_exist(i) .and. age(i)>1 ) then
         
         !Distance from gap center
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:rm
!         d1 = (x - (bole_x(i)-real(Max_loc)) )**2 + (y - (bole_y(i)-real(Max_loc)) )**2
!         d2 = (x - (bole_x(i)              ) )**2 + (y - (bole_y(i)-real(Max_loc)) )**2
!         d3 = (x - (bole_x(i)+real(Max_loc)) )**2 + (y - (bole_y(i)-real(Max_loc)) )**2
!         d4 = (x - (bole_x(i)-real(Max_loc)) )**2 + (y - (bole_y(i)              ) )**2
!         d5 = (x - (bole_x(i)              ) )**2 + (y - (bole_y(i)              ) )**2
!         d6 = (x - (bole_x(i)+real(Max_loc)) )**2 + (y - (bole_y(i)              ) )**2
!         d7 = (x - (bole_x(i)-real(Max_loc)) )**2 + (y - (bole_y(i)+real(Max_loc)) )**2
!         d8 = (x - (bole_x(i)              ) )**2 + (y - (bole_y(i)+real(Max_loc)) )**2
!         d9 = (x - (bole_x(i)+real(Max_loc)) )**2 + (y - (bole_y(i)+real(Max_loc)) )**2
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:rm
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
         d1 = (x - (bole_x(i)-real(GRID%Max_x)) )**2 + (y - (bole_y(i)-real(GRID%Max_y)) )**2
         d2 = (x - (bole_x(i)                 ) )**2 + (y - (bole_y(i)-real(GRID%Max_y)) )**2
         d3 = (x - (bole_x(i)+real(GRID%Max_x)) )**2 + (y - (bole_y(i)-real(GRID%Max_y)) )**2
         d4 = (x - (bole_x(i)-real(GRID%Max_x)) )**2 + (y - (bole_y(i)                 ) )**2
         d5 = (x - (bole_x(i)                 ) )**2 + (y - (bole_y(i)                 ) )**2
         d6 = (x - (bole_x(i)+real(GRID%Max_x)) )**2 + (y - (bole_y(i)                 ) )**2
         d7 = (x - (bole_x(i)-real(GRID%Max_x)) )**2 + (y - (bole_y(i)+real(GRID%Max_y)) )**2
         d8 = (x - (bole_x(i)                 ) )**2 + (y - (bole_y(i)+real(GRID%Max_y)) )**2
         d9 = (x - (bole_x(i)+real(GRID%Max_x)) )**2 + (y - (bole_y(i)+real(GRID%Max_y)) )**2
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
         
         dist = min(d1,d2,d3,d4,d5,d6,d7,d8,d9)
         dist = sqrt(dist)
         
         if ( dist<Dist_GapEntangled ) then
         if ( Frac_GapEntangled(id_layer(no)) > randf() ) then
            death_list_gap(i)=.true.
         endif
         endif
         
      endif
      enddo
      
   endif
   
END DO
   
!_____________ Kill assigned trees
   do no=1, Max_no 
   if (death_list(no) .or. death_list_gap(no)) then
      
      !Increase pool_litter
      flux_litter_trunk = flux_litter_trunk + mass_trunk(no)+ mass_available(no)+ mass_stock(no)
      flux_litter_leaf  = flux_litter_leaf  + mass_leaf (no)
      flux_litter_root  = flux_litter_root  + mass_root (no)
      
      !Initialize tree
      tree_exist     (no) = .false.
      pft            (no) = 1  
      age            (no) = 0  
      height         (no) = 0  
      bole           (no) = 0  
      dbh_sapwood    (no) = 0.0
      dbh_heartwood  (no) = 0.0
      crown_area     (no) = 0.0
      crown_diameter (no) = 0.0
      la             (no) = 0.0
      mort_regu1     (no) = 0.0
      mort_regu2     (no) = 0.0
      mass_leaf      (no) = 0.0
      mass_root      (no) = 0.0
      mass_trunk     (no) = 0.0
      mass_stock     (no) = 0.0
      mass_available (no) = 0.0
      
      height_limit   (no) = 0
      radius_limit   (no) = 0.0
      
      npp_crowntop        (no)   = 0.0
      npp_crownbottom     (no,:) = 0.0
      
   endif
   end do
   
!_____________ Increament of tree age
   do no=1, Max_no 
   if ( tree_exist(no) ) then !when tree exist
       age(no) = age(no) + 1  !add tree age
   end if
   end do
   
END SUBROUTINE mortality

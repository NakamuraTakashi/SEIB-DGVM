!**************************************************************************************************
! Calculation of photosynthesis parameters (for each PFT)
!**************************************************************************************************
SUBROUTINE photosynthesis_condition (tmp_air)

!_____________ Set variables
!Namespace
   USE data_structure
   USE time_counter
   USE vegi_status_current1
   USE vegi_status_current2
   USE grid_status_current1
   USE grid_status_current2
   implicit none
   
!Define Auguments
   real,intent(IN)::tmp_air
   
!Define local variables
   real,dimension(PFT_no,MaxParClass)::lue_today    !light use efficiency (mol CO2 mol photon-1)
   real,dimension(PFT_no,MaxParClass)::co2cmp_today !CO2 compensation point (ppmv)
   real,dimension(PFT_no,MaxParClass)::psat_today   !light-saturated photosynthesis rate (micro mol CO2 m-2 s-1)
   
   real   p_ave    !average photosynthetic rate
   real   gs_ave   !daily mean stomatal conductance (mol H2O m-2 s-1)
   real   co2cell  !monthly intercellular CO2 concentration (ppmv)
   real   topt     !adjusted optimum temperature (deg C)
   real   ce_tmp   !temperature effect on photosynthesis
   real   ce_co2   !stomatal effect on photosynthesis via intercellular CO2 concentration
   real,dimension(PFT_no)::ce_water !non-stomatal effect on photosynthesis
   
   real     a1, a2, x
   integer  loop, p, i, no
   real     par_current
   real     par_ave           !average PAR intensity receiving for each tree
   
!_____________ Condition that plants can survive by Koppen (1936) (optional)
!   a1 = sum(prec_RunningRecord(1:Day_in_Year))                  !annual precipitation (mm/year)
!   a2 = sum(tmp_air_RunningRecord(1:Day_in_Year)) / Day_in_Year !annual mean temperature (Celcius)
!   if ( a1 < max(100.0, 10.0*a2) ) then
!      lue          (:) = 0.0
!      co2cmp       (:) = 0.0
!      psat         (:) = 0.0
!      lue_grass    (:) = 0.0
!      co2cmp_grass (:) = 0.0
!      psat_grass   (:) = 0.0
!      return
!   endif
   
!_____________ Preparation of PFT specific varables
   Do p=1, PFT_no
      !eK: adjusted light-attenuation-coefficient
      eK(p) = EK0(p) / min(1.0, max(0.3, sin(0.86 * sl_hgt(doy) * DtoR)))
      
      !ce_water: limitation on photosynthesis via soil water
      x = min( 1.0 , max(0.001, stat_water(p)) )
     !ce_water(p) = x ** 0.5           !Method for SEIB V1.00-V2.61
      ce_water(p) = 2.0 * x - x ** 2.0 !Method for SEIB V2.62-
   End do
   
!_____________ Compute photosynthesis hparameters for each PFT and each relative light intensity
psat_today  (:,:) = 0.0 !light-saturated photosynthesis rate (micro mol CO2 m-2 s-1)
lue_today   (:,:) = 0.0 !light use efficiency (mol CO2 mol photon-1)
co2cmp_today(:,:) = 0.0 !CO2 compensation point (ppmv)

DO p = 1, PFT_no
   if ( .not. pft_exist(p) ) cycle
!  if ( .not. phenology(p) ) cycle
   
   Do i = 1, MaxParClass !for each relative PAR intensity
      par_current = 0.5 * par * (i/real(MaxParClass)) !Daily mean PAR intensity for this loop 
      
      !initialize
      psat_today(p,i) = Pmax(p) * 0.7 !initial psat, initial light_saturated photosynthesis rate
      co2cell         = co2atm  * 0.7 !initial intercellular CO2 concentration
      
      !loop for parameter convergence
      do loop = 0, 5
         !------------------------------------------------------------!
         !lue: light use efficiency (mol CO2 mon photon-1)
         !------------------------------------------------------------!
         a1             = (52.0-tmp_air) / (3.5+0.75*(52.0-tmp_air)) !temperature dependence factor
         a2             = co2cell / (90.0+0.6*co2cell)               !CO2 dependence factor
         
                              lue_today(p,i) = Lue0(p) * a1 * a2 
         if (Life_type(p)==4) lue_today(p,i) = Lue0(p)           !temperature and CO2 insensitive for C4 grass
         
         !------------------------------------------------------------!
         ! ce_tmp: temperature effect on photosynthesis
         !------------------------------------------------------------!
         !topt: optimum photosynthesis temperature (Celcius)
         if (Life_type(p)==3)     then
            !C3 grass
            topt = sum( tmp_ave_GrassGrowth_RR(1:20) ) / 20.0
            topt = max(10.0, topt)
            topt = min(30.0, topt)
         elseif (Life_type(p)==4) then
            !C4 grass
            topt = sum( tmp_ave_GrassGrowth_RR(1:20) ) / 20.0
            topt = max(20.0, topt)
            topt = min(40.0, topt)
         else
            !Woody PFTs
            topt = Topt0(p) + 0.01 * co2cell
         endif
         
         !ce_tmp: temperature effect on photosynthesis
         x = (tmp_air-Tmax(p)) * (tmp_air-Tmin(p)) - (tmp_air-topt)**2
         if ( (abs(x)<0.0001) .or. (tmp_air<=Tmin(p)) .or. (tmp_air>=Tmax(p)) ) then
            ce_tmp=0.0
         else
            ce_tmp =  (tmp_air-Tmax(p)) * (tmp_air-Tmin(p)) / x
            ce_tmp = min(1.0, max(0.0, ce_tmp ))
         endif
         
         !------------------------------------------------------------!
         ! ce_co2: stomatal limitation on photosynthesis via co2cell
         !------------------------------------------------------------!
         ! compcd: CO2 compensation point
         if (Life_type(p) .ne. 4) then
            a1                = 1.00 + 0.0451*(tmp_air-20.0) + 0.000347*((tmp_air-20.0)**2)
            co2cmp_today(p,i) = CO2cmp0(p) * max(0.0, a1) 
            ce_co2            = 0.30 + 0.70 * (co2cell-co2cmp_today(p,i)) / (KM(p)+co2cell)
         else
            co2cmp_today(p,i) = CO2cmp0(p)
            ce_co2            = 0.50 + 0.50 * (co2cell-co2cmp_today(p,i)) / (KM(p)+co2cell)
         endif
         ce_co2            = min(1.0, max(0.0, ce_co2))
         
         !------------------------------------------------------------!
         !psat: light-saturated photosynthesis rate (micro mol CO2 m-2 s-1)
         !------------------------------------------------------------!
         psat_today(p,i)= Pmax(p) * ce_tmp * ce_co2 * ce_water(p)
         
         !------------------------------------------------------------!
         !p_ave: daily mean photosynthetic rate (micro mol CO2 m-2 s-1)
         !------------------------------------------------------------!
         a1 = psat_today(p,i) * lue_today(p,i) * par_current
         a2 = psat_today(p,i) + lue_today(p,i) * par_current
         
         if (a2 > 0.0) then; p_ave = a1/a2
         else              ; p_ave = 0.0  
         endif                            
         
         !------------------------------------------------------------!
         !gs_ave: daily mean stomatal conductance (mol H2O m-2 s-1)
         !        Leuning (1995)
         !------------------------------------------------------------!
         gs_ave = GS_b1(p) + ( GS_b2(p)*p_ave ) / ( (co2atm-co2cmp_today(p,i))*(1.0+vpd/GS_b3(p)) )
         
         !------------------------------------------------------------!
         !co2cell: intercellular CO2 concentration (ppmv)
         !------------------------------------------------------------!
         a1      = gs_ave / 1.56       !1.56: conversion from H2O to CO2 conductance
         co2cell = co2atm - p_ave / a1 !Leuning (1995)
         co2cell = max(0.0, min(co2atm, co2cell))
      end do
   End Do
END DO

!_____________ Give photosynthesis condition for each tree
!              accroding to their PFT and mean daytime PAR
   
   !x: fraction of direct PAR (0.0~1.0)
   x = par_direct / max(0.0001, par_direct+par_diffuse)
   x = max(x, 0.0)
   x = min(x, 1.0)
   
!$omp parallel private(no)
!$omp do private(p,par_ave,i)
   DO no=1, Max_no
      p = pft(no) ! p: PFT number
      
      if ( .not. tree_exist(no)) cycle
      if ( .not. phenology(p)  ) cycle
      if ( la(no) <=0.0        ) cycle
      
      !par_ave: Mean daytime PAR of each tree (micro mol photon m-2 s-1)
      par_ave = 0.0 
      do i = height(no), bole(no)+1, -1
         par_ave = par_ave  + x * par_direct_rel(no,i) + (1.0-x) * par_diffuse_rel(i)
      end do
      
      !i: MaxParClass * relative par intensity of this tree (1~MaxParClass)
      i = int( real(MaxParClass) * par_ave/(height(no)-bole(no)) ) 
      i = min(max(i,1),MaxParClass)
      
      psat  (no) = psat_today   (p,i)
      lue   (no) = lue_today    (p,i)
      co2cmp(no) = co2cmp_today (p,i)
      
   END DO
!$omp end do
!$omp end parallel
   
!_____________ Give photosynthesis condition for each grass cell
   if (pft_exist(C3g_no)) then
      p=C3g_no
   else
      p=C4g_no
   endif
   
   DO i=1, MaxParClass
      psat_grass  (i) = psat_today   (p,i)
      lue_grass   (i) = lue_today    (p,i)
      co2cmp_grass(i) = co2cmp_today (p,i)
   END DO
   
END SUBROUTINE photosynthesis_condition



!**************************************************************************************************
! Daily GPP and Canopy_conductance
!**************************************************************************************************
SUBROUTINE photosynthesis ()

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
   
!Define local variables
   real    const0, const1, const2, const3, const4, const5
   real    gpp_ind
   real    x, a1, a2, a3, a4 !for general usage
   integer i, j, k, no, p        !for loop counter
   integer count                 !counter
   
!Define Local parameter
   real Unit_converter
        Unit_converter = 3600.0 * 12.0 / 1000000.0 / C_in_drymass
   ! constant for unit conversion
   ! from (micro mol CO2 m-2 s-1) to (g dm m-2 hour-1)
   ! micro mol  -> mol       :* 10^-6
   ! mol        -> g(Carbon) :* 12
   ! secound^-1 -> hour^-1   :* 3600
   
!Initialize canopy conductance
   canopy_cond = 0.0
   
!_____________ Woody PFTs
count = 0
!$omp parallel private(no)
!$omp do private(i,gpp_ind,p,const0,const4,const5,a1,a2,a4)
DO no=1, Max_no
   
   !initialize gpp of a tree
   gpp_ind = 0.000
   p = pft(no)
   
   if ( .not. tree_exist(no) ) cycle
   if ( .not. phenology(p)   ) cycle
   if ( la(no)  <=0.0        ) cycle
   if ( psat(no) <= 0.0      ) cycle
   
   const0 = GS_b2(p) * psat(no) / max(0.01, co2atm-co2cmp(no)) / (1.0 + vpd/GS_b3(p))
   const4 = Unit_converter * dlen(doy) * psat(no)  !
   const5 = la(no) / (height(no)-bole(no))    !leaf area per crown layer (m2/step)
   do i = height(no), bole(no)+1, -1
      !Photosynthesis
      a1      = par_diffuse*par_diffuse_rel(i) + par_direct*par_direct_rel(no,i) !PAR income
      a2      = 1.0 - (1.0+a1*lue(no)/psat(no))**(-0.5)
      a4      = const4*const5*a2                             !GPP (g dm/layer/day)
      gpp_ind = gpp_ind + a4
      
      !Canopy conductance
      canopy_cond = canopy_cond + const5 * ( GS_b1(p) + const0 * a2 )
      
      !Production recorder for each foliage disk
      if (i-bole(no)<=10) npp_crownbottom(no, i-bole(no)) = & 
                          npp_crownbottom(no, i-bole(no)) + a4
      if (i==height(no) ) npp_crowntop(no)                = & 
                          npp_crowntop(no)                + a4
      if (i==bole(no)+1 ) npp_crownbottom_daily(no)       = & 
                          npp_crownbottom_daily(no)       + a4
   end do
   
   !Increase mass_available & c_uptake
   mass_available (no) = mass_available(no)  + gpp_ind
   gpp(p)              = gpp(p)              + gpp_ind
   npp(p)              = npp(p)              + gpp_ind
   gpp_daily_ind(no)   = gpp_daily_ind(no)   + gpp_ind
   
   mort_regu1(no)      = mort_regu1(no)      + gpp_ind
   flux_c_uptake_RR(1) = flux_c_uptake_RR(1) + gpp_ind
   
END DO
!$omp end do
!$omp end parallel

!_____________ Herbaceous PFTs
   !set PFT number, p
   if (pft_exist(C3g_no)) then
      p=C3g_no
   else
      p=C4g_no
   endif
   
   !obtain GPP per unit area
!$omp parallel private(i,j)
!$omp do private(k,const0,const1,const2,a2,a3,const3,x)
!DO i=1, DivedG !!!>>>>>>>>>>>>TN:rm
!DO j=1, DivedG !!!>>>>>>>>>>>>TN:rm
DO i=1, GRID%N_x !!!<<<<<<<<<<<<TN:add
DO j=1, GRID%N_y !!!<<<<<<<<<<<<TN:add
   
   k = int( MaxParClass * par_grass_rel(i,j) ) !PAR intensity class
   if ( k <= 0 )               cycle
   if ( psat_grass(k) <= 0.0 ) cycle
   
   const0 = GS_b2(p) * psat_grass(k) / max(0.01, co2atm-co2cmp_grass(k)) / (1.0 + vpd/GS_b3(p))
   const1 = Unit_converter * 2.0 * dlen(doy) * psat_grass(k) / eK(p)
   const2 = eK(p) * lue_grass(k) / psat_grass(k)
   
   !obtain daily gpp per unit area
   a2  = 1.0 + sqrt( 1.0 + (par * par_grass_rel(i,j)) * const2 )
   a3  = 1.0 + sqrt( 1.0 + (par * par_grass_rel(i,j)) * const2 * exp(-1.0*eK(p)*lai_grass(i,j)) )
   const3 = log(a2/a3)
   x      = const1 * const3          !gpp (g dm / m2   / day)
!   x      = x * ((Max_loc/DivedG)**2) !gpp (g dm / cell / day) !!!>>>>>>>>>>>>TN:rm
   x      = x * (real(GRID%Area)/real(GRID%N_tot)) !gpp (g dm / cell / day) !!!<<<<<<<<<<<<TN:add
   
   !increase gmass_available (g dm / cell / day)
   gmass_available(i,j) = gmass_available(i,j) + x
   
   !carbon balance
   gpp(p)              = gpp(p)              + x
   npp(p)              = npp(p)              + x
   flux_c_uptake_RR(1) = flux_c_uptake_RR(1) + x
   
   !canopy conductance
   canopy_cond = canopy_cond + &
!                 ( GS_b1(p)*lai_grass(i,j) + const0*(2.0/eK(p))*const3 ) * ((Max_loc/DivedG)**2) !!!>>>>>>>>>>>>TN:rm
                 ( GS_b1(p)*lai_grass(i,j) + const0*(2.0/eK(p))*const3 ) * (real(GRID%Area)/real(GRID%N_tot)) !!!<<<<<<<<<<<<TN:add
   
END DO
END DO
!$omp end do
!$omp end parallel

!Adjust unit of canopy conductance
!canopy_cond = canopy_cond / Max_loc / Max_loc !!!>>>>>>>>>>>>TN:rm
canopy_cond = canopy_cond / real(GRID%Area) !!!<<<<<<<<<<<<TN:add

END SUBROUTINE photosynthesis



!**************************************************************************************************
! Optimum Leaf Area Index: daily computation using running mean of PAR
!**************************************************************************************************
SUBROUTINE lai_optimum (tmp_air, tmp_soil)

!_____________ Set variables
!Namespace
   USE data_structure
   USE time_counter
   USE vegi_status_current1
   USE vegi_status_current2
   USE grid_status_current2
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none
   
!Argments
   real   ,intent(IN)::tmp_air   !2m air temperature (Celcius)
   real   ,intent(IN)::tmp_soil  !average soil temperature for 50cm  depth (Celcius)
   
!Local variables
!   real,dimension(DivedG,DivedG):: lai_opt_grass   !optimum LAI for grass (m2/m2) !!!>>>>>>>>>>>>TN:rm
   real,dimension(GRID%N_x,GRID%N_y):: lai_opt_grass   !optimum LAI for grass (m2/m2) !!!<<<<<<<<<<<<TN:add
   
   real    qt                    !QT10 of maintenance respiration
   real    tmp_sensibility_air   !temperature sensibility multiplier for foliage and stem
   real    tmp_sensibility_soil  !temperature sensibility multiplier for root
   real    root_per_foliage      !realized root per foliage fraction
   real    const1, const2        !temporal usage
   real    a0, a1, a2, a3        !for general usage
   integer i, j, k, p            !loop counter
   real,dimension(10)::lai_opt
   
!Reset
   lai_opt_grass(:,:) = 0.0
   
!p: PFT number of dominant grass type
   if (pft_exist(C3g_no)) then
      p=C3g_no
   else
      p=C4g_no
   endif
   
!Return when further calculation is unnecessary
!   if ( dlen(doy)  <= 0.1   ) return
!   if ( eK(p) <= 0.001 ) return
   
!Calculate temperature sensibility
   qt                   = 2.0 * exp( -0.009 * (tmp_air - 15.0) )
   tmp_sensibility_air  = exp( (tmp_air  - 15.0) * log(qt) / 10.0 )
   qt                   = 2.0 * exp( -0.009 * (tmp_soil - 15.0) )
   tmp_sensibility_soil = exp( (tmp_soil - 15.0) * log(qt) / 10.0 )
   
!Calculate realized root per foliage fraction
   root_per_foliage = max(1.0 , sum(gmass_root(:,:))) / &
                      max(1.0 , sum(gmass_leaf(:,:)))    
   root_per_foliage = min(root_per_foliage , 1.0/FR_ratio(p))
   
!Calculate daily maintenance cost
   !a0: maintenance respiration rate of leaf (gDM / gDM / day)
   a0 = RM(p) * PN_f(p) * tmp_sensibility_air
   
   !a1: turnover cost of leaf (g DM g DM-1 day-1)
   a1 = (TO_f(p)/Day_in_Year) * ( RG_f(p) - RG_f_suck(p) )
   
   !a2: maintenance respiration rate of root per unit leaf-mass (g DM g DM-1 day-1)
   !a2 = RM(p) * PN_r(p) * tmp_sensibility_soil / FR_ratio(p) !Previous method
   a2 = RM(p) * PN_r(p) * tmp_sensibility_soil * root_per_foliage
   
   !a3: turnover cost of root (g DM g DM-1 day-1)
   !a3 = 0.0
   !a3 = (TO_r(p)/Day_in_Year) * RG_r(p) / FR_ratio(p) !Previous method
   a3 = (TO_r(p)/Day_in_Year) * RG_r(p) * root_per_foliage
   
   !const1: daily maintenance cost of leaf (gDM / m2 / day)
   const1 = (a0+a1+a2+a3)/SLA(p)
   const1 = max(0.0, const1)
   
!optimum leaf index
lai_opt(:) = 0.0
   DO k=1, 10 !For each PAR intensity class
      if ( psat_grass(k) <= 0.001 ) cycle
      if ( lue_grass (k) <= 0.001 ) cycle
      if ( par           <= 1.00  ) cycle
      
      const2 = 1.0 - const1 / ( 0.09093 * dlen(doy) * psat_grass(k) )
      if (const2 .ne. 0.0) then
         const2 = 1.0 / (const2*const2) - 1.0
         const2 = (psat_grass(k)/lue_grass(k)) * const2
      endif
      
      if ( const2>0.0 ) then
         lai_opt(k) = ( log( (par*k/10.0) / const2) ) / eK(p)
         lai_opt(k) = max(0.0, lai_opt(k))
      endif
      
   END DO
   
!   DO i=1, DivedG !!!>>>>>>>>>>>>TN:rm
!   DO j=1, DivedG !!!>>>>>>>>>>>>TN:rm
!$omp parallel private(i,j)
!$omp do
   DO i=1, GRID%N_x !!!<<<<<<<<<<<<TN:add
   DO j=1, GRID%N_y !!!<<<<<<<<<<<<TN:add
      lai_opt_grass(i,j) = lai_opt( max(1, int( 10*par_grass_rel(i,j) )) )
      
      lai_opt_grass_RunningRecord(2:20,i,j) = lai_opt_grass_RunningRecord(1:19,i,j)
      lai_opt_grass_RunningRecord(1,i,j)    = lai_opt_grass(i,j)
   END DO
   END DO
!$omp end do
!$omp end parallel
   
END SUBROUTINE lai_optimum



!**************************************************************************************************
! Leaf phenology (daily computation)
!**************************************************************************************************
SUBROUTINE leaf_season (LAT, tmp_soil)

!_____________ Namespace
   USE data_structure
   USE time_counter
   USE vegi_status_current1
   USE vegi_status_current2
   USE grid_status_current1
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none
   
!_____________ Local parameters
!Minimum stock mass for dormance
   real,parameter::Stockmass_min_tree  = 10.0 !For woody PFTs (g dm/individual)
   real,parameter::Stockmass_min_grass = 50.0 !FOr grass PFTs (g dm/m2)
   
!Boundary conditions for leaf-pehnology (in days)
   integer,parameter::Days_dormance_min  = 60  !Minimum Dormance  Phase
   integer,parameter::Days_foliation_min = 60  !Minimum Foliation duration
   integer,parameter::Days_foliation_max = 300 !Maximum Foliation duration for Deciduous woody PFTs
   
   integer,parameter::Days_leaf_shed      = 14 !day length required for full leaf drop
   integer,parameter::Days_release_default= 30 !days required for full release of stock energy
   integer,parameter::Days_release_larch  = 14 !days required for full release of stock energy for larch

!_____________ Arguments
   real,intent(IN)::LAT       !latitude  (degree)
   real,dimension(NumSoil),intent(IN)::tmp_soil !Average soil temperature for each layer (C)
   
!_____________ Local variables
   real    tmp_air_10d_ave    !1- days running mean of air temperature (Celcius)
   real    mass_combust       !biomass combust     (g DM/day/stand)
   integer day_length_release !days required for full release of stocked energy
   
   logical                    flag_dry     !flag for available water
   logical                    flag_season  !flag for seasonal cycle
   logical,dimension(PFT_no)::flag         !flag for phenology change
   
   integer i, j, p, no     !for loop counter
   integer k, l, m, n      !for general usage
   real    x, y            !for general usage
   
   real pheno_gdd0 !GDD on 0 degree base (reset: 1 Jan for N-hemisphere, 1 July for S-hemisphere)
   real pheno_gdd5 !GDD on 5 degree base (reset: 1 Jan for N-hemisphere, 1 July for S-hemisphere)
   
   integer pheno_ncd !Number of chilling day, defined as number of days on which 
                     !the daily mean air temperature is below 5 Cecius degree
                     !(reset: 1 Nov for N-hemisphere, 1 May for S-hemisphere)
   
   integer day_until_bare  !day until full dormance
   
   integer doy_now
   
!_____________ Initialize local variables
   mass_combust = 0.000
   
!_____________ Prepare related variables
pheno_gdd0 = 0.0
pheno_gdd5 = 0.0
pheno_ncd  = 1

do i= 1, Day_in_Year !from 364 days ago to today
   
   !day of year for current repeat
   doy_now = doy - Day_in_Year + i
   if (doy_now <= 0) doy_now = doy_now + Day_in_Year
   
   !air temperature for current repeat
   x = tmp_air_RunningRecord (Day_in_Year-i+1)
   
   !pheno_gdd
   if ( (doy_now==1 .and. LAT>=0.0) .or. (doy_now==182 .and. LAT<0.0) ) then
   !reset on 1 Jan for N-hemisphere, 1 July for S-hemisphere
        pheno_gdd0 = 0.0 
        pheno_gdd5 = 0.0
   endif
   pheno_gdd0 = pheno_gdd0 + max(0.0, x    )
   pheno_gdd5 = pheno_gdd5 + max(0.0, x-5.0)
   
   !pheno_ncd
   if ( (doy_now==304 .and. LAT>=0.0) .or. (doy_now==121 .and. LAT<0.0) ) then
        pheno_ncd = 0 !reset on 1 Nov for N-hemisphere, 1 May for S-hemisphere
   endif
   
   if ( x < 5.0 ) then
        pheno_ncd = pheno_ncd + 1
   end if
   
enddo

!Prepare climatic factor 1
   tmp_air_10d_ave = sum( tmp_air_RunningRecord(1:10) ) / 10.0
   
!Prepare climatic factor 2
   m = 0 ; n = 0 ; l = 7
   do j = 0, 4
      x = sum( tmp_air_RunningRecord(1+j*l    :l+j*l    ) )
      y = sum( tmp_air_RunningRecord(1+(j+1)*l:l+(j+1)*l) )
      if ( x > y ) then
         m = m + 1
      else
         n = n + 1
      endif
   enddo
   
   if (m>n) then
      flag_season = .false.  !temperature increasing
   else
      flag_season = .true.   !temperature decreasing
   endif
   
!Prepare climatic factor 3
!dry_month: Number of dry month based on Priestley-Taylor model
   flag_dry = .false.
   x = 0.0
   y = 0.0
   do i=1, 20
      x = x + prec_RunningRecord  (i) !x: precipitation                 (mm/month)
      y = y + ev_pot_RunningRecord(i) !y: potential evapotranspirationn (mm/month)
   enddo
   if (x<y) flag_dry = .true.
   
!__________________ Checker (Foliage -> Domance) __________________
flag(:) = .false.
DO p = 1, PFT_no
!IF ( .not. pft_exist(p)   ) cycle
IF ( .not. phenology(p)   ) cycle
IF ( Phenology_type(p)==0 ) cycle
   
   x = tmp_air_10d_ave
   select case ( Phenology_type(p))
     
     case (1) !BoNS
       if (x < 7.0)                                          flag(p)=.true.
       
     case (2) !BoBS
       if ( sum(tmp_soil(1:5)) / 5.0 < 2.0 )                 flag(p)=.true.
       
     case (3) !TeBS
       if ( (x<9.0) .or. (x<10.0+tmp_coldest_20yr_ave) )     flag(p)=.true.
       
     case (4) !water-controlling PFTs 
       !if ( sum(stat_water_RunningRecord(1:10,p))/10.0<0.7 ) flag(p)=.true. !(Botta et al. 2000)
       
       !Modified for African vegetationz
       if ( sum(stat_water_RunningRecord(1:10,p))/10.0<0.0 ) flag(p)=.true.
       
     case (5:6) !grass
        flag(p)=.true.
        do k=1, 7
!          y = sum(lai_opt_grass_RunningRecord(k,:,:)) / DivedG / DivedG !!!>>>>>>>>>>>>TN:rm
          y = sum(lai_opt_grass_RunningRecord(k,:,:)) / real(GRID%N_tot) !!!<<<<<<<<<<<<TN:add
          if (y > 0.01)                                      flag(p)=.false.
        enddo
       
   end select
   
   !Adjustment using boundary conditions
   if ( dfl_leaf_onset(p) < Days_foliation_min ) then
       flag(p)=.false.
   endif 
   
   if ( dfl_leaf_onset(p) > Days_foliation_max ) then
   if ( Life_type(P) .ne. 3 )                    then
   if ( Life_type(P) .ne. 4 )                    then
      flag(p)=.true.
   endif
   endif
   endif
   
END DO

!__________________ Procedure (Foliage -> Dormance) __________________
DO p = 1, PFT_no
IF ( .not. flag(p) ) cycle
   !if (Logging==1) write (*,*)  'Foliage phase ---> Domance phase', p, doy
   phenology(p)     = .false.
   dfl_leaf_shed(p) = 0
END DO

!__________________ Gradual defoliation __________________
DO p = 1, PFT_no
day_until_bare = Days_leaf_shed - dfl_leaf_shed(p)
IF (day_until_bare < 1) cycle
   
   !for woody PFTs
   if ( (p.ne.C3g_no) .and. (p.ne.C4g_no) ) then
      do no=1, Max_no
      if (tree_exist(no) .and. pft(no)==p) then
          
          x                    = mass_leaf(no) / day_until_bare !leaf mass for shedding (g DM)
          mass_leaf     (no)   = mass_leaf(no) - x
          la            (no)   = mass_leaf(no) * SLA(P)
          mass_available(no)   = mass_available(no) + x * RG_f_suck(p)
          flux_litter_leaf     = flux_litter_leaf   + x * (1.0-RG_f_suck(p))
          
          !give minimum stock mass on the last leaf fall day
          if (day_until_bare==1 .and. mass_stock(no)<Stockmass_min_tree) then
             flux_litter_leaf  = flux_litter_leaf + mass_stock(no) - Stockmass_min_tree
             mass_stock(no)    = Stockmass_min_tree
          end if
          
      end if
      end do
   
   !for Grass PFTs
   else
      
!$omp parallel private(i,j)
!$omp do private(x)
!      do i=1, DivedG !!!>>>>>>>>>>>>TN:rm
!      do j=1, DivedG !!!>>>>>>>>>>>>TN:rm
      do i=1, GRID%N_x !!!<<<<<<<<<<<<TN:add
      do j=1, GRID%N_y !!!<<<<<<<<<<<<TN:add
         x                    = gmass_leaf(i,j) / day_until_bare !leaf mass for shedding (gDM/cell)
         gmass_available(i,j) = gmass_available(i,j) + x * RG_f_suck(p)
         flux_litter_ag       = flux_litter_ag   + x * (1.0- RG_f_suck(p))
         gmass_leaf(i,j)      = gmass_leaf(i,j)  - x
!         lai_grass(i,j)       = gmass_leaf(i,j) * SLA(p) / ((Max_loc/DivedG)**2) !!!>>>>>>>>>>>>TN:rm
         lai_grass(i,j)       = gmass_leaf(i,j) * SLA(p) / (real(GRID%Area)/real(GRID%N_tot)) !!!<<<<<<<<<<<<TN:add
         
         !give minimum stock mass
!         x = Stockmass_min_grass * ((Max_loc/DivedG)**2) !Minimum stock mass for grass layer (gDM/cell) !!!>>>>>>>>>>>>TN:rm
         x = Stockmass_min_grass * (real(GRID%Area)/real(GRID%N_tot)) !Minimum stock mass for grass layer (gDM/cell) !!!<<<<<<<<<<<<TN:add
         if ( day_until_bare==1 .and. gmass_stock(i,j)<x ) then
            
            if ( pool_litter_ag >= x-gmass_stock(i,j) ) then
               flux_litter_ag = flux_litter_ag - (x-gmass_stock(i,j))
            else
               flux_litter_ag = flux_litter_ag + pool_litter_ag
               mass_combust = mass_combust - (x-gmass_stock(i,j)-pool_litter_ag)
            endif
            gmass_stock(i,j) = x
            
         end if
      end do
      end do
!$omp end do
!$omp end parallel
   end if
   
END DO

!__________________ Checker (Dormance -> Foliage) __________________
flag(:) = .false.
DO p = 1, PFT_no
!IF (.not. pft_exist(p)   ) cycle
IF (phenology(p)         ) cycle
IF (Phenology_type(p)==0 ) cycle
   
   !Checker (Dormance -> Foliage)
   select case ( Phenology_type(p) )
     
     case (1) !boreal deciduous needle-leaved-woods [by Picard et al 2005]
       !x = 0.0
       !do j = 1, doy
       !   x = x + max(0.0, tmp_air_RunningRecord(j)-4.1)
       !end do
       !if ( x >= 65.0 )                                        flag(p) = .true. 
       
       ![by Yamazaki et al. 2007, Hydrological Processes: 21(15)]
       if ( pheno_gdd0>100.0 .and. tmp_soil(1)>5.0 ) flag(p) = .true. 
       
!      if ( LAT>=0.0 .and. (doy<100 .or. doy>190) )            flag(p) = .false.
!      if ( LAT< 0.0 .and. (doy<212-LAT .or. doy>312-LAT) )    flag(p) = .false.
       
     case (2:3) !other temperature-controlling PFTs [Botta et al 2000 procedure]
       if ( pheno_gdd5 > -68.0 + 638.0 * exp(-0.01*pheno_ncd) ) flag(p) = .true. 
       if ( LAT>=0.0 .and. (doy< 30+LAT .or. doy>130+LAT) )     flag(p) = .false.
       if ( LAT< 0.0 .and. (doy<212-LAT .or. doy>312-LAT) )     flag(p) = .false.
       
     case (4) !water-controlling deciduous
       if ( sum(stat_water_RunningRecord(1:10,p)) / 10.0 > 0.0 )    flag(p) = .true.
       
     case (5:6) !grass
       flag(p)=.true.
       do k=1, 7
!          y = sum(lai_opt_grass_RunningRecord(k,:,:)) / DivedG / DivedG !!!>>>>>>>>>>>>TN:rm
          y = sum(lai_opt_grass_RunningRecord(k,:,:)) / real(GRID%N_tot) !!!<<<<<<<<<<<<TN:add
          if (y <= 0.01)                                    flag(p) = .false.
       enddo
       
       !if ( LAT>=0.0 .and. (doy< 30+LAT .or. doy>130+LAT) ) flag(p) = .false.
       !if ( LAT< 0.0 .and. (doy<212-LAT .or. doy>312-LAT) ) flag(p) = .false.
       
   end select
   if ( dfl_leaf_shed(p) < Days_dormance_min ) flag(p)=.false.
   
END DO

!__________________ Procedure (Dormance -> Foliage) __________________
DO p = 1, PFT_no
IF (.not. flag(p)) cycle
   !if (Logging==1) write (*,*) 'Domance phase ---> Foliage phase', p, doy
   phenology(p) = .true.
   dfl_leaf_onset(p) = 0
END DO

!__________________ Gradual release of stock energy __________________
DO p = 1, PFT_no
!IF (.not. pft_exist(p)  ) cycle
IF (.not. phenology(p)  ) cycle
IF (Phenology_type(p)==0) cycle
   
   if ( Life_type(p).ne.2 ) then
      day_length_release = Days_release_default
   else
      day_length_release = Days_release_larch
   endif
   
IF (dfl_leaf_onset(p) >= day_length_release) cycle
   
   !For woody PFTs
   If ( (p.ne.C3g_no) .and. (p.ne.C4g_no) ) then
       do no=1, Max_no
       if (tree_exist(no) .and. pft(no)==p) then
          !x: stock mass for releasing in the day (g)
          x = mass_stock(no)/(day_length_release - dfl_leaf_onset(p))
          
          !release process
          mass_available (no) = mass_available(no) + x / RG_stock_out
          mass_stock     (no) = mass_stock    (no) - x
          resp_trunk     (no) = resp_trunk    (no) + (x - x / RG_stock_out)
          mass_combust        = mass_combust       + (x - x / RG_stock_out)
          mort_regu1     (no) = mort_regu1    (no) - (x - x / RG_stock_out)
          npp            (p)  = npp(p)             - (x - x / RG_stock_out)
       endif
       enddo
   
   !For Grass species
   Else
!$omp parallel private(i,j)
!$omp do private(x)
!      do i=1, DivedG !!!>>>>>>>>>>>>TN:rm
!      do j=1, DivedG !!!>>>>>>>>>>>>TN:rm
      do i=1, GRID%N_x !!!<<<<<<<<<<<<TN:add
      do j=1, GRID%N_y !!!<<<<<<<<<<<<TN:add
         !x: stock mass consumed within a day (g/cell)
         x = gmass_stock(i,j) / (day_length_release - dfl_leaf_onset(p))
         
         !release process
         gmass_available(i,j) = gmass_available(i,j) + x / RG_stock_out
         gmass_stock(i,j)     = gmass_stock(i,j)     - x
         resp_grass_bg        = resp_grass_bg   + (x - x / RG_stock_out)
         mass_combust         = mass_combust    + (x - x / RG_stock_out)
         npp (p)              = npp(p)          - (x - x / RG_stock_out)
      end do
      end do
!$omp end do
!$omp end parallel
   Endif
   
END DO

!_____________ Total carbon emission by maintenance respiration
   flux_c_gro_RR(1) = flux_c_gro_RR(1) + mass_combust
   
END SUBROUTINE leaf_season



!**************************************************************************************************
! Maintenance respiration (daily computation)
!**************************************************************************************************
SUBROUTINE maintenance_resp (tmp_air, tmp_soil)

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
   
!Set Local Parameter
   !daily proportion biomass remove when resource for maintenance respiration resource absent
   real,parameter::Frac_organ_remove=0.01
   
!Define arguments
   real,intent(IN)::tmp_air  !air  temperature (Cecius)
   real,intent(IN)::tmp_soil !average soil temperature of 50cm depth (Cecius)
   
!Define local variables
   real qt                    !QT
   real tmp_sensibility_air   !temperature sensibility multiplier for foliage and stem
   real tmp_sensibility_soil  !temperature sensibility multiplier for root
   real mass_required         !resource required for maintenance respiration (g)
   real mass_resp             !resource usage for mainatenance respiration (g)
   real f_sapwood             !sapwood_mass / stem_mass
   
   real mass_combust          !biomass combust (g DM/day/stand)
   
   real mass_required_leaf  !carbon required for respiration for leaf               (gDM)
   real mass_required_trunc !carbon required for respiration for truc               (gDM)
   real mass_required_root  !carbon required for respiration for root               (gDM)
   real mass_required_ag    !carbon required for respiration for above ground grass (gDM)
   real mass_required_bg    !carbon required for respiration for below ground grass (gDM)
   
   integer no, i, j, p     !loop counter
   real    a1, a2, a3      !for general usage
   real    x, y            !for general usage
   
!_____________ Set initial values and prepare related variables
!dry-mass to combust (g dm / stand)
   mass_combust = 0.000
   
!temperature sensibility, for foliage and sapwood
   !for aboveground organ
   qt                   = 2.0 * exp( -0.009 * (tmp_air-15.0)  )
   tmp_sensibility_air  = exp( (tmp_air  - 15.0) * log(qt) / 10.0 )
   !for belowground organ
   qt                   = 2.0 * exp( -0.009 * (tmp_soil-15.0) )
   tmp_sensibility_soil = exp( (tmp_soil - 15.0) * log(qt) / 10.0 )
   
!_____________ Woody PFTs
!$omp parallel private(no)
!$omp do private(p,f_sapwood,mass_required_trunc,mass_required_root,mass_required,mass_resp,x,y,a1,a2,i)
DO no=1, Max_no
   !set PFT number
   p = pft(no)
   
   If ( .not. tree_exist(no) ) cycle
   
   !Calculate daily resource required for maintenance respiration
   f_sapwood           = 1.0-(dbh_heartwood(no)**2)/( (dbh_sapwood(no)+dbh_heartwood(no))**2 )
   mass_required_leaf  = RM(p)* PN_f(p)* mass_leaf (no)           * tmp_sensibility_air
   mass_required_trunc = RM(p)* PN_s(p)* mass_trunk(no)* f_sapwood* tmp_sensibility_air
   mass_required_root  = RM(p)* PN_r(p)* mass_root (no)           * tmp_sensibility_soil
   mass_required       = mass_required_leaf + mass_required_trunc + mass_required_root
   
   !phase 1: use mass_available
   if (mass_required > 0.001) then
      if (mass_required < mass_available(no)) then
         resp_leaf (no)    = resp_leaf (no)     + mass_required_leaf  
         resp_trunk(no)    = resp_trunk(no)     + mass_required_trunc 
         resp_root (no)    = resp_root (no)     + mass_required_root  
         mort_regu1(no)    = mort_regu1(no)     - mass_required
         npp(p)            = npp(p)             - mass_required
         mass_combust      = mass_combust       + mass_required
         
         mass_available(no)= mass_available(no) - mass_required
         mass_required     = 0.0
      else
         resp_leaf (no)    = resp_leaf (no)     + mass_available(no) * (mass_required_leaf  / mass_required)
         resp_trunk(no)    = resp_trunk(no)     + mass_available(no) * (mass_required_trunc / mass_required)
         resp_root (no)    = resp_root (no)     + mass_available(no) * (mass_required_root  / mass_required)
         mort_regu1(no)    = mort_regu1(no)     - mass_available(no)
         npp(p)            = npp(p)             - mass_available(no)
         mass_combust      = mass_combust       + mass_available(no)
         
         mass_required     = mass_required      - mass_available(no)
         mass_available(no)= 0.0
      endif
   endif
   
   !phase 2: use mass_stock
   if (mass_required > 0.001) then
      mass_resp      = min(mass_required, mass_stock(no)/ RG_stock_out)
      
      mass_stock(no) = mass_stock(no) - mass_resp * RG_stock_out
      
      resp_trunk(no) = resp_trunk(no) + mass_resp * RG_stock_out
      mort_regu1(no) = mort_regu1(no) - mass_resp * RG_stock_out
      npp(p)         = npp(p)         - mass_resp * RG_stock_out
      mass_combust   = mass_combust   + mass_resp * RG_stock_out
      mass_required  = mass_required  - mass_resp
   end if
   
   !phase 3: when resource absent, a portion of leaves drop
   if (mass_required > 0.001) then
      a1 = Frac_organ_remove * mass_leaf(no)   !leaf mass that will purged (gDM)
      a2 = Frac_organ_remove * mass_root(no)   !root mass that will purged (gDM)
      a3 = dbh_sapwood(no) * Frac_organ_remove !length of sapwood to heartwood (m)
      
      mass_leaf    (no) = mass_leaf    (no) - a1
      mass_root    (no) = mass_root    (no) - a2
      dbh_sapwood  (no) = dbh_sapwood  (no) - a3
      dbh_heartwood(no) = dbh_heartwood(no) + a3
      la           (no) = mass_leaf(no) * SLA(p)
      
      flux_litter_leaf = flux_litter_leaf + a1
      flux_litter_root = flux_litter_leaf + a2
   end if
   
   !record crown-bottom-layer stat_leaf
   x  = mass_leaf(no) / ( height(no)-bole(no) ) !leaf mass per crown disk(gDM/step)
   y  = x / FR_ratio(p)                         !required fine root for each crown disk(gDM/step)
   
   !a1: leaf maintenance cost (g/day/step)
   a1 = x * ( PN_f(p)*RM(p)*tmp_sensibility_air  + (TO_f(p)/Day_in_Year)*(RG_f(p)-RG_f_suck(p)) )
   
   !a2: root maintenance cost (g/day/step)
   a2 = y * ( PN_r(p)*RM(p)*tmp_sensibility_soil + (TO_r(p)/Day_in_Year)*RG_r(p)                )
   
   do i = 1, 10
   npp_crownbottom(no,i)     = npp_crownbottom(no,i)      - a1 - a2
   end do
   
   npp_crowntop(no)          = npp_crowntop(no)           - a1 - a2
   npp_crownbottom_daily(no) = npp_crownbottom_daily (no) - a1 - a2
   
End Do
!$omp end do
!$omp end parallel
   
!_____________ Herbaceous PFTs
   !set pft number, p
   if (pft_exist(C3g_no)) then
      p=C3g_no
   else
      p=C4g_no
   endif
   
if (phenology(p)) then
!$omp parallel private(i,j)
!$omp do private(mass_resp,mass_required_ag,mass_required_bg,mass_required,mass_resp)
!do i=1, DivedG !!!>>>>>>>>>>>>TN:rm
!do j=1, DivedG !!!>>>>>>>>>>>>TN:rm
do i=1, GRID%N_x !!!<<<<<<<<<<<<TN:add
do j=1, GRID%N_y !!!<<<<<<<<<<<<TN:add
   
   !reset mass respiration 
   mass_resp  = 0.000
   
   !calculate daily mass_required for maintenance respiration (g/cell)
   mass_required_ag = RM(p) * PN_f(p) * tmp_sensibility_air  * gmass_leaf(i,j)
   mass_required_bg = RM(p) * PN_r(p) * tmp_sensibility_soil * gmass_root(i,j)
   mass_required    = mass_required_ag + mass_required_bg
   
   !Source 1: mass_available
   if (mass_required > 0.001) then
      if (mass_required < gmass_available(i,j)) then
         resp_grass_ag  = resp_grass_ag + mass_required * (mass_required_ag / mass_required)
         resp_grass_bg  = resp_grass_bg + mass_required * (mass_required_bg / mass_required)
         mass_combust   = mass_combust  + mass_required
         npp(p)         = npp(p)        - mass_required
         
         gmass_available(i,j) = gmass_available(i,j) - mass_required
         mass_required        = 0.0
      else
         resp_grass_ag  = resp_grass_ag + gmass_available(i,j) * (mass_required_ag / mass_required)
         resp_grass_bg  = resp_grass_bg + gmass_available(i,j) * (mass_required_bg / mass_required)
         mass_combust   = mass_combust  + gmass_available(i,j)
         npp(p)         = npp(p)        - gmass_available(i,j)
         
         mass_required        = mass_required - gmass_available(i,j)
         gmass_available(i,j) = 0.0
      endif
   endif
   
   !Source 2: mass_stock
   if (mass_required > 0.001) then
      mass_resp        = min(mass_required, gmass_stock(i,j)/RG_stock_out)
      
      gmass_stock(i,j) = gmass_stock(i,j) - mass_resp * RG_stock_out
      
      resp_grass_bg    = resp_grass_bg    + mass_resp * RG_stock_out
      mass_combust     = mass_combust     + mass_resp * RG_stock_out
      npp(p)           = npp(p)           - mass_resp * RG_stock_out
      mass_required    = mass_required    - mass_resp
   endif
   
   !Source 3: when resource absent
   if (mass_required > 0.001) then
      a1 = Frac_organ_remove * gmass_leaf(i,j)
      a2 = Frac_organ_remove * gmass_root(i,j)
      
      gmass_leaf(i,j)  = gmass_leaf(i,j)  - a1
      gmass_root(i,j)  = gmass_root(i,j)  - a2
      
      flux_litter_ag   = flux_litter_ag   + a1
      flux_litter_bg   = flux_litter_bg   + a2
      
!      lai_grass(i,j)   = gmass_leaf(i,j) * SLA(p) / ((Max_loc/DivedG)**2) !!!>>>>>>>>>>>>TN:rm
      lai_grass(i,j)   = gmass_leaf(i,j) * SLA(p) / (real(GRID%Area)/real(GRID%N_tot)) !!!<<<<<<<<<<<<TN:add
   endif
   
end do
end do
!$omp end do
!$omp end parallel
end if

!_____________ Total carbon emission by maintenance respiration
   flux_c_mnt_RR(1) = flux_c_mnt_RR(1) + mass_combust
   
END SUBROUTINE maintenance_resp



!**************************************************************************************************
! Turnover (daily computation)
!**************************************************************************************************
SUBROUTINE turnover ()
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
   
!Local variables
   real    t_rate, t_rate_ag, t_rate_bg !daily turnover rates
   real    t_mass      !daily turnover mass
   integer p           !for PFT number
   integer no, i, j    !for loop counter
   real    x, y, z     !for general usage
   
!_____________ Woody PFTs
!$omp parallel private(no)
!$omp do private(p,t_rate,t_mass,x,y,z)
DO no=1, Max_no 
if (tree_exist(no)) then
   !set pft number
   p = pft(no)
   
   !foliage
   t_rate = TO_f(p) / Day_in_Year
   
   t_mass             = mass_leaf(no) * t_rate !leaf mass to be drop (g DM)
   mass_leaf(no)      = mass_leaf(no)      - t_mass 
   mass_available(no) = mass_available(no) + t_mass * RG_f_suck(p)
   la(no)             = mass_leaf(no) * SLA(p)
   
   flux_litter_leaf = flux_litter_leaf + t_mass * (1.0 - RG_f_suck(p))
   
   !root
   t_rate = TO_r(p) / Day_in_Year
   
   t_mass           = mass_root(no) *    t_rate !root mass to be drop (g DM)
   mass_root(no)    = mass_root(no)    - t_mass
   
   flux_litter_root = flux_litter_root + t_mass
   
   !sapwood
  if   (Life_type(p)==1 .or. Life_type(p)==5) then
      !Tropical rain trees
      x = ALM5(p)       * (dbh_sapwood(no) + dbh_heartwood(no)) !new sapwood diameter   (m)
      y = (1.0-ALM5(p)) * (dbh_sapwood(no) + dbh_heartwood(no)) !new heartwood diameter (m)
      
      dbh_sapwood(no)   = x
      dbh_heartwood(no) = y
      
   elseif (Life_type(p)==2 .or. Life_type(p)==6) then
      !Boreal deciduous trees
      x                 = dbh_sapwood(no) + dbh_heartwood(no) !DBH(m)
      dbh_sapwood(no)   = min(x, 0.0188)
      dbh_heartwood(no) = x - dbh_sapwood(no)
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add Mangrove feature ***tentative
  else if (Life_type(p)==7 .or. Life_type(p)==8) then
      !Tropical rain trees
      x = ALM5(p)       * (dbh_sapwood(no) + dbh_heartwood(no)) !new sapwood diameter   (m)
      y = (1.0-ALM5(p)) * (dbh_sapwood(no) + dbh_heartwood(no)) !new heartwood diameter (m)
      
      dbh_sapwood(no)   = x
      dbh_heartwood(no) = y
      
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   else
      x                 = dbh_sapwood(no)
      y                 = dbh_heartwood(no)
      z                 = TO_s(pft(no)) / Day_in_Year
      t_rate            = sqrt( y*y + z*( x*x + 2*x*y ) ) - y
      dbh_sapwood(no)   = dbh_sapwood(no)   - t_rate
      dbh_heartwood(no) = dbh_heartwood(no) + t_rate
   endif
   
end if
END DO
!$omp end do
!$omp end parallel

!_____________ Herbaceous PFTs
   !set pft number, p
   if (pft_exist(C3g_no)) then
      p=C3g_no
   else
      p=C4g_no
   endif
   
   !turnover rates
   t_rate_ag = TO_f(p)/Day_in_Year !above ground turnover rate
   t_rate_bg = TO_r(p)/Day_in_Year !below ground turnover rate
   
!$omp parallel private(i,j)
!$omp do private(t_mass)
!DO i=1, DivedG !!!>>>>>>>>>>>>TN:rm
!DO j=1, DivedG !!!>>>>>>>>>>>>TN:rm
DO i=1, GRID%N_x !!!<<<<<<<<<<<<TN:add
DO j=1, GRID%N_y !!!<<<<<<<<<<<<TN:add
   
   !above ground turnover
   t_mass = gmass_leaf(i,j) * t_rate_ag  !leaf mass to drop (g/cell)
   
   gmass_leaf     (i,j) = gmass_leaf     (i,j) - t_mass
   gmass_available(i,j) = gmass_available(i,j) + t_mass * RG_f_suck(p)
!   lai_grass      (i,j) = gmass_leaf     (i,j) * SLA(p) / ((Max_loc/DivedG)**2) !!!>>>>>>>>>>>>TN:rm
   lai_grass      (i,j) = gmass_leaf     (i,j) * SLA(p) / (real(GRID%Area)/real(GRID%N_tot)) !!!<<<<<<<<<<<<TN:add
   
   flux_litter_ag   = flux_litter_ag   + t_mass * (1.0-RG_f_suck(p))
   
   !below ground turnover
   t_mass          = gmass_root(i,j) * t_rate_bg  !root mass to drop (g/m2)
   
   gmass_root(i,j) = gmass_root(i,j) - t_mass
   flux_litter_bg  = flux_litter_bg  + t_mass
   
END DO
END DO

END SUBROUTINE turnover



!**************************************************************************************************
! growth woody PFTs (daily computation)
!**************************************************************************************************
SUBROUTINE growth_wood ()
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
   
!Local parameters for tropical rain forest from Huth & Ditzer (2000)
   ! bottom height of forest layer no 2,3,4,5 (m) 
   !(top    height of forest layer no 1,2,3,4 (m) )
   real,dimension(4),parameter::Layer_top = (/1.3, 15.0, 25.0, 36.0/)
   
   !maximum LAI for each layer (m2)
   real,parameter::LAI_max_cohort = 2.0
   
   !maximum LA that can increase within a day (fraction to crown area)
   real,parameter::LAgrow_max_day = 0.1
   
!_____________ Define local variables
!   real   ,dimension((int(0.99999+Max_loc/20.0))**2, 5)::cohort_lai     !!!>>>>>>>>>>>>TN:rm
!   logical,dimension((int(0.99999+Max_loc/20.0))**2, 5)::cohort_crowded !!!>>>>>>>>>>>>TN:rm
   real   ,dimension(int(0.99999+GRID%Max_x/20.0)*int(0.99999+GRID%Max_y/20.0), 5)::cohort_lai      !!!<<<<<<<<<<<<TN:add これで良いかわからない。要確認
   logical,dimension(int(0.99999+GRID%Max_x/20.0)*int(0.99999+GRID%Max_y/20.0), 5)::cohort_crowded  !!!<<<<<<<<<<<<TN:add これで良いかわからない。要確認
   integer,dimension(Max_no)::id_location
   integer,dimension(Max_no)::id_layer
   
   integer p                      !for PFT number
   integer i, j, no               !for loop counter
   real    x, a1, a2, a3, a4, a6  !for temporal usage
   real    mass_combust           !biomass combust (g DM/day/stand)
   
!_____________ Initialize
   mass_combust = 0.000
   
!_____________ Determine cohort structure
   cohort_lai      (:,:) = 0.0     !(location,layer) crown area for each cohort (m2)
   cohort_crowded  (:,:) = .false. !(location,layer) crowded flag for each cohort
   id_location     (:)   = 0       !(tree_number)    location number that each tree belongs
   id_layer        (:)   = 0       !(tree_number)    layer number that each tree belongs
   
!$omp parallel private(no)
!$omp do private(x)
   Do no=1, Max_no 
   If ( tree_exist(no) ) then
      
!      id_location(no) = 1+int(bole_x(no)/20) + int(0.99999+Max_loc/20.0) * int(bole_y(no)/20)  !!!>>>>>>>>>>>>TN:rm
      id_location(no) = 1+int(bole_x(no)/20) + int(0.99999+GRID%Max_y/20.0) * int(bole_y(no)/20)  !!!<<<<<<<<<<<<TN:add これで良いかわからない。要確認
      
      x = height(no)*STEP+1.3 !x: tree height (m)
      if     (x<Layer_top(1)) then ;id_layer(no) = 1
      elseif (x<Layer_top(2)) then ;id_layer(no) = 2
      elseif (x<Layer_top(3)) then ;id_layer(no) = 3
      elseif (x<Layer_top(4)) then ;id_layer(no) = 4
      else                         ;id_layer(no) = 5
      endif
      
      cohort_lai(id_location(no),id_layer(no)) = &
      cohort_lai(id_location(no),id_layer(no)) + la(no)*SLA(pft(no))
      
   End if
   End do
!$omp end do
!$omp end parallel
   
!   Do i=1, (int(0.99999+Max_loc/20.0))**2 !for each location !!!>>>>>>>>>>>>TN:rm
   Do i=1, int(0.99999+GRID%Max_x/20.0)*int(0.99999+GRID%Max_y/20.0) !for each location !!!<<<<<<<<<<<<TN:add これで良いかわからない。要確認
   Do j=1, 5                              !for each layer
      if ( cohort_lai(i,j)/400 > LAI_max_cohort ) cohort_crowded(i,j)=.true.
   End do
   End do
   
!_____________ Each individual of woody PFT procedure
!$omp parallel private(no)
!$omp do private(p,a1,a2,a3,a4,a6)
DO no=1, Max_no 
IF ( .not. tree_exist(no)     ) cycle
IF ( .not. phenology(pft(no)) ) cycle
    
    !SET PFT NUMBER
    p = pft(no)
    
    !ROOT MASS SUPPLEMENT
    !a1: required root mass (g)
    a1 = max( 0.00000, mass_leaf(no)/FR_ratio(p) - mass_root(no) )
    !a2: maximum root mass increament allowded by available resource (g)
    a2 = mass_available(no) / RG_r(p)
    !a3 = actual increament of Root mass with considering resource_available
    a3 = min( a1, a2 )
    !root mass supplement precedure
    mass_root(no)      = mass_root(no)      +  a3                
    mass_available(no) = mass_available(no) -  a3 * RG_r(p)      
    
    resp_root (no)     = resp_root (no)     + (a3 * RG_r(p) - a3)
    mort_regu1(no)     = mort_regu1(no)     - (a3 * RG_r(p) - a3)
    mass_combust       = mass_combust       + (a3 * RG_r(p) - a3)
    npp(p)             = npp(p)             - (a3 * RG_r(p) - a3)
    
    !STOCK RESOURCE SUPPLEMENT
    if (dfl_leaf_onset(p) > 30) then !after more than 30 days from foliation
       !a1 = resource requirement for stock (g)
       a1 = max( 0.00000, mass_leaf(no) - mass_stock(no) )
       !a2 = maximum increment of mass_stock (g)
       a2 = mass_available(no) / RG_stock_in
       !a3:actual increment of mass_stock (g/m2)
       a3 = min( a1, a2 )
       
       !stock resource supplement precedure
       mass_stock    (no) = mass_stock    (no) +  a3                    
       mass_available(no) = mass_available(no) -  a3 * RG_stock_in      
       
       resp_root     (no) = resp_root      (no)+ (a3 * RG_stock_in - a3)
       mass_combust       = mass_combust       + (a3 * RG_stock_in - a3)
       npp(p)             = npp(p)             - (a3 * RG_stock_in - a3)
       mort_regu1(no)     = mort_regu1(no)     - (a3 * RG_stock_in - a3)
    end if
    
    !FOLIAGE MASS SUPPLEMENT
    !a1: deficit of foliage mass until fullfill existing crown size (g)
    if ( Life_type(p)==1 .or. Life_type(p)==5 ) then
    !For TrBE
       a1 = crown_area(no)
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add Mangrove feature ***tentative
    else if ( Life_type(p)==7 .or. Life_type(p)==8 ) then
    !For Mangrove
       a1 = crown_area(no)
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
    else
    !For other woody PFTs
       a1 = crown_area(no) + PI * crown_diameter(no) * real(height(no)-bole(no)) * STEP
    endif
    a1 = a1 * LA_max(p) / SLA(p)
    a1 = max(0.0, a1-mass_leaf(no) )
    
    !a2: deficit of foliage mass until fullfill existing Sapwood (g)
    if     ( Life_type(p)==1 .or. Life_type(p)==5 ) then
    !For TrBE and TrBS
       a2 = 100000.0 !This large value means that no constrain from this factor
    elseif ( Life_type(p)==2 ) then
    !For BoNS
!      a2 = 14.8 * ( (100*dbh_sapwood(no)+100*dbh_heartwood(no))**1.68 ) !Kajimoto et al (2006) 
       a2 = 330 + 50300 * ( (dbh_sapwood(no)+dbh_heartwood(no))**2 )      !Schulze et al (1995)
       a2 = max(0.0, a2-mass_leaf(no))
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add Mangrove feature ***tentative
    else if ( Life_type(p)==7 .or. Life_type(p)==8 ) then
    !For TrBE and TrBS
       a2 = 100000.0 !This large value means that no constrain from this factor
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
    else 
    !For other woody PFTs
       a2 = ALM1(p) * (1/SLA(p)) * (PI/4) * (dbh_sapwood(no)**2 + 2*dbh_sapwood(no)*dbh_heartwood(no))
       a2 = max(0.0, a2-mass_leaf(no) )
    endif
    
    !a3: maximum foliage mass with current resource available (g)
    a3 = mass_available(no) / RG_f(p)
    
    !a4: maximum foliage mass that can increase within a day (g)
    a4 = LAgrow_max_day * crown_area(no) / SLA(p)
    
    !a6: acutual increament of foliage mass (g)
    a6 = min(a1, a2, a3, a4)
    
    if ( Life_type(p)==1 .or. Life_type(p)==5 ) then
    !For TrBE (Tropical broadleaved evergreeen trees)
      !if ( cohort_crowded(id_location(no),id_layer(no)) ) a6=0.0
       if ( npp_crownbottom_daily(no) < 0.0        ) a6=0.0
    elseif ( Life_type(p)==2 ) then
    !For BoNS (larch)
       if ( npp_crownbottom_daily(no) < 0.0        ) a6=0.0
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add Mangrove feature ***tentative
    elseif ( Life_type(p)==7 .or. Life_type(p)==8 ) then
    !For TrBE (Tropical broadleaved evergreeen trees)
      !if ( cohort_crowded(id_location(no),id_layer(no)) ) a6=0.0
       if ( npp_crownbottom_daily(no) < 0.0        ) a6=0.0
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
    endif
    
    !reset daily leaf area stat on bottom of crown layer
    npp_crownbottom_daily(no) = 0.0
    
    !a6: actual increament (g)
    if ( a6 > 0.0) then
        mass_leaf(no)      = mass_leaf(no)       +  a6                
        mass_available(no) = mass_available(no)  -  a6 * RG_f(p)      
        la(no)             = mass_leaf(no) * SLA(p)
        
        resp_leaf(no)      = resp_leaf(no)       + (a6 * RG_f(p) - a6)
        mass_combust       = mass_combust        + (a6 * RG_f(p) - a6)
        npp(p)             = npp(p)              - (a6 * RG_f(p) - a6)
        mort_regu1(no)     = mort_regu1(no)      - (a6 * RG_f(p) - a6)
    end if
    
    !update mort_regu2 (annumal mean of leaf area)
    mort_regu2(no) = mort_regu2(no) + la(no)/Day_in_Year
    
END DO
!$omp end do
!$omp end parallel
   
!_____________ carbon emission
   flux_c_gro_RR(1) = flux_c_gro_RR(1) + mass_combust
   
End subroutine growth_wood



!**************************************************************************************************
! growth_grass (daily computation)
!**************************************************************************************************
SUBROUTINE growth_grass ()
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
   !Grass seeds biomass that always present (g/m2), This value should be very small
   real,parameter::Grass_seeds = 0.0!!!>>>>>>>>>>>>TN:changed 0.1 -> 0.0 for no grass establishment
   
   !Fraction of stock resource to aboveground biomass
!  real,parameter::Frac_Stock = 0.45
!  real,parameter::Frac_Stock = 0.71
   real,parameter::Frac_Stock = 1.3
   
   !Minimum stock resource for reproduction (gDM m-2)
   real   ,parameter::Repro_Stock_min = 100.0 
   
   !No reproduction days from onset date (day)
   integer,parameter::NoRepro_Phase  = 30    
   
   !No saving days from onset date (day)
   integer,parameter::NoStock_Phase  = 30    
   
   !Period for averaging LAI optimum when regulating leaf growth (day)
   integer,parameter::LaiOpt_period_mean = 7
   
!Local variables
   real    Unit_conv     !Unit_conv: m2 -> grass_cell
   real    lai_opt       !optimum LAI of grass layer (m2/m2)
   real    mass_combust  !biomass combust (g dm/day/stand)
   real    x, a1, a2, a3 !for temporal usage
   integer p             !PFT number
   integer i, j          !counter
   
!_____________ Initialize local variables
   !initialize combustion biomass
   mass_combust = 0.000 !
   
   !p: PFT number of the dominant grass type
                          p=C4g_no
   if (pft_exist(C3g_no)) p=C3g_no
   
   !Unit_conv: m2 -> grass_cell
!   Unit_conv = (Max_loc/DivedG)**2 !!!>>>>>>>>>>>>TN:rm
   Unit_conv = real(GRID%Area)/real(GRID%N_tot) !!!<<<<<<<<<<<<TN:add
   
!_____________ each grass type procedure
IF ( phenology(p) ) then
!$omp parallel private(i,j)
!$omp do private(lai_opt,a1,a2,a3)
!DO i=1, DivedG !!!>>>>>>>>>>>>TN:rm
!DO j=1, DivedG !!!>>>>>>>>>>>>TN:rm
DO i=1, GRID%N_x !!!<<<<<<<<<<<<TN:add
DO j=1, GRID%N_y !!!<<<<<<<<<<<<TN:add
    
    !set lai_opt, running mean of optimum leaf area index
    lai_opt = sum(lai_opt_grass_RunningRecord(1:LaiOpt_period_mean,i,j)) / real(LaiOpt_period_mean)
    
    !excess foliage senescence procedure (optional)
    ! if (lai_grass(i,j) > lai_opt) then
    !    !a1: excess gmass_leaf (g drymatter/cell)
    !    a1             = (lai_grass(i,j)-lai_opt) / SLA(p) * Unit_conv
    !    !a2: gmass_leaf that senescence (g drymatter/cell)
    !    a2               = a1 * 0.020
    !    gmass_leaf(i,j)  = gmass_leaf(i,j)  - a2
    !    flux_litter_ag   = flux_litter_ag   + a2
    !    lai_grass(i,j)   = gmass_leaf(i,j) * SLA(p) / Unit_conv
    ! end if
    
    !root mass supplement
       !a1: required root mass (g dm / cell)
       !a1 = max(0.0, gmass_leaf(i,j)/FR_ratio(p) - gmass_root(i,j)) !Previous method
       a1 = max(0.0, gmass_leaf(i,j)/FR_ratio(p) - gmass_root(i,j) - gmass_stock(i,j))
       
       !a2: maximum root mass increament (g dm / cell)
       a2 = gmass_available(i,j) / RG_r(p)
       
       !a3: actual increamet of root mass (g dm / grass_cell)
       a3 = min(a1, a2)
       
       !root mass supplement procedure
       gmass_available(i,j) = gmass_available(i,j) - a3 * RG_r(p)
       gmass_root     (i,j) = gmass_root     (i,j) + a3          
       
       resp_grass_bg   = resp_grass_bg   + a3 * (RG_r(p)-1.0)
       mass_combust    = mass_combust    + a3 * (RG_r(p)-1.0)
       npp(p)          = npp(p)          - a3 * (RG_r(p)-1.0)
    
    !stock resource supplement
    if ( dfl_leaf_onset(p) > NoStock_Phase ) then !after more than 30 days from foliation
       !a1:resource requirement for stock (g DM/cell)
!       a1 = max( 0.0, gmass_leaf(i,j) - gmass_stock(i,j) )
        a1 = max( 0.0, gmass_leaf(i,j)*Frac_Stock - gmass_stock(i,j) )
       !a2:maximum increment of mass_stock with resource available (g DM/cell)
        a2 = gmass_available(i,j) / RG_stock_in
       !a3:actual increment of mass_stock (g DM/cell)
        a3 = min( a1, a2 )
       !stock mass increment precedure
        gmass_available(i,j) = gmass_available(i,j) - a3 * RG_stock_in 
        gmass_stock    (i,j) = gmass_stock    (i,j) + a3               
        
        resp_grass_bg   = resp_grass_bg   + a3 * (RG_stock_in-1.0)
        mass_combust    = mass_combust    + a3 * (RG_stock_in-1.0)
        npp (p)         = npp (p)         - a3 * (RG_stock_in-1.0)
    endif
    
    !foliage growth
       !a1: maximum foliage mass requirement (g DM/cell)
        a1 = Unit_conv * (lai_opt - lai_grass(i,j)) / SLA(p)
        a1 = max(a1, 0.0)
       !a2: maximum increament of folige mass that available resource permits (g/cell)
        a2 = gmass_available(i,j) / RG_f(p)
       !a3: actual increament of foliage mass (g DM/cell)
        a3 = min(a1, a2)
       
       !foliage mass increament procedure
        gmass_available(i,j) = gmass_available(i,j) - a3 * RG_f(p)
        gmass_leaf     (i,j) = gmass_leaf     (i,j) + a3          
        lai_grass      (i,j) = gmass_leaf     (i,j) * SLA(p) / Unit_conv
        
        resp_grass_ag   = resp_grass_ag   + a3 * (RG_f(p)-1.0)
        mass_combust    = mass_combust    + a3 * (RG_f(p)-1.0)
        npp (p)         = npp (p)         - a3 * (RG_f(p)-1.0)
    
    !reproduction: excess resource will be added to litter
       if (dfl_leaf_onset(p)<= NoRepro_Phase ) cycle !No reproduction phase
       if (gmass_stock(i,j) <= Repro_Stock_min * Unit_conv ) cycle !No enough stock mass
       
       flux_litter_ag      = flux_litter_ag   + gmass_available(i,j)
       gmass_available(i,j)= 0.000000
    
END DO
END DO
!$omp end do
!$omp end parallel
END IF
   
!_____________ Give minimum resource as seeds
IF (doy==1) then
!$omp parallel private(i,j)
!$omp do private(x)
!DO i=1, DivedG !!!>>>>>>>>>>>>TN:rm
!DO j=1, DivedG !!!>>>>>>>>>>>>TN:rm
DO i=1, GRID%N_x !!!<<<<<<<<<<<<TN:add
DO j=1, GRID%N_y !!!<<<<<<<<<<<<TN:add
   x = gmass_leaf(i,j)+gmass_root(i,j)+gmass_available(i,j)+gmass_stock(i,j)
   if ( x < Grass_seeds * Unit_conv ) then
      gmass_stock(i,j) = Grass_seeds * Unit_conv
   end if
END DO
END DO
!$omp end do
!$omp end parallel
END IF
   
!_____________ carbon emission (gDM / day / stand)
   flux_c_gro_RR(1) = flux_c_gro_RR(1) + mass_combust
   
End subroutine growth_grass



!**************************************************************************************************
! Crown depth adjustment
!**************************************************************************************************
SUBROUTINE crown_adjust ()
!_____________ Set variables
!Namespace
   USE data_structure
   USE vegi_status_current1
   USE grid_status_current1
   implicit none
   
!Local Parameter (minimum crown depth in STEP)
   integer,parameter::Crown_depth_min = 10
   
!Local variables
   real    hgt          !height of tree (m)
   real    bol          !height of bole (m)
   integer count        !number of crown layers to purge
   integer p            !PFT number
   integer no,i         !for general usage
   real    x            !for general usage
   
!_____________ Determine number of crown disks that will be purged
!$omp parallel private(no)
!$omp do private(p,hgt,bol,count,x,i)
DO no=1, Max_no 
if ( .not. tree_exist(no) )                   cycle
if ( height(no)-bole(no) <= Crown_depth_min ) cycle
   
   p     = pft(no)               !PFT number
   hgt   = height(no)*STEP + 1.3 !tree height (m)
   bol   = bole  (no)*STEP + 1.3 !bole height (m)
   count = 0                     !number of crown disk that will be purged
   Select case ( Life_type(p) )
   Case (1)     !Trees of tropical rain forest
      !if ( hgt > 34.0 ) then
      !   x = 0.18
      !else
      !   x = 0.0002 * hgt * hgt - 0.0182 * hgt + 0.6
      !endif
      !x = x * hgt                    !new crown depth (m)
      
      x = ALM6(p) * hgt              !new crown height (m)
      i = int( height(no) - x/STEP ) !new bole height (STEP)
      count = max(0, i - bole(no))
      
   Case (5)     !Trees of African trees
      x = ALM6(p) * hgt              !new crown height (m)
      i = int( height(no) - x/STEP ) !new bole height (STEP)
      count = max(0, i - bole(no))
      
   Case (2,6)     !Siberian trees
      x = max( 0.0, height(no)*Step+1.3 - 10.0) !minimum bole height (m)
      i = int( (x-1.3) / STEP )                 !minimum bole height (STEP)
      
      count = max(count, i-bole(no) )
      
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add Mangrove feature ***tentative
   Case (7,8)     !Trees of tropical rain forest
      !if ( hgt > 34.0 ) then
      !   x = 0.18
      !else
      !   x = 0.0002 * hgt * hgt - 0.0182 * hgt + 0.6
      !endif
      !x = x * hgt                    !new crown depth (m)
      
      x = ALM6(p) * hgt              !new crown height (m)
      i = int( height(no) - x/STEP ) !new bole height (STEP)
      count = max(0, i - bole(no))
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
      
   Case default !default standard
      do i = 1, min(10, height(no)-bole(no)-1)
         if ( sum(npp_crownbottom(no,1:i))/i < npp_crowntop(no)*ALM4(p) ) count = i
      end do
   
   End select
   count = min(count, height(no)-bole(no)-Crown_depth_min)
   
    !when maximum_crown_depth was assumed
    ! max_depth = 30
    ! if ( height(no)-bole(no)-max_depth > 0 ) then
    !    count = max(count, height(no)-bole(no)-max_depth )
    ! end if
    
    !crown adjust procedures
    if (count > 0) then
    !  la_step          = la(no) / real( height(no)-bole(no) )
    !  la(no)           = la(no)      - la_step * count
    !  mass_leaf(no)    = mass_leaf(no)    - la_step * count / SLA(pft(no))
    !  flux_litter_leaf = flux_litter_leaf + la_step * count / SLA(pft(no))
       
       bole(no)         = bole(no)    + count
    endif
   
END DO
!$omp end do
!$omp end parallel
   
!_____________ Reset bole-height-adjust control variables
   npp_crownbottom(:,:) = 0.0
   npp_crowntop(:)      = 0.0
   
END SUBROUTINE crown_adjust



!**************************************************************************************************
! Monthly growth and reproduction of woody PFT (each individual_tree procedure)
!**************************************************************************************************
SUBROUTINE growth_trunk ()

!_____________ Set variables
!Namespace
   USE data_structure
   USE time_counter
   USE vegi_status_current1
   USE vegi_status_current2
   USE grid_status_current1
   implicit none
   
!Local parameter
   !fraction of available resource allocates to reproduction 
   real,parameter::Frac_mass_reproduction = 0.10
   
   !delay of stem growth and reproduction process after foliation (day)
   integer,parameter::Delay_from_foliation = 0
   
   !Threshold biomass that initiate reproduction (DM g)
   real,parameter::Thres_mass_reproduct = 10000.0
   
   !Buffer for available resource (DM/tree)
   real,parameter::Mass_buffer = 10.0
   
!Local variables
   real    diameter_new   !new diameter (m)
   real    mass_trunk_tmp !stem biomass for temporal usage (dm g)
   integer height_new     !new tree height (STEP)
   real    dn, md, up     !for convergence loop
   real    x, y           !for general usage
   integer no, p, count   !
   
real mass_combust
real mass_reserved

mass_combust = 0.0
mass_reserved = 0.0

!_____________ Growth process for each tree
!$omp parallel private(no)
!$omp do private(p,mass_reserved,x,y,up,dn,count,md,diameter_new,height_new,mass_trunk_tmp)
DO no=1, Max_no 

p = pft(no)

IF ( .not. tree_exist(no)                       ) cycle !when tree does not exist
IF ( .not. phenology(p)                         ) cycle !when dormance phase
IF ( dfl_leaf_onset(p)  <= Delay_from_foliation ) cycle !before specific day from foliation
IF ( flag_suppress(no) .ne. 0                   ) cycle !if suppressed
IF ( mass_available(no) <= Mass_buffer          ) cycle !when mass available doed not exist
   
   !keep some available resource
   mass_reserved     = min(mass_available(no), Mass_buffer)
   mass_available(no)= mass_available(no) - mass_reserved
   
   !Reproduction
   x = mass_leaf(no)+mass_trunk(no)+mass_root(no)+mass_stock(no)+mass_available(no)
   If (x > Thres_mass_reproduct) then
      y = mass_available(no)*Frac_mass_reproduction
      mass_available(no) = mass_available(no) - y
      flux_litter_trunk  = flux_litter_trunk  + y
   Endif
   
   !Diminish mass available
   if (Life_type(p)==1) then
      !For tropical rain trees
      x = ( (dbh_heartwood(no)+dbh_sapwood(no)) / DBH_limit(p) )**2
      x = min(1.0, x)      !diminishing factor due to size limitation
      y = mass_available(no) / max(0.01, gpp_daily_ind(no))
      y = min(1.0, y)      !diminishing factor due to NPP/GPP balance
      
      mass_combust      = mass_combust        + mass_available(no) * x * y 
      resp_leaf     (no)= resp_leaf       (no)+ mass_available(no) * x * y 
      npp            (p)= npp              (p)- mass_available(no) * x * y 
      mort_regu1    (no)= mort_regu1      (no)- mass_available(no) * x * y  !Annual NPP (gDM)
      mass_available(no)= mass_available  (no)- mass_available(no) * x * y 
      
   else if (Life_type(p)==2 .or. Life_type(p)==6) then
      !For Siberian woody PFTs
      x = ( (dbh_heartwood(no)+dbh_sapwood(no)) / DBH_limit(p) )
      x = min(1.0, x)      !diminishing factor due to size limitation
      
      mass_combust       = mass_combust       + mass_available(no) * x
      resp_leaf     (no) = resp_leaf     (no) + mass_available(no) * x
      npp            (p) = npp            (p) - mass_available(no) * x
      mort_regu1    (no) = mort_regu1    (no) - mass_available(no) * x !Annual NPP (gDM)
      mass_available(no) = mass_available(no) - mass_available(no) * x

!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add Mangrove feature ***tentative
   else if (Life_type(p)==7 .or. Life_type(p)==8) then
      !For tropical rain trees
      x = ( (dbh_heartwood(no)+dbh_sapwood(no)) / DBH_limit(p) )**2
      x = min(1.0, x)      !diminishing factor due to size limitation
      y = mass_available(no) / max(0.01, gpp_daily_ind(no))
      y = min(1.0, y)      !diminishing factor due to NPP/GPP balance
      
      mass_combust      = mass_combust        + mass_available(no) * x * y 
      resp_leaf     (no)= resp_leaf       (no)+ mass_available(no) * x * y 
      npp            (p)= npp              (p)- mass_available(no) * x * y 
      mort_regu1    (no)= mort_regu1      (no)- mass_available(no) * x * y  !Annual NPP (gDM)
      mass_available(no)= mass_available  (no)- mass_available(no) * x * y 
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
      
   end if
   
   !Allocate all remaining mass_available to sapwood production
   x = mass_available(no) !x: allocated biomass       (gDM/tree)
   y = x  / RG_s(p)       !y: increased trunk biomass (gDM/tree)
   
   resp_trunk    (no) = resp_trunk(no) + (x-y)
   npp            (p) = npp        (p) - (x-y)
   mort_regu1    (no) = mort_regu1(no) - (x-y)
   mass_combust       = mass_combust   + (x-y)
   
   mass_trunk(no)     = mass_trunk(no) +  y  
   mass_available(no) = mass_reserved        !take back reserved resource
   
!x: gap in trunk drymass (gDM)
x = mass_trunk(no) - stem_weight(dbh_heartwood(no) + dbh_sapwood(no), height(no), p)
if (x<1.0) cycle

   !determine DBH increament
  !up = 0.1000000; dn = 0.0000000 !lengh to examine convergence
   up = 0.0020000; dn = 0.0000000 !lengh to examine convergence
   count = 0
   Do
      count = count + 1
      !current sapwood diameter increament (m)
      md             = (up+dn) / 2.0
      !current diameter (m)
      diameter_new   = dbh_heartwood(no) + dbh_sapwood(no) + md
      !current height (STEP)
      height_new     = max(height(no), min(height_limit(no), h(diameter_new, p) ))
      !current trunk mass (g)
      mass_trunk_tmp = stem_weight(diameter_new, height_new, p)
      !x: difference in stem mass (g)
      x = mass_trunk_tmp - mass_trunk(no)
      !examine convergence
      if    (x < 0.0000000) then ; dn = md
      else                       ; up = md
      end if                              
      
      !Error should be less than 0.1 mm
      !This threshold value should be very low as it significantly affects woody biomass
      if (up-dn < 0.0001) exit 
      
      !Avoid infinite loop
      if (count > 100)    exit 
   End do
   
   !adjust trunc diameter and height
   height(no)      = height_new
   dbh_sapwood(no) = dbh_sapwood(no) + md 
   mort_regu4(no)  = mort_regu4(no)  + md !annual increament of DBH
   
   !adjust crown diameter ( x=crown_area(m2) )
   if     (Life_type(p)==1 .and. diameter_new>0.2) then
      !tropical evergreen trees (large)
      x = ( 25*diameter_new     )**2 * PI / 4.0
   elseif (Life_type(p)==1              )          then
      !tropical evergreen trees (small)
      x = ( (50-125*diameter_new)*diameter_new )**2 * PI / 4.0
   elseif (Life_type(p)==2              )          then
      !boreal deciduous needle leaf trees
      x = 80*diameter_new
   elseif (Life_type(p)==5              )          then
      !Tropical trees in Africa
!     x = PI * ((46.4587*diameter_new)**2)
      x = PI * ( 0.37 * (height_new*STEP+1.3) )**2
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add Mangrove feature ***tentative
   else if ((Life_type(p)==7 .or. Life_type(p)==8) .and. diameter_new>0.2) then
      !tropical evergreen trees (large)
      x = ( 25*diameter_new     )**2 * PI / 4.0
   else if (Life_type(p)==7 .or. Life_type(p)==8)   then
      !tropical evergreen trees (small)
      x = ( (50-125*diameter_new)*diameter_new )**2 * PI / 4.0
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   else
      !other trees
      x = ALM2(p) * (diameter_new**1.6)
   endif
   
   x = min( x, PI*CD_max(p)*CD_max(p)/4.0           )
   x = min( x, PI*radius_limit(no)*radius_limit(no) )
   x = max( x, 0.001                                )
   
   crown_area    (no) = x
   crown_diameter(no) = 2 * sqrt(x/PI)
   
END DO
!$omp end do
!$omp end parallel

!_____________ Carbon emission
flux_c_gro_RR(1) = flux_c_gro_RR(1) + mass_combust

END SUBROUTINE growth_trunk



!**************************************************************************************************
! Decomposition process for litter and Soil Organic Matter (SOM)
!**************************************************************************************************
SUBROUTINE decomposition (W_fi)
   
   USE data_structure
   USE grid_status_current1
   implicit none
   
!_____________ Define arguments
   real,intent(IN)   ::W_fi !field capacity   (m3/m3, 0.0 -> 1.0)

!_____________ Define local variable
   real effect_tmp      !temperature coefficient for decomposition
   real effect_mois     !moisture    coefficient for decomposition
   real d_rate_litter   !decomposition rate of litter   (month-1)
   real d_rate_som_int  !decomposition rate of som_int  (month-1)
   real d_rate_som_slow !decomposition rate of som_slow (month-1)
   real mass_decomp     !decomposed organic matter (g dm month-1 stand-1)
   real mass_combust    !biomass combust (g DM/day/stand)
   real x               !for general usage
   
mass_combust = 0.0
   
!_____________ Calculatet decomposition coefficients 
  !x: 30days running average of tmp_soil for 50cm depth (Celcius)
  !x = sum(tmp_soil_RunningRecord(1:30, 1:5)) /30.0 /5.0
  
  !x: Daily Average Soil Temperature for 0-50cm Depth (Celcius)
  !x = sum(tmp_soil_RunningRecord(1, 1:5)) /5.0
  
  !x: 30days running average of tmp_soil for the top layer (Celcius)
  !x = sum(tmp_soil_RunningRecord(1:30, 1)) /30.0
  
  !x: Daily Average Soil Temperature of the top soil layer (Celcius)
  x = tmp_soil_RunningRecord(1,1)
  
  !x: Daily Air Temperature (Celcius)
  !x = tmp_air_RunningRecord(1)
  
  !effect_tmp = exp( 0.230259 * (x-20.0) ) !Q10 method
   effect_tmp = exp( 308.56 * ( 1.0/66.02 - 1.0/( max(-45.0, x)+46.02) ) )
   
   x = (sum(pool_w(1:5))/5.0) / W_fi / Depth
   effect_mois = 0.25 + 0.75 * x
   
!_____________ Calculate decomposition rates
!  aet = sum(flux_ic_RunningRecord(1:Day_in_Year)) + &
!        sum(flux_ev_RunningRecord(1:Day_in_Year)) + &
!        sum(flux_tr_RunningRecord(1:Day_in_Year))    
!  
!  d_rate_litter   = max(0.0 , 0.000597*aet-0.0131369 ) !Meentemeyer(1978)
!  d_rate_litter   = 10**(0.0014175*aet-1.4553)         !Foley(1995)
!  d_rate_litter   = d_rate_litter / 12.0               !Annualy to Monthly
   
   d_rate_litter   = (effect_tmp * effect_mois) * (TO_litter/ 12.0)
   d_rate_som_int  = (effect_tmp * effect_mois) * (TO_fast  / 12.0)
   d_rate_som_slow = (effect_tmp * effect_mois) * (TO_slow  / 12.0)
   
   d_rate_litter   = min(1.0, d_rate_litter  )/30.0 !monthly to daily
   d_rate_som_int  = min(1.0, d_rate_som_int )/30.0 !monthly to daily
   d_rate_som_slow = min(1.0, d_rate_som_slow)/30.0 !monthly to daily
   
!_____________ Decomposition procedure1, litter
   mass_decomp = (pool_litter_trunk+pool_litter_leaf+pool_litter_root+pool_litter_ag+pool_litter_bg) * d_rate_litter
   
   pool_litter_trunk = pool_litter_trunk * (1.0 - d_rate_litter)
   pool_litter_leaf  = pool_litter_leaf  * (1.0 - d_rate_litter)
   pool_litter_root  = pool_litter_root  * (1.0 - d_rate_litter)
   pool_litter_ag    = pool_litter_ag    * (1.0 - d_rate_litter)
   pool_litter_bg    = pool_litter_bg    * (1.0 - d_rate_litter)
   
   pool_som_int     = pool_som_int     + mass_decomp * (1.00-F_air) * F_inter
   pool_som_slow    = pool_som_slow    + mass_decomp * (1.00-F_air) * (1.00 - F_inter)
   mass_combust     = mass_combust     + mass_decomp * F_air
   
!_____________ Decomposition procedure2, s.o.m. with intermediated composition rate
   mass_decomp  = pool_som_int * d_rate_som_int
   pool_som_int = pool_som_int - mass_decomp
   mass_combust = mass_combust + mass_decomp
   
!_____________ Decomposition procedure3, s.o.m. with slow composition rate
   mass_decomp   = pool_som_slow * d_rate_som_slow
   pool_som_slow = pool_som_slow - mass_decomp
   mass_combust  = mass_combust + mass_decomp
   
!_____________ Decomposition procedure4, Decay of fire-ignition-fuel 
!              an extension for African simulation
   pool_fuel_standT = pool_fuel_standT * (1.0 - 0.00075)
   pool_fuel_standG = pool_fuel_standG * (1.0 - 0.00075)
   
!_____________ Record carebon emission
   flux_c_htr_RR(1) = flux_c_htr_RR(1) + mass_combust
   
END SUBROUTINE decomposition
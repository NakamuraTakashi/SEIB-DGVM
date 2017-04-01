!**************************************************************************************************
! Air properties
!**************************************************************************************************
SUBROUTINE air (tmp_air, humid, ALT)

!_____________ Set variables
!Namespace
   USE data_structure
   USE grid_status_current2
   implicit none
   
!Augments
   real,intent(IN) ::tmp_air   !2m air temperature (Celcius)
   real,intent(IN) ::humid     !SPecific humidity (kg kg-1)
   real,intent(IN) ::ALT       !altitude (m above MSL)
   
!Local variables
   real a1, a2  !for temporal memory
   
!_____________ Main part
!ap: air pressure (hPa)
   ap = 1013.25 * exp ( -0.2838472*ALT /(8.3144*(tmp_air+ZAT)) )
   
!vp: vapour pressure (hPa)
   vp = ap * humid / ( 0.622 + 0.378 * humid )
   
!vp_sat: saturated vapour pressure (hPa)
   if(tmp_air > 0.0)then   !@water surface
      vp_sat = 6.1078*(10.0**((7.5*tmp_air)/(237.3+tmp_air)))
   else                    !@ice surface
      vp_sat = 6.1078*(10.0**((9.5*tmp_air)/(265.3+tmp_air)))
   endif
   vp_sat = max(0.0, vp_sat)
   
!vpd: vapour pressure deficit (hPa)
   vpd = max(0.0, vp_sat - vp)
   
!slope_vps: slope of saturated vapour pressure as a function of temperature
   if(tmp_air > 0.0)then  !  at  water  surface
      a1 = 6.1078 * (2500.0-2.4*tmp_air)
      a2 = ( 10.0**((7.5*tmp_air)/(237.3+tmp_air)) )
   else                 !  at  ice  surface
      a1 = 6.1078 * 2834.0
      a2 = (10.0**((9.5*tmp_air)/(265.3+tmp_air)))
   endif
   slope_vps = ( a1 / (0.4615*(ZAT+tmp_air)*(ZAT+tmp_air)) ) * a2
   
!dnsa: air density
   dnsa = 1.293 * (ZAT / (tmp_air+ZAT)) * (ap / 1013.25) * (1.0 - 0.378 * vp / ap)
   
END SUBROUTINE air



!**************************************************************************************************
! Calculate seasonal chages in radiation-related-variables as a function of latitude
!**************************************************************************************************
SUBROUTINE radiation_seasonal_change (LAT)

!_____________ Set variables
!Namespace
   USE data_structure
   USE vegi_status_current1
   USE grid_status_current2
   implicit none
   
!Augments
   real   ,intent(IN)::LAT        !latitude  (degree)
   
!Local variables
   real    sl_dec     !solar  declination (degree)
   real    ha         !angle from sun-rise to median passage (degree)
   integer doy        !day of year (1~Day_in_Year)
   real    x, a1, a2  !for general usage
   
!_____________ Main part
Do doy=1, Day_in_Year
   !sl_dec: solar declination (degree)
      sl_dec = 23.45 * sin(DtoR * 360.0 * ( real(doy) - 81.0 ) / Day_in_Year)
   
   !sl_hgt: solar hight at midday (degree)
      x           = sin(LAT * DtoR) * sin(sl_dec * DtoR) + &
                    cos(LAT * DtoR) * cos(sl_dec * DtoR)
      x           = min(1.0, max(-1.0, x))
      sl_hgt(doy) = asin(x) * RtoD
   
   !dlen: day length (hour)
      if (sl_hgt(doy) <= 0.1) then
         dlen(doy)=0.0
      else
         x           = -tan(LAT * DtoR) * tan(sl_dec * DtoR)
         ha          = RtoD * acos( min(1.0, max(-1.0, x)) ) !angle from sun-rise to median passage
         dlen (doy)  = 2.0 * (ha / 15.0)
      endif
      dlen (doy)  = max(min(dlen(doy), 24.0), 0.0)
   
   !rad_stratosphere: shortwave radiation at the atmosphere-top (W/m2) 
      a1 = DtoR * 360.0 * (real(doy)/Day_in_Year) !seasonal angle of the earth's orbit (radian)
      a2 = 1.00011 + 0.034221 * cos(a1) + 0.00128 * sin(a1) + &
           0.000719 * cos(2.0*a1) + 0.000077 * sin(2.0*a1)
      rad_stratosphere (doy) = max(0.0, 1367.0 * sin(sl_hgt(doy) * DtoR) * a2)
   
End Do

END SUBROUTINE radiation_seasonal_change



!**************************************************************************************************
! Properties for shortwave radiation
!**************************************************************************************************
SUBROUTINE radiation (LAT, cloud)

!_____________ Set variables
!Namespace
   USE data_structure
   USE time_counter
   USE vegi_status_current1
   USE grid_status_current2
   implicit none
   
!Augments
   real   ,intent(IN) ::LAT       !latitude  (degree)
   real   ,intent(IN) ::cloud     !total cloudness (fraction)
   
!Local variables
   real rad_diffuse   !diffused radiation in rad (W/m2)
   real rad_direct    !direct radiation in rad (W/m2)
   real x             !for general usage
   
!_____________ Main part
!trap for high latitude regions
if (sl_hgt(doy)<=0.1) then
   rad         = 0.0
   par         = 0.0
   par_diffuse = 0.0
   par_direct  = 0.0
   return
endif
    
!rad: Shortwave radiation at midday at the ground surface (W/m2)
   x = 0.803 - 0.34*cloud - 0.458*cloud*cloud  !Black's
   x = 0.8964 - 0.5392 * cloud   !new regression based on NCEP/NCAR data
   x = min(1.0,  max(0.0, x))
   rad = rad_stratosphere(doy) * x
   
!par: photosynthetically active radiation in mid_day (micro mol photon m-2 s-1)
!   based  on  the  empirical  Tooming's  equation
   if(rad_stratosphere(doy) /= 0.0) then
      rad_diffuse = max( 0.0 , rad * (0.958 - 0.982*(rad/rad_stratosphere(doy))) )
      rad_direct  = max( 0.0 , rad - rad_diffuse )
      
      par_diffuse = 4.2 * 0.57 * rad_diffuse
      par_direct  = 4.6 * 0.43 * rad_direct 
      par         = par_diffuse + par_direct
   else
      par         = 0.0
      par_diffuse = 0.0
      par_direct  = 0.0
   endif
   
END SUBROUTINE radiation



!**************************************************************************************************
! Albedo
!**************************************************************************************************
SUBROUTINE albedo_calc (Albedo_soil0)

!_____________ Set variables
!Namespace
!USE data_structure
   USE vegi_status_current1
!USE vegi_status_current2
   USE grid_status_current1
   USE grid_status_current2
   implicit none
   
!Augments
   real,intent(IN)::Albedo_soil0 !albedo, default
   
!Local variables
   real    albedo_leaf0 !leaf albedo @ snow free condition
   
!_____________ Main part
!albedo_soil: soil surface albedo
   albedo_soil = Albedo_soil0 + (0.7-Albedo_soil0) / (1+exp(-0.05*(pool_snow-70.0)))
   albedo_soil = max(0.05, min(0.70, albedo_soil))
   
!albedo_leaf: leaf surface albedo (temporal)
   select case (biome)
      case (3:8)   !default forests (Jones 1992)
         albedo_leaf0 = 0.15
      case (9)     !larch forests (Fukuda personal observation)
         albedo_leaf0 = 0.107
      case default !grass land, desert, etc (Jones 1992)
         albedo_leaf0 = 0.24
   end select
   
   albedo_leaf = albedo_leaf0 + (0.7-albedo_leaf0) / (1+exp(-0.05*(pool_snow-70.0)))
   albedo_leaf = max(0.05, min(0.70, albedo_leaf))
   
!albedo_mean
   albedo_mean = albedo_leaf * (1.0-ir_tree*ir_grass) + albedo_soil * (ir_tree*ir_grass)
   albedo_mean = max(0.0, min(999.9, albedo_mean))
   
END SUBROUTINE albedo_calc



!**************************************************************************************************
! Net radiation
!**************************************************************************************************
SUBROUTINE net_radiation (tmp_air, cloud)

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
   real,intent(IN)::tmp_air      !2m air temperature (Celcius)
   real,intent(IN)::cloud        !total cloudness, fraction
   
!_____________ Main part
!radnet_long: net longwave radiation (+:Upward direction)
   radnet_long = (5.6703/100000000.0) * ((tmp_air + ZAT)**4.0) &
                  * (0.39 + 0.058/(vp+1.0)) * (1.0 - 0.65*cloud)
   
!radnet_veg: net radiation of plant canopy (W m-2, day time mean)
!            0.5-> convertion coefficient from midday to daytime mean value
   radnet_veg = 0.5 * rad * (1.0-albedo_leaf) * (1.0 - ir_tree*ir_grass) &
              - radnet_long                   * (1.0 - ir_tree*ir_grass)  
   
!radnet_soil: net radiation of soil surface (W m-2, whole day mean)
!            0.5-> convertion coefficient from midday to daytime mean value
   radnet_soil = 0.5 * rad * (1.0-albedo_soil) * (ir_tree*ir_grass) * dlen(doy) / 24.0 &
               - radnet_long                   * (ir_tree*ir_grass)
   
END SUBROUTINE net_radiation



!**************************************************************************************************
! Water budget module, daily computation
!**************************************************************************************************
SUBROUTINE waterbudget (W_fi, W_wilt, prec, wind, tmp_air, tmp_soil)

!_____________ Set variables
!Namespace
   USE data_structure
   USE time_counter
   USE vegi_status_current1
   USE vegi_status_current2
   USE grid_status_current1
   USE grid_status_current2
   implicit none
   
!Local parameters
   !Soil Conductance Coefficient (determined by just only parametarization)
   real,parameter::C_soil_coefficient = 0.00150
!  real,parameter::C_soil_coefficient = 0.00300
   
!Augments
   real,intent(IN)::tmp_air           !air temperature (Celcius)
   real,intent(IN)::tmp_soil(NumSoil) !soil temperature for each layer (Celcius)
   real,intent(IN)::prec              !precipitation,  mm day-1
   real,intent(IN)::wind              !wind velocity (m s-1)
   
   real,intent(IN)::W_fi     !filed capacity   (m3/m3, 0.0 -> 1.0)
   real,intent(IN)::W_wilt   !wilting point    (m3/m3, 0.0 -> 1.0)
   
!Local variables
   real flux_rain  !rain fall (mm/day)
   real flux_snow  !snow fall (mm/day)
   real flux_tw    !snow thaw water (mm/day)
   real flux_sfc   !available water on soil surface (mm/day)
   
   real flux_ic    !intercepted rain (mm/day)
   real flux_ev    !evaporation flux (mm/day)
   real flux_tr    !transpiration flux (mm/day)
   real flux_ro    !runoff flux (mm/day)
   
   real flux_ev_pm !potential evaporation rate  (mm/day)
   real flux_tr_pm !potential transpiration rate  (mm/day)
   
   real rain_no    !expected number of rain event (day-1)
   
   real c_aero     !aerodynamic conductance (s-1 m)
   real c_soil     !soil conductance        (s-1 m)
   real c_leaf     !canopy conductance      (s-1 m)
   
   real lh         !Latent heat of water vaporization (MJ/kg H2O)
   real gamma      !Psychrometer constant (hPa k-1), around 0.667
   
   real pool_total          !for checker (total water available in this virtual field)
   real pool_total_reserved !for checker (total water available in this virtual field)
   
   real    water_available
   real    x, y, a, b, c    !for general usage
   integer layer, no, p     !for counter
   
!_____________ initialize daily hydrological variables
   flux_ro = 0.0
   flux_ic = 0.0
   flux_ev = 0.0
   flux_tr = 0.0
   
!_____________ lh: heat of water vaporization at 1atm (J/g H2O)
   lh = 2500.25 - 2.365 * tmp_air   ! Fritschen and Gay (1979)
   
!_____________ Psychrometer Coefficient (hPa k-1), around 0.667
   ! 1004.: Specific-Heat at Constant-Pressure of Dry-Air (J/kg/K)
   ! 0.622: fraction of 
   !        18.015 [molecular mass of water vapor (g/mol)] to
   !        28.966 [molecular mass of     dry air (g/mol)].  
   gamma = 1004. * ap / (lh*1000) / 0.622
   
!_____________ Water-leakage checker#0 (record total water available on the virtual forest)
   pool_total = prec + pool_snow + sum(pool_w(:))
   pool_total_reserved = pool_total
   
!_____________ Calculate Potential Eveporation & Transpiration
   
   !*** Preparetion for Evapotranspiration
   !c_aero: aerodynamic conductance (s m-1)
   c_aero = (1+0.537*wind) / 250.08
   
   !c_soil: surface soil conductance (s m-1)
   x = sum(pool_w(1:5)) / (5*W_fi*Depth) ; x = min(1.0, x)
   c_soil = C_soil_coefficient * (x**2.0) !Coefficient -> Just Only Parametarization
   
   !c_leaf: canopy conductance (m s-1)
   c_leaf = 0.0224 * canopy_cond !1 mol gas = 0.0224 m3 at 0C and 1atm
   
   !*** flux_tr_pm: Potential transpiration rate
   if (c_leaf>0.0001) then
      !a: sum of daily netradiation (J m-2 day-1)
      a = radnet_veg           * dlen(doy) * 3600 / (lh*1000)
      b = dnsa * 1012 * vpd * c_aero * dlen(doy) * 3600 / (lh*1000)
      c = slope_vps + gamma * (1.0+c_aero/c_leaf)
      
      flux_tr_pm = ( a*slope_vps + b ) / c
      !     3600   (s hr-1) = secound -> hour
      !     1000   (g mm-2) = 1mm equivalent water for evaporation
   else
      flux_tr_pm = 0.0
   endif
   
   !*** flux_ev_pm: Potential evaporation rate
   if (c_soil>0.0001) then
      !a: sum of daily netradiation (J m-2 day-1)
      a = radnet_soil                * 24 * 60 * 60 / (lh*1000)
      b = dnsa * 1012 * vpd * c_aero * 24 * 60 * 60 / (lh*1000)
      c = slope_vps + gamma * (1.0+c_aero/c_soil)
      
      flux_ev_pm = ( a*slope_vps + b ) / c
   else	
      flux_ev_pm = 0.0
   endif
   
  !*** ev_pot_RunningRecord: Potential evapotranspiration rate with Priestley-Taylor model
   ev_pot_RunningRecord(2:Day_in_Year) = ev_pot_RunningRecord(1:Day_in_Year-1)
   
   
   a = (radnet_veg*dlen(doy) + radnet_soil*24)*60*60 !a: sum of daily netradiation (J m-2 day-1)
   ev_pot_RunningRecord (1) = 1.26 * (slope_vps/(slope_vps+gamma)) * a / (lh*1000)
   !1.26            = an empirical coefficient of Priestley and Taylor
   !60              = secound -> min, min -> hour
   !1000   (g mm-2) = 1mm equivalent water for evaporation
   
   !Other method
   !ev_pot_RunningRecord (1) = flux_tr_pm + flux_ev_pm
   
!_____________ Fluxus of snow, rain, and thaw water
   !x: proportion of snow in precipitation
   x = 1.0 / ( 1.0 + exp(0.75*tmp_air - 1.5) )
   
   !y: proportion of snow melt in this day
   y = 1.0 / ( 1.0 + exp( -0.3*(tmp_soil(1)-10.0) ) )
   
   !water balance
   flux_snow   = prec * x
   flux_rain   = prec * (1.0-x)
   pool_snow   = pool_snow + flux_snow
   
   flux_tw     = pool_snow  * y
   pool_snow   = pool_snow - flux_tw
   
!   !Water-leakage checker #1
!   y = flux_rain + flux_tw + pool_snow + sum(pool_w(:))
!   if ( Logging==1 .and. abs(pool_total-y)>0.1 ) then
!      write(*,'(a,2f7.1)') 'error@waterbudget_1', pool_total, y
!   endif
   
!_____________ Canopy Interception
   !determine fraction interception
   rain_no = 0.1 * flux_rain
   if ( flux_tr_pm + flux_ev_pm <= 50.0 * 12.0 / Day_in_Year ) then
      rain_no = min( rain_no,  5.0 * 12.0 / Day_in_Year )
   else
      rain_no = min( rain_no, 10.0 * 12.0 / Day_in_Year )
   endif
   
   flux_ic = 3.0 * rain_no * ( 1.0 - exp(-1.0*lai) )
   flux_ic = min(flux_rain, flux_ic)          
   
   !water balance
   !flux_ic: intercepted water (mm)
   flux_sfc = flux_rain - flux_ic + flux_tw 
   
!   !Water-leakage cheker #2
!   pool_total = pool_total - flux_ic
!   y = flux_sfc + pool_snow + sum(pool_w(:))
!   if ( Logging==1 .and. abs(pool_total-y)>0.1 )  then
!      write(*,'(a,2f7.1)') 'error@waterbudget_2', pool_total, y
!   endif
   
!_____________ Percolation
   Do layer=1, NumSoil
      if (tmp_soil(layer) < 0.0) exit
      x = min(flux_sfc, W_fi*Depth - pool_w(layer)) !x: water that add this layer (mm)
      pool_w(layer) = pool_w(layer) + x
      flux_sfc      = flux_sfc      - x
   End do
   
   flux_ro  = flux_sfc !recorder
   flux_sfc = 0.0
   
!   !Water-leakage cheker #4
!   pool_total = pool_total - flux_ro !checker
!   y = pool_snow + sum(pool_w(:))
!   if ( Logging==1 .and. abs(pool_total-y)>0.1 )  then
!      write(*,'(a,2f7.1)') 'error@waterbudget_4', pool_total, y
!   endif
   
!_____________ Actual transpiration (flux_tr)
   flux_tr_pm = max(0.0, flux_tr_pm - flux_ic)
   
   !Count layers where transpiration occurs
   no = 1
   do p=1, PFT_no
      if (pft_exist(p)) then
         no = max(no, RootDepth(p))
      endif
   end do
   
   !water_available: Available water for transpiration (mm)
   water_available = 0.00
   do layer=1, no
   if ( tmp_soil(layer) >= 0.0 ) then
      x = pool_w(layer) - Depth*W_wilt
      water_available = water_available + max(x, 0.0)
   endif
   enddo
   
   !Compute actual transpiration
   a       = 0.85 !a = 0.1
   b       = water_available + flux_tr_pm
   c       = water_available * flux_tr_pm
   flux_tr = (b - sqrt( b*b - 4.0*a*c )) / (2.0*a)
   flux_tr = max(flux_tr, 0.0)
   flux_tr = min(flux_tr, water_available)
   
   !** water balance (transpiration water will be paied from the soil-top-layer)
   x = flux_tr
   do layer=1, no
   if ( tmp_soil(layer) >= 0.0 ) then
      y = max( 0.0, pool_w(layer) - Depth*W_wilt ) !water availbale at the present layer
      y = min (x, y)                               !actual water transpirated from this layer
      pool_w(layer) = pool_w(layer) - y
      x             = x             - y
   endif
   enddo
   flux_tr = flux_tr - x
   
!   !Water-leakage cheker #5
!   pool_total = pool_total - flux_tr
!   y = pool_snow + sum(pool_w(:))
!   if ( Logging==1 .and. abs(pool_total-y)>0.1 )  then
!      write(*,'(a,2f7.1)') 'error@waterbudget_5', pool_total, y
!   endif
   
!_____________ Actual evaporation
   water_available = 0.00
   do layer=1, 5
   if ( tmp_soil(layer) >= 0.0 ) then
      water_available = water_available + max(pool_w(layer), 0.0)
   endif
   enddo
   
   a = 0.85 !a = 0.1
!  b = pool_w(1) + flux_ev_pm
!  c = pool_w(1) * flux_ev_pm
   b = water_available + flux_ev_pm
   c = water_available * flux_ev_pm
   flux_ev = (b - sqrt( b*b - 4.0*a*c )) / (2.0*a)
   flux_ev = max(flux_ev, 0.0      )
!  flux_ev = min(flux_ev, pool_w(1))
   flux_ev = min(flux_ev, water_available)
!  pool_w(1) = pool_w(1) - flux_ev
   
   x = flux_ev
   do layer=1, 5
   if ( tmp_soil(layer) >= 0.0 ) then
      y = min (x, pool_w(layer))      !actual water transpirated from this layer
      pool_w(layer) = pool_w(layer) - y
      x             = x             - y
   endif
   enddo
   flux_ev = flux_ev - x
   
!   !Water-leakage cheker #6
!   pool_total = pool_total - flux_ev
!   y = pool_snow + sum(pool_w(:))
!   if ( Logging==1 .and. abs(pool_total-y)>0.1 )  then
!      write(*,'(a,2f7.1)') 'error@waterbudget_6', pool_total, y
!   endif
   
!_____________ Update running records
   flux_ro_RunningRecord(2:Day_in_Year) = flux_ro_RunningRecord(1:Day_in_Year-1)
   flux_ic_RunningRecord(2:Day_in_Year) = flux_ic_RunningRecord(1:Day_in_Year-1)
   flux_ev_RunningRecord(2:Day_in_Year) = flux_ev_RunningRecord(1:Day_in_Year-1)
   flux_tr_RunningRecord(2:Day_in_Year) = flux_tr_RunningRecord(1:Day_in_Year-1)
   
   flux_ro_RunningRecord(1) = flux_ro
   flux_ic_RunningRecord(1) = flux_ic
   flux_ev_RunningRecord(1) = flux_ev
   flux_tr_RunningRecord(1) = flux_tr
   
   pool_w1_RunningRecord(2:Day_in_Year) = pool_w1_RunningRecord(1:Day_in_Year-1)
   pool_w1_RunningRecord(1)             = sum(pool_w(1:5)) / 5.0
   
!_____________ Water-leakage checker #7 (Final check)
   x = pool_total_reserved - (pool_snow + sum(pool_w(:))) &
     - (flux_ro + flux_ic + flux_ev + flux_tr)
   if ( abs(x)>0.01 )  then
      write(*,'(a,2f8.2)') 'error@waterbudget water leakage occured!', pool_total_reserved, x
   endif
   
END SUBROUTINE waterbudget



!*************************************************************************************************
! light interruption index (daily computation)
!*************************************************************************************************
SUBROUTINE ir_index ()

!_____________ Set variables
!Namespace
   USE data_structure
   USE vegi_status_current1
   USE vegi_status_current2
   USE grid_status_current2
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none
   
!Local variables
   integer no, p
   real    x, y
   
!_____________ Main part
   ir_tree        = 0.00
   ir_grass       = 0.00
   
IF ( par > 1.0 ) then
   !ir_tree: radiation-interruption-coefficient by tree crown
   x = 0.0 ; y = 0.0 
   do no = 1, Max_no
      if (tree_exist(no)) then
         p = pft(no)
         if (phenology(p)) then
            x = x + la(no) * max(EK0(p),eK(p)) !sum of light attenuation index (direct radiation)
            y = y + la(no) * EK0(p)            !sum of light attenuation index (diffused radiation)
         endif
      endif
   enddo
!   x = x / Max_loc / Max_loc !!!>>>>>>>>>>>>TN:rm
!   y = y / Max_loc / Max_loc !!!>>>>>>>>>>>>TN:rm
   x = x / real(GRID%Area) !!!<<<<<<<<<<<<TN:add
   y = y / real(GRID%Area) !!!<<<<<<<<<<<<TN:add
   
   ir_tree = (par_direct/par)*exp(-1.0*x) + (par_diffuse/par)*exp(-1.0*y)
   ir_tree = min(1.0, max(0.0, ir_tree) )
   
   !ir_grass: radiation-interruption-coefficient by grass leaf
   if (pft_exist(C3g_no)) then
      p=C3g_no
   else
      p=C4g_no
   endif
   
!   x = sum(lai_grass(:,:)) / DivedG / DivedG !!!>>>>>>>>>>>>TN:rm
   x = sum(lai_grass(:,:)) / real(GRID%N_tot) !!!<<<<<<<<<<<<TN:add
   
   ir_grass = (par_direct  / par) * exp(-1.0 * x *  eK(p)) + &
              (par_diffuse / par) * exp(-1.0 * x * EK0(p))    
   ir_grass = min(1.0, max(0.0, ir_grass) )
   
ENDIF

END SUBROUTINE ir_index

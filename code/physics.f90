!**************************************************************************************************
! Air properties
!**************************************************************************************************
SUBROUTINE air (tmp_air, rh, ALT)

!_____________ Set variables
!Namespace
   USE data_structure
   USE grid_status_current2
   implicit none
   
!Augments
   real,intent(IN)  ::tmp_air   !2m air temperature (Celcius)
   real,intent(IN)  ::rh        !relative humidity (kg  kg-1)
   real,intent(IN)  ::ALT       !altitude (m above MSL)
   
!_____________ Main part
!ap: air pressure (hPa)
   ap = 1013.25 * exp ( -0.2838472*ALT /(8.3144*(tmp_air+ABS_ZERO)) )
   
!vp_sat: saturated vapour pressure (hPa)
   vp_sat = Sat_vp (tmp_air)
   
!vp: vapour pressure (hPa)
   vp = vp_sat * rh / 100.0
   
!vpd: vapour pressure deficit (hPa)
   vpd = vp_sat - vp
   
!slope_vps: slope of saturated vapour pressure as a function of temperature
   slope_vps = Delta_sat_vp (tmp_air)
   
!dnsa: air density
   dnsa = 1.293 * (ABS_ZERO / (tmp_air+ABS_ZERO)) * (ap / 1013.25) * (1.0 - 0.378 * vp / ap)
   
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
   real    ha         !angle from sun-rise to median passage (degree)
   integer doy        !day of year (1~Day_in_Year)
   real    x, a1, a2  !for general usage
   
!_____________ Main part
Do doy=1, Day_in_Year
   !sl_dec: solar declination (degree)
      sl_dec(doy) = 23.45 * sin(DtoR * 360.0 * ( real(doy) - 81.0 ) / Day_in_Year)
   
   !sl_hgt: solar hight at midday (degree)
      x           = sin(LAT * DtoR) * sin(sl_dec(doy) * DtoR) + &
                    cos(LAT * DtoR) * cos(sl_dec(doy) * DtoR)
      x           = min(1.0, max(-1.0, x))
      sl_hgt(doy) = asin(x) * RtoD
   
   !dlen: day length (hour)
      if (sl_hgt(doy) <= 0.1) then
         dlen(doy)=0.0
      else
         x           = -tan(LAT * DtoR) * tan(sl_dec(doy) * DtoR)
         ha          = RtoD * acos( min(1.0, max(-1.0, x)) ) !angle from sun-rise to median passage
         dlen (doy)  = 2.0 * (ha / 15.0)
      endif
   
   !rad_stratosphere: shortwave radiation at the atmosphere-top (W/m2) 
      a1 = DtoR * 360.0 * (real(doy)/Day_in_Year) !seasonal angle of the earth's orbit (radian)
      a2 = 1.00011 + 0.034221 * cos(a1) + 0.00128 * sin(a1) + &
           0.000719 * cos(2.0*a1) + 0.000077 * sin(2.0*a1)
      rad_stratosphere (doy) = max(0.0, 1367.0 * sin(sl_hgt(doy) * DtoR) * a2)
   
End Do

END SUBROUTINE radiation_seasonal_change



!**************************************************************************************************
! Radiation properties
!**************************************************************************************************
SUBROUTINE radiation (LAT, rad_short, cloud)

!_____________ Set variables
!Namespace
   USE data_structure
   USE time_counter
   USE vegi_status_current1
   USE grid_status_current2
   implicit none
   
!Augments
   real   ,intent(IN) ::LAT       !latitude  (degree)
   real   ,intent(IN) ::rad_short !shortwave radiation @ midday (W/m2)
   real   ,intent(OUT)::cloud     !total cloudness (fraction)
   
!Local variables
   real    rad_diffuse   !diffused radiation in rad (W/m2)
   real    rad_direct    !direct radiation in rad (W/m2)
   real    x             !for general usage
   integer dom_mid       !midday of the month
   
!_____________ Main part
!trap for high latitude regions
if (sl_hgt(doy)<=0.1) then
   par         = 0.0
   par_diffuse = 0.0
   par_direct  = 0.0
   return
endif
    
!rad_short: Shortwave radiation at midday at the ground surface (W/m2)
   !x = 0.803 - 0.34*cloud - 0.458*cloud*cloud  !Black's
   !x = 0.8964 - 0.5392 * cloud   !new regression based on NCEP/NCAR data
   !x = min(1.0,  max(0.0, x))
   !rad_short = rad_stratosphere(doy) * x
   
   dom_mid = sum( Day_in_month(1:Month(doy)) )-15 !Mid day of the month
   x       = rad_short / max(0.01, rad_stratosphere(dom_mid))
   cloud   = (0.8964 - x) / 0.5392
   cloud   = max(0.0, min(cloud, 1.0))
   
!par: photosynthetically active radiation in mid_day (micro mol photon m-2 s-1)
!   based  on  the  empirical  Tooming's  equation
   if(rad_stratosphere(doy) /= 0.0) then
      rad_diffuse = max( 0.0 , rad_short * (0.958 - 0.982*(rad_short/max(0.01,rad_stratosphere(dom_mid)))) )
      rad_direct  = max( 0.0 , rad_short - rad_diffuse )
      
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
SUBROUTINE net_radiation (tmp_air, rad_short, rad_long, cloud)

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
   real,intent(IN)::rad_short    !downward shortwave radiation (W/m2)
   real,intent(IN)::rad_long     !downward longwave  radiation (W/m2)
   real,intent(IN)::cloud        !total cloudness, fraction
   
!_____________ Main part
!radlong_up: longwave radiation (+:Upward direction) PREVIOUS 
   radlong_up = (5.6703/100000000.0) * ((tmp_air + ABS_ZERO)**4.0) - rad_long
  !Previous way of calculation
  !radlong_up = (5.6703/100000000.0) * ((tmp_air + ABS_ZERO)**4.0) &
  !           * (0.39 + 0.058/(vp+1.0)) * (1.0 - 0.65*cloud)
   
!radnet_veg: net radiation of plant canopy (W m-2, day time mean)
!            0.5-> convertion coefficient from midday to daytime mean value
   radnet_veg = 0.5 * rad_short * (1.0-albedo_leaf) * (1.0 - ir_tree*ir_grass) * dlen(doy) / 24.0 &
              - radlong_up                    * (1.0 - ir_tree*ir_grass)  
   
!radnet_soil: net radiation of soil surface (W m-2, whole day mean)
!            0.5-> convertion coefficient from midday to daytime mean value
   radnet_soil = 0.5 * rad_short * (1.0-albedo_soil) * (ir_tree*ir_grass) * dlen(doy) / 24.0 &
               - radlong_up                    * (ir_tree*ir_grass)
   
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
   lh = Latent_heat (tmp_air) / (1000.0)
   
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
   
   x       = sum(lai_RunningRecord(1,:)) !LAI for all PFTs
   flux_ic = 3.0 * rain_no * ( 1.0 - exp(-1.0*x) )
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
   flux_ro1_RunningRecord(2:Day_in_Year) = flux_ro1_RunningRecord(1:Day_in_Year-1)
   flux_ro2_RunningRecord(2:Day_in_Year) = flux_ro2_RunningRecord(1:Day_in_Year-1)
   flux_ic_RunningRecord (2:Day_in_Year) = flux_ic_RunningRecord (1:Day_in_Year-1)
   flux_ev_RunningRecord (2:Day_in_Year) = flux_ev_RunningRecord (1:Day_in_Year-1)
   flux_tr_RunningRecord (2:Day_in_Year) = flux_tr_RunningRecord (1:Day_in_Year-1)
   flux_sl_RunningRecord (2:Day_in_Year) = flux_sl_RunningRecord (1:Day_in_Year-1)
   flux_tw_RunningRecord (2:Day_in_Year) = flux_tw_RunningRecord (1:Day_in_Year-1)
   flux_sn_RunningRecord (2:Day_in_Year) = flux_sn_RunningRecord (1:Day_in_Year-1)
   
   flux_ro1_RunningRecord(1) = flux_ro
   flux_ro2_RunningRecord(1) = 0.0
   flux_ic_RunningRecord (1) = flux_ic
   flux_ev_RunningRecord (1) = flux_ev
   flux_tr_RunningRecord (1) = flux_tr
   flux_sl_RunningRecord (1) = 0.0
   flux_tw_RunningRecord (1) = flux_tw
   flux_sn_RunningRecord (1) = flux_snow
   
   pool_w1_RunningRecord(2:Day_in_Year) = pool_w1_RunningRecord(1:Day_in_Year-1)
   pool_w1_RunningRecord(1)             = sum(pool_w(1:5)) / 5.0
   
!_____________ Water-leakage checker #7 (Final check)
   x = pool_total_reserved - (pool_snow + sum(pool_w(:))) &
     - (flux_ro + flux_ic + flux_ev + flux_tr)
   if ( abs(x)>0.01 )  then
      write(*,'(a,2f8.2)') 'error@waterbudget water leakage occured!', pool_total_reserved, x
   endif
   
! !TMP
! 1 format(E10.4e1,a)  !Output format (Sample: "0.284E+3"+",")
! 2 format(E10.4e1  )  !Output format (Sample: "0.284E+3"    )
! write(*, '(i3,a)', advance='no') year, ","
! write(*, '(i3,a)', advance='no') doy, ","
! write(*, 1, advance='no') prec     , ","
! write(*, 1, advance='no') flux_ro  , ","
! write(*, 1, advance='no') flux_ic  , ","
! write(*, 1, advance='no') flux_ev  , ","
! write(*, 1, advance='no') flux_tr  , ","
! write(*, 1, advance='no') flux_rain, ","
! write(*, 1, advance='no') flux_snow, ","
! write(*, 2, advance='no') flux_tw       
! write(*,*)

END SUBROUTINE waterbudget

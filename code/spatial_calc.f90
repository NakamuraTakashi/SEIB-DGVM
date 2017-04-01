!**************************************************************************************************
! Relative intensity of diffused PAR for each forest layer
!**************************************************************************************************
SUBROUTINE diffused_radiation ()

!_____________ Set variables
!Namespace
   USE data_structure
   USE vegi_status_current1
   USE grid_status_current2
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none
   
!Local variables
   real,dimension(PFT_no,Max_hgt)::lai_each  !LAI for each PFT for each forest layer (m2/m2)
   real,dimension(Max_hgt)     ::attenuation_idx !attenuation index for each forest layer
   integer no, i, p !loop counter
   real    x        !for general usage
   
!_____________ Leaf Area Index for each layer for each PFT (lai_each)
   lai_each(:,:) = 0.0
   DO no = 1, Max_no
   if ( .not. tree_exist(no)     ) cycle
   if ( .not. phenology(pft(no)) ) cycle
      
      x = la(no) / real(height(no)-bole(no)) !leaf area per each foliage step (m2/STEP)
      
      do i = bole(no)+1, height(no)
         lai_each(pft(no), i) = lai_each(pft(no), i) + x
      end do
   END DO
!   lai_each(:,:) = lai_each(:,:) / real(Max_loc) / real(Max_loc)  !!!>>>>>>>>>>>>TN:rm
   lai_each(:,:) = lai_each(:,:) / real(GRID%Area)                 !!!<<<<<<<<<<<<TN:add
   
!_____________ Light attenuation coefficent for each forest layer
   DO i = 1, Max_hgt
      x = 0.0
      do p = 1, PFT_no
         x = x + lai_each(p, i) * EK0(p)
      end do
      attenuation_idx(i) = x
   END DO
   
!_____________ Relative intensity of diffuse PAR for each forest layer (par_diffuse_rel)
   par_diffuse_rel(Max_hgt) = 1.0 !top layer
   DO i = Max_hgt, 2, -1
      x = par_diffuse_rel(i) * exp( -1.0 * attenuation_idx(i) )
      par_diffuse_rel(i-1) = x
   END DO
   
END SUBROUTINE diffused_radiation



!**************************************************************************************************
! Fraction of tree crown coverage 
!**************************************************************************************************
SUBROUTINE crown_coverage ()

!_____________ Set variables
!Namespace
   USE data_structure
   USE vegi_status_current1
   USE grid_status_current2
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none
   
!Local variables
   real    x
   integer no
   real    d1, d2, d3, d4, d5, d6, d7, d8, d9, proxy, y
   integer count, i, j
   logical flag_coverage
   
!_____________ Main(Method1)
   count = 0
!   Do i=1, Max_loc    !!!>>>>>>>>>>>>TN:rm
!   Do j=1, Max_loc    !!!>>>>>>>>>>>>TN:rm
   Do i=1, GRID%Max_x  !!!<<<<<<<<<<<<TN:add
   Do j=1, GRID%Max_y  !!!<<<<<<<<<<<<TN:add
      flag_coverage = .false.
      
      do no = 1, Max_no
      if ( flag_coverage       ) cycle
      if (.not. tree_exist(no) ) cycle
      if ( mass_leaf(no)==0.0  ) cycle
         
         !d1~9: Square distance between centers for each "mirror location"
         x = crown_x(no) - (real(i)-0.5)
         y = crown_y(no) - (real(j)-0.5)
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:rm
!         d1 = (x - real(Max_loc))**2 + (y - real(Max_loc))**2
!         d2 = (x                )**2 + (y - real(Max_loc))**2
!         d3 = (x + real(Max_loc))**2 + (y - real(Max_loc))**2
!         d4 = (x - real(Max_loc))**2 + (y                )**2
!         d5 = (x                )**2 + (y                )**2
!         d6 = (x + real(Max_loc))**2 + (y                )**2
!         d7 = (x - real(Max_loc))**2 + (y + real(Max_loc))**2
!         d8 = (x                )**2 + (y + real(Max_loc))**2
!         d9 = (x + real(Max_loc))**2 + (y + real(Max_loc))**2
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:rm
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
         d1 = (x - real(GRID%Max_x))**2 + (y - real(GRID%Max_y))**2
         d2 = (x                   )**2 + (y - real(GRID%Max_y))**2
         d3 = (x + real(GRID%Max_x))**2 + (y - real(GRID%Max_y))**2
         d4 = (x - real(GRID%Max_x))**2 + (y                   )**2
         d5 = (x                   )**2 + (y                   )**2
         d6 = (x + real(GRID%Max_x))**2 + (y                   )**2
         d7 = (x - real(GRID%Max_x))**2 + (y + real(GRID%Max_y))**2
         d8 = (x                   )**2 + (y + real(GRID%Max_y))**2
         d9 = (x + real(GRID%Max_x))**2 + (y + real(GRID%Max_y))**2
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
         
         !proxy: distance between grid-center and crown-center
         proxy = min(d1,d2,d3,d4,d5,d6,d7,d8,d9)
         proxy = sqrt(proxy)
         
         if ( 0.5*crown_diameter(no) > proxy) then
            flag_coverage = .true.
         endif
         
      end do
      
      if (flag_coverage) count = count+1
      
   End Do
   End Do
   
!   frac_crown_coverage = real(count) / real(Max_loc*Max_loc)   !!!>>>>>>>>>>>>TN:rm
   frac_crown_coverage = real(count) / real(GRID%Area)          !!!<<<<<<<<<<<<TN:add
   
!_____________ Main (Method2)
!   x = 0.0
!   do no=1, Max_no
!   if (.not. tree_exist(no)) cycle
!      x = x + la(no)
!   enddo
!   x = x / Max_loc / Max_loc
!   
!   x = min(5.0,x)/5.0
!   x = x*(2-x)
!   frac_crown_coverage = x
   
END SUBROUTINE crown_coverage



!**************************************************************************************************
! Relative intensity of PAR on each cell of forest floor
!**************************************************************************************************
SUBROUTINE floor_radiation ()

!_____________ Set variables
!Namespace
   USE data_structure
   USE vegi_status_current1
   USE grid_status_current2
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none
   
!Local variables
   real    frac_intercepted, cell_area, cell_length
   real    x, y, dist_x, dist_y, d1, d2, d3, proxy
   integer count, loop_no, loop_no_max
   integer,dimension(Max_no):: tree_no
   integer no, i, j
!   logical,dimension(Dived,Dived)::flag !!!>>>>>>>>>>>>TN:rm
   logical,dimension(GRID%N_x, GRID%N_y)::flag !!!<<<<<<<<<<<<TN:add
   
!_____________ Initialize return variables
   par_grass_rel(:,:) = 1.0
   par_floor_rel(:,:) = 1.0
   
!_____________ count loop number
   loop_no = 0
   Do no = 1, Max_no
      if (.not. tree_exist(no)) cycle
      loop_no          = loop_no + 1
      tree_no(loop_no) = no
   End Do
   loop_no_max = loop_no
   
   !Return if thre are no treess
   IF (loop_no_max == 0) return
   
!_____________ par_grass_rel
!   cell_area   = real(Max_loc**2) / real(DivedG**2) !Area for each grass cell (m2)        !!!>>>>>>>>>>>>TN:rm
!   cell_length = real(Max_loc   ) / real(DivedG   ) !Side length for each grass cell (m)  !!!>>>>>>>>>>>>TN:rm
   cell_area   = 0.25 !Area for each grass cell (m2)         !!!<<<<<<<<<<<<TN:add
   cell_length = 0.5  !Side length for each grass cell (m)   !!!<<<<<<<<<<<<TN:add
   
   DO loop_no = 1, loop_no_max
      no        = tree_no(loop_no)
      count     = 0
      flag(:,:) = .false.
      
!      do i=1, DivedG !!!>>>>>>>>>>>>TN:rm
      do i=1, GRID%N_x !!!<<<<<<<<<<<<TN:add
         !x: location of the center of this grid cells on x-axis
         x = (real(i)-0.5)*cell_length
         
         !dn: distances between centers of crown and grid cell on x-axis
         d1 = crown_x(no) -  x
!         d2 = crown_x(no) - (x - real(Max_loc))!!!>>>>>>>>>>>>TN:rm
!         d3 = crown_x(no) - (x + real(Max_loc))!!!>>>>>>>>>>>>TN:rm
         d2 = crown_x(no) - (x - real(GRID%Max_x))!!!<<<<<<<<<<<<TN:add
         d3 = crown_x(no) - (x + real(GRID%Max_x))!!!<<<<<<<<<<<<TN:add
         
         !dist_x: square of minimum distance between centers of crown and grid cell on x-axis
         dist_x = min(d1**2, d2**2, d3**2)
         
         if ( 0.5*crown_diameter(no) + 0.5*cell_length < sqrt(dist_x) ) cycle
!         do j=1, DivedG !!!>>>>>>>>>>>>TN:rm
         do j=1, GRID%N_y !!!<<<<<<<<<<<<TN:add
            !y: location of the center of this grid cells on j-axis
            y = (real(j)-0.5)*cell_length
            
            !dn: distances between centers of crown and grid cell on y-axis
            d1 = crown_y(no) -  y
!            d2 = crown_y(no) - (y - real(Max_loc))!!!>>>>>>>>>>>>TN:rm
!            d3 = crown_y(no) - (y + real(Max_loc))!!!>>>>>>>>>>>>TN:rm
            d2 = crown_y(no) - (y - real(GRID%Max_y))!!!<<<<<<<<<<<<TN:add
            d3 = crown_y(no) - (y + real(GRID%Max_y))!!!<<<<<<<<<<<<TN:add
            
            !dist_y: square of minimum distance between centers of crown and grid cell on x-axis
            dist_y = min(d1**2, d2**2, d3**2)
            
            !proxy: distance between centers of crown and grid cell on xy-plate
            proxy = sqrt(dist_x + dist_y)
            
            !sumup par_grass_rel, if this cell is under the tree canopy
            if ( 0.5*crown_diameter(no) + 0.5*cell_length < proxy ) cycle
            count     = count+1
            flag(i,j) = .true.
            
         end do
      end do
      if (count==0) cycle
      
      frac_intercepted = exp( -1.0 * EK0(pft(no)) * (la(no)/real(count)/cell_area) )
!      do i=1, DivedG !!!>>>>>>>>>>>>TN:rm
!      do j=1, DivedG !!!>>>>>>>>>>>>TN:rm
      do i=1, GRID%N_x !!!<<<<<<<<<<<<TN:add
      do j=1, GRID%N_y !!!<<<<<<<<<<<<TN:add
         if (.not. flag(i,j)) cycle
         par_grass_rel(i,j) = par_grass_rel(i,j) * frac_intercepted
      end do
      end do
      
   END DO
   
!_____________ Need further computation?
   !When grid layouts of grass cell and establishment patch cell are identical,
   !no need for detailed computation
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:rm
!   IF (Dived==DivedG) then
!      do i=1, Dived
!      do j=i, Dived
!         par_floor_rel(i,j) = par_grass_rel(i,j)
!      enddo
!      enddo
!      return
!   ENDIF
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:rm
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
    do i=1, GRID%N_x
    do j=i, GRID%N_y
       par_floor_rel(i,j) = par_grass_rel(i,j)
    enddo
    enddo
    return
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   
!_____________ par_floor_rel
!   cell_area   = real(Max_loc**2) / real(Dived**2) !Area for each establishment cell (m2)!!!>>>>>>>>>>>>TN:rm
!   cell_length = real(Max_loc   ) / real(Dived   ) !Side length for each establishment cell (m)!!!>>>>>>>>>>>>TN:rm
   cell_area   = 0.25 !Area for each establishment cell (m2)!!!<<<<<<<<<<<<TN:add
   cell_length = 0.5  !Side length for each establishment cell (m)!!!<<<<<<<<<<<<TN:add
   
   DO loop_no = 1, loop_no_max
      no = tree_no(loop_no)
      count = 0
      flag(:,:) = .false.
      
!      do i=1, Dived !!!>>>>>>>>>>>>TN:rm
      do i=1, GRID%N_x !!!<<<<<<<<<<<<TN:add
         !x: location of the center of this grid cells on x-axis
         x = (real(i)-0.5)*cell_length
         
         !dn: distances between centers of crown and grid cell on x-axis
         d1 = crown_x(no) -  x
!         d2 = crown_x(no) - (x - real(Max_loc))!!!>>>>>>>>>>>>TN:rm
!         d3 = crown_x(no) - (x + real(Max_loc))!!!>>>>>>>>>>>>TN:rm
         d2 = crown_x(no) - (x - real(GRID%Max_x))!!!<<<<<<<<<<<<TN:add
         d3 = crown_x(no) - (x + real(GRID%Max_x))!!!<<<<<<<<<<<<TN:add
         
         !dist_x: square of minimum distance between centers of crown and grid cell on x-axis
         dist_x = min(d1**2, d2**2, d3**2)
         
         if ( 0.5*crown_diameter(no) + 0.5*cell_length < sqrt(dist_x) ) cycle
!         do j=1, Dived !!!>>>>>>>>>>>>TN:rm
         do j=1, GRID%N_y !!!<<<<<<<<<<<<TN:add
            !y: location of the center of this grid cells on j-axis
            y = (real(j)-0.5)*cell_length
            
            !dn: distances between centers of crown and grid cell on y-axis
            d1 = crown_y(no) -  y
!            d2 = crown_y(no) - (y - real(Max_loc))!!!>>>>>>>>>>>>TN:rm
!            d3 = crown_y(no) - (y + real(Max_loc))!!!>>>>>>>>>>>>TN:rm
            d2 = crown_y(no) - (y - real(GRID%Max_y))!!!<<<<<<<<<<<<TN:add
            d3 = crown_y(no) - (y + real(GRID%Max_y))!!!<<<<<<<<<<<<TN:add
            
            !dist_y: square of minimum distance between centers of crown and grid cell on x-axis
            dist_y = min(d1**2, d2**2, d3**2)
            
            !proxy: distance between centers of crown and grid cell on xy-plate
            proxy = sqrt(dist_x + dist_y)
            
            !sumup par_grass_rel, if this cell is under the tree canopy
            if ( 0.5*crown_diameter(no) + 0.5*cell_length < proxy ) cycle
            count     = count+1
            flag(i,j) = .true.
            
         end do
      end do
      if (count==0) cycle
      
      frac_intercepted = exp( -1.0 * EK0(pft(no)) * (la(no)/real(count)/cell_area) )
!      do i=1, Dived !!!>>>>>>>>>>>>TN:rm
!      do j=1, Dived !!!>>>>>>>>>>>>TN:rm
      do i=1, GRID%N_x !!!<<<<<<<<<<<<TN:add
      do j=1, GRID%N_y !!!<<<<<<<<<<<<TN:add
         if (.not. flag(i,j)) cycle
         par_floor_rel(i,j) = par_floor_rel(i,j) * frac_intercepted
      end do
      end do
      
   END DO
   
END SUBROUTINE floor_radiation



!*************************************************************************************************
! Calculation of Direct PAR intensity for each crown disk of each tree

! Target variable -> par_direct_rel(individual_id, layer_no in step)
!                    Relative intensity of direct radiation for each crown disk of each tree
!
! Note: y_dimension -> West - East
!       x_dimension -> North - South
!**************************************************************************************************
SUBROUTINE direct_radiation ()

!_____________ Set variables
!Namespace
   USE data_structure
   USE time_counter
   USE vegi_status_current1
   USE grid_status_current2
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none
   
!_____________ Define local parameters
   real    Sl_hgt_adj
   integer Hgt_thres 
   real    Dx        
   
!_____________ Define local variables
   real                        overlap_sum !
   real   ,dimension(Max_hgt)::overlap     !fraction overlap (shading by proximate trees)
   real   ,dimension(Max_no) ::attenu      !light attenuation efficiency
   
   integer,dimension(Max_no       )::save_value1
   integer,dimension(Max_no       )::save_value2
   integer,dimension(Max_no*Max_no)::save_value3
   integer,dimension(Max_no*Max_no)::save_value4
   integer loop_no1 !Required loop number of comuputation: for self-shading
   integer loop_no3 !Required loop number of comuputation: for neigbor-shading
   integer loop_no2 !Required loop number of comuputation: for simple computation
   
   !for geometric calculation
   real    dist, dist_x, dist_y           !distance (m)
   real    r, r1, r2                      !ridius (m)
   real    y1, y2                         !
   real    x1_me,  x2_me,  y1_me,  y2_me  !coodination (m)
   real    x1_you, x2_you, y1_you, y2_you !coodination (m)
   real    y_shift, xband_me, xband_you   !coodination (m)
   real    cosine, cosine1, cosine2       !cosine
   
   !for general usage
   integer no, me, you, loop               !loop counter
   integer l1, l2                          
   integer i, j                            
   real    x, y, a1, a2, a3, c1, c2
   integer flag                            
   
!_____________ Define local parameters
   !Angle of virtual pipe (unit=degree)
   !(this angle equally divides daily sum of sun radiation into upper and lower parts)
   Sl_hgt_adj = 0.86* sl_hgt(doy)
   
   !Threshold height, below which simple computation will be applied for par_direct_rel
   !(unit=STEP, range=0~Max_hgt)
   Hgt_thres  = 0
   
   !Shade lengh for single height step (unit=m/STEP)
   Dx = STEP / tan(Sl_hgt_adj*DtoR)
   
!_____________ Following procedures required only if sunlight exist
   IF (Sl_hgt_adj <= 10.0) return
   
!_____________ Initialize variables
   loop_no1 = 0
   loop_no2 = 0
   loop_no3 = 0
   par_direct_rel(:,:) = 1.00000 !(individual_no, layer_no in step), default values
   
!_____________ Select trees that require computations
   DO no = 1, Max_no
   IF ( .not. tree_exist(no)      ) cycle
   IF ( .not. phenology(pft(no))  ) cycle
   IF ( crown_diameter(no) == 0.0 ) cycle
   IF ( la(no) == 0.0             ) cycle
      
      !Trees that require intensive computaiton
      if ( height(no) >= Hgt_thres+1 ) then
         loop_no1              = loop_no1 + 1 !tree number for intensive computation
         save_value1(loop_no1) = no           !save tree ID for each loop number
         !x: Leaf Area Density for each tree (m2/m2/STEP)
         x = la(no) / max(0.001,crown_area(no)) / real( max(1,height(no)-bole(no)) )
         !Light attenuation coefficient for each tree (EK0 should be used, not eK)
         attenu(no) = -1.000 * EK0(pft(no)) * x
      endif
      
      !Trees that require simple computaiton
      if ( bole(no)+1 <= Hgt_thres ) then
         loop_no2              = loop_no2 + 1 !tree number for simple computation
         save_value2(loop_no2) = no           !save tree ID for each loop number
      endif
      
   END DO
   
!_____________ Simple computation for light attnuation
   If (loop_no2 >= 1) then
   Do loop = 1, loop_no2      !For each tree
      me = save_value2(loop)  !Recall ID number of the current tree
      
      do i = bole(me)+1, min(height(me),Hgt_thres)
         par_direct_rel(me,i) = par_diffuse_rel(i)
      end do
      
   End Do
   End If
   
!_____________ Intensive computation of self-shading
   IF (loop_no1 == 0) return
   DO loop = 1, loop_no1             !For each tree
      !Prepare tree specific variables
      me      = save_value1(loop)               !ID number of the current tree
      r       = crown_diameter(me) / 2.0        !crown radius of the current tree (m)
      
      !cancel trees without computation requirement
      if ( height(me)-bole(me) < 2   ) cycle !less than 2 crown layers exist
      if ( height(me) == Hgt_thres+1 ) cycle !less than 2 crown layers exist above Hgt_thres
      
      !Reset sum of shade fraction between disks
      overlap_sum = 0.0
      
      !Compute self shading for each crown layer (iは相対高度、jは絶対高度)
      Do i = 1, height(me)-max(bole(me),Hgt_thres)-1 !substration of 1 is due to top layer
         !j:height of mensioned crown disk (in STEP)
         j = height(me) - i
         
         !dist: x-axis distance between shade and current disk
         dist = Dx * real(i)
!         dist = min(dist, abs(dist-Max_loc))!!!>>>>>>>>>>>>TN:rm
         dist = min(dist, abs(dist-real(GRID%Max_x)))!!!<<<<<<<<<<<<TN:add
         
         !sumup fraction of overlap area
         if ( dist < 2.0*r ) then
            cosine = dist / (2.0*r)
            cosine = min(1.0,max(-1.0,cosine))
            overlap_sum = overlap_sum + (2/PI)*( acos(cosine) - cosine * sqrt(1.0 - cosine**2) )
         endif
         
         !Give fraction available light from top layer
         par_direct_rel(me, j) = exp( attenu(me) * overlap_sum )
         
      End do
      
   END DO
   
!_____________ List up tree-pairs that interfare
   IF (loop_no1 < 2) return
   Do i = 1, loop_no1 !tree 'me' 
   Do j = 1, loop_no1 !tree 'you'
   
   !Recall tree ID
      me  = save_value1(i)
      you = save_value1(j)
      
   !Cancel no-interferebce trees 1
      if ( bole(me)+1 >= height(you) ) cycle
      if ( i == j )                    cycle
      
   !Cancel no-interferebce trees 2
      !check whether crowns of 'me' and 'you' share some range of y axis
      a1 = abs( crown_y(me) - crown_y(you)           )
!      a2 = abs( crown_y(me) - crown_y(you) + Max_loc )!!!>>>>>>>>>>>>TN:rm
!      a3 = abs( crown_y(me) - crown_y(you) - Max_loc )!!!>>>>>>>>>>>>TN:rm
      a2 = abs( crown_y(me) - crown_y(you) + real(GRID%Max_y) )!!!<<<<<<<<<<<<TN:add
      a3 = abs( crown_y(me) - crown_y(you) - real(GRID%Max_y) )!!!<<<<<<<<<<<<TN:add
      if (min(a1,a2,a3) > (crown_diameter(me) + crown_diameter(you))/2 ) cycle
      
   !Cancel no-interferebce trees 3
      !check whether crowns 'me' interferes with crown 'you' on z-x plane
      !xband_me: x-axis range of crown 'me' that can interfere with crown 'you'
      y1_me  = crown_y(me)  + crown_diameter(me)  / 2.0
      y2_me  = crown_y(me)  - crown_diameter(me)  / 2.0
      y1_you = crown_y(you) + crown_diameter(you) / 2.0
      y2_you = crown_y(you) - crown_diameter(you) / 2.0
      
      !y_shift: y-axis shift length of crown 'you' for adjusting repeat world
      a1 = abs( y2_me - (y2_you        ) )
!      a2 = abs( y2_me - (y2_you+Max_loc) )!!!>>>>>>>>>>>>TN:rm
!      a3 = abs( y2_me - (y2_you-Max_loc) )!!!>>>>>>>>>>>>TN:rm
      a2 = abs( y2_me - (y2_you+real(GRID%Max_y)) )!!!<<<<<<<<<<<<TN:add
      a3 = abs( y2_me - (y2_you-real(GRID%Max_y)) )!!!<<<<<<<<<<<<TN:add
      
      y_shift = 0.0
!      if ( a2<=a1 .and. a2<=a3 ) y_shift =        real(Max_loc)!!!>>>>>>>>>>>>TN:rm
!      if ( a3<=a1 .and. a3<=a2 ) y_shift = -1.0 * real(Max_loc)!!!>>>>>>>>>>>>TN:rm
      if ( a2<=a1 .and. a2<=a3 ) y_shift =        real(GRID%Max_y)!!!<<<<<<<<<<<<TN:add
      if ( a3<=a1 .and. a3<=a2 ) y_shift = -1.0 * real(GRID%Max_y)!!!<<<<<<<<<<<<TN:add
      
      !y1: larger  y-location of interfering range of crowns
      !y2: smaller y-location of interfering range of crowns
      y1 = min(y1_me, y1_you + y_shift)
      y2 = max(y2_me, y2_you + y_shift)
      
      if ( (crown_y(me)<y1) .and. (crown_y(me)>y2) ) then
         xband_me = crown_diameter(me)
      else
         x = 0.25 * ( crown_diameter(me)**2 )
         y = min( abs(crown_y(me)-y1) , abs(crown_y(me)-y2) ); y = y**2
         xband_me = 2.0 * sqrt(abs(x-y))
      endif
      
      !xband_you: x-axis range of crown 'you' that can interfere with crown 'me'
      if ( (crown_y(you)+y_shift<y1) .and. (crown_y(you)+y_shift>y2) ) then
         xband_you = crown_diameter(you)
      else
         x = 0.25 * ( crown_diameter(you)**2 )
         y = min( abs(crown_y(you)+y_shift-y1) , abs(crown_y(you)+y_shift-y2) ); y = y**2
         xband_you = 2.0 * sqrt(abs(x-y))
      endif
      
      !Set coordination
      x1_me  = crown_x(me)  + xband_me  / 2.0 !left  edge of crown me  on x-axixs
      x2_me  = crown_x(me)  - xband_me  / 2.0 !right edge of crown me  on x-axixs
      x1_you = crown_x(you) + xband_you / 2.0 !left  edge of crown you on x-axixs
      x2_you = crown_x(you) - xband_you / 2.0 !right edge of crown you on x-axixs
      
      flag = 1
      if (Sl_hgt_adj<90.0) then
         c1 = x1_me + Dx * real( height(you) - max(bole(me)+1, Hgt_thres+1) )
         c2 = x2_me + Dx * real( max( height(me) , bole(you)+1 ) - height(me) )
         if ( c1>x2_you         .and. c2<x1_you         ) flag= 0
!         if ( c1>x2_you+Max_loc .and. c2<x1_you+Max_loc ) flag= 0!!!>>>>>>>>>>>>TN:rm
         if ( c1>x2_you+real(GRID%Max_x) .and. c2<x1_you+real(GRID%Max_x) ) flag= 0!!!<<<<<<<<<<<<TN:add
      else
         c1 = x2_me - Dx * real( height(you) - max(bole(me)+1, Hgt_thres+1) )
         c2 = x1_me - Dx * real( max( height(me) , bole(you)+1 ) - height(me) )
         if ( c1<x1_you         .and. c2>x2_you         ) flag= 0
!         if ( c1<x1_you-Max_loc .and. c2>x2_you-Max_loc ) flag= 0!!!>>>>>>>>>>>>TN:rm
         if ( c1<x1_you-real(GRID%Max_x) .and. c2>x2_you-real(GRID%Max_x) ) flag= 0!!!<<<<<<<<<<<<TN:add
      endif
      if (flag == 1) cycle
      
   !Save values
      loop_no3              = loop_no3 + 1
      save_value3(loop_no3) = me          
      save_value4(loop_no3) = you         
      
   End do
   End do
   
!_____________ Compute among trees shading
!CDIR LOOPCNT = 15000
DO loop = 1, loop_no3 !for each combination of individuals
   
   !Reset shade fraction between disks
   overlap(:) = 0.0
   
   !Recall variables for this comparison
   me     = save_value3(loop) !ID of current mentioned tree
   you    = save_value4(loop) !ID of interfering tree
   
   !dist_y: absolute y-axis distance between circle centers (adjusted for mirror world)
   dist_y = crown_y(me)-crown_y(you)
!   dist_y = min( abs(dist_y), abs(dist_y-Max_loc), abs(dist_y+Max_loc) )!!!>>>>>>>>>>>>TN:rm
   dist_y = min( abs(dist_y), abs(dist_y-real(GRID%Max_y)), abs(dist_y+real(GRID%Max_y)) )!!!<<<<<<<<<<<<TN:add
   
   !Compute fraction of shading from above crown disk
   do i = max(1,bole(you)+1-height(me)) , height(you) - max(bole(me)+1,Hgt_thres+1)
      !i: Height difference between crown disks in STEP, Relative value
      
      !dist_x: absolute x-axis distance between circle centers (adjusted for mirror world)
      if (Sl_hgt_adj<=90.0) then
         dist_x = crown_x(me) - (crown_x(you)-Dx*real(i))
      else
         dist_x = crown_x(me) - (crown_x(you)+Dx*real(i))
      endif
!      dist_x = min( abs(dist_x), abs(dist_x-Max_loc), abs(dist_x+Max_loc) )!!!>>>>>>>>>>>>TN:rm
      dist_x = min( abs(dist_x), abs(dist_x-real(GRID%Max_x)), abs(dist_x+real(GRID%Max_x)) )!!!<<<<<<<<<<<<TN:add
      
      !dist: x-y-plane distance between crown and shade center
      dist = sqrt(dist_x**2 + dist_y**2)
      
      !x: overlapping area (in m^2)
      r1 = min(crown_diameter(me),crown_diameter(you)) / 2.0 !radius of smaller circle (m)
      r2 = max(crown_diameter(me),crown_diameter(you)) / 2.0 !radius of larger  circle (m)
      if     (r1+r2 <= dist) then
         x = 0.000
         
      elseif (dist + r1 <= r2) then
         x = PI * r1 * r1
         
      elseif (r2 < dist ) then
         cosine1 = (r1**2 - r2**2 + dist**2) / (2.0*r1*dist); cosine1 = min(1.0,max(-1.0,cosine1))
         cosine2 = (r2**2 - r1**2 + dist**2) / (2.0*r2*dist); cosine2 = min(1.0,max(-1.0,cosine2))
         x       = r1*r1*( acos(cosine1) - cosine1 * sqrt(1.0 - cosine1**2) ) &
                 + r2*r2*( acos(cosine2) - cosine2 * sqrt(1.0 - cosine2**2) )  
         
      else
         cosine1 = (r1**2 - r2**2 + dist**2) / (2.0*r1*dist); cosine1 = min(1.0,max(-1.0,cosine1))
         cosine2 = (r2**2 - r1**2 + dist**2) / (2.0*r2*dist); cosine2 = min(1.0,max(-1.0,cosine2))
         x       = r1*r1*acos(cosine1) + r2*r2*acos(cosine2)   &
                 - dist * r2 * sqrt(1.0 - cosine2**2)  
         
      end if
      overlap(i) = min(1.000, x / crown_area(me)) !convert to fraction overlap
   end do
   
   !reflect shadeing to par_direct_rel for each crown disk
   do i = max(Hgt_thres+1, bole(me)+1), min(height(me), height(you)-1)
      !i: height of current mentioned crown, Absolute value
      l1 = max(1,   bole(you) + 1 - i) !lower, Relatige value
      l2 = max(1, height(you)     - i) !upper, Relatige value
      par_direct_rel(me,i) = par_direct_rel(me,i) * exp( attenu(you) * sum(overlap(l1:l2)) )
   end do
   
END DO

END SUBROUTINE direct_radiation



!**************************************************************************************************
! obtain ground vacant where new sapling can establish
!**************************************************************************************************
SUBROUTINE ground_vacant ()

!_____________ Set variables
!Namespace
   USE data_structure
   USE vegi_status_current1
   USE vegi_status_current2
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none
   
!Local parameters
   real   ,parameter::Crown_diameter_new = 1.0 !maximum crown diameter of new tree (m)
   integer,parameter::Height_new         = 20  !maximum height of new tree      (step)
   
!Local variables
   real    distance                           !distance (m)
   real    d1, d2, d3, d4, d5, d6, d7, d8, d9 !distance (m)
   real    x, y                               !for temporal usage
   integer no, i, j                           !for loop counter
   
!_____________ Main part
!Reset output variables
!   patch_vacant(:,:) = .true.!!!>>>>>>>>>>>>TN:rm
   patch_vacant(:,:) = GRID%mask(:,:) !!!<<<<<<<<<<<<TN:add
   
!Calculate safe-site in relation to proximity of previous trees
Do no=1, Max_no
if ( .not. tree_exist(no)   ) cycle
   
   !two trees cannot exist at the same establishment location
!   x = real(Dived)/real(Max_loc)!!!>>>>>>>>>>>>TN:rm
   x = 2.0 !!!<<<<<<<<<<<<TN:add
   patch_vacant( int(bole_x(no)*x)+1 , int(bole_y(no)*x)+1 ) = .false.
   
if ( bole(no) >= Height_new ) cycle
   
   !examine for proximate area for each grid of establishment coodinate
!   do i = 1, Dived !!!>>>>>>>>>>>>TN:rm
!   do j = 1, Dived !!!>>>>>>>>>>>>TN:rm
   do i = 1, GRID%N_x !!!<<<<<<<<<<<<TN:add
   do j = 1, GRID%N_y !!!<<<<<<<<<<<<TN:add
      
      !location of this tree under nomal coodinate
!      x = real(i) * real(Max_loc) / real(Dived) - 0.5 !x location (m)!!!>>>>>>>>>>>>TN:rm
!      y = real(j) * real(Max_loc) / real(Dived) - 0.5 !y location (m)!!!>>>>>>>>>>>>TN:rm
      x = (real(i) - 0.5) * 0.5  !x location (m)!!!<<<<<<<<<<<<TN:add?????bug?????要確認
      y = (real(j) - 0.5) * 0.5  !y location (m)!!!<<<<<<<<<<<<TN:add?????bug?????要確認
      
      !calculate distance between my place and target place
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:rm
!      d1 = (crown_x(no)-Max_loc-x)**2 + (crown_y(no)-Max_loc-y)**2
!      d2 = (crown_x(no)        -x)**2 + (crown_y(no)-Max_loc-y)**2
!      d3 = (crown_x(no)+Max_loc-x)**2 + (crown_y(no)-Max_loc-y)**2
!      d4 = (crown_x(no)-Max_loc-x)**2 + (crown_y(no)        -y)**2
!      d5 = (crown_x(no)        -x)**2 + (crown_y(no)        -y)**2
!      d6 = (crown_x(no)+Max_loc-x)**2 + (crown_y(no)        -y)**2
!      d7 = (crown_x(no)-Max_loc-x)**2 + (crown_y(no)+Max_loc-y)**2
!      d8 = (crown_x(no)        -x)**2 + (crown_y(no)+Max_loc-y)**2
!      d9 = (crown_x(no)+Max_loc-x)**2 + (crown_y(no)+Max_loc-y)**2
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:rm
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
      d1 = (crown_x(no)-real(GRID%Max_x)-x)**2 + (crown_y(no)-real(GRID%Max_y)-y)**2
      d2 = (crown_x(no)                 -x)**2 + (crown_y(no)-real(GRID%Max_y)-y)**2
      d3 = (crown_x(no)+real(GRID%Max_x)-x)**2 + (crown_y(no)-real(GRID%Max_y)-y)**2
      d4 = (crown_x(no)-real(GRID%Max_x)-x)**2 + (crown_y(no)                 -y)**2
      d5 = (crown_x(no)                 -x)**2 + (crown_y(no)                 -y)**2
      d6 = (crown_x(no)+real(GRID%Max_x)-x)**2 + (crown_y(no)                 -y)**2
      d7 = (crown_x(no)-real(GRID%Max_x)-x)**2 + (crown_y(no)+real(GRID%Max_y)-y)**2
      d8 = (crown_x(no)                 -x)**2 + (crown_y(no)+real(GRID%Max_y)-y)**2
      d9 = (crown_x(no)+real(GRID%Max_x)-x)**2 + (crown_y(no)+real(GRID%Max_y)-y)**2
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
      
      distance = min(d1,d2,d3,d4,d5,d6,d7,d8,d9)
      distance = sqrt(distance)
      
      !determine whether current site is safe or not
      if ( distance < 0.5 * (crown_diameter(no)+Crown_diameter_new) ) then
         patch_vacant(i,j) = .false.
      endif
   end do
   end do
   
End do
   
END SUBROUTINE ground_vacant



!**************************************************************************************************
! Calculate spacial limitation on growth (for woody PFTs, Monthely computation)
!**************************************************************************************************
SUBROUTINE spatial_limitation ()

!_____________ Set variables
!Namespace
   USE data_structure
   USE vegi_status_current1
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none
   
!Local parameters
   real,parameter::Margine=0.1  !minimum proximate distance between crowns (m)
   

!Local variable
   real    proxy, x, y
   integer me, you
   logical flag1, flag2
   
!_____________ for each combination of trees
DO me=1, Max_no
if ( .not. tree_exist(me)     ) cycle !case no tree
if ( .not. phenology(pft(me)) ) cycle !case dormant phase
   
   !Give potential maximum
   height_limit(me) = Max_hgt
   radius_limit(me) = CD_max(pft(me)) / 2.0
   
   !Examine limitations of proxymate trees
   DO you=1, Max_no
   if ( .not. tree_exist(you)     ) cycle !case no tree
   if ( me         == you         ) cycle !case same tree
   if ( bole  (me) >= height(you) ) cycle !case no interference
      
      !Obtain most proximate distance among 9 duplicated-target-trees
      x = crown_x(you) - crown_x(me) !difference in x-location
!      x = min(x**2, (x+real(Max_loc))**2 , (x-real(Max_loc))**2 )!!!>>>>>>>>>>>>TN:rm
      x = min(x**2, (x+real(GRID%Max_x))**2 , (x-real(GRID%Max_x))**2 )!!!<<<<<<<<<<<<TN:add
      
      y = crown_y(you) - crown_y(me) !difference in x-location
!      y = min(y**2, (y+real(Max_loc))**2 , (y-real(Max_loc))**2 )!!!>>>>>>>>>>>>TN:rm
      y = min(y**2, (y+real(GRID%Max_y))**2 , (y-real(GRID%Max_y))**2 )!!!<<<<<<<<<<<<TN:add
      
      proxy = sqrt(x+y)
      
      !Give height or radius growth limit
      flag1 = ( bole(you) >= height(me) )
      flag2 = ( proxy < crown_diameter(me)/2 + crown_diameter(you)/2 + Margine )
      
      if ( flag1 ) then
         !Height growth limit
         if ( flag2 ) then
         height_limit(me) = min( height_limit(me), bole(you) )
         endif
      else
         !Radius growth limit
         x = proxy - crown_diameter(you)/2 - Margine
         radius_limit(me) = min( radius_limit(me), x )
      endif
      
   END DO
   
END DO

END SUBROUTINE spatial_limitation


!**************************************************************************************************
! Crown movement
!**************************************************************************************************
SUBROUTINE crown_shake ()

!_____________ Set variables
!Namespace
   USE data_structure
   USE vegi_status_current1
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE mod_grid
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none
   
!Local parameters
   integer,parameter::Direct_resolusion = 36  !resolution of movement direction   (n)
   real   ,parameter::Crown_plasticity  = 0.5 !maximum fraction of crown movement (fraction)
   real   ,parameter::Crown_max_move    = 0.2 !maximum rate of crown movement     (m/year)
   real   ,parameter::Search_radius     = 3.0 ! (m)
   integer,parameter::Max_proxy         = 100 !maximum number of proximate trees to examine (n)
   
!Local variable
   !maximum lengh of crown movement for each direction (m)
   real,dimension(Direct_resolusion)::lenght_max1 !by constraint due to Crown_plasticity
   real,dimension(Direct_resolusion)::lenght_max2 !by constraint due to surrounding trees
   
   !other varibles
   real,dimension(Max_proxy, 3)     ::proxy_crown !(no_you, your r/crown_x/crown_y)
   
   !other variables
   real    angle, angle0, angle1, angle2 !angle (radius)
   real    x, y, z                       !for temporal usage
   real    dx, dy, r, d                  !for temporal usage
   real    a1, a2, a3                    !for temporal usage
   real    distance                      !
   integer direction                     !
   integer count1, count2                !for counter
   integer me, you, i                    !for temporal counter
   
!_____________ Start of loop
DO me=1, Max_no
IF ( .not. tree_exist(me) )  cycle
   
!_____________ Determine maximum plasticity of crown movement by considering self-tree constraint
   !output --> lenght_max1 (Direct_resolusion)
   
   !deviations between bole_center and crown_center (m)
   dx = crown_x(me) - bole_x(me) !on x-axis
   dy = crown_y(me) - bole_y(me) !on y-axis
   d  = sqrt( dx**2 + dy**2 )    !on x-y-plane
   
   !maximum deviations between bole_center and crown_center on x-y-plane (m)
   r = 0.5 * crown_diameter(me) * Crown_plasticity
   
   IF (d==0.0) then !when bole_center == crown_center
      
      lenght_max1(:) = r
      
   ELSE
      angle0 = acos( -1.0 * dy / d )
      
      !for each direction of the movement
      DO i=1, Direct_resolusion
         !obtain angle2
         angle1 = 2.0 * PI * ( real(i) / real(Direct_resolusion) )
         
         if ( dx >= 0.0 ) then
            angle2 = angle0 + angle1
         else
            angle2 = abs(angle0 - angle1)
         endif
         
         do while (angle2 > 2.0*PI); angle2 = angle2 - 2.0*PI; enddo
         
         if ( angle2 > PI ) angle2 = 2.0*PI - angle2
         
         !give maximum movement range
         z = d*d*cos(angle2)*cos(angle2) + r*r - d*d
         if (z < 0.0) then
            lenght_max1(i) = 0.0
         else
            lenght_max1(i) = d*cos(angle2) + sqrt(z)
         endif
         
      END DO
      
   ENDIF
   
!_____________ Find crowns that may cause interference
!   output --> proxy_crown (Max_proxy, 3), count1
   
   !initialize
   count1 = 0
   proxy_crown(:,:) = 0.0
   
   !for each tree
   DO you = 1, Max_no
   if ( count1==Max_proxy         ) exit                       
   if ( .not. tree_exist(you)     ) cycle !case tree not exist 
   if ( bole  (me) >  height(you) ) cycle !case no interference
   if ( height(me) <= bole  (you) ) cycle !case no interference
   if ( me         == you         ) cycle !case same tree      
      
      !common variables in the following loop
      a1 = crown_x(you) - crown_x(me)
      a2 = crown_y(you) - crown_y(me)
      a3 = (crown_diameter(me)+crown_diameter(you))/2.0 + max(Search_radius, Crown_max_move)
      
      DO i=1, 9
         
         !dx, dy: deviations of crown centers on x and y axis, respectively
         select case (i)
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:rm
!           case(1); dx = a1 - Max_loc ; dy = a2 - Max_loc
!           case(2); dx = a1           ; dy = a2 - Max_loc
!           case(3); dx = a1 + Max_loc ; dy = a2 - Max_loc
!           case(4); dx = a1 - Max_loc ; dy = a2          
!           case(5); dx = a1           ; dy = a2          
!           case(6); dx = a1 + Max_loc ; dy = a2          
!           case(7); dx = a1 - Max_loc ; dy = a2 + Max_loc
!           case(8); dx = a1           ; dy = a2 + Max_loc
!           case(9); dx = a1 + Max_loc ; dy = a2 + Max_loc
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:rm
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
           case(1); dx = a1 - real(GRID%Max_x) ; dy = a2 - real(GRID%Max_y)
           case(2); dx = a1                    ; dy = a2 - real(GRID%Max_y)
           case(3); dx = a1 + real(GRID%Max_x) ; dy = a2 - real(GRID%Max_y)
           case(4); dx = a1 - real(GRID%Max_x) ; dy = a2                   
           case(5); dx = a1                    ; dy = a2                   
           case(6); dx = a1 + real(GRID%Max_x) ; dy = a2                   
           case(7); dx = a1 - real(GRID%Max_x) ; dy = a2 + real(GRID%Max_y)
           case(8); dx = a1                    ; dy = a2 + real(GRID%Max_y)
           case(9); dx = a1 + real(GRID%Max_x) ; dy = a2 + real(GRID%Max_y)
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
         end select
         
         !list up for all interference crown
         if ( count1==Max_proxy       ) exit
         if ( sqrt(dx*dx+dy*dy) >= a3 ) cycle
         if ( dx==0.0 .and. dy==0.0   ) cycle
         
         !in this cas, this tree can interfare, and thus its information is saved
         count1 = count1 + 1
         proxy_crown(count1, 1) = 0.5*crown_diameter(you) !radius of the target tree (m)
         proxy_crown(count1, 2) = dx !deviation of crown centers on x-axis (m)
         proxy_crown(count1, 3) = dy !deviation of crown centers on y-axis (m)
         
      END DO
   END DO
   
!_____________ Adjust maximum movement of crown center by proxymate trees
!   output --> lenght_max2 (Direct_resolusion)
   
   !initialize
   lenght_max2(:) = Search_radius
   
   !for each proxymate tree
   DO count2=1, count1
      
      r        = proxy_crown(count2, 1) !radius (m)
      dx       = proxy_crown(count2, 2) !relative x location (m)
      dy       = proxy_crown(count2, 3) !relative y location (m)
      distance = sqrt( dx*dx + dy*dy )
      
      if ( dx >= 0.0 ) then
         angle0 =          acos( dy / distance )
      else
         angle0 = 2.0*PI - acos( dy / distance )
      endif
      
      !for each direction
      DO i=1, Direct_resolusion
         angle1 = 2.0 * PI * ( real(i) / real(Direct_resolusion) )
         angle2 = angle0 - angle1
         
         !covert target crown location on the new coodinate system
         x = distance * sin(angle2) !converted x (m)
         y = distance * cos(angle2) !converted y (m)
         
         !select interference situation
         if ( y <= 0.0                       ) cycle
         if (  0.5*crown_diameter(me) <= x-r ) cycle
         if ( -0.5*crown_diameter(me) >= x+r ) cycle
         
         !adjust maximum movement, z (m)
         z = y - sqrt( (0.5*crown_diameter(me)+r)**2 - x**2 )
         z = max(0.0, z)
         lenght_max2(i) = min( lenght_max2(i), z )
         
      END DO
   END DO
   
!_____________ Execute crown movement
!   output --> crown_x(me), crown_y(me)
   !In case of no interference, following procedure will be omitted
   IF ( minval(lenght_max2(:)) .ne. maxval(lenght_max2(:)) ) then
      
      !determine direction and distance of crown movement
      direction = 1
      distance  = 0.0
      DO i=1, Direct_resolusion  !for each direction
         if ( lenght_max2(i)>distance .and. lenght_max1(i)>0.0 ) then
            direction = i
            distance  = lenght_max2(i)
         endif
      END DO
      angle = 2.0 * PI * ( real(direction) / real(Direct_resolusion) )
      distance  = min( lenght_max2(direction), lenght_max1(direction), Crown_max_move )
      
      !execute crown movement
      crown_x(me) = crown_x(me) + distance * sin(angle)
      crown_x(me) = max(crown_x(me),           0.0)
!      crown_x(me) = min(crown_x(me), real(Max_loc))!!!>>>>>>>>>>>>TN:rm
      crown_x(me) = min(crown_x(me), real(GRID%Max_x))!!!<<<<<<<<<<<<TN:add
      
      crown_y(me) = crown_y(me) + distance * cos(angle)
      crown_y(me) = max(crown_y(me),           0.0)
!      crown_y(me) = min(crown_y(me), real(Max_loc))!!!>>>>>>>>>>>>TN:rm
      crown_y(me) = min(crown_y(me), real(GRID%Max_y))!!!<<<<<<<<<<<<TN:add
      
   ENDIF
   
!_____________ End of loop
END DO

END SUBROUTINE crown_shake

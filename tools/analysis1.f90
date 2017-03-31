!**************************************************************************
! Tool for SEIB-DGVM, 1 : Forest demography analysis
! (C)2005 Hisashi SATO (FRCGC/JAMSTEC)
! 
! Input  -> Forest structure data outputted by SEIB-DGVM
! Output -> Frequency, Growth_rate, and Mortality for each tree size class
!           (size class is defined by Diameter at Breast Height)
!**************************************************************************

PROGRAM analysis1

implicit none

!Set simulation properties (should be same as 'parameter.txt')
   real   ,parameter::Max_loc         = 30.0 !width lenght of forest stand (m)
   integer,parameter::Max_no          = 900  !maximum individual number in a plot
!   integer,parameter::Simulation_year = 1000  !simulation year (yr)

!Set I/O file names
   character*128,parameter::Fn1='output_forest.txt' !input data name
   character*128,parameter::Fn2='analysis1.txt'     !output file name

!Set parameters: Input file definition
   integer,parameter::Interval1 = 5  !comparison length for mortality calculation (year)
   integer,parameter::Interval2 = 50 !interval for analysis repeat (year)

!Set parameters: Analysis definition
   integer,parameter            ::Sizeclass_num = 10 !maximum number of size class
   real,dimension(Sizeclass_num)::Size_set           !sizeclass structure (m)
   real,dimension(Sizeclass_num,3)::output !tree density & mortality for each sizeclass
    data Size_set / 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50 /
!   data Size_set / 0.02, 0.04, 0.06, 0.08, 0.10, 0.12 /
!   data Size_set / 0.05, 0.10, 0.20, 0.30, 0.40, 0.50 /

!Set variables
   real   ,dimension(1000, Max_no, 9)::data_input !inpur data
   integer,dimension(1000)           ::tree_num  !tree number for each year
   
   integer  Repeat
   integer  now, prev
   integer  size
   real     dbh
   integer  i, j, k, r, year, no


   integer   Simulation_year  !simulation year (yr)
   read(*,*) Simulation_year

!Initiate valiables
   tree_num(:)       = 0  
   data_input(:,:,:) = 0.0
   output(:,:)       = 0.0

!Set repeat number of analysis
   Repeat = int(real(Simulation_year) / real(Interval2)) !number of comparison repeat

!Omit this program?
   !case1
   if ( Repeat < 1 ) then
      write(*,*) 'Aborted!: data record is not enough'
      stop
   end if
   
   !case2
   if ( Interval1 > Interval2 ) then
      write(*,*) 'Aborted!: invalid parameter definition'
      stop
   end if

!Read data file
   open (1, file = Fn1)
   do year=1, Simulation_year
      
      !tree_num(year)        -> number of tree for each year
      read (1,*) tree_num(year)
      
      !data_input(year,no,:) -> each individual data
      do no=1, tree_num(year) ; read (1,*) data_input(year,no,:) ; end do
      
   end do
   close(1)

!Analysis
DO r=1, Repeat !For each repeat comparison
   prev = r*Interval2-Interval1 !year 0: comparison base year
   now  = r*Interval2           !year 1: current year
   
DO no=1, tree_num(prev)
   
   !determine size class of trees in previous year
   dbh  = data_input(prev,no,7) !diameter at breast height (m)
   size = 0                     !reset size class
   do i=1, Sizeclass_num
      if ( dbh >= Size_set(i) ) size = size + 1
   enddo
   
   !if size_class is below the minimum threshold, omit this individual
   if (size==0) cycle
   
   !count tree number for each size_class in the 'previous' year
   output(size,1) = output(size,1) + 1
   
   !comparison between previous and current years
   do j=1, tree_num(now)
      !screen alived tree in the current year
      if ( data_input(prev,no,1) .ne. data_input(now,j,1) ) cycle !loc_x
      if ( data_input(prev,no,2) .ne. data_input(now,j,2) ) cycle !loc_y
      if ( data_input(prev,no,5) >    data_input(now,j,5) ) cycle !bole height
      if ( data_input(prev,no,7) >    data_input(now,j,7) ) cycle !dbh
      if ( data_input(prev,no,9) .ne. data_input(now,j,9) ) cycle !pft
      
      !add total increament of DBH of survived trees (mm)
      output(size,2) = output(size,2) + 1000.0 * (data_input(now,j,7) - dbh)
      
      !add number of survived trees
      output(size,3) = output(size,3) + 1
      
      exit
   end do
END DO
END DO
   
!Process output variables
   do i=1, Sizeclass_num
   if ( output(i,1) > 0.0 ) then
      !Growth rate: (total increament of DBH) / (survived tree tree)
      output(i,2) = output(i,2) / output(i,3) / real(Interval1)
      
      !Death rate per year: 1 - (survived tree) / (previous tree)
      output(i,3) = ( 1.00-real(output(i,3))/real(output(i,1)) ) / real(Interval1)
      
      !tree density of previous year
      output(i,1) = real(output(i,1)) / (Max_loc*Max_loc*0.0001) / real(Repeat)
   endif
   enddo

!Write on file
open (1, file=Fn2)
   !write setting parameters
   write(1,*) 'comparison length for mortality calculation (year)', Interval1
   write(1,*) 
   write(1,*) 'interval year for analysis repeart (year)         ', Interval2
   write(1,*) 
   write(1,*) 'sizeclass setting (m)'
   do i = 1, Sizeclass_num; write(1,'(f7.3)') Size_set(i); enddo
   write(1,*) 
   
   !write population density
   write(1,*) 'RESULT1: population density (ha-1)'
   do i = 1, Sizeclass_num; write(1,'(f7.2)') output(i,1); enddo
   write(1,*)
   
   !write growth rate
   write(1,*) 'RESULT2: DBH growth rate (mm year-1)'
   do i = 1, Sizeclass_num; write(1,'(f7.2)') output(i,2); enddo
   write(1,*)
   
   !write mortality
   write(1,*) 'RESULT3: mortality (year-1)'
   do i = 1, Sizeclass_num; write(1,'(f7.3)') output(i,3); enddo
close(1)

END program analysis1

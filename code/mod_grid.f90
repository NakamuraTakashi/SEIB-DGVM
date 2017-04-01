
!!!=== ver 2017/03/31   Copyright (c) 2017 Takashi NAKAMURA  =====

!!!**** MODULE OF NETCDF READ & WRITE ****************

  module mod_grid
  
    USE netcdf
    
    implicit none

    TYPE T_GRID
    
      integer, pointer :: N_x, N_y       ! Number of cells == Dived
      integer, pointer :: Max_x, Max_y   ! Lenght of computational domain == Max_loc
      integer, pointer :: Area           ! Area of domain
      integer, pointer :: N_tot          ! Total number of cells
      real,    pointer :: h(:,:)         ! depth (m)
      logical, pointer :: mask(:,:)      ! Land mask: land -> .false., wet -> .true.
      real,    pointer :: sal_ave(:,:)   ! Mean salinity (psu)
      real,    pointer :: sal_max(:,:)   ! Maximum salinity (psu)
      real,    pointer :: sal_min(:,:)   ! Minimum salinity (psu)

    END TYPE T_GRID

    TYPE (T_GRID) :: GRID

  CONTAINS

!**** Read ROMS GRID & HIS NetCDF file ************************************

      SUBROUTINE read_ROMS_files( GRID_FILE, HIS_FILE, resol)
      
      character(len=*), intent( in) :: GRID_FILE
      character(len=*), intent( in) :: HIS_FILE
      integer, intent( in) :: resol

      real, allocatable :: roms_h(:,:) 
      real, allocatable :: roms_mask(:,:)
      real, allocatable :: roms_sal(:,:)   ! salinity (psu)
      real, allocatable :: roms_sal_ave(:,:)   ! salinity (psu)
      real, allocatable :: roms_sal_max(:,:)   ! salinity (psu)
      real, allocatable :: roms_sal_min(:,:)   ! salinity (psu)
      
      integer :: N_xi_rho, N_eta_rho, N_time
      integer :: ncid, var_id
      integer :: i,j
      
!---- Read ROMS grid netCDF file --------------------------------

      write(*,*) "OPEN: ", GRID_FILE
      
      ! Open NetCDF grid file
      call check( nf90_open(GRID_FILE, nf90_nowrite, ncid) )
      ! Get dimension data
      call get_dimension(ncid, 'xi_rho',  N_xi_rho)
      call get_dimension(ncid, 'eta_rho', N_eta_rho)
      allocate( roms_h   (N_xi_rho, N_eta_rho)     )
      allocate( roms_mask(N_xi_rho, N_eta_rho)     )
      allocate( roms_sal (N_xi_rho, N_eta_rho)     )
      allocate( roms_sal_ave(N_xi_rho, N_eta_rho)     )
      allocate( roms_sal_max(N_xi_rho, N_eta_rho)     )
      allocate( roms_sal_min(N_xi_rho, N_eta_rho)     )
      allocate( GRID%N_x,   GRID%N_y   )
      allocate( GRID%Max_x, GRID%Max_y )
      allocate( GRID%Area,  GRID%N_tot)
      
      GRID%N_x   = N_xi_rho *resol*2
      GRID%N_y   = N_eta_rho*resol*2
      GRID%Max_x = N_xi_rho *resol
      GRID%Max_y = N_eta_rho*resol
      GRID%Area  = GRID%Max_x * GRID%Max_y
      GRID%N_tot = GRID%N_x * GRID%N_y
      
      allocate( GRID%h      (GRID%N_x, GRID%N_y) )
      allocate( GRID%mask   (GRID%N_x, GRID%N_y) )
      allocate( GRID%sal_ave(GRID%N_x, GRID%N_y) )
      allocate( GRID%sal_max(GRID%N_x, GRID%N_y) )
      allocate( GRID%sal_min(GRID%N_x, GRID%N_y) )
      
      GRID%N_x   = N_xi_rho *resol*2
      GRID%N_y   = N_eta_rho*resol*2
      GRID%Max_x = N_xi_rho *resol
      GRID%Max_y = N_eta_rho*resol
      GRID%Area  = GRID%Max_x * GRID%Max_y
      GRID%N_tot = GRID%N_x * GRID%N_y
      
      ! Get variable id
      call check( nf90_inq_varid(ncid, 'h', var_id     ) ) 
      call check( nf90_get_var(ncid, var_id, roms_h    ) )
      call check( nf90_inq_varid(ncid, 'mask_rho', var_id  ) ) 
      call check( nf90_get_var(ncid, var_id, roms_mask ) )
      ! Close NetCDF file
      call check( nf90_close(ncid) )
      
      write(*,*) "CLOSE: ", GRID_FILE
      
!---- Read ROMS his netCDF file --------------------------------

      write(*,*) "OPEN: ", HIS_FILE
      
      ! Open NetCDF grid file
      call check( nf90_open(HIS_FILE, nf90_nowrite, ncid) )
      
      ! Get variable id
!      call check( nf90_inq_varid(ncid, 'salt', var_id    ) ) 
!      call check( nf90_get_var(ncid, var_id, roms_sal    ) )
      ! Close NetCDF file
      call check( nf90_close(ncid) )
      
      write(*,*) "CLOSE: ", HIS_FILE
      
      
      CALL ex_grid1(N_xi_rho, N_eta_rho, resol, roms_h,    GRID%h   )
      CALL ex_grid2(N_xi_rho, N_eta_rho, resol, roms_mask, GRID%mask)
      
      Do i=1, N_eta_rho
        write(99,*) GRID%mask(:,i)
      End do
      
      END SUBROUTINE read_ROMS_files
      
      
!**** Expand grid resolution from ROMS to SEIV-DGVM ***************
      
      SUBROUTINE ex_grid1(Nx, Ny, resol, grd1, grd2)
      
      integer, intent( in) :: Nx,Ny
      integer, intent( in) :: resol
      real,    intent( in) :: grd1(Nx,Ny)
      real,    intent(out) :: grd2(Nx*resol,Ny*resol)
      
      integer :: i,j
      integer :: i2,j2
      
      Do j=1, Ny*resol
        Do i=1, Nx*resol
          i2 = int((i-1)/resol)+1
          j2 = int((j-1)/resol)+1
          grd2(i,j) = grd1(i2,j2)
        End do
      End do
      
      END SUBROUTINE ex_grid1
!-------------------------------------------------------------------
      SUBROUTINE ex_grid2(Nx, Ny, resol, grd1, grd2)
      
      integer, intent( in) :: Nx,Ny
      integer, intent( in) :: resol
      real,    intent( in) :: grd1(Nx,Ny)
      logical, intent(out) :: grd2(Nx*resol,Ny*resol)
      
      integer :: i,j
      integer :: i2,j2
      
      Do j=1, Ny*resol
        Do i=1, Nx*resol
          i2 = int((i-1)/resol)+1
          j2 = int((j-1)/resol)+1
          if(grd1(i2,j2)==1.0)then
            grd2(i,j) = .true.
          else
            grd2(i,j) = .false.
          end if
        End do
      End do
      
      END SUBROUTINE ex_grid2
      
!**** NetCDF utility **********************************************
      
      SUBROUTINE get_dimension(ncid, name, dim)
      
      integer,           intent( in) :: ncid
      character(len=*),  intent( in) :: name
      integer,           intent(out) :: dim

      integer :: dimid
      call check( nf90_inq_dimid(ncid, name, dimid) )
      call check( nf90_inquire_dimension(ncid, dimid, len=dim) )
      
      END SUBROUTINE get_dimension
! -------------------------------------------------------------------------
      SUBROUTINE  get_dimension2(ncid, name, dim, dims)
      integer,           intent( in) :: ncid
      character(len=*),  intent( in) :: name
      integer,           intent(out) :: dim
      real, allocatable, intent(out) :: dims(:)

      integer :: varid, dimid
      integer :: err

      call check( nf90_inq_dimid(ncid, name, dimid) )
      call check( nf90_inquire_dimension(ncid, dimid, len=dim) )
      allocate(dims(dim), stat=err)
      if (err /= 0) print *, name, ": Allocation request denied"
      call check( nf90_inq_varid(ncid, name, varid) )
      call check( nf90_get_var(ncid, varid, dims) )
      END SUBROUTINE get_dimension2

! -------------------------------------------------------------------------

      SUBROUTINE check(status)
      
      integer, intent(in) :: status

      if (status /= nf90_noerr) then 
          print *, trim(nf90_strerror(status))
          stop "Stopped"
      end if
      
      END SUBROUTINE check
      
! -------------------------------------------------------------------------
      
      SUBROUTINE check2(status, err_flag)
      
      integer, intent( in) :: status
      integer, intent(out) :: err_flag
      
      err_flag = 0

      if (status /= nf90_noerr) then 
          print *, trim(nf90_strerror(status))
          err_flag = 1
!          stop "Stopped"
      end if
      
      END SUBROUTINE check2

  END MODULE mod_grid


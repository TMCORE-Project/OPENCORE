module parameters_mod
  use constants_mod
  implicit none
  
  ! Namelist parameters
  ! Domain
  real    :: dx     !  grid-spacing in the x-direction
  real    :: dy     !  grid-spacing in the y-direction
  
  integer, parameter :: xhalo = 6 !  halo number of x-diretion
  integer, parameter :: yhalo = 6 !  halo number of y-diretion
  
  ! Index parameter
  integer :: ids      ! The starting index in the x-direction (Physical domain)
  integer :: ide      ! The ending index in the x-direction  (Physical domain)
  integer :: jds      ! The starting index in the y-direction  (Physical domain)
  integer :: jde      ! The ending index in the y-direction  (Physical domain)
  
  integer :: ips      ! The starting index in the x-direction (PV domain)
  integer :: ipe      ! The ending index in the x-direction  (PV domain)
  integer :: jps      ! The starting index in the y-direction  (PV domain)
  integer :: jpe      ! The ending index in the y-direction  (PV domain)
  
  integer :: ifs      ! The starting index of patch(face)
  integer :: ife      ! The ending index of patch(face)
                      
  integer :: Nx       ! Element numbers in the x-direction
  integer :: Ny       ! Element numbers in the y-direction
  
  integer, parameter :: Nf = 6           ! Number of cube faces
  
  !real, parameter :: x_min = -45.   !  start location of x-direction
  !real, parameter :: x_max =  45.   !  end location of x-direction
  !real, parameter :: y_min = -45.   !  start location of y-direction
  !real, parameter :: y_max =  45.   !  end location of y-direction
  
  real, parameter :: x_min = 0.   !  start location of x-direction
  real, parameter :: x_max = 90.  !  end location of x-direction
  real, parameter :: y_min = 0.   !  start location of y-direction
  real, parameter :: y_max = 90.  !  end location of y-direction
  
  integer :: nPVHalo
  
  namelist /domain/ dx,dy
  
  contains
  
  subroutine readNamelist
    
    open(1, file = 'namelist.input',status='old')
    read(1, nml  = domain       )
    close(1)
    
  end subroutine readNamelist
  
  subroutine initParameters
    
    ! Setting default values
    dx    = 2.
    dy    = 2.
    
    ! Read namelist
    call readNamelist
    
    ! Check if dx and dy are avaliable to achieve the integral element number
    if( (x_max - x_min)/dx - int((x_max - x_min)/dx) /= 0. )then
      stop '90 divide dx must be integer, choose another dx'
    end if
    
    if( (y_max - y_min)/dy - int((y_max - y_min)/dy) /= 0. )then
        stop '90 divide dy must be integer, choose another dy'
    end if
    
    ! Calculate element numbers on x/y direction
    Nx = int((x_max - x_min)/dx)
    Ny = int((y_max - y_min)/dy)
    Nx = 2 * Nx + 1
    Ny = 2 * Ny + 1
    
    ! Calculate starting and ending index for physical domain
    ids  = 1
    ide  = nx
    jds  = 1
    jde  = ny
    
    ! Calculate starting and ending index for memory array
    ips  = ids - xhalo
    ipe  = ide + xhalo
    jps  = jds - yhalo
    jpe  = jde + yhalo
    
    nPVHalo = xhalo
    
    ! Setting the starting patch index and ending patch index
    ifs = 1
    ife = Nf
    
    ! Convert Degree to map coordinate
    dx = dx * D2R
    dy = dy * D2R
    
  end subroutine initParameters
  
end module parameters_mod
    
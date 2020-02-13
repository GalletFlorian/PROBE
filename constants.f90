module constants

implicit none

integer, parameter :: dp = selected_real_kind(15, 307)
integer, parameter :: qp = selected_real_kind(33, 4931)



real(dp), parameter :: G         =  6.6742367d-11      !m^3.kg^-1.s^-2
real(dp), parameter :: AU        =  1.49598d11         !m
real(dp), parameter :: yr        =  365.25*24*3.6d3    !s
real(dp), parameter :: day       =  24*3.6d3           !s
real(dp), parameter :: hr        =  3.6d3              !s
real(dp), parameter ::  pi = 3.14159

real(dp), parameter :: Msun      =  1.98892d30               !kg
real(dp), parameter :: Mjup      =  1.8986112d27	          !kg
real(dp), parameter :: Mearth    =  3.d-6 * Msun             !kg
real(dp), parameter :: Rjup      =  69173.d3                 !m
real(dp), parameter :: Rearth    =  6371.0d3                 !m
real(dp), parameter :: Rsun      =  6.96d8                   !m

real(dp), parameter :: spinsun   =  2.87d-6                   !s-1

!moment of inertia (I/mR^2)
real(dp), parameter :: rg2_sun    =  5.9d-2
real(dp), parameter :: rg2_dM     =  2.0d-1
real(dp), parameter :: rg2_jup    =  2.54d-1
real(dp), parameter :: rg2_earth  =  3.308d-1


!_______Control braking law_______ 
  integer :: brklaw = 0   !brklaw = braking law for magnbrak.f, brklaw=3 = Matt 2015
  real :: Kbr=1.0
  real :: mmp=0.22
  real :: K1MP = 1.8d0
  real :: K2MP = 0.0506d0
  real :: K3 = 2.8
  real :: K4 = 0.2
  real :: mvic = 0.11
  real(dp) :: chi = 10.0
  real(dp) :: p = 2.0
  real(dp) ::  taudecinit = 6.0d6
  real(dp) :: tdiskparam = 2.4d6
  real(dp), parameter :: mstar = 1.1d0
  real(dp), parameter :: initper = 1.0d0
!_________________________________ 



end module constants

module monin_obukhov_functions_mod

#define _PURE
!#include <fms_platform.h>

use, intrinsic :: ieee_arithmetic
use integrate_mod, only : integrate_romberg_trapezoid, integrate_romberg_midpoint, integrate_romberg_midpoint_inv
use rsl_functions_mod, only : rsl_functions_T

implicit none
private

public :: most_functions_T
public :: make_most1_functions, make_most2_functions, make_brutsaert_functions, &
          make_neutral_functions

! representation of Monin-Obukhov Similarity Theory (MOST) stability correction functions
type, abstract :: most_functions_T
  logical :: neutral = .FALSE. ! only true for neutral stability functions
  real :: rich_crit ! it is here because it is used in Monin-Obukhov solver for all stability options,
                    ! and in some stability functions
  class(rsl_functions_T), pointer :: rsl => NULL () ! pointer to RSL functions
  ! lookup tables for RSL integrals Im and It
  logical :: use_RSL_lookup = .TRUE.
  real, allocatable :: a(:)    ! coordinates along axis a
  real, allocatable :: b(:)    ! coordinates along axis b
  real, allocatable :: Im(:,:) ! values of integral Im
  real, allocatable :: It(:,:) ! values of integral It
contains
  procedure(most_derivative_function), deferred :: derivative_m ! stability correction for momentum
  procedure(most_derivative_function), deferred :: derivative_t ! stability correction for heat and tracers
  procedure(most_integral_function),   deferred :: integral_m   ! integral stability correction for momentum
  procedure(most_integral_function),   deferred :: integral_t   ! integral stability correction for heat
  procedure(most_integral_function),   deferred :: integral_q   ! integral stability correction for tracers

  procedure(most_stable_mix),          deferred :: stable_mix

  ! procedures related to roughness sublayer (RSL)
  procedure :: set_rsl_functions   ! assign RSL functions and possibly do preliminary calculations (e.g. tabulate additive part of RSL integrals)
  procedure :: add_rsl_integral_m => add_rsl_integral_m ! add RSL integral stability term -- and optionally its derivative -- for momentum
  procedure :: add_rsl_integral_t => add_rsl_integral_t ! add RSL integral stability term -- and optionally its derivative -- for heat
  procedure :: add_rsl_integral_q => add_rsl_integral_t ! add RSL integral stability term -- and optionally its derivative -- for tracers
                            ! currently the tracer term is the same as for heat
end type most_functions_T

abstract interface
  _PURE subroutine most_derivative_function(this,n,mask,zeta,phi,ier)
     import :: most_functions_T
     class(most_functions_T), intent(in)   :: this
     integer, intent(in   )                :: n
     logical, intent(in   ), dimension(n)  :: mask
     real   , intent(in   ), dimension(n)  :: zeta
     real   , intent(inout), dimension(n)  :: phi
     integer, intent(  out)                :: ier
  end subroutine most_derivative_function
  _PURE subroutine most_integral_function(this, n, mask, zeta, zeta_0, ln_z_z0, F, ier)
     import :: most_functions_T
     class(most_functions_T), intent(in)   :: this
     integer, intent(in   )                :: n
     logical, intent(in   ), dimension(n)  :: mask
     real   , intent(in   ), dimension(n)  :: zeta, zeta_0, ln_z_z0
     real   , intent(inout), dimension(n)  :: F
     integer, intent(  out)                :: ier
  end subroutine most_integral_function
  _PURE subroutine most_stable_mix(this, n, rich, mix, ier)
     import :: most_functions_T
     class(most_functions_T), intent(in)   :: this
     integer, intent(in   )                :: n
     real   , intent(in   ), dimension(n)  :: rich
     real   , intent(  out), dimension(n)  :: mix
     integer, intent(  out)                :: ier
  end subroutine most_stable_mix
end interface

type, extends(most_functions_T) :: neutral_functions_T
! stable option 1
contains
  procedure :: derivative_m => neutral_deriv_m
  procedure :: derivative_t => neutral_deriv_t
  procedure :: integral_m   => neutral_integral_m
  procedure :: integral_t   => neutral_integral_tq
  procedure :: integral_q   => neutral_integral_tq
  procedure :: stable_mix   => neutral_stable_mix
end type neutral_functions_T

type, extends(most_functions_T) :: most1_functions_T
! stable option 1
contains
  procedure :: derivative_m => most1_deriv_m
  procedure :: derivative_t => most1_deriv_t
  procedure :: integral_m   => most1_integral_m
  procedure :: integral_t   => most1_integral_tq
  procedure :: integral_q   => most1_integral_tq
  procedure :: stable_mix   => most1_stable_mix
end type most1_functions_T

type, extends(most_functions_T) :: most2_functions_T
! stable option 1
  real :: zeta_trans
contains
  procedure :: derivative_m => most2_deriv_m
  procedure :: derivative_t => most2_deriv_t
  procedure :: integral_m   => most2_integral_m
  procedure :: integral_t   => most2_integral_tq
  procedure :: integral_q   => most2_integral_tq
  procedure :: stable_mix   => most2_stable_mix
end type most2_functions_T

type, extends(most_functions_T) :: brutsaert_functions_T
  real :: a_s = 6.1,  b_s = 2.5 ! parameters of stable regime
  real :: a_u = 0.33, b_u = 0.41, c_u = 0.33, d_u = 0.057, n_u = 0.78 ! parameters of unstable regime
contains
  procedure :: derivative_m => brutsaert_deriv_m
  procedure :: derivative_t => brutsaert_deriv_t
  procedure :: integral_m   => brutsaert_integral_m
  procedure :: integral_t   => brutsaert_integral_tq
  procedure :: integral_q   => brutsaert_integral_tq
  procedure :: stable_mix   => brutsaert_stable_mix
end type brutsaert_functions_T

real, parameter :: &
  RSL_UPPER_LIMIT = HUGE(1.0)*1e-4, & ! actual upper limit when integrating RSL corrections to "infinity"
  RSL_RTOL        = 1e-8              ! relative tolerance for RSL integrals; setting it to 1e-7 or below results in significant artifacts

contains ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!>/brief Set up roughness sublayer (RSL) parameterization
!!
!! given a pointer to RSL object, stores with the Monin-Obukhov stability correction
!! functions, and calculates the look-up table for the RSL integrals
subroutine set_rsl_functions(this,rsl, use_RSL_lookup, a_min, a_max, a_nsteps, b_min, b_max, b_nsteps)
  class(most_functions_T), intent(inout) :: this
  class(rsl_functions_T),  pointer       :: rsl !> pointer to RSL function object
  logical, intent(in) :: use_RSL_lookup !> use loookup tables to compute RSL integrals; otherwise
                                   !! calculate integrals directly: this can be used for, say,
                                   !! testing of quality of lookup in a single point runs, but would
                                   !! very likely be prohibitively slow in global simulations
  real,    intent(in) :: a_min     !> lower lookup table limit for parameter a of I_m and I_h RSL integrals: a_min > 0.
  real,    intent(in) :: a_max     !> upper lookup table limit for parameter a of I_m and I_h RSL integrals: a_max > a_min > 0.
  integer, intent(in) :: a_nsteps  !> number of lookup table steps along the axis a.
  real,    intent(in) :: b_min     !> lower lookup table limit for parameter b of I_m and I_h RSL integrals.
  real,    intent(in) :: b_max     !> upper lookup table limit for parameter b of I_m and I_h RSL integrals
  integer, intent(in) :: b_nsteps  !> number of lookup table steps along the axis b.

  integer :: i,j
  integer :: ierr
  real    :: x0,x1,x, y0,y1,y

  this%rsl => rsl
  if (.not.associated(this%rsl)) return ! don't do anything further

  this%use_RSL_lookup = use_RSL_lookup
  if (this%use_RSL_lookup) then
     if (allocated(this%a))  deallocate(this%a)
     if (allocated(this%b))  deallocate(this%b)
     if (allocated(this%Im)) deallocate(this%Im)
     if (allocated(this%It)) deallocate(this%It)

     allocate(this%a(a_nsteps+1),             &
              this%b(b_nsteps+1),             &
              this%Im(a_nsteps+1,b_nsteps+1), &
              this%It(a_nsteps+1,b_nsteps+1))

     x0 = sqrt(a_min); x1 = sqrt(a_max)
     do i = 1,a_nsteps+1
        x = x0+(x1-x0)/a_nsteps*(i-1)
        this%a(i) = x**2
     enddo

     ! NOTE: sign (a, b) returns the absolute value of a times the sign of b
     y0 = sign(abs(b_min)**(1./3.),b_min); y1 = sign(abs(b_max)**(1./3.),b_max)
     do j = 1,b_nsteps+1
        y = y0+(y1-y0)/b_nsteps*(j-1)
        this%b(j) = y**3
     enddo

     do i = 1,a_nsteps+1
     do j = 1,b_nsteps+1
        call RSL_integral_I_m(this,this%a(i),this%b(j),this%Im(i,j),ierr)
        call RSL_integral_I_t(this,this%a(i),this%b(j),this%It(i,j),ierr)
     enddo
     enddo
!      write(*,*) 'a ='
!      write(*,'(10g14.5)')this%a
!      write(*,*) 'b ='
!      write(*,'(10g14.5)')this%b
  endif
end subroutine set_rsl_functions

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! add the value of integral stability function roughness sublayer correction for momentum
_PURE subroutine add_rsl_integral_m(this, n, mask, l_inv, z1, z2, zR, F, df, ierr)
   class(most_functions_T), intent(in) :: this
   integer, intent(in)    :: n          ! size of the input/output arrays
   logical, intent(in)    :: mask(n)    ! don't do calculations where this mask is FALSE
   real,    intent(in)    :: l_inv(n)   ! 1/L, reciprocal of Monin-Obukhov length
   real,    intent(in)    :: z1(n)      ! lower limit of the RSL integral, m
   real,    intent(in)    :: z2(n)      ! upper limit of the RSL integral, m
   real,    intent(in)    :: zR(n)      ! roughness sublayer length scale, m
   ! the following inout arguments are updated (incremented) by this subroutine
   real,    intent(inout), optional :: F (n) ! value of the integral function
   real,    intent(inout), optional :: df(n) ! derivative of the integral function w.r.t. zeta
                                        ! where zeta is assumed to be z2/L
   integer, intent(out),   optional :: ierr  ! error code

   real, parameter :: delta_l_inv = 0.01 ! small increment of 1/L for derivative calculation

   integer :: i
   real :: R0 ! value of RSL integral for given parameters
   real :: R1 ! value of RSL integral with small zeta increment, for derivative calculations

   if (.not.associated(this%rsl)) return ! don't do anything if ther is no RSL

   do i = 1, n
      if (.not.mask(i))    cycle ! skip maske-out points
      if (.not.zR(i)>0) cycle ! skip points without roughness sublayer

      if (this%use_RSL_lookup) then
         call lookup_R_rsl(this, z1(i),z2(i),zR(i), l_inv(i), this%Im, R0, ierr)
      else
         call integralR_m_rsl(this, z1(i),z2(i),zR(i), l_inv(i), R0, ierr)
      endif
      if (present(F)) F(i) = F(i) - R0
      ! derivative of RSL correction w.r.t has to be calculated numerically
      if (present(df)) then
         if (this%use_RSL_lookup) then
            call lookup_R_rsl(this, z1(i),z2(i),zR(i), l_inv(i)+delta_l_inv, this%Im, R1, ierr)
         else
            call integralR_m_rsl(this, z1(i),z2(i),zR(i), l_inv(i)+delta_l_inv, R1, ierr)
         endif
         dF(i) = dF(i) - (R1-R0)/(delta_l_inv*z2(i))
      endif
   enddo
end subroutine add_rsl_integral_m

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! add the value of integral stability function roughness sublayer correction for heat
_PURE subroutine add_rsl_integral_t(this, n, mask, l_inv, z1, z2, zR, F, df, ierr)
   class(most_functions_T), intent(in) :: this
   integer, intent(in)    :: n          ! size of the input/output arrays
   logical, intent(in)    :: mask(n)    ! don't do calculations where this mask is FALSE
   real,    intent(in)    :: l_inv(n)   ! 1/L, reciprocal of Monin-Obukhov length
   real,    intent(in)    :: z1(n)      ! lower limit of the RSL integral, m
   real,    intent(in)    :: z2(n)      ! upper limit of the RSL integral, m
   real,    intent(in)    :: zR(n)      ! roughness sublayer length scale, m
   ! the following inout arguments are updated (incremented) by this subroutine
   real,    intent(inout), optional :: F (n) ! value of the integral function
   real,    intent(inout), optional :: df(n) ! derivative of the integral function w.r.t. zeta
                                        ! where zeta is assumed to be z2/L
   integer, intent(out),   optional :: ierr  ! error code

   real, parameter :: delta_l_inv = 0.01 ! small increment of 1/L for derivative calculation

   integer :: i
   real :: R0 ! value of RSL integral for given parameters
   real :: R1 ! value of RSL integral with small zeta increment, for derivative calculations

   if (.not.associated(this%rsl)) return ! don't do anything if ther is no RSL

   do i = 1, n
      if (.not.mask(i))    cycle ! skip maske-out points
      if (.not.zR(i)>0) cycle ! skip points without roughness sublayer

      if (this%use_RSL_lookup) then
         call lookup_R_rsl(this, z1(i),z2(i),zR(i), l_inv(i), this%It, R0, ierr)
      else
         call integralR_t_rsl(this, z1(i), z2(i), zR(i), l_inv(i), R0, ierr)
      endif
      if (present(F)) F(i) = F(i) - R0
      ! derivative of RSL correction w.r.t has to be calculated numerically
      if (present(df)) then
         if (this%use_RSL_lookup) then
            call lookup_R_rsl(this, z1(i), z2(i), zR(i), l_inv(i)+delta_l_inv, this%It, R1, ierr)
         else
            call integralR_t_rsl(this, z1(i), z2(i), zR(i), l_inv(i)+delta_l_inv, R1, ierr)
         endif
         dF(i) = dF(i) - (R1-R0)/(delta_l_inv*z2(i))
      endif
   enddo
end subroutine add_rsl_integral_t


_PURE subroutine lookup_R_rsl(most, z1, z2, z_rsl, l_inv, table, s, ierr)
  class(most_functions_T), intent(in) :: most
  real,    intent(in)  :: z1     !< lower limit of the integral R, m
  real,    intent(in)  :: z2     !< upper limit of the integral R, m
  real,    intent(in)  :: z_rsl  !< roughness sublayer length scale, m
  real,    intent(in)  :: l_inv  !< reciprocal of Monin-Obukhov length, 1/m
  real,    intent(in)  :: table(:,:) !< lookup table, Im or It
  real,    intent(out) :: s      !< value of the integral
  integer, intent(out) :: ierr   !< error code

  real :: a1, a2, b, s1, s2
  s  = ieee_value( s, ieee_signaling_nan )

  a1 = z1/z_rsl
  a2 = z2/z_rsl
  b  = z_rsl*l_inv

  call lookup_I_rsl(most,a1,b,table,s1,ierr); if (ierr.ne.0) return
  call lookup_I_rsl(most,a2,b,table,s2,ierr); if (ierr.ne.0) return
  s = s1 - s2
end subroutine lookup_R_rsl

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
_PURE subroutine lookup_I_rsl(most,a,b,table,s,ierr)
  class(most_functions_T), intent(in) :: most
  real,    intent(in)  :: a      !< parameter of the integral, z_1/z_R
  real,    intent(in)  :: b      !< parameter of the integral, z_R/L
  real,    intent(in)  :: table(:,:) !< lookup table, Im or It
  real,    intent(out) :: s      !< value of the integral
  integer, intent(out) :: ierr   !< error code, 0 = no error

  integer :: i,j
  real    :: da,db,f1,f2

  s    = ieee_value( s, ieee_signaling_nan )
  ierr = 1
  i = bisect(most%a,a)
  if (i<1.or.i>=size(most%a)) then
      ! bisect did not find appropriate interval for interpolation
      ! write(*,*) 'a out of bounds :',a
      return
  endif
  j = bisect(most%b,b,extrapolate_high=.TRUE.)
  if (j<1.or.j>=size(most%b)) then
      ! bisect did not find appropriate interval for interpolation
      ! write(*,*) 'b out of bounds :',b
      return
  endif

  da = (a-most%a(i))/(most%a(i+1)-most%a(i))
  f1 = table(i,j  )*(1-da)+table(i+1,j  )*da
  f2 = table(i,j+1)*(1-da)+table(i+1,j+1)*da

  db = (b-most%b(j))/(most%b(j+1)-most%b(j))
  s  = f1*(1-db) + f2*db
  ierr = 0
end subroutine lookup_I_rsl

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Given an array xx(n) in ascending order, and coordinate x1, finds
! the interval within xx that x1 belongs to, returning index i such that xx(i)<=x1<=xx(i+1)
!
! Normally this function returns -1 if the coordinate is out of range of the input array,
! but optional parameters extrapolate_high and extrapolate_low allow to control that
! behavior. If extrapolate_high is TRUE and x1 is higher than xx(n) [the upper boundary
! of input array values, since it is assumed to be in ascending order] (n-1) is
! returned. Similarly, if extrapolate_low is TRUE and the coordinate is less than the
! lower boundary of xx, 1 is returned.
_PURE integer function bisect(xx, x1, extrapolate_high, extrapolate_low)
  real, intent(in) :: xx(:) ! array of boundaries
  real, intent(in) :: x1    ! point to locate
  logical, intent(in), optional :: extrapolate_high ! if TRUE, extrapolate at the high end
  logical, intent(in), optional :: extrapolate_low  ! if TRUE, extrapolate at the low end

   ! ---- local vars
  real    :: x              ! duplicate of input value
  integer :: low, high, mid
  integer :: n              ! size of the input array
  logical :: extr_h, extr_l

  extr_h = .FALSE.; if (present(extrapolate_high)) extr_h = extrapolate_high
  extr_l = .FALSE.; if (present(extrapolate_low))  extr_l = extrapolate_low

  n = size(xx)
  x = x1

  ! find the coordinates
  if (x >= xx(1).and.x<=xx(n)) then
     low = 1; high = n
     do while (high-low > 1)
        mid = (low+high)/2
        if (xx(mid) <= x) then
           low = mid
        else
           high = mid
        endif
     enddo
     bisect = low
  else if (x>xx(n).and.extr_h) then
     bisect = n-1
  else if (x<xx(1).and.extr_l) then
     bisect = 1
  else
     bisect = -1
  endif
end function bisect

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
_PURE subroutine integralR_m_rsl(most, z1, z2, z_rsl, l_inv, s, ierr)
  class(most_functions_T), intent(in) :: most
  real,    intent(in)  :: z1, z2 ! lower and upper limits of the integral, m
  real,    intent(in)  :: z_rsl  ! roughness sublayer length scale, m
  real,    intent(in)  :: l_inv  ! reciprocal of Monin-Obukhov length, 1/m
  real,    intent(out) :: s      ! value of the integral
  integer, intent(out) :: ierr   ! error code

!   call integrate_romberg_trapezoid(f1 ,a,RSL_UPPER_LIMIT,RSL_RTOL,s,ierr)
  call integrate_romberg_midpoint(f1,z1,z2,RSL_RTOL,s,ierr)

contains
  ! internal function that returns the integrand
  _PURE real function f1(x)
     real, intent(in) :: x

     logical :: mask_1(1)
     real    :: phi_1(1), l_inv_1(1)
     real    :: rsl
     integer :: ierr_ignored

     ! calculate stability correction function
     mask_1  = .TRUE.; l_inv_1 = l_inv
     call most%derivative_m(1,mask_1,x*l_inv_1,phi_1,ierr_ignored)
     ! calculate roughness sublayer correction function
     rsl = most%rsl%rsl_m(x/z_rsl)
     ! finally, function under the integral
     f1 = phi_1(1)*(1-rsl)/x
  end function f1
end subroutine integralR_m_rsl

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
_PURE subroutine integralR_t_rsl(most,z1, z2, z_rsl, l_inv, s, ierr)
  class(most_functions_T), intent(in) :: most
  real,    intent(in)  :: z1, z2 ! lower and upper limits of the integral, m
  real,    intent(in)  :: z_rsl  ! roughness sublayer length scale, m
  real,    intent(in)  :: l_inv  ! reciprocal of Monin-Obukhov length, 1/m
  real,    intent(out) :: s      ! value of the integral
  integer, intent(out) :: ierr   ! error code

!   call integrate_romberg_trapezoid(f1 ,a,RSL_UPPER_LIMIT,RSL_RTOL,s,ierr)
  call integrate_romberg_midpoint(f1,z1,z2,RSL_RTOL,s,ierr)

contains
  ! internal function that returns the integrand
  _PURE real function f1(x)
     real, intent(in) :: x

     logical :: mask_1(1)
     real    :: phi_1(1), l_inv_1(1)
     real    :: rsl
     integer :: ierr_ignored

     ! calculate stability correction function
     mask_1  = .TRUE.; l_inv_1 = l_inv
     call most%derivative_t(1,mask_1,x*l_inv_1,phi_1,ierr_ignored)
     ! calculate roughness sublayer correction function
     rsl = most%rsl%rsl_t(x/z_rsl)
     ! finally, function under the integral
     f1 = phi_1(1)*(1-rsl)/x
  end function f1
end subroutine integralR_t_rsl

! ==== neutral stability option =========================================================
function make_neutral_functions(rich_crit) result(ptr)
   class(neutral_functions_T), pointer :: ptr
   real, intent(in) :: rich_crit

   allocate(ptr)
   ptr%rich_crit = rich_crit
   ptr%neutral   = .TRUE.
end function make_neutral_functions

! neutral stability functions are not really used: instead, a simplified non-iterative
! special case solution is employed by Monin-Obukhov kernel module
_PURE subroutine neutral_deriv_m(this,n,mask,zeta,phi,ier)
  class(neutral_functions_T), intent(in) :: this
  integer, intent(in   )                :: n
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(in   ), dimension(n)  :: zeta
  real   , intent(inout), dimension(n)  :: phi
  integer, intent(  out)                :: ier

  ier = 0
  phi = 1.0
end subroutine neutral_deriv_m

_PURE subroutine neutral_deriv_t(this,n,mask,zeta,phi,ier)
  class(neutral_functions_T), intent(in) :: this
  integer, intent(in   )                :: n
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(in   ), dimension(n)  :: zeta
  real   , intent(inout), dimension(n)  :: phi
  integer, intent(  out)                :: ier

  ier = 0
  phi = 1.0
end subroutine neutral_deriv_t

_PURE subroutine neutral_integral_m(this, n, mask, zeta, zeta_0, ln_z_z0, F, ier)
  class(neutral_functions_T), intent(in) :: this
  integer, intent(in   )                :: n
  real   , intent(in   ), dimension(n)  :: zeta, zeta_0, ln_z_z0
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(inout), dimension(n)  :: F
  integer, intent(  out)                :: ier

  ier = 0
  F = ln_z_z0
end subroutine neutral_integral_m

_PURE subroutine neutral_integral_tq(this, n, mask, zeta, zeta_0, ln_z_z0, F, ier)
  class(neutral_functions_T), intent(in) :: this
  integer, intent(in   )                :: n
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(in   ), dimension(n)  :: zeta, zeta_0, ln_z_z0
  real   , intent(inout), dimension(n)  :: F
  integer, intent(  out)                :: ier

  ier = 0
  F = ln_z_z0
end subroutine neutral_integral_tq

_PURE subroutine neutral_stable_mix(this, n, rich, mix, ier)
  class(neutral_functions_T), intent(in) :: this
  integer, intent(in   )                :: n
  real   , intent(in   ), dimension(n)  :: rich
  real   , intent(  out), dimension(n)  :: mix
  integer, intent(  out)                :: ier

  ier = 0

  mix = 0.0
  where (rich > 0.0) mix = 1.0
end subroutine neutral_stable_mix

! ==== first stability option ===========================================================
function make_most1_functions(rich_crit) result(ptr)
   class(most1_functions_T), pointer :: ptr
   real, intent(in) :: rich_crit

   allocate(ptr)
   ptr%rich_crit = rich_crit
end function make_most1_functions

_PURE subroutine most1_deriv_m(this,n,mask,zeta,phi,ier)
  class(most1_functions_T), intent(in)    :: this
  integer, intent(in   )                :: n
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(in   ), dimension(n)  :: zeta
  real   , intent(inout), dimension(n)  :: phi
  integer, intent(  out)                :: ier

  logical, dimension(n) :: stable, unstable
  real   , dimension(n) :: x
  real                  :: b_stab

  ier = 0
  b_stab   = 1.0/this%rich_crit

  stable   = mask .and. zeta >= 0.0
  unstable = mask .and. zeta <  0.0

  where (unstable)
     x     = (1 - 16.0*zeta  )**(-0.5)
     phi = sqrt(x)  ! phi = (1 - 16.0*zeta)**(-0.25)
  end where
  where (stable)
     phi = 1.0 + zeta  *(5.0 + b_stab*zeta)/(1.0 + zeta)
  end where
end subroutine most1_deriv_m

_PURE subroutine most1_deriv_t(this,n,mask,zeta,phi,ier)
  class(most1_functions_T), intent(in)    :: this
  integer, intent(in   )                :: n
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(in   ), dimension(n)  :: zeta
  real   , intent(inout), dimension(n)  :: phi
  integer, intent(  out)                :: ier

  logical, dimension(n) :: stable, unstable
  real                  :: b_stab

  ier = 0
  b_stab     = 1.0/this%rich_crit

  stable   = mask .and. zeta >= 0.0
  unstable = mask .and. zeta <  0.0

  where (unstable)
     phi = (1 - 16.0*zeta)**(-0.5)
  end where
  where (stable)
     phi = 1.0 + zeta*(5.0 + b_stab*zeta)/(1.0 + zeta)
  end where
end subroutine most1_deriv_t

_PURE subroutine most1_integral_m(this, n, mask, zeta, zeta_0, ln_z_z0, F, ier)
  class(most1_functions_T), intent(in)     :: this
  integer, intent(in   )                :: n
  real   , intent(in   ), dimension(n)  :: zeta, zeta_0, ln_z_z0
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(inout), dimension(n)  :: F
  integer, intent(  out)                :: ier

  real                   :: b_stab
  real, dimension(n) :: x, x_0, x1, x1_0, num, denom, y
  logical, dimension(n) :: stable, unstable

  ier = 0

  b_stab     = 1.0/this%rich_crit

  stable   = mask .and. zeta >= 0.0
  unstable = mask .and. zeta <  0.0

  where(unstable)
     x     = sqrt(1 - 16.0*zeta)
     x_0   = sqrt(1 - 16.0*zeta_0)

     x      = sqrt(x)
     x_0    = sqrt(x_0)

     x1     = 1.0 + x
     x1_0   = 1.0 + x_0

     num    = x1*x1*(1.0 + x*x)
     denom  = x1_0*x1_0*(1.0 + x_0*x_0)
     y      = atan(x) - atan(x_0)
     F  = ln_z_z0 - log(num/denom) + 2*y
  end where

  where (stable)
     F = ln_z_z0 + (5.0 - b_stab)*log((1.0 + zeta)/(1.0 + zeta_0)) &
          + b_stab*(zeta - zeta_0)
  end where
end subroutine most1_integral_m

_PURE subroutine most1_integral_tq(this, n, mask, zeta, zeta_0, ln_z_z0, F, ier)
  class(most1_functions_T), intent(in)    :: this
  integer, intent(in   )                :: n
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(in   ), dimension(n)  :: zeta, zeta_0, ln_z_z0
  real   , intent(inout), dimension(n)  :: F
  integer, intent(  out)                :: ier

  real, dimension(n)     :: x, x_t
  logical, dimension(n)  :: stable, unstable
  real                   :: b_stab

  ier = 0

  b_stab     = 1.0/this%rich_crit

  stable   = mask .and. zeta >= 0.0
  unstable = mask .and. zeta <  0.0

  where(unstable)
    x     = sqrt(1 - 16.0*zeta)
    x_t   = sqrt(1 - 16.0*zeta_0)

    F = ln_z_z0 - 2.0*log( (1.0 + x)/(1.0 + x_t) )
  end where
  where (stable)
    F = ln_z_z0 + (5.0 - b_stab)*log((1.0 + zeta)/(1.0 + zeta_0)) &
       + b_stab*(zeta - zeta_0)
  end where
end subroutine most1_integral_tq

_PURE subroutine most1_stable_mix(this, n, rich, mix, ier)
  class(most1_functions_T), intent(in)  :: this
  integer, intent(in   )                :: n
  real   , intent(in   ), dimension(n)  :: rich
  real   , intent(  out), dimension(n)  :: mix
  integer, intent(  out)                :: ier

  real    :: r, a, b, c, zeta, phi
  real    :: b_stab
  integer :: i

  ier = 0

  mix = 0.0
  b_stab     = 1.0/this%rich_crit

  c = - 1.0
  do i = 1, n
     if(rich(i) > 0.0 .and. rich(i) < this%rich_crit) then
        r = 1.0/rich(i)
        a = r - b_stab
        b = r - (1.0 + 5.0)
        zeta = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a)
        phi = 1.0 + b_stab*zeta + (5.0 - b_stab)*zeta/(1.0 + zeta)
        mix(i) = 1./(phi*phi)
     endif
  end do
end subroutine most1_stable_mix

! ==== second stability option ===========================================================
function make_most2_functions(rich_crit, zeta_trans) result(ptr)
   class(most2_functions_T), pointer :: ptr
   real, intent(in) :: rich_crit, zeta_trans

   allocate(ptr)
   ptr%rich_crit = rich_crit
   ptr%zeta_trans = zeta_trans
end function make_most2_functions

_PURE subroutine most2_deriv_m(this,n,mask,zeta,phi,ier)
  class(most2_functions_T), intent(in)    :: this
  integer, intent(in   )                :: n
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(in   ), dimension(n)  :: zeta
  real   , intent(inout), dimension(n)  :: phi
  integer, intent(  out)                :: ier

  logical, dimension(n) :: stable, unstable
  real   , dimension(n) :: x
  real                  :: b_stab, lambda

  ier = 0
  b_stab     = 1.0/this%rich_crit

  stable   = mask .and. zeta >= 0.0
  unstable = mask .and. zeta <  0.0

  where (unstable)
     x     = (1 - 16.0*zeta  )**(-0.5)
     phi = sqrt(x)  ! phi = (1 - 16.0*zeta)**(-0.25)
  end where

  lambda = 1.0 + (5.0 - b_stab)*this%zeta_trans

  where (stable .and. zeta < this%zeta_trans)
     phi = 1 + 5.0*zeta
  end where
  where (stable .and. zeta >= this%zeta_trans)
     phi = lambda + b_stab*zeta
  end where
end subroutine most2_deriv_m

_PURE subroutine most2_deriv_t(this,n,mask,zeta,phi,ier)
  class(most2_functions_T), intent(in)    :: this
  integer, intent(in   )                :: n
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(in   ), dimension(n)  :: zeta
  real   , intent(inout), dimension(n)  :: phi
  integer, intent(  out)                :: ier
  logical, dimension(n) :: stable, unstable
  real                  :: b_stab, lambda

  ier = 0
  b_stab     = 1.0/this%rich_crit

  stable   = mask .and. zeta >= 0.0
  unstable = mask .and. zeta <  0.0

  where (unstable)
     phi = (1 - 16.0*zeta)**(-0.5)
  end where

  lambda = 1.0 + (5.0 - b_stab)*this%zeta_trans

  where (stable .and. zeta < this%zeta_trans)
     phi = 1 + 5.0*zeta
  end where
  where (stable .and. zeta >= this%zeta_trans)
     phi = lambda + b_stab*zeta
  end where
end subroutine most2_deriv_t

_PURE subroutine most2_integral_m(this, n, mask, zeta, zeta_0, ln_z_z0, F, ier)
  class(most2_functions_T), intent(in)     :: this
  integer, intent(in   )                :: n
  real   , intent(in   ), dimension(n)  :: zeta, zeta_0, ln_z_z0
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(inout), dimension(n)  :: F
  integer, intent(  out)                :: ier

  real                   :: b_stab, lambda
  real, dimension(n) :: x, x_0, x1, x1_0, num, denom, y
  logical, dimension(n) :: stable, unstable, &
       weakly_stable, strongly_stable

  ier = 0

  b_stab     = 1.0/this%rich_crit

  stable   = mask .and. zeta >= 0.0
  unstable = mask .and. zeta <  0.0

  where(unstable)
     x     = sqrt(1 - 16.0*zeta)
     x_0   = sqrt(1 - 16.0*zeta_0)

     x      = sqrt(x)
     x_0    = sqrt(x_0)

     x1     = 1.0 + x
     x1_0   = 1.0 + x_0

     num    = x1*x1*(1.0 + x*x)
     denom  = x1_0*x1_0*(1.0 + x_0*x_0)
     y      = atan(x) - atan(x_0)
     F  = ln_z_z0 - log(num/denom) + 2*y
  end where

  lambda = 1.0 + (5.0 - b_stab)*this%zeta_trans

  weakly_stable   = stable .and. zeta <= this%zeta_trans
  strongly_stable = stable .and. zeta >  this%zeta_trans

  where (weakly_stable)
     F = ln_z_z0 + 5.0*(zeta - zeta_0)
  end where

  where(strongly_stable)
     x = (lambda - 1.0)*log(zeta/this%zeta_trans) + b_stab*(zeta - this%zeta_trans)
  end where

  where (strongly_stable .and. zeta_0 <= this%zeta_trans)
     F = ln_z_z0 + x + 5.0*(this%zeta_trans - zeta_0)
  end where
  where (strongly_stable .and. zeta_0 > this%zeta_trans)
     F = lambda*ln_z_z0 + b_stab*(zeta  - zeta_0)
  end where
end subroutine most2_integral_m

_PURE subroutine most2_integral_tq(this, n, mask, zeta, zeta_0, ln_z_z0, F, ier)
  class(most2_functions_T), intent(in)    :: this
  integer, intent(in   )                :: n
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(in   ), dimension(n)  :: zeta, zeta_0, ln_z_z0
  real   , intent(inout), dimension(n)  :: F
  integer, intent(  out)                :: ier

  real, dimension(n)     :: x, x_t
  logical, dimension(n)  :: stable, unstable, &
                             weakly_stable, strongly_stable
  real                   :: b_stab, lambda

  ier = 0

  b_stab     = 1.0/this%rich_crit

  stable   = mask .and. zeta >= 0.0
  unstable = mask .and. zeta <  0.0

  where(unstable)
    x     = sqrt(1 - 16.0*zeta)
    x_t   = sqrt(1 - 16.0*zeta_0)

    F = ln_z_z0 - 2.0*log( (1.0 + x)/(1.0 + x_t) )
  end where

  lambda = 1.0 + (5.0 - b_stab)*this%zeta_trans

  weakly_stable   = stable .and. zeta <= this%zeta_trans
  strongly_stable = stable .and. zeta >  this%zeta_trans

  where (weakly_stable)
    F = ln_z_z0 + 5.0*(zeta - zeta_0)
  end where

  where(strongly_stable)
    x = (lambda - 1.0)*log(zeta/this%zeta_trans) + b_stab*(zeta - this%zeta_trans)
  end where

  where (strongly_stable .and. zeta_0 <= this%zeta_trans)
    F = ln_z_z0 + x + 5.0*(this%zeta_trans - zeta_0)
  end where
  where (strongly_stable .and. zeta_0 > this%zeta_trans)
    F = lambda*ln_z_z0 + b_stab*(zeta  - zeta_0)
  end where
end subroutine most2_integral_tq

_PURE subroutine most2_stable_mix(this, n, rich, mix, ier)
  class(most2_functions_T), intent(in)  :: this
  integer, intent(in   )                :: n
  real   , intent(in   ), dimension(n)  :: rich
  real   , intent(  out), dimension(n)  :: mix
  integer, intent(  out)                :: ier

  real    :: b_stab, rich_trans, lambda

  ier = 0

  mix = 0.0
  b_stab     = 1.0/this%rich_crit
  rich_trans = this%zeta_trans/(1.0 + 5.0*this%zeta_trans)

  lambda = 1.0 + (5.0 - b_stab)*this%zeta_trans

  where(rich > 0.0 .and. rich <= rich_trans)
    mix = (1.0 - 5.0*rich)**2
  end where
  where(rich > rich_trans .and. rich < this%rich_crit)
    mix = ((1.0 - b_stab*rich)/lambda)**2
  end where
end subroutine most2_stable_mix

! ==== Brutsaert stability option ===========================================================
function make_brutsaert_functions(rich_crit) result(ptr)
   class(brutsaert_functions_T), pointer :: ptr
   real, intent(in) :: rich_crit

   allocate(ptr)
   ptr%rich_crit = rich_crit
end function make_brutsaert_functions

_PURE subroutine brutsaert_deriv_m(this,n,mask,zeta,phi,ier)
  class(brutsaert_functions_T), intent(in)    :: this
  integer, intent(in   )                :: n
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(in   ), dimension(n)  :: zeta
  real   , intent(inout), dimension(n)  :: phi
  integer, intent(  out)                :: ier

  integer :: i
  real    :: y

  ier = 0

  do i = 1, n
     if (.not.mask(i)) cycle ! skip the points that are masked out
     if (zeta(i) < 0.0) then ! unstable
        y = -zeta(i)
        if (y <= 1.0/this%b_u**3) then
           phi(i) = (this%a_u + this%b_u * y**(4.0/3.0))/(this%a_u + y)
        else
           phi(i) = 1.0
        endif
     else ! zeta(i) >= 0.0, stable
        phi(i) = 1 + this%a_s * (zeta(i) + zeta(i)**this%b_s*(1+zeta(i)**this%b_s)**(1.0/this%b_s-1))/ &
                                (zeta(i) + (1+zeta(i)**this%b_s)**(1.0/this%b_s))
     endif
  enddo
end subroutine brutsaert_deriv_m

_PURE subroutine brutsaert_deriv_t(this,n,mask,zeta,phi,ier)
  class(brutsaert_functions_T), intent(in)    :: this
  integer, intent(in   )                :: n
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(in   ), dimension(n)  :: zeta
  real   , intent(inout), dimension(n)  :: phi
  integer, intent(  out)                :: ier

  logical, dimension(n) :: stable, unstable
  real   , dimension(n) :: y

  ier = 0

  stable   = mask .and. zeta >= 0.0
  unstable = mask .and. zeta <  0.0

  where (unstable)
     y = -zeta
     phi = (this%c_u + this%d_u * y**this%n_u)/(this%c_u + y**this%n_u)
  end where

  where (stable)
     phi = 1 + this%a_s * (zeta + zeta**this%b_s*(1+zeta**this%b_s)**(1.0/this%b_s-1))/ &
                          (zeta + (1+zeta**this%b_s)**(1.0/this%b_s))
  end where
end subroutine brutsaert_deriv_t

elemental real function brutsaert_psi_m(this, zeta) result(psi_m)
  class(brutsaert_functions_T), intent(in) :: this
  real,                       intent(in) :: zeta

  real, parameter :: r3 = 1.0/3.0
  real, parameter :: s3 = sqrt(3.0)
  real, parameter :: pi = 3.1415926535
  real :: x, y, psi_0

  if (zeta >= 0.0) then
     ! stable
     psi_m = -this%a_s*log(zeta+(1+zeta**this%b_s)**(1.0/this%b_s))
  else
     ! unstable
     y = min(-zeta,1.0/this%b_u**3)
     x = (y/this%a_u)**r3
     psi_0 = -log(this%a_u) + s3*this%b_u*this%a_u**r3*pi/6 ! can be precomputed
     psi_m = log(this%a_u+y) - 3*this%b_u*y**r3 &
           + this%b_u * this%a_u**r3/2 * log((1+x)**2/(1-x+x**2)) &
           + s3*this%b_u*this%a_u**r3*atan((2*x-1)/s3) &
           + psi_0
  endif
end function brutsaert_psi_m

_PURE subroutine brutsaert_integral_m(this, n, mask, zeta, zeta_0, ln_z_z0, F, ier)
  class(brutsaert_functions_T), intent(in)    :: this
  integer, intent(in   )                :: n
  real   , intent(in   ), dimension(n)  :: zeta, zeta_0, ln_z_z0
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(inout), dimension(n)  :: F
  integer, intent(  out)                :: ier

  ier = 0

  where (mask)
     F = ln_z_z0 - brutsaert_psi_m(this, zeta) + brutsaert_psi_m(this, zeta_0)
  end where
end subroutine brutsaert_integral_m

elemental real function brutsaert_psi_h(this, zeta) result(psi_h)
  class(brutsaert_functions_T), intent(in) :: this
  real,                       intent(in) :: zeta

  real :: y

  if(zeta>=0) then
     !stable
     psi_h = -this%a_s*log(zeta+(1+zeta**this%b_s)**(1.0/this%b_s))
  else
     ! unstable
     y = -zeta
     psi_h = (1-this%d_u)/this%n_u * log((this%c_u+y**this%n_u)/this%c_u)
  endif
end function brutsaert_psi_h

_PURE subroutine brutsaert_integral_tq(this, n, mask, zeta, zeta_0, ln_z_z0, F, ier)
  class(brutsaert_functions_T), intent(in)    :: this
  integer, intent(in   )                :: n
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(in   ), dimension(n)  :: zeta, zeta_0, ln_z_z0
  real   , intent(inout), dimension(n)  :: F
  integer, intent(  out)                :: ier

  where (mask)
     F = ln_z_z0 - brutsaert_psi_h(this, zeta) + brutsaert_psi_h(this, zeta_0)
  end where
end subroutine brutsaert_integral_tq

_PURE subroutine brutsaert_stable_mix(this, n, rich, mix, ier)
  class(brutsaert_functions_T), intent(in) :: this
  integer, intent(in   )                :: n
  real   , intent(in   ), dimension(n)  :: rich
  real   , intent(  out), dimension(n)  :: mix
  integer, intent(  out)                :: ier

  ! NOT IMPLEMENTED: solving zeta for given Ri seems complicated because of the
  ! more complex formulation of stable phi_m and phi_h; possibly solve numerically,
  ! tabulate, and use lookup table?

  ier = 1
end subroutine brutsaert_stable_mix

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!> \brief Calculate RSL integral I for momentum
!!
!! given parameters a = z_atm/z_RRL and b = z_R/L, calculates the RSL integral I_m
!! for momentum
!! I_m = \int_a^\infty \phi_m (b z) (1-phi_{RSL,m}(z))/z dz
subroutine RSL_integral_I_m(most,a,b,s,ierr)
  class(most_functions_T), intent(in) :: most
  real, intent(in)     :: a    !< parameter of the integral, z_1/z_R
  real, intent(in)     :: b    !< parameter of the integral, z_R/L
  real, intent(out)    :: s    !< value of the integral
  integer, intent(out) :: ierr !< error code

  if (.not.associated(most%rsl)) then
     s = 0.0; ierr = 0
     return
  endif
  call integrate_romberg_midpoint_inv(f,a,RSL_UPPER_LIMIT,RSL_RTOL,s,ierr)
contains
  _PURE real function f(x)
     real, intent(in) :: x

     real,    dimension(1) :: zeta,phi,rsl
     logical, dimension(1) :: mask
     integer :: ierr_ignored

     mask=.TRUE.
     zeta = x*b
     call most%derivative_m(1,mask,zeta,phi,ierr_ignored)
     rsl = most%rsl%rsl_m(x)
     f = phi(1)*(1-rsl(1))/x
  end function f
end subroutine RSL_integral_I_m

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!> \brief Calculate RSL integral I for heat
!!
!! given parameters a = z_atm/z_RRL and b = z_R/L, calculates the RSL integral I_t
!! for heat
!! I_t = \int_a^\infty \phi_t (b z) (1-phi_{RSL,t}(z))/z dz
subroutine RSL_integral_I_t(most,a,b,s,ierr)
  class(most_functions_T), intent(in) :: most
  real, intent(in)     :: a    !< parameter of the integral, z_1/z_R
  real, intent(in)     :: b    !< parameter of the integral, z_R/L
  real, intent(out)    :: s    !< value of the integral
  integer, intent(out) :: ierr !< error code

  if (.not.associated(most%rsl)) then
     s = 0.0; ierr = 0
     return
  endif
  call integrate_romberg_midpoint_inv(f,a,RSL_UPPER_LIMIT,RSL_RTOL,s,ierr)
contains
  _PURE real function f(x)
     real, intent(in) :: x

     real,    dimension(1) :: zeta,phi,rsl
     logical, dimension(1) :: mask
     integer :: ierr_ignored

     mask=.TRUE.
     zeta = x*b
     call most%derivative_t(1,mask,zeta,phi,ierr_ignored)
     rsl = most%rsl%rsl_t(x)
     f = phi(1)*(1-rsl(1))/x
  end function f
end subroutine RSL_integral_I_t

end module monin_obukhov_functions_mod

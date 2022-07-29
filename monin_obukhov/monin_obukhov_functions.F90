module monin_obukhov_functions_mod
#include <fms_platform.h>

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
contains
  procedure(most_derivative_function), deferred :: derivative_m ! stability correction for momentum
  procedure(most_derivative_function), deferred :: derivative_t ! stability correction for heat and tracers
  procedure(most_integral_function),   deferred :: integral_m   ! integral stability correction for momentum
  procedure(most_integral_function),   deferred :: integral_t   ! integral stability correction for heat
  procedure(most_integral_function),   deferred :: integral_q   ! integral stability correction for tracers

  procedure(most_stable_mix),          deferred :: stable_mix

  ! procedures related to roughness sublayer (RSL)
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
! set up rougness sublayer parameterization
subroutine set_rsl_functions(this,rsl)
  class(most_functions_T), intent(inout) :: this
  class(rsl_functions_T),  pointer       :: rsl
  this%rsl => rsl
end subroutine set_rsl_functions

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! add the value of integral stability function roughness sublayer correction for momentum
pure subroutine add_rsl_integral_m(this, n, mask, zeta, z0, zA, zR, F, df, ierr)
   class(most_functions_T), intent(in) :: this
   integer, intent(in)    :: n          ! size of the input/output arrays
   logical, intent(in)    :: mask(n)    ! don't do calculations where this mask is FALSE
   real,    intent(in)    :: zeta(n)    ! zA/L, ratio of MO layer height to MO length scale
   real,    intent(in)    :: z0(n)      ! roughness length, m
   real,    intent(in)    :: zA(n)      ! MO layer thickness (height of the lowest atmos layer), m
   real,    intent(in)    :: zR(n)      ! roughness sublayer length scale, m
   real,    intent(inout), optional :: F (n) ! value of the integral function
   real,    intent(inout), optional :: df(n) ! derivative of the integral function w.r.t. zeta
   integer, intent(out),   optional :: ierr  ! error code

   real, parameter :: d_zeta = 0.01 ! small increment of zeta for derivative calculation

   integer :: i
   real :: R0 ! value of RSL integral for given parameters
   real :: R1 ! value of RSL integral with small zeta increment, for derivative calculations

   if (.not.associated(this%rsl)) return ! don't do anything if ther is no RSL

   do i = 1, n
      if (.not.mask(i))    cycle ! skip maske-out points
      if (.not.zR(i)>0) cycle ! skip points without roughness sublayer

      call integralR_m_rsl(this, z0(i),zA(i),zR(i), zeta(i), R0, ierr)
      if (present(F)) F(i) = F(i) - R0
      ! derivative of RSL correction w.r.t has to be calculated numerically
      if (present(df)) then
         call integralR_m_rsl(this, z0(i),zA(i),zR(i), zeta(i)+d_zeta, R1, ierr)
         dF(i) = dF(i) - (R1-R0)/d_zeta
      endif
   enddo
end subroutine add_rsl_integral_m

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! add the value of integral stability function roughness sublayer correction for momentum
pure subroutine add_rsl_integral_t(this, n, mask, zeta, z0, zA, zR, F, df, ierr)
   class(most_functions_T), intent(in) :: this
   integer, intent(in)    :: n          ! size of the input/output arrays
   logical, intent(in)    :: mask(n)    ! don't do calculations where this mask is FALSE
   real,    intent(in)    :: zeta(n)    ! zA/L, ratio of MO layer height to MO length scale
   real,    intent(in)    :: z0(n)      ! roughness length, m
   real,    intent(in)    :: zA(n)      ! MO layer thickness (height of the lowest atmos layer), m
   real,    intent(in)    :: zR(n)      ! roughness sublayer length scale, m
   real,    intent(inout), optional :: F (n) ! value of the integral function
   real,    intent(inout), optional :: df(n) ! derivative of the integral function w.r.t. zeta
   integer, intent(out),   optional :: ierr  ! error code

   real, parameter :: d_zeta = 0.01 ! small increment of zeta for derivative calculation

   integer :: i
   real :: R0 ! value of RSL integral for given parameters
   real :: R1 ! value of RSL integral with small zeta increment, for derivative calculations

   if (.not.associated(this%rsl)) return ! don't do anything if ther is no RSL

   do i = 1, n
      if (.not.mask(i))    cycle ! skip maske-out points
      if (.not.zR(i)>0) cycle ! skip points without roughness sublayer

      call integralR_t_rsl(this, z0(i),zA(i),zR(i), zeta(i), R0, ierr)
      if (present(F)) F(i) = F(i) - R0
      ! derivative of RSL correction w.r.t has to be calculated numerically
      if (present(df)) then
         call integralR_t_rsl(this, z0(i),zA(i),zR(i), zeta(i)+d_zeta, R1, ierr)
         dF(i) = dF(i) - (R1-R0)/d_zeta
      endif
   enddo
end subroutine add_rsl_integral_t

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pure subroutine integralR_m_rsl(most,z0,za,z_rsl,zeta, s, ierr)
  class(most_functions_T), intent(in) :: most
  real, intent(in)     :: z0, za, z_rsl, zeta
  real, intent(out)    :: s
  integer, intent(out) :: ierr

!   call integrate_romberg_trapezoid(f1 ,a,RSL_UPPER_LIMIT,RSL_RTOL,s,ierr)
  call integrate_romberg_midpoint(f1,z0,za,RSL_RTOL,s,ierr)

contains
  ! internal function that returns the integrand
  pure real function f1(x)
     real, intent(in) :: x

     logical :: mask(1)
     real    :: zeta1(1),phi(1)
     real    :: rsl
     integer :: ierr_ignored

     mask  = .TRUE.
     zeta1 = zeta
     call most%derivative_m(1,mask,x*zeta1/za,phi,ierr_ignored)
     rsl = most%rsl%rsl_m(x/z_rsl)
     f1 = phi(1)*(1-rsl)/x
  end function f1
end subroutine integralR_m_rsl

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pure subroutine integralR_t_rsl(most,z0,za,z_rsl,zeta, s, ierr)
  class(most_functions_T), intent(in) :: most
  real, intent(in) :: z0, za, z_rsl, zeta
  real, intent(out) :: s
  integer, intent(out) :: ierr

!   call integrate_romberg_trapezoid(f1 ,a,RSL_UPPER_LIMIT,RSL_RTOL,s,ierr)
  call integrate_romberg_midpoint(f1,z0,za,RSL_RTOL,s,ierr)

contains
  ! internal function that returns the integrand
  pure real function f1(x)
     real, intent(in) :: x

     logical :: mask(1)
     real    :: zeta1(1),phi(1)
     real    :: rsl
     integer :: ierr_ignored

     mask  = .TRUE.
     zeta1 = zeta
     call most%derivative_t(1,mask,x*zeta1/za,phi,ierr_ignored)
     rsl = most%rsl%rsl_t(x/z_rsl)
     f1 = phi(1)*(1-rsl)/x
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

  real, dimension(n)     :: x, x_t, x_q
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

end module monin_obukhov_functions_mod

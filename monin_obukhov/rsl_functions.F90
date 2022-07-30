module rsl_functions_mod

implicit none
private

public :: rsl_functions_T
public :: make_ridder2010_rsl_functions, make_ghannam2022_rsl_functions

! representation of Roughness SubLayer (RSL) corrections
! - RSL formulation should be independent from MOST formulation -- that is, we should
!   be able to combine any RSL functional dependence with any MOST functional dependence
! - solution of MOST for given sfc and atm state depends on RSL corrections
! - special case of no RSL correction should revert to original MOST solution precisely
! - neutral case can still use RSL correction, and therefore it could depart from log
!   profile assumed in current MOST code; perhaps a special case must be made for MOST
!   without RSL ("pure neutral" solution)
type, abstract :: rsl_functions_T
contains
  procedure(most_rsl_function), deferred :: rsl_m ! RSL correction for momentum
  procedure(most_rsl_function), deferred :: rsl_t ! RSL correction for heat and tracers
end type rsl_functions_T

abstract interface
  elemental real function most_rsl_function(this, zeta_r)
     import :: rsl_functions_T
     class(rsl_functions_T), intent(in) :: this
     real   , intent(in)  :: zeta_r  ! z/z_R, ratio of height to roughness sublayer height
  end function most_rsl_function
end interface



type, extends(rsl_functions_T) :: ridder2010_rsl_functions_T
  real :: mu_m = 2.59 ! exponent of RSL momentum correction
  real :: mu_t = 0.95 ! exponent of RSL heat and tracers correction
contains
  procedure :: rsl_m => ridder2010_rsl_m
  procedure :: rsl_t => ridder2010_rsl_t
end type ridder2010_rsl_functions_T

type, extends(rsl_functions_T) :: ghannam2022_rsl_functions_T
  ! 1 - mu_1*exp(-mu z/z_R)
  real :: mu_1 = 0.67 ! scaling of exponential correction
  real :: mu_m = 2.59 ! exponent of RSL momentum correction
  real :: mu_t = 0.95 ! exponent of RSL heat and tracers correction
contains
  procedure :: rsl_m => ghannam2022_rsl_m
  procedure :: rsl_t => ghannam2022_rsl_t
end type ghannam2022_rsl_functions_T

real, parameter :: RSL_UPPER_LIMIT = HUGE(1.0)*1e-4 ! actual upper limit when integrating RSL corrections to "infinity"
real, parameter :: RSL_RTOL        = 1e-8           ! relative tolerance for RSL integrals; setting it to 1e-7 or below results in significant artifacts

contains

! ---- Ridder (2010) RSL correction
function make_ridder2010_rsl_functions(mu_m,mu_t) result(ptr)
  class(ridder2010_rsl_functions_T), pointer :: ptr
  real, intent(in) :: mu_m, mu_t ! parameters or RSL correction for momentum and heat respectively

  allocate(ptr)
  ptr%mu_m = mu_m
  ptr%mu_t = mu_t
end function make_ridder2010_rsl_functions

elemental real function ridder2010_rsl_m(this, zeta_r)
   class(ridder2010_rsl_functions_T), intent(in)   :: this
   real, intent(in   )  :: zeta_r
   ridder2010_rsl_m = 1.0 - exp(-this%mu_m*zeta_r)
end function ridder2010_rsl_m

elemental real function ridder2010_rsl_t(this, zeta_r)
   class(ridder2010_rsl_functions_T), intent(in)   :: this
   real, intent(in   )  :: zeta_r
   ridder2010_rsl_t = 1.0 - exp(-this%mu_t*zeta_r)
end function ridder2010_rsl_t

! ---- Ghannam (2022) RSL corrections
function make_ghannam2022_rsl_functions(mu_1, mu_m, mu_t) result(ptr)
  class(ghannam2022_rsl_functions_T), pointer :: ptr
  real, intent(in) :: mu_1, mu_m, mu_t ! parameters or RSL correction for momentum and heat respectively

  allocate(ptr)
  ptr%mu_1 = mu_1
  ptr%mu_m = mu_m
  ptr%mu_t = mu_t
end function make_ghannam2022_rsl_functions

elemental real function ghannam2022_rsl_m(this, zeta_r)
   class(ghannam2022_rsl_functions_T), intent(in)   :: this
   real   , intent(in   )  :: zeta_r
   ghannam2022_rsl_m = 1.0 - this%mu_1 * exp(-this%mu_m*zeta_r)
end function ghannam2022_rsl_m

elemental real function ghannam2022_rsl_t(this, zeta_r)
   class(ghannam2022_rsl_functions_T), intent(in)   :: this
   real   , intent(in   )  :: zeta_r
   ghannam2022_rsl_t = 1.0 - this%mu_1 * exp(-this%mu_t*zeta_r)
end function ghannam2022_rsl_t

end module rsl_functions_mod

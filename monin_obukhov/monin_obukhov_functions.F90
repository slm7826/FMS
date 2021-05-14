module monin_obukhov_functions_mod
#include <fms_platform.h>

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
contains
  procedure(most_derivative_function), deferred :: derivative_m
  procedure(most_derivative_function), deferred :: derivative_t
  procedure(most_integral_m),          deferred :: integral_m
  procedure(most_integral_tq),         deferred :: integral_tq
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
  _PURE subroutine most_integral_m(this, n, mask, zeta, zeta_0, ln_z_z0, F_m, ier)
     import :: most_functions_T
     class(most_functions_T), intent(in)   :: this
     integer, intent(in   )                :: n
     logical, intent(in   ), dimension(n)  :: mask
     real   , intent(in   ), dimension(n)  :: zeta, zeta_0, ln_z_z0
     real   , intent(inout), dimension(n)  :: F_m
     integer, intent(  out)                :: ier
  end subroutine most_integral_m
  _PURE subroutine most_integral_tq(this, n, mask, zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq, F_t, F_q, ier)
     import :: most_functions_T
     class(most_functions_T), intent(in)   :: this
     integer, intent(in   )                :: n
     real   , intent(in   ), dimension(n)  :: zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq
     logical, intent(in   ), dimension(n)  :: mask
     real   , intent(inout), dimension(n)  :: F_t, F_q
     integer, intent(  out)                :: ier
  end subroutine most_integral_tq
end interface

type, extends(most_functions_T) :: neutral_functions_T
! stable option 1
contains
  procedure :: derivative_m => neutral_deriv_m
  procedure :: derivative_t => neutral_deriv_t
  procedure :: integral_m   => neutral_integral_m
  procedure :: integral_tq  => neutral_integral_tq
end type neutral_functions_T

type, extends(most_functions_T) :: most1_functions_T
! stable option 1
contains
  procedure :: derivative_m => most1_deriv_m
  procedure :: derivative_t => most1_deriv_t
  procedure :: integral_m   => most1_integral_m
  procedure :: integral_tq  => most1_integral_tq
end type most1_functions_T

type, extends(most_functions_T) :: most2_functions_T
! stable option 1
  real :: zeta_trans
contains
  procedure :: derivative_m => most2_deriv_m
  procedure :: derivative_t => most2_deriv_t
  procedure :: integral_m   => most2_integral_m
  procedure :: integral_tq  => most2_integral_tq
end type most2_functions_T

type, extends(most_functions_T) :: brutsaert_functions_T
  real :: a_s = 6.1,  b_s = 2.5 ! parameters of stable regime
  real :: a_u = 0.33, b_u = 0.41, c_u = 0.33, d_u = 0.057, n_u = 0.78 ! parameters of unstable regime
contains
  procedure :: derivative_m => brutsaert_deriv_m
  procedure :: derivative_t => brutsaert_deriv_t
  procedure :: integral_m   => brutsaert_integral_m
  procedure :: integral_tq  => brutsaert_integral_tq
end type brutsaert_functions_T

contains ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

_PURE subroutine neutral_integral_m(this, n, mask, zeta, zeta_0, ln_z_z0, F_m, ier)
  class(neutral_functions_T), intent(in) :: this
  integer, intent(in   )                :: n
  real   , intent(in   ), dimension(n)  :: zeta, zeta_0, ln_z_z0
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(inout), dimension(n)  :: F_m
  integer, intent(  out)                :: ier

  ier = 0
  F_m = ln_z_z0
end subroutine neutral_integral_m

_PURE subroutine neutral_integral_tq(this, n, mask, zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq, F_t, F_q, ier)
  class(neutral_functions_T), intent(in) :: this
  integer, intent(in   )                :: n
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(in   ), dimension(n)  :: zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq
  real   , intent(inout), dimension(n)  :: F_t, F_q
  integer, intent(  out)                :: ier

  ier = 0
  F_t = ln_z_zt
  F_q = ln_z_zq
end subroutine neutral_integral_tq

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

_PURE subroutine most1_integral_m(this, n, mask, zeta, zeta_0, ln_z_z0, F_m, ier)
  class(most1_functions_T), intent(in)     :: this
  integer, intent(in   )                :: n
  real   , intent(in   ), dimension(n)  :: zeta, zeta_0, ln_z_z0
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(inout), dimension(n)  :: F_m
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
     F_m  = ln_z_z0 - log(num/denom) + 2*y
  end where

  where (stable)
     F_m = ln_z_z0 + (5.0 - b_stab)*log((1.0 + zeta)/(1.0 + zeta_0)) &
          + b_stab*(zeta - zeta_0)
  end where
end subroutine most1_integral_m

_PURE subroutine most1_integral_tq(this, n, mask, zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq, F_t, F_q, ier)
  class(most1_functions_T), intent(in)    :: this
  integer, intent(in   )                :: n
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(in   ), dimension(n)  :: zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq
  real   , intent(inout), dimension(n)  :: F_t, F_q
  integer, intent(  out)                :: ier

  real, dimension(n)     :: x, x_t, x_q
  logical, dimension(n)  :: stable, unstable
  real                   :: b_stab

  ier = 0

  b_stab     = 1.0/this%rich_crit

  stable   = mask .and. zeta >= 0.0
  unstable = mask .and. zeta <  0.0

  where(unstable)
    x     = sqrt(1 - 16.0*zeta)
    x_t   = sqrt(1 - 16.0*zeta_t)
    x_q   = sqrt(1 - 16.0*zeta_q)

    F_t = ln_z_zt - 2.0*log( (1.0 + x)/(1.0 + x_t) )
    F_q = ln_z_zq - 2.0*log( (1.0 + x)/(1.0 + x_q) )
  end where
  where (stable)
    F_t = ln_z_zt + (5.0 - b_stab)*log((1.0 + zeta)/(1.0 + zeta_t)) &
       + b_stab*(zeta - zeta_t)
    F_q = ln_z_zq + (5.0 - b_stab)*log((1.0 + zeta)/(1.0 + zeta_q)) &
       + b_stab*(zeta - zeta_q)
  end where

end subroutine most1_integral_tq

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

_PURE subroutine most2_integral_m(this, n, mask, zeta, zeta_0, ln_z_z0, F_m, ier)
  class(most2_functions_T), intent(in)     :: this
  integer, intent(in   )                :: n
  real   , intent(in   ), dimension(n)  :: zeta, zeta_0, ln_z_z0
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(inout), dimension(n)  :: F_m
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
     F_m  = ln_z_z0 - log(num/denom) + 2*y
  end where

  lambda = 1.0 + (5.0 - b_stab)*this%zeta_trans

  weakly_stable   = stable .and. zeta <= this%zeta_trans
  strongly_stable = stable .and. zeta >  this%zeta_trans

  where (weakly_stable)
     F_m = ln_z_z0 + 5.0*(zeta - zeta_0)
  end where

  where(strongly_stable)
     x = (lambda - 1.0)*log(zeta/this%zeta_trans) + b_stab*(zeta - this%zeta_trans)
  end where

  where (strongly_stable .and. zeta_0 <= this%zeta_trans)
     F_m = ln_z_z0 + x + 5.0*(this%zeta_trans - zeta_0)
  end where
  where (strongly_stable .and. zeta_0 > this%zeta_trans)
     F_m = lambda*ln_z_z0 + b_stab*(zeta  - zeta_0)
  end where
end subroutine most2_integral_m

_PURE subroutine most2_integral_tq(this, n, mask, zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq, F_t, F_q, ier)
  class(most2_functions_T), intent(in)    :: this
  integer, intent(in   )                :: n
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(in   ), dimension(n)  :: zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq
  real   , intent(inout), dimension(n)  :: F_t, F_q
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
    x_t   = sqrt(1 - 16.0*zeta_t)
    x_q   = sqrt(1 - 16.0*zeta_q)

    F_t = ln_z_zt - 2.0*log( (1.0 + x)/(1.0 + x_t) )
    F_q = ln_z_zq - 2.0*log( (1.0 + x)/(1.0 + x_q) )
  end where

  lambda = 1.0 + (5.0 - b_stab)*this%zeta_trans

  weakly_stable   = stable .and. zeta <= this%zeta_trans
  strongly_stable = stable .and. zeta >  this%zeta_trans

  where (weakly_stable)
    F_t = ln_z_zt + 5.0*(zeta - zeta_t)
    F_q = ln_z_zq + 5.0*(zeta - zeta_q)
  end where

  where(strongly_stable)
    x = (lambda - 1.0)*log(zeta/this%zeta_trans) + b_stab*(zeta - this%zeta_trans)
  end where

  where (strongly_stable .and. zeta_t <= this%zeta_trans)
    F_t = ln_z_zt + x + 5.0*(this%zeta_trans - zeta_t)
  end where
  where (strongly_stable .and. zeta_t > this%zeta_trans)
    F_t = lambda*ln_z_zt + b_stab*(zeta  - zeta_t)
  end where

  where (strongly_stable .and. zeta_q <= this%zeta_trans)
    F_q = ln_z_zq + x + 5.0*(this%zeta_trans - zeta_q)
  end where
  where (strongly_stable .and. zeta_q > this%zeta_trans)
    F_q = lambda*ln_z_zq + b_stab*(zeta  - zeta_q)
  end where
end subroutine most2_integral_tq

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

  logical, dimension(n) :: stable, unstable
  real   , dimension(n) :: y

  ier = 0

  stable   = mask .and. zeta >= 0.0
  unstable = mask .and. zeta <  0.0

  where (unstable.and.y <= 1.0/this%b_u**3)
     y = -zeta
     phi = (this%a_u + this%b_u * y**(4.0/3.0))/(this%a_u + y)
  end where
  where (unstable.and.y > 1.0/this%b_u**3)
     phi = 1.0
  end where

  where (stable)
     phi = 1 + this%a_s * (zeta + zeta**this%b_s*(1+zeta**this%b_s)**(1.0/this%b_s-1))/ &
                          (zeta + (1+zeta**this%b_s)**(1.0/this%b_s))
  end where
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

_PURE subroutine brutsaert_integral_m(this, n, mask, zeta, zeta_0, ln_z_z0, F_m, ier)
  class(brutsaert_functions_T), intent(in)    :: this
  integer, intent(in   )                :: n
  real   , intent(in   ), dimension(n)  :: zeta, zeta_0, ln_z_z0
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(inout), dimension(n)  :: F_m
  integer, intent(  out)                :: ier

  ier = 0

  where (mask)
     F_m = ln_z_z0 - brutsaert_psi_m(this, zeta) + brutsaert_psi_m(this, zeta_0)
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

_PURE subroutine brutsaert_integral_tq(this, n, mask, zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq, F_t, F_q, ier)
  class(brutsaert_functions_T), intent(in)    :: this
  integer, intent(in   )                :: n
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(in   ), dimension(n)  :: zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq
  real   , intent(inout), dimension(n)  :: F_t, F_q
  integer, intent(  out)                :: ier

  where (mask)
     F_t = ln_z_zt - brutsaert_psi_h(this, zeta) + brutsaert_psi_h(this, zeta_t)
     F_q = ln_z_zq - brutsaert_psi_h(this, zeta) + brutsaert_psi_h(this, zeta_q)
  end where
end subroutine brutsaert_integral_tq

end module monin_obukhov_functions_mod

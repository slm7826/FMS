!> \file
!! \brief `integrate_mod` provides procedures for numerical calculation of definite integrals

module integrate_mod

#define _PURE

use, intrinsic :: ieee_arithmetic

implicit none
private

public :: integrand
public :: integrate_trapezoid
public :: integrate_midpoint
public :: integrate_simpson
public :: integrate_romberg_trapezoid
public :: integrate_romberg_midpoint
public :: integrate_romberg_midpoint_inv
public :: integrate_romberg_midpoint_exp


abstract interface
   _PURE real function integrand(x)
      real, intent(in) :: x
   end function integrand
   _PURE subroutine refiner(f,a,b,s,n)
      import :: integrand
      procedure(integrand)   :: f   ! function to integrate
      real   , intent(in)    :: a,b ! limits of integration
      real   , intent(inout) :: s   ! estimate of integral to be refined
      integer, intent(in)    :: n   ! number of integration steps
   end subroutine refiner
end interface

integer, parameter, public :: &
  NO_ERROR         = 0, &  ! success indicator
  RTOL_NOT_CHECKED = 1, &  ! relative tolerance of integral values was not checked (e.g. refinement degree too low)
  RTOL_NOT_REACHED = 2, &  ! relative tolerance of integral values checked, but not reached
  EVAL_ERROR       = 3     ! Error in evaluation of integral, e.g. function returned NaN or Inf

contains ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! ---------------------------------------------------------------------------------------
!> \brief Refine integral using trapezoidal rule
!!
!! Given a function f, limits of integration a,b, estimate of integral from the previous
!! stage of refinement s, and stage of refinement n, updates estimate of integral with higher
!! accuracy.
!!
!! Consecutive calls to `refine_trapezoid` will improve accuracy by adding 2**(n-2) interior
!! points.
!!
!! Note that this requires that the function is defined at the boundaries iof interval [a,b]
_PURE subroutine refine_trapezoid(f,a,b,s,n)
  procedure(integrand)   :: f   !< function to integrate
  real   , intent(in)    :: a   !< lower limit of integration
  real   , intent(in)    :: b   !< upper limit of integration
  real   , intent(inout) :: s   !< estimate of integral from the previous stage of refinement
  integer, intent(in)    :: n   !< stage of refinement. If n==1, value of s passed in this
                                !! subroutine is not used, and an estimate of integral
                                !! using trapezoid between a and b is returned in s.

  integer :: it ! number of extra steps
  integer :: i
  real    :: x, dx, sum
  if (n==1) then
     s = 0.5*(b-a)*(f(a)+f(b))
  else
     it = 2**(n-2) ! number of additional interior points
     dx = (b-a)/it
     sum = 0
     do i = 1,it
        x = a+(i-0.5)*dx
        sum = sum + f(x)
     enddo
     s = 0.5*(s+(b-a)*sum/it)
  endif
end subroutine refine_trapezoid

! ---------------------------------------------------------------------------------------
!> \brief Refine integral using midpoint rule
!!
!! Given a function f, limits of integration [a,b] an estimate of integral from the previous
!! stage of refinement s, and stage of refinement n, updates estimate of integral with higher
!! accuracy.
!!
!! Consecutive calls to `refine_midpoint` will improve accuracy by adding 3**(n-2) interior
!! points.
_PURE subroutine refine_midpoint(f,a,b,s,n)
  procedure(integrand)   :: f   !< function to integrate
  real   , intent(in)    :: a   !< lower limit of integration
  real   , intent(in)    :: b   !< upper limit of integration
  real   , intent(inout) :: s   !< estimate of integral from the previous stage of refinement
  integer, intent(in)    :: n   !< stage of refinement. If n==1, value of s passed in this
                                !! subroutine is not used, and an estimate of integral
                                !! using midpoint between a and b is returned in s.

  integer :: i, it
  real    :: x, dx, sum

  if (n==1) then
     s = (b-a)*f(0.5*(b+a))
  else
     it = 3**(n-2)
     dx = (b-a)/(3.0*it)
     x  = a+0.5*dx
     sum = 0.0
     do i = 1, it
       sum = sum+f(x)+f(x+2*dx)
       x   = x+3*dx
     enddo
     s = (s+(b-a)*sum/it)/3.0
  endif
end subroutine refine_midpoint

! ---------------------------------------------------------------------------------------
!> \brief Refine integral using midpoint rule with 1/x variable substitution
!!
!! Given a function f, limits of integration [a,b], an estimate of integral from the previous
!! stage of refinement s, and stage of refinement n, update estimate of integral with higher
!! accuracy. The function is evaluated at equally spaced points in 1/x rather than x. This
!! allows the upper limit to be as large and positive as the computer allows, or lower limit
!! to be large and negative (but not both). Both limits of integration must have the same sign.
!!
!! Consecutive calls to `refine_midpoint_inv` will improve accuracy by adding 3**(n-2) interior
!! points.
_PURE subroutine refine_midpoint_inv(ff,aa,bb,s,n)
  procedure(integrand)   :: ff  !< function to integrate
  real   , intent(in)    :: aa  !< lower limit of integration
  real   , intent(in)    :: bb  !< upper limit of integration
  real   , intent(inout) :: s   !< estimate of integral from the previous stage of refinement
  integer, intent(in)    :: n   !< stage of refinement. If n==1, value of s passed in this
                                !! subroutine is not used, and an estimate of integral
                                !! using midpoint between a and b is returned in s.

  integer :: i, it
  real    :: a, b, x, dx, sum

!   f(x)=ff(1.0/x)/x**2

  b = 1.0/aa
  a = 1.0/bb
  if (n==1) then
     s = (b-a)*f(0.5*(b+a))
  else
     it = 3**(n-2)
     dx = (b-a)/(3.0*it)
     x  = a+0.5*dx
     sum = 0.0
     do i = 1, it
       sum = sum+f(x)+f(x+2*dx)
       x   = x+3*dx
     enddo
     s = (s+(b-a)*sum/it)/3.0
  endif
!   write(*,*)'refine_midpoint_inv: n=',n,'s=',s
contains
  ! This function effects the change of variable.
  _PURE real function f(x)
     real, intent(in) :: x
     f = ff(1.0/x)/x**2
  end function f
end subroutine refine_midpoint_inv

! ---------------------------------------------------------------------------------------
!> \brief Refine integral using midpoint rule with exp(-x) variable substitution
!!
!! Given a function f, limits of integration [a,b], an estimate of integral from the previous
!! stage of refinement s, and stage of refinement n, update estimate of integral with higher
!! accuracy. The function is evaluated at equally spaced points in exp(-x) rather than x. This
!! allows the upper limit to be as large and positive as the computer allows. Both limits
!! must have the same sign.
!!
!! Consecutive calls to `refine_midpoint_exp` will improve accuracy by adding 3**(n-2) interior
!! points.
_PURE subroutine refine_midpoint_exp(ff,aa,bb,s,n)
  procedure(integrand)   :: ff  !< function to integrate
  real   , intent(in)    :: aa  !< lower limit of integration
  real   , intent(in)    :: bb  !< upper limit of integration
  real   , intent(inout) :: s   !< estimate of integral from the previous stage of refinement
  integer, intent(in)    :: n   !< stage of refinement. If n==1, value of s passed in this
                                !! subroutine is not used, and an estimate of integral
                                !! using midpoint between a and b is returned in s.

  integer :: i, it
  real    :: a, b, x, dx, sum

  b = exp(-aa)
  a = exp(-bb)
  if (n==1) then
     s = (b-a)*f(0.5*(b+a))
  else
     it = 3**(n-2)
     dx = (b-a)/(3.0*it)
     x  = a+0.5*dx
     sum = 0.0
     do i = 1, it
       sum = sum+f(x)+f(x+2*dx)
       x   = x+3*dx
     enddo
     s = (s+(b-a)*sum/it)/3.0
  endif
contains
  ! This function effects the change of variable.
  _PURE real function f(x)
     real, intent(in) :: x
     f = ff(-log(x))/x
  end function f
end subroutine refine_midpoint_exp


! ---------------------------------------------------------------------------------------
!> \brief Integrate function numerically within a given interval using trapezoidal rule
!!
subroutine integrate_trapezoid(f,a,b, rtol, s, ierr, maxD, nDeg, nSteps)
  procedure(integrand)   :: f    !< function to integrate
  real   , intent(in)    :: a    !< lower limit of integration
  real   , intent(in)    :: b    !< upper limit of integration
  real   , intent(in)    :: rtol !< relative tolerance of the integral
  real   , intent(out)   :: s    !< estimate of integral
  integer, intent(out)   :: ierr !< indication of error (0 means no error)
  integer, intent(in),  optional :: maxD    !< maximum degree of refinement
                                 !! corresponding maximum number of steps is 2**(maxD-1)
  integer, intent(out), optional :: nDeg    !< degree of refinement reached
  integer, intent(out), optional :: nSteps  !< number of steps reached

  integer, parameter :: minD = 5 ! minimum degree of refinement

  real :: s0
  integer :: i
  integer :: maxD_

  maxD_ = 20
  if (present(maxD)) maxD_ = maxD

  ierr = RTOL_NOT_CHECKED
  s0 = -HUGE(1.0)
  do i = 1, maxD_
     call refine_trapezoid(f,a,b,s,i)
     !! Note that for degree of refinement <= minD the convergence to relative tolerance
     !! is never checked, so that the returned ierr is always RTOL_NOT_REACHED
     if(i > minD) then
        ierr = RTOL_NOT_REACHED
        if (abs(s-s0)<rtol*abs(s0).or.(s==0.and.s0==0)) then
           ierr = NO_ERROR; exit
        endif
     endif
     s0 = s
  enddo
  if (present(nDeg)) nDeg = i
  if (present(nSteps)) nSteps = 2**(i-1)
end subroutine integrate_trapezoid

! ---------------------------------------------------------------------------------------
!> \brief Integrate function numerically within a given interval using midpoint rule
subroutine integrate_midpoint(f,a,b, rtol, s, ierr, maxD, nDeg, nSteps)
  procedure(integrand)   :: f    !< function to integrate
  real   , intent(in)    :: a    !< lower limit of integration
  real   , intent(in)    :: b    !< upper limit of integration
  real   , intent(in)    :: rtol !< relative tolerance of the integral
  real   , intent(out)   :: s    !< estimate of integral
  integer, intent(out)   :: ierr !< indication of error (0 means no error)
  integer, intent(in), optional :: maxD !< maximum degree of refinement;
                                 !! corresponding maximum number of steps = 3**(maxD-1)
  integer, intent(out), optional :: nDeg    !< degree of refinement reached
  integer, intent(out), optional :: nSteps  !< number of steps reached

  integer, parameter :: minD = 5 ! minimum degree of refinement

  real :: s0
  integer :: i
  integer :: maxD_

  maxD_ = 20
  if (present(maxD)) maxD_ = maxD

  ierr = RTOL_NOT_CHECKED
  s0 = -HUGE(1.0)
  do i = 1, maxD_
     call refine_midpoint(f,a,b,s,i)
     if(i > minD) then
        ierr = RTOL_NOT_REACHED
        if (abs(s-s0)<rtol*abs(s0).or.(s==0.and.s0==0)) then
           ierr = NO_ERROR; exit
        endif
     endif
     s0 = s
  enddo
  if (present(nDeg)) nDeg = i
  if (present(nSteps)) nSteps = 3**(i-1)
end subroutine integrate_midpoint

! ---------------------------------------------------------------------------------------
!> \brief Integrate function numerically within a given interval using Simpson rule
subroutine integrate_simpson(f,a,b, rtol, s, ierr, maxD, nDeg, nSteps)
  procedure(integrand)   :: f    !< function to integrate
  real   , intent(in)    :: a    !< lower limit of integration
  real   , intent(in)    :: b    !< upper limit of integration
  real   , intent(in)    :: rtol !< relative tolerance
  real   , intent(out)   :: s    !< estimate of integral
  integer, intent(out)   :: ierr !< indication of error (0 means no error)
  integer, intent(in), optional :: maxD !< maximum degree of refinement;
                                 !! corresponding maximum number of steps = 2**(maxd-1)
  integer, intent(out), optional :: nDeg    !< degree of refinement reached
  integer, intent(out), optional :: nSteps  !< number of steps reached

  integer, parameter :: minD = 5 ! minimum degree of refinement

  real :: os, ost, st
  integer :: i
  integer :: maxD_

  maxD_ = 20
  if (present(maxD)) maxD_ = maxD

  ierr = RTOL_NOT_CHECKED
  ost = -HUGE(1.0)
  os  = -HUGE(1.0)
  do i = 1, maxD_
     call refine_trapezoid(f,a,b,st,i)
     s=(4.0*st-ost)/3.0
     if(i > minD) then
        ierr = RTOL_NOT_REACHED
        if (abs(s-os)<rtol*abs(os).or.(s==0.and.os==0)) then
           ierr = NO_ERROR; exit
        endif
     endif
     os  = s
     ost = st
  enddo
  if (present(nDeg)) nDeg = i
  if (present(nSteps)) nSteps = 2**(i-1)
end subroutine integrate_simpson

! ---------------------------------------------------------------------------------------
!> Calculate integral of a given function using Romberg procedure with trapezoidal integration rule
subroutine integrate_romberg_trapezoid(f,a,b, rtol, ss, ierr, maxD, nDeg, nSteps)
  procedure(integrand)   :: f    !< function to integrate
  real   , intent(in)    :: a    !< lower limit of integration
  real   , intent(in)    :: b    !< upper limit of integration
  real   , intent(in)    :: rtol !<relative tolerance
  real   , intent(out)   :: ss   !< estimate of integral
  integer, intent(out)   :: ierr !< indication of error (0 means no error)
  integer, intent(in), optional :: maxD !< maximum degree of refinement;
                                 !! corresponding maximum number of steps is 2**(maxd-1)
  integer, intent(out), optional :: nDeg    !< degree of refinement reached
  integer, intent(out), optional :: nSteps  !< number of steps reached

  integer, parameter :: K = 5, KM = K-1

  integer :: i, maxD_
  real :: dss
  real, allocatable :: h(:),s(:)

  ierr = RTOL_NOT_CHECKED
  maxD_ = 20
  if (present(maxD)) maxD_ = maxD
  allocate(h(maxD_+1),s(maxD_+1))
  h(1) = 1.0
  do i = 1, maxD_
     call refine_trapezoid(f,a,b,s(i),i)
     if(i >= K) then
        ierr = RTOL_NOT_REACHED
        call polint(h(i-KM:),s(i-KM:),K,0.0,ss,dss)
        if (abs(dss) <= rtol*abs(ss)) then
           ierr = NO_ERROR; exit
        endif
     endif
     s(i+1) = s(i)
     h(i+1) = h(i)/4
  enddo
  deallocate(h,s)
  if (present(nDeg)) nDeg = i
  if (present(nSteps)) nSteps = 2**(i-1)
end subroutine integrate_romberg_trapezoid

! ---------------------------------------------------------------------------------------
!> Calculate integral of a given function using Romberg procedure with midpoint integration rule
_PURE subroutine integrate_romberg3(f, a, b, refine, rtol, ss, ierr, maxD, nDeg, nSteps)
  procedure(integrand)   :: f       !< function to integrate
  real   , intent(in)    :: a       !< lower limit of integration
  real   , intent(in)    :: b       !< upper limit of integration
  procedure(refiner)     :: refine  !< subroutine that refines integral value, one of refine_midpoint, refine_midpoint_inv, refine_midpoint_exp
  real   , intent(in)    :: rtol    !< relative tolerance
  real   , intent(out)   :: ss      !< estimate of integral
  integer, intent(out)   :: ierr    !< indication of error (0 means no error)
  integer, intent(in), optional :: maxD !< maximum degree of refinement;
                                 !! corresponding maximum number of steps = 3**(maxD-1)
  integer, intent(out), optional :: nDeg    !< degree of refinement reached
  integer, intent(out), optional :: nSteps  !< number of steps reached

  integer, parameter :: K = 5, KM = K-1

  integer :: i, maxD_
  real :: dss
  real, allocatable :: h(:),s(:)

  ierr = RTOL_NOT_CHECKED
  maxD_ = 20
  if (present(maxD)) maxD_ = maxD
  allocate(h(maxD_+1),s(maxD_+1))
  h(1) = 1.0; s(1) = 0.0
  do i = 1, maxD_
     call refine(f,a,b,s(i),i)
!      write(*,*) i,s(i)
     if(.not.ieee_is_finite(s(i)))then
        ierr = EVAL_ERROR; exit
     endif
     if(i >= K) then
        ierr = RTOL_NOT_REACHED
        call polint(h(i-KM:),s(i-KM:),K,0.0,ss,dss)
        if (abs(dss) <= rtol*abs(ss)) then
           ierr = NO_ERROR; exit
        endif
     endif
     s(i+1) = s(i)
     h(i+1) = h(i)/9 ! this is where the assumption of step tripling is used
  enddo
  deallocate(h,s)
  if (present(nDeg)) nDeg = i
  if (present(nSteps)) nSteps = 3**(i-1)
end subroutine integrate_romberg3

! ---------------------------------------------------------------------------------------
!> Calculate integral of a given function using Romberg procedure with midpoint integration rule
_PURE subroutine integrate_romberg_midpoint(f, a, b, rtol, s, ierr, maxD, nDeg, nSteps)
  procedure(integrand)   :: f    !< function to integrate
  real   , intent(in)    :: a    !< lower limit of integration
  real   , intent(in)    :: b    !< upper limit of integration
  real   , intent(in)    :: rtol !< relative tolerance
  real   , intent(out)   :: s    !< estimate of integral
  integer, intent(out)   :: ierr !< indication of error (0 means no error)
  integer, intent(in), optional :: maxD !< maximum degree of refinement;
                                 !! corresponding maximum number of steps = 3**(maxD-1)
  integer, intent(out), optional :: nDeg    !< degree of refinement reached
  integer, intent(out), optional :: nSteps  !< number of steps reached

  call integrate_romberg3(f, a, b, refine_midpoint, rtol, s, ierr, maxD, nDeg, nSteps)
end subroutine integrate_romberg_midpoint

! ---------------------------------------------------------------------------------------
!> Calculate integral of a given function using Romberg procedure with midpoint integration rule
_PURE subroutine integrate_romberg_midpoint_inv(f, a, b, rtol, s, ierr, maxD, nDeg, nSteps)
  procedure(integrand)   :: f    !< function to integrate
  real   , intent(in)    :: a    !< lower limit of integration
  real   , intent(in)    :: b    !< upper limit of integration
  real   , intent(in)    :: rtol !< relative tolerance
  real   , intent(out)   :: s    !< estimate of integral
  integer, intent(out)   :: ierr !< indication of error (0 means no error)
  integer, intent(in), optional :: maxD !< maximum degree of refinement;
                                 !! corresponding maximum number of steps = 3**(maxD-1)
  integer, intent(out), optional :: nDeg    !< degree of refinement reached
  integer, intent(out), optional :: nSteps  !< number of steps reached

  call integrate_romberg3(f, a, b, refine_midpoint_inv, rtol, s, ierr, maxD, nDeg, nSteps)
end subroutine integrate_romberg_midpoint_inv

! ---------------------------------------------------------------------------------------
!> Calculate integral of a given function using Romberg procedure with midpoint integration rule
subroutine integrate_romberg_midpoint_exp(f, a, b, rtol, s, ierr, maxD, nDeg, nSteps)
  procedure(integrand)   :: f    !< function to integrate
  real   , intent(in)    :: a    !< lower limit of integration
  real   , intent(in)    :: b    !< upper limit of integration
  real   , intent(in)    :: rtol !< relative tolerance
  real   , intent(out)   :: s    !< estimate of integral
  integer, intent(out)   :: ierr !< indication of error (0 means no error)
  integer, intent(in), optional :: maxD !< maximum degree of refinement;
                                 !! corresponding maximum number of steps = 3**(maxD-1)
  integer, intent(out), optional :: nDeg    !< degree of refinement reached
  integer, intent(out), optional :: nSteps  !< number of steps reached

  call integrate_romberg3(f, a, b, refine_midpoint_exp, rtol, s, ierr, maxD, nDeg, nSteps)
end subroutine integrate_romberg_midpoint_exp

! ---------------------------------------------------------------------------------------
!> \brief Polinomial interpolation
!!
!! Given arrays xa and ya, each of length n, and given a value x, this routine returns
!! a value y, and an error estimate dy. If P (x) is the polynomial of degree N − 1 such
!! that P(xa(i)) = ya(i),i = 1,...,n, then the returned value y = P(x)
_PURE subroutine polint(xa,ya,n,x,y,dy)
  real,    intent(in)  :: xa(:) !< array of x-coordinates
  real,    intent(in)  :: ya(:) !< array of y-coordinates
  integer, intent(in)  :: n     !< size of points in the arrays, degree of interpolation
  real,    intent(in)  :: x     !< point to interpolate t
  real,    intent(out) :: y     !< result of interpolation
  real,    intent(out) :: dy    !< error estimate

  integer, parameter :: NMAX=10 ! Largest anticipated value of n.

  integer :: i,m,ns
  real    :: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

  ns=1
  dif=abs(x-xa(1))

  ! Here we find the index ns of the closest table entry,
  ! and initialize the tableau of c’s and d’s.
  do i=1,n
     dift=abs(x-xa(i))
     if (dift<dif) then
        ns=i
        dif=dift
    endif
    c(i)=ya(i)
    d(i)=ya(i)
  enddo
  y=ya(ns) ! This is the initial approximation to y.
  ns=ns-1
  ! For each column of the tableau,
  ! we loop over the current c’s and d’s and update them.
  do m=1,n-1
     do i=1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        ! if(den==0.0) pause
        ! This error can occur only if two input xa’s are (to within roundoff) identical.
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
     enddo
     if (2*ns < n-m) then
        dy=c(ns+1)
     else
        dy=d(ns)
        ns=ns-1
     endif
     y=y+dy
  enddo

end subroutine polint


end module integrate_mod

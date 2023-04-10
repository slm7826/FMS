!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

module monin_obukhov_mod


!=======================================================================
!
!                         MONIN-OBUKHOV MODULE
!
!          Routines for computing surface drag coefficients
!                 from data at the lowest model level
!              and for computing the profile of fields
!           between the lowest model level and the ground
!                  using Monin-Obukhov scaling
!
!=======================================================================


use constants_mod, only: grav, vonkarm
use mpp_mod,       only: input_nml_file
use fms_mod,       only: error_mesg, FATAL, NOTE, file_exist,   &
                         check_nml_error, open_namelist_file,      &
                         mpp_pe, mpp_root_pe, close_file, stdlog, &
                         write_version_number, lowercase, string
use monin_obukhov_functions_mod, only: most_functions_T, &
                         make_most1_functions, make_most2_functions, &
                         make_brutsaert_functions, make_neutral_functions
use rsl_functions_mod, only: rsl_functions_T, &
                         make_ridder2010_rsl_functions, make_ghannam2022_rsl_functions
use monin_obukhov_kernel, only: monin_obukhov_diff, monin_obukhov_drag_1d, &
                         monin_obukhov_solve_zeta, monin_obukhov_profile_1d

implicit none
private

!=======================================================================
 public monin_obukhov_init
 public monin_obukhov_end
 public mo_drag
 public mo_profile
 public mo_diff
 public stable_mix
!=======================================================================

interface mo_drag
    module procedure  mo_drag_0d, mo_drag_1d, mo_drag_2d
end interface


interface mo_profile
    module procedure  mo_profile_0d,   mo_profile_1d,   mo_profile_2d, &
                      mo_profile_0d_n, mo_profile_1d_n, mo_profile_2d_n
end interface

interface mo_diff
    module procedure  mo_diff_0d_n, mo_diff_0d_1, &
                      mo_diff_1d_n, mo_diff_1d_1, &
                      mo_diff_2d_n, mo_diff_2d_1
end interface

!--------------------- version number ---------------------------------
! Include variable "version" to be written to log file.
#include <file_version.h>

!=======================================================================
!  DEFAULT VALUES OF NAMELIST PARAMETERS:

real    :: rich_crit      = 2.0
real    :: drag_min_heat  = 1.e-05
real    :: drag_min_moist = 1.e-05
real    :: drag_min_mom   = 1.e-05
character(32) :: stable_option  = '1' ! valid options: "1", "2", "brutsaert", or "neutral"
real    :: zeta_trans     = 0.5

character(32) :: rsl_option = 'none'
real    :: rsl_mu_1 = 0.67 ! parameter of RSL correction
real    :: rsl_mu_m = 2.59 ! parameter of RSL momentum correction
real    :: rsl_mu_t = 0.95 ! parameter of RSL heat and tracer correction
! parameters of RSL integrals Im and It lookup tables
logical :: use_RSL_lookup = .TRUE. !> use loookup tables to compute RSL integrals; otherwise
                              !! calculate integrals directly: this can be used for, say,
                              !! testing of quality of lookup in a single point runs, but would
                              !! very likely be prohibitively slow in global simulations
real    :: a_min    = 1e-5    !> lower lookup table limit for parameter a of I_m and I_h RSL integrals: a_min > 0.
real    :: a_max    = 100     !> upper lookup table limit for parameter a of I_m and I_h RSL integrals: a_max > a_min > 0.
integer :: a_nsteps = 100     !> number of lookup table steps along the axis a.
real    :: b_min    = -10.0   !> lower lookup table limit for parameter b of I_m and I_h RSL integrals.
real    :: b_max    =  1000.0 !> upper lookup table limit for parameter b of I_m and I_h RSL integrals
integer :: b_nsteps = 100     !> number of lookup table steps along the axis b.

namelist /monin_obukhov_nml/ rich_crit, drag_min_heat, drag_min_moist, drag_min_mom, &
                             stable_option, zeta_trans, & !miz
                             rsl_option, rsl_mu_1, rsl_mu_m, rsl_mu_t, &
                             use_RSL_lookup, a_min, a_max, a_nsteps, b_min, b_max, b_nsteps


!=======================================================================
!  MODULE VARIABLES

logical            :: module_is_initialized = .false.
class(most_functions_T), pointer :: most => NULL() ! pointer to an object representing
                   ! stability correction functions

contains

!=======================================================================

subroutine monin_obukhov_init

integer :: unit, ierr, io, logunit
class(rsl_functions_T),  pointer :: rsl => NULL()  ! pointer to roughness sublayer (RSL)
                   ! correction functions

!------------------- read namelist input -------------------------------

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=monin_obukhov_nml, iostat=io)
      ierr = check_nml_error(io,"monin_obukhov_nml")
#else
      if (file_exist('input.nml')) then
         unit = open_namelist_file ()
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=monin_obukhov_nml, iostat=io, end=10)
            ierr = check_nml_error(io,'monin_obukhov_nml')
         enddo
  10     call close_file (unit)
      endif
#endif

!---------- output namelist to log-------------------------------------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number('monin_obukhov_nml', version)
           logunit = stdlog()
           write (logunit, nml=monin_obukhov_nml)
      endif

!----------------------------------------------------------------------

if(rich_crit.le.0.25)  call error_mesg( &
        'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', &
        'rich_crit in monin_obukhov_mod must be > 0.25', FATAL)

if(drag_min_heat.le.0.0)  call error_mesg( &
        'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', &
        'drag_min_heat in monin_obukhov_mod must be >= 0.0', FATAL)

if(drag_min_moist.le.0.0)  call error_mesg( &
        'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', &
        'drag_min_moist in monin_obukhov_mod must be >= 0.0', FATAL)

if(drag_min_mom.le.0.0)  call error_mesg( &
        'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', &
        'drag_min_mom in monin_obukhov_mod must be >= 0.0', FATAL)

select case(trim(lowercase(stable_option)))
case('1')
   most => make_most1_functions(rich_crit)
case('2')
   if(zeta_trans < 0) call error_mesg( &
        'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', &
        'zeta_trans must be positive', FATAL)
   most => make_most2_functions(rich_crit, zeta_trans)
case('brutsaert')
   most => make_brutsaert_functions(rich_crit)
case('neutral')
   most => make_neutral_functions(rich_crit)
case default
   call error_mesg( &
        'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', &
        'stable_option = "'//trim(stable_option)//'" is incorrect, use "1", "2", "brutsaert", or "neutral"', FATAL)
end select

! set up roughness sublayer (RSL) corrections
select case(trim(lowercase(rsl_option)))
case('none')
   rsl=>NULL()
case('ridder2010')
   rsl=>make_ridder2010_rsl_functions(rsl_mu_m,rsl_mu_t)
case('ghannam2022')
   rsl=>make_ghannam2022_rsl_functions(rsl_mu_1,rsl_mu_m,rsl_mu_t)
case default
   call error_mesg( &
      'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', &
      'rsl_option = "'//trim(rsl_option)//'" is incorrect, use "none", "ghannam2022", or "ridder2010"', FATAL)
end select

call error_mesg('MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', 'Will set up RSL functions', NOTE)
call most%set_rsl_functions(rsl,use_RSL_lookup,a_min,a_max,a_nsteps,b_min,b_max,b_nsteps)
call error_mesg('MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', 'Did set up RSL functions', NOTE)

module_is_initialized = .true.

end subroutine monin_obukhov_init

!=======================================================================
subroutine monin_obukhov_end

if (associated(most)) deallocate(most)
module_is_initialized = .false.

end subroutine monin_obukhov_end

!=======================================================================
subroutine mo_drag_1d &
         (pt, pt0, z, z0, zt, zq, zR, speed, drag_m, drag_t, drag_q, &
          u_star, b_star, avail)
  real, intent(in)   , dimension(:) :: pt, pt0, z, z0, zt, zq, speed
  real, intent(in)   , dimension(:) :: zR ! roughness sublayer scale
  real, intent(inout), dimension(:) :: drag_m, drag_t, drag_q, u_star, b_star
  logical, intent(in), optional, dimension(:) :: avail

  integer            :: n, ier

  integer, parameter :: max_iter = 20
  real   , parameter :: error=1.e-04, zeta_min=1.e-06, small=1.e-04

  real   , dimension(size(pt)) :: rich, zeta

  if(.not.module_is_initialized) call error_mesg('mo_drag_1d in monin_obukhov_mod', &
       'monin_obukhov_init has not been called', FATAL)

  if (present(avail)) then
     if (count(avail) .eq. 0) return
  endif
  n = size(pt)
  call monin_obukhov_drag_1d(most, grav, vonkarm,                  &
       & error, zeta_min, max_iter, small,                         &
       & drag_min_heat, drag_min_moist, drag_min_mom,              &
       & n, pt, pt0, z, z0, zt, zq, zR, speed, drag_m, drag_t,         &
       & drag_q, u_star, b_star, rich, zeta, ier, avail)

end subroutine mo_drag_1d

!=======================================================================
subroutine mo_profile_1d(zref, zref_t, z, z0, zt, zq, zR, u_star, b_star, q_star, &
                         del_m, del_t, del_q, avail)
  real,    intent(in)                :: zref, zref_t
  real,    intent(in) , dimension(:) :: z, z0, zt, zq, u_star, b_star, q_star
  real,    intent(in) , dimension(:) :: zR ! roughness sublayer length scale
  real,    intent(out), dimension(:) :: del_m, del_t, del_q
  logical, intent(in) , optional, dimension(:) :: avail

  integer                            :: n, ier

  if(.not.module_is_initialized) call error_mesg('mo_profile_1d in monin_obukhov_mod', &
       'monin_obukhov_init has not been called', FATAL)

  if (present(avail)) then
     if (count(avail) .eq. 0) return
  endif
  n = size(z)
  call monin_obukhov_profile_1d(most, vonkarm, &
       & n, zref, zref_t, z, z0, zt, zq, zR, u_star, b_star, q_star, &
       & del_m, del_t, del_q, ier, avail)

end subroutine mo_profile_1d

!=======================================================================
subroutine stable_mix(rich, z_ag, rsl_scale, mix)

real, intent(in) , dimension(:,:,:)  :: rich      ! Richardson number
real, intent(in) , dimension(:,:,:)  :: z_ag      ! height above ground, m
real, intent(in) , dimension(:,:)    :: rsl_scale ! roughness sublayer scale, m
real, intent(out), dimension(:,:,:)  :: mix       !

integer :: i,j,k,n
integer :: ier ! error code returned by most%stable_mix
real, dimension(size(rich,1),size(rich,2),size(rich,3)) :: &
    pm2, & ! RSL correction for momentum squared
    Ri     ! Richardson number scaled with RSL corrections
real :: pt ! RSL correction for heat

if(.not.module_is_initialized) call error_mesg('stable_mix_3d in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

if(size(rich,3).ne.size(z_ag,3)) call error_mesg('stable_mix_3d in monin_obukhov_mod', &
     'vertical sizes of "rich" ('//string(size(rich,3))//') and "z_ag" (' &
     //string(size(z_ag,3))//') are inconsistent', FATAL)

n = size(rich,1)*size(rich,2)*size(rich,3)

if(associated(most%rsl)) then
  ! scale Richardson number with roughness sublayer corrections, where necessary
  do j = 1,size(rich,2)
  do i = 1,size(rich,1)
    if (rsl_scale(i,j)>0.0) then
      do k = 1,size(rich,3)
        pm2(i,j,k) = most%rsl%rsl_m(z_ag(i,j,k)/rsl_scale(i,j))**2
        pt         = most%rsl%rsl_t(z_ag(i,j,k)/rsl_scale(i,j))
        Ri (i,j,k) = rich(i,j,k) * pm2(i,j,k)/pt
      enddo
    else
      pm2(i,j,:) = 1.0
      Ri (i,j,:) = rich(i,j,:)
    endif
  enddo
  enddo

  call most%stable_mix(n, Ri, mix, ier)

  ! scale mixing factor with roughness sublayer correction
  mix(:,:,:) = mix(:,:,:)/pm2(:,:,:)
else
  call most%stable_mix(n, rich, mix, ier)
endif

if (ier.ne.0) call error_mesg('stable_mix_3d in monin_obukhov_mod', &
     'stable_mix calculations for stable_option "'//trim(stable_option)//'" returned an error', FATAL)

end subroutine stable_mix

!=======================================================================
subroutine mo_diff_2d_n(z, u_star, b_star, zR, k_m, k_h)

real, intent(in),  dimension(:,:,:) :: z
real, intent(in),  dimension(:,:)   :: u_star, b_star, zR
real, intent(out), dimension(:,:,:) :: k_m, k_h

integer            :: ni, nj, nk, ier
real, parameter    :: ustar_min = 1.e-10

if(.not.module_is_initialized) call error_mesg('mo_diff_2d_n in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

ni = size(z, 1); nj = size(z, 2); nk = size(z, 3)
call monin_obukhov_diff(most, vonkarm,                           &
          & ustar_min,                                     &
          & ni, nj, nk, z, u_star, b_star, zR, k_m, k_h, ier)

end subroutine mo_diff_2d_n

!=======================================================================
! The following routines allow the public interfaces to be used
! with different dimensions of the input and output
!
!=======================================================================
subroutine mo_drag_2d &
    (pt, pt0, z, z0, zt, zq, zR, speed, drag_m, drag_t, drag_q, u_star, b_star)

real, intent(in)   , dimension(:,:) :: z, speed, pt, pt0, z0, zt, zq, zR
real, intent(out)  , dimension(:,:) :: drag_m, drag_t, drag_q
real, intent(inout), dimension(:,:) :: u_star, b_star

integer :: j

do j = 1, size(pt,2)
  call mo_drag_1d (pt(:,j), pt0(:,j), z(:,j), z0(:,j), zt(:,j), zq(:,j), zR(:,j),&
                   speed(:,j), drag_m(:,j), drag_t(:,j), drag_q(:,j), &
                   u_star(:,j), b_star(:,j))
end do

end subroutine mo_drag_2d

!=======================================================================
subroutine mo_drag_0d &
    (pt, pt0, z, z0, zt, zq, zR, speed, drag_m, drag_t, drag_q, u_star, b_star)

real, intent(in)    :: z, speed, pt, pt0, z0, zt, zq, zR
real, intent(out)   :: drag_m, drag_t, drag_q, u_star, b_star

real, dimension(1) :: pt_1, pt0_1, z_1, z0_1, zt_1, zq_1, zR_1, speed_1, &
                      drag_m_1, drag_t_1, drag_q_1, u_star_1, b_star_1

pt_1   (1) = pt
pt0_1  (1) = pt0
z_1    (1) = z
z0_1   (1) = z0
zt_1   (1) = zt
zq_1   (1) = zq
zR_1   (1) = zR
speed_1(1) = speed

call mo_drag_1d (pt_1, pt0_1, z_1, z0_1, zt_1, zq_1, zR_1, speed_1, &
                 drag_m_1, drag_t_1, drag_q_1, u_star_1, b_star_1)

drag_m = drag_m_1(1)
drag_t = drag_t_1(1)
drag_q = drag_q_1(1)
u_star = u_star_1(1)
b_star = b_star_1(1)

end subroutine mo_drag_0d

!=======================================================================
subroutine mo_profile_2d(zref, zref_t, z, z0, zt, zq, zR, u_star, b_star, q_star, &
                         del_m, del_h, del_q)

real, intent(in)                  :: zref, zref_t
real, intent(in) , dimension(:,:) :: z, z0, zt, zq, u_star, b_star, q_star
real, intent(in) , dimension(:,:) :: zR ! roughness sublayer length scale
real, intent(out), dimension(:,:) :: del_m, del_h, del_q

integer :: j

do j = 1, size(z,2)
  call mo_profile_1d (zref, zref_t, z(:,j), z0(:,j), zt(:,j),         &
                      zq(:,j), zR(:,j), u_star(:,j), b_star(:,j), q_star(:,j), &
                      del_m(:,j), del_h (:,j), del_q (:,j))
enddo

end subroutine mo_profile_2d

!=======================================================================
subroutine mo_profile_0d(zref, zref_t, z, z0, zt, zq, zR, u_star, b_star, q_star, &
                         del_m, del_h, del_q)

real, intent(in)  :: zref, zref_t
real, intent(in)  :: z, z0, zt, zq, u_star, b_star, q_star
real, intent(in)  :: zR ! roughness sublayer length scale
real, intent(out) :: del_m, del_h, del_q

real, dimension(1) :: z_1, z0_1, zt_1, zq_1, zR_1, u_star_1, b_star_1, q_star_1, &
                      del_m_1, del_h_1, del_q_1

z_1     (1) = z
z0_1    (1) = z0
zt_1    (1) = zt
zq_1    (1) = zq
zR_1    (1) = zR
u_star_1(1) = u_star
b_star_1(1) = b_star
q_star_1(1) = q_star

call mo_profile_1d (zref, zref_t, z_1, z0_1, zt_1, zq_1, zR_1, &
                    u_star_1, b_star_1, q_star_1,              &
                    del_m_1, del_h_1, del_q_1)

del_m = del_m_1(1)
del_h = del_h_1(1)
del_q = del_q_1(1)

end subroutine mo_profile_0d

!=======================================================================
subroutine mo_profile_1d_n(zref, z, z0, zt, zq, zR, u_star, b_star, q_star, &
                         del_m, del_t, del_q, avail)
  real,    intent(in),  dimension(:)   :: zref
  real,    intent(in) , dimension(:)   :: z, z0, zt, zq, u_star, b_star, q_star
  real,    intent(in) , dimension(:)   :: zR ! roughness sublayer length scale
  real,    intent(out), dimension(:,:) :: del_m, del_t, del_q
  logical, intent(in) , optional, dimension(:) :: avail

  integer :: k

  do k = 1, size(zref(:))
     call mo_profile_1d (zref(k), zref(k), z, z0, zt, zq, zR, &
        u_star, b_star, q_star, del_m(:,k), del_t(:,k), del_q(:,k), avail)
  enddo
end subroutine mo_profile_1d_n

!=======================================================================
subroutine mo_profile_0d_n(zref, z, z0, zt, zq, zR, u_star, b_star, q_star, &
                         del_m, del_t, del_q)

real,    intent(in),  dimension(:) :: zref
real,    intent(in)                :: z, z0, zt, zq, u_star, b_star, q_star
real,    intent(in)                :: zR ! roughness sublayer length scale
real,    intent(out), dimension(:) :: del_m, del_t, del_q

integer :: k

do k = 1, size(zref(:))
  call mo_profile_0d (zref(k), zref(k), z, z0, zt, zq, zR, &
       u_star, b_star, q_star, del_m(k), del_t(k), del_q(k))
enddo

end subroutine mo_profile_0d_n

!=======================================================================
subroutine mo_profile_2d_n(zref, z, z0, zt, zq, zR, u_star, b_star, q_star, &
                         del_m, del_t, del_q)

real,    intent(in),  dimension(:)     :: zref
real,    intent(in),  dimension(:,:)   :: z, z0, zt, zq, u_star, b_star, q_star
real,    intent(in) , dimension(:,:)   :: zR ! roughness sublayer length scale
real,    intent(out), dimension(:,:,:) :: del_m, del_t, del_q

integer :: k

do k = 1, size(zref(:))
  call mo_profile_2d (zref(k), zref(k), z, z0, zt, zq, zR, &
       u_star, b_star, q_star, del_m(:,:,k), del_t(:,:,k), del_q(:,:,k))
enddo

end subroutine mo_profile_2d_n

!=======================================================================
subroutine mo_diff_2d_1(z, u_star, b_star, zR, k_m, k_h)

real, intent(in),  dimension(:,:) :: z, u_star, b_star, zR
real, intent(out), dimension(:,:) :: k_m, k_h

real   , dimension(size(z,1),size(z,2),1) :: z_n, k_m_n, k_h_n

z_n(:,:,1) = z

call mo_diff_2d_n(z_n, u_star, b_star, zR, k_m_n, k_h_n)

k_m = k_m_n(:,:,1)
k_h = k_h_n(:,:,1)

end subroutine mo_diff_2d_1

!=======================================================================
subroutine mo_diff_1d_1(z, u_star, b_star, zR, k_m, k_h)

real, intent(in),  dimension(:) :: z, u_star, b_star, zR
real, intent(out), dimension(:) :: k_m, k_h

real, dimension(size(z),1,1) :: z_n, k_m_n, k_h_n
real, dimension(size(z),1)   :: u_star_n, b_star_n, zR_n

z_n   (:,1,1) = z
u_star_n(:,1) = u_star
b_star_n(:,1) = b_star
zR_n(:,1)     = zR

call mo_diff_2d_n(z_n, u_star_n, b_star_n, zR_n, k_m_n, k_h_n)

k_m = k_m_n(:,1,1)
k_h = k_h_n(:,1,1)

end subroutine mo_diff_1d_1

!=======================================================================
subroutine mo_diff_1d_n(z, u_star, b_star, zR, k_m, k_h)

real, intent(in),  dimension(:,:) :: z
real, intent(in),  dimension(:)   :: u_star, b_star, zR
real, intent(out), dimension(:,:) :: k_m, k_h

real, dimension(size(z,1),1)            :: u_star2, b_star2, zR2
real, dimension(size(z,1),1, size(z,2)) :: z2, k_m2, k_h2

integer :: n

do n = 1, size(z,2)
  z2   (:,1,n) = z(:,n)
enddo
u_star2(:,1) = u_star
b_star2(:,1) = b_star
zR2    (:,1) = zR

call mo_diff_2d_n(z2, u_star2, b_star2, zR2, k_m2, k_h2)

do n = 1, size(z,2)
  k_m(:,n) = k_m2(:,1,n)
  k_h(:,n) = k_h2(:,1,n)
enddo

end subroutine mo_diff_1d_n

!=======================================================================
subroutine mo_diff_0d_1(z, u_star, b_star, zR, k_m, k_h)

real, intent(in)  :: z, u_star, b_star, zR
real, intent(out) :: k_m, k_h

integer            :: ni, nj, nk, ier
real, parameter    :: ustar_min = 1.e-10

real, dimension(1,1,1) :: z_, k_m_, k_h_
real, dimension(1,1)   :: u_star_, b_star_, zR_

if(.not.module_is_initialized) call error_mesg('mo_diff_0d_1 in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

ni = 1; nj = 1; nk = 1
z_(1,1,1) = z; u_star_ = u_star; b_star_ = b_star; zR_ = zR
call monin_obukhov_diff(most, vonkarm,                     &
          & ustar_min,                                     &
          & ni, nj, nk, z_, u_star_, b_star_, zR_, k_m_, k_h_, ier)
k_m = k_m_(1,1,1); k_h = k_h_(1,1,1)

end subroutine mo_diff_0d_1

!=======================================================================
subroutine mo_diff_0d_n(z, u_star, b_star, zR, k_m, k_h)

real, intent(in),  dimension(:) :: z
real, intent(in)                :: u_star, b_star, zR
real, intent(out), dimension(:) :: k_m, k_h

integer            :: ni, nj, nk, ier
real, parameter    :: ustar_min = 1.e-10
real, dimension(1,1,size(z)) :: z_, k_m_, k_h_
real, dimension(1,1)         :: u_star_, b_star_, zR_

if(.not.module_is_initialized) call error_mesg('mo_diff_0d_n in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

ni = 1; nj = 1; nk = size(z(:))
z_(1,1,:) = z; u_star_ = u_star; b_star_ = b_star; zR_ = zR
call monin_obukhov_diff(most, vonkarm,                     &
          & ustar_min,                                     &
          & ni, nj, nk, z_, u_star_, b_star_, zR_, k_m_, k_h_, ier)
k_m(:) = k_m_(1,1,:); k_h(:) = k_h_(1,1,:)

end subroutine mo_diff_0d_n

end module monin_obukhov_mod


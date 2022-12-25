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
!> @defgroup monin_obukhov_inter monin_obukhov_inter
!> @ingroup monin_obukhov
!> @brief Utility routines to be used in @ref monin_obukhov_mod

!> @addtogroup monin_obukhov_inter
!> @{
module monin_obukhov_inter

use monin_obukhov_functions_mod, only: most_functions_T

implicit none
private


public :: monin_obukhov_diff
public :: monin_obukhov_drag_1d
public :: monin_obukhov_solve_zeta
public :: monin_obukhov_profile_1d


contains

pure subroutine monin_obukhov_diff(most, vonkarm,     &
     & ustar_min,                                     &
     & ni, nj, nk, z, u_star, b_star, k_m, k_h, ier)
  class(most_functions_T), intent(in) :: most
  real   , intent(in   )                        :: vonkarm
  real   , intent(in   )                        :: ustar_min ! = 1.e-10
  integer, intent(in   )                        :: ni, nj, nk
  real   , intent(in   ), dimension(ni, nj, nk) :: z
  real   , intent(in   ), dimension(ni, nj)     :: u_star, b_star
  real   , intent(  out), dimension(ni, nj, nk) :: k_m, k_h
  integer, intent(  out)                        :: ier

  real , dimension(ni, nj) :: phi_m, phi_h, zeta, uss
  integer :: j, k
  logical, dimension(ni) :: mask

  ier = 0

  mask = .true.
  uss = max(u_star, ustar_min)

  if(most%neutral) then
     do k = 1, size(z,3)
        k_m(:,:,k) = vonkarm *uss*z(:,:,k)
        k_h(:,:,k) = k_m(:,:,k)
     end do
  else
     do k = 1, size(z,3)
        zeta = - vonkarm * b_star*z(:,:,k)/(uss*uss)
        do j = 1, size(z,2)
           call most%derivative_m(ni, mask, zeta(:,j), phi_m(:,j), ier)
           call most%derivative_t(ni, mask, zeta(:,j), phi_h(:,j), ier)
        enddo
        k_m(:,:,k) = vonkarm * uss*z(:,:,k)/phi_m
        k_h(:,:,k) = vonkarm * uss*z(:,:,k)/phi_h
     end do
  endif

end subroutine monin_obukhov_diff


pure subroutine monin_obukhov_drag_1d(most, grav, vonkarm,       &
     & error, zeta_min, max_iter, small,                         &
     & drag_min_heat, drag_min_moist, drag_min_mom,              &
     & n, pt, pt0, z, z0, zt, zq, zR, speed, drag_m, drag_t,     &
     & drag_q, u_star, b_star, rich, zeta, ier, avail)

  class(most_functions_T), intent(in)   :: most ! set of stability functions
  real   , intent(in   )                :: grav
  real   , intent(in   )                :: vonkarm
  real   , intent(in   )                :: error    ! = 1.e-04
  real   , intent(in   )                :: zeta_min ! = 1.e-06
  integer, intent(in   )                :: max_iter ! = 20
  real   , intent(in   )                :: small    ! = 1.e-04
  real   , intent(in   )                :: drag_min_heat, drag_min_moist, drag_min_mom
  integer, intent(in   )                :: n
  real   , intent(in   ), dimension(n)  :: pt, pt0
  real   , intent(in   ), dimension(n)  :: z ! top of the Monin-Obukhov layer (that is, lowest atmos layer height), m
  real   , intent(in   ), dimension(n)  :: z0, zt, zq ! roughness lengths for momentum, heat, and tracers, respectively, m
  real   , intent(in   ), dimension(n)  :: zR ! roughness sublayer length scale, m
  real   , intent(in   ), dimension(n)  :: speed
  real   , intent(inout), dimension(n)  :: drag_m, drag_t, drag_q, u_star, b_star, zeta, rich
  integer, intent(out  )                :: ier
  logical, intent(in   ), dimension(n), optional :: avail  ! provided mask

  real   , dimension(n) :: fm, ft, fq, zz
  logical, dimension(n) :: mask, mask_1, mask_2
  real   , dimension(n) :: delta_b !!, us, bs, qs
  real                  :: r_crit, sqrt_drag_min_heat
  real                  :: sqrt_drag_min_moist, sqrt_drag_min_mom
  real                  :: us, bs, qs
  integer               :: i

  ier = 0
  r_crit = 0.95*most%rich_crit  ! convergence can get slow if one is
                           ! close to rich_crit
  sqrt_drag_min_heat = 0.0
  if(drag_min_heat.ne.0.0) sqrt_drag_min_heat = sqrt(drag_min_heat)
  sqrt_drag_min_moist = 0.0
  if(drag_min_moist.ne.0.0) sqrt_drag_min_moist = sqrt(drag_min_moist)
  sqrt_drag_min_mom = 0.0
  if(drag_min_mom.ne.0.0) sqrt_drag_min_mom = sqrt(drag_min_mom)

  if(present(avail)) then
     mask = avail
  else
     mask = .true.
  endif

  where(mask)
     delta_b = grav*(pt0 - pt)/pt0
     rich    = - z*delta_b/(speed*speed + small)
     zz      = max(z,z0,zt,zq)
  elsewhere
     rich = 0.0
  end where

  if(most%neutral) then

     do i = 1, n
        if(mask(i)) then
           fm(i)   = log(zz(i)/z0(i))
           ft(i)   = log(zz(i)/zt(i))
           fq(i)   = log(zz(i)/zq(i))
           us   = vonkarm/fm(i)
           bs   = vonkarm/ft(i)
           qs   = vonkarm/fq(i)
           drag_m(i)    = us*us
           drag_t(i)    = us*bs
           drag_q(i)    = us*qs
           u_star(i) = us*speed(i)
           b_star(i) = bs*delta_b(i)
        end if
     enddo

  else

     mask_1 = mask .and. rich <  r_crit
     mask_2 = mask .and. rich >= r_crit

     do i = 1, n
        if(mask_2(i)) then
           drag_m(i)   = drag_min_mom
           drag_t(i)   = drag_min_heat
           drag_q(i)   = drag_min_moist
           us       = sqrt_drag_min_mom
           bs       = sqrt_drag_min_heat
           u_star(i)   = us*speed(i)
           b_star(i)   = bs*delta_b(i)
        end if
     enddo

     call monin_obukhov_solve_zeta (most, error, zeta_min, max_iter, small, &
          & n, rich, zz, z0, zt, zq, zR, fm, ft, fq, zeta, mask_1, ier)

     do i = 1, n
        if(mask_1(i)) then
           us   = max(vonkarm/fm(i), sqrt_drag_min_mom)
           bs   = max(vonkarm/ft(i), sqrt_drag_min_heat)
           qs   = max(vonkarm/fq(i), sqrt_drag_min_moist)
           drag_m(i)   = us*us
           drag_t(i)   = us*bs
           drag_q(i)   = us*qs
           u_star(i)   = us*speed(i)
           b_star(i)   = bs*delta_b(i)
        endif
     enddo

  end if

end subroutine monin_obukhov_drag_1d


pure subroutine monin_obukhov_solve_zeta(most, error, zeta_min, max_iter, small,  &
     & n, rich, z, z0, zt, zq, zR, f_m, f_t, f_q, zeta, mask, ier)
  class(most_functions_T), intent(in)     :: most
  real   , intent(in   )                :: error    ! = 1.e-04, solution tolerance
  real   , intent(in   )                :: zeta_min ! = 1.e-06, for zeta < zeta_min solution is assumed neutral
  integer, intent(in   )                :: max_iter ! = 20 maximum number of iteration steps
  real   , intent(in   )                :: small    ! = 1.e-04
  integer, intent(in   )                :: n
  real   , intent(in   ), dimension(n)  :: rich ! bulk Richardson number
  real   , intent(in   ), dimension(n)  :: z ! top of the MO layer (that is, lowest atmos layer height), m
  real   , intent(in   ), dimension(n)  :: z0, zt, zq ! roughness length for momentum, heat, and tracers, respectively m
  real   , intent(in   ), dimension(n)  :: zR ! roughness sublayer length scale, m
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(  out), dimension(n)  :: f_m, f_t, f_q ! final values of integral stability correction functions for momentum, heat, and tracers respectively
  real   , intent(  out), dimension(n)  :: zeta ! solution for zeta (z/L)
  integer, intent(  out)                :: ier

  real    :: max_cor
  integer :: iter
  real, dimension(n) ::   &
       d_rich, rich_1, correction, corr, z_z0, z_zt, z_zq, &
       ln_z_z0, ln_z_zt, ln_z_zq,                          &
       phi_m, phi_m_0, phi_t, phi_t_0, rzeta,              &
       zeta_0, zeta_t, zeta_q, df_m, df_t, l_inv
  logical, dimension(n) :: mask_1, mask_n

!   integer :: i

  ier = 0

  z_z0 = z/z0
  z_zt = z/zt
  z_zq = z/zq
  ln_z_z0 = log(z_z0)
  ln_z_zt = log(z_zt)
  ln_z_zq = log(z_zq)

  corr = 0.0
  mask_1 = mask

  ! initial guess

  zeta = 0.0
  where(mask_1)
     zeta = rich*ln_z_z0*ln_z_z0/ln_z_zt
  end where

  where (mask_1 .and. rich >= 0.0)
     zeta = zeta/(1.0 - rich/most%rich_crit)
  end where

  iter_loop: do iter = 1, max_iter

     ! handle points in neutral or near-neutral condition. Note that with RSL the profile
     ! is only logarithmic where zR == 0
     mask_n = mask_1 .and. (abs(zeta)<zeta_min)

     where (mask_n)
        zeta = 0.0
        f_m = ln_z_z0
        f_t = ln_z_zt
        f_q = ln_z_zq
     end where
     ! add roughness sublayer corrections
     where (mask_n) l_inv = zeta/z
     call most%add_rsl_integral_m(n, mask_n, l_inv, z0, z, zR, f_m, ierr=ier)
     call most%add_rsl_integral_t(n, mask_n, l_inv, zt, z, zR, f_t, ierr=ier)
     call most%add_rsl_integral_q(n, mask_n, l_inv, zq, z, zR, f_q, ierr=ier)
     ! do not do any more calculations at these points
     where (mask_n) mask_1 = .false.

     zeta_0 = 0.0
     zeta_t = 0.0
     zeta_q = 0.0
     where (mask_1)
        rzeta  = 1.0/zeta
        zeta_0 = zeta/z_z0
        zeta_t = zeta/z_zt
        zeta_q = zeta/z_zq
     end where

     call most%derivative_m(n, mask_1, zeta,   phi_m,   ier)
     call most%derivative_m(n, mask_1, zeta_0, phi_m_0, ier)
     call most%derivative_t(n, mask_1, zeta,   phi_t  , ier)
     call most%derivative_t(n, mask_1, zeta_t, phi_t_0, ier)

     where (mask_1)
        df_m  = (phi_m - phi_m_0)*rzeta
        df_t  = (phi_t - phi_t_0)*rzeta
     endwhere

     call most%integral_m(n, mask_1, zeta, zeta_0, ln_z_z0, f_m, ier)
     call most%integral_t(n, mask_1, zeta, zeta_t, ln_z_zt, f_t, ier)
     call most%integral_q(n, mask_1, zeta, zeta_q, ln_z_zq, f_q, ier)

     ! add roughness sublaye corrections
     where (mask_1) l_inv = zeta/z
     call most%add_rsl_integral_m(n, mask_1, l_inv, z0, z, zR, f_m, df_m, ierr=ier)
     call most%add_rsl_integral_t(n, mask_1, l_inv, zt, z, zR, f_t, df_t, ierr=ier)
     ! we need the value of f_q to return to the calling subroutine, but it is not used
     ! in the solver
     call most%add_rsl_integral_q(n, mask_1, l_inv, zq, z, zR, f_q,       ierr=ier)

     where (mask_1)
        rich_1 = zeta*f_t/(f_m*f_m)
        d_rich = rich_1*( rzeta +  df_t/f_t - 2.0 *df_m/f_m)
        correction = (rich - rich_1)/d_rich
        corr = min(abs(correction),abs(correction/zeta))
        ! the criterion corr < error seems to work ok, but is a bit arbitrary
        !  when zeta is small the tolerance is reduced
     end where

     max_cor= maxval(corr)

     if(max_cor > error) then
        mask_1 = mask_1 .and. (corr > error)
        ! change the mask so computation proceeds only on non-converged points
        where(mask_1)
           zeta = zeta + correction
        end where
        cycle iter_loop
     else
        return
     end if

  end do iter_loop

  ier = 1 ! surface drag iteration did not converge

end subroutine monin_obukhov_solve_zeta


pure subroutine monin_obukhov_profile_1d(most, &
     vonkarm, &
     & n, zref, zref_t, z, z0, zt, zq, zR, u_star, b_star, q_star, &
     & del_m, del_t, del_q, ier, avail)

  class(most_functions_T), intent(in) :: most
  real   , intent(in   )                :: vonkarm
  integer, intent(in   )                :: n
  real,    intent(in   )                :: zref, zref_t
  real,    intent(in   ), dimension(n)  :: z, z0, zt, zq, u_star, b_star, q_star
  real,    intent(in   ), dimension(n)  :: zR ! roughness sublayer length scale
  real,    intent(  out), dimension(n)  :: del_m, del_t, del_q
  integer, intent(out  )                :: ier
  logical, intent(in   ), dimension(n), optional :: avail ! provided mask

  real, dimension(n) :: zeta, zeta_0, zeta_t, zeta_q, zeta_ref, zeta_ref_t, &
       ln_z_z0, ln_z_zt, ln_z_zq, ln_z_zref, ln_z_zref_t,  &
       f_m_ref, f_m, f_t_ref, f_t, f_q_ref, f_q,           &
       mo_length_inv, zref_n

  logical, dimension(n) :: mask

  ier = 0

  if (present(avail)) then
     mask = avail
  else
     mask = .true.
  endif

  del_m = 0.0  ! zero output arrays
  del_t = 0.0
  del_q = 0.0

  where(mask)
     ln_z_z0     = log(z/z0)
     ln_z_zt     = log(z/zt)
     ln_z_zq     = log(z/zq)
     ln_z_zref   = log(z/zref)
     ln_z_zref_t = log(z/zref_t)
  endwhere

  if(most%neutral) then

     where(mask)
        del_m = 1.0 - ln_z_zref  /ln_z_z0
        del_t = 1.0 - ln_z_zref_t/ln_z_zt
        del_q = 1.0 - ln_z_zref_t/ln_z_zq
     endwhere

  else

     where(mask .and. u_star > 0.0)
        mo_length_inv = - vonkarm * b_star/(u_star*u_star)
        zeta       = z     *mo_length_inv
        zeta_0     = z0    *mo_length_inv
        zeta_t     = zt    *mo_length_inv
        zeta_q     = zq    *mo_length_inv
        zeta_ref   = zref  *mo_length_inv
        zeta_ref_t = zref_t*mo_length_inv
     endwhere

     call most%integral_m(n, mask, zeta, zeta_0,     ln_z_z0,     f_m,     ier)
     call most%integral_m(n, mask, zeta, zeta_ref,   ln_z_zref,   f_m_ref, ier)

     call most%integral_t(n, mask, zeta, zeta_t,     ln_z_zt,     f_t,     ier)
     call most%integral_t(n, mask, zeta, zeta_ref_t, ln_z_zref_t, f_t_ref, ier)

     call most%integral_q(n, mask, zeta, zeta_q,     ln_z_zq,     f_q,     ier)
     call most%integral_q(n, mask, zeta, zeta_ref_t, ln_z_zref_t, f_q_ref, ier)

     ! add roughness sublayer corrections
     zref_n(:) = zref
     call most%add_rsl_integral_m(n, mask, mo_length_inv, z0,     z, zR, f_m,     ierr=ier)
     call most%add_rsl_integral_m(n, mask, mo_length_inv, zref_n, z, zR, f_m_ref, ierr=ier)

     zref_n(:) = zref_t
     call most%add_rsl_integral_t(n, mask, mo_length_inv, zt,     z, zR, f_t,     ierr=ier)
     call most%add_rsl_integral_t(n, mask, mo_length_inv, zref_n, z, zR, f_t_ref, ierr=ier)

     call most%add_rsl_integral_q(n, mask, mo_length_inv, zq,     z, zR, f_q,     ierr=ier)
     call most%add_rsl_integral_q(n, mask, mo_length_inv, zref_n, z, zR, f_q_ref, ierr=ier)

     where(mask)
        del_m = 1.0 - f_m_ref/f_m
        del_t = 1.0 - f_t_ref/f_t
        del_q = 1.0 - f_q_ref/f_q
     endwhere

  end if
end subroutine monin_obukhov_profile_1d

end module monin_obukhov_inter
!> @}
! close documentation grouping

!-----------------------------------------------------------------------
!  os-axion-primakoff.f03
!
!  Forward:  gamma -> a   (Primakoff conversion in screened Coulomb field)
!  Inverse:  a -> gamma
!
!  Notes:
!   - Uses Poisson macro-events with variance reduction (w_axion, lambda_cap)
!   - Depletes emitter by reducing its macro weight q(i)
!   - Uses deposited ion charge density VDF rho_ion (BH-style) as target density:
!       rho_norm ~ (Z * n_i / n0)
!       n_i [cm^-3] = |rho_norm| * n0 / Z_ion
!   - omega_keV = |p| * 511 keV (massless)
!-----------------------------------------------------------------------
module m_axion_primakoff
  use m_system,         only: mpi_node
  use m_parameters,     only: p_double, p_k_part, p_x_dim
  use m_vdf_define,     only: t_vdf
  use m_photons_define, only: t_photons
  use m_axions_define,  only: t_axions
  use m_axion_primakoff_info, only: t_axion_primakoff_info

  implicit none
  private

  public :: primakoff_photon_to_axion
  public :: primakoff_axion_to_photon

  ! physics helpers (kept private unless you want them exposed)
  private :: primakoff_sigma_cm2, primakoff_rate_sinv, kappa_debye_keV
  private :: poisson_knuth, rho_ion_ngp
  private :: maybe_init_photons_planck, sample_planck_w, sample_planck_n

  real(p_double), parameter :: pi   = 3.1415926535897932384626433832795_p_double
  real(p_double), parameter :: alpha= 7.2973525693e-3_p_double
  real(p_double), parameter :: c_cms= 2.99792458e10_p_double
  real(p_double), parameter :: hbarc_keV_cm = 1.973269804e-8_p_double  ! (ħc) in keV·cm
  real(p_double), parameter :: mec2_keV = 510.99895_p_double


! --- Primakoff balance diagnostics (module state) ---
integer, save :: prim_bal_step = 0
logical, save :: prim_seen_fwd = .false., prim_seen_inv = .false.
real(p_double), save :: prim_qrem_fwd = 0.0_p_double, prim_qrem_inv = 0.0_p_double
real(p_double), save :: prim_qcre_fwd = 0.0_p_double, prim_qcre_inv = 0.0_p_double
real(p_double), save :: prim_qexp_fwd = 0.0_p_double, prim_qexp_inv = 0.0_p_double
real(p_double), save :: prim_Quse_fwd = 0.0_p_double, prim_Quse_inv = 0.0_p_double
real(p_double), save :: prim_omega_q_fwd = -1.0_p_double, prim_omega_q_inv = -1.0_p_double
real(p_double), save :: prim_kappa_q_fwd = -1.0_p_double, prim_kappa_q_inv = -1.0_p_double
real(p_double), save :: prim_omega_qp_fwd = -1.0_p_double, prim_omega_qp_inv = -1.0_p_double
real(p_double), save :: prim_kappa_qp_fwd = -1.0_p_double, prim_kappa_qp_inv = -1.0_p_double

contains

  pure function primakoff_sigma_cm2(omega_keV, kappa_keV, Z, g_agamma_GeVinv) result(sig_cm2)
    real(p_double), intent(in) :: omega_keV, kappa_keV, Z, g_agamma_GeVinv
    real(p_double) :: sig_cm2
    real(p_double) :: s, F, g_keVinv, sig_keVminus2

    sig_cm2 = 0.0_p_double
    if (omega_keV <= 0.0_p_double) return
    if (kappa_keV <= 0.0_p_double) return
    if (Z <= 0.0_p_double) return
    if (g_agamma_GeVinv <= 0.0_p_double) return

    g_keVinv = g_agamma_GeVinv / 1.0e6_p_double  ! 1 GeV = 1e6 keV
    s = 4.0_p_double * omega_keV * omega_keV / (kappa_keV * kappa_keV)

    if (s <= 0.0_p_double) return
    F = (1.0_p_double + 1.0_p_double/s) * log(1.0_p_double + s) - 1.0_p_double
    if (F < 0.0_p_double) F = 0.0_p_double

    sig_keVminus2 = (g_keVinv*g_keVinv) * alpha * (Z*Z) * 0.125_p_double * F  ! 1/8
    sig_cm2 = sig_keVminus2 * (hbarc_keV_cm*hbarc_keV_cm)
  end function primakoff_sigma_cm2

  pure function primakoff_rate_sinv(omega_keV, kappa_keV, n_target_cm3, Z, g_agamma_GeVinv) result(gamma_sinv)
    real(p_double), intent(in) :: omega_keV, kappa_keV, n_target_cm3, Z, g_agamma_GeVinv
    real(p_double) :: gamma_sinv
    real(p_double) :: sig

    gamma_sinv = 0.0_p_double
    if (n_target_cm3 <= 0.0_p_double) return
    sig = primakoff_sigma_cm2(omega_keV, kappa_keV, Z, g_agamma_GeVinv)
    gamma_sinv = n_target_cm3 * sig * c_cms
  end function primakoff_rate_sinv

  pure function kappa_debye_keV(Te_keV, ne_cm3, sum_Z2_ni_cm3) result(kappa_keV)
    real(p_double), intent(in) :: Te_keV, ne_cm3, sum_Z2_ni_cm3
    real(p_double) :: kappa_keV
    real(p_double) :: n_eff_cm3, n_eff_keV3, k2

    kappa_keV = 0.0_p_double
    if (Te_keV <= 0.0_p_double) return

    n_eff_cm3 = max(0.0_p_double, ne_cm3) + max(0.0_p_double, sum_Z2_ni_cm3)
    if (n_eff_cm3 <= 0.0_p_double) return

    n_eff_keV3 = n_eff_cm3 * (hbarc_keV_cm*hbarc_keV_cm*hbarc_keV_cm)
    k2 = 4.0_p_double * pi * alpha * (n_eff_keV3 / Te_keV)
    if (k2 <= 0.0_p_double) return
    kappa_keV = sqrt(k2)
  end function kappa_debye_keV

  integer function poisson_knuth(lambda) result(k)
    real(p_double), intent(in) :: lambda
    real(p_double) :: L, p, u
    integer :: kk

    if (lambda <= 0.0_p_double) then
      k = 0
      return
    end if

    L = exp(-lambda)
    p = 1.0_p_double
    kk = 0
    do
      kk = kk + 1
      call random_number(u)
      if (u <= 0.0_p_double) u = 1.0e-16_p_double
      p = p * u
      if (p <= L) exit
    end do
    k = kk - 1
  end function poisson_knuth

  !---------------------------------------------------------------------------------------
  ! Stochastic rounding depletion for float32 particle weights.
  !
  ! Motivation: photons%q (and often axions%q) are stored in p_k_part (float32) at the OSIRIS
  ! core level. For Primakoff, per-step fractional depletion can be far below the float32 ULP,
  ! so direct assignment q := real(q_old - dq) rounds to q_old (no depletion).
  !
  ! This routine applies an *unbiased* depletion in float32 by stepping to the next representable
  ! float below q with probability dq/ULP when dq < ULP (and iterating if dq spans multiple ULPs).
  !---------------------------------------------------------------------------------------
  subroutine deplete_q_stoch(qf, dq)
    real(p_k_part), intent(inout) :: qf
    real(p_double), intent(in)    :: dq
    real(p_double) :: qd_old, qd_target
    real(p_double) :: q0d, step, prob, u
    real(p_k_part) :: q0, q_down, q_up

    if (dq <= 0.0_p_double) return
    if (qf <= 0.0_p_k_part) return

    ! Work in double for the intended update.
    qd_old    = real(qf, p_double)
    qd_target = qd_old - dq
    if (qd_target <= 0.0_p_double) then
      qf = 0.0_p_k_part
      return
    end if

    ! Stochastically round the *target* value to float32 in an unbiased way.
    ! This avoids any systematic drift that can arise from deterministic rounding.
    q0  = real(qd_target, p_k_part)
    q0d = real(q0, p_double)

    if (q0d == qd_target) then
      qf = q0
      return
    end if

    if (q0d > qd_target) then
      ! q0 rounded above the target: choose between q0 and the next representable below.
      q_down = nearest(q0, -1.0_p_k_part)
      step   = real(q0 - q_down, p_double)
      if (step > 0.0_p_double) then
        prob = (q0d - qd_target) / step
        prob = min(1.0_p_double, max(0.0_p_double, prob))
        call random_number(u)
        if (u < prob) q0 = q_down
      end if
    else
      ! q0 rounded below the target: choose between q0 and the next representable above.
      q_up  = nearest(q0,  1.0_p_k_part)
      step  = real(q_up - q0, p_double)
      if (step > 0.0_p_double) then
        prob = (qd_target - q0d) / step
        prob = min(1.0_p_double, max(0.0_p_double, prob))
        call random_number(u)
        if (u < prob) q0 = q_up
      end if
    end if

    if (q0 < 0.0_p_k_part) q0 = 0.0_p_k_part
    qf = q0
  end subroutine deplete_q_stoch



integer function sample_planck_n(nmax) result(n)
  integer, intent(in) :: nmax
  real(p_double) :: u, s, w
  integer :: k

  n = 1
  if (nmax <= 1) return

  ! draw n with P(n) ∝ 1/n^3 truncated at nmax
  call random_number(u)
  if (u <= 0.0_p_double) u = 1.0e-16_p_double

  s = 0.0_p_double
  do k = 1, nmax
    s = s + 1.0_p_double / real(k, p_double)**3
  end do

  w = 0.0_p_double
  do k = 1, nmax
    w = w + (1.0_p_double / real(k, p_double)**3) / s
    if (u <= w) then
      n = k
      return
    end if
  end do
  n = nmax
end function sample_planck_n


real(p_double) function sample_planck_w(nmax) result(w)
  integer, intent(in) :: nmax
  real(p_double) :: u1, u2, u3
  integer :: n

  ! Draw dimensionless w = ω/T with pdf ∝ w^2/(exp(w)-1)
  ! Series method: choose n with P(n) ∝ 1/n^3, then w = (x1+x2+x3)/n with xi = -ln(Ui)
  n = sample_planck_n(nmax)

  call random_number(u1); if (u1 <= 0.0_p_double) u1 = 1.0e-16_p_double
  call random_number(u2); if (u2 <= 0.0_p_double) u2 = 1.0e-16_p_double
  call random_number(u3); if (u3 <= 0.0_p_double) u3 = 1.0e-16_p_double

  w = (-log(u1) - log(u2) - log(u3)) / real(n, p_double)
end function sample_planck_w


subroutine maybe_init_photons_planck(photons)
  class(t_photons), intent(inout) :: photons
  integer :: i, npar, nmax
  real(p_double) :: Te_keV, w, omega_keV, omega_dimless
  real(p_double) :: u, v, phi, costh, sinth
  real(p_k_part) :: px, py, pz

  associate(info => photons%axion_primakoff_info)

    if (info%photon_init_done) return
    if (.not. info%if_primakoff) return
    if (trim(info%photon_init) /= 'planck') then
      info%photon_init_done = .true.
      return
    end if

    Te_keV = info%Te_keV
    if (info%photon_Te_keV > 0.0_p_double) Te_keV = info%photon_Te_keV
    if (Te_keV <= 0.0_p_double) then
      if (mpi_node() == 0) write(*,*) "[Primakoff] photon_init='planck' but Te_keV<=0; skipping."
      info%photon_init_done = .true.
      return
    end if

    nmax = info%photon_nmax
    if (nmax < 1) nmax = 1

    npar = photons%num_par
    if (npar <= 0) then
      info%photon_init_done = .true.
      return
    end if

    if (mpi_node() == 0) then
      write(*,*) "[Primakoff] Resampling photon momenta to Planck distribution: T_keV=", Te_keV, &
                 " quiet_pairs=", info%photon_quiet_pairs, " nmax=", nmax
    end if

    if (info%photon_quiet_pairs) then
      i = 1
      do while (i <= npar)
        w = sample_planck_w(nmax)
        omega_keV = w * Te_keV
        omega_dimless = omega_keV / mec2_keV

        call random_number(u); if (u <= 0.0_p_double) u = 1.0e-16_p_double
        call random_number(v); if (v <= 0.0_p_double) v = 1.0e-16_p_double

        costh = 1.0_p_double - 2.0_p_double * u
        sinth = sqrt(max(0.0_p_double, 1.0_p_double - costh*costh))
        phi   = 2.0_p_double * pi * v

        px = real(omega_dimless * sinth * cos(phi), p_k_part)
        py = real(omega_dimless * sinth * sin(phi), p_k_part)
        pz = real(omega_dimless * costh, p_k_part)

        photons%p(1,i) = px
        photons%p(2,i) = py
        photons%p(3,i) = pz

        if (i+1 <= npar) then
          photons%p(1,i+1) = -px
          photons%p(2,i+1) = -py
          photons%p(3,i+1) = -pz
        end if

        i = i + 2
      end do
    else
      do i = 1, npar
        w = sample_planck_w(nmax)
        omega_keV = w * Te_keV
        omega_dimless = omega_keV / mec2_keV

        call random_number(u); if (u <= 0.0_p_double) u = 1.0e-16_p_double
        call random_number(v); if (v <= 0.0_p_double) v = 1.0e-16_p_double

        costh = 1.0_p_double - 2.0_p_double * u
        sinth = sqrt(max(0.0_p_double, 1.0_p_double - costh*costh))
        phi   = 2.0_p_double * pi * v

        photons%p(1,i) = real(omega_dimless * sinth * cos(phi), p_k_part)
        photons%p(2,i) = real(omega_dimless * sinth * sin(phi), p_k_part)
        photons%p(3,i) = real(omega_dimless * costh, p_k_part)
      end do
    end if

    info%photon_init_done = .true.
  end associate
end subroutine maybe_init_photons_planck

  pure function rho_ion_ngp(rho_ion, ix) result(rho)
    type(t_vdf), intent(in) :: rho_ion
    integer, intent(in) :: ix(3)
    real(p_double) :: rho

    select case ( p_x_dim )
    case ( 1 )
        if ( ix(1) < lbound(rho_ion%f1, 2) .or. ix(1) > ubound(rho_ion%f1, 2) .or. &
             ix(1) < lbound(rho_ion%f1, 2) .or. ix(1) > ubound(rho_ion%f1, 2) ) then
          rho = 0.0_p_double
        else
          rho = real( rho_ion%f1(1, ix(1)), p_double )
        end if
    case ( 2 )
        if ( ix(1) < lbound(rho_ion%f2, 2) .or. ix(1) > ubound(rho_ion%f2, 2) .or. &
             ix(2) < lbound(rho_ion%f2, 3) .or. ix(2) > ubound(rho_ion%f2, 3) ) then
          rho = 0.0_p_double
        else
          rho = real( rho_ion%f2(1, ix(1), ix(2)), p_double )
        end if
    case ( 3 )
        if ( ix(1) < lbound(rho_ion%f3, 2) .or. ix(1) > ubound(rho_ion%f3, 2) .or. &
             ix(2) < lbound(rho_ion%f3, 3) .or. ix(2) > ubound(rho_ion%f3, 3) .or. &
             ix(3) < lbound(rho_ion%f3, 4) .or. ix(3) > ubound(rho_ion%f3, 4) ) then
          rho = 0.0_p_double
        else
          rho = real( rho_ion%f3(1, ix(1), ix(2), ix(3)), p_double )
        end if
    case default
        rho = 0.0_p_double
    end select
  end function rho_ion_ngp

  subroutine primakoff_photon_to_axion(photons, axions, dt)
    class(t_photons), intent(inout) :: photons
    class(t_axions),  intent(inout) :: axions
    real(p_double),   intent(in)    :: dt

    integer :: i, npar0, N, Nfull
    integer :: ix3(3)
    integer :: newix(p_x_dim)
    real(p_k_part) :: newx(p_x_dim)
    real(p_k_part) :: newp(3)

    real(p_double) :: q_emit_base, q_emit, q_ph, q_total, rem
    real(p_double) :: omega_dimless, omega_keV
    real(p_double) :: rho_norm, n0_cm3, ni_cm3
    real(p_double) :: Te_keV, kappa_keV, ne_cm3, sumZ2ni_cm3
    real(p_double) :: gamma_sinv, a, pconv
    real(p_double) :: dt_phys
    real(p_double) :: lambda, lambda_eff, scale
    real(p_double) :: Z

    integer, save :: dbg_done = 0
    real(p_double) :: lambda_dbg



    ! Debug accounting (weight conservation).
    real(p_double) :: q_removed_step, q_created_step, q_rem_written_step, q_expect_step
    real(p_double) :: q_sum, qomega_sum, qkappa_sum
    real(p_double) :: qp_sum, qpomega_sum, qpkappa_sum
    real(p_double) :: q_old_d, q_new_d, dq_done

    if (dt <= 0.0_p_double) return

    associate(info => photons%axion_primakoff_info)
    q_removed_step = 0.0_p_double
    q_created_step = 0.0_p_double
    q_expect_step = 0.0_p_double
    q_rem_written_step = 0.0_p_double
    q_sum = 0.0_p_double
    qomega_sum = 0.0_p_double
    qkappa_sum = 0.0_p_double
    qp_sum = 0.0_p_double
    qpomega_sum = 0.0_p_double
    qpkappa_sum = 0.0_p_double

    if (.not. info%if_primakoff) return
    if (info%g_agamma <= 0.0_p_double) return
    if (info%w_axion <= 0.0_p_double) return

    call maybe_init_photons_planck(photons)

    n0_cm3 = real(photons%n0, p_double)
    if (n0_cm3 <= 0.0_p_double) return

    if (real(photons%omega_p0, p_double) <= 0.0_p_double) return
    dt_phys = dt / real(photons%omega_p0, p_double)

    Z = real(info%Z_ion, p_double)
    if (Z < 1.0_p_double) Z = 1.0_p_double
    Te_keV = info%Te_keV

    npar0 = photons%num_par
    do i = 1, npar0

      q_ph = real(photons%q(i), p_double)
      q_old_d = real(photons%q(i), p_double)
      if (q_ph <= 0.0_p_double) cycle

      omega_dimless = sqrt( real(photons%p(1,i),p_double)**2 + real(photons%p(2,i),p_double)**2 + real(photons%p(3,i),p_double)**2 )
      if (omega_dimless <= 0.0_p_double) cycle
      omega_keV = omega_dimless * mec2_keV

      ix3 = 1
      ix3(1) = photons%ix(1,i)
      if (p_x_dim >= 2) ix3(2) = photons%ix(2,i)
      if (p_x_dim >= 3) ix3(3) = photons%ix(3,i)

      rho_norm = abs( rho_ion_ngp(photons%axion_primakoff_info%rho_ion, ix3) )

      ni_cm3 = rho_norm * n0_cm3 / Z
      if (ni_cm3 <= 0.0_p_double) cycle

      if (info%use_screening) then
        if (info%kappa_keV > 0.0_p_double) then
          kappa_keV = info%kappa_keV
        else
          sumZ2ni_cm3 = (Z*Z) * ni_cm3
          if (info%include_ee) then
            if (info%assume_neutral) then
              ne_cm3 = info%ne_factor * Z * ni_cm3
            else
              ne_cm3 = 0.0_p_double
            end if
          else
            ne_cm3 = 0.0_p_double
          end if
          kappa_keV = kappa_debye_keV(Te_keV, ne_cm3, sumZ2ni_cm3)
        end if
      else
        kappa_keV = 1.0e-9_p_double
      end if
      if (kappa_keV <= 0.0_p_double) cycle

      gamma_sinv = primakoff_rate_sinv(omega_keV, kappa_keV, ni_cm3, Z, info%g_agamma)
      if (gamma_sinv <= 0.0_p_double) cycle

      a = gamma_sinv * dt_phys
      if (a <= 0.0_p_double) cycle
      if (a < 1.0e-8_p_double) then
        pconv = a
      else
        pconv = 1.0_p_double - exp(-a)
      end if
      if (pconv <= 0.0_p_double) cycle
      q_expect_step = q_expect_step + (q_ph * pconv)
      q_sum = q_sum + q_ph
      qomega_sum = qomega_sum + q_ph * omega_keV
      qkappa_sum = qkappa_sum + q_ph * kappa_keV
      qp_sum = qp_sum + q_ph * pconv
      qpomega_sum = qpomega_sum + q_ph * pconv * omega_keV
      qpkappa_sum = qpkappa_sum + q_ph * pconv * kappa_keV

      if (dbg_done == 0) then
        lambda_dbg = (q_ph * pconv) / info%w_axion
        if (mpi_node() == 0) then
          write(*,*) "[AXION PRIM ONEPH] E_keV=", omega_keV, " kappa_keV=", kappa_keV
          write(*,*) "[AXION PRIM ONEPH] rate_sinv=", gamma_sinv, " dt_s=", dt_phys, " a=", a, " pconv=", pconv
          write(*,*) "[AXION PRIM ONEPH] qph=", q_ph, " w_axion=", info%w_axion, " lambda=", lambda_dbg
        end if
        dbg_done = 1
      end if


      lambda = (q_ph * pconv) / info%w_axion
      if (lambda <= 0.0_p_double) cycle

      lambda_eff = min(lambda, info%lambda_cap)
      if (lambda_eff <= 0.0_p_double) cycle

      scale = lambda / lambda_eff
      q_emit_base = info%w_axion
      q_emit = q_emit_base * scale

      N = poisson_knuth(lambda_eff)
      if (N <= 0) cycle

      q_total = q_emit * real(N, p_double)

      if (q_total >= q_ph) then
        Nfull = int(q_ph / q_emit)
        rem   = q_ph - real(Nfull, p_double) * q_emit
        q_total = q_ph
      else
        Nfull = N
        rem   = 0.0_p_double
      end if

      q_removed_step = q_removed_step + q_total
      q_created_step = q_created_step + ( real(Nfull, p_double)*q_emit + rem )

      if (Nfull <= 0 .and. rem <= 0.0_p_double) cycle

      newix = photons%ix(1:p_x_dim,i)
      newx  = photons%x(1:p_x_dim,i)
      newp  = photons%p(:,i)

      if (Nfull > 0) then
        do N = 1, Nfull
          call axions%create_particle( newix, newx, newp, real(q_emit, p_k_part) )
        end do
      end if
      if (rem > 0.0_p_double) then
        call axions%create_particle( newix, newx, newp, real(rem, p_k_part) )
      end if

      q_ph = q_ph - q_total
      if (q_ph < 0.0_p_double) q_ph = 0.0_p_double
      call deplete_q_stoch(photons%q(i), q_total)
      q_new_d = real(photons%q(i), p_double)
      dq_done = q_old_d - q_new_d
      q_rem_written_step = q_rem_written_step + (q_total - dq_done)



    end do
    if (info%debug) then
      if (mpi_node() == 0) then
        write(*,*) "[AXION PRIM STEP] photon->axion: q_removed=", q_removed_step, &
                   " q_created=", q_created_step, " q_round_loss=", q_rem_written_step
      end if
    end if



    ! --- Primakoff balance diagnostics: record one line per full (fwd+inv) step when info%debug ---
    if (info%debug) then
      if (mpi_node() == 0) then
        prim_Quse_fwd = q_sum
        if (q_sum > 0.0_p_double) then
          prim_omega_q_fwd = qomega_sum / q_sum
          prim_kappa_q_fwd = qkappa_sum / q_sum
        else
          prim_omega_q_fwd = -1.0_p_double
          prim_kappa_q_fwd = -1.0_p_double
        end if
        if (qp_sum > 0.0_p_double) then
          prim_omega_qp_fwd = qpomega_sum / qp_sum
          prim_kappa_qp_fwd = qpkappa_sum / qp_sum
        else
          prim_omega_qp_fwd = -1.0_p_double
          prim_kappa_qp_fwd = -1.0_p_double
        end if
        prim_qrem_fwd = q_removed_step
        prim_qcre_fwd = q_created_step
        prim_qexp_fwd = q_expect_step
        prim_seen_fwd = .true.
        if (prim_seen_fwd .and. prim_seen_inv) then
          prim_bal_step = prim_bal_step + 1
          ! totals after both conversions (whichever call arrives second writes)
          ! NOTE: these are FULL species totals (not RAW sampling).
          call write_prim_balance_line( prim_bal_step, dt, dt_phys, photons, axions )
          prim_seen_fwd = .false.; prim_seen_inv = .false.
        end if
      end if
    end if

    end associate

end subroutine primakoff_photon_to_axion

  subroutine primakoff_axion_to_photon(axions, photons, dt)
    class(t_axions),  intent(inout) :: axions
    class(t_photons), intent(inout) :: photons
    real(p_double),   intent(in)    :: dt

    integer :: i, npar0, N, Nfull
    integer :: ix3(3)
    integer :: newix(p_x_dim)
    real(p_k_part) :: newx(p_x_dim)
    real(p_k_part) :: newp(3)

    real(p_double) :: q_emit_base, q_emit, q_ax, q_total, rem
    real(p_double) :: omega_dimless, omega_keV
    real(p_double) :: rho_norm, n0_cm3, ni_cm3
    real(p_double) :: Te_keV, kappa_keV, ne_cm3, sumZ2ni_cm3
    real(p_double) :: gamma_sinv, a, pconv
    real(p_double) :: dt_phys
    real(p_double) :: lambda, lambda_eff, scale
    real(p_double) :: Z


    ! Debug accounting (weight conservation).
    real(p_double) :: q_removed_step, q_created_step, q_rem_written_step, q_expect_step
    real(p_double) :: q_sum, qomega_sum, qkappa_sum
    real(p_double) :: qp_sum, qpomega_sum, qpkappa_sum
    real(p_double) :: q_old_d, q_new_d, dq_done

    if (dt <= 0.0_p_double) return


    if (.not. photons%axion_primakoff_info%if_inverse) return

    associate(info => photons%axion_primakoff_info)
    q_removed_step = 0.0_p_double
    q_created_step = 0.0_p_double
    q_expect_step = 0.0_p_double
    q_rem_written_step = 0.0_p_double
    q_sum = 0.0_p_double
    qomega_sum = 0.0_p_double
    qkappa_sum = 0.0_p_double
    qp_sum = 0.0_p_double
    qpomega_sum = 0.0_p_double
    qpkappa_sum = 0.0_p_double

    if (info%g_agamma <= 0.0_p_double) return
    if (info%w_axion <= 0.0_p_double) return

    call maybe_init_photons_planck(photons)

    n0_cm3 = real(photons%n0, p_double)
    if (n0_cm3 <= 0.0_p_double) return

    if (real(photons%omega_p0, p_double) <= 0.0_p_double) return
    dt_phys = dt / real(photons%omega_p0, p_double)

    Z = real(info%Z_ion, p_double)
    if (Z < 1.0_p_double) Z = 1.0_p_double
    Te_keV = info%Te_keV

    npar0 = axions%num_par
    do i = 1, npar0

      q_ax = real(axions%q(i), p_double)
      q_old_d = real(axions%q(i), p_double)
      if (q_ax <= 0.0_p_double) cycle

      omega_dimless = sqrt( real(axions%p(1,i),p_double)**2 + real(axions%p(2,i),p_double)**2 + real(axions%p(3,i),p_double)**2 )
      if (omega_dimless <= 0.0_p_double) cycle
      omega_keV = omega_dimless * mec2_keV

      ix3 = 1
      ix3(1) = axions%ix(1,i)
      if (p_x_dim >= 2) ix3(2) = axions%ix(2,i)
      if (p_x_dim >= 3) ix3(3) = axions%ix(3,i)

      rho_norm = abs( rho_ion_ngp(photons%axion_primakoff_info%rho_ion, ix3) )

      ni_cm3 = rho_norm * n0_cm3 / Z
      if (ni_cm3 <= 0.0_p_double) cycle

      if (info%use_screening) then
        if (info%kappa_keV > 0.0_p_double) then
          kappa_keV = info%kappa_keV
        else
          sumZ2ni_cm3 = (Z*Z) * ni_cm3
          if (info%include_ee) then
            if (info%assume_neutral) then
              ne_cm3 = info%ne_factor * Z * ni_cm3
            else
              ne_cm3 = 0.0_p_double
            end if
          else
            ne_cm3 = 0.0_p_double
          end if
          kappa_keV = kappa_debye_keV(Te_keV, ne_cm3, sumZ2ni_cm3)
        end if
      else
        kappa_keV = 1.0e-9_p_double
      end if
      if (kappa_keV <= 0.0_p_double) cycle

      gamma_sinv = primakoff_rate_sinv(omega_keV, kappa_keV, ni_cm3, Z, info%g_agamma)
      if (gamma_sinv <= 0.0_p_double) cycle

      a = gamma_sinv * dt_phys
      if (a <= 0.0_p_double) cycle
      if (a < 1.0e-8_p_double) then
        pconv = a
      else
        pconv = 1.0_p_double - exp(-a)
      end if
      if (pconv <= 0.0_p_double) cycle
      q_expect_step = q_expect_step + (q_ax * pconv)
      q_sum = q_sum + q_ax
      qomega_sum = qomega_sum + q_ax * omega_keV
      qkappa_sum = qkappa_sum + q_ax * kappa_keV
      qp_sum = qp_sum + q_ax * pconv
      qpomega_sum = qpomega_sum + q_ax * pconv * omega_keV
      qpkappa_sum = qpkappa_sum + q_ax * pconv * kappa_keV

      lambda = (q_ax * pconv) / info%w_axion
      if (lambda <= 0.0_p_double) cycle

      lambda_eff = min(lambda, info%lambda_cap)
      if (lambda_eff <= 0.0_p_double) cycle

      scale = lambda / lambda_eff
      q_emit_base = info%w_axion
      q_emit = q_emit_base * scale

      N = poisson_knuth(lambda_eff)
      if (N <= 0) cycle

      q_total = q_emit * real(N, p_double)

      if (q_total >= q_ax) then
        Nfull = int(q_ax / q_emit)
        rem   = q_ax - real(Nfull, p_double) * q_emit
        q_total = q_ax
      else
        Nfull = N
        rem   = 0.0_p_double
      end if

      q_removed_step = q_removed_step + q_total
      q_created_step = q_created_step + ( real(Nfull, p_double)*q_emit + rem )

      if (Nfull <= 0 .and. rem <= 0.0_p_double) cycle

      newix = axions%ix(1:p_x_dim,i)
      newx  = axions%x(1:p_x_dim,i)
      newp  = axions%p(:,i)

      if (Nfull > 0) then
        do N = 1, Nfull
          call photons%create_particle( newix, newx, newp, real(q_emit, p_k_part) )
        end do
      end if
      if (rem > 0.0_p_double) then
        call photons%create_particle( newix, newx, newp, real(rem, p_k_part) )
      end if

      q_ax = q_ax - q_total
      if (q_ax < 0.0_p_double) q_ax = 0.0_p_double
      call deplete_q_stoch(axions%q(i), q_total)
      q_new_d = real(axions%q(i), p_double)
      dq_done = q_old_d - q_new_d
      q_rem_written_step = q_rem_written_step + (q_total - dq_done)



    end do
    if (info%debug) then
      if (mpi_node() == 0) then
        write(*,*) "[AXION PRIM STEP] axion->photon: q_removed=", q_removed_step, &
                   " q_created=", q_created_step, " q_round_loss=", q_rem_written_step
      end if
    end if



    ! --- Primakoff balance diagnostics: record one line per full (fwd+inv) step when info%debug ---
    if (info%debug) then
      if (mpi_node() == 0) then
        prim_Quse_inv = q_sum
        if (q_sum > 0.0_p_double) then
          prim_omega_q_inv = qomega_sum / q_sum
          prim_kappa_q_inv = qkappa_sum / q_sum
        else
          prim_omega_q_inv = -1.0_p_double
          prim_kappa_q_inv = -1.0_p_double
        end if
        if (qp_sum > 0.0_p_double) then
          prim_omega_qp_inv = qpomega_sum / qp_sum
          prim_kappa_qp_inv = qpkappa_sum / qp_sum
        else
          prim_omega_qp_inv = -1.0_p_double
          prim_kappa_qp_inv = -1.0_p_double
        end if
        prim_qrem_inv = q_removed_step
        prim_qcre_inv = q_created_step
        prim_qexp_inv = q_expect_step
        prim_seen_inv = .true.
        if (prim_seen_fwd .and. prim_seen_inv) then
          prim_bal_step = prim_bal_step + 1
          ! totals after both conversions (whichever call arrives second writes)
          ! NOTE: these are FULL species totals (not RAW sampling).
          call write_prim_balance_line( prim_bal_step, dt, dt_phys, photons, axions )
          prim_seen_fwd = .false.; prim_seen_inv = .false.
        end if
      end if
    end if

    end associate

end subroutine primakoff_axion_to_photon


  subroutine write_prim_balance_line(step, dt_norm, dt_phys, photons, axions)
    integer, intent(in) :: step
    real(p_double), intent(in) :: dt_norm, dt_phys
    class(t_photons), intent(in) :: photons
    class(t_axions),  intent(in) :: axions
    integer :: iu, i
    real(p_double) :: Qph, Qax
    Qph = 0.0_p_double
    do i = 1, photons%num_par
      Qph = Qph + real(photons%q(i), p_double)
    end do
    Qax = 0.0_p_double
    do i = 1, axions%num_par
      Qax = Qax + real(axions%q(i), p_double)
    end do
    open(newunit=iu, file='primakoff_balance.dat', status='unknown', position='append', action='write')
    write(iu,'(I8,1X,ES12.4,1X,ES12.4,1X,ES16.8,1X,ES16.8,1X,ES16.8,1X,ES16.8,1X,ES16.8,1X,ES16.8)') &
      step, dt_norm, dt_phys, Qph, Qax, prim_qexp_fwd, prim_qexp_inv, prim_qrem_fwd, prim_qrem_inv
    close(iu)
    open(newunit=iu, file='primakoff_balance_extra.dat', status='unknown', position='append', action='write')
    if (step == 1) then
      write(iu,'(A)') '# step dt_norm dt_phys Qph Qax qexp_fwd qexp_inv qrem_fwd qrem_inv Quse_fwd Quse_inv omega_q_fwd omega_q_inv kappa_q_fwd kappa_q_inv omega_qp_fwd omega_qp_inv kappa_qp_fwd kappa_qp_inv'
    end if
    write(iu,*) step, dt_norm, dt_phys, Qph, Qax, prim_qexp_fwd, prim_qexp_inv, prim_qrem_fwd, prim_qrem_inv, &
    & prim_Quse_fwd, prim_Quse_inv, prim_omega_q_fwd, prim_omega_q_inv, prim_kappa_q_fwd, prim_kappa_q_inv, &
    & prim_omega_qp_fwd, prim_omega_qp_inv, prim_kappa_qp_fwd, prim_kappa_qp_inv
    close(iu)
  end subroutine write_prim_balance_line

end module m_axion_primakoff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     QED particle species class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_push_qed

#include "memory/memory.h"

   use m_system
   use m_parameters
   use m_math
   use m_species_define
   use m_species_define_qed
   use m_species_current
   use m_vdf_define
   use m_emf_define, only: t_emf
   use m_qed
   use m_qed_coulomb

   implicit none

   private

   interface dudt_boris_qed
      module procedure dudt_boris_qed
   end interface

   interface dudt_exact_qed
      module procedure dudt_exact_qed
   end interface

   public :: dudt_boris_qed, dudt_exact_qed

   integer, parameter :: p_iter_max = 10

contains

!-----------------------------------------------------------------------------------------
   subroutine dudt_boris_qed( this, emf, dt, i0, i1, energy, t )
!-----------------------------------------------------------------------------------------
! Advance velocities and calculate time-centered total energy
! process particles in the range [i0,i1]
! This routine is meant to be used for testing simd code and is written to use the same numerics
!-----------------------------------------------------------------------------------------

      implicit none

      class( t_species ), intent(inout) :: this
      class( t_emf ), intent( in )  ::  emf
      real(p_double),                 intent(in) :: dt
      integer, intent(in) :: i0, i1
      real(p_double), intent(inout) :: energy
      real(p_double), intent(in) :: t

      ! dummy variables
      real(p_k_part) :: tem
      integer :: i, ptrcur, np, pp
      real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep, utemp
      real(p_k_part), dimension(p_cache_size) ::  gam_tem, otsq
      real(p_k_part) :: gamma, u2
      real(p_double) :: loc_ene

      ! qed specific variables
      real(p_double), dimension(p_cache_size) :: gl, k_damp
      real(p_k_part), dimension(p_p_dim,p_cache_size) :: bpm, epm, um

      ! executable statements

      select type( this )
       class is ( t_species_qed )

         ! get factor that includes timestep and charge-to-mass ratio
         ! note that the charge-to-mass ratio is used since the
         ! momentum p is not the total momentum of a particles but
         ! the momentum per unit restmass (which is the electron mass)

         tem = real( 0.5_p_double * (dt / this%rqm), p_k_part )

         !loop through all particles
         do ptrcur = i0, i1, p_cache_size

            ! check if last copy of table and set np
            if( ptrcur + p_cache_size > i1 ) then
               np = i1 - ptrcur + 1
            else
               np = p_cache_size
            endif

            ! modify bp & ep to include timestep and charge-to-mass ratio
            ! and perform half the electric field acceleration.
            ! Result is stored in UTEMP.

            call this % get_emf( emf, bp, ep, np, ptrcur, t=t )

            call emit_particles( this, ep, bp, dt, gl, k_damp, ptrcur, np )

            if (this%if_damp_classical) then
               ! copy the momentum at t = n-1/2
               pp = ptrcur
               do i=1,np
                  um(:,i) = this%p(:,pp)
                  pp = pp+1
               enddo

               ! copy the electric and magnetic fields
               do i=1, np
                  epm(:,i) = ep(:,i)
                  bpm(:,i) = bp(:,i)
               end do
            endif

            do i=1, np
               ep(1,i) = ep(1,i) * tem
               ep(2,i) = ep(2,i) * tem
               ep(3,i) = ep(3,i) * tem
            end do

            loc_ene = 0

            pp = ptrcur
            do i=1,np
               utemp(1,i) = this%p(1,pp) + ep(1,i)
               utemp(2,i) = this%p(2,pp) + ep(2,i)
               utemp(3,i) = this%p(3,pp) + ep(3,i)

               ! Get time centered gamma
               u2 =  (utemp(1,i)**2 + utemp(2,i)**2) + utemp(3,i)**2
               gamma = sqrt( u2 + 1 )

               ! accumulate time centered energy
               ! this is done in double precicion always
               loc_ene = loc_ene + this%q(pp) * u2 / (gamma + 1.0_p_double)

               gam_tem(i)= tem / gamma

               pp = pp + 1
            enddo

            ! accumulate global energy
            energy = energy + loc_ene

            do i=1,np
               bp(1,i) = bp(1,i)*gam_tem(i)
               bp(2,i) = bp(2,i)*gam_tem(i)
               bp(3,i) = bp(3,i)*gam_tem(i)
            end do

            pp = ptrcur
            do i=1,np
               this%p(1,pp) = utemp(1,i) + utemp(2,i) * bp(3,i)
               this%p(2,pp) = utemp(2,i) + utemp(3,i) * bp(1,i)
               this%p(3,pp) = utemp(3,i) + utemp(1,i) * bp(2,i)
               pp = pp + 1
            end do

            pp = ptrcur
            do i=1,np
               this%p(1,pp) = this%p(1,pp) - utemp(3,i) * bp(2,i)
               this%p(2,pp) = this%p(2,pp) - utemp(1,i) * bp(3,i)
               this%p(3,pp) = this%p(3,pp) - utemp(2,i) * bp(1,i)
               pp = pp + 1
            end do

            do i=1,np
               otsq(i) = 2.0_p_k_part / ( ((1.0_p_k_part + bp(1,i)**2) + bp(2,i)**2) + bp(3,i)**2)
            end do

            do i=1,np
               bp(1,i) = bp(1,i) * otsq(i)
               bp(2,i) = bp(2,i) * otsq(i)
               bp(3,i) = bp(3,i) * otsq(i)
            end do

            pp = ptrcur
            do i=1,np
               utemp(1,i) = utemp(1,i) + this%p(2,pp) * bp(3,i)
               utemp(2,i) = utemp(2,i) + this%p(3,pp) * bp(1,i)
               utemp(3,i) = utemp(3,i) + this%p(1,pp) * bp(2,i)
               pp = pp + 1
            end do

            pp = ptrcur
            do i=1,np
               utemp(1,i) = utemp(1,i) - this%p(3,pp) * bp(2,i)
               utemp(2,i) = utemp(2,i) - this%p(1,pp) * bp(3,i)
               utemp(3,i) = utemp(3,i) - this%p(2,pp) * bp(1,i)
               pp = pp + 1
            end do

            ! Perform second half of electric field acceleration.
            pp = ptrcur
            do i=1,np
               this%p(1,pp) = utemp(1,i) + ep(1,i)
               this%p(2,pp) = utemp(2,i) + ep(2,i)
               this%p(3,pp) = utemp(3,i) + ep(3,i)
               pp = pp + 1
            end do

            if (this%if_damp_qed) call damp_qed( this, k_damp, ptrcur, np )
            if (this%if_damp_classical) call damp_classical( this, epm, bpm, dt, gl, um, k_damp, ptrcur, np )

            ! bremsstrahlung photon emission
            if (this % bremsstrahlung_info % if_bremsstrahlung) then
               call bremsstrahlung_emission(this, dt, k_damp, ptrcur, i0, i1, np)
            endif

            ! bremsstrahlung axion emission
            if (this % axion_bremsstrahlung_info % if_bremsstrahlung) then
               call axion_bremsstrahlung_emission(this, dt, k_damp, ptrcur, i0, i1, np)
            endif

            ! Compton-like axion emission
            if (this % axion_compton_info % if_compton) then
               call axion_compton_emission(this, dt, k_damp, ptrcur, i0, i1, np)
            endif

         enddo

       class default
         ! this must never happen
         call abort_program( p_err_invalid )
      end select

   end subroutine dudt_boris_qed
!-----------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Emits particles according to selected pair production model
!---------------------------------------------------------------------------------------------------
   subroutine emit_particles( this, ep, bp, dt, gl, k_damp, ptrcur, np )
      !---------------------------------------------------------------------------------------------------

      use m_random
      use m_species_diagnostics_qed

      implicit none

      class( t_species_qed ),    intent(inout) :: this
      real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: ep, bp
      real(p_double), intent(in) :: dt
      real(p_double), dimension(p_cache_size), intent(inout) :: gl, k_damp
      integer, intent(in) :: ptrcur, np

      ! process pair production model
      select case ( this%pairprod )
       case ( p_pairprod_qed )
         call emit_particles_qed( this, ep, bp, dt, gl, k_damp, ptrcur, np )
       case ( p_pairprod_gthr )
         call emit_particles_gthr( this, ep, bp, dt, gl, k_damp, ptrcur, np )
      end select

   end subroutine emit_particles
   !---------------------------------------------------------------------------------------------------

   !---------------------------------------------------------------------------------------------------
   ! Emits particles (photons) according to QED pair production model
   !---------------------------------------------------------------------------
   subroutine emit_particles_qed( this, ep, bp, dt, gl, k_damp, ptrcur, np )
!---------------------------------------------------------------------------------------------------

      use m_random
      use m_species_diagnostics_qed

      implicit none

      class( t_species_qed ),    intent(inout) :: this
      real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: ep, bp
      real(p_double), intent(in) :: dt
      real(p_double), dimension(p_cache_size), intent(inout) :: gl, k_damp
      integer, intent(in) :: ptrcur, np

      ! dummy variables
      integer :: i, pp
      real(p_double) :: eta, eta_s, chi_g, coef_QED, norm_schw, p2, pdotE2, &
         gamE_plus_pcrossB2, varand, Wrad, interp

      real(p_double), dimension(2) :: energy_rad

      real(p_k_part) :: p_frac, p_g, q_g
      real(p_k_part), dimension(p_p_dim) :: p_g_vec
      integer, dimension(p_x_dim) :: ix_g

      real( p_double ), parameter :: M_PI  = 3.14159265358979323846264338327950288d0 ! pi

      ! Values from the NIST reference on Constants, Units and Uncertainty
      ! http://physics.nist.gov/cuu/Constants/index.html
      real( p_double ), parameter :: hbar  = 1.054571800d-34 ! Planck constant over 2 pi
      real( p_double ), parameter :: cl    = 299792458.0d0   ! speed of light in vacuum [m/s]
      real( p_double ), parameter :: me    = 9.10938356d-31  ! electron rest mass [kg]
      real( p_double ), parameter :: alpha = 7.2973525664d-3 ! fine-structure constant

      ! executable statements

      ! calculate the quantum coefficient
      norm_schw = hbar*this%omega_p0/(me*cl**2.d0)
      coef_QED = sqrt(3.d0)*alpha/(2.d0*M_PI*norm_schw)

      ! reset energy radiation diagnostic
      energy_rad = 0.0_p_double

      ! reset damping factor
      k_damp = 1.0_p_double

      pp = ptrcur
      do i = 1, np
         ! compute lorentz factor
         p2 = this%p(1,pp)**2+this%p(2,pp)**2+this%p(3,pp)**2
         gl(i) = sqrt(1.0d0+p2)

         ! compute eta parameter terms
         pdotE2 = (this%p(1,pp)*ep(1,i)+this%p(2,pp)*ep(2,i)+this%p(3,pp)*ep(3,i))**2
         gamE_plus_pcrossB2 = (-gl(i)*ep(1,i)-bp(3,i)*this%p(2,pp)+bp(2,i)*this%p(3,pp))**2 &
            + (-gl(i)*ep(2,i)+bp(3,i)*this%p(1,pp)-bp(1,i)*this%p(3,pp))**2 &
            + (-gl(i)*ep(3,i)-bp(2,i)*this%p(1,pp)+bp(1,i)*this%p(2,pp))**2


         ! compute the quantum parameter eta
         eta = sqrt(abs(gamE_plus_pcrossB2-pdotE2))*norm_schw
         eta_s = eta * this%qed_zeta

         if ( (eta_s > 1.0d-4) .and. (eta_s < 5.0d4) .and. (gl(i) > this%qed_g_cutoff) ) then

            if (eta_s > 300.d0) then
               Wrad = (1.46d0*alpha/norm_schw)*(eta_s**(2.0d0/3.0d0))/gl(i)
            else
               interp = interp_spec_QED(eta_s)
               Wrad = interp*coef_QED*eta_s/gl(i)
            endif

            call random_number( varand )

            if ( varand < Wrad * dt ) then
               ! compute photon chi
               chi_g = find_photon_chi(eta_s)

               ! compute momentum fraction lost to emitted photon
               p_frac = 2.0d0*chi_g/eta_s
               k_damp(i) = real(1.0 - p_frac, p_k_part)

               ! compute photon momentum
               p_g = p_frac * gl(i)

               ! compute photon momentum vector
               p_g_vec = p_frac*this%p(:,pp)

               ! compute photon weight
               q_g = abs(this%q(pp))
               ix_g = this%ix(1:p_x_dim,pp)

               ! Diagnostics
               if (this%ndump_fac_radspect > 0) call radspect_hist( this, p_g, q_g )
               if (this%ndump_fac_radanglespect > 0) call radanglespect_hist( this, p_g_vec, q_g )
               if (this%ndump_fac_radsphere > 0) call radsphere_hist( this, p_g_vec, q_g )
               if (this%ndump_fac_raddetector>0) call detector_hist( this, p_g_vec, q_g )
               if (this%ndump_fac_chi > 0) call chi_hist( this, ix_g, q_g, real(eta_s, p_k_part) )
               if (this%ndump_fac_chi_emit > 0) call emission_chi_hist( this, real(eta_s, p_k_part), real(chi_g, p_k_part), real(q_g, p_k_part) )

               ! condition for photons that can emit a pair
               if (p_g > this%p_emit_cutoff) then

                  call this%photons%create_particle( &
                     this%ix(1:p_x_dim,pp), &
                     this%x(1:p_x_dim,pp), &
                     p_g_vec, q_g &
                     )

                  if (this%ndump_fac_n_emit > 0) call n_emit_hist( this, ix_g )

                  ! add photon energy to QED tracked bin
                  energy_rad(2) = energy_rad(2) +  p_g * q_g
               else
                  ! add photon energy to QED non-tracked bin
                  energy_rad(1) = energy_rad(1) +  p_g * q_g
               endif

            endif
         endif

         pp = pp + 1
      end do

      ! accumulate radiated energy
      this % energy_rad(2) = this % energy_rad(2) + energy_rad(1)
      this % energy_rad(3) = this % energy_rad(3) + energy_rad(2)

   end subroutine emit_particles_qed
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Emits particles (electrons + positrons) according to in-place,
! energy threshold pair production model
!---------------------------------------------------------------------------------------------------
   subroutine emit_particles_gthr( this, ep, bp, dt, gl, k_damp, ptrcur, np )
      !---------------------------------------------------------------------------------------------------

      use m_random
      use m_species_diagnostics_qed

      implicit none

      class( t_species_qed ),    intent(inout) :: this
      real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: ep, bp
      real(p_double), intent(in) :: dt
      real(p_double), dimension(p_cache_size), intent(inout) :: gl, k_damp
      integer, intent(in) :: ptrcur, np

      ! dummy variables
      integer :: i, pp, n_emit
      real(p_double) :: p2, p_frac
      real(p_k_part), dimension(p_p_dim) :: p_pair_vec

      ! executable statements
      ! do nothing to the energy radiation diagnostic
      ! reset damping factor
      k_damp = 1.0_p_double

      pp = ptrcur
      do i = 1, np
         ! compute lorentz factor
         p2 = this%p(1,pp)**2+this%p(2,pp)**2+this%p(3,pp)**2
         gl(i) = sqrt(1.0d0+p2)

         if (gl(i) > this%pairprod_gthr) then
            ! compute momentum fraction lost to emitted photon
            p_frac = this%pairprod_gpair / gl(i)

            p_pair_vec = 0.5 * p_frac * this%p(:,pp)

            call this%create_particle( &
               this%ix(1:p_x_dim,pp), &
               this%x(1:p_x_dim,pp), &
               p_pair_vec, &
               this%q(pp) &
               )

            call this%recip_spec%create_particle( &
               this%ix(1:p_x_dim,pp), &
               this%x(1:p_x_dim,pp), &
               p_pair_vec, &
               -this%q(pp) &
               )

            ! register emission
            if (this%ndump_fac_n_emit > 0) then
               select case ( p_x_dim )
                case ( 1 )
                  n_emit = this%n_emit%f1(1, this%ix(1,pp))
                  this%n_emit%f1(1, this%ix(1,pp)) = n_emit + 1
                case ( 2 )
                  n_emit = this%n_emit%f2(1, this%ix(1,pp), this%ix(2,pp))
                  this%n_emit%f2(1, this%ix(1,pp), this%ix(2,pp)) = n_emit + 1
                case ( 3 )
                  n_emit = this%n_emit%f3(1, this%ix(1,pp), this%ix(2,pp), this%ix(3,pp))
                  this%n_emit%f3(1, this%ix(1,pp), this%ix(2,pp), this%ix(3,pp)) = n_emit + 1
               end select
            endif

            k_damp(i) = real(1.0 - p_frac, p_k_part)
         endif
         pp = pp + 1
      end do

   end subroutine emit_particles_gthr
   !---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Damps particle momentum due to qed emission
!---------------------------------------------------------------------------------------------------
   subroutine damp_qed( this, k_damp, ptrcur, np )
!---------------------------------------------------------------------------------------------------

      implicit none

      class( t_species_qed ),    intent(inout) :: this
      real(p_double), dimension(p_cache_size), intent(in) :: k_damp
      integer, intent(in) :: ptrcur, np

      ! dummy variables
      integer :: i, pp

      ! QED radiation damping
      pp = ptrcur
      do i=1,np
         this%p(:,pp) = this%p(:,pp)*k_damp(i)
         pp = pp + 1
      end do

   end subroutine damp_qed
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Damps particle momentum due to classical radiation reaction
!---------------------------------------------------------------------------------------------------
   subroutine damp_classical( this, ep, bp, dt, gl, um, k_damp, ptrcur, np )
!---------------------------------------------------------------------------------------------------

      use m_random
      use m_species_diagnostics_qed

      implicit none

      class( t_species_qed ),    intent(inout) :: this
      real(p_double), intent(in) :: dt
      real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: ep, bp, um
      real(p_double), dimension(p_cache_size), intent(in) :: gl, k_damp
      integer, intent(in) :: ptrcur, np

      ! dummy variables
      integer :: i, pp
      real(p_double) :: norm_schw, p2, gf, energy_rad
      real(p_k_part), dimension(p_p_dim,p_cache_size) :: un, frad, frad2, frad3, &
         frad2_aux, frad3_aux
      real(p_k_part), dimension(p_cache_size) :: uE, uB, B2, gamn, inv_gamn, gi
      real(p_k_part), dimension(p_p_dim) :: p_g_vec
      real(p_k_part) :: q_g

      ! Values from the NIST reference on Constants, Units and Uncertainty
      ! http://physics.nist.gov/cuu/Constants/index.html
      real( p_double ), parameter :: hbar  = 1.054571800d-34 ! Planck constant over 2 pi
      real( p_double ), parameter :: cl    = 299792458.0d0   ! speed of light in vacuum [m/s]
      real( p_double ), parameter :: me    = 9.10938356d-31  ! electron rest mass [kg]
      real( p_double ), parameter :: alpha = 7.2973525664d-3 ! fine-structure constant

      ! executable statements

      ! calculate the quantum coefficient
      norm_schw = hbar*this%omega_p0/(me*cl**2.d0)

      if (this%ndump_fac_rad > 0) then
         pp = ptrcur
         do i=1,np
            if (gl(i) < this%qed_g_cutoff) then
               gi(i) = sqrt(1.0_p_k_part + this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2)
            endif
            pp = pp + 1
         end do
      endif

      ! get time centered momentum and Lorentz factor
      pp = ptrcur
      do i=1,np
         if (gl(i) < this%qed_g_cutoff) then
            un(:,i) = 0.5_p_k_part*(um(:,i) + this%p(:,pp))
            gamn(i) = sqrt(1.0_p_k_part + un(1,i)**2 + un(2,i)**2 + un(3,i)**2)
            inv_gamn(i)= 1.0_p_k_part / gamn(i)
         endif
         pp = pp + 1
      end do

      ! compute auxiliary dot products
      do i=1,np
         if (gl(i) < this%qed_g_cutoff) then
            uE(i) = dot_product(un(:,i), ep(:,i))
            uB(i) = dot_product(un(:,i), bp(:,i))
            B2(i) = dot_product(bp(:,i), bp(:,i))
         endif
      end do

      do i=1,np
         if (gl(i) < this%qed_g_cutoff) then
            frad2_aux(:,i) = bp(:,i)*uB(i)-un(:,i)*B2(i)+ep(:,i)*uE(i)
            frad3_aux(:,i) = (gamn(i)*ep(1,i)+un(2,i)*bp(3,i)-un(3,i)*bp(2,i))**2
         endif
      end do

      do i=1,np
         if (gl(i) < this%qed_g_cutoff) then
            frad3_aux(:,i) = frad3_aux(:,i) + (gamn(i)*ep(2,i)+un(3,i)*bp(1,i)-un(1,i)*bp(3,i))**2
         endif
      end do

      do i=1,np
         if (gl(i) < this%qed_g_cutoff) then
            frad3_aux(:,i) = frad3_aux(:,i) + (gamn(i)*ep(3,i)+un(1,i)*bp(2,i)-un(2,i)*bp(1,i))**2
         endif
      end do

      do i=1,np
         if (gl(i) < this%qed_g_cutoff) then
            frad3_aux(:,i) = frad3_aux(:,i) - uE(i)**2
         endif
      end do

      do i=1,np
         if (gl(i) < this%qed_g_cutoff) then
            frad3_aux(:,i) = frad3_aux(:,i) * inv_gamn(i)
         endif
      end do

      ! calculate full radiation force from L&L
      pp = ptrcur
      do i=1,np
         if (gl(i) < this%qed_g_cutoff) then
            ! we neglect the first term coming from space time gradients of E and B
            ! second term
            frad2(1,i) = -ep(2,i)*bp(3,i)+ep(3,i)*bp(2,i)- &
               inv_gamn(i)*frad2_aux(1,i)

            frad2(2,i) = -ep(3,i)*bp(1,i)+ep(1,i)*bp(3,i)- &
               inv_gamn(i)*frad2_aux(2,i)

            frad2(3,i) = -ep(1,i)*bp(2,i)+ep(2,i)*bp(1,i)- &
               inv_gamn(i)*frad2_aux(3,i)

            ! third term
            frad3(:,i) = un(:,i)*frad3_aux(:,i)

            frad(:,i) = (0.6666666666666666_p_k_part*alpha*norm_schw)*(frad2(:,i)+frad3(:,i))

            ! calculate classical recoil from L&L
            this%p(:,pp) = this%p(:,pp) - dt*frad(:,i)
         endif
         pp = pp + 1
      end do

      if (this%ndump_fac_rad > 0) then

         energy_rad = 0.0_p_double

         pp = ptrcur
         do i=1,np
            if (gl(i) < this%qed_g_cutoff) then
               ! compute final lorentz factor
               p2 = this%p(1,pp)**2+this%p(2,pp)**2+this%p(3,pp)**2
               gf = sqrt(1.0d0+p2)

               energy_rad = energy_rad + (gi(i) - gf) * this%q(pp)
               if (this%if_add_classical_radsphere) then
                  if (this%ndump_fac_radsphere > 0) then
                     p_g_vec = (this%p(1:3,pp))*(gi(i) - gf) / p2
                     q_g=abs(this%q(pp))
                     call radsphere_hist( this, p_g_vec, q_g )
                  endif
               endif
            endif
            pp = pp + 1
         end do

         ! accumulate radiated energy
         this%energy_rad(1) = this%energy_rad(1) + energy_rad

      endif

   end subroutine damp_classical
!---------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
   subroutine dudt_exact_qed( this, emf, dt, i0, i1, energy, time )
!-----------------------------------------------------------------------------------------
! Advance velocities and calculate time-centered total energy
! process particles in the range [i0,i1]
! This pusher follows the ref. https://arxiv.org/abs/2007.07556
!-----------------------------------------------------------------------------------------

      implicit none

      class( t_species ), intent(inout) :: this
      class( t_emf ), intent( in ) :: emf
      real(p_double), intent(in) :: dt
      integer, intent(in) :: i0, i1
      real(p_double), intent(inout) :: energy
      real(p_double), intent(in) :: time

      integer :: i, j, ptrcur, np, pp
      real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep
      real(p_k_part), dimension(4) :: u0, u, fu0, f2u0, f3u0, uk0, fuk0, uk, uo0, fuo0, uo
      real(p_k_part) :: t, cn_o, snc_o, cnc_o, ch_k, shc_k, chc_k, tmp, k2o2, c1, c2
      real(p_k_part) :: inv1, inv2, e2, b2, fac_rr, qm
      real(p_double) :: loc_ene
      real(p_k_part), dimension(p_cache_size) :: k, o, k2, o2
      integer, dimension(p_cache_size) :: ix1, ix2
      integer :: n1, n2

      real(p_double), save :: rel_tol = 1.0d-12
      real(p_double), save :: weak_fld_tol = 1.0d-12

      ! qed specific variables
      real(p_double), dimension(p_cache_size) :: gl, k_damp
      real(p_double), dimension(p_cache_size) :: g_elec
      real(p_double) :: g_elec_tmp
      real(p_double), dimension(2) :: energy_rad

      select type( this )
       class is ( t_species_qed )

         ! executable statements
         qm = 1.0_p_k_part / this%rqm
         fac_rr = this%k_rr * this%q_real * dt * 0.5_p_k_part * qm

         !loop through all particles
         do ptrcur = i0, i1, p_cache_size

            ! check if last copy of table and set np
            if( ptrcur + p_cache_size > i1 ) then
               np = i1 - ptrcur + 1
            else
               np = p_cache_size
            endif

            call this % get_emf( emf, bp, ep, np, ptrcur, t=time )
            call emit_particles( this, ep, bp, dt, gl, k_damp, ptrcur, np )

            loc_ene = 0.0_p_double
            energy_rad = 0.0_p_double

            pp = ptrcur
            do i = 1, np
               ep(1,i) = ep(1,i) * qm
               ep(2,i) = ep(2,i) * qm
               ep(3,i) = ep(3,i) * qm
               bp(1,i) = bp(1,i) * qm
               bp(2,i) = bp(2,i) * qm
               bp(3,i) = bp(3,i) * qm
               pp = pp + 1
            enddo

            ! -------------------------------------------------------------------------
            ! half-dt push of radiation reaction
            ! -------------------------------------------------------------------------
            if ( this%rad_react ) then
               pp = ptrcur
               do i = 1, np

                  u0(1:3) = this%p(1:3,pp)
                  u0(4) = sqrt( 1.0_p_k_part + u0(1)*u0(1) + u0(2)*u0(2) + u0(3)*u0(3) )

                  fu0  = get_fu( ep(:,i), bp(:,i), u0 )
                  f2u0 = get_fu( ep(:,i), bp(:,i), fu0 )
                  tmp  = dot4( u0, f2u0 )

                  u0(4) = fac_rr / u0(4)

                  ! Energy of the emitted photon
                  g_elec(i) = sqrt( 1.0_p_k_part + u0(1)*u0(1) + u0(2)*u0(2) + u0(3)*u0(3) )

                  ! If the photon energy is lower than the threshold for its emission
                  ! we keep the damping of the exact pusher
                  if(g_elec(i) < this%qed_g_cutoff) then

                     this%p(1,pp) = u0(1) + u0(4) * ( f2u0(1) - tmp * u0(1) )
                     this%p(2,pp) = u0(2) + u0(4) * ( f2u0(2) - tmp * u0(2) )
                     this%p(3,pp) = u0(3) + u0(4) * ( f2u0(3) - tmp * u0(3) )

                  endif

                  pp = pp + 1
               enddo
            endif

            ! -------------------------------------------------------------------------
            ! categorize particles
            ! -------------------------------------------------------------------------
            n1 = 0; ix1 = 0
            n2 = 0; ix2 = 0
            do i = 1, np

               e2 = ep(1,i) * ep(1,i) + ep(2,i) * ep(2,i) + ep(3,i) * ep(3,i)
               b2 = bp(1,i) * bp(1,i) + bp(2,i) * bp(2,i) + bp(3,i) * bp(3,i)

               if ( e2 + b2 < weak_fld_tol ) then

                  k2(i) = 0.0_p_k_part; k(i) = 0.0_p_k_part
                  o2(i) = 0.0_p_k_part; o(i) = 0.0_p_k_part

               else

                  ! Lorentz invariants
                  inv1 = e2 - b2
                  inv2 = ep(1,i) * bp(1,i) + ep(2,i) * bp(2,i) + ep(3,i) * bp(3,i)

                  ! calculate eigenvalues
                  tmp = sqrt( inv1 * inv1 + 4.0_p_k_part * inv2 * inv2 )
                  k2(i) = 0.5_p_k_part * ( inv1 + tmp ); k(i) = sqrt(k2(i))
                  o2(i) = 0.5_p_k_part * (-inv1 + tmp ); o(i) = sqrt(o2(i))

                  if ( k2(i) < rel_tol .and. o2(i) < rel_tol ) then
                     ! light-like
                     n1 = n1 + 1; ix1(n1) = i
                  else
                     ! general case
                     n2 = n2 + 1; ix2(n2) = i
                  endif

               endif
            enddo

            ! -------------------------------------------------------------------------
            ! advance particles in light-like field (category 1)
            ! -------------------------------------------------------------------------
            do j = 1, n1

               i = ix1(j); pp = ptrcur + i - 1

               u0(1:3) = this%p(:,pp)
               u0(4) = sqrt( 1.0_p_k_part + u0(1)*u0(1) + u0(2)*u0(2) + u0(3)*u0(3) )

               fu0  = get_fu( ep(:,i), bp(:,i), u0 )
               f2u0 = get_fu( ep(:,i), bp(:,i), fu0 )
               f3u0 = get_fu( ep(:,i), bp(:,i), f2u0 )

               ! solve for the proper time step according to dt
               call get_t_light( dt, u0(4), fu0(4), f2u0(4), f3u0(4), c1, c2, t, this%iter_tol )
               u = u0 + t * fu0 + c1 * f2u0 + c2 * f3u0

               this%p(:,pp) = u(1:3)

            enddo

            ! -------------------------------------------------------------------------
            ! advance particles in arbitrary field (category 2)
            ! -------------------------------------------------------------------------
            do j = 1, n2

               i = ix2(j); pp = ptrcur + i - 1

               u0(1:3) = this%p(:,pp)
               u0(4) = sqrt( 1.0_p_k_part + u0(1)*u0(1) + u0(2)*u0(2) + u0(3)*u0(3) )

               ! subspace decomposition
               k2o2 = k2(i) + o2(i)
               fu0  = get_fu( ep(:,i), bp(:,i), u0 )
               f2u0 = get_fu( ep(:,i), bp(:,i), fu0 )
               uk0  = ( f2u0 + o2(i) * u0 ) / k2o2
               uo0  = u0 - uk0

               fuk0 = get_fu( ep(:,i), bp(:,i), uk0 )
               fuo0 = fu0 - fuk0

               ! solve for the proper time step according to dt
               call get_t( dt, k(i), uk0(4), fuk0(4), o(i), uo0(4), fuo0(4), &
                  ch_k, shc_k, chc_k, cn_o, snc_o, cnc_o, t, this%iter_tol )

               uk = uk0 * ch_k + fuk0 * shc_k * t
               uo = uo0 * cn_o + fuo0 * snc_o * t
               u = uk + uo

               this%p(1:3,pp) = u(1:3)

            enddo

            ! -------------------------------------------------------------------------
            ! half-dt push of radiation reaction
            ! -------------------------------------------------------------------------
            if ( this%rad_react ) then
               pp = ptrcur
               do i = 1, np
                  u0(1:3) = this%p(1:3,pp)
                  u0(4) = sqrt( 1.0_p_k_part + u0(1) * u0(1) + u0(2) * u0(2) + u0(3) * u0(3) )

                  fu0  = get_fu( ep(:,i), bp(:,i), u0 )
                  f2u0 = get_fu( ep(:,i), bp(:,i), fu0 )
                  tmp  = dot4( u0, f2u0 )

                  u0(4) = fac_rr / u0(4)

                  ! If the photon energy is lower than the threshold for its emission
                  ! we keep the damping of the exact pusher
                  if(g_elec(i) < this%qed_g_cutoff) then
                     this%p(1,pp) = this%p(1,pp) + u0(4) * ( f2u0(1) - tmp * u0(1) )
                     this%p(2,pp) = this%p(2,pp) + u0(4) * ( f2u0(2) - tmp * u0(2) )
                     this%p(3,pp) = this%p(3,pp) + u0(4) * ( f2u0(3) - tmp * u0(3) )

                     ! Diagnostic of energy radiated classically
                     if(this%ndump_fac_rad>0) then
                        g_elec_tmp = sqrt( 1.0_p_k_part + this%p(1,pp) * this%p(1,pp) + this%p(2,pp) * this%p(2,pp) + this%p(3,pp) * this%p(3,pp) )
                        energy_rad(1) = energy_rad(1) + (g_elec(i) - g_elec_tmp) * this%q(pp)
                     endif

                  endif

                  pp = pp + 1
               enddo
            endif

            if(( this%rad_react ) .and. ( this%if_damp_qed ) ) call damp_qed( this, k_damp, ptrcur, np )

            pp = ptrcur
            do i = 1, np
               tmp = sqrt(1.0 + this%p(1,pp)*this%p(1,pp) + this%p(2,pp)*this%p(2,pp) + this%p(3,pp)*this%p(3,pp))
               loc_ene = loc_ene + this%q(pp) * (tmp - 1.0)
               pp = pp + 1
            enddo

            energy = energy + loc_ene

            if( ( this%ndump_fac_rad>0 ) .and. ( this%if_damp_classical ) ) this%energy_rad(1) = this%energy_rad(1) + energy_rad(1)

         enddo

       class default
         ! this must never happen
         call abort_program( p_err_invalid )
      end select

   end subroutine dudt_exact_qed
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
   subroutine get_t_light( dt, u0, fu0, f2u0, f3u0, c1, c2, t, tol )
!-----------------------------------------------------------------------------------------
! Using Newton-Raphson method to calculate proper time step for the light-like case
!-----------------------------------------------------------------------------------------

      implicit none

      real(p_double), intent(in) :: dt, tol
      real(p_k_part), intent(in) :: u0, fu0, f2u0, f3u0
      real(p_k_part), intent(inout) :: t, c1, c2

      real(p_k_part) :: t1, t2, t3, fun, dfun
      integer :: iter

      t1 = dt / u0
      iter = 0

      do

         t2 = t1 * t1
         t3 = t2 * t1
         c1 = 0.5_p_double * t2
         c2 = 0.166666666666667_p_double * t3

         fun  = t1 * u0 + c1 * fu0 + c2 * f2u0 + 0.041666666666667_p_double * t3 * t1 * f3u0 - dt
         dfun = u0 + t * fu0 + c1 * f2u0 + c2 * f3u0

         t = t1 - fun / dfun
         iter = iter + 1

         if ( abs(t - t1) / t1 < tol .or. iter > p_iter_max ) then
            exit
         else
            t1 = t
         endif

      enddo

   end subroutine get_t_light
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
   subroutine get_t( dt, k, uk0, fuk0, o, uo0, fuo0, ch_k, shc_k, chc_k, cn_o, snc_o, cnc_o, t, tol )
!-----------------------------------------------------------------------------------------
! Using Newton-Raphson method to calculate proper time step for general case
!-----------------------------------------------------------------------------------------

      implicit none

      real(p_double), intent(in) :: dt, tol
      real(p_k_part), intent(in) :: k, uk0, fuk0
      real(p_k_part), intent(in) :: o, uo0, fuo0
      real(p_k_part), intent(inout) :: ch_k, shc_k, chc_k, cn_o, snc_o, cnc_o, t

      real(p_k_part) :: fun, dfun, t0
      integer :: iter

      t0 = dt / ( uk0 + uo0 )
      iter = 0

      do

         call  trig_func( o * t0, cn_o, snc_o, cnc_o )
         call hyper_func( k * t0, ch_k, shc_k, chc_k )

         fun = ( uk0 * shc_k + fuk0 * chc_k * t0 ) * t0 + &
            ( uo0 * snc_o - fuo0 * cnc_o * t0 ) * t0 - dt

         dfun = uk0 * ch_k + fuk0 * shc_k * t0 + &
            uo0 * cn_o + fuo0 * snc_o * t0

         t = t0 - fun / dfun
         iter = iter + 1

         if ( abs(t - t0) / t0 < tol .or. iter > p_iter_max ) then
            exit
         else
            t0 = t
         endif

      enddo

   end subroutine get_t
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
   subroutine trig_func( x, cn, snc, cnc, sncc )
!-----------------------------------------------------------------------------------------
! Fast evaluation of the trigonometric functions cos(x), sinc(x), (cos(x)-1)/x^2 and
! (sinc(x)-1)/x^2. When x is smaller than a threshold value, they are evaluated using
! Taylor expansion, otherwise using tan(x/2) for evaluation.
!-----------------------------------------------------------------------------------------
      implicit none
      real(p_k_part), intent(in) :: x
      real(p_k_part), intent(out) :: cn, snc, cnc
      real(p_k_part), intent(out), optional :: sncc

      real(p_double) :: xd, x2, x4, tmp1, tmp2, tmp3, cn_d, snc_d
      ! threshold making Taylor expansion accurate to machine precision
      real(p_double), save :: eps = 0.007367081148666_p_double

      xd = real(x,p_double)
      x2 = xd * xd

      if ( abs(xd) < eps ) then
         x4 = x2 * x2
         cn  = real( 1.0_p_double - 0.500000000000000_p_double * x2 + 0.041666666666667_p_double * x4, p_k_part )
         snc = real( 1.0_p_double - 0.166666666666667_p_double * x2 + 0.008333333333333_p_double * x4, p_k_part )
         cnc =real( -0.5_p_double + 0.041666666666667_p_double * x2 - 0.001388888888889_p_double * x4, p_k_part )
         if ( present(sncc) ) then
            sncc = real( -0.166666666666667_p_double + 0.008333333333333_p_double * x2 - 1.984126984126984d-4 * x4, p_k_part )
         endif
      else
         tmp1 = tan( 0.5_p_double * xd )
         tmp2 = tmp1 * tmp1
         tmp3 = 1.0_p_double / (1.0_p_double + tmp2)
         cn_d = ( 1.0_p_double - tmp2 ) * tmp3
         cn = real( cn_d, p_k_part )
         snc_d = 2.0_p_double * tmp1 * tmp3 / xd
         snc = real( snc_d, p_k_part )
         cnc = real( ( cn_d - 1.0_p_double ) / x2, p_k_part )
         if ( present(sncc) ) then
            sncc = real( ( snc_d - 1.0_p_double ) / x2, p_k_part )
         endif
      endif

   end subroutine trig_func

!-----------------------------------------------------------------------------------------
   subroutine hyper_func( x, ch, shc, chc, shcc )
!-----------------------------------------------------------------------------------------
! Fast evaluation of the hyperbolic functions cosh(x), sinhc(x), (cosh(x)-1)/x^2 and
! (sinhc(x)-1)/x^2. When x is smaller than a threshold value, they are evaluated using
! Taylor expansion, otherwise using tanh(x/2) for evaluation.
!-----------------------------------------------------------------------------------------
      implicit none
      real(p_k_part), intent(in) :: x
      real(p_k_part), intent(out) :: ch, shc, chc
      real(p_k_part), intent(out), optional :: shcc

      real(p_double) :: xd, x2, x4, tmp1, tmp2, tmp3, ch_d, shc_d
      ! threshold making Taylor expansion accurate to machine precision
      real(p_double), save :: eps = 0.007367081148666_p_double

      xd = real(x,p_double)
      x2 = xd * xd

      if ( abs(xd) < eps ) then
         x4 = x2 * x2
         ch  = real( 1.0_p_double + 0.500000000000000_p_double * x2 + 0.041666666666667_p_double * x4, p_k_part )
         shc = real( 1.0_p_double + 0.166666666666667_p_double * x2 + 0.008333333333333_p_double * x4, p_k_part )
         chc = real( 0.5_p_double + 0.041666666666667_p_double * x2 + 0.001388888888889_p_double * x4, p_k_part )
         if ( present(shcc) ) then
            shcc = real( 0.166666666666667_p_double + 0.008333333333333_p_double * x2 + 1.984126984126984d-4 * x4, p_k_part )
         endif
      else
         tmp1 = tanh( 0.5_p_double * xd )
         tmp2 = tmp1 * tmp1
         tmp3 = 1.0_p_double / (1.0_p_double - tmp2)
         ch_d = ( 1.0_p_double + tmp2 ) * tmp3
         ch = real( ch_d, p_k_part )
         shc_d = 2.0_p_double * tmp1 * tmp3 / xd
         shc = real( shc_d, p_k_part )
         chc = real( ( ch_d - 1.0_p_double ) / x2, p_k_part )
         if ( present(shcc) ) then
            shcc = real( ( shc_d - 1.0_p_double ) / x2, p_k_part )
         endif
      endif

   end subroutine hyper_func

!-----------------------------------------------------------------------------------------
   function dot4( u, v )

      implicit none
      real(p_k_part), intent(in), dimension(4) :: u, v
      real(p_k_part) :: dot4

      dot4 = u(4) * v(4) - u(1) * v(1) - u(2) * v(2) - u(3) * v(3)

   end function dot4

!-----------------------------------------------------------------------------------------
   function get_fu( ep, bp, u )
!-----------------------------------------------------------------------------------------
! Calculate the product of field tensor F and a four-vector u
!-----------------------------------------------------------------------------------------
      implicit none
      real(p_k_part), intent(in), dimension(3) :: ep, bp
      real(p_k_part), intent(in), dimension(4) :: u
      real(p_k_part), dimension(4) :: get_fu

      get_fu(4) =                ep(1) * u(1) + ep(2) * u(2) + ep(3) * u(3)
      get_fu(1) = ep(1) * u(4) +                bp(3) * u(2) - bp(2) * u(3)
      get_fu(2) = ep(2) * u(4) - bp(3) * u(1) +                bp(1) * u(3)
      get_fu(3) = ep(3) * u(4) + bp(2) * u(1) - bp(1) * u(2)
   end function get_fu

   subroutine bremsstrahlung_emission(this, dt, k_damp, ptrcur, i0, i1, np)
      use m_parameters

      implicit none

      ! input parameters
      class( t_species_qed ),    intent(inout) :: this
      real(p_double), intent(in) :: dt
      real(p_double), dimension(p_cache_size), intent(in) ::  k_damp

      integer, intent(in) :: ptrcur, np, i0, i1

      ! dummy variables
      integer :: i, pp
      real(p_double) :: proba, varand, ge, gam_k, g_p, w_p, erad
      real(p_k_part), dimension(p_p_dim) :: g_p_vec

      ! Bremsstrahlung variables
      real(p_double) :: factor_mult, omega_p, n_p, ni, momentum_damped, normp
      
      ! Diagnostics
      real(p_k_part), dimension(2) :: energy_rad

      integer :: Z_ion, cs_model
      logical :: energy_damp

      ! Bremsstrahlung variables of interest
      factor_mult = this % bremsstrahlung_info % proba_mult
      omega_p = this%omega_p0
      n_p = this%n0
      Z_ion = int(this % bremsstrahlung_info % Z_ion)
      cs_model = int(this % bremsstrahlung_info % cs_model)
      energy_damp = this % bremsstrahlung_info % energy_damp
      energy_rad = 0.0_p_k_part

      ! loop similar to os-photons
      pp = ptrcur
      do i=1, np

         ge = sqrt(1.0d0 + this%p(1,pp)**2+this%p(2,pp)**2+this%p(3,pp)**2)
         normp = sqrt(this%p(1,pp)**2+this%p(2,pp)**2+this%p(3,pp)**2)

         if (ge .gt. 1.0d0) then

            ! Ion density at the position of the photon
            select case ( p_x_dim )
            case ( 1 )
               ni = this % bremsstrahlung_info % rho_ion % f1(1, this%ix(1, pp))
            case ( 2 )
               ni = this % bremsstrahlung_info % rho_ion % f2(1, this%ix(1, pp), this%ix(2, pp))
            case ( 3 )
               ni = this % bremsstrahlung_info % rho_ion % f3(1, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp))
            end select

            call br_proba(Z_ion, ge, ni, dt, omega_p, n_p, factor_mult, proba)

            call random_number( varand )

            if (varand < proba) then

               call br_spec(Z_ion, ge, gam_k)

               ! photon momentum
               g_p_vec = gam_k * this%p(:,pp) / normp

               ! photon weight
               w_p = abs(this%q(pp)) / factor_mult

               ! Add the photon and update diagnostics
               if (gam_k > this%p_emit_cutoff) then
                  call this%photons%create_particle( &
                     this%ix(1:p_x_dim,pp), &
                     this%x(1:p_x_dim,pp), &
                     g_p_vec, w_p &
                     )
                  energy_rad(2) = energy_rad(2) + gam_k * w_p
               else
                  energy_rad(1) = energy_rad(1) + gam_k * w_p
               endif

               ! Damp the momentum of the electron, via energy balance
               if (energy_damp) then
                  momentum_damped = sqrt((ge - gam_k / factor_mult)**2 - 1.0_p_k_part)
                  this%p(:,pp) = momentum_damped * this%p(:,pp) / normp
               end if

            endif

         endif
         pp = pp + 1
      enddo

      ! accumulate radiated energy
      this % energy_rad(4) = this % energy_rad(4) + energy_rad(1)
      this % energy_rad(5) = this % energy_rad(5) + energy_rad(2)

   end subroutine bremsstrahlung_emission
   subroutine axion_bremsstrahlung_emission(this, dt, k_damp, ptrcur, i0, i1, np)
      use m_parameters

      implicit none

      ! input parameters
      class( t_species_qed ),    intent(inout) :: this
      real(p_double), intent(in) :: dt
      real(p_double), dimension(p_cache_size), intent(in) ::  k_damp

      integer, intent(in) :: ptrcur, np, i0, i1

      ! dummy variables
      integer :: i, pp, k_emit, j_emit
      integer :: ix1, ix2, ix3
      real(p_double) :: proba, ge, ge_new, gam_k, w_p, sum_gam
      real(p_k_part), dimension(p_p_dim) :: g_p_vec

      ! template (photon-clone) variables
      real(p_double) :: factor_mult, factor_mult_eff, omega_p, n_p, ni, momentum_damped, normp
      integer :: Z_ion, cs_model, ax_model
      logical :: energy_damp

      ! model-2 thermal variables (Option A)
      real(p_double) :: g_ae, cs_scale, Te_keV, Te_me, w_axion, X_H, A_ion, ne_factor, lambda_cap
      logical :: include_ee, assume_neutral, use_screening

      ! model-2 physics (Raffelt–Weiss style nondegenerate brems; cgs)
      real(p_double) :: dt_phys, T8, T7, rho_gcc, rho6, rho2, alpha_ae, alpha_26
      real(p_double) :: comp_fac, comp_ei, comp_ee, screen_ei, screen_ee, eps_ei, eps_ee
      real(p_double) :: eps_a, Q_a, kT_erg, meanE_erg, ndot_cm3_s, ne_cm3
      real(p_double) :: rate_per_e_s, lambda_macro, lambda_eff, w_scale
      real(p_double) :: energy_loss_per_e

      ! RNG / angles
      real(p_double) :: u1, u2, u3, u4, mu, phi, st

      ! constants (cgs)
      real(p_double), parameter :: kB_eV_per_K = 8.617333262145d-5
      real(p_double), parameter :: erg_per_eV  = 1.602176634d-12
      real(p_double), parameter :: amu_g       = 1.66053906660d-24
      real(p_double), parameter :: sqrt2       = 1.4142135623730950488d0
      real(p_double), parameter :: tinyu       = 1.0d-300

      ! spectrum proxy parameters for w = E/T (Gamma distribution)
      ! Tuned so mean ~ 1.74 and mode ~ 0.82 (roughly solar free-free shape).
      real(p_double), parameter :: w_shape = 1.891304347826087d0
      real(p_double), parameter :: w_theta = 0.92d0
      real(p_double), parameter :: w_mean  = w_shape*w_theta

      real(p_double) :: pi

      logical, save :: printed_model2 = .false.

      ! Diagnostics
      real(p_k_part), dimension(2) :: energy_rad

      factor_mult = this % axion_bremsstrahlung_info % proba_mult
      omega_p     = this%omega_p0
      n_p         = this%n0
      Z_ion       = int(this % axion_bremsstrahlung_info % Z_ion)
      cs_model    = int(this % axion_bremsstrahlung_info % cs_model)
      energy_damp = this % axion_bremsstrahlung_info % energy_damp

      ax_model    = this % axion_bremsstrahlung_info % ax_model
      g_ae        = this % axion_bremsstrahlung_info % g_ae
      cs_scale    = this % axion_bremsstrahlung_info % cs_scale
      Te_keV      = this % axion_bremsstrahlung_info % Te_keV
      include_ee  = this % axion_bremsstrahlung_info % include_ee
      assume_neutral = this % axion_bremsstrahlung_info % assume_neutral
      use_screening  = this % axion_bremsstrahlung_info % use_screening
      w_axion        = this % axion_bremsstrahlung_info % w_axion

      ! new fields (defaulted in the input reader if not present)
      X_H        = this % axion_bremsstrahlung_info % X_H
      A_ion      = this % axion_bremsstrahlung_info % A_ion
      ne_factor  = this % axion_bremsstrahlung_info % ne_factor
      lambda_cap = this % axion_bremsstrahlung_info % lambda_cap

      ! Feedback scaffolding: allocate per-cell accumulators if requested.
      if ( this % axion_bremsstrahlung_info % apply_feedback .or. &
           this % axion_bremsstrahlung_info % evolve_Te ) then
         call this % axion_bremsstrahlung_info % ensure_feedback()
      end if

      ! If we are doing conservative feedback, do NOT also damp the emitting electron.
      if ( this % axion_bremsstrahlung_info % apply_feedback ) energy_damp = .false.

      if (lambda_cap <= 0.0_p_double) lambda_cap = 20.0_p_double
      if (A_ion <= 0.0_p_double) A_ion = 1.0_p_double
      if (ne_factor <= 0.0_p_double) ne_factor = 1.0_p_double

      energy_rad = 0.0_p_k_part
      pi   = acos(-1.0_p_double)

      ! loop over electrons
      pp = ptrcur
      do i=1, np

         ! Cache local cell indices (used for Te lookup and/or deposits)
         ix1 = this%ix(1, pp)
         if (p_x_dim >= 2) ix2 = this%ix(2, pp)
         if (p_x_dim >= 3) ix3 = this%ix(3, pp)

         ge = sqrt(1.0d0 + this%p(1,pp)**2+this%p(2,pp)**2+this%p(3,pp)**2)
         normp = sqrt(this%p(1,pp)**2+this%p(2,pp)**2+this%p(3,pp)**2)

         ! Ion density at the position of the particle (normalized density in OSIRIS units)
         select case ( p_x_dim )
         case ( 1 )
            ni = this % axion_bremsstrahlung_info % rho_ion % f1(1, this%ix(1, pp))
         case ( 2 )
            ni = this % axion_bremsstrahlung_info % rho_ion % f2(1, this%ix(1, pp), this%ix(2, pp))
         case ( 3 )
            ni = this % axion_bremsstrahlung_info % rho_ion % f3(1, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp))
         end select
         if (ni < 0.0_p_double) ni = 0.0_p_double

         ! -------------------------------------------------------------------
         ! Model 2: Thermal plasma axion bremsstrahlung (Option A)
         !   - Uses a literature-style nondegenerate brems emissivity (Raffelt–Weiss form)
         !   - Converts emissivity -> number rate per electron
         !   - Samples energies ~ O(kT)
         !   - Uses importance sampling with unbiased weight compensation
         ! -------------------------------------------------------------------
         if (ax_model == 2) then

            ! Optional evolving temperature: use local Te_cell (keV) instead of constant Te_keV.
            if ( this % axion_bremsstrahlung_info % evolve_Te ) then
               select case ( p_x_dim )
               case ( 1 )
                  Te_keV = this % axion_bremsstrahlung_info % Te_cell % f1(1, ix1)
               case ( 2 )
                  Te_keV = this % axion_bremsstrahlung_info % Te_cell % f2(1, ix1, ix2)
               case ( 3 )
                  Te_keV = this % axion_bremsstrahlung_info % Te_cell % f3(1, ix1, ix2, ix3)
               end select
            end if

            if (omega_p > 0.0_p_double) then
               dt_phys = dt / omega_p
            else
               dt_phys = 0.0_p_double
            end if

            if ( (dt_phys <= 0.0_p_double) .or. (ni <= 0.0_p_double) .or. &
                 (Te_keV <= 0.0_p_double) .or. (g_ae <= 0.0_p_double) ) then
               k_emit = 0
            else

               ! kT/(mec^2)
               Te_me = (Te_keV * 1.0d3) / 5.11d5

               ! T8 = T / 1e8 K
               T8 = (Te_keV * 1.0d3) / (kB_eV_per_K * 1.0d8)

               ! Local mass density (g/cm^3) from ni*n0 (cm^-3)
               rho_gcc = (ni * n_p) * (A_ion * amu_g)
               rho6    = rho_gcc / 1.0d6

               ! Composition factors (H mass fraction X_H)
               comp_ei  = (1.0_p_double + X_H)
               comp_ee  = (comp_ei*comp_ei) / (2.0_p_double * sqrt2)
               comp_fac = comp_ei
               if (include_ee) comp_fac = comp_fac + comp_ee

               ! alpha_ae = g_ae^2/(4*pi); alpha_26 = alpha_ae / 1e-26
               alpha_ae = (g_ae*g_ae) / (4.0_p_double*pi)
               alpha_26 = alpha_ae / 1.0d-26

               ! Debye–Hückel screening correction (Raffelt 1986, Eqs. 43 & 50)
               screen_ei = 1.0_p_double
               screen_ee = 1.0_p_double
               if (use_screening) then
                  T7   = 10.0_p_double * T8
                  rho2 = rho_gcc / 1.0d2
                  if (T7 > 0.0_p_double) then
                     screen_ei = 1.0_p_double - 6.98d-2 * (3.0_p_double + X_H) * rho2 / (T7*T7)
                     screen_ee = 1.0_p_double - 1.40d-1 * (3.0_p_double + X_H) * rho2 / (T7*T7)
                  end if
                  if (screen_ei < 0.0_p_double) screen_ei = 0.0_p_double
                  if (screen_ee < 0.0_p_double) screen_ee = 0.0_p_double
               end if

               ! Nondegenerate brems energy-loss rate per mass (erg g^-1 s^-1)
               eps_ei = 297.0_p_double * alpha_26 * (T8**2.5_p_double) * rho6 * comp_ei * screen_ei
               eps_ee = 0.0_p_double
               if (include_ee) eps_ee = 297.0_p_double * alpha_26 * (T8**2.5_p_double) * rho6 * comp_ee * screen_ee
               eps_a = eps_ei + eps_ee


! Volumetric emissivity (erg cm^-3 s^-1)
               Q_a = rho_gcc * eps_a

               ! Mean axion energy for number rate conversion: <E> = <w> kT
               kT_erg    = (Te_keV * 1.0d3) * erg_per_eV
               meanE_erg = w_mean * kT_erg
               ndot_cm3_s = Q_a / max(meanE_erg, 1.0d-99)

               ! Local electron density (cm^-3): neutrality by default
               if (assume_neutral) then
                  ne_cm3 = (dble(Z_ion) * ni * n_p) * ne_factor
               else
                  ! fallback: still use neutrality (we can add a true ne deposit later)
                  ne_cm3 = (dble(Z_ion) * ni * n_p) * ne_factor
               end if

               rate_per_e_s = cs_scale * ndot_cm3_s / max(ne_cm3, 1.0d-99)
               if (rate_per_e_s < 0.0_p_double) rate_per_e_s = 0.0_p_double

               ! Choose macro-axion weight (variance control):
               !  - if w_axion > 0 : each emitted macroparticle carries weight w_axion (real axions)
               !  - else : legacy path uses proba_mult with weight |q|/proba_mult
               if (w_axion > 0.0_p_double) then
                  w_p = w_axion
                  lambda_macro = abs(this%q(pp)) * rate_per_e_s * dt_phys / max(w_axion, 1.0d-99)
               else
                  w_p = abs(this%q(pp)) / max(factor_mult, 1.0_p_double)
                  lambda_macro = max(factor_mult, 1.0_p_double) * rate_per_e_s * dt_phys
               end if

               ! Safety cap: keep macro emission count reasonable; compensate weight to preserve expectation
               lambda_eff = lambda_macro
               if (lambda_eff > lambda_cap) then
                  w_scale = lambda_eff / lambda_cap
                  lambda_eff = lambda_cap
                  w_p = w_p * w_scale
               end if

               if (.not. printed_model2) then
                  printed_model2 = .true.
                  write(0,'(a,1pe12.4)') 'ax_model2(RW): omega_p0 [1/s]   = ', omega_p
                  write(0,'(a,1pe12.4)') 'ax_model2(RW): n0 [cm^-3]       = ', n_p
                  write(0,'(a,1pe12.4)') 'ax_model2(RW): ni*n0 [cm^-3]    = ', ni*n_p
                  write(0,'(a,1pe12.4)') 'ax_model2(RW): rho [g/cc]       = ', rho_gcc
                  write(0,'(a,1pe12.4)') 'ax_model2(RW): Te_keV           = ', Te_keV
                  write(0,'(a,1pe12.4)') 'ax_model2(RW): screen_ei        = ', screen_ei
                  write(0,'(a,1pe12.4)') 'ax_model2(RW): screen_ee        = ', screen_ee
                  write(0,'(a,1pe12.4)') 'ax_model2(RW): rate/e [1/s]     = ', rate_per_e_s
                  write(0,'(a,1pe12.4)') 'ax_model2(RW): dt_phys [s]      = ', dt_phys
                  write(0,'(a,1pe12.4)') 'ax_model2(RW): lambda_macro(e0) = ', lambda_macro
                  write(0,'(a,1pe12.4)') 'ax_model2(RW): lambda_cap       = ', lambda_cap
                  write(0,'(a,1pe12.4)') 'ax_model2(RW): w_macro          = ', w_p
               end if

               k_emit = poisson_sample(lambda_eff)

            end if

            if (k_emit > 0) then
               sum_gam = 0.0_p_double

               do j_emit = 1, k_emit

                  ! Sample w = E/T from tuned Gamma proxy; then E/(mec^2) = w * Te_me
                  gam_k = gamma_sample(w_shape, w_theta) * Te_me
                  if (gam_k <= 0.0_p_double) cycle

                  sum_gam = sum_gam + gam_k

                  ! isotropic direction
                  call random_number(u3)
                  call random_number(u4)
                  mu  = 2.0_p_double*u3 - 1.0_p_double
                  phi = 2.0_p_double*pi*u4
                  st  = sqrt(max(0.0_p_double, 1.0_p_double - mu*mu))

                  g_p_vec(1) = gam_k * st * cos(phi)
                  g_p_vec(2) = gam_k * st * sin(phi)
                  g_p_vec(3) = gam_k * mu

                  ! Layer-0 bookkeeping: deposit plasma recoil into per-cell accumulators.
                  ! Convention: emission removes energy/momentum from plasma.
                  if ( this % axion_bremsstrahlung_info % apply_feedback ) then
                     select case ( p_x_dim )
                     case ( 1 )
                        this % axion_bremsstrahlung_info % depE % f1(1, ix1) = &
                             this % axion_bremsstrahlung_info % depE % f1(1, ix1) - w_p * gam_k
                        this % axion_bremsstrahlung_info % depP % f1(1, ix1) = &
                             this % axion_bremsstrahlung_info % depP % f1(1, ix1) - w_p * g_p_vec(1)
                        this % axion_bremsstrahlung_info % depP % f1(2, ix1) = &
                             this % axion_bremsstrahlung_info % depP % f1(2, ix1) - w_p * g_p_vec(2)
                        this % axion_bremsstrahlung_info % depP % f1(3, ix1) = &
                             this % axion_bremsstrahlung_info % depP % f1(3, ix1) - w_p * g_p_vec(3)
                     case ( 2 )
                        this % axion_bremsstrahlung_info % depE % f2(1, ix1, ix2) = &
                             this % axion_bremsstrahlung_info % depE % f2(1, ix1, ix2) - w_p * gam_k
                        this % axion_bremsstrahlung_info % depP % f2(1, ix1, ix2) = &
                             this % axion_bremsstrahlung_info % depP % f2(1, ix1, ix2) - w_p * g_p_vec(1)
                        this % axion_bremsstrahlung_info % depP % f2(2, ix1, ix2) = &
                             this % axion_bremsstrahlung_info % depP % f2(2, ix1, ix2) - w_p * g_p_vec(2)
                        this % axion_bremsstrahlung_info % depP % f2(3, ix1, ix2) = &
                             this % axion_bremsstrahlung_info % depP % f2(3, ix1, ix2) - w_p * g_p_vec(3)
                     case ( 3 )
                        this % axion_bremsstrahlung_info % depE % f3(1, ix1, ix2, ix3) = &
                             this % axion_bremsstrahlung_info % depE % f3(1, ix1, ix2, ix3) - w_p * gam_k
                        this % axion_bremsstrahlung_info % depP % f3(1, ix1, ix2, ix3) = &
                             this % axion_bremsstrahlung_info % depP % f3(1, ix1, ix2, ix3) - w_p * g_p_vec(1)
                        this % axion_bremsstrahlung_info % depP % f3(2, ix1, ix2, ix3) = &
                             this % axion_bremsstrahlung_info % depP % f3(2, ix1, ix2, ix3) - w_p * g_p_vec(2)
                        this % axion_bremsstrahlung_info % depP % f3(3, ix1, ix2, ix3) = &
                             this % axion_bremsstrahlung_info % depP % f3(3, ix1, ix2, ix3) - w_p * g_p_vec(3)
                     end select

                     this % axion_bremsstrahlung_info % fb_dE_step = this % axion_bremsstrahlung_info % fb_dE_step - w_p * gam_k
                     this % axion_bremsstrahlung_info % fb_dP_step(1) = this % axion_bremsstrahlung_info % fb_dP_step(1) - w_p * g_p_vec(1)
                     this % axion_bremsstrahlung_info % fb_dP_step(2) = this % axion_bremsstrahlung_info % fb_dP_step(2) - w_p * g_p_vec(2)
                     this % axion_bremsstrahlung_info % fb_dP_step(3) = this % axion_bremsstrahlung_info % fb_dP_step(3) - w_p * g_p_vec(3)
                     this % axion_bremsstrahlung_info % fb_n_emit_step = this % axion_bremsstrahlung_info % fb_n_emit_step + 1
                  end if

                  if (gam_k > this%p_emit_cutoff) then
                     call this%axions%create_particle( &
                        this%ix(1:p_x_dim,pp), &
                        this%x(1:p_x_dim,pp), &
                        g_p_vec, w_p &
                        )
                     energy_rad(2) = energy_rad(2) + gam_k * w_p
                  else
                     energy_rad(1) = energy_rad(1) + gam_k * w_p
                  endif

               enddo

               ! Optional energy damping of the emitting electron macroparticle
               if (energy_damp .and. (sum_gam > 0.0_p_double)) then
                  energy_loss_per_e = (sum_gam * w_p) / max(abs(this%q(pp)), 1.0d-99)
                  ge_new = ge - energy_loss_per_e
                  if (ge_new < 1.0_p_double) ge_new = 1.0_p_double
                  if (normp > 0.0_p_double) then
                     momentum_damped = sqrt(ge_new**2 - 1.0_p_double)
                     this%p(:,pp) = momentum_damped * this%p(:,pp) / normp
                  end if
               end if

            endif

         else
            ! -------------------------------------------------------------------
            ! Model 0/1: photon-clone template using Coulomb brems tables
            ! -------------------------------------------------------------------
            factor_mult_eff = factor_mult
            if (ax_model == 1) then
               factor_mult_eff = factor_mult * cs_scale * (g_ae**2)
            end if

            call br_proba(Z_ion, ge, ni, dt, omega_p, n_p, factor_mult_eff, proba)

            k_emit = poisson_sample(proba)

            if (k_emit > 0) then

               sum_gam = 0.0_p_double

               do j_emit = 1, k_emit

                  call br_spec(Z_ion, ge, gam_k)

                  ! Prevent unphysical energy loss if the boosted event count is very large
                  if (energy_damp) then
                     if ( (sum_gam + gam_k) / max(factor_mult_eff,1.0_p_double) > (ge - 1.0_p_double) ) then
                        gam_k = (ge - 1.0_p_double) * max(factor_mult_eff,1.0_p_double) - sum_gam
                     end if
                  end if
                  if (gam_k <= 0.0_p_double) exit

                  sum_gam = sum_gam + gam_k

                  ! axion momentum (template: collinear with lepton momentum)
                  if (normp > 0.0_p_double) then
                     g_p_vec = gam_k * this%p(:,pp) / normp
                  else
                     g_p_vec = 0.0_p_k_part
                  end if

                  ! axion weight (template)
                  w_p = abs(this%q(pp)) / max(factor_mult_eff, 1.0_p_double)

                  if (gam_k > this%p_emit_cutoff) then
                     call this%axions%create_particle( &
                        this%ix(1:p_x_dim,pp), &
                        this%x(1:p_x_dim,pp), &
                        g_p_vec, w_p &
                        )
                     energy_rad(2) = energy_rad(2) + gam_k * w_p
                  else
                     energy_rad(1) = energy_rad(1) + gam_k * w_p
                  endif

               enddo

               ! Damp the momentum of the electron, via energy balance (use total emitted energy)
               if (energy_damp .and. (sum_gam > 0.0_p_double)) then
                  ge_new = ge - sum_gam / max(factor_mult_eff,1.0_p_double)
                  if (ge_new < 1.0_p_double) ge_new = 1.0_p_double
                  if (normp > 0.0_p_double) then
                     momentum_damped = sqrt(ge_new**2 - 1.0_p_double)
                     this%p(:,pp) = momentum_damped * this%p(:,pp) / normp
                  end if
               end if

            endif

         endif

         pp = pp + 1
      enddo

      ! accumulate radiated energy
      this % energy_rad(4) = this % energy_rad(4) + energy_rad(1)
      this % energy_rad(5) = this % energy_rad(5) + energy_rad(2)

   contains

      real(p_double) function randn_boxmuller()
         real(p_double) :: a, b
         call random_number(a)
         call random_number(b)
         randn_boxmuller = sqrt(-2.0_p_double * log(max(a, tinyu))) * cos(2.0_p_double*pi*b)
      end function randn_boxmuller

      real(p_double) function gamma_mt(k, theta)
         ! Marsaglia & Tsang method for k >= 1
         real(p_double), intent(in) :: k, theta
         real(p_double) :: d, c, x, v, u

         d = k - 1.0_p_double/3.0_p_double
         c = 1.0_p_double / sqrt(9.0_p_double*d)

         do
            x = randn_boxmuller()
            v = 1.0_p_double + c*x
            if (v <= 0.0_p_double) cycle
            v = v**3
            call random_number(u)
            if (u < 1.0_p_double - 0.0331_p_double*(x**4)) exit
            if (log(max(u, tinyu)) < 0.5_p_double*x*x + d*(1.0_p_double - v + log(v))) exit
         end do

         gamma_mt = theta * d * v
      end function gamma_mt

      real(p_double) function gamma_sample(k, theta)
         ! Gamma(k,theta) sampler, valid for all k>0
         real(p_double), intent(in) :: k, theta
         real(p_double) :: u

         if (k <= 0.0_p_double .or. theta <= 0.0_p_double) then
            gamma_sample = 0.0_p_double
            return
         end if

         if (k >= 1.0_p_double) then
            gamma_sample = gamma_mt(k, theta)
            return
         end if

         ! 0 < k < 1: use boost trick
         call random_number(u)
         gamma_sample = gamma_mt(k + 1.0_p_double, theta) * u**(1.0_p_double/k)

      end function gamma_sample

      integer function poisson_sample(lambda) result(k)
         ! Simple, robust Poisson sampler:
         !  - tiny lambda: Bernoulli shortcut
         !  - moderate lambda: Knuth product method
         !  - large lambda: normal approximation
         real(p_double), intent(in) :: lambda
         real(p_double) :: L, p, u1, u2, z, two_pi, tinyu

         if (lambda <= 0.0_p_double) then
            k = 0
            return
         end if

         if (lambda < 1.0e-3_p_double) then
            call random_number(u1)
            if (u1 < lambda) then
               k = 1
            else
               k = 0
            end if
            return
         end if

         if (lambda < 30.0_p_double) then
            L = exp(-lambda)
            p = 1.0_p_double
            k = -1
            do
               k = k + 1
               call random_number(u1)
               p = p * u1
               if (p <= L) exit
            end do
            return
         end if

         ! Normal approximation for large lambda
         two_pi = 6.28318530717958647692_p_double
         tinyu  = 1.0e-300_p_double
         call random_number(u1)
         call random_number(u2)
         z = sqrt(-2.0_p_double * log(max(u1, tinyu))) * cos(two_pi * u2)

         k = int(lambda + sqrt(lambda) * z + 0.5_p_double)
         if (k < 0) k = 0

      end function poisson_sample

   end subroutine axion_bremsstrahlung_emission
   subroutine axion_compton_emission(this, dt, k_damp, ptrcur, i0, i1, np)

      use m_system
      use m_parameters
      use m_species_define_qed, only : t_species_qed

      implicit none

      class(t_species_qed), intent(inout) :: this
      real(p_double), intent(in) :: dt
      real(p_double), dimension(p_cache_size), intent(in) :: k_damp
      integer, intent(in) :: ptrcur, np, i0, i1

      integer :: pp, k_emit, kk
      integer :: ip
      real(p_double) :: omega_p, dt_phys, n0_cm3
      real(p_double) :: g_ae, cs_scale, Te_keV, X_H, A_ion, ne_factor
      real(p_double) :: proba_mult, w_axion, lambda_cap
      logical :: assume_neutral

      ! constants (cgs)
      real(p_double), parameter :: kB_eV_per_K   = 8.617333262145d-5
      real(p_double), parameter :: erg_per_eV    = 1.602176634d-12
      real(p_double), parameter :: amu_g         = 1.66053906660d-24

      ! spectrum proxy (Gamma in w=E/T): mean 4.9, mode 3.9 -> k=4.9, theta=1
      real(p_double), parameter :: w_shape = 4.9d0
      real(p_double), parameter :: w_theta = 1.0d0
      real(p_double), parameter :: w_mean  = 4.9d0

      real(p_double) :: pi, T7, kT_erg, meanE_erg
      real(p_double) :: ni, rho_gcc, eps_c, Q_a, ndot_cm3_s, ne_cm3
      real(p_double) :: alpha_ae
      real(p_double) :: rate_per_e_s
      real(p_double) :: lambda_macro, lambda_eff, w_scale, w_p

      ! sampled axion kinematics
      real(p_double) :: w, E_keV, gam_k, Te_me
      real(p_double) :: phi, costh, sinth
      real(p_double) :: g_p_vec(3)
      real(p_double) :: tinyu

      integer :: ix1, ix2, ix3

      if (.not. this%axion_compton_info%if_compton) return

      g_ae       = this%axion_compton_info%g_ae
      cs_scale   = this%axion_compton_info%cs_scale
      Te_keV     = this%axion_compton_info%Te_keV
      X_H        = this%axion_compton_info%X_H
      A_ion      = this%axion_compton_info%A_ion
      ne_factor  = this%axion_compton_info%ne_factor
      assume_neutral = this%axion_compton_info%assume_neutral

      proba_mult = this%axion_compton_info%proba_mult
      w_axion    = this%axion_compton_info%w_axion
      lambda_cap = this%axion_compton_info%lambda_cap
      if (lambda_cap <= 0.0_p_double) lambda_cap = 20.0_p_double

      if (Te_keV <= 0.0_p_double) return
      if (g_ae   <= 0.0_p_double) return
      if (cs_scale <= 0.0_p_double) return

      omega_p = this%omega_p0
      if (omega_p <= 0.0_p_double) return
      dt_phys = dt / omega_p
      if (dt_phys <= 0.0_p_double) return

      n0_cm3 = this%n0
      pi = acos(-1.0_p_double)
      tinyu = 1.0e-300_p_double

      ! alpha_ae = g_ae^2/(4*pi)
      alpha_ae = (g_ae*g_ae) / (4.0_p_double*pi)

      ! convert Te_keV -> T7 (T / 1e7 K)
      T7 = (Te_keV * 1.0d3) / (kB_eV_per_K * 1.0d7)

      kT_erg    = (Te_keV * 1.0d3) * erg_per_eV
      meanE_erg = w_mean * kT_erg

      ! dimensionless temperature in mec^2 units for momentum conversion
      Te_me = (Te_keV * 1.0d3) / 5.11d5

      if ( this%axion_compton_info%apply_feedback .or. this%axion_compton_info%evolve_Te ) then
         call this%axion_compton_info%ensure_feedback()
      end if
      

      pp = ptrcur
      

      do ip = 1, np
         ix1 = this%ix(1, pp)
         if (p_x_dim >= 2) ix2 = this%ix(2, pp)
         if (p_x_dim >= 3) ix3 = this%ix(3, pp)

         ! local ion density proxy (normalized)
         select case ( p_x_dim )
         case ( 1 )
            ni = this % axion_compton_info % rho_ion % f1(1, ix1)
         case ( 2 )
            ni = this % axion_compton_info % rho_ion % f2(1, ix1, ix2)
         case ( 3 )
            ni = this % axion_compton_info % rho_ion % f3(1, ix1, ix2, ix3)
         end select
         if (ni <= 0.0_p_double) cycle

         ! evolving Te option
         if (this%axion_compton_info%evolve_Te) then
            select case ( p_x_dim )
            case ( 1 )
               Te_keV = this%axion_compton_info%Te_cell%f1(1, ix1)
            case ( 2 )
               Te_keV = this%axion_compton_info%Te_cell%f2(1, ix1, ix2)
            case ( 3 )
               Te_keV = this%axion_compton_info%Te_cell%f3(1, ix1, ix2, ix3)
            end select
            if (Te_keV <= 0.0_p_double) cycle
            T7 = (Te_keV * 1.0d3) / (kB_eV_per_K * 1.0d7)
            kT_erg    = (Te_keV * 1.0d3) * erg_per_eV
            meanE_erg = w_mean * kT_erg
            Te_me = (Te_keV * 1.0d3) / 5.11d5
         end if

         ! mass density rho (g/cm^3) from ion proxy
         rho_gcc = (ni * n0_cm3) * (A_ion * amu_g)

         ! Raffelt Compton energy-loss rate per mass (erg g^-1 s^-1)
         ! eps_c = 2.67e-2 * (1+X_H) * T7^6 * (alpha_ae / 1.60e-23)
         eps_c = 2.67d-2 * (1.0_p_double + X_H) * (T7**6) * (alpha_ae / 1.60d-23)

         Q_a = rho_gcc * eps_c
         ndot_cm3_s = Q_a / max(meanE_erg, 1.0d-99)

         ! electron density proxy
         ne_cm3 = this%axion_compton_info%Z_ion * (ni * n0_cm3) * ne_factor
         if (ne_cm3 <= 0.0_p_double) cycle

         rate_per_e_s = cs_scale * ndot_cm3_s / ne_cm3
         if (rate_per_e_s <= 0.0_p_double) cycle

         ! Unbiased variance-control weight scheme (matches axion brems):
         if (w_axion > 0.0_p_double) then
            w_p = w_axion
            lambda_macro = abs(this%q(pp)) * rate_per_e_s * dt_phys / max(w_p, 1.0d-99)
         else
            w_p = abs(this%q(pp)) / max(proba_mult, 1.0d-99)
            lambda_macro = proba_mult * rate_per_e_s * dt_phys
         end if

         if (lambda_macro <= 0.0_p_double) cycle

         ! cap Poisson mean but keep expectation by reweighting
         lambda_eff = lambda_macro
         w_scale = 1.0_p_double
         if (lambda_cap > 0.0_p_double .and. lambda_eff > lambda_cap) then
            w_scale = lambda_eff / lambda_cap
            lambda_eff = lambda_cap
         end if

         k_emit = poisson_sample(lambda_eff)
         
         
         if (k_emit <= 0) cycle

         w_p = w_p * w_scale

         do kk = 1, k_emit

            w = gamma_sample(w_shape, w_theta)
            if (w <= 0.0_p_double) cycle

            ! E_keV = w * Te_keV ; gamma_k = E/(mec^2) = E_keV/511
            E_keV = w * Te_keV
            gam_k = (E_keV) / 5.11d2

            ! isotropic emission in electron rest frame
            call random_number(phi)
            call random_number(costh)
            phi = 2.0_p_double*pi*phi
            costh = 2.0_p_double*costh - 1.0_p_double
            sinth = sqrt(max(0.0_p_double, 1.0_p_double - costh*costh))

            g_p_vec(1) = gam_k * sinth * cos(phi)
            g_p_vec(2) = gam_k * sinth * sin(phi)
            g_p_vec(3) = gam_k * costh

            ! create axion (stored in this%axions species)
            call this%axions%create_particle( this%ix(:,pp), this%x(:,pp), g_p_vec, w_p )

            ! feedback deposits (negative to electrons)
            if (this%axion_compton_info%apply_feedback) then
               select case ( p_x_dim )
               case ( 1 )
                  this%axion_compton_info%depE%f1(1, ix1) = this%axion_compton_info%depE%f1(1, ix1) - w_p * gam_k
                  this%axion_compton_info%depP%f1(1, ix1) = this%axion_compton_info%depP%f1(1, ix1) - w_p * g_p_vec(1)
                  this%axion_compton_info%depP%f1(2, ix1) = this%axion_compton_info%depP%f1(2, ix1) - w_p * g_p_vec(2)
                  this%axion_compton_info%depP%f1(3, ix1) = this%axion_compton_info%depP%f1(3, ix1) - w_p * g_p_vec(3)
               case ( 2 )
                  this%axion_compton_info%depE%f2(1, ix1, ix2) = this%axion_compton_info%depE%f2(1, ix1, ix2) - w_p * gam_k
                  this%axion_compton_info%depP%f2(1, ix1, ix2) = this%axion_compton_info%depP%f2(1, ix1, ix2) - w_p * g_p_vec(1)
                  this%axion_compton_info%depP%f2(2, ix1, ix2) = this%axion_compton_info%depP%f2(2, ix1, ix2) - w_p * g_p_vec(2)
                  this%axion_compton_info%depP%f2(3, ix1, ix2) = this%axion_compton_info%depP%f2(3, ix1, ix2) - w_p * g_p_vec(3)
               case ( 3 )
                  this%axion_compton_info%depE%f3(1, ix1, ix2, ix3) = this%axion_compton_info%depE%f3(1, ix1, ix2, ix3) - w_p * gam_k
                  this%axion_compton_info%depP%f3(1, ix1, ix2, ix3) = this%axion_compton_info%depP%f3(1, ix1, ix2, ix3) - w_p * g_p_vec(1)
                  this%axion_compton_info%depP%f3(2, ix1, ix2, ix3) = this%axion_compton_info%depP%f3(2, ix1, ix2, ix3) - w_p * g_p_vec(2)
                  this%axion_compton_info%depP%f3(3, ix1, ix2, ix3) = this%axion_compton_info%depP%f3(3, ix1, ix2, ix3) - w_p * g_p_vec(3)
               end select

               this%axion_compton_info%fb_dE_step = this%axion_compton_info%fb_dE_step - w_p * gam_k
               this%axion_compton_info%fb_dP_step(1) = this%axion_compton_info%fb_dP_step(1) - w_p * g_p_vec(1)
               this%axion_compton_info%fb_dP_step(2) = this%axion_compton_info%fb_dP_step(2) - w_p * g_p_vec(2)
               this%axion_compton_info%fb_dP_step(3) = this%axion_compton_info%fb_dP_step(3) - w_p * g_p_vec(3)
               this%axion_compton_info%fb_n_emit_step = this%axion_compton_info%fb_n_emit_step + 1
            end if

         end do

      pp = pp + 1

      end do
   contains

      real(p_double) function randn_boxmuller()
         real(p_double) :: a, b
         call random_number(a)
         call random_number(b)
         randn_boxmuller = sqrt(-2.0_p_double * log(max(a, tinyu))) * cos(2.0_p_double*pi*b)
      end function randn_boxmuller

      real(p_double) function gamma_mt(k, theta)
         real(p_double), intent(in) :: k, theta
         real(p_double) :: d, c, x, v, u

         d = k - 1.0_p_double/3.0_p_double
         c = 1.0_p_double / sqrt(9.0_p_double*d)

         do
            x = randn_boxmuller()
            v = 1.0_p_double + c*x
            if (v <= 0.0_p_double) cycle
            v = v**3
            call random_number(u)
            if (u < 1.0_p_double - 0.0331_p_double*(x**4)) exit
            if (log(max(u, tinyu)) < 0.5_p_double*x*x + d*(1.0_p_double - v + log(v))) exit
         end do

         gamma_mt = theta * d * v
      end function gamma_mt

      real(p_double) function gamma_sample(k, theta)
         real(p_double), intent(in) :: k, theta
         real(p_double) :: u

         if (k <= 0.0_p_double .or. theta <= 0.0_p_double) then
            gamma_sample = 0.0_p_double
            return
         end if

         if (k >= 1.0_p_double) then
            gamma_sample = gamma_mt(k, theta)
            return
         end if

         call random_number(u)
         gamma_sample = gamma_mt(k + 1.0_p_double, theta) * u**(1.0_p_double/k)

      end function gamma_sample

      integer function poisson_sample(lambda) result(k)
         real(p_double), intent(in) :: lambda
         real(p_double) :: L, p, u1, u2, z, two_pi

         if (lambda <= 0.0_p_double) then
            k = 0
            return
         end if

         if (lambda < 1.0e-3_p_double) then
            call random_number(u1)
            if (u1 < lambda) then
               k = 1
            else
               k = 0
            end if
            return
         end if

         if (lambda < 30.0_p_double) then
            L = exp(-lambda)
            p = 1.0_p_double
            k = -1
            do
               k = k + 1
               call random_number(u1)
               p = p * u1
               if (p <= L) exit
            end do
            return
         end if

         two_pi = 6.28318530717958647692_p_double
         call random_number(u1)
         call random_number(u2)
         z = sqrt(-2.0_p_double * log(max(u1, tinyu))) * cos(two_pi * u2)

         k = int(lambda + sqrt(lambda) * z + 0.5_p_double)
         if (k < 0) k = 0

      end function poisson_sample

   end subroutine axion_compton_emission

end module m_species_push_qed

!-----------------------------------------------------------------------------------------
! Set the correct dudt function for qed
!-----------------------------------------------------------------------------------------
subroutine set_dudt_spec_qed( this  )

   use m_parameters
   use m_species_define_qed, only : t_species_qed
   use m_species_push_qed, only : dudt_boris_qed, dudt_exact_qed

   implicit none

   class( t_species_qed ), intent(inout) :: this

   ! Set pointer to dudt function
   select case ( this % push_type )
    case(p_std)
      this % dudt => dudt_boris_qed
    case(p_exact)
      this % dudt => dudt_exact_qed
    case default
      this % push_type = p_std
      this % dudt => dudt_boris_qed
   end select

end subroutine set_dudt_spec_qed
!-----------------------------------------------------------------------------------------

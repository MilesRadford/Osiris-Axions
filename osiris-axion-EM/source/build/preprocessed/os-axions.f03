# 1 "qed/os-axions.f03"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "qed/os-axions.f03"
# 1 "./os-config.h" 1
! Configuration file for osiris

! ----------------------------------------------------------------------------------------
! Algorithm options
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! System options
! ----------------------------------------------------------------------------------------

! MPI supports MPI_IN_PLACE in global operations
!#define __HAS_MPI_IN_PLACE__

! Use OpenMP
!#define _OPENMP

! Use 1 optimized code
!#define 1
!#define 1
!#define SIMD_AVX
!#define SIMD_BGQ
!#define SIMD_MIC

! Use SION for checkpointing (not available on all systems)
!#define __RST_IO__ = __RST_SION__

! Use log files
!#define __USE_LOG__

! Use MPE for logging and profiling
!#define __USE_MPE__

! Compiler does not fully support the sizeof intrinsic which is a Fortran 2003 feature
!#define __NO_SIZEOF__

! Use the PAPI library for profiling
!#define __USE_PAPI__


! ----------------------------------------------------------------------------------------
! Distribution Options
! ----------------------------------------------------------------------------------------
!
! Optional modules to be removed for distribution
! These are all turned on by default unless the __DISTRO__ preprocessor macro is defined
!
! ----------------------------------------------------------------------------------------



! Include Ionization module


! Include binary collisions module


! Include particle tracking module


! Include perfectly matched layers boundary conditions for EMF


! Include spin advance module which is disabled by default
! #define __HAS_SPIN__

! Include a debug flag for cylindrical modes simulations
! #define __CYL_MODES_DEBUG__

! Compile functions for reporting B fields in the RaDiO module
! #define __HAS_RAD_BFLD__
# 2 "qed/os-axions.f03" 2
# 1 "./os-preprocess.fpp" 1
!
! File: os-preprocess.fpp
!
! A set of preprocessing macros for osiris, and the fpp/cpp preprocessor
!
!



! Macros for the IBM CPP/GNU preprocessor




! Assertion macro



! Debug functions
# 33 "./os-preprocess.fpp"
! SCR hostname functions







! LOG functions




! ERROR functions






! WARNING functions




! functions for restart io
# 3 "qed/os-axions.f03" 2

module m_axions

# 1 "./memory/memory.h" 1
! Include file for the memory module
! Files must include this file using #include "memory/memory.h" rather than just using the module
! the #include must be placed where the module use statement would usually be i.e.
!
! module module2
!
! use module1
! #include "memory/memory.h"
!
# 19 "./memory/memory.h"
use memory
# 7 "qed/os-axions.f03" 2

use m_system
use m_parameters
use m_math

use m_species_define, only: p_cell_near
use m_axions_define, only: t_axions
use m_species_define, only : p_ene_recalc
use m_emf_define, only: t_emf
use m_qed
use m_vdf_define
use m_qed_coulomb
use m_species_charge, only: deposit_rho

implicit none

private

public advance_axions_1d, advance_axions_2d, advance_axions_3d, advance_axions_cyl_2d
public create_pairs, ntrim_qed, create_pairs_betheheitler
public axion_bremsstrahlung_absorb
public axion_compton_absorb
public apply_axion_compton_feedback
public apply_axion_brems_feedback

contains

!---------------------------------------------------------------------------------------------------
subroutine advance_axions_1d( this, gdt, i0, i1)
!---------------------------------------------------------------------------------------------------

    implicit none

    integer, parameter :: rank = 2

    ! dummy variables
    class( t_axions ), intent(inout) :: this
    real(p_double), intent(in) :: gdt
    integer, intent(in) :: i0, i1

    real(p_k_part) :: rdx
    integer :: i, pp, np, ptrcur

    real(p_k_part), dimension(rank,p_cache_size) :: xbuf
    integer, dimension(rank,p_cache_size) :: dxi
    real(p_k_part), dimension(p_cache_size) :: rp, rpp
    real(p_k_part) :: dt

    ! executable statements
    rdx = real( 1.0_p_double/this%dx(1), p_k_part )
    dt = real( gdt, p_k_part )

    ! advance position of photons i0 to i1 in chunks of p_cache_size
    do ptrcur = i0, i1, p_cache_size

        ! check if last copy of table and set np
        if( ptrcur + p_cache_size > i1 ) then
            np = i1 - ptrcur + 1
        else
            np = p_cache_size
        endif

        ! this type of loop is actually faster than
        ! using a forall construct
        pp = ptrcur

        do i=1,np
            rp(i) = 1.0_p_k_part / &
            sqrt( this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2 )
            rpp(i) = dt * rp(i)
        pp = pp + 1
        end do

        ! reference implementation
        pp = ptrcur
        do i=1,np
            xbuf(1,i) = this%x(1,pp) + this%p(1,pp) * rpp(i) * rdx
            dxi(1,i) = ntrim_qed( xbuf(1,i) )

            pp = pp + 1
        end do


        ! copy data from buffer to species data trimming positions
        pp = ptrcur
        do i = 1, np
            this%x(1,pp) = xbuf(1,i) - dxi(1,i)
            this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)

            pp = pp + 1
        end do

    enddo

end subroutine advance_axions_1d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine advance_axions_2d( this, gdt, i0, i1)
!---------------------------------------------------------------------------------------------------

    implicit none

    integer, parameter :: rank = 2

    ! dummy variables
    class( t_axions ), intent(inout) :: this
    real(p_double), intent(in) :: gdt
    integer, intent(in) :: i0, i1

    real(p_k_part), dimension(p_x_dim) :: rdx
    integer :: i, pp, np, ptrcur

    real(p_k_part), dimension(rank,p_cache_size) :: xbuf
    integer, dimension(rank,p_cache_size) :: dxi
    real(p_k_part), dimension(p_cache_size) :: rp, rpp
    real(p_k_part) :: dt

    ! executable statements
    rdx(1) = real( 1.0_p_double/this%dx(1), p_k_part )
    rdx(2) = real( 1.0_p_double/this%dx(2), p_k_part )
    dt = real( gdt, p_k_part )

    ! advance position of photons i0 to i1 in chunks of p_cache_size
    do ptrcur = i0, i1, p_cache_size

        ! check if last copy of table and set np
        if( ptrcur + p_cache_size > i1 ) then
            np = i1 - ptrcur + 1
        else
            np = p_cache_size
        endif

        ! this type of loop is actually faster than
        ! using a forall construct
        pp = ptrcur

        do i=1,np
            rp(i) = 1.0_p_k_part / &
            sqrt( this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2 )
            rpp(i) = dt * rp(i)
        pp = pp + 1
        end do

        ! reference implementation
        pp = ptrcur
        do i=1,np
            xbuf(1,i) = this%x(1,pp) + this%p(1,pp) * rpp(i) * rdx(1)
            xbuf(2,i) = this%x(2,pp) + this%p(2,pp) * rpp(i) * rdx(2)
            dxi(1,i) = ntrim_qed( xbuf(1,i) )
            dxi(2,i) = ntrim_qed( xbuf(2,i) )

            pp = pp + 1
        end do


        ! copy data from buffer to species data trimming positions
        pp = ptrcur
        do i = 1, np
            this%x(1,pp) = xbuf(1,i) - dxi(1,i)
            this%x(2,pp) = xbuf(2,i) - dxi(2,i)
            this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
            this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

            pp = pp + 1
        end do

    enddo

end subroutine advance_axions_2d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine advance_axions_3d( this, gdt, i0, i1)
!---------------------------------------------------------------------------------------------------

    implicit none

    integer, parameter :: rank = 3

    ! dummy variables

    class( t_axions ), intent(inout) :: this
    real(p_double), intent(in) :: gdt
    integer, intent(in) :: i0, i1

    real(p_k_part), dimension(p_x_dim) :: rdx
    integer :: i, pp, np, ptrcur

    real(p_k_part), dimension(rank,p_cache_size) :: xbuf
    integer, dimension(rank,p_cache_size) :: dxi
    real(p_k_part), dimension(p_cache_size) :: rp, rpp
    real(p_k_part) :: dt

    ! executable statements
    rdx(1) = real( 1.0_p_double/this%dx(1), p_k_part )
    rdx(2) = real( 1.0_p_double/this%dx(2), p_k_part )
    rdx(3) = real( 1.0_p_double/this%dx(3), p_k_part )
    dt = real( gdt, p_k_part )

    ! advance position of photons i0 to i1 in chunks of p_cache_size
    do ptrcur = i0, i1, p_cache_size

        ! check if last copy of table and set np
        if( ptrcur + p_cache_size > i1 ) then
            np = i1 - ptrcur + 1
        else
            np = p_cache_size
        endif

        ! this type of loop is actually faster than
        ! using a forall construct
        pp = ptrcur

        do i=1,np
            rp(i) = 1.0_p_k_part / &
            sqrt( this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2 )
            rpp(i) = dt * rp(i)
            pp = pp + 1
        end do

        ! reference implementation
        pp = ptrcur
        do i=1,np
            xbuf(1,i) = this%x(1,pp) + this%p(1,pp) * rpp(i) * rdx(1)
            xbuf(2,i) = this%x(2,pp) + this%p(2,pp) * rpp(i) * rdx(2)
            xbuf(3,i) = this%x(3,pp) + this%p(3,pp) * rpp(i) * rdx(3)
            dxi(1,i) = ntrim_qed( xbuf(1,i) )
            dxi(2,i) = ntrim_qed( xbuf(2,i) )
            dxi(3,i) = ntrim_qed( xbuf(3,i) )

            pp = pp + 1
        end do


        ! copy data from buffer to species data trimming positions
        pp = ptrcur
        do i = 1, np
            this%x(1,pp) = xbuf(1,i) - dxi(1,i)
            this%x(2,pp) = xbuf(2,i) - dxi(2,i)
            this%x(3,pp) = xbuf(3,i) - dxi(3,i)
            this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
            this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)
            this%ix(3,pp) = this%ix(3,pp) + dxi(3,i)

            pp = pp + 1
        end do
    enddo

end subroutine advance_axions_3d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine advance_axions_cyl_2d( this, gdt, i0, i1 )
!---------------------------------------------------------------------------------------------------

    implicit none

    integer, parameter :: rank = 2

    class( t_axions ), intent(inout) :: this
    real(p_double), intent(in) :: gdt
    integer, intent(in) :: i0, i1

    integer :: i, pp, np, ptrcur
    real(p_double), dimension(p_x_dim) :: xmin_g
    real(p_k_part), dimension(rank,p_cache_size) :: xbuf
    integer, dimension(rank,p_cache_size) :: dxi
    real(p_k_part), dimension(p_cache_size) :: rgamma

    real(p_k_part) :: dt_dx1
    integer :: gix2, shift_ix2
    real(p_double) :: dr, rdr, tmp

    real(p_double) :: x2_new, x3_new, r_old, r_new
    real(p_double) :: r_shift

    ! executable statements

    xmin_g(1) = this%g_box( p_lower, 1 )
    xmin_g(2) = this%g_box( p_lower, 2 )

    shift_ix2 = this%my_nx_p(p_lower, 2) - 2
    dr = this%dx(p_r_dim)
    rdr = 1.0_p_double/dr

    dt_dx1 = real( gdt / this%dx(1), p_k_part )

    ! Since r = 0 is at the center of cell ix2 = 1 the radial position of particles will be
    ! - Positions defined with regard to the center of the cell (odd interpolation)
    ! r = ( (gix2-1) + x2 ) * dr
    ! - Positions defined with regard to the corner of the cell (even interpolation)
    ! r = ( (gix2-1) + x2 - 0.5 ) * dr

    if ( this%pos_type == p_cell_near ) then
        r_shift = 0.5_p_double
    else
        r_shift = 0.0_p_double
    endif

    ! advance position
    do ptrcur = i0, i1, p_cache_size

        ! check if last copy of table and set np
        if( ptrcur + p_cache_size > i1 ) then
            np = i1 - ptrcur + 1
        else
            np = p_cache_size
        endif

        ! this type of loop is actually faster than
        ! using a forall construct
        pp = ptrcur

        do i=1,np
            rgamma(i) = 1.0_p_k_part / &
            sqrt( ( ( this%p(1,pp)**2 ) + this%p(2,pp)**2 ) + this%p(3,pp)**2 )

            pp = pp + 1
        end do

        ! advance particle position
        pp = ptrcur

        do i=1,np
            xbuf(1,i) = this%x(1,pp) + ( this%p(1,pp) * rgamma(i) ) * dt_dx1

            ! Convert radial "cell" position to "box" position in double precision
            gix2 = this%ix(2,pp) + shift_ix2
            r_old = ( this%x(2,pp) + (gix2 - r_shift) ) * dr

            ! Particle is pushed in 3D like fashion
            x2_new = r_old + ( this%p(2,pp) * rgamma(i) ) * gdt
            x3_new = ( this%p(3,pp) * rgamma(i) ) * gdt

            r_new = sqrt( x2_new**2 + x3_new**2 )

            ! Convert new position to "cell" type position

            ! this is a protection against roundoff for cold plasmas
            if ( r_old == r_new ) then
            xbuf(2,i) = this%x(2,pp)
            else
            xbuf(2,i) = real( r_new * rdr - ( gix2 - r_shift ) , p_k_part)
            endif

            ! Correct p_r and p_\theta to conserve angular momentum
            ! there is a potential division by zero here
            tmp = 1.0_p_double / r_new
            this%p(2,pp) = real( ( this%p(2,pp) * x2_new + this%p(3,pp)*x3_new ) * tmp, p_k_part )
            this%p(3,pp) = real( this%p(3,pp) * ( r_old * tmp ), p_k_part )

            dxi(1,i) = ntrim_qed( xbuf(1,i) )
            dxi(2,i) = ntrim_qed( xbuf(2,i) )

            pp = pp + 1
        end do

        ! copy data from buffer to species data trimming positions
        pp = ptrcur
        do i = 1, np
        this%x(1,pp) = xbuf(1,i) - dxi(1,i)
        this%x(2,pp) = xbuf(2,i) - dxi(2,i)
        this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
        this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

        pp = pp + 1
        end do

    enddo


end subroutine advance_axions_cyl_2d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function ntrim_qed(x)
!---------------------------------------------------------------------------------------------------
! Returns the integer shift (-1, 0 or +1) so that the coordinate remains in the [-0.5, 0.5[
! range. This is the fastest implementation (twice as fast as a sequence of ifs) because
! the two if structures compile as conditional moves and can be processed independently.
! This has no precision problem and is only 12% slower than the previous "int(x+1.5)-1"
! routine that would break for x = nearest( 0.5, -1.0 )
!---------------------------------------------------------------------------------------------------
    implicit none

    real(p_k_part), intent(in) :: x
    integer :: ntrim_qed, a, b

    if ( x < -.5 ) then
    a = -1
    else
    a = 0
    endif

    if ( x >= .5 ) then
    b = +1
    else
    b = 0
    endif

    ntrim_qed = a+b

end function ntrim_qed
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine create_pairs( this, electrons, positrons, emf, dt, i0, i1)
!---------------------------------------------------------------------------------------------------

    use m_species_define_qed, only : t_species_qed
    use m_species_comm, only : t_part_idx, remove_particles
    use m_axions_diagnostics
    !use m_cross_section

    implicit none

    integer, parameter :: rank = 2

    ! dummy variables

    class( t_axions ), intent(inout) :: this
    class( t_species_qed ), intent(inout) :: electrons
    class( t_species_qed ), intent(inout) :: positrons
    class( t_emf ), intent( in ) :: emf
    real(p_double), intent(in) :: dt
    integer, intent(in) :: i0, i1

    real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep
    integer :: i, pp, np, ptrcur

    ! attributes for new particles
    real(p_k_part) , dimension(p_p_dim) :: p_new

    ! index of photons to remove
    type(t_part_idx) :: local_idx

    ! QED local variables
    real(p_double) :: chi_g, chi_g_s
    real(p_double) :: coef_QED, omega_p, varand, Wpair, eta, interp, norm_schw, &
    enephot, pdotE2, enephotE_plus_pcrossB2, gam_e, gam_p

    ! dummy variables for chi parameter diagnostic
    real(p_k_part) :: q_g
    integer, dimension(p_x_dim) :: ix_g

    ! Physical constants taken from NIST reference data
    real(p_double), parameter :: hbar = 1.054571726d-34 ! J s
    real(p_double), parameter :: me = 9.10938291d-31 ! kg
    real(p_double), parameter :: cl = 299792458.0 ! m/s
    real(p_double), parameter :: alpha = 7.2973525698d-3 ! (no units)

    real(p_double), parameter :: PI = 3.141592653589793d0

    ! executable statements
    omega_p = electrons%omega_p0
    norm_schw = hbar * omega_p/(me * cl**2)
    coef_QED = alpha/(sqrt(3.d0)*PI*norm_schw)

    ! advance position of photons i0 to i1 in chunks of p_cache_size
    do ptrcur = i0, i1, p_cache_size

        ! check if last copy of table and set np
        if( ptrcur + p_cache_size > i1 ) then
            np = i1 - ptrcur + 1
        else
            np = p_cache_size
        endif

        call this % get_emf( emf, bp, ep, np, ptrcur )

        ! this type of loop is actually faster than
        ! using a forall construct
        pp = ptrcur
        do i=1, np
            ! calculate the quantum parameter chi for all photons
            pdotE2=((this%p(1,pp)*ep(1,i))+(this%p(2,pp)*ep(2,i))+(this%p(3,pp)*ep(3,i)))**2

            enephot = sqrt(this%p(1,pp)**2+this%p(2,pp)**2+this%p(3,pp)**2)

            enephotE_plus_pcrossB2=(-enephot*ep(1,i)-bp(3,i)*this%p(2,pp)+bp(2,i)*this%p(3,pp))**2 &
            + (-enephot*ep(2,i)+bp(3,i)*this%p(1,pp)-bp(1,i)*this%p(3,pp))**2 &
            + (-enephot*ep(3,i)-bp(2,i)*this%p(1,pp)+bp(1,i)*this%p(2,pp))**2

            chi_g = sqrt(abs(enephotE_plus_pcrossB2-pdotE2))*norm_schw
            chi_g_s = chi_g * this%qed_zeta

            ! condition to radiate significant number of pairs
            if ( (chi_g_s > this%pp_chi_cutoff) .and. (chi_g_s < 5.0d4) .and.(enephot > 2.0d0) ) then

                interp = interp_spec_pairQED(chi_g_s)
                Wpair = coef_QED*interp/enephot

                call random_number( varand )

                if ( varand < Wpair * dt ) then
                    ! check to be sure that the pair created have at least each an energy > mc^2
                    gam_e = 0.0d0
                    gam_p = 0.0d0

                    eta = specpair_fast(chi_g_s)
                    gam_e = enephot*eta/chi_g_s
                    gam_p = enephot - gam_e

                    q_g = this%q(pp)
                    ix_g = this%ix(1:p_x_dim,pp)

                    ! Diagnostics
                    if (this%ndump_fac_chi > 0) call chi_hist( this, ix_g, real(q_g, p_k_part), real(chi_g_s, p_k_part) )

                    if ( (gam_e.gt.1.0d0) .and. (gam_p.gt.1.0d0) ) then

                        this%num_pairs = this%num_pairs + 1
                        this%weight_pairs = this%weight_pairs + q_g
                        this%energy_pairs = this%energy_pairs + q_g*(gam_e + gam_p)

                        ! create electron
                        p_new = real( sqrt((gam_e)**2-1)*this%p(:,pp)/enephot, p_k_part )
                        call electrons%create_particle( this%ix(:,pp), this%x(:,pp), p_new, -1.0_p_k_part*q_g )

                        ! create positron
                        p_new = real( sqrt((gam_p)**2-1)*this%p(:,pp)/enephot, p_k_part )
                        call positrons%create_particle( this%ix(:,pp), this%x(:,pp), p_new, +1.0_p_k_part*q_g )

                        ! register emission
                        if (this%ndump_fac_n_emit > 0) call n_emit_hist( this, this%ix(1:p_x_dim,pp) )

                    endif

                    ! mark photon for deletion
                    call add_particle_index( local_idx, pp )

                endif
            endif
            pp = pp + 1
        enddo
    enddo

    ! delete photons
    call remove_particles( this%t_species , local_idx )

    ! cleanup the local_idx buffer for deleted photons
    if ( local_idx%buf_size > 0 ) call freemem(local_idx%idx,"qed/os-axions.f03",550)
    local_idx%start = 1
    local_idx%nidx = 0
    local_idx%buf_size = 0

    contains

    !---------------------------------------------------------------------------------------------------
    subroutine add_particle_index( list, i )

        use m_species_define

        implicit none

        type(t_part_idx), intent(inout) :: list
        integer, intent(in) :: i

        integer, dimension(:), pointer :: tmp_idx => null()

        ! grow buffer if necessary
        if ( list%nidx + 1 > list%buf_size ) then
            list%buf_size = list%buf_size + p_spec_buf_block
            call alloc(tmp_idx, (/ list%buf_size /),"qed/os-axions.f03",572)
            if ( list%nidx > 0 ) call memcpy( tmp_idx, list%idx, list%nidx )
            if (associated(list%idx)) call freemem(list%idx,"qed/os-axions.f03",574)
            list%idx => tmp_idx
        endif

        ! add particle index to send index list
        list%nidx = list%nidx + 1
        list%idx( list%nidx ) = i

    end subroutine add_particle_index
    !---------------------------------------------------------------------------------------------------

end subroutine create_pairs
!---------------------------------------------------------------------------------------------------

  subroutine create_pairs_betheheitler( this, electrons, positrons, dt, i0, i1)

    use m_species_define_qed, only : t_species_qed
    use m_species_define, only : t_species
    use m_species_comm, only : t_part_idx, remove_particles
    use m_axions_diagnostics

    implicit none

    ! dummy variables

    class( t_axions ), intent(inout) :: this
    class( t_species_qed ), intent(inout) :: electrons
    class( t_species_qed ), intent(inout) :: positrons
    real(p_double), intent(in) :: dt
    integer, intent(in) :: i0, i1

    integer :: i, pp, np, ptrcur

    ! attributes for new particles
    real(p_k_part) , dimension(p_p_dim) :: p_new

    ! index of photons to remove
    type(t_part_idx) :: local_idx
    real(p_double) :: proba, varand, enephot, gam_e, gam_p

    ! Bethe Heitler variables
    real(p_double) :: factor_mult, omega_p, n_p, ni
    real(p_k_part) :: q_pair

    integer :: Z_ion, cs_model
    logical :: decay_photons

    ! We introduce shorter names for all Bethe Heitler variables of interest
    factor_mult = this % betheheitler_info % proba_mult
    omega_p = this%omega_p0
    n_p = this%n0
    Z_ion = int(this % betheheitler_info % Z_ion)
    cs_model = int(this % betheheitler_info % cs_model)
    decay_photons = this % betheheitler_info % decay_photons

    ! advance position of photons i0 to i1 in chunks of p_cache_size
    do ptrcur = i0, i1, p_cache_size

        ! check if last copy of table and set np
        if( ptrcur + p_cache_size > i1 ) then
            np = i1 - ptrcur + 1
        else
            np = p_cache_size
        endif

        ! this type of loop is actually faster than
        ! using a forall construct
        pp = ptrcur
        do i=1, np

                ! Photon energy
                enephot = sqrt(this%p(1,pp)**2+this%p(2,pp)**2+this%p(3,pp)**2)

                ! The photon energy must be at least 2 mc2 for pair production
                if (enephot .gt. 2.0d0) then

                    ! Ion density at the position of the photon
                    select case ( p_x_dim )
                    case ( 1 )
                        ni = this % betheheitler_info % rho_ion % f1(1, this%ix(1, pp))
                    case ( 2 )
                        ni = this % betheheitler_info % rho_ion % f2(1, this%ix(1, pp), this%ix(2, pp))
                    case ( 3 )
                        ni = this % betheheitler_info % rho_ion % f3(1, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp))
                    end select

                    call bh_proba(Z_ion, enephot, ni, dt, omega_p, n_p, factor_mult, proba)
                    call random_number( varand )

                    if ( varand < proba ) then
                        ! check to be sure that the pair created have at least each an energy > mc^2
                        gam_e = 0.0d0
                        gam_p = 0.0d0

                        ! Sample positron energy and deduce electron energy
                        call bh_spec(Z_ion, enephot, gam_p)
                        gam_e = enephot - gam_p
                        q_pair = this%q(pp) / factor_mult

                        ! Diagnostics
                        if ( (gam_e.gt.1.0d0) .and. (gam_p.gt.1.0d0) ) then

                            this%num_pairs_bh = this%num_pairs_bh + 1
                            this%weight_pairs_bh = this%weight_pairs_bh + q_pair
                            this%energy_pairs_bh = this%energy_pairs_bh + q_pair*(gam_e + gam_p)

                            ! create electron
                            p_new = real( sqrt((gam_e)**2-1)*this%p(:,pp)/enephot, p_k_part )
                            call electrons%create_particle_single_cell_p( this%ix(:,pp), this%x(:,pp), p_new, -1.0_p_k_part*q_pair )

                            ! create positron
                            p_new = real( sqrt((gam_p)**2-1)*this%p(:,pp)/enephot, p_k_part )
                            call positrons%create_particle_single_cell_p( this%ix(:,pp), this%x(:,pp), p_new, 1.0_p_k_part*q_pair )

                        endif

                        ! Decay of the photon (if we increased proba, we don't decay the photon, but just decrease its weight) 
                        if (decay_photons) then
                            if (factor_mult .eq. 1.0_p_double) then
                                call add_particle_index( local_idx, pp )
                            else
                                this%q(pp) = this%q(pp) * (1.0_p_k_part - (gam_e + gam_p - 2.)/(factor_mult*enephot))
                            endif
                        endif

                    endif
                endif
                pp = pp + 1
        enddo
    enddo

    ! delete photons
    if (decay_photons) call remove_particles( this%t_species , local_idx )

    ! cleanup the local_idx buffer for deleted photons
    if ( local_idx%buf_size > 0 ) call freemem(local_idx%idx,"qed/os-axions.f03",709)
    local_idx%start = 1
    local_idx%nidx = 0
    local_idx%buf_size = 0

    contains

    !---------------------------------------------------------------------------------------------------
    subroutine add_particle_index( list, i )

        use m_species_define

        implicit none

        type(t_part_idx), intent(inout) :: list
        integer, intent(in) :: i

        integer, dimension(:), pointer :: tmp_idx => null()

        ! grow buffer if necessary
        if ( list%nidx + 1 > list%buf_size ) then
            list%buf_size = list%buf_size + p_spec_buf_block
            call alloc(tmp_idx, (/ list%buf_size /),"qed/os-axions.f03",731)
            if ( list%nidx > 0 ) call memcpy( tmp_idx, list%idx, list%nidx )
            if (associated(list%idx)) call freemem(list%idx,"qed/os-axions.f03",733)
            list%idx => tmp_idx
        endif

        ! add particle index to send index list
        list%nidx = list%nidx + 1
        list%idx( list%nidx ) = i

    end subroutine add_particle_index
    !---------------------------------------------------------------------------------------------------


  end subroutine create_pairs_betheheitler

!---------------------------------------------------------------------------------------------------
! Inverse axion bremsstrahlung (absorption) for thermal plasma model (ax_model = 2)
!
! We compute an absorption rate per axion energy using detailed balance:
! d\dot n_a/dE = n_e * (d\Gamma_em/dE)
! (dn_a^{eq}/dE) * \Gamma_abs(E) = cs_scale * d\dot n_a/dE
! so that in equilibrium the net production vanishes bin-by-bin in energy.
!
! The differential emissivity is constructed from the same local RW-style
! emissivity model used in the emission routine, and the same proxy spectrum
! pdf in w = E/T (Gamma distribution) used to sample axion energies.
!---------------------------------------------------------------------------------------------------
subroutine axion_bremsstrahlung_absorb( this, electrons, gdt, i0, i1 )

    use m_system
    use m_parameters
    use m_axions_define, only : t_axions
    use m_species_define_qed, only : t_species_qed
    use m_species_comm, only : t_part_idx, remove_particles

    implicit none

    class( t_axions ), intent(inout) :: this
    class( t_species_qed ), intent(inout) :: electrons
    real(p_double), intent(in) :: gdt
    integer, intent(in) :: i0, i1

    integer :: pp, i1_use
    real(p_double) :: omega_p, dt_phys, n0_cm3

    ! local plasma / model parameters
    real(p_double) :: Te_keV, g_ae, cs_scale
    real(p_double) :: X_H, A_ion
    logical :: include_ee, use_screening, use_bose
    logical :: apply_feedback, evolve_Te

    ! constants (cgs)
    real(p_double), parameter :: kB_eV_per_K = 8.617333262145d-5
    real(p_double), parameter :: erg_per_eV = 1.602176634d-12
    real(p_double), parameter :: amu_g = 1.66053906660d-24
    real(p_double), parameter :: hbar_erg_s = 1.054571817d-27
    real(p_double), parameter :: c_cms = 2.99792458d10
    real(p_double), parameter :: hbarc_erg_cm = hbar_erg_s * c_cms

    ! spectrum proxy parameters for w = E/T (Gamma distribution)
    real(p_double), parameter :: w_shape = 1.891304347826087d0
    real(p_double), parameter :: w_theta = 0.92d0
    real(p_double), parameter :: w_mean = w_shape*w_theta
    real(p_double), parameter :: w_norm = 1.2211041230563731d0 ! 1/(Gamma(k) theta^k)

    real(p_double) :: pi
    real(p_double) :: T8, T7, rho_gcc, rho6, rho2
    real(p_double) :: alpha_ae, alpha_26
    real(p_double) :: comp_ei, comp_ee, screen_ei, screen_ee
    real(p_double) :: eps_ei, eps_ee, eps_a, Q_a
    real(p_double) :: kT_erg, meanE_erg, ndot_cm3_s
    real(p_double) :: Te_use, T8_use, T7_use, kT_use, meanE_use

    ! particle energy / spectrum
    real(p_double) :: p2, E_me, w, pdf_w
    real(p_double) :: E_erg
    real(p_double) :: dn_eq_dE, dndot_dE, gamma_abs, tau, p_abs, u
    real(p_double) :: pref_ne

    ! local ion density (normalized)
    real(p_double) :: ni

    type( t_part_idx ) :: local_idx

    ! Defensive initialization (some compilers/builds may not guarantee default init)
    local_idx%nidx = 0
    local_idx%buf_size = 0
    nullify(local_idx%idx)

    if (.not. electrons%axion_bremsstrahlung_info%if_bremsstrahlung) return
    if (electrons%axion_bremsstrahlung_info%ax_model /= 2) return
    if (.not. electrons%axion_bremsstrahlung_info%if_inverse) return

    Te_keV = electrons%axion_bremsstrahlung_info%Te_keV
    g_ae = electrons%axion_bremsstrahlung_info%g_ae
    cs_scale = electrons%axion_bremsstrahlung_info%cs_scale
    include_ee = electrons%axion_bremsstrahlung_info%include_ee
    use_screening = electrons%axion_bremsstrahlung_info%use_screening
    use_bose = electrons%axion_bremsstrahlung_info%use_bose
    X_H = electrons%axion_bremsstrahlung_info%X_H
    A_ion = electrons%axion_bremsstrahlung_info%A_ion

    apply_feedback = electrons%axion_bremsstrahlung_info%apply_feedback
    evolve_Te = electrons%axion_bremsstrahlung_info%evolve_Te

    ! Ensure feedback/evolving-Te storage is allocated (does nothing if already allocated)
    if (apply_feedback) call electrons%axion_bremsstrahlung_info%ensure_feedback()

    if (Te_keV <= 0.0_p_double) return
    if (g_ae <= 0.0_p_double) return
    if (cs_scale <= 0.0_p_double) return

    omega_p = electrons%omega_p0
    if (omega_p <= 0.0_p_double) return
    dt_phys = gdt / omega_p
    if (dt_phys <= 0.0_p_double) return

    n0_cm3 = electrons%n0

    pi = acos(-1.0_p_double)

    ! T8 = T / 1e8 K
    T8 = (Te_keV * 1.0d3) / (kB_eV_per_K * 1.0d8)
    T7 = 10.0_p_double * T8

    ! mean energy conversion constants
    kT_erg = (Te_keV * 1.0d3) * erg_per_eV
    meanE_erg = w_mean * kT_erg

    ! alpha_ae = g_ae^2/(4*pi); alpha_26 = alpha_ae / 1e-26
    alpha_ae = (g_ae*g_ae) / (4.0_p_double*pi)
    alpha_26 = alpha_ae / 1.0d-26

    ! equilibrium dn/dE prefactor: E^2 / [ 2 pi^2 (hbar c)^3 ]
    pref_ne = 1.0_p_double / (2.0_p_double * pi*pi * (hbarc_erg_cm**3))

    ! init delete list
    local_idx%start = 1
    local_idx%nidx = 0
    local_idx%buf_size = 0
    local_idx%idx => null()

    ! Defensive clamp: this routine may be called with stale i1 if other operators
    ! removed axions earlier in the same timestep.
    i1_use = min(i1, this%num_par)
    if (i0 > i1_use) return

    do pp = i0, i1_use

        ! local ion density at axion location (OSIRIS normalized units)
        select case ( p_x_dim )
        case ( 1 )
            ni = electrons % axion_bremsstrahlung_info % rho_ion % f1(1, this%ix(1, pp))
        case ( 2 )
            ni = electrons % axion_bremsstrahlung_info % rho_ion % f2(1, this%ix(1, pp), this%ix(2, pp))
        case ( 3 )
            ni = electrons % axion_bremsstrahlung_info % rho_ion % f3(1, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp))
        end select
        if (ni < 0.0_p_double) ni = 0.0_p_double
        if (ni <= 0.0_p_double) cycle

        ! local temperature (keV) optionally from evolving field
        Te_use = Te_keV
        T8_use = T8
        T7_use = T7
        kT_use = kT_erg
        meanE_use = meanE_erg

        if (evolve_Te) then
            select case ( p_x_dim )
            case ( 1 )
                if (associated(electrons%axion_bremsstrahlung_info%Te_cell%f1)) then
                    Te_use = electrons%axion_bremsstrahlung_info%Te_cell%f1(1, this%ix(1, pp))
                end if
            case ( 2 )
                if (associated(electrons%axion_bremsstrahlung_info%Te_cell%f2)) then
                    Te_use = electrons%axion_bremsstrahlung_info%Te_cell%f2(1, this%ix(1, pp), this%ix(2, pp))
                end if
            case ( 3 )
                if (associated(electrons%axion_bremsstrahlung_info%Te_cell%f3)) then
                    Te_use = electrons%axion_bremsstrahlung_info%Te_cell%f3(1, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp))
                end if
            end select
            if (Te_use <= 0.0_p_double) cycle

            T8_use = (Te_use * 1.0d3) / (kB_eV_per_K * 1.0d8)
            T7_use = 10.0_p_double * T8_use
            kT_use = (Te_use * 1.0d3) * erg_per_eV
            meanE_use = w_mean * kT_use
        end if

        ! axion energy (dimensionless E/(mec^2) = |p|)
        p2 = this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2
        if (p2 <= 0.0_p_double) cycle
        E_me = sqrt(p2)

        ! w = E/T with E in keV, T in keV; E = E_me * 511 keV
        w = (E_me * 5.11d2) / Te_use
        if (w <= 0.0_p_double) cycle

        ! proxy emission spectrum pdf in w
        if (w < 1.0d-300) cycle
        pdf_w = w_norm * (w**(w_shape - 1.0_p_double)) * exp( - w / w_theta )

        ! local mass density (g/cm^3)
        rho_gcc = (ni * n0_cm3) * (A_ion * amu_g)
        rho6 = rho_gcc / 1.0d6

        ! composition factors (H mass fraction X_H)
        comp_ei = (1.0_p_double + X_H)
        comp_ee = (comp_ei*comp_ei) / (2.0_p_double * 1.4142135623730950488d0)

        ! Debye-Huckel screening correction
        screen_ei = 1.0_p_double
        screen_ee = 1.0_p_double
        if (use_screening) then
            rho2 = rho_gcc / 1.0d2
            if (T7_use > 0.0_p_double) then
                screen_ei = 1.0_p_double - 6.98d-2 * (3.0_p_double + X_H) * rho2 / (T7_use*T7_use)
                screen_ee = 1.0_p_double - 1.40d-1 * (3.0_p_double + X_H) * rho2 / (T7_use*T7_use)
            end if
            if (screen_ei < 0.0_p_double) screen_ei = 0.0_p_double
            if (screen_ee < 0.0_p_double) screen_ee = 0.0_p_double
        end if

        ! RW-style nondegenerate brems energy-loss rate per mass (erg g^-1 s^-1)
        eps_ei = 297.0_p_double * alpha_26 * (T8_use**2.5_p_double) * rho6 * comp_ei * screen_ei
        eps_ee = 0.0_p_double
        if (include_ee) eps_ee = 297.0_p_double * alpha_26 * (T8_use**2.5_p_double) * rho6 * comp_ee * screen_ee
        eps_a = eps_ei + eps_ee

        ! volumetric emissivity and total number emissivity
        Q_a = rho_gcc * eps_a
        ndot_cm3_s = Q_a / max(meanE_use, 1.0d-99)

        ! differential number emissivity per energy: d\dot n/dE = (ndot) * pdf(w) / (kT)
        dndot_dE = ndot_cm3_s * pdf_w / max(kT_use, 1.0d-99)

        ! equilibrium axion number density per energy
        E_erg = w * kT_use
        if (use_bose) then
            if (w > 700.0_p_double) then
                dn_eq_dE = 0.0_p_double
            else
                dn_eq_dE = pref_ne * (E_erg*E_erg) / max(exp(w) - 1.0_p_double, 1.0d-99)
            end if
        else
            if (w > 700.0_p_double) then
                dn_eq_dE = 0.0_p_double
            else
                dn_eq_dE = pref_ne * (E_erg*E_erg) * exp(-w)
            end if
        end if

        if (dn_eq_dE <= 0.0_p_double) then
            p_abs = 1.0_p_double
        else
            gamma_abs = cs_scale * dndot_dE / dn_eq_dE
            if (gamma_abs < 0.0_p_double) gamma_abs = 0.0_p_double
            tau = gamma_abs * dt_phys
            if (tau > 50.0_p_double) then
                p_abs = 1.0_p_double
            else
                p_abs = 1.0_p_double - exp(-tau)
            end if
        end if

        call random_number(u)
        if (u < p_abs) then
            if (apply_feedback) then
                ! Deposit energy/momentum back to plasma (units: mec^2 and mec)
                ! Guard against out-of-range ix indices (axion species may carry guard-cell indices).
                select case ( p_x_dim )
                case ( 1 )
                    if (associated(electrons%axion_bremsstrahlung_info%depE%f1)) then
                        if ( this%ix(1,pp) >= lbound(electrons%axion_bremsstrahlung_info%depE%f1,2) .and. &
                             this%ix(1,pp) <= ubound(electrons%axion_bremsstrahlung_info%depE%f1,2) ) then
                            electrons%axion_bremsstrahlung_info%depE%f1(1, this%ix(1, pp)) = &
                                electrons%axion_bremsstrahlung_info%depE%f1(1, this%ix(1, pp)) + abs(this%q(pp)) * E_me
                            electrons%axion_bremsstrahlung_info%depP%f1(1, this%ix(1, pp)) = &
                                electrons%axion_bremsstrahlung_info%depP%f1(1, this%ix(1, pp)) + abs(this%q(pp)) * this%p(1, pp)
                            electrons%axion_bremsstrahlung_info%depP%f1(2, this%ix(1, pp)) = &
                                electrons%axion_bremsstrahlung_info%depP%f1(2, this%ix(1, pp)) + abs(this%q(pp)) * this%p(2, pp)
                            electrons%axion_bremsstrahlung_info%depP%f1(3, this%ix(1, pp)) = &
                                electrons%axion_bremsstrahlung_info%depP%f1(3, this%ix(1, pp)) + abs(this%q(pp)) * this%p(3, pp)
                        end if
                    end if

                case ( 2 )
                    if (associated(electrons%axion_bremsstrahlung_info%depE%f2)) then
                        if ( this%ix(1,pp) >= lbound(electrons%axion_bremsstrahlung_info%depE%f2,2) .and. &
                             this%ix(1,pp) <= ubound(electrons%axion_bremsstrahlung_info%depE%f2,2) .and. &
                             this%ix(2,pp) >= lbound(electrons%axion_bremsstrahlung_info%depE%f2,3) .and. &
                             this%ix(2,pp) <= ubound(electrons%axion_bremsstrahlung_info%depE%f2,3) ) then
                            electrons%axion_bremsstrahlung_info%depE%f2(1, this%ix(1, pp), this%ix(2, pp)) = &
                                electrons%axion_bremsstrahlung_info%depE%f2(1, this%ix(1, pp), this%ix(2, pp)) + abs(this%q(pp)) * E_me
                            electrons%axion_bremsstrahlung_info%depP%f2(1, this%ix(1, pp), this%ix(2, pp)) = &
                                electrons%axion_bremsstrahlung_info%depP%f2(1, this%ix(1, pp), this%ix(2, pp)) + abs(this%q(pp)) * this%p(1, pp)
                            electrons%axion_bremsstrahlung_info%depP%f2(2, this%ix(1, pp), this%ix(2, pp)) = &
                                electrons%axion_bremsstrahlung_info%depP%f2(2, this%ix(1, pp), this%ix(2, pp)) + abs(this%q(pp)) * this%p(2, pp)
                            electrons%axion_bremsstrahlung_info%depP%f2(3, this%ix(1, pp), this%ix(2, pp)) = &
                                electrons%axion_bremsstrahlung_info%depP%f2(3, this%ix(1, pp), this%ix(2, pp)) + abs(this%q(pp)) * this%p(3, pp)
                        end if
                    end if

                case ( 3 )
                    if (associated(electrons%axion_bremsstrahlung_info%depE%f3)) then
                        if ( this%ix(1,pp) >= lbound(electrons%axion_bremsstrahlung_info%depE%f3,2) .and. &
                             this%ix(1,pp) <= ubound(electrons%axion_bremsstrahlung_info%depE%f3,2) .and. &
                             this%ix(2,pp) >= lbound(electrons%axion_bremsstrahlung_info%depE%f3,3) .and. &
                             this%ix(2,pp) <= ubound(electrons%axion_bremsstrahlung_info%depE%f3,3) .and. &
                             this%ix(3,pp) >= lbound(electrons%axion_bremsstrahlung_info%depE%f3,4) .and. &
                             this%ix(3,pp) <= ubound(electrons%axion_bremsstrahlung_info%depE%f3,4) ) then
                            electrons%axion_bremsstrahlung_info%depE%f3(1, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp)) = &
                                electrons%axion_bremsstrahlung_info%depE%f3(1, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp)) + abs(this%q(pp)) * E_me
                            electrons%axion_bremsstrahlung_info%depP%f3(1, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp)) = &
                                electrons%axion_bremsstrahlung_info%depP%f3(1, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp)) + abs(this%q(pp)) * this%p(1, pp)
                            electrons%axion_bremsstrahlung_info%depP%f3(2, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp)) = &
                                electrons%axion_bremsstrahlung_info%depP%f3(2, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp)) + abs(this%q(pp)) * this%p(2, pp)
                            electrons%axion_bremsstrahlung_info%depP%f3(3, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp)) = &
                                electrons%axion_bremsstrahlung_info%depP%f3(3, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp)) + abs(this%q(pp)) * this%p(3, pp)
                        end if
                    end if
                end select

                electrons%axion_bremsstrahlung_info%fb_dE_step = electrons%axion_bremsstrahlung_info%fb_dE_step + abs(this%q(pp)) * E_me
                electrons%axion_bremsstrahlung_info%fb_dP_step(1) = electrons%axion_bremsstrahlung_info%fb_dP_step(1) + abs(this%q(pp)) * this%p(1, pp)
                electrons%axion_bremsstrahlung_info%fb_dP_step(2) = electrons%axion_bremsstrahlung_info%fb_dP_step(2) + abs(this%q(pp)) * this%p(2, pp)
                electrons%axion_bremsstrahlung_info%fb_dP_step(3) = electrons%axion_bremsstrahlung_info%fb_dP_step(3) + abs(this%q(pp)) * this%p(3, pp)
            end if

            call add_particle_index(local_idx, pp)
        end if

    end do

    if (local_idx%nidx > 0) then
        call remove_particles( this%t_species, local_idx )
    end if

    if (associated(local_idx%idx)) call freemem(local_idx%idx,"qed/os-axions.f03",1072)

contains

    !---------------------------------------------------------------------------------------------------
    subroutine add_particle_index( list, i )

        use m_species_define, only: p_spec_buf_block

        type( t_part_idx ), intent(inout) :: list
        integer, intent(in) :: i

        integer, dimension(:), pointer :: tmp_idx => null()

        ! grow buffer if necessary
        if ( list%nidx + 1 > list%buf_size ) then
            list%buf_size = list%buf_size + p_spec_buf_block
            call alloc(tmp_idx, (/ list%buf_size /),"qed/os-axions.f03",1089)
            if ( list%nidx > 0 ) call memcpy( tmp_idx, list%idx, list%nidx )
            if (associated(list%idx)) call freemem(list%idx,"qed/os-axions.f03",1091)
            list%idx => tmp_idx
        endif

        list%nidx = list%nidx + 1
        list%idx( list%nidx ) = i

    end subroutine add_particle_index
    !---------------------------------------------------------------------------------------------------

end subroutine axion_bremsstrahlung_absorb
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Apply per-cell energy/momentum deposits from axion bremsstrahlung emission/absorption to the plasma.
!
! This implements a conservative, cell-local affine mapping of particle momenta
! p -> s p + b
! such that each cell receives the requested total momentum change dP and (as much as
! possible) the requested kinetic energy change dE, without biasing the mean.
!
! Notes:
! - Uses a non-relativistic energy model E = sum w * p^2 / (2 mu), where mu = m_s/m_e.
! - Requires depE/depP to be accumulated during the timestep (Layer 0).
! - If evolve_Te is enabled, updates Te_cell from the post-map thermal variance.
!---------------------------------------------------------------------------------------------------
subroutine apply_axion_brems_feedback( electrons )

    use m_species_define, only : t_species
    use m_species_define_qed, only : t_species_qed

    implicit none

    class( t_species_qed ), intent(inout) :: electrons
    class( t_species ), pointer :: ions

    real(p_double) :: fe, one_minus_fe
    real(p_double) :: Zabs, mu_i
    logical :: do_update_Te

    if (.not. electrons%axion_bremsstrahlung_info%if_bremsstrahlung) return
    if (.not. electrons%axion_bremsstrahlung_info%apply_feedback) return

    ! nothing to do if feedback buffers are not allocated
    select case ( p_x_dim )
    case ( 1 )
        if (.not. associated(electrons%axion_bremsstrahlung_info%depE%f1)) return
    case ( 2 )
        if (.not. associated(electrons%axion_bremsstrahlung_info%depE%f2)) return
    case ( 3 )
        if (.not. associated(electrons%axion_bremsstrahlung_info%depE%f3)) return
    end select

    fe = electrons%axion_bremsstrahlung_info%fe_electrons
    if (fe < 0.0_p_double) fe = 0.0_p_double
    if (fe > 1.0_p_double) fe = 1.0_p_double
    one_minus_fe = 1.0_p_double - fe

    do_update_Te = electrons%axion_bremsstrahlung_info%evolve_Te

    ! electrons: Z=1, mu=1
    if (fe > 0.0_p_double) then
        call apply_to_species( electrons, 1.0_p_double, 1.0_p_double, fe, do_update_Te )
    else if (do_update_Te) then
        ! even if fe=0, we can still update Te_cell from the current distribution
        call apply_to_species( electrons, 1.0_p_double, 1.0_p_double, 0.0_p_double, .true. )
    end if

    ! ions: use pointer stored in the brems info
    ions => electrons%axion_bremsstrahlung_info%ion
    if (one_minus_fe > 0.0_p_double) then
        if (associated(ions)) then
            Zabs = abs(real(electrons%axion_bremsstrahlung_info%Z_ion, p_double))
            if (Zabs <= 0.0_p_double) Zabs = 1.0_p_double
            mu_i = max(abs(ions%rqm), 1.0d-30) / Zabs
            call apply_to_species( ions, Zabs, mu_i, one_minus_fe, .false. )
        end if
    end if

    ! clear deposits for next step
    call electrons%axion_bremsstrahlung_info%clear_feedback()

contains

    !-----------------------------------------------------------------------------------------------
    subroutine apply_to_species( sp, Zcharge_abs, mu, frac, update_Te )

        class( t_species ), intent(inout) :: sp
        real(p_double), intent(in) :: Zcharge_abs
        real(p_double), intent(in) :: mu
        real(p_double), intent(in) :: frac
        logical, intent(in) :: update_Te

        integer :: pp
        integer :: ix, iy, iz
        integer :: nx, ny, nz
        integer :: lx, ux, ly, uy, lz, uz
        real(p_double) :: wnum
        real(p_double) :: p1, p2, p3
        real(p_double) :: Wc, Px, Py, Pz, U2, meanx, meany, meanz
        real(p_double) :: dE, dPx, dPy, dPz, dx, dy, dz
        real(p_double) :: var, var_new, s, bx, by, bz
        real(p_double) :: dE_bulk, dE_rem
        real(p_double) :: Te_est, Te_old, Te_new, relax, floorT, ceilT
        real(p_double), parameter :: tinyW = 1.0d-99
        real(p_double), parameter :: tinyVar = 1.0d-99
        real(p_double), parameter :: mec2_keV = 5.11d2

        ! fast exit: if frac==0 and no Te update requested
        if (frac <= 0.0_p_double .and. .not. update_Te) return

        relax = electrons%axion_bremsstrahlung_info%Te_relax
        if (relax < 0.0_p_double) relax = 0.0_p_double
        if (relax > 1.0_p_double) relax = 1.0_p_double
        floorT = electrons%axion_bremsstrahlung_info%Te_floor_keV
        ceilT = electrons%axion_bremsstrahlung_info%Te_ceil_keV
        if (floorT <= 0.0_p_double) floorT = 0.0_p_double
        if (ceilT <= 0.0_p_double) ceilT = 1.0d99

        select case ( p_x_dim )

        case ( 1 )
            lx = lbound(electrons%axion_bremsstrahlung_info%depE%f1, 2)
            ux = ubound(electrons%axion_bremsstrahlung_info%depE%f1, 2)

            block

            ! cell moments
            real(p_double), allocatable :: W(:), Px1(:), Py1(:), Pz1(:), U2_1(:)
            real(p_double), allocatable :: s1(:), bx1(:), by1(:), bz1(:)
            allocate(W(lx:ux), Px1(lx:ux), Py1(lx:ux), Pz1(lx:ux), U2_1(lx:ux), s1(lx:ux), bx1(lx:ux), by1(lx:ux), bz1(lx:ux))
            W = 0.0_p_double; Px1 = 0.0_p_double; Py1 = 0.0_p_double; Pz1 = 0.0_p_double; U2_1 = 0.0_p_double
            s1 = 1.0_p_double; bx1 = 0.0_p_double; by1 = 0.0_p_double; bz1 = 0.0_p_double

            do pp = 1, sp%num_par
                ix = sp%ix(1, pp)
                if (ix < lx .or. ix > ux) cycle
                wnum = abs(real(sp%q(pp), p_double)) / max(Zcharge_abs, 1.0d-99)
                p1 = real(sp%p(1,pp), p_double)
                p2 = real(sp%p(2,pp), p_double)
                p3 = real(sp%p(3,pp), p_double)
                W(ix) = W(ix) + wnum
                Px1(ix) = Px1(ix) + wnum * p1
                Py1(ix) = Py1(ix) + wnum * p2
                Pz1(ix) = Pz1(ix) + wnum * p3
                U2_1(ix)= U2_1(ix)+ wnum * (p1*p1 + p2*p2 + p3*p3)
            end do

            do ix = lx, ux
                Wc = W(ix)
                if (Wc <= 0.0_p_double) cycle

                dE = frac * electrons%axion_bremsstrahlung_info%depE%f1(1, ix)
                dPx = frac * electrons%axion_bremsstrahlung_info%depP%f1(1, ix)
                dPy = frac * electrons%axion_bremsstrahlung_info%depP%f1(2, ix)
                dPz = frac * electrons%axion_bremsstrahlung_info%depP%f1(3, ix)

                Px = Px1(ix); Py = Py1(ix); Pz = Pz1(ix); U2 = U2_1(ix)

                meanx = Px / Wc
                meany = Py / Wc
                meanz = Pz / Wc

                var = U2 - (Px*Px + Py*Py + Pz*Pz) / max(Wc, tinyW)
                if (var < 0.0_p_double) var = 0.0_p_double

                ! momentum shift per particle
                dx = dPx / Wc
                dy = dPy / Wc
                dz = dPz / Wc

                ! energy change due to bulk shift
                dE_bulk = 0.5_p_double/mu * ( 2.0_p_double*(dx*Px + dy*Py + dz*Pz) + Wc*(dx*dx + dy*dy + dz*dz) )

                dE_rem = dE - dE_bulk
                var_new = var + 2.0_p_double*mu*dE_rem
                if (var_new < 0.0_p_double) var_new = 0.0_p_double

                if (var > tinyVar) then
                    s = sqrt( var_new / var )
                else
                    ! cannot change thermal variance if it's effectively zero
                    s = 1.0_p_double
                    var_new = var
                end if

                bx = (1.0_p_double - s)*meanx + dx
                by = (1.0_p_double - s)*meany + dy
                bz = (1.0_p_double - s)*meanz + dz

                s1(ix) = s
                bx1(ix) = bx
                by1(ix) = by
                bz1(ix) = bz

                if (update_Te) then
                    if (associated(electrons%axion_bremsstrahlung_info%Te_cell%f1)) then
                        Te_est = mec2_keV * (var_new / max(Wc, tinyW)) / 3.0_p_double
                        Te_old = electrons%axion_bremsstrahlung_info%Te_cell%f1(1, ix)
                        Te_new = (1.0_p_double - relax)*Te_old + relax*Te_est
                        if (Te_new < floorT) Te_new = floorT
                        if (Te_new > ceilT) Te_new = ceilT
                        electrons%axion_bremsstrahlung_info%Te_cell%f1(1, ix) = Te_new
                    end if
                end if
            end do

            do pp = 1, sp%num_par
                ix = sp%ix(1, pp)
                if (ix < lx .or. ix > ux) cycle
                sp%p(1,pp) = real( s1(ix) * real(sp%p(1,pp), p_double) + bx1(ix), p_k_part )
                sp%p(2,pp) = real( s1(ix) * real(sp%p(2,pp), p_double) + by1(ix), p_k_part )
                sp%p(3,pp) = real( s1(ix) * real(sp%p(3,pp), p_double) + bz1(ix), p_k_part )
            end do

            deallocate(W, Px1, Py1, Pz1, U2_1, s1, bx1, by1, bz1)

            end block

        case ( 2 )
            lx = lbound(electrons%axion_bremsstrahlung_info%depE%f2, 2)
            ux = ubound(electrons%axion_bremsstrahlung_info%depE%f2, 2)
            ly = lbound(electrons%axion_bremsstrahlung_info%depE%f2, 3)
            uy = ubound(electrons%axion_bremsstrahlung_info%depE%f2, 3)

            block

            real(p_double), allocatable :: W(:,:), Px2(:,:), Py2(:,:), Pz2(:,:), U2_2(:,:)
            real(p_double), allocatable :: s2(:,:), bx2(:,:), by2(:,:), bz2(:,:)
            allocate(W(lx:ux,ly:uy), Px2(lx:ux,ly:uy), Py2(lx:ux,ly:uy), Pz2(lx:ux,ly:uy), U2_2(lx:ux,ly:uy), s2(lx:ux,ly:uy), bx2(lx:ux,ly:uy), by2(lx:ux,ly:uy), bz2(lx:ux,ly:uy))
            W = 0.0_p_double; Px2 = 0.0_p_double; Py2 = 0.0_p_double; Pz2 = 0.0_p_double; U2_2 = 0.0_p_double
            s2 = 1.0_p_double; bx2 = 0.0_p_double; by2 = 0.0_p_double; bz2 = 0.0_p_double

            do pp = 1, sp%num_par
                ix = sp%ix(1, pp)
                iy = sp%ix(2, pp)
                if (ix < lx .or. ix > ux .or. iy < ly .or. iy > uy) cycle
                wnum = abs(real(sp%q(pp), p_double)) / max(Zcharge_abs, 1.0d-99)
                p1 = real(sp%p(1,pp), p_double)
                p2 = real(sp%p(2,pp), p_double)
                p3 = real(sp%p(3,pp), p_double)
                W(ix,iy) = W(ix,iy) + wnum
                Px2(ix,iy) = Px2(ix,iy) + wnum * p1
                Py2(ix,iy) = Py2(ix,iy) + wnum * p2
                Pz2(ix,iy) = Pz2(ix,iy) + wnum * p3
                U2_2(ix,iy)= U2_2(ix,iy)+ wnum * (p1*p1 + p2*p2 + p3*p3)
            end do

            ! Loop over the actual allocated cell bounds (including guard cells)
            do iy = ly, uy
                do ix = lx, ux
                    Wc = W(ix,iy)
                    if (Wc <= 0.0_p_double) cycle

                    dE = frac * electrons%axion_bremsstrahlung_info%depE%f2(1, ix, iy)
                    dPx = frac * electrons%axion_bremsstrahlung_info%depP%f2(1, ix, iy)
                    dPy = frac * electrons%axion_bremsstrahlung_info%depP%f2(2, ix, iy)
                    dPz = frac * electrons%axion_bremsstrahlung_info%depP%f2(3, ix, iy)

                    Px = Px2(ix,iy); Py = Py2(ix,iy); Pz = Pz2(ix,iy); U2 = U2_2(ix,iy)
                    meanx = Px / Wc
                    meany = Py / Wc
                    meanz = Pz / Wc

                    var = U2 - (Px*Px + Py*Py + Pz*Pz) / max(Wc, tinyW)
                    if (var < 0.0_p_double) var = 0.0_p_double

                    dx = dPx / Wc
                    dy = dPy / Wc
                    dz = dPz / Wc

                    dE_bulk = 0.5_p_double/mu * ( 2.0_p_double*(dx*Px + dy*Py + dz*Pz) + Wc*(dx*dx + dy*dy + dz*dz) )
                    dE_rem = dE - dE_bulk
                    var_new = var + 2.0_p_double*mu*dE_rem
                    if (var_new < 0.0_p_double) var_new = 0.0_p_double

                    if (var > tinyVar) then
                        s = sqrt( var_new / var )
                    else
                        s = 1.0_p_double
                        var_new = var
                    end if

                    bx = (1.0_p_double - s)*meanx + dx
                    by = (1.0_p_double - s)*meany + dy
                    bz = (1.0_p_double - s)*meanz + dz

                    s2(ix,iy) = s
                    bx2(ix,iy) = bx
                    by2(ix,iy) = by
                    bz2(ix,iy) = bz

                    if (update_Te) then
                        if (associated(electrons%axion_bremsstrahlung_info%Te_cell%f2)) then
                            Te_est = mec2_keV * (var_new / max(Wc, tinyW)) / 3.0_p_double
                            Te_old = electrons%axion_bremsstrahlung_info%Te_cell%f2(1, ix, iy)
                            Te_new = (1.0_p_double - relax)*Te_old + relax*Te_est
                            if (Te_new < floorT) Te_new = floorT
                            if (Te_new > ceilT) Te_new = ceilT
                            electrons%axion_bremsstrahlung_info%Te_cell%f2(1, ix, iy) = Te_new
                        end if
                    end if
                end do
            end do

            do pp = 1, sp%num_par
                ix = sp%ix(1, pp)
                iy = sp%ix(2, pp)
                if (ix < lx .or. ix > ux .or. iy < ly .or. iy > uy) cycle
                sp%p(1,pp) = real( s2(ix,iy) * real(sp%p(1,pp), p_double) + bx2(ix,iy), p_k_part )
                sp%p(2,pp) = real( s2(ix,iy) * real(sp%p(2,pp), p_double) + by2(ix,iy), p_k_part )
                sp%p(3,pp) = real( s2(ix,iy) * real(sp%p(3,pp), p_double) + bz2(ix,iy), p_k_part )
            end do

            deallocate(W, Px2, Py2, Pz2, U2_2, s2, bx2, by2, bz2)

            end block

        case ( 3 )
            lx = lbound(electrons%axion_bremsstrahlung_info%depE%f3, 2)
            ux = ubound(electrons%axion_bremsstrahlung_info%depE%f3, 2)
            ly = lbound(electrons%axion_bremsstrahlung_info%depE%f3, 3)
            uy = ubound(electrons%axion_bremsstrahlung_info%depE%f3, 3)
            lz = lbound(electrons%axion_bremsstrahlung_info%depE%f3, 4)
            uz = ubound(electrons%axion_bremsstrahlung_info%depE%f3, 4)

            block

            real(p_double), allocatable :: W(:,:,:), Px3(:,:,:), Py3(:,:,:), Pz3(:,:,:), U2_3(:,:,:)
            real(p_double), allocatable :: s3(:,:,:), bx3(:,:,:), by3(:,:,:), bz3(:,:,:)
            allocate(W(lx:ux,ly:uy,lz:uz), Px3(lx:ux,ly:uy,lz:uz), Py3(lx:ux,ly:uy,lz:uz), Pz3(lx:ux,ly:uy,lz:uz), U2_3(lx:ux,ly:uy,lz:uz), s3(lx:ux,ly:uy,lz:uz), bx3(lx:ux,ly:uy,lz:uz), by3(lx:ux,ly:uy,lz:uz), bz3(lx:ux,ly:uy,lz:uz))
            W = 0.0_p_double; Px3 = 0.0_p_double; Py3 = 0.0_p_double; Pz3 = 0.0_p_double; U2_3 = 0.0_p_double
            s3 = 1.0_p_double; bx3 = 0.0_p_double; by3 = 0.0_p_double; bz3 = 0.0_p_double

            do pp = 1, sp%num_par
                ix = sp%ix(1, pp)
                iy = sp%ix(2, pp)
                iz = sp%ix(3, pp)
                if (ix < lx .or. ix > ux .or. iy < ly .or. iy > uy .or. iz < lz .or. iz > uz) cycle
                wnum = abs(real(sp%q(pp), p_double)) / max(Zcharge_abs, 1.0d-99)
                p1 = real(sp%p(1,pp), p_double)
                p2 = real(sp%p(2,pp), p_double)
                p3 = real(sp%p(3,pp), p_double)
                W(ix,iy,iz) = W(ix,iy,iz) + wnum
                Px3(ix,iy,iz) = Px3(ix,iy,iz) + wnum * p1
                Py3(ix,iy,iz) = Py3(ix,iy,iz) + wnum * p2
                Pz3(ix,iy,iz) = Pz3(ix,iy,iz) + wnum * p3
                U2_3(ix,iy,iz)= U2_3(ix,iy,iz)+ wnum * (p1*p1 + p2*p2 + p3*p3)
            end do

            ! Loop over the actual allocated cell bounds (including guard cells)
            do iz = lz, uz
                do iy = ly, uy
                    do ix = lx, ux
                        Wc = W(ix,iy,iz)
                        if (Wc <= 0.0_p_double) cycle

                        dE = frac * electrons%axion_bremsstrahlung_info%depE%f3(1, ix, iy, iz)
                        dPx = frac * electrons%axion_bremsstrahlung_info%depP%f3(1, ix, iy, iz)
                        dPy = frac * electrons%axion_bremsstrahlung_info%depP%f3(2, ix, iy, iz)
                        dPz = frac * electrons%axion_bremsstrahlung_info%depP%f3(3, ix, iy, iz)

                        Px = Px3(ix,iy,iz); Py = Py3(ix,iy,iz); Pz = Pz3(ix,iy,iz); U2 = U2_3(ix,iy,iz)
                        meanx = Px / Wc
                        meany = Py / Wc
                        meanz = Pz / Wc

                        var = U2 - (Px*Px + Py*Py + Pz*Pz) / max(Wc, tinyW)
                        if (var < 0.0_p_double) var = 0.0_p_double

                        dx = dPx / Wc
                        dy = dPy / Wc
                        dz = dPz / Wc

                        dE_bulk = 0.5_p_double/mu * ( 2.0_p_double*(dx*Px + dy*Py + dz*Pz) + Wc*(dx*dx + dy*dy + dz*dz) )
                        dE_rem = dE - dE_bulk
                        var_new = var + 2.0_p_double*mu*dE_rem
                        if (var_new < 0.0_p_double) var_new = 0.0_p_double

                        if (var > tinyVar) then
                            s = sqrt( var_new / var )
                        else
                            s = 1.0_p_double
                            var_new = var
                        end if

                        bx = (1.0_p_double - s)*meanx + dx
                        by = (1.0_p_double - s)*meany + dy
                        bz = (1.0_p_double - s)*meanz + dz

                        s3(ix,iy,iz) = s
                        bx3(ix,iy,iz) = bx
                        by3(ix,iy,iz) = by
                        bz3(ix,iy,iz) = bz

                        if (update_Te) then
                            if (associated(electrons%axion_bremsstrahlung_info%Te_cell%f3)) then
                                Te_est = mec2_keV * (var_new / max(Wc, tinyW)) / 3.0_p_double
                                Te_old = electrons%axion_bremsstrahlung_info%Te_cell%f3(1, ix, iy, iz)
                                Te_new = (1.0_p_double - relax)*Te_old + relax*Te_est
                                if (Te_new < floorT) Te_new = floorT
                                if (Te_new > ceilT) Te_new = ceilT
                                electrons%axion_bremsstrahlung_info%Te_cell%f3(1, ix, iy, iz) = Te_new
                            end if
                        end if
                    end do
                end do
            end do

            do pp = 1, sp%num_par
                ix = sp%ix(1, pp)
                iy = sp%ix(2, pp)
                iz = sp%ix(3, pp)
                if (ix < lx .or. ix > ux .or. iy < ly .or. iy > uy .or. iz < lz .or. iz > uz) cycle
                sp%p(1,pp) = real( s3(ix,iy,iz) * real(sp%p(1,pp), p_double) + bx3(ix,iy,iz), p_k_part )
                sp%p(2,pp) = real( s3(ix,iy,iz) * real(sp%p(2,pp), p_double) + by3(ix,iy,iz), p_k_part )
                sp%p(3,pp) = real( s3(ix,iy,iz) * real(sp%p(3,pp), p_double) + bz3(ix,iy,iz), p_k_part )
            end do

            deallocate(W, Px3, Py3, Pz3, U2_3, s3, bx3, by3, bz3)

            end block
        end select

    end subroutine apply_to_species

end subroutine apply_axion_brems_feedback
!---------------------------------------------------------------------------------------------------

subroutine axion_compton_absorb( this, electrons, gdt, i0, i1 )

    use m_system
    use m_parameters
    use m_axions_define, only : t_axions
    use m_species_define_qed, only : t_species_qed
    use m_species_comm, only : t_part_idx, remove_particles

    implicit none

    class( t_axions ), intent(inout) :: this
    class( t_species_qed ), intent(inout) :: electrons
    real(p_double), intent(in) :: gdt
    integer, intent(in) :: i0, i1

    integer :: pp, i1_use
    real(p_double) :: omega_p, dt_phys, n0_cm3

    ! local plasma / model parameters
    real(p_double) :: Te_keV, g_ae, cs_scale
    real(p_double) :: X_H, A_ion
    logical :: use_bose, apply_feedback, evolve_Te

    ! constants (cgs)
    real(p_double), parameter :: kB_eV_per_K = 8.617333262145d-5
    real(p_double), parameter :: erg_per_eV = 1.602176634d-12
    real(p_double), parameter :: amu_g = 1.66053906660d-24
    real(p_double), parameter :: hbar_erg_s = 1.054571817d-27
    real(p_double), parameter :: c_cms = 2.99792458d10
    real(p_double), parameter :: hbarc_erg_cm = hbar_erg_s * c_cms

    ! spectrum proxy parameters for w = E/T (Gamma distribution): k=4.9, theta=1
    real(p_double), parameter :: w_shape = 4.9d0
    real(p_double), parameter :: w_theta = 1.0d0
    real(p_double), parameter :: w_mean = 4.9d0
    real(p_double), parameter :: w_norm = 4.8385412738965787d-2 ! 1/Gamma(4.9)

    real(p_double) :: pi
    real(p_double) :: T7, rho_gcc
    real(p_double) :: alpha_ae
    real(p_double) :: eps_c, Q_a
    real(p_double) :: kT_erg, meanE_erg, ndot_cm3_s

    ! particle energy / spectrum
    real(p_double) :: p2, E_me, w, pdf_w
    real(p_double) :: E_erg
    real(p_double) :: dn_eq_dE, dndot_dE, gamma_abs, tau, p_abs, u
    real(p_double) :: pref_ne

    ! local ion density (normalized)
    real(p_double) :: ni

    type( t_part_idx ) :: local_idx

    ! Defensive initialization
    local_idx%nidx = 0
    local_idx%buf_size = 0
    nullify(local_idx%idx)

    if (.not. electrons%axion_compton_info%if_compton) return
    if (.not. electrons%axion_compton_info%if_inverse) return

    Te_keV = electrons%axion_compton_info%Te_keV
    g_ae = electrons%axion_compton_info%g_ae
    cs_scale = electrons%axion_compton_info%cs_scale
    X_H = electrons%axion_compton_info%X_H
    A_ion = electrons%axion_compton_info%A_ion
    use_bose = electrons%axion_compton_info%use_bose

    apply_feedback = electrons%axion_compton_info%apply_feedback
    evolve_Te = electrons%axion_compton_info%evolve_Te

    if (apply_feedback) call electrons%axion_compton_info%ensure_feedback()

    if (Te_keV <= 0.0_p_double) return
    if (g_ae <= 0.0_p_double) return
    if (cs_scale <= 0.0_p_double) return

    omega_p = electrons%omega_p0
    if (omega_p <= 0.0_p_double) return
    dt_phys = gdt / omega_p
    if (dt_phys <= 0.0_p_double) return

    n0_cm3 = electrons%n0
    pi = acos(-1.0_p_double)

    ! alpha_ae = g_ae^2/(4*pi)
    alpha_ae = (g_ae*g_ae) / (4.0_p_double*pi)

    ! equilibrium dn/dE prefactor: E^2 / [ 2 pi^2 (hbar c)^3 ]
    pref_ne = 1.0_p_double / (2.0_p_double * pi*pi * (hbarc_erg_cm**3))

    ! init delete list
    local_idx%start = 1
    local_idx%nidx = 0
    local_idx%buf_size = 0
    local_idx%idx => null()

    ! Defensive clamp: this routine may be called with stale i1 if other operators
    ! removed axions earlier in the same timestep.
    i1_use = min(i1, this%num_par)
    if (i0 > i1_use) return

    do pp = i0, i1_use

        ! local ion density at axion location (normalized units)
        select case ( p_x_dim )
        case ( 1 )
            ni = electrons % axion_compton_info % rho_ion % f1(1, this%ix(1, pp))
        case ( 2 )
            ni = electrons % axion_compton_info % rho_ion % f2(1, this%ix(1, pp), this%ix(2, pp))
        case ( 3 )
            ni = electrons % axion_compton_info % rho_ion % f3(1, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp))
        end select
        if (ni < 0.0_p_double) ni = 0.0_p_double
        if (ni <= 0.0_p_double) cycle

        ! local temperature (keV) optionally from evolving field
        if (evolve_Te) then
            select case ( p_x_dim )
            case ( 1 )
                if (associated(electrons%axion_compton_info%Te_cell%f1)) then
                    Te_keV = electrons%axion_compton_info%Te_cell%f1(1, this%ix(1, pp))
                end if
            case ( 2 )
                if (associated(electrons%axion_compton_info%Te_cell%f2)) then
                    Te_keV = electrons%axion_compton_info%Te_cell%f2(1, this%ix(1, pp), this%ix(2, pp))
                end if
            case ( 3 )
                if (associated(electrons%axion_compton_info%Te_cell%f3)) then
                    Te_keV = electrons%axion_compton_info%Te_cell%f3(1, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp))
                end if
            end select
            if (Te_keV <= 0.0_p_double) cycle
        end if

        ! convert Te_keV -> T7
        T7 = (Te_keV * 1.0d3) / (kB_eV_per_K * 1.0d7)
        kT_erg = (Te_keV * 1.0d3) * erg_per_eV
        meanE_erg = w_mean * kT_erg

        ! axion energy (dimensionless E/(mec^2) = |p|)
        p2 = this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2
        if (p2 <= 0.0_p_double) cycle
        E_me = sqrt(p2)

        ! w = E/T with E in keV, T in keV; E = E_me * 511 keV
        w = (E_me * 5.11d2) / Te_keV
        if (w <= 0.0_p_double) cycle

        ! proxy emission spectrum pdf in w
        pdf_w = w_norm * (w**(w_shape - 1.0_p_double)) * exp( - w / w_theta )

        ! local mass density (g/cm^3)
        rho_gcc = (ni * n0_cm3) * (A_ion * amu_g)

        ! Raffelt Compton emissivity per mass
        eps_c = 2.67d-2 * (1.0_p_double + X_H) * (T7**6) * (alpha_ae / 1.60d-23)

        Q_a = rho_gcc * eps_c
        ndot_cm3_s = Q_a / max(meanE_erg, 1.0d-99)

        ! differential number emissivity per energy
        dndot_dE = ndot_cm3_s * pdf_w / max(kT_erg, 1.0d-99)

        ! equilibrium axion number density per energy
        E_erg = w * kT_erg
        if (use_bose) then
            if (w > 700.0_p_double) then
                dn_eq_dE = 0.0_p_double
            else
                dn_eq_dE = pref_ne * (E_erg*E_erg) / max(exp(w) - 1.0_p_double, 1.0d-99)
            end if
        else
            if (w > 700.0_p_double) then
                dn_eq_dE = 0.0_p_double
            else
                dn_eq_dE = pref_ne * (E_erg*E_erg) * exp(-w)
            end if
        end if

        if (dn_eq_dE <= 0.0_p_double) then
            p_abs = 1.0_p_double
        else
            gamma_abs = cs_scale * dndot_dE / dn_eq_dE
            if (gamma_abs < 0.0_p_double) gamma_abs = 0.0_p_double
            tau = gamma_abs * dt_phys
            if (tau > 50.0_p_double) then
                p_abs = 1.0_p_double
            else
                p_abs = 1.0_p_double - exp(-tau)
            end if
        end if

        call random_number(u)
        if (u < p_abs) then

            if (apply_feedback) then
                ! Deposit energy/momentum back to plasma (units: mec^2 and mec)
                select case ( p_x_dim )
                case ( 1 )
                    if (associated(electrons%axion_compton_info%depE%f1)) then
                        if ( this%ix(1,pp) >= lbound(electrons%axion_compton_info%depE%f1,2) .and. &
                             this%ix(1,pp) <= ubound(electrons%axion_compton_info%depE%f1,2) ) then
                            electrons%axion_compton_info%depE%f1(1, this%ix(1, pp)) = &
                                electrons%axion_compton_info%depE%f1(1, this%ix(1, pp)) + abs(this%q(pp)) * E_me
                            electrons%axion_compton_info%depP%f1(1, this%ix(1, pp)) = &
                                electrons%axion_compton_info%depP%f1(1, this%ix(1, pp)) + abs(this%q(pp)) * this%p(1, pp)
                            electrons%axion_compton_info%depP%f1(2, this%ix(1, pp)) = &
                                electrons%axion_compton_info%depP%f1(2, this%ix(1, pp)) + abs(this%q(pp)) * this%p(2, pp)
                            electrons%axion_compton_info%depP%f1(3, this%ix(1, pp)) = &
                                electrons%axion_compton_info%depP%f1(3, this%ix(1, pp)) + abs(this%q(pp)) * this%p(3, pp)
                        end if
                    end if

                case ( 2 )
                    if (associated(electrons%axion_compton_info%depE%f2)) then
                        if ( this%ix(1,pp) >= lbound(electrons%axion_compton_info%depE%f2,2) .and. &
                             this%ix(1,pp) <= ubound(electrons%axion_compton_info%depE%f2,2) .and. &
                             this%ix(2,pp) >= lbound(electrons%axion_compton_info%depE%f2,3) .and. &
                             this%ix(2,pp) <= ubound(electrons%axion_compton_info%depE%f2,3) ) then
                            electrons%axion_compton_info%depE%f2(1, this%ix(1, pp), this%ix(2, pp)) = &
                                electrons%axion_compton_info%depE%f2(1, this%ix(1, pp), this%ix(2, pp)) + abs(this%q(pp)) * E_me
                            electrons%axion_compton_info%depP%f2(1, this%ix(1, pp), this%ix(2, pp)) = &
                                electrons%axion_compton_info%depP%f2(1, this%ix(1, pp), this%ix(2, pp)) + abs(this%q(pp)) * this%p(1, pp)
                            electrons%axion_compton_info%depP%f2(2, this%ix(1, pp), this%ix(2, pp)) = &
                                electrons%axion_compton_info%depP%f2(2, this%ix(1, pp), this%ix(2, pp)) + abs(this%q(pp)) * this%p(2, pp)
                            electrons%axion_compton_info%depP%f2(3, this%ix(1, pp), this%ix(2, pp)) = &
                                electrons%axion_compton_info%depP%f2(3, this%ix(1, pp), this%ix(2, pp)) + abs(this%q(pp)) * this%p(3, pp)
                        end if
                    end if

                case ( 3 )
                    if (associated(electrons%axion_compton_info%depE%f3)) then
                        if ( this%ix(1,pp) >= lbound(electrons%axion_compton_info%depE%f3,2) .and. &
                             this%ix(1,pp) <= ubound(electrons%axion_compton_info%depE%f3,2) .and. &
                             this%ix(2,pp) >= lbound(electrons%axion_compton_info%depE%f3,3) .and. &
                             this%ix(2,pp) <= ubound(electrons%axion_compton_info%depE%f3,3) .and. &
                             this%ix(3,pp) >= lbound(electrons%axion_compton_info%depE%f3,4) .and. &
                             this%ix(3,pp) <= ubound(electrons%axion_compton_info%depE%f3,4) ) then
                            electrons%axion_compton_info%depE%f3(1, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp)) = &
                                electrons%axion_compton_info%depE%f3(1, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp)) + abs(this%q(pp)) * E_me
                            electrons%axion_compton_info%depP%f3(1, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp)) = &
                                electrons%axion_compton_info%depP%f3(1, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp)) + abs(this%q(pp)) * this%p(1, pp)
                            electrons%axion_compton_info%depP%f3(2, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp)) = &
                                electrons%axion_compton_info%depP%f3(2, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp)) + abs(this%q(pp)) * this%p(2, pp)
                            electrons%axion_compton_info%depP%f3(3, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp)) = &
                                electrons%axion_compton_info%depP%f3(3, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp)) + abs(this%q(pp)) * this%p(3, pp)
                        end if
                    end if
                end select

                electrons%axion_compton_info%fb_dE_step = electrons%axion_compton_info%fb_dE_step + abs(this%q(pp)) * E_me
                electrons%axion_compton_info%fb_dP_step(1) = electrons%axion_compton_info%fb_dP_step(1) + abs(this%q(pp)) * this%p(1, pp)
                electrons%axion_compton_info%fb_dP_step(2) = electrons%axion_compton_info%fb_dP_step(2) + abs(this%q(pp)) * this%p(2, pp)
                electrons%axion_compton_info%fb_dP_step(3) = electrons%axion_compton_info%fb_dP_step(3) + abs(this%q(pp)) * this%p(3, pp)
            end if

            call add_particle_index(local_idx, pp)
        end if

    end do

    if (local_idx%nidx > 0) then
        call remove_particles( this%t_species, local_idx )
    end if

    if (associated(local_idx%idx)) call freemem(local_idx%idx,"qed/os-axions.f03",1787)

contains

    subroutine add_particle_index( list, i )

        use m_species_define, only: p_spec_buf_block

        type( t_part_idx ), intent(inout) :: list
        integer, intent(in) :: i

        integer, dimension(:), pointer :: tmp_idx => null()

        if ( list%nidx + 1 > list%buf_size ) then
            list%buf_size = list%buf_size + p_spec_buf_block
            call alloc(tmp_idx, (/ list%buf_size /),"qed/os-axions.f03",1802)
            if ( list%nidx > 0 ) call memcpy( tmp_idx, list%idx, list%nidx )
            if (associated(list%idx)) call freemem(list%idx,"qed/os-axions.f03",1804)
            list%idx => tmp_idx
        endif

        list%nidx = list%nidx + 1
        list%idx( list%nidx ) = i

    end subroutine add_particle_index

end subroutine axion_compton_absorb
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Apply per-cell energy/momentum deposits from axion Compton emission/absorption to the plasma.
! Identical operator to apply_axion_brems_feedback but acting on axion_compton_info buffers.
!---------------------------------------------------------------------------------------------------
subroutine apply_axion_compton_feedback( electrons )

    use m_species_define, only : t_species
    use m_species_define_qed, only : t_species_qed

    implicit none

    class( t_species_qed ), intent(inout) :: electrons
    class( t_species ), pointer :: ions

    real(p_double) :: fe, one_minus_fe
    real(p_double) :: Zabs, mu_i
    logical :: do_update_Te

    if (.not. electrons%axion_compton_info%if_compton) return
    if (.not. electrons%axion_compton_info%apply_feedback) return

    ! nothing to do if feedback buffers are not allocated
    select case ( p_x_dim )
    case ( 1 )
        if (.not. associated(electrons%axion_compton_info%depE%f1)) return
    case ( 2 )
        if (.not. associated(electrons%axion_compton_info%depE%f2)) return
    case ( 3 )
        if (.not. associated(electrons%axion_compton_info%depE%f3)) return
    end select

    fe = electrons%axion_compton_info%fe_electrons
    if (fe < 0.0_p_double) fe = 0.0_p_double
    if (fe > 1.0_p_double) fe = 1.0_p_double
    one_minus_fe = 1.0_p_double - fe

    do_update_Te = electrons%axion_compton_info%evolve_Te

    ! electrons: Z=1, mu=1
    if (fe > 0.0_p_double) then
        call apply_to_species( electrons, 1.0_p_double, 1.0_p_double, fe, do_update_Te )
    else if (do_update_Te) then
        ! even if fe=0, we can still update Te_cell from the current distribution
        call apply_to_species( electrons, 1.0_p_double, 1.0_p_double, 0.0_p_double, .true. )
    end if

    ! ions: use pointer stored in the brems info
    ions => electrons%axion_compton_info%ion
    if (one_minus_fe > 0.0_p_double) then
        if (associated(ions)) then
            Zabs = abs(real(electrons%axion_compton_info%Z_ion, p_double))
            if (Zabs <= 0.0_p_double) Zabs = 1.0_p_double
            mu_i = max(abs(ions%rqm), 1.0d-30) / Zabs
            call apply_to_species( ions, Zabs, mu_i, one_minus_fe, .false. )
        end if
    end if

    ! clear deposits for next step
    call electrons%axion_compton_info%clear_feedback()

contains

    !-----------------------------------------------------------------------------------------------
    subroutine apply_to_species( sp, Zcharge_abs, mu, frac, update_Te )

        class( t_species ), intent(inout) :: sp
        real(p_double), intent(in) :: Zcharge_abs
        real(p_double), intent(in) :: mu
        real(p_double), intent(in) :: frac
        logical, intent(in) :: update_Te

        integer :: pp
        integer :: ix, iy, iz
        integer :: nx, ny, nz
        integer :: lx, ux, ly, uy, lz, uz
        real(p_double) :: wnum
        real(p_double) :: p1, p2, p3
        real(p_double) :: Wc, Px, Py, Pz, U2, meanx, meany, meanz
        real(p_double) :: dE, dPx, dPy, dPz, dx, dy, dz
        real(p_double) :: var, var_new, s, bx, by, bz
        real(p_double) :: dE_bulk, dE_rem
        real(p_double) :: Te_est, Te_old, Te_new, relax, floorT, ceilT
        real(p_double), parameter :: tinyW = 1.0d-99
        real(p_double), parameter :: tinyVar = 1.0d-99
        real(p_double), parameter :: mec2_keV = 5.11d2

        ! fast exit: if frac==0 and no Te update requested
        if (frac <= 0.0_p_double .and. .not. update_Te) return

        relax = electrons%axion_compton_info%Te_relax
        if (relax < 0.0_p_double) relax = 0.0_p_double
        if (relax > 1.0_p_double) relax = 1.0_p_double
        floorT = electrons%axion_compton_info%Te_floor_keV
        ceilT = electrons%axion_compton_info%Te_ceil_keV
        if (floorT <= 0.0_p_double) floorT = 0.0_p_double
        if (ceilT <= 0.0_p_double) ceilT = 1.0d99

        select case ( p_x_dim )

        case ( 1 )
            lx = lbound(electrons%axion_compton_info%depE%f1, 2)
            ux = ubound(electrons%axion_compton_info%depE%f1, 2)

            block

            ! cell moments
            real(p_double), allocatable :: W(:), Px1(:), Py1(:), Pz1(:), U2_1(:)
            real(p_double), allocatable :: s1(:), bx1(:), by1(:), bz1(:)
            allocate(W(lx:ux), Px1(lx:ux), Py1(lx:ux), Pz1(lx:ux), U2_1(lx:ux), s1(lx:ux), bx1(lx:ux), by1(lx:ux), bz1(lx:ux))
            W = 0.0_p_double; Px1 = 0.0_p_double; Py1 = 0.0_p_double; Pz1 = 0.0_p_double; U2_1 = 0.0_p_double
            s1 = 1.0_p_double; bx1 = 0.0_p_double; by1 = 0.0_p_double; bz1 = 0.0_p_double

            do pp = 1, sp%num_par
                ix = sp%ix(1, pp)
                if (ix < lx .or. ix > ux) cycle
                wnum = abs(real(sp%q(pp), p_double)) / max(Zcharge_abs, 1.0d-99)
                p1 = real(sp%p(1,pp), p_double)
                p2 = real(sp%p(2,pp), p_double)
                p3 = real(sp%p(3,pp), p_double)
                W(ix) = W(ix) + wnum
                Px1(ix) = Px1(ix) + wnum * p1
                Py1(ix) = Py1(ix) + wnum * p2
                Pz1(ix) = Pz1(ix) + wnum * p3
                U2_1(ix)= U2_1(ix)+ wnum * (p1*p1 + p2*p2 + p3*p3)
            end do

            do ix = lx, ux
                Wc = W(ix)
                if (Wc <= 0.0_p_double) cycle

                dE = frac * electrons%axion_compton_info%depE%f1(1, ix)
                dPx = frac * electrons%axion_compton_info%depP%f1(1, ix)
                dPy = frac * electrons%axion_compton_info%depP%f1(2, ix)
                dPz = frac * electrons%axion_compton_info%depP%f1(3, ix)

                Px = Px1(ix); Py = Py1(ix); Pz = Pz1(ix); U2 = U2_1(ix)

                meanx = Px / Wc
                meany = Py / Wc
                meanz = Pz / Wc

                var = U2 - (Px*Px + Py*Py + Pz*Pz) / max(Wc, tinyW)
                if (var < 0.0_p_double) var = 0.0_p_double

                ! momentum shift per particle
                dx = dPx / Wc
                dy = dPy / Wc
                dz = dPz / Wc

                ! energy change due to bulk shift
                dE_bulk = 0.5_p_double/mu * ( 2.0_p_double*(dx*Px + dy*Py + dz*Pz) + Wc*(dx*dx + dy*dy + dz*dz) )

                dE_rem = dE - dE_bulk
                var_new = var + 2.0_p_double*mu*dE_rem
                if (var_new < 0.0_p_double) var_new = 0.0_p_double

                if (var > tinyVar) then
                    s = sqrt( var_new / var )
                else
                    ! cannot change thermal variance if it's effectively zero
                    s = 1.0_p_double
                    var_new = var
                end if

                bx = (1.0_p_double - s)*meanx + dx
                by = (1.0_p_double - s)*meany + dy
                bz = (1.0_p_double - s)*meanz + dz

                s1(ix) = s
                bx1(ix) = bx
                by1(ix) = by
                bz1(ix) = bz

                if (update_Te) then
                    if (associated(electrons%axion_compton_info%Te_cell%f1)) then
                        Te_est = mec2_keV * (var_new / max(Wc, tinyW)) / 3.0_p_double
                        Te_old = electrons%axion_compton_info%Te_cell%f1(1, ix)
                        Te_new = (1.0_p_double - relax)*Te_old + relax*Te_est
                        if (Te_new < floorT) Te_new = floorT
                        if (Te_new > ceilT) Te_new = ceilT
                        electrons%axion_compton_info%Te_cell%f1(1, ix) = Te_new
                    end if
                end if
            end do

            do pp = 1, sp%num_par
                ix = sp%ix(1, pp)
                if (ix < lx .or. ix > ux) cycle
                sp%p(1,pp) = real( s1(ix) * real(sp%p(1,pp), p_double) + bx1(ix), p_k_part )
                sp%p(2,pp) = real( s1(ix) * real(sp%p(2,pp), p_double) + by1(ix), p_k_part )
                sp%p(3,pp) = real( s1(ix) * real(sp%p(3,pp), p_double) + bz1(ix), p_k_part )
            end do

            deallocate(W, Px1, Py1, Pz1, U2_1, s1, bx1, by1, bz1)

            end block

        case ( 2 )
            lx = lbound(electrons%axion_compton_info%depE%f2, 2)
            ux = ubound(electrons%axion_compton_info%depE%f2, 2)
            ly = lbound(electrons%axion_compton_info%depE%f2, 3)
            uy = ubound(electrons%axion_compton_info%depE%f2, 3)

            block

            real(p_double), allocatable :: W(:,:), Px2(:,:), Py2(:,:), Pz2(:,:), U2_2(:,:)
            real(p_double), allocatable :: s2(:,:), bx2(:,:), by2(:,:), bz2(:,:)
            allocate(W(lx:ux,ly:uy), Px2(lx:ux,ly:uy), Py2(lx:ux,ly:uy), Pz2(lx:ux,ly:uy), U2_2(lx:ux,ly:uy), s2(lx:ux,ly:uy), bx2(lx:ux,ly:uy), by2(lx:ux,ly:uy), bz2(lx:ux,ly:uy))
            W = 0.0_p_double; Px2 = 0.0_p_double; Py2 = 0.0_p_double; Pz2 = 0.0_p_double; U2_2 = 0.0_p_double
            s2 = 1.0_p_double; bx2 = 0.0_p_double; by2 = 0.0_p_double; bz2 = 0.0_p_double

            do pp = 1, sp%num_par
                ix = sp%ix(1, pp)
                iy = sp%ix(2, pp)
                if (ix < lx .or. ix > ux .or. iy < ly .or. iy > uy) cycle
                wnum = abs(real(sp%q(pp), p_double)) / max(Zcharge_abs, 1.0d-99)
                p1 = real(sp%p(1,pp), p_double)
                p2 = real(sp%p(2,pp), p_double)
                p3 = real(sp%p(3,pp), p_double)
                W(ix,iy) = W(ix,iy) + wnum
                Px2(ix,iy) = Px2(ix,iy) + wnum * p1
                Py2(ix,iy) = Py2(ix,iy) + wnum * p2
                Pz2(ix,iy) = Pz2(ix,iy) + wnum * p3
                U2_2(ix,iy)= U2_2(ix,iy)+ wnum * (p1*p1 + p2*p2 + p3*p3)
            end do

            ! Loop over the actual allocated cell bounds (including guard cells)
            do iy = ly, uy
                do ix = lx, ux
                    Wc = W(ix,iy)
                    if (Wc <= 0.0_p_double) cycle

                    dE = frac * electrons%axion_compton_info%depE%f2(1, ix, iy)
                    dPx = frac * electrons%axion_compton_info%depP%f2(1, ix, iy)
                    dPy = frac * electrons%axion_compton_info%depP%f2(2, ix, iy)
                    dPz = frac * electrons%axion_compton_info%depP%f2(3, ix, iy)

                    Px = Px2(ix,iy); Py = Py2(ix,iy); Pz = Pz2(ix,iy); U2 = U2_2(ix,iy)
                    meanx = Px / Wc
                    meany = Py / Wc
                    meanz = Pz / Wc

                    var = U2 - (Px*Px + Py*Py + Pz*Pz) / max(Wc, tinyW)
                    if (var < 0.0_p_double) var = 0.0_p_double

                    dx = dPx / Wc
                    dy = dPy / Wc
                    dz = dPz / Wc

                    dE_bulk = 0.5_p_double/mu * ( 2.0_p_double*(dx*Px + dy*Py + dz*Pz) + Wc*(dx*dx + dy*dy + dz*dz) )
                    dE_rem = dE - dE_bulk
                    var_new = var + 2.0_p_double*mu*dE_rem
                    if (var_new < 0.0_p_double) var_new = 0.0_p_double

                    if (var > tinyVar) then
                        s = sqrt( var_new / var )
                    else
                        s = 1.0_p_double
                        var_new = var
                    end if

                    bx = (1.0_p_double - s)*meanx + dx
                    by = (1.0_p_double - s)*meany + dy
                    bz = (1.0_p_double - s)*meanz + dz

                    s2(ix,iy) = s
                    bx2(ix,iy) = bx
                    by2(ix,iy) = by
                    bz2(ix,iy) = bz

                    if (update_Te) then
                        if (associated(electrons%axion_compton_info%Te_cell%f2)) then
                            Te_est = mec2_keV * (var_new / max(Wc, tinyW)) / 3.0_p_double
                            Te_old = electrons%axion_compton_info%Te_cell%f2(1, ix, iy)
                            Te_new = (1.0_p_double - relax)*Te_old + relax*Te_est
                            if (Te_new < floorT) Te_new = floorT
                            if (Te_new > ceilT) Te_new = ceilT
                            electrons%axion_compton_info%Te_cell%f2(1, ix, iy) = Te_new
                        end if
                    end if
                end do
            end do

            do pp = 1, sp%num_par
                ix = sp%ix(1, pp)
                iy = sp%ix(2, pp)
                if (ix < lx .or. ix > ux .or. iy < ly .or. iy > uy) cycle
                sp%p(1,pp) = real( s2(ix,iy) * real(sp%p(1,pp), p_double) + bx2(ix,iy), p_k_part )
                sp%p(2,pp) = real( s2(ix,iy) * real(sp%p(2,pp), p_double) + by2(ix,iy), p_k_part )
                sp%p(3,pp) = real( s2(ix,iy) * real(sp%p(3,pp), p_double) + bz2(ix,iy), p_k_part )
            end do

            deallocate(W, Px2, Py2, Pz2, U2_2, s2, bx2, by2, bz2)

            end block

        case ( 3 )
            lx = lbound(electrons%axion_compton_info%depE%f3, 2)
            ux = ubound(electrons%axion_compton_info%depE%f3, 2)
            ly = lbound(electrons%axion_compton_info%depE%f3, 3)
            uy = ubound(electrons%axion_compton_info%depE%f3, 3)
            lz = lbound(electrons%axion_compton_info%depE%f3, 4)
            uz = ubound(electrons%axion_compton_info%depE%f3, 4)

            block

            real(p_double), allocatable :: W(:,:,:), Px3(:,:,:), Py3(:,:,:), Pz3(:,:,:), U2_3(:,:,:)
            real(p_double), allocatable :: s3(:,:,:), bx3(:,:,:), by3(:,:,:), bz3(:,:,:)
            allocate(W(lx:ux,ly:uy,lz:uz), Px3(lx:ux,ly:uy,lz:uz), Py3(lx:ux,ly:uy,lz:uz), Pz3(lx:ux,ly:uy,lz:uz), U2_3(lx:ux,ly:uy,lz:uz), s3(lx:ux,ly:uy,lz:uz), bx3(lx:ux,ly:uy,lz:uz), by3(lx:ux,ly:uy,lz:uz), bz3(lx:ux,ly:uy,lz:uz))
            W = 0.0_p_double; Px3 = 0.0_p_double; Py3 = 0.0_p_double; Pz3 = 0.0_p_double; U2_3 = 0.0_p_double
            s3 = 1.0_p_double; bx3 = 0.0_p_double; by3 = 0.0_p_double; bz3 = 0.0_p_double

            do pp = 1, sp%num_par
                ix = sp%ix(1, pp)
                iy = sp%ix(2, pp)
                iz = sp%ix(3, pp)
                if (ix < lx .or. ix > ux .or. iy < ly .or. iy > uy .or. iz < lz .or. iz > uz) cycle
                wnum = abs(real(sp%q(pp), p_double)) / max(Zcharge_abs, 1.0d-99)
                p1 = real(sp%p(1,pp), p_double)
                p2 = real(sp%p(2,pp), p_double)
                p3 = real(sp%p(3,pp), p_double)
                W(ix,iy,iz) = W(ix,iy,iz) + wnum
                Px3(ix,iy,iz) = Px3(ix,iy,iz) + wnum * p1
                Py3(ix,iy,iz) = Py3(ix,iy,iz) + wnum * p2
                Pz3(ix,iy,iz) = Pz3(ix,iy,iz) + wnum * p3
                U2_3(ix,iy,iz)= U2_3(ix,iy,iz)+ wnum * (p1*p1 + p2*p2 + p3*p3)
            end do

            ! Loop over the actual allocated cell bounds (including guard cells)
            do iz = lz, uz
                do iy = ly, uy
                    do ix = lx, ux
                        Wc = W(ix,iy,iz)
                        if (Wc <= 0.0_p_double) cycle

                        dE = frac * electrons%axion_compton_info%depE%f3(1, ix, iy, iz)
                        dPx = frac * electrons%axion_compton_info%depP%f3(1, ix, iy, iz)
                        dPy = frac * electrons%axion_compton_info%depP%f3(2, ix, iy, iz)
                        dPz = frac * electrons%axion_compton_info%depP%f3(3, ix, iy, iz)

                        Px = Px3(ix,iy,iz); Py = Py3(ix,iy,iz); Pz = Pz3(ix,iy,iz); U2 = U2_3(ix,iy,iz)
                        meanx = Px / Wc
                        meany = Py / Wc
                        meanz = Pz / Wc

                        var = U2 - (Px*Px + Py*Py + Pz*Pz) / max(Wc, tinyW)
                        if (var < 0.0_p_double) var = 0.0_p_double

                        dx = dPx / Wc
                        dy = dPy / Wc
                        dz = dPz / Wc

                        dE_bulk = 0.5_p_double/mu * ( 2.0_p_double*(dx*Px + dy*Py + dz*Pz) + Wc*(dx*dx + dy*dy + dz*dz) )
                        dE_rem = dE - dE_bulk
                        var_new = var + 2.0_p_double*mu*dE_rem
                        if (var_new < 0.0_p_double) var_new = 0.0_p_double

                        if (var > tinyVar) then
                            s = sqrt( var_new / var )
                        else
                            s = 1.0_p_double
                            var_new = var
                        end if

                        bx = (1.0_p_double - s)*meanx + dx
                        by = (1.0_p_double - s)*meany + dy
                        bz = (1.0_p_double - s)*meanz + dz

                        s3(ix,iy,iz) = s
                        bx3(ix,iy,iz) = bx
                        by3(ix,iy,iz) = by
                        bz3(ix,iy,iz) = bz

                        if (update_Te) then
                            if (associated(electrons%axion_compton_info%Te_cell%f3)) then
                                Te_est = mec2_keV * (var_new / max(Wc, tinyW)) / 3.0_p_double
                                Te_old = electrons%axion_compton_info%Te_cell%f3(1, ix, iy, iz)
                                Te_new = (1.0_p_double - relax)*Te_old + relax*Te_est
                                if (Te_new < floorT) Te_new = floorT
                                if (Te_new > ceilT) Te_new = ceilT
                                electrons%axion_compton_info%Te_cell%f3(1, ix, iy, iz) = Te_new
                            end if
                        end if
                    end do
                end do
            end do

            do pp = 1, sp%num_par
                ix = sp%ix(1, pp)
                iy = sp%ix(2, pp)
                iz = sp%ix(3, pp)
                if (ix < lx .or. ix > ux .or. iy < ly .or. iy > uy .or. iz < lz .or. iz > uz) cycle
                sp%p(1,pp) = real( s3(ix,iy,iz) * real(sp%p(1,pp), p_double) + bx3(ix,iy,iz), p_k_part )
                sp%p(2,pp) = real( s3(ix,iy,iz) * real(sp%p(2,pp), p_double) + by3(ix,iy,iz), p_k_part )
                sp%p(3,pp) = real( s3(ix,iy,iz) * real(sp%p(3,pp), p_double) + bz3(ix,iy,iz), p_k_part )
            end do

            deallocate(W, Px3, Py3, Pz3, U2_3, s3, bx3, by3, bz3)

            end block
        end select

    end subroutine apply_to_species

end subroutine apply_axion_compton_feedback
end module m_axions

!---------------------------------------------------------------------------------------------------
! advance photons
!---------------------------------------------------------------------------------------------------
subroutine advance_axions( this, electrons, positrons, emf, gdt, tid, n_threads)

    use m_system
    use m_parameters
    use m_axions_define, only : t_axions
    use m_species_define, only : p_ene_recalc
    use m_species_define_qed, only : t_species_qed
    use m_emf_define, only: t_emf
    use m_axions
    use m_species_define, only : p_ene_recalc
    use m_species, only : sort

    implicit none

    class( t_axions ), intent(inout) :: this
    class( t_species_qed ), intent(inout) :: electrons
    class( t_species_qed ), intent(inout) :: positrons
    class( t_emf ), intent( inout ) :: emf
    real(p_double), intent(in) :: gdt
    integer, intent(in) :: tid ! local thread id
    integer, intent(in) :: n_threads ! total number of threads

    integer :: chunk, ip0, ip1

    ! sanity check
    ! Code does not yet support OpenMP parallelism (mainly because of pair production)
    if ( n_threads > 1 ) then
      write(0,*) 'QED photon pusher is not yet supported with multiple threads per node'
      call abort_program()
    endif

    ! Time centered energy diagnostic is not yet implemented, force recalculation later
    ! when get_energy() is called
    ! (the push should not change photon energy)
    this%energy = p_ene_recalc

    ! range of photons for each thread
    chunk = ( this%num_par + n_threads - 1 ) / n_threads
    ip0 = tid * chunk + 1
    ip1 = min( (tid+1) * chunk, this%num_par )

    ! If there are no axion macroparticles, nothing to advance (avoid invalid ranges)
    if ( this%num_par == 0 ) return

    ! Create pairs and absorb photons
    if (this%if_pairprod) call create_pairs(this, electrons, positrons, emf, gdt, ip0, ip1)

    ! Warning : we want to avoid looping on photons that already decayed into pairs
    ip1 = min( (tid+1) * chunk, this%num_par )

    ! If there are no axion macroparticles, nothing to advance (avoid invalid ranges)
    if ( this%num_par == 0 ) return

    ! call bethe heitler pair creation
    if (this % betheheitler_info % if_betheheitler) then
        call create_pairs_betheheitler(this, electrons, positrons, gdt, ip0, ip1)
    endif

    ! Push photons. Boundary crossings will be checked at update_boundary
    select case ( this%coordinates )
        case default
            select case ( p_x_dim )
                case (1)
                    call advance_axions_1d( this, gdt, ip0, ip1 )

                case (2)
                    call advance_axions_2d( this, gdt, ip0, ip1 )

                case (3)
                    call advance_axions_3d( this, gdt, ip0, ip1 )

                case default
                    write(err_buf__,*) 'Not implemented for x_dim = ',p_x_dim;call err__("qed/os-axions.f03",2298)
                    call abort_program(p_err_invalid)

            end select
        case ( p_cylindrical_b )
            call advance_axions_cyl_2d( this, gdt, ip0, ip1 )
    end select

    ! Inverse process (absorption) for thermal-plasma axion bremsstrahlung
    if (electrons%axion_bremsstrahlung_info%if_bremsstrahlung) then
        if (electrons%axion_bremsstrahlung_info%if_inverse) then
            call axion_bremsstrahlung_absorb( this, electrons, gdt, ip0, ip1 )
        end if
     if (electrons%axion_bremsstrahlung_info%apply_feedback) then
         call apply_axion_brems_feedback( electrons )
     end if
    end if

    ! Inverse process (absorption) for thermal-bath axion Compton-like channel
    if (electrons%axion_compton_info%if_compton) then
        if (electrons%axion_compton_info%if_inverse) then
            call axion_compton_absorb( this, electrons, gdt, ip0, ip1 )
        end if
        if (electrons%axion_compton_info%apply_feedback) then
            call apply_axion_compton_feedback( electrons )
        end if
    end if

end subroutine advance_axions
!---------------------------------------------------------------------------------------------------

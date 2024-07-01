#include "os-config.h"
#include "os-preprocess.fpp"

module m_axions

#include "memory/memory.h"

use m_system
use m_parameters
use m_math

use m_species_define, only: p_cell_near
use m_axions_define, only: t_axions
use m_species_define, only : p_ene_recalc
use m_photons, only : ntrim_qed
use m_emf_define, only: t_emf
use m_qed
use m_vdf_define
use m_qed_coulomb
use m_species_charge, only: deposit_rho

implicit none

private

public advance_axions_1d, advance_axions_2d, advance_axions_3d, advance_axions_cyl_2d
!public ----photondecay----ntrim_qed

contains

!---------------------------------------------------------------------------------------------------
subroutine advance_axions_1d( this, gdt, i0, i1)
!---------------------------------------------------------------------------------------------------

implicit none

    integer, parameter :: rank = 2

    ! dummy variables
    class( t_axions ),    intent(inout) :: this
    real(p_double),   intent(in) :: gdt
    integer,   intent(in) :: i0, i1

    real(p_k_part) :: rdx
    integer :: i, pp, np, ptrcur

    real(p_k_part), dimension(rank,p_cache_size) :: xbuf
    integer, dimension(rank,p_cache_size)        :: dxi
    real(p_k_part), dimension(p_cache_size)      :: rp, rpp
    real(p_k_part) :: dt

    ! executable statements
    rdx = real( 1.0_p_double/this%dx(1), p_k_part )
    dt = real( gdt, p_k_part )
    
    ! advance position of axions i0 to i1 in chunks of p_cache_size
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
            this%x(1,pp)  = xbuf(1,i) - dxi(1,i)
            this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)

            pp = pp + 1
        end do

    enddo

end subroutine advance_photons_1d

!---------------------------------------------------------------------------------------------------
subroutine advance_axions_2d( this, gdt, i0, i1)
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine advance_axions_3d( this, gdt, i0, i1)
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine advance_axions_cyl_2d( this, gdt, i0, i1 )
!---------------------------------------------------------------------------------------------------  

!---------------------------------------------------------------------------------------------------
subroutine create_photonpair( this, photon, emf, dt, i0, i1)
!---------------------------------------------------------------------------------------------------

    use m_species_define_qed, only : t_species_qed
    use m_species_comm, only : t_part_idx, remove_particles
    use m_photons_diagnostics
    
    implicit none
    
    integer, parameter :: rank = 2
    
    ! dummy variables
    
    class( t_photons ), intent(inout) :: photons
    class( t_axions ), intent(inout) :: this
    class( t_emf), intent(in) :: emf
    real(p_double),   intent(in) :: dt
    integer,   intent(in) :: i0, i1
    
    real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep
    integer :: i, pp, np, ptrcur

    ! attributes for new particles
    real(p_k_part)   , dimension(p_p_dim) :: p_new

    ! index of photons to remove
    type(t_part_idx) :: local_idx
    
    ! dummy variables for chi parameter diagnostic
    real(p_k_part) :: q_g
    integer, dimension(p_x_dim) :: ix_g
    
    ! Physical constants taken from NIST reference data
    real(p_double), parameter :: hbar  = 1.054571726d-34     ! J s
    real(p_double), parameter :: me    = 9.10938291d-31      ! kg
    real(p_double), parameter :: cl    = 299792458.0         ! m/s
    real(p_double), parameter :: alpha = 7.2973525698d-3     ! (no units)

    real(p_double), parameter :: PI    = 3.141592653589793d0
    
    ! see notes a + emf -- gamma (x2?) 

    

!---------------------------------------------------------------------------------------------------
subroutine advance_axions( this, photons, emf, gdt, tid, n_threads)
!---------------------------------------------------------------------------------------------------

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
    class( t_photons ), intent(inout) :: photons
    class( t_emf ), intent(inout) :: emf
    real(p_double), intent(in) :: gdt
    integer, intent(in) :: tid    ! local thread id
    integer, intent(in) :: n_threads  ! total number of threads
    
    integer :: chunk, ip0, ip1
    
    ! sanity check
    ! Code does not yet support OpenMP parallelism (mainly because of pair production)
    if ( n_threads > 1 ) then
      write(0,*) 'QED axion pusher is not yet supported with multiple threads per node'
      call abort_program()
    endif
    
    ! Time centered energy diagnostic is not yet implemented, force recalculation later
    ! when get_energy() is called
    ! (the push should not change axion energy)
    this%energy = p_ene_recalc
    
    ! range of axions for each thread
    chunk = ( this%num_par + n_threads - 1 ) / n_threads
    ip0    = tid * chunk + 1
    ip1    = min( (tid+1) * chunk, this%num_par )
    
    ! Create photon pairs and absorb axions
    
    
    
    ! Push axions. Boundary crossings will be checked at update_boundary
    
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
                    ERROR('Not implemented for x_dim = ',p_x_dim)
                    call abort_program(p_err_invalid)
        
            end select
        case ( p_cylindrical_b )
            call advance_axions_cyl_2d( this, gdt, ip0, ip1 )
    end select
  
end subroutine advance_photons




















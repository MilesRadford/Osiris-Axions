#include "os-config.h"
#include "os-preprocess.fpp"

module m_axions

#include "memory/memory.h"

use m_system
use m_parameters
use m_math

! GitHUb seems pretty rubbish

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
public ntrim_qed

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

end subroutine advance_axions_1d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine advance_axions_2d( this, gdt, i0, i1)
!---------------------------------------------------------------------------------------------------

    implicit none

    integer, parameter :: rank = 2

    ! dummy variables
    class( t_axions ),    intent(inout) :: this
    real(p_double),   intent(in) :: gdt
    integer,   intent(in) :: i0, i1

    real(p_k_part), dimension(p_x_dim) :: rdx
    integer :: i, pp, np, ptrcur

    real(p_k_part), dimension(rank,p_cache_size) :: xbuf
    integer, dimension(rank,p_cache_size)        :: dxi
    real(p_k_part), dimension(p_cache_size)      :: rp, rpp
    real(p_k_part) :: dt

    ! executable statements
    rdx(1) = real( 1.0_p_double/this%dx(1), p_k_part )
    rdx(2) = real( 1.0_p_double/this%dx(2), p_k_part )
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
            xbuf(1,i) = this%x(1,pp) + this%p(1,pp) * rpp(i) * rdx(1)
            xbuf(2,i) = this%x(2,pp) + this%p(2,pp) * rpp(i) * rdx(2)
            dxi(1,i) = ntrim_qed( xbuf(1,i) )
            dxi(2,i) = ntrim_qed( xbuf(2,i) )

            pp = pp + 1
        end do


        ! copy data from buffer to species data trimming positions
        pp = ptrcur
        do i = 1, np
            this%x(1,pp)  = xbuf(1,i) - dxi(1,i)
            this%x(2,pp)  = xbuf(2,i) - dxi(2,i)
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

    class( t_axions ),    intent(inout) :: this
    real(p_double),   intent(in) :: gdt
    integer,   intent(in) :: i0, i1

    real(p_k_part), dimension(p_x_dim) :: rdx
    integer :: i, pp, np, ptrcur

    real(p_k_part), dimension(rank,p_cache_size) :: xbuf
    integer, dimension(rank,p_cache_size)        :: dxi
    real(p_k_part), dimension(p_cache_size)      :: rp, rpp
    real(p_k_part) :: dt

    ! executable statements
    rdx(1) = real( 1.0_p_double/this%dx(1), p_k_part )
    rdx(2) = real( 1.0_p_double/this%dx(2), p_k_part )
    rdx(3) = real( 1.0_p_double/this%dx(3), p_k_part )
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
            this%x(1,pp)  = xbuf(1,i) - dxi(1,i)
            this%x(2,pp)  = xbuf(2,i) - dxi(2,i)
            this%x(3,pp)  = xbuf(3,i) - dxi(3,i)
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
    real(p_double), dimension(p_x_dim)   :: xmin_g
    real(p_k_part), dimension(rank,p_cache_size) :: xbuf
    integer, dimension(rank,p_cache_size)        :: dxi
    real(p_k_part), dimension(p_cache_size)      :: rgamma

    real(p_k_part) :: dt_dx1
    integer :: gix2, shift_ix2
    real(p_double) :: dr, rdr, tmp

    real(p_double) :: x2_new, x3_new, r_old, r_new
    real(p_double) :: r_shift

    ! executable statements

    xmin_g(1) = this%g_box( p_lower, 1 )
    xmin_g(2) = this%g_box( p_lower, 2 )

    shift_ix2 = this%my_nx_p(p_lower, 2) - 2
    dr  = this%dx(p_r_dim)
    rdr = 1.0_p_double/dr

    dt_dx1 = real( gdt / this%dx(1), p_k_part )

    ! Since r = 0 is at the center of cell ix2 = 1 the radial position of particles will be
    ! - Positions defined with regard to the center of the cell (odd interpolation)
    !       r = ( (gix2-1) + x2 ) * dr
    ! - Positions defined with regard to the corner of the cell (even interpolation)
    !       r = ( (gix2-1) + x2 - 0.5 ) * dr

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
            gix2  = this%ix(2,pp)  + shift_ix2
            r_old = ( this%x(2,pp) + (gix2 - r_shift) ) * dr

            ! Particle is pushed in 3D like fashion
            x2_new    = r_old + ( this%p(2,pp) * rgamma(i) ) * gdt
            x3_new    =         ( this%p(3,pp) * rgamma(i) ) * gdt

            r_new     = sqrt( x2_new**2 + x3_new**2 )

            ! Convert new position to "cell" type position

            ! this is a protection against roundoff for cold plasmas
            if ( r_old == r_new ) then
            xbuf(2,i) = this%x(2,pp)
            else
            xbuf(2,i) = real( r_new * rdr - ( gix2 - r_shift ) , p_k_part)
            endif

            ! Correct p_r and p_\theta to conserve angular momentum
            ! there is a potential division by zero here
            tmp          = 1.0_p_double / r_new
            this%p(2,pp) = real( ( this%p(2,pp) * x2_new + this%p(3,pp)*x3_new ) * tmp, p_k_part )
            this%p(3,pp) = real( this%p(3,pp) * ( r_old * tmp ), p_k_part )

            dxi(1,i) = ntrim_qed( xbuf(1,i) )
            dxi(2,i) = ntrim_qed( xbuf(2,i) )

            pp = pp + 1
        end do

        ! copy data from buffer to species data trimming positions
        pp = ptrcur
        do i = 1, np
        this%x(1,pp)  = xbuf(1,i) - dxi(1,i)
        this%x(2,pp)  = xbuf(2,i) - dxi(2,i)
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


end module m_axions

!---------------------------------------------------------------------------------------------------
! advance axions
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
    class( t_emf ), intent( inout )  ::  emf
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
  
end subroutine advance_axions
!---------------------------------------------------------------------------------------------------

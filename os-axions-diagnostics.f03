#include "os-config.h"
#include "os-preprocess.fpp"

module m_axions_diagnostics

#include "memory/memory.h"

use m_system
use m_parameters
use m_math

use m_node_conf
use m_time_step
use m_space, only : t_space
use m_grid_define, only : t_grid
use m_vdf_define
use m_diagfile, only : p_time_length

use m_axions_define, only: t_axions

implicit none

private

interface n_emit_hist
    module procedure n_emit_hist
end interface

interface report_n_emit
    module procedure report_n_emit
end interface 

interface chi_hist
    module procedure chi_hist
end interface

interface report_chi
    module procedure report_chi
end interface 

public :: n_emit_hist, report_n_emit

contains
subroutine n_emit_hist(this, ix)

    implicit none
  
    class(t_axions), intent(inout) :: this
    integer, dimension(p_x_dim), intent(in) :: ix
    integer :: n_emit
  
    select case ( p_x_dim )
    case ( 1 )
        n_emit = this%n_emit%f1(1, ix(1))
        this%n_emit%f1(1, ix(1)) = n_emit + 1
    case ( 2 )
        n_emit = this%n_emit%f2(1, ix(1), ix(2))
        this%n_emit%f2(1, ix(1), ix(2)) = n_emit + 1
    case ( 3 )
        n_emit = this%n_emit%f3(1, ix(1), ix(2), ix(3))
        this%n_emit%f3(1, ix(1), ix(2), ix(3)) = n_emit + 1
    end select
    
end subroutine n_emit_hist

!------------------------------------------------------------------------------------------

subroutine report_n_emit( this, g_space, grid, no_co, tstep, t )

    use m_diagnostic_utilities
    use stringutil, only: idx_string
  
    implicit none
  
    class(t_axions), intent(inout) :: this
    type( t_space ),       intent(in) :: g_space
    class( t_grid ),       intent(in) :: grid
    class( t_node_conf ),  intent(in) :: no_co
    type( t_time_step ),   intent(in) :: tstep
    real(p_double),        intent(in) :: t

    type( t_vdf_report ) :: info
    integer :: ns, nd

    ! radiated energy spectrum diagnostic
    if (test_if_report( tstep, this%ndump_fac_n_emit )) then
        
        ns = n(tstep)
        nd = ndump(tstep) * this%ndump_fac_n_emit

        ! Prepare the metadata
        info%name = 'n_emit'
        info%ndump = nd

        info%xname  = (/'x1', 'x2', 'x3'/)

        info%xlabel = (/'x_1', 'x_2', 'x_3'/)
        info%xunits = (/'c / \omega_a', 'c / \omega_a', 'c / \omega_a'/)

        info%time_units = '1 / \omega_a'
        info%dt         = dt(tstep)

        info%label = 'N'
        info%units = ''

        info%n = ns
        info%t = t

        info%fileLabel = ''
        info%path  = trim(path_mass) // 'QED' // p_dir_sep // 'n_emit' // &
                     p_dir_sep // trim(this%name) // p_dir_sep
        info%filename = 'n_emit_a-' // trim(this%name) // '-' // idx_string( ns/nd, p_time_length )

        ! write the file
        call this % n_emit % write( info, 1, g_space, grid, no_co )
        
        call this % n_emit % zero()
    endif
end subroutine report_n_emit
















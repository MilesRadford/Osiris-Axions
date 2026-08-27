module m_axion_field_diag_define

#include "os-config.h"
#include "os-preprocess.fpp"

  use m_parameters
  use m_input_file, only : t_input_file
  use m_vdf_define, only : t_vdf_report
  use m_diagnostic_utilities

  implicit none
  private

  integer, public :: diag_axion_ev = 0

  character(len=10), dimension(12), parameter, public :: p_axion_report_quants = &
    (/ 'phi       ' , 'phidot    ' , 'aold      ' , 'anew      ' , &
       'ae1       ' , 'ae2       ' , 'ae3       ' , 'ab1       ' , 'ab2       ' , 'ab3       ' , &
       'edotb     ' , 'emsrcmag  ' /)

  integer, parameter, public :: p_phi=1, p_phidot=2, p_aold=3, p_anew=4, p_ae1=5, p_ae2=6, p_ae3=7, &
                                p_ab1=8, p_ab2=9, p_ab3=10, p_edotb=11, p_emsrcmag=12

  integer, parameter :: p_max_diag_len = 10

  type, public :: t_axion_diag
    character( len = p_max_diag_len ), dimension(:), pointer :: report_quants => null()
    type(t_vdf_report), pointer :: reports => null()
  contains
    procedure :: avail_report_quants => avail_report_quants_axion
    procedure :: init_report_quants  => init_report_quants_axion
    procedure :: read_input          => read_input_diag_axion
    procedure :: init                => init_diag_axion
    procedure :: cleanup             => cleanup_diag_axion
  end type t_axion_diag

  interface
    ! Forward the diagnostic-quantity count query to the compiled implementation module.
    function avail_report_quants_axion( this ) result(nq)
      import t_axion_diag
      class(t_axion_diag), intent(in) :: this
      integer :: nq
    end function avail_report_quants_axion
  end interface

  interface
    ! Forward diagnostic-quantity initialization to the compiled implementation module.
    subroutine init_report_quants_axion( this )
      import t_axion_diag
      class(t_axion_diag), intent(inout) :: this
    end subroutine init_report_quants_axion
  end interface

  interface
    ! Forward axion diagnostic input parsing to the compiled implementation module.
    subroutine read_input_diag_axion( this, input_file )
      import t_axion_diag, t_input_file
      class(t_axion_diag), intent(inout) :: this
      class(t_input_file), intent(inout) :: input_file
    end subroutine read_input_diag_axion
  end interface

  interface
    ! Forward axion diagnostic initialization to the compiled implementation module.
    subroutine init_diag_axion( this )
      import t_axion_diag
      class(t_axion_diag), intent(inout) :: this
    end subroutine init_diag_axion
  end interface

  interface
    ! Forward axion diagnostic cleanup to the compiled implementation module.
    subroutine cleanup_diag_axion( this )
      import t_axion_diag
      class(t_axion_diag), intent(inout) :: this
    end subroutine cleanup_diag_axion
  end interface

end module m_axion_field_diag_define

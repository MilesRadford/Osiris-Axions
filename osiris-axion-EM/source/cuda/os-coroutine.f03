module m_cosimulation

  use m_system
  use m_parameters
  use m_simulation , only : t_simulation
  use iso_c_binding

  implicit none

  private

  interface
    subroutine init_cosimulation() bind(C)
      use iso_c_binding
    end subroutine
  end interface

  interface
    subroutine handoff_context_cpp(path)  bind(c)
      use iso_c_binding
      character(kind=c_char), intent(in) :: path(*)
    end subroutine handoff_context_cpp
  end interface

  interface
    subroutine start_sim_itteration_cpp()  bind(c)
      use iso_c_binding
    end subroutine start_sim_itteration_cpp
  end interface

  interface
    subroutine test_str_cpp(path)  bind(c)
      use iso_c_binding
      character(kind=c_char), intent(in) :: path(*)
    end subroutine test_str_cpp
  end interface


  interface
    subroutine execute_sim1_iter_cpp()  bind(c)
      use iso_c_binding
    end subroutine execute_sim1_iter_cpp
  end interface

  interface
    subroutine execute_sim2_iter_cpp()  bind(c)
      use iso_c_binding
    end subroutine execute_sim2_iter_cpp
  end interface

  integer :: num_contexts
  integer :: current_context
  character(len=50) :: context_return_value

  public :: setup_cosimluation_magic, start_sim_itteration, handoff_context, test_str, context_return_value
  public :: execute_sim1_iter, execute_sim2_iter

  contains

  subroutine setup_cosimluation_magic( sim  )
    implicit none

    class(t_simulation) :: sim
    #ifdef ENABLE_COROUTINE
      print *, "Starting setup_cosimluation_magic (F03 subroutine)"
      call init_cosimulation()
      print *, "At end of setup_cosimluation_magic (F03 subroutine)"

      num_contexts = 2
      current_context = -1
    #endif
    !call handoff_context_cpp("COOLEST_CONTEXT")

  end subroutine

  subroutine start_sim_itteration( )
    implicit none

    #ifdef ENABLE_COROUTINE
      call start_sim_itteration_cpp()
    #endif

  end subroutine

  subroutine handoff_context(context_name)
    implicit none
    character(len=*), intent(in) :: context_name

    character(len=len_trim(context_name)+1) :: lcontext_name
    #ifdef ENABLE_COROUTINE
      ! turn Fortan string into one that C/C++ likes (zero terminated)
      lcontext_name = trim(context_name)//char(0)
      call handoff_context_cpp(lcontext_name)
    #endif

  end subroutine handoff_context

  subroutine test_str()

    implicit none

    #ifdef ENABLE_COROUTINE
      call test_str_cpp(context_return_value)
    #endif

  end subroutine test_str
  
  subroutine execute_sim1_iter

    implicit none

    #ifdef ENABLE_COROUTINE
      call execute_sim1_iter_cpp
      call test_str()
    #endif

  end subroutine

  subroutine execute_sim2_iter

    implicit none

    #ifdef ENABLE_COROUTINE
      !print *, "^^^^^^^^^^^^^^^^^ 1"
      call execute_sim2_iter_cpp
      !print *, "^^^^^^^^^^^^^^^^^ 2"
      call test_str()
      !print *, "^^^^^^^^^^^^^^^^^ 3"
    #endif
  
  end subroutine

end module
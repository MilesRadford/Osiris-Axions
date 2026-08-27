module m_handoff
  use m_system
  use m_parameters

  use iso_c_binding

  implicit none

  private

  interface handoff_context_new
    module procedure handoff_context_new1
    module procedure handoff_context_new2
  end interface

  !------------------------------------------------------------------------------
  !------ Interface(s) to make C++ function(s) easily callable by Fortran -------
  !------------------------------------------------------------------------------
  interface
    subroutine handoff_context_cpp(path)  bind(c)
      use iso_c_binding
      character(kind=c_char), intent(in) :: path(*)
    end subroutine handoff_context_cpp
  end interface
  !------------------------------------------------------------------------------

  
  integer :: num_contexts
  integer :: current_context
  character(len=50) :: context_return_value

  public :: handoff_context_new, get_trailing_int_in_str

  contains

  subroutine handoff_context_new1(context_name)
    implicit none
    character(len=*), intent(in) :: context_name

    character(len=len_trim(context_name)+1) :: lcontext_name

    #ifdef ENABLE_COROUTINE
      ! turn Fortan string into one that C/C++ likes (zero terminated)
      lcontext_name = trim(context_name)//char(0)
      call handoff_context_cpp(lcontext_name)
    #endif

  end subroutine handoff_context_new1


  subroutine handoff_context_new2(context_name, idx)
    implicit none
    character(len=*), intent(in) :: context_name
    integer :: idx

    character(len=len_trim(context_name)+100) :: lcontext_name

    #ifdef COSIM
      ! Make a context description that looks like:
      !     context_name##idx
      ! Example:
      !     context_name == 'ADVANCED_DEPOSIT', idx = 3 would give
      !     a context description of: ADVANCED_DEPOSIT##3
      lcontext_name = trim(context_name)//'##'//trim(int2str(idx))//char(0)
      call handoff_context_cpp(lcontext_name)
    #endif
  end subroutine handoff_context_new2


  !-----------------------------------------------------------------------------------------
  ! Get an integer tag at the end of a sting. The tag demlimter is '##'.
  !  The following example will result in val1==77:
  !     integer :: val1
  !     character(len=50) :: str
  !     str = 'ADV_DEPOSIT##77'
  !     val1 = get_trailing_int_in_str( str)
  !
  ! Optionally:
  ! We can an optional index of where in the strng to stop searching. On return, this value 
  !  is overwriten with the index of where to stop searching the string such that
  !  the '##' tag just read is excluded. This value can be passed back in to subsquent calls 
  !  to 'get_int_tag' to get multiple tags from the string.
  !
  !   Example: 
  !     The following will result in val1==77,val2==343:
  !
  !       integer :: idx, val1, val2
  !       character(len=50) :: str
  !       str = 'ADV_DEPOSIT##343##77'
  !       idx = len(str)
  !       val1 = get_trailing_int_in_str( str,str_idx)
  !       val2 = get_trailing_int_in_str( str,str_idx)
  !
  !-----------------------------------------------------------------------------------------
  function get_trailing_int_in_str(str,start_index_,none_value_)
    implicit none
    character(len=*),intent(in) :: str
    integer,optional            :: start_index_
    integer,optional            :: none_value_

    integer                     :: get_trailing_int_in_str

    integer :: start_index, none_value, stridx, int_val

    if(present(start_index_)) then
      start_index = start_index_
    else
      start_index = len(str)
    endif
    if(present(none_value_)) then
      none_value = none_value_
    else
      none_value = -1
    endif

    int_val = none_value

    stridx = INDEX( str(0:start_index), '##', .true.)
    !print *, 'stridx', str(0:start_index), stridx
    if(stridx > 0) then
      call str2int( trim(str(stridx+1:start_index)), int_val )
   endif
   get_trailing_int_in_str = int_val

   if(present(start_index_) .and. stridx > 0) then
      start_index_ = stridx-2
   endif

  end function

  !-----------------------------------------------------------------------------------------
  !  Return a string representing the integer 'int'
  !-----------------------------------------------------------------------------------------  
  function int2str(int)
    implicit none
    integer :: int
    character(len=20) :: int2str

    write (int2str, *) int
    int2str = adjustl(int2str)
  end function

  !-----------------------------------------------------------------------------------------
  !  If 'str' is a string representing an integer, return that integer.
  !-----------------------------------------------------------------------------------------
  subroutine str2int(str,int)
    implicit none
    character(len=*),intent(in) :: str
    integer,intent(inout)       :: int
    integer :: stat

    read(str,*,iostat=stat)  int
  end subroutine

end module
module m_axion_fps_context

#include "os-config.h"
#include "os-preprocess.fpp"

  use m_vdf_define, only : t_vdf

  implicit none
  private

  logical, save :: fps_active = .false.
  logical, save :: fps_part_grid_center = .false.
  type(t_vdf), pointer, save :: fps_e1 => null()
  type(t_vdf), pointer, save :: fps_b1 => null()

  public :: set_axion_fps_context
  public :: clear_axion_fps_context
  public :: axion_fps_is_active
  public :: get_axion_fps_fields
  public :: axion_fps_part_grid_center

contains

subroutine set_axion_fps_context( e1, b1, part_grid_center )
  class(t_vdf), target, intent(in) :: e1, b1
  logical, intent(in), optional :: part_grid_center

  fps_e1 => e1
  fps_b1 => b1
  fps_active = .true.
  if ( present(part_grid_center) ) then
    fps_part_grid_center = part_grid_center
  else
    fps_part_grid_center = .false.
  endif
end subroutine set_axion_fps_context

subroutine clear_axion_fps_context()
  fps_active = .false.
  fps_part_grid_center = .false.
  nullify( fps_e1 )
  nullify( fps_b1 )
end subroutine clear_axion_fps_context

logical function axion_fps_is_active()
  axion_fps_is_active = fps_active .and. associated(fps_e1) .and. associated(fps_b1)
end function axion_fps_is_active

subroutine get_axion_fps_fields( e1, b1 )
  type(t_vdf), pointer, intent(out) :: e1, b1
  e1 => fps_e1
  b1 => fps_b1
end subroutine get_axion_fps_fields

logical function axion_fps_part_grid_center()
  axion_fps_part_grid_center = fps_part_grid_center
end function axion_fps_part_grid_center

end module m_axion_fps_context

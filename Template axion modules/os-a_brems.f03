#include "os-preprocess.fpp"
#include "os-config.h"

module m_a_brems

#include "memory/memory.h"
  !use m_random_class

  use m_parameters
  use m_species_define_qed
  use m_system
  use m_node_conf
  use m_random
  use m_vdf_define
  use m_vdf
  use m_species_define

implicit none

private

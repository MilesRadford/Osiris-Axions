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

public :: n_emit_hist, report_n_emit, chi_hist, report_chi

contains


















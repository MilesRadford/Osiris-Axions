! this module provides a data structure which abstracts
! boosted grid diagnostics. specifically it handles the conversion
! of lab and boosted coordinates


#include "os-config.h"
#include "os-preprocess.fpp"

module m_boosted_diag

#include "memory/memory.h"

use m_system
use m_parameters

use stringutil
use m_vdf_define
use m_vdf_memory
use m_diagnostic_utilities
use m_grid_define
use m_space
use m_time_step
use m_utilities
use m_boosted_grid_mgr
use m_input_file
use m_parameters
use m_restart

implicit none


character(len=*), parameter :: p_boosted_diag_rst_id = "boosted diag rst data - 0x0000"

! configurations for boosted diagnostics to work with
! these change how the data is stored, updated, and lorentz-transformed
integer, parameter, public :: p_boosted_fields = 0
integer, parameter, public :: p_boosted_4vector = 1


! note: in all the code below, the 1 and 2 indexing makes life a lot easier.
! just remember that the simulation frame is frame 1, and the lab frame
! which we're transforming to is frame 2.


!-------------------------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------------------------
type :: t_boosted_diag

    ! shared object to keep track of coordinates for each calculation
    type(t_boosted_grid_mgr), pointer :: grid_mgr => null()

    ! spatial dimensionality of the data
    integer :: x_dim_

    ! current options: p_boosted_4vector, p_boosted_fields
    integer :: diag_type

    integer :: ndump_fac

    ! previous and current data
    ! the implementation could be generalized by instead keeping a slot
    ! for generic data here (e.g. type(t_vdf) :: cur_data, prev_data),
    ! but the current implementation makes life easier since there aren't any other boosted
    ! quantities we would want to look at.

    ! type(t_vdf), pointer :: cur_e, cur_b
    type(t_vdf), pointer :: prev_e => null(), prev_b => null()
    type(t_vdf), pointer :: prev_4vector => null()

    ! buffers storing partial data of a time slice in frame 2
    ! for each timeslice
    type(t_vdf), pointer :: buffered_e => null(), buffered_b => null()
    type(t_vdf), pointer :: buffered_4vector => null()

    ! current and previous leftmost global cell index (needed for interp)
    real(p_double) :: xmin1_prev

    ! size of data buffer, calculated to fit max amount of data given
    integer     :: data_buf_size

    ! number of data currently stored in buffer
    integer :: num_data_in_buf

    ! index of grid % boosted_cell_coords that buffer index 0 corresponds to
    integer :: buf_start_idx

    ! directory where this data will be stored
    character(len=1024) :: basepath

    ! number of components to write to disk
    integer :: num_diag_components

    ! indices of the specific components to write to disk (size = ndiags)
    integer, dimension(6) :: diag_components

    integer :: num_diag_components_to_buffer
    integer, dimension(6) :: diag_components_to_buffer

    ! names of all components
    character(len=128), dimension(6) :: component_names

    contains

    procedure :: restart_write => restart_write_boosted_diag
    procedure :: restart_read => restart_read_boosted_diag

    procedure :: init => init_boosted_diag
    procedure :: compute_data_buf_size => compute_data_buf_size_boosted_diag
    procedure :: init_grid_mgr => init_grid_mgr_boosted_diag
    procedure :: cleanup => cleanup_boosted_diag
    procedure :: read_input => read_input_boosted_diag

    procedure :: parse_reports_string => parse_reports_string_boosted_diag
    procedure :: init_components_to_buffer => init_components_to_buffer_boosted_diag
    procedure :: get_dir_fname_label => get_dir_fname_label_boosted_diag

    procedure :: create_files => create_files_boosted_diag
    procedure :: create_files_component => create_files_component_boosted_diag

    procedure :: report => report_boosted_diag
    procedure :: report_component => report_component_boosted_diag

    procedure :: extract => extract_boosted_diag

    procedure :: allocate_objs => allocate_objs_boosted_diag

end type
!-------------------------------------------------------------------------------------------------


public :: t_boosted_diag


contains



!-------------------------------------------------------------------------------------------------
subroutine restart_write_boosted_diag( this, restart_handle )

    use m_restart

    implicit none

    class( t_boosted_diag ), intent(in) :: this
    type( t_restart_handle ), intent(inout) :: restart_handle

    character(len=*), parameter :: err_msg = 'error writing restart data for boosted_diag object.'
    integer :: ierr

    restart_io_wr( p_boosted_diag_rst_id, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%xmin1_prev, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%buf_start_idx, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%num_data_in_buf, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    select case( this%diag_type )

    case( p_boosted_4vector )

        call this % buffered_4vector % write_checkpoint( restart_handle )
        call this % prev_4vector % write_checkpoint( restart_handle )

    case( p_boosted_fields )

        call this % buffered_e % write_checkpoint( restart_handle )
        call this % buffered_b % write_checkpoint( restart_handle )

        call this % prev_e % write_checkpoint( restart_handle )
        call this % prev_b % write_checkpoint( restart_handle )

    end select

end subroutine restart_write_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine restart_read_boosted_diag( this, restart_handle )

    use m_restart

    implicit none

    class( t_boosted_diag ), intent(inout) :: this
    type( t_restart_handle ), intent(in) :: restart_handle

    character(len=*), parameter :: err_msg = 'error writing restart data for boosted_diag object.'
    integer :: ierr
    character(len=len(p_boosted_diag_rst_id)) :: rst_id

    integer, dimension(3) :: data_buf_nx


    restart_io_rd( rst_id, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd )

    ! check if restart file is compatible
    if ( rst_id /= p_boosted_diag_rst_id ) then
        ERROR('Corrupted restart file, or restart file from incompatible binary (emf)')
        call abort_program(p_err_rstrd)
    endif

    restart_io_rd( this%xmin1_prev, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_rd( this%buf_start_idx, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_rd( this%num_data_in_buf, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    select case( this%diag_type )

    case( p_boosted_4vector )

        call this % buffered_4vector % read_checkpoint( restart_handle )

        data_buf_nx = this%buffered_4vector%nx_
        data_buf_nx(1) = this%data_buf_size

        call reshape_array_vdf( this % buffered_4vector, data_buf_nx )

        call this % prev_4vector % read_checkpoint( restart_handle )

    case( p_boosted_fields )

        call this % buffered_e % read_checkpoint( restart_handle )
        call this % buffered_b % read_checkpoint( restart_handle )

        data_buf_nx = this%buffered_e%nx_
        data_buf_nx(1) = this%data_buf_size

        call reshape_array_vdf( this % buffered_e, data_buf_nx )
        call reshape_array_vdf( this % buffered_b, data_buf_nx )

        call this % prev_e % read_checkpoint( restart_handle )
        call this % prev_b % read_checkpoint( restart_handle )

    end select

end subroutine restart_read_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine init_grid_mgr_boosted_diag( this, grid, g_space, dt1, tmax1 )

    implicit none

    class(t_boosted_diag), intent(inout) :: this

    class( t_grid ), intent(in) :: grid
    type( t_space ),     intent(in) :: g_space
    real(p_double), intent(in) :: dt1, tmax1

    call this % grid_mgr % set_space1_params( grid, g_space, dt1, tmax1 )
    call this % grid_mgr % init()

end subroutine init_grid_mgr_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine init_boosted_diag( this, &
                    grid, g_space, dt1, tmax1, &
                    ndump_fac_global, basepath, ref_vdf, &
                    gc_num, restart, restart_handle, init_grid_mgr_ )

    implicit none

    class(t_boosted_diag), intent(inout) :: this

    class( t_grid ), intent(in) :: grid
    type( t_space ),     intent(in) :: g_space
    real(p_double), intent(in) :: dt1, tmax1
    integer, intent(in) :: ndump_fac_global
    character(len=*), intent(in) :: basepath
    class(t_vdf), intent(in)           :: ref_vdf
    integer, dimension (2, p_x_dim), intent(in) :: gc_num
    logical, intent(in) :: restart
    type( t_restart_handle ), intent(in) :: restart_handle
    logical, intent(in), optional :: init_grid_mgr_

    integer :: ndump_total
    integer, dimension(3) :: data_buf_nx
    logical :: init_grid_mgr

    if( present( init_grid_mgr_ )) then
        init_grid_mgr = init_grid_mgr_
    else
        init_grid_mgr = .true.
    endif

    ndump_total = this%ndump_fac * ndump_fac_global

    ! this is an input because in quasi-3D we only call it for one of the 2m+1 diags
    if( init_grid_mgr ) then
        call this % init_grid_mgr( grid, g_space, dt1, tmax1 )

        ! calculate the size to be used for the data buffer
        ! this is the most data that a local node would ever need to store
        call this % compute_data_buf_size( ndump_total )
    endif

    ! allocate the appropriate vdfs -- needs to happen before either restart
    ! or initialization later in this function
    call this % allocate_objs()

    ! non-restart initialization
    this%basepath = basepath
    this%x_dim_ = ref_vdf%x_dim_
    call this % init_components_to_buffer()

    ! the rest of the initialization doesn't apply if we're restarting
    if( restart ) then

        call this % restart_read( restart_handle )

    else

        this%xmin1_prev = this%grid_mgr%xmin1

        ! init vars not passed in
        this % num_data_in_buf = 0
        this % buf_start_idx = 1

        data_buf_nx(1) = this%data_buf_size
        data_buf_nx(2:3) = ref_vdf%nx_(2:3)

        select case( this%diag_type )

        case( p_boosted_4vector )

            call this%prev_4vector    %new( ref_vdf%x_dim_, 4, ref_vdf%nx_, gc_num, ref_vdf%dx_, .true. )
            call this%buffered_4vector%new( ref_vdf%x_dim_, 4, data_buf_nx, gc_num, ref_vdf%dx_, .true. )

        case( p_boosted_fields )

            call this%prev_e%new( ref_vdf%x_dim_, 3, ref_vdf%nx_, gc_num, ref_vdf%dx_, .true. )
            call this%prev_b%new( ref_vdf%x_dim_, 3, ref_vdf%nx_, gc_num, ref_vdf%dx_, .true. )

            call this%buffered_e%new( ref_vdf%x_dim_, 3, data_buf_nx, gc_num, ref_vdf%dx_, .true. )
            call this%buffered_b%new( ref_vdf%x_dim_, 3, data_buf_nx, gc_num, ref_vdf%dx_, .true. )

        case default
            SCR_ROOT( "ERROR in os-boosted-diag: unrecognized boosted diagnostic type request.d")
            stop

        end select

    endif

end subroutine init_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine compute_data_buf_size_boosted_diag( this, ndump_total )

    implicit none

    class(t_boosted_diag), intent(inout) :: this
    integer, intent(in) :: ndump_total

    integer :: i1, tmp_buf_size, start_idx, final_idx

    this%data_buf_size = 1

    do i1 = 1, ceiling( this%grid_mgr%nt1 * 1.0 / ndump_total )

        start_idx = ( i1 - 1 ) * ndump_total + 1
        final_idx = min( i1 * ndump_total, this%grid_mgr%nt1 )

        tmp_buf_size = sum( this%grid_mgr%num_t2_points_at_t1(start_idx : final_idx))
        this%data_buf_size = max( this%data_buf_size, tmp_buf_size )
    enddo

end subroutine
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
! allocate all objects except our grid manager
! the grid manager is allocated separately because in quasi-3D we only
! want to allocate grid manager which is shared by the different diags.
!-------------------------------------------------------------------------------------------------
subroutine allocate_objs_boosted_diag( this )

    implicit none

    class(t_boosted_diag), intent(inout) :: this

    select case( this%diag_type )

    case( p_boosted_4vector )

        call alloc( this%prev_4vector )
        call alloc( this%buffered_4vector )

    case( p_boosted_fields )

        call alloc( this%prev_e )
        call alloc( this%prev_b )
        call alloc( this%buffered_e )
        call alloc( this%buffered_b )

    end select

end subroutine allocate_objs_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine cleanup_boosted_diag( this )

    implicit none

    class(t_boosted_diag), intent(inout) :: this

    if( associated(this%grid_mgr) ) then
        call this%grid_mgr%cleanup()
        deallocate( this%grid_mgr )
    endif

    if (this%ndump_fac > 0) then

        select case( this%diag_type )
        case( p_boosted_4vector )

            call this%prev_4vector%cleanup()
            call this%buffered_4vector%cleanup()

            call freemem( this%prev_4vector )
            call freemem( this%buffered_4vector )

        case( p_boosted_fields )

            call this%prev_e%cleanup()
            call this%prev_b%cleanup()
            call this%buffered_e%cleanup()
            call this%buffered_b%cleanup()

            call freemem( this%prev_e )
            call freemem( this%prev_b )
            call freemem( this%buffered_e )
            call freemem( this%buffered_b )

        end select

    endif

end subroutine cleanup_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine read_input_boosted_diag( this, input_file, quant_list, sim_gamma, &
                                     diag_type, if_use_boosted_diag )

    implicit none

    class(t_boosted_diag), intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file
    character(len=*), dimension(:), intent(in) :: quant_list
    real(p_double), intent(in) :: sim_gamma
    integer, intent(in) :: diag_type
    logical, intent(out) :: if_use_boosted_diag

    integer        :: ndump_fac
    real(p_double) :: boost_signed_gamma
    real(p_double) :: boost_tmin, boost_tmax, boost_dt
    real(p_double) :: boost_xmin, boost_xmax
    integer        :: boost_nx
    logical        :: boost_if_move
    real(p_double) :: boost_move_vel

    integer :: ierr

    character( len = p_max_reports_len ), dimension( 6 ) :: boosted_reports

    namelist /nl_boosted_diag/ ndump_fac, boost_signed_gamma, &
                           boost_tmin, boost_tmax, boost_dt, &
                           boost_xmin, boost_xmax, boost_nx, &
                           boost_if_move, boost_move_vel, &
                           boosted_reports

    ndump_fac = 0
    boost_signed_gamma = -sim_gamma
    boost_tmin = 0.
    boost_tmax = 0.
    boost_dt = 0.
    boost_xmin = 0.
    boost_xmax = 0.
    boost_nx = 0
    boost_if_move = .true.
    boost_move_vel = 1.
    boosted_reports = "-"

    this%diag_type = diag_type

    ! disable boosted diags by default 
    this%ndump_fac = 0 
    
    ! return and say we didn't find a boosted diag if there's no boosted_diag namelist 
    if_use_boosted_diag = .false.
    if (trim(input_file % get_section_name()) /= 'nl_boosted_diag') return

    ! Get namelist text from input file
    call get_namelist( input_file, "nl_boosted_diag", ierr )

    if (ierr == 0) then
        read (input_file%nml_text, nml = nl_boosted_diag, iostat = ierr)
        if (ierr /= 0) then
            SCR_ROOT( "Error reading boosted_diag parameters" )
            SCR_ROOT( "aborting..." )
            stop
        endif
    else
        if (ierr < 0) then
            SCR_ROOT( "Error reading boosted_diag parameters" )
            SCR_ROOT( "aborting..." )
            stop
        else
            if (disp_out(input_file)) then
                SCR_ROOT(" - no diagnostics specified")
            endif
        endif
    endif

    this % ndump_fac = ndump_fac
    if_use_boosted_diag = ( ndump_fac > 0 )

    if (.not. if_use_boosted_diag) return

    if( .not. associated( this % grid_mgr ) ) then
        allocate( this%grid_mgr )
    else
        SCR_ROOT( 'ERROR in os-boosted-diag: grid manager was already allocated in read_input.')
        SCR_ROOT( 'This should never happen.')
        stop
    endif


    call this % parse_reports_string( boosted_reports, quant_list )

    call this % grid_mgr % set_space2_params( boost_signed_gamma, &
                                         boost_tmin, boost_tmax, boost_dt, &
                                         boost_xmin, boost_xmax, boost_nx, &
                                         boost_if_move, boost_move_vel )

end subroutine read_input_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
! boosted_reports_str: list of reports to track. only these quantities will be
! buffered. e.g. ( "e1", "b3" )
! quant_list: list of ALL options for components for this diagnostic. example:
! ( "e1", "e2", "e3", "b1", "b2", "b3" )
!-------------------------------------------------------------------------------------------------
subroutine parse_reports_string_boosted_diag( this, boosted_reports_requested, quant_list )

    implicit none

    class(t_boosted_diag), intent(inout) :: this
    character(len=*), dimension(:), intent(in) :: boosted_reports_requested
    character(len=*), dimension(:), intent(in) :: quant_list

    integer :: component
    character(len=64) :: requested_report

    integer :: i, j

    this%num_diag_components = 0

    ! SCR_ROOT( 'requested: ', boosted_reports_requested )
    ! SCR_ROOT( 'quant_list: ', quant_list )

    ! loop thru all requested quants
    do i = 1, size( boosted_reports_requested )

        requested_report = boosted_reports_requested(i)

        if ( trim(requested_report ) == '-' ) exit

        component = -1

        ! loop through all available components
        ! make sure we match one, and if so store the component
        do j = 1, size(quant_list)

            if ( trim(requested_report) == trim(quant_list(j)) ) then
                component = j
                exit
            endif

        enddo

        if ( component == -1 ) then
            SCR_ROOT( "ERROR in os-boosted-diag: invalid boosted quantity requested: ", requested_report )
            stop
        endif

        this%diag_components(i) = component
        this%num_diag_components = this%num_diag_components + 1

    enddo

    this%component_names(1:size(quant_list)) = quant_list

end subroutine parse_reports_string_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine init_components_to_buffer_boosted_diag( this )

    implicit none

    class(t_boosted_diag), intent(inout) :: this

    integer :: i

    this%num_diag_components_to_buffer = 0

    ! compute necessary diag components here
    select case( this%diag_type )

    case( p_boosted_4vector)

        i = 1

        if( any( this%diag_components(1:this%num_diag_components) .eq. 1 ) &
            .or. any( this%diag_components(1:this%num_diag_components) .eq. 2 ) ) then

            this%diag_components_to_buffer(i) = 1
            i = i + 1

            this%diag_components_to_buffer(i) = 2
            i = i + 1

            this%num_diag_components_to_buffer = this%num_diag_components_to_buffer + 2

        endif

        if( any( this%diag_components(1:this%num_diag_components) .eq. 3 ) ) then
            this%diag_components_to_buffer(i) = 3
            this%num_diag_components_to_buffer = this%num_diag_components_to_buffer + 1
            i = i + 1
        endif

        if( any( this%diag_components(1:this%num_diag_components) .eq. 4 ) ) then
            this%diag_components_to_buffer(i) = 4
            this%num_diag_components_to_buffer = this%num_diag_components_to_buffer + 1
            i = i + 1
        endif

    end select

end subroutine init_components_to_buffer_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine get_dir_fname_label_boosted_diag( this, component, t2_idx, dir, fname, label )

    implicit none

    class(t_boosted_diag), intent(in) :: this
    integer, intent(in)              :: t2_idx
    integer, intent(in)              :: component
    character(len=1024 ), intent(out) :: dir
    character(len=128 ), intent(out) :: fname
    character(len=128 ), intent(out) :: label

    character(len=128) :: name

    name = this%component_names( component )

    dir = trim(this%basepath) // p_dir_sep // trim(name) // p_dir_sep

    label = trim(name)

    fname = trim( name ) // "-" // idx_string( t2_idx, 6 )

end subroutine get_dir_fname_label_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine create_files_boosted_diag( this, grid )

    implicit none

    class(t_boosted_diag), intent(in) :: this
    class( t_grid ),      intent(in) :: grid

    integer :: i

    if( this%ndump_fac .le. 0 ) return

    do i = 1, this%num_diag_components
        call this % create_files_component( this%diag_components(i), grid )
    enddo

end subroutine create_files_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine create_files_component_boosted_diag( this, component, grid )

    implicit none

    class(t_boosted_diag), intent(in) :: this
    integer, intent(in) :: component
    class( t_grid ),      intent(in) :: grid

    character(len=1024) :: dir
    character(len=128 ) :: filename
    character(len=128 ) :: label

    integer :: i2, j
    real(p_double) :: t2, xmin_current2, xmax_current2
    integer, dimension(3) :: dataset_dims

    class( t_diag_file ), allocatable   :: diagfile
    type( t_diag_dataset ) :: dset
    character(len=12), dimension(3) :: xnames, xlabels, xunits
    character(len=8) :: label_
    integer, parameter :: izero = ichar('0')

    xnames  = (/'x1', 'x2', 'x3'/)
    xlabels = (/'x_1', 'x_2', 'x_3'/)
    xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)

    do i2 = 0, this%grid_mgr%nt2 - 1

        ! get the time and spatial bounds for this slice
        t2 = this%grid_mgr%tmin2 + i2 * this%grid_mgr%dt2

        if( this%grid_mgr%if_move2 ) then
            xmin_current2 = this%grid_mgr%xmin2 + t2
        else
            xmin_current2 = this%grid_mgr%xmin2
        endif

        xmax_current2 = xmin_current2 + ( this%grid_mgr%xmax2 - this%grid_mgr%xmin2 )

        ! this allocates the diagfile -- must be deallocated later
        call create_diag_file( diagfile )

        call this % get_dir_fname_label( component, i2, dir, filename, label )

        diagfile % ftype    = p_diag_grid

        diagfile % filepath = dir
        diagfile % filename = filename
        diagfile % name     = label

        ! Iteration info
        diagfile%iter%n         = i2
        diagfile%iter%t         = t2

        ! Grid info
        diagfile%grid%ndims = this % x_dim_


        select case( this%diag_type )

        case( p_boosted_4vector )

            diagfile%grid%units  = 'n_0'

            select case( component )

            case(1)
                label_ = '\rho'

                if ( p_x_dim == 1 ) then
                    diagfile%grid%units  = 'e \omega_p / c'
                else
                    diagfile%grid%units  = 'e \omega_p^'//char(izero+p_x_dim)// &
                            '/ c^'//char(izero+p_x_dim)
                endif

            case(2)
                label_ = 'j_1'
            case(3)
                label_ = 'j_2'
            case(4)
                label_ = 'j_3'
            end select

            ! diagfile%units = 'm_e c \omega_p e^{-1}'


        case( p_boosted_fields )

            diagfile%grid%units = 'm_e c \omega_p e^{-1}'

            select case( component )
            case(1)
                label_ = 'E_1'
            case(2)
                label_ = 'E_2'
            case(3)
                label_ = 'E_3'
            case(4)
                label_ = 'B_1'
            case(5)
                label_ = 'B_2'
            case(6)
                label_ = 'B_3'
            end select

        end select

        diagFile%grid%label = trim( label_ )
        diagfile%grid%name  = diagfile%name

        ! handle i=1 separately
        diagfile%grid%count(1)     = this%grid_mgr%nx2
        diagfile%grid%axis(1)%type = diag_axis_linear
        diagfile%grid%axis(1)%min  = xmin_current2
        diagfile%grid%axis(1)%max  = xmax_current2

        diagfile%grid%axis(1)%name  = xnames(1)
        diagfile%grid%axis(1)%label = xlabels(1)
        diagfile%grid%axis(1)%units = xunits(1)

        ! handle transverse directions
        do j = 2, this % x_dim_
            diagfile%grid%count(j)     = this%grid_mgr%gnxperp(j-1)
            diagfile%grid%axis(j)%type = diag_axis_linear
            diagfile%grid%axis(j)%min  = this%grid_mgr%gxperp_min(j-1)
            diagfile%grid%axis(j)%max  = this%grid_mgr%gxperp_max(j-1)

            diagfile%grid%axis(j)%name  = xnames(j)
            diagfile%grid%axis(j)%label = xlabels(j)
            diagfile%grid%axis(j)%units = xunits(j)
        enddo

        call diagFile % open( p_diag_create, grid%io%comm )

        dataset_dims(1) = this%grid_mgr%nx2
        dataset_dims(2:3) = this%grid_mgr%gnxperp(1:2)

        ! don't chunk right now, use contiguous dataset
        call diagFile % start_cdset( label, this%x_dim_, &
                dataset_dims, &
                diag_type_real( p_single ), dset )

        call diagFile % close()

        deallocate( diagFile )

    enddo

end subroutine create_files_component_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine report_boosted_diag( this, tstep, grid )

    implicit none

    class(t_boosted_diag), intent(inout) :: this
    type( t_time_step ), intent(in) :: tstep
    class( t_grid ), intent(in) :: grid

    integer :: i, n1

    n1 = n( tstep )

    if( n1 .eq. 0 ) then

        SCR_ROOT('Creating boosted files...')
        call this % create_files( grid )
        SCR_ROOT('Done.')

    else

        if( .not. test_if_report( tstep, this%ndump_fac ) ) return

        do i = 1, this%num_diag_components
            call this % report_component( this%diag_components(i), grid )
        enddo

        ! adjust position of all buffers
        this%buf_start_idx = this%buf_start_idx + this%num_data_in_buf
        this%num_data_in_buf = 0

    endif

end subroutine report_boosted_diag
!-------------------------------------------------------------------------------------------------




!-------------------------------------------------------------------------------------------------
! name: prefix for the data
! component: component of the vdf where the data is stored
! i indexes time, j indexes space
!-------------------------------------------------------------------------------------------------
subroutine report_component_boosted_diag( this, component, grid )

    implicit none

    class(t_boosted_diag), intent(inout) :: this
    class( t_grid ),      intent(in) :: grid

    character(len=1024) :: dir
    character(len=128 ) :: filename
    character(len=128 ) :: label

    class( t_diag_file ), allocatable   :: diagfile
    type( t_diag_dataset ) :: dset
    type( t_diag_chunk ) :: chunk

    ! type(t_vdf) :: write_buf
    type(t_vdf), pointer :: data_buf

    real(p_single), dimension(:),     pointer :: write_buf_1d
    real(p_single), dimension(:,:),   pointer :: write_buf_2d
    real(p_single), dimension(:,:,:), pointer :: write_buf_3d


    integer :: component, i2, j, idx, cur_i2, cur_j2, data_idx
    integer :: start_idx, num_data, cur_idx, min_j2, buff_idx
    integer :: cur_idx_within_databuffer

    ! i2 values obtained in the unsorted buffer
    integer, dimension(:), pointer :: i2_vals, num_data_at_i2

    ! integer, dimension(:,:), pointer :: buffer_indices_for_each_i2

    ! sorted indices of those i2 indices
    integer, dimension(:), pointer :: sorted_i2_indices

    ! set pointers to null
    write_buf_1d => null(); write_buf_2d => null(); write_buf_3d => null()
    i2_vals => null(); num_data_at_i2 => null(); sorted_i2_indices => null()

    select case( this%diag_type )

    case( p_boosted_fields )

        ! components 1 to 3 store e
        if( component <= 3 ) then
            data_buf => this%buffered_e
            data_idx = component

        ! components 4 to 6 store b
        else
           data_buf => this%buffered_b
           data_idx = component - 3
       endif

    case( p_boosted_4vector )

       data_buf => this%buffered_4vector
       data_idx = component

    end select

    start_idx = 0

    call alloc( num_data_at_i2, (/ this%grid_mgr%nt2 /) )
    num_data_at_i2(:) = 0

    ! first compute sorted_i2_indices
    if( this%num_data_in_buf > 0 ) then

        call alloc( i2_vals, (/ this%num_data_in_buf /) )
        call alloc( sorted_i2_indices, (/ this%num_data_in_buf /) )

        do j=1, this%num_data_in_buf
            idx = this%buf_start_idx + (j-1)
            i2 = this%grid_mgr%boosted_cell_coords(idx)%t2_idx
            i2_vals(j) = i2
            num_data_at_i2( i2 ) = num_data_at_i2( i2 ) + 1
        enddo

        call heapsort( i2_vals, sorted_i2_indices )

    endif

    call create_diag_file( diagfile )


    ! loop through all lab time slices to reconstruct
    do i2 = 1, this%grid_mgr%nt2

        num_data = num_data_at_i2(i2)
        min_j2 = 1 ! default value in case num_data is 0

        ! we still need to proceed if num_data == 0 because other processes
        ! may have data; in this case, we need to tell the mpi process
        ! doing the writing that we don't have any data, otherwise it
        ! will get stuck waiting for us. however, if we don't have any data
        ! then we don't need to do anything other than give this info to
        ! the writing process.

        if( num_data > 0 ) then

            ! allocate the appropriate dimension buffer. this will hold
            ! all the sorted data for a given i2 until we write it.
            ! the dimensions have to exactly match the ones we eventually
            ! write, which is why we need to keep reallocating this array.
            select case( this%x_dim_ )
            case(1)
                call alloc( write_buf_1d, (/ num_data /) )
                chunk % data = c_loc( write_buf_1d )
            case(2)
                call alloc( write_buf_2d, (/ num_data, data_buf%nx_(2) /) )
                chunk % data = c_loc( write_buf_2d )
            case(3)
                call alloc( write_buf_3d, (/ num_data, data_buf%nx_(2), data_buf%nx_(3) /) )
                chunk % data = c_loc( write_buf_3d )
            end select

        endif

        ! find the minimum j2 attained for all the available data (the heapsort
        ! doesn't respect the j2 ordering posessed by the data prior to sorting )
        do j = 1, num_data
            ! cur_idx_within_databuffer  = sorted_i2_indices( start_idx + j )
            cur_idx  = this%buf_start_idx + sorted_i2_indices( start_idx + j ) - 1
            cur_j2 = this%grid_mgr%boosted_cell_coords(cur_idx)%x2_idx

            if( j .eq. 1 ) then
                min_j2 = cur_j2
            else
                min_j2 = min( min_j2, cur_j2 )
            endif

        enddo

        ! count all the data in the buffer which matches this i2
        do j = 1, num_data

            cur_idx_within_databuffer  = sorted_i2_indices( start_idx + j )
            cur_idx  = this%buf_start_idx + sorted_i2_indices( start_idx + j ) - 1

            cur_i2 = this%grid_mgr%boosted_cell_coords(cur_idx)%t2_idx

            if( i2 .ne. cur_i2 ) then
                print *,  "Fatal error: wrong i2 in os-boosted-diag.f03"
                stop
            endif

            cur_j2 = this%grid_mgr%boosted_cell_coords(cur_idx)%x2_idx

            buff_idx = cur_j2 - min_j2 + 1

            ! finally load the write buffer
            select case( this%x_dim_ )

            case(1)
                write_buf_1d(buff_idx)       = data_buf%f1( data_idx, cur_idx_within_databuffer )
            case(2)
                write_buf_2d(buff_idx, 1:data_buf%nx_(2))    = data_buf%f2( data_idx, cur_idx_within_databuffer, 1:data_buf%nx_(2) )
            case(3)
                write_buf_3d(buff_idx, 1:data_buf%nx_(2), 1:data_buf%nx_(3)) = data_buf%f3( &
                                data_idx, cur_idx_within_databuffer, 1:data_buf%nx_(2), 1:data_buf%nx_(3) )

            end select

            ! ! useful debug info -- keeping this here for now
            ! if( num_data > 0 .and. component .eq. 1 .and. i2 .eq. 5 .and. this%diag_type .eq. p_boosted_4vector ) then
            !     print *, 'debug ', mpi_node(), cur_j2, min_j2, num_data, write_buf_2d(buff_idx, 128)
            ! endif

        enddo

        ! update start index to account for progress through the
        ! cell coordinates we've looped through
        start_idx = start_idx + num_data

        call this % get_dir_fname_label( component, i2 - 1, dir, filename, label )

        diagfile % filename = filename
        diagfile % filepath = dir
        diagfile % name = label
        dset % name = label
        dset % ndims = this%x_dim_

        chunk % count(1) = num_data
        chunk % count(2:3) = this%grid_mgr%nxperp(1:2)

        chunk % start(1) = min_j2 - 1
        chunk % start(2:3) = this%grid_mgr%nxperp_min(1:2) - 1

        chunk % stride(1:3) = 1

        ! does not currently work
        ! call diagfile % open( p_diag_update, grid%io%comm, DIAG_MPIIO_INDEPENDENT )

        ! call diagfile % open( p_diag_update, grid%io%comm, DIAG_MPI )
        call diagfile % open( p_diag_update, grid%io%comm )

        call diagfile % open_cdset( dset )

        ! this is necessary because this metadata is only read by
        ! the one process that opens it
        dset % count(1) = this%grid_mgr%nx2
        dset % count(2:3) = this%grid_mgr%gnxperp(1:2)
        dset % data_type = p_diag_float32

        call diagfile % write_par_cdset( dset, chunk )

        call diagfile % end_cdset( dset )

        call diagfile % close()

        if( num_data > 0 ) then
            call freemem( write_buf_1d )
            call freemem( write_buf_2d )
            call freemem( write_buf_3d )
        endif

    enddo ! i2 loop

    deallocate( diagFile )

    call freemem( num_data_at_i2 )

    if( this%num_data_in_buf > 0 ) then

        call freemem( i2_vals )
        call freemem( sorted_i2_indices )

    endif

end subroutine report_component_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
! this is called at every timestep to extract all data out of the current and
! previous timestep that corresponds to a lab frame timeslice
!-------------------------------------------------------------------------------------------------
subroutine extract_boosted_diag( this, tstep, t1, grid, g_space, vdf1, vdf2 )

    implicit none

    class(t_boosted_diag), intent(inout) :: this
    type(t_time_step), intent(in) :: tstep
    real(p_double), intent(in) :: t1
    class( t_grid ), intent(in) :: grid
    type( t_space ),     intent(in) :: g_space
    class(t_vdf), intent(in)           :: vdf1
    class(t_vdf), intent(in), optional :: vdf2

    integer :: i1, j, buff_idx, cell_idx !, n
    integer :: x1_idx, x1_idx_prev
    real(p_double) :: x1_cell, t1_cell, xmin1
    real(p_double) :: dx1, dt1
    real(p_double) :: t1_p, t1_m, x1_p, x1_m
    real(p_double) :: coef_tp, coef_tm, coef_xp, coef_xm
    real(p_double) :: coef_pp, coef_pm, coef_mp, coef_mm
    real(p_double) :: gamma, beta

    real(p_k_fld), dimension(3) :: e_pp, e_pm, e_mp, e_mm
    real(p_k_fld), dimension(3) :: b_pp, b_pm, b_mp, b_mm
    real(p_k_fld), dimension(3) :: buffered_e, buffered_b
    real(p_k_fld), dimension(4) :: x_pp, x_pm, x_mp, x_mm
    real(p_k_fld), dimension(4) :: buffered_4vector

    integer :: i2, i3

    if( this%ndump_fac .le. 0 ) return

    gamma = this%grid_mgr%gamma
    beta  = this%grid_mgr%beta

    ! i1 has the property t = (i1 - 1) * dt1 + tmin
    i1 = n( tstep ) + 1

    buff_idx = this%num_data_in_buf + 1

    cell_idx = this%grid_mgr%starting_cell_idx_at_t1(i1)

    xmin1 = xmin( g_space, 1 )  + this%grid_mgr%dx1 * (grid%my_nx(p_lower,1) - 1)

    dx1 = this%grid_mgr%dx1
    dt1 = this%grid_mgr%dt1

    t1_p = t1
    t1_m = t1 - dt1

    do j = 1, this % grid_mgr % num_t2_points_at_t1(i1)

        t1_cell = this%grid_mgr%boosted_cell_coords(cell_idx)%t1
        x1_cell = this%grid_mgr%boosted_cell_coords(cell_idx)%x1

        ! find the spatial cells to use for this interpolation
        x1_idx = floor( ( x1_cell - xmin1 ) / dx1 )
        x1_idx_prev = floor( ( x1_cell - this%xmin1_prev ) / dx1 )

        ! gridpoints of the interpolation
        x1_m = xmin1 + x1_idx * dx1
        x1_p = x1_m + dx1

! for serious debugging
#if 0
        if( t1_cell > t1 .or. t1_cell < t1_m ) then
            print *, 'ERROR: os-boosted-diag: invalid cell coordinate', t1_cell, t1, t1_m
            stop
        endif

        if( min( x1_idx, x1_idx_prev ) < 0 &
            .or. max( x1_idx, x1_idx_prev ) > grid%my_nx(3,1) ) then

             print *, 'ERROR: os-boosted-diag: spatial interpolation index is out of bounds ', &
             x1_cell, &
             xmin1, &!  xmax1, &
             ! t1, dt1, &
             this%grid_mgr%xmin1, this%grid_mgr%xmax1, &
             x1_idx, x1_idx_prev, grid%my_nx(3,1), &
             dx1, xmin( g_space, 1 ), xmax( g_space, 1 ), g_space%total_i_moved(1)

            stop
        endif
#endif

        ! m, p --> minus, plus
        coef_tp = (t1_p - t1_cell) / dt1
        coef_xp = (x1_p - x1_cell) / dx1
        coef_tm = (t1_cell - t1_m) / dt1
        coef_xm = (x1_cell - x1_m) / dx1

        coef_pp = coef_tp * coef_xp
        coef_pm = coef_tp * coef_xm
        coef_mp = coef_tm * coef_xp
        coef_mm = coef_tm * coef_xm

        select case( this%diag_type )

        case( p_boosted_fields )

            select case( this%x_dim_ )

            case(1)

                e_pp = vdf1%f1(:,x1_idx + 1)
                e_pm = vdf1%f1(:,x1_idx )
                e_mp = this%prev_e%f1(:,x1_idx_prev + 1)
                e_mm = this%prev_e%f1(:,x1_idx_prev )

                b_pp = vdf2%f1(:,x1_idx + 1)
                b_pm = vdf2%f1(:,x1_idx )
                b_mp = this%prev_b%f1(:,x1_idx_prev + 1)
                b_mm = this%prev_b%f1(:,x1_idx_prev )

                call extract_fields_helper( e_pp, e_pm, e_mp, e_mm, &
                                            b_pp, b_pm, b_mp, b_mm, &
                                            coef_pp, coef_pm, coef_mp, coef_mm, &
                                            gamma, beta, &
                                            buffered_e, buffered_b )

                this%buffered_e%f1(:,buff_idx) = buffered_e
                this%buffered_b%f1(:,buff_idx) = buffered_b

            case(2)

                do i2 = 1, this%grid_mgr%nxperp(1)
                    ! SCR_ROOT('i2', i2)
                    e_pp = vdf1%f2(:, x1_idx + 1, i2 )
                    e_pm = vdf1%f2(:, x1_idx, i2 )
                    e_mp = this%prev_e%f2(:, x1_idx_prev + 1, i2 )
                    e_mm = this%prev_e%f2(:, x1_idx_prev, i2 )

                    b_pp = vdf2%f2(:, x1_idx + 1, i2 )
                    b_pm = vdf2%f2(:, x1_idx, i2 )
                    b_mp = this%prev_b%f2(:, x1_idx_prev + 1, i2)
                    b_mm = this%prev_b%f2(:, x1_idx_prev, i2 )

                    call extract_fields_helper( e_pp, e_pm, e_mp, e_mm, &
                                                b_pp, b_pm, b_mp, b_mm, &
                                                coef_pp, coef_pm, coef_mp, coef_mm, &
                                                gamma, beta, &
                                                buffered_e, buffered_b )

                    this%buffered_e%f2(:, buff_idx, i2) = buffered_e
                    this%buffered_b%f2(:, buff_idx, i2) = buffered_b
                enddo

            case(3)

                do i2 = 1, this%grid_mgr%nxperp(1)
                    do i3 = 1, this%grid_mgr%nxperp(2)

                        e_pp = vdf1%f3(:, x1_idx + 1, i2, i3)
                        e_pm = vdf1%f3(:, x1_idx, i2, i3)
                        e_mp = this%prev_e%f3(:, x1_idx_prev + 1, i2, i3)
                        e_mm = this%prev_e%f3(:, x1_idx_prev, i2, i3)

                        b_pp = vdf2%f3(:, x1_idx + 1, i2, i3)
                        b_pm = vdf2%f3(:, x1_idx, i2, i3)
                        b_mp = this%prev_b%f3(:, x1_idx_prev + 1, i2, i3)
                        b_mm = this%prev_b%f3(:, x1_idx_prev, i2, i3)

                        call extract_fields_helper( e_pp, e_pm, e_mp, e_mm, &
                                                    b_pp, b_pm, b_mp, b_mm, &
                                                    coef_pp, coef_pm, coef_mp, coef_mm, &
                                                    gamma, beta, &
                                                    buffered_e, buffered_b )

                        this%buffered_e%f3(:, buff_idx, i2, i3) = buffered_e
                        this%buffered_b%f3(:, buff_idx, i2, i3) = buffered_b

                    enddo
                enddo

            end select


        case( p_boosted_4vector )

            select case( this%x_dim_ )

            case(1)

                x_pp = vdf1%f1(:,x1_idx + 1)
                x_pm = vdf1%f1(:,x1_idx )
                x_mp = this%prev_4vector%f1(:,x1_idx_prev + 1)
                x_mm = this%prev_4vector%f1(:,x1_idx_prev )

                call extract_4vector_helper( x_pp, x_pm, x_mp, x_mm, &
                                             coef_pp, coef_pm, coef_mp, coef_mm, &
                                             gamma, beta, &
                                             buffered_4vector )

                this%buffered_4vector%f1(:,buff_idx) = buffered_4vector

            case(2)

                do i2 = 1, this%grid_mgr%nxperp(1)

                    x_pp = vdf1%f2(:,x1_idx + 1, i2)
                    x_pm = vdf1%f2(:,x1_idx, i2 )
                    x_mp = this%prev_4vector%f2(:,x1_idx_prev + 1, i2)
                    x_mm = this%prev_4vector%f2(:,x1_idx_prev, i2 )

                    call extract_4vector_helper( x_pp, x_pm, x_mp, x_mm, &
                                                 coef_pp, coef_pm, coef_mp, coef_mm, &
                                                 gamma, beta, &
                                                 buffered_4vector )

                    this%buffered_4vector%f2(:,buff_idx,i2) = buffered_4vector

                enddo

            case(3)

                do i2 = 1, this%grid_mgr%nxperp(1)
                    do i3 = 1, this%grid_mgr%nxperp(2)

                        x_pp = vdf1%f3(:,x1_idx + 1, i2, i3)
                        x_pm = vdf1%f3(:,x1_idx, i2, i3 )
                        x_mp = this%prev_4vector%f3(:,x1_idx_prev + 1, i2, i3)
                        x_mm = this%prev_4vector%f3(:,x1_idx_prev, i2, i3 )

                        call extract_4vector_helper( x_pp, x_pm, x_mp, x_mm, &
                                                     coef_pp, coef_pm, coef_mp, coef_mm, &
                                                     gamma, beta, &
                                                     buffered_4vector )

                        this%buffered_4vector%f3(:,buff_idx,i2, i3) = buffered_4vector

                    enddo
                enddo

            end select

        end select

        cell_idx = cell_idx + 1
        buff_idx = buff_idx + 1

    enddo

    ! account for new stored data
    this%num_data_in_buf = this%num_data_in_buf + this%grid_mgr%num_t2_points_at_t1(i1)

    this%xmin1_prev = xmin1

    select case( this%diag_type )

    case( p_boosted_fields )

        call check_vdf_compatability( this%prev_e, vdf1 )
        call check_vdf_compatability( this%prev_b, vdf2 )

        call this%prev_e%copy( vdf1 )
        call this%prev_b%copy( vdf2 )

    case( p_boosted_4vector )

        call check_vdf_compatability( this%prev_4vector, vdf1 )

        call this%prev_4vector%copy( vdf1 )

    end select

end subroutine extract_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
! don't mind me, just a nice little helper function at your service
! this function / usage of it is not optimized
!-------------------------------------------------------------------------------------------------
subroutine extract_fields_helper(    e_pp, e_pm, e_mp, e_mm, &
                                     b_pp, b_pm, b_mp, b_mm, &
                                     coef_pp, coef_pm, coef_mp, coef_mm, &
                                     gamma, beta, &
                                     buffered_e, buffered_b )

    implicit none

    real(p_k_fld),  dimension(3), intent(in) :: e_pp, e_pm, e_mp, e_mm
    real(p_k_fld),  dimension(3), intent(in) :: b_pp, b_pm, b_mp, b_mm
    real(p_double), intent(in) :: coef_pp, coef_pm, coef_mp, coef_mm
    real(p_double), intent(in)               :: gamma, beta
    real(p_k_fld),  dimension(3), intent(out) :: buffered_e, buffered_b

    real(p_k_fld),  dimension(3) :: e_interp, b_interp

    ! first index is time, second is space
    e_interp = real( coef_pp * e_mm + coef_pm * e_mp &
                   + coef_mp * e_pm + coef_mm * e_pp, p_k_fld )

    b_interp = real( coef_pp * b_mm + coef_pm * b_mp &
                   + coef_mp * b_pm + coef_mm * b_pp, p_k_fld )

    ! now populate the buffers
    ! Eperp ' = gamma[ E_perp - v \times B ] for boost into lab
    buffered_e(1) = e_interp(1)
    buffered_e(2) = real( gamma * ( e_interp(2) - beta * b_interp(3) ), p_k_fld )
    buffered_e(3) = real( gamma * ( e_interp(3) + beta * b_interp(2) ), p_k_fld )

    buffered_b(1) = b_interp(1)
    buffered_b(2) = real( gamma * ( b_interp(2) + beta * e_interp(3) ), p_k_fld )
    buffered_b(3) = real( gamma * ( b_interp(3) - beta * e_interp(2) ), p_k_fld )

end subroutine extract_fields_helper
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine extract_4vector_helper(   x_pp, x_pm, x_mp, x_mm, &
                                     coef_pp, coef_pm, coef_mp, coef_mm, &
                                     gamma, beta, &
                                     buffered_4vector )

    implicit none

    real(p_k_fld),  dimension(4), intent(in) :: x_pp, x_pm, x_mp, x_mm
    real(p_double), intent(in) :: coef_pp, coef_pm, coef_mp, coef_mm
    real(p_double), intent(in)               :: gamma, beta

    real(p_k_fld),  dimension(4), intent(out) :: buffered_4vector

    real(p_k_fld),  dimension(4) :: x_interp

    ! first index is time, second is space
    x_interp = real( coef_pp * x_mm + coef_pm * x_mp &
                   + coef_mp * x_pm + coef_mm * x_pp, p_k_fld )

    buffered_4vector(1) = real( gamma * ( x_interp(1) - beta * x_interp(2) ), p_k_fld )
    buffered_4vector(2) = real( gamma * ( x_interp(2) - beta * x_interp(1) ), p_k_fld )
    buffered_4vector(3:4) = x_interp(3:4)

end subroutine extract_4vector_helper
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine check_vdf_compatability( vdf1, vdf2 )

    implicit none

    type(t_vdf), intent(in)  :: vdf1
    type(t_vdf), intent(in)  :: vdf2

    logical :: success
    integer :: x_dim

    x_dim = vdf1%x_dim_

    success = ( (vdf1%x_dim_ .eq. vdf2%x_dim_) &
            .and. all( vdf1%nx_(1:x_dim) .eq. vdf2%nx_(1:x_dim) ) &
            .and. all( vdf1%gc_num_(:,1:x_dim) .eq. vdf2%gc_num_(:,1:x_dim) ) )

    if( .not. success ) then
        ERROR( 'ERROR: incompatible vdf shapes in os-boosted-diag' )
        ERROR( 'x_dim_: ', vdf1%x_dim_, vdf2%x_dim_ )
        ERROR( 'nx_:', vdf1%nx_(1:x_dim), vdf2%nx_(1:x_dim) )
        ERROR( 'gc_num_:', vdf1%gc_num_(:,1:x_dim), vdf2%gc_num_(:,1:x_dim) )
        call abort_program(p_err_invalid)
    endif

end subroutine check_vdf_compatability
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine reshape_array_vdf( vdf, nx_new )

    implicit none

    type(t_vdf), intent(inout), pointer  :: vdf
    integer, dimension(3), intent(in) :: nx_new

    type(t_vdf) :: vdf_tmp
    integer, dimension(3) :: nx_old

    nx_old(:) = vdf%nx_(:)

    call vdf_tmp % new( vdf, copy = .true. )

    call vdf % cleanup()
    call vdf % new( vdf_tmp%x_dim_, vdf_tmp%f_dim_, nx_new, vdf_tmp%gc_num_, vdf_tmp%dx_, .true. )

    ! copy old data into the fresh vdf
    select case( vdf%x_dim_ )
    case(1)
        vdf%f1( :, 1:nx_old(1) ) = vdf_tmp%f1( :, 1:nx_old(1) )
    case(2)
        vdf%f2( :, 1:nx_old(1), 1:nx_old(2) ) = vdf_tmp%f2( :, 1:nx_old(1), 1:nx_old(2) )
    case(3)
        vdf%f3( :, 1:nx_old(1), 1:nx_old(2), 1:nx_old(3) ) = vdf_tmp%f3( :, 1:nx_old(1), 1:nx_old(2), 1:nx_old(3) )
    end select

end subroutine reshape_array_vdf
!-------------------------------------------------------------------------------------------------




end module m_boosted_diag
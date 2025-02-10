!> A non-functioning template of the public structures provided through MARBL interface
module marbl_interface_public_types

    use marbl_logging, only : marbl_log_type

    implicit none
    private ! Only want a few types to be public

    !> A non-functioning template of MARBL diagnostic type
    type :: marbl_single_diagnostic_type
        character(len=0)                  :: long_name  !< dummy name
        character(len=0)                  :: short_name  !< dummy name
        character(len=0)                  :: units  !< dummy units
        character(len=0)                  :: vertical_grid  !< dummy grid
        logical                           :: compute_now  !< dummy flag
        logical                           :: ltruncated_vertical_extent  !< dummy flag
        integer                           :: ref_depth  !< dummy depth
        real, allocatable, dimension(:)   :: field_2d  !< dummy field
        real, allocatable, dimension(:,:) :: field_3d  !< dummy field
    end type marbl_single_diagnostic_type

    !> A non-functioning template of MARBL diagnostic type
    type, public :: marbl_diagnostics_type
        type(marbl_single_diagnostic_type), dimension(:), pointer :: diags => NULL()  !< dummy point
    end type marbl_diagnostics_type

    !> A non-functioning template of MARBL saved state type
    type :: marbl_single_saved_state_type
        integer :: rank  !< dummy rank
        character(len=0) :: short_name  !< dummy name
        character(len=0) :: units  !< dummy units
        character(len=0) :: vertical_grid  !< dummy grid
        real, allocatable :: field_2d(:)  !< dummy field
        real, allocatable :: field_3d(:,:)  !< dummy field
    end type marbl_single_saved_state_type

    !> A non-functioning template of MARBL saved state type
    type, public :: marbl_saved_state_type
        integer :: saved_state_cnt  !< dummy counter
        type(marbl_single_saved_state_type), dimension(:), pointer :: state => NULL()  !< dummy pointer
    end type marbl_saved_state_type

    !> A non-functioning template of MARBL forcing metadata type
    type :: marbl_forcing_fields_metadata_type
        character(len=0) :: varname  !< dummy name
    end type marbl_forcing_fields_metadata_type

    !> A non-functioning template of MARBL forcing type
    type, public :: marbl_forcing_fields_type
        type(marbl_forcing_fields_metadata_type) :: metadata  !< dummy metadata
        real, pointer :: field_0d(:)   => NULL()  !< dummy pointer
        real, pointer :: field_1d(:,:) => NULL()  !< dummy pointer
    end type marbl_forcing_fields_type

    !> A non-functioning template of MARBL tracer metadata type
    type, public :: marbl_tracer_metadata_type
        character(len=0) :: short_name  !< dummy name
        character(len=0) :: long_name  !< dummy name
        character(len=0) :: units  !< dummy units
    end type marbl_tracer_metadata_type

    !> A non-functioning template of MARBL domain type
    type, public :: marbl_domain_type
        integer :: kmt  !< dummy index
        integer :: km   !< dummy index
        real, allocatable :: zt(:)  !< dummy depths
        real, allocatable :: zw(:)  !< dummy depths
        real, allocatable :: delta_z(:)  !< dummy thicknesses
    end type marbl_domain_type

    !> A non-functioning template of MARBL single output type
    type, public :: marbl_single_output_type
        ! marbl_single_output :
        ! a private type, this contains both the metadata and
        ! the actual data for a single field computed in either
        ! surface_flux_compute() or interior_tendency_compute()
        ! that needs to be passed to the GCM / flux coupler.
        ! Data must be accessed via the marbl_output_for_GCM_type
        ! data structure.
        character(len=0) :: short_name  !< dummy name
        real, allocatable, dimension(:)   :: forcing_field_0d  !< dummy forcing_field_0d
        real, allocatable, dimension(:,:) :: forcing_field_1d  !< forcing_field_1d
    end type marbl_single_output_type

    !> A non-functioning template of MARBL output for GCM type
    type, public :: marbl_output_for_GCM_type
        type(marbl_single_output_type), dimension(:), pointer :: outputs_for_GCM => NULL()  !< dummy outputs_for_GCM
    end type marbl_output_for_GCM_type

end module marbl_interface_public_types
!> Contains wrapper routines that call the ufs tracing routines
module mom_cap_profiling

#ifdef UFS_TRACING
  use ufs_trace_mod, only: ufs_trace_init, ufs_trace, ufs_trace_finalize
#endif

  implicit none

  private

  public cap_profiling_init
  public cap_profiling
  public cap_profiling_finalize

contains

!> Wrapper routine that calls ufs_trace_init
  subroutine cap_profiling_init()
#ifdef UFS_TRACING
    call ufs_trace_init()
#endif
    return
  end subroutine cap_profiling_init

!> Wrapper routine that calls ufs_trace
  subroutine cap_profiling(component, routine, ph)
    character(len=*), intent(in) :: component !< Name of the component, 'mom'
    character(len=*), intent(in) :: routine   !< Name of the profiled subroutine
    character(len=*), intent(in) :: ph        !< Duration event phase type. 'B' or 'E' for begin/end
#ifdef UFS_TRACING
    call ufs_trace(component, routine, ph)
#endif
    return
  end subroutine cap_profiling

!> Wrapper routine that calls ufs_trace_finalize
  subroutine cap_profiling_finalize()
#ifdef UFS_TRACING
    call ufs_trace_finalize()
#endif
    return
  end subroutine cap_profiling_finalize

end module mom_cap_profiling

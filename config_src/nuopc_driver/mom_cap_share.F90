module mom_cap_share
  ! Temporary module for sharing ccp defs and other settings
  ! betwen NEMS and CMEPS

  use ESMF         , only: ESMF_GeomType_Flag
  use ESMF         , only: ESMF_GEOMTYPE_MESH, ESMF_GEOMTYPE_GRID
#ifdef CESMCOUPLED
  use shr_file_mod , only: shr_file_setLogUnit, shr_file_getLogUnit
#endif

  implicit none
  public

  integer :: shrlogUnit

#ifdef CESMCOUPLED
  logical :: cesm_coupled = .true.
  type(ESMF_GeomType_Flag) :: geomtype = ESMF_GEOMTYPE_MESH
#else
  logical :: cesm_coupled = .false.
  type(ESMF_GeomType_Flag) :: geomtype = ESMF_GEOMTYPE_GRID
#endif

contains

#ifndef CESMCOUPLED
  subroutine shr_file_setLogUnit(nunit)
    integer, intent(in) :: nunit
    ! do nothing for this stub - its just here to replace
    ! having cppdefs in the main program
  end subroutine shr_file_setLogUnit
#endif

end module mom_cap_share

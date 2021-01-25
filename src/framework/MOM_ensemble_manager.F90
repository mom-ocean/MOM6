!> A simple (very thin) wrapper for managing ensemble member layout information
module MOM_ensemble_manager

! This file is part of MOM6. See LICENSE.md for the license.

use ensemble_manager_mod, only : get_ensemble_id, get_ensemble_size
use ensemble_manager_mod, only : get_ensemble_pelist, get_ensemble_filter_pelist

implicit none ; private

public get_ensemble_id, get_ensemble_size, get_ensemble_pelist, get_ensemble_filter_pelist


end module MOM_ensemble_manager

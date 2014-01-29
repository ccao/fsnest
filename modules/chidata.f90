MODULE chidata
  !
  use constants
  !
  implicit none
  !
  complex(dp), allocatable :: chi_loc(:, :, :, :)
  complex(dp), allocatable :: chi_bare(:, :, :, :)
  complex(dp), allocatable :: chi_rpa(:, :, :, :)
  complex(dp), allocatable :: u_mat(:, :, :, :)
  complex(dp), allocatable :: chi_tmp(:, :, :, :)
  complex(dp), allocatable :: u_chi(:, :, :, :)
  !
CONTAINS
  !
  SUBROUTINE finalize_chi
    !
    implicit none
    !
    if(allocated(chi_loc)) deallocate(chi_loc)
    if(allocated(chi_bare)) deallocate(chi_bare)
    if(allocated(chi_rpa)) deallocate(chi_rpa)
    if(allocated(u_mat)) deallocate(u_mat)
    if(allocated(chi_tmp)) deallocate(chi_tmp)
    if(allocated(u_chi)) deallocate(u_chi)
  END SUBROUTINE
  !
END MODULE


MODULE banddata
  !
  use constants
  !
  implicit none
  !
  real(dp) ef
  real(dp), allocatable :: kvec(:, :)
  real(dp), allocatable :: fs(:)   ! # of fermi surface sheets at ik
  integer nbnd
  integer nkx, nky, nkz
  integer nkpt  ! = nkx*nky*nkz = dimension of kvec
  !
CONTAINS

SUBROUTINE init_band
  !
  implicit none
  !
  integer ikx, iky, ikz, ik
  !
  allocate(fs(1:nkpt))
  allocate(kvec(1:3, 1:nkpt))
  !
  fs(:)=0
  !
  do ikx=1, nkx
    do iky=1, nky
      do ikz=1, nkz
        ik=(ikz-1)*nkx*nky+(iky-1)*nkx+ikx
        kvec(1, ik)=(ikx-1.d0)/nkx
        kvec(2, ik)=(iky-1.d0)/nky
        kvec(3, ik)=(ikz-1.d0)/nkz
      enddo
    enddo
  enddo
  !
END SUBROUTINE

SUBROUTINE finalize_band()
  !
  implicit none
  !
  if(allocated(fs)) deallocate(fs)
  if(allocated(kvec)) deallocate(kvec)
  !
END SUBROUTINE

END MODULE


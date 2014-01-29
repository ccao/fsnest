MODULE para
  !
#if defined __MPI
  use mpi
#endif
  !
  implicit none
  !
  integer inode, nnode
  integer first_k, last_k
  !
INTERFACE para_sync
  MODULE PROCEDURE para_sync_int0, para_sync_real0, para_sync_real1, para_sync_real2, para_sync_cmplx3
END INTERFACE
  !
INTERFACE para_merge
  MODULE PROCEDURE para_merge_int1, para_merge_real0, para_merge_real1, para_merge_real2, para_merge_cmplx0, para_merge_cmplx1, para_merge_cmplx3, para_merge_cmplx4
END INTERFACE
  !
CONTAINS
  !
SUBROUTINE distribute_k()
  !
  use banddata, only : nkpt
  !
  implicit none
  !
#if defined __MPI
  first_k=inode*nkpt/nnode+1
  last_k=(inode+1)*nkpt/nnode
#else
  first_k=1
  last_k=nkpt
#endif
  !
END SUBROUTINE

SUBROUTINE init_para()
  !
  use constants
  !
  implicit none
  !
  integer ierr
  !
#if defined __MPI
  CALL mpi_init(ierr)
  !
  CALL mpi_comm_rank(mpi_comm_world, inode, ierr)
  CALL mpi_comm_size(mpi_comm_world, nnode, ierr)
  !
  if (inode.eq.0) write(stdout, *) "# WannChi running on ", nnode, " nodes..."
#else
  inode=0
  nnode=1
  write(stdout, *) "# WannChi serial ..."
#endif
  !
END SUBROUTINE

SUBROUTINE para_sync_int0(dat)
  !
  implicit none
  !
  integer :: dat
  !
  integer ierr
  !
#if defined __MPI
  CALL mpi_bcast(dat, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE

SUBROUTINE para_sync_real0(dat)
  !
  use constants, only: dp
  !
  implicit none
  !
  real(dp) :: dat
  !
  integer ierr
  !
#if defined __MPI
  CALL mpi_bcast(dat, 1, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE

SUBROUTINE para_sync_real1(dat, dat_size)
  !
  use constants, only : dp
  !
  implicit none
  !
  integer :: dat_size
  real(dp) :: dat(1:dat_size)
  !
  integer ierr
  !
#if defined __MPI
  CALL mpi_bcast(dat, dat_size, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE

SUBROUTINE para_sync_real2(dat, size1, size2)
  !
  use constants, only: dp
  !
  implicit none
  !
  integer :: size1, size2
  real(dp) :: dat(1:size1, 1:size2)
  !
  integer ierr, dat_size
  !
#if defined __MPI
  dat_size=size1*size2
  CALL mpi_bcast(dat, dat_size, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE

SUBROUTINE para_sync_cmplx3(dat, size1, size2, size3)
  !
  use constants, only : dp
  !
  implicit none
  !
  integer :: size1, size2, size3
  complex(dp) :: dat(1:size1, 1:size2, 1:size3)
  !
  integer ierr, dat_size
  !
#if defined __MPI
  dat_size=size1*size2*size3
  CALL mpi_bcast(dat, dat_size, MPI_DOUBLE_COMPLEX, 0, mpi_comm_world, ierr)
#endif
  !
END SUBROUTINE

SUBROUTINE para_merge_int1(dat, dat_size)
  !
  implicit none
  !
  integer :: dat_size
  integer :: dat(1:dat_size)
  integer, allocatable :: buf(:)
  !
  integer ierr
  !
#if defined __MPI
  allocate(buf(1:dat_size))
  CALL mpi_allreduce(dat, buf, dat_size, MPI_INTEGER, MPI_SUM, mpi_comm_world, ierr)
  dat(:)=buf(:)
  deallocate(buf)
#endif
  !
END SUBROUTINE

SUBROUTINE para_merge_real0(dat)
  !
  use constants, only :dp
  !
  implicit none
  !
  real(dp) dat, tot
  !
  integer ierr
  !
#if defined __MPI
  CALL mpi_allreduce(dat, tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
  dat=tot
#endif
  !
END SUBROUTINE

SUBROUTINE para_merge_real1(dat, dat_size)
  !
  use constants, only :dp
  !
  implicit none
  !
  integer :: dat_size
  real(dp) :: dat(1:dat_size)
  real(dp), allocatable :: buf(:)
  !
  integer ierr
  !
#if defined __MPI
  allocate(buf(1:dat_size))
  CALL mpi_allreduce(dat, buf, dat_size, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
  dat(:)=buf(:)
  deallocate(buf)
#endif
  !
END SUBROUTINE

SUBROUTINE para_merge_real2(dat, size1, size2)
  !
  use constants, only :dp
  !
  implicit none
  !
  integer :: size1, size2
  real(dp) :: dat(1:size1, size2)
  !
  integer ierr
  !
  real(dp), allocatable :: buf(:, :)
  !
#if defined __MPI
  allocate(buf(1:size1, 1:size2))
  CALL mpi_allreduce(dat, buf, size1*size2, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
  dat(:,:)=buf(:,:)
  deallocate(buf)
#endif
  !
END SUBROUTINE

SUBROUTINE para_merge_cmplx0(dat)
  !
  use constants, only : dp
  !
  implicit none
  !
  complex(dp) :: dat, tot
  integer ierr
  !
#if defined __MPI
  CALL mpi_allreduce(dat, tot, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_world, ierr)
  dat=tot
#endif
  !
END SUBROUTINE

SUBROUTINE para_merge_cmplx1(dat, dat_size)
  !
  use constants, only : dp
  !
  implicit none
  !
  integer :: dat_size
  complex(dp) :: dat(1:dat_size)
  !
  integer ierr
  !
  complex(dp), allocatable :: buf(:)
  !
#if defined __MPI
  allocate(buf(1:dat_size))
  CALL mpi_allreduce(dat, buf, dat_size, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_world, ierr)
  dat(:)=buf(:)
  deallocate(buf)
#endif
  !
END SUBROUTINE

SUBROUTINE para_merge_cmplx3(dat, size1, size2, size3)
  !
  use constants, only : dp
  !
  implicit none
  !
  integer :: size1, size2, size3
  complex(dp) :: dat(1:size1, 1:size2, 1:size3)
  !
  complex(dp), allocatable :: buf(:, :, :)
  !
  integer ierr
  !
#if defined __MPI
  allocate(buf(1:size1, 1:size2, 1:size3))
  CALL mpi_allreduce(dat, buf, size1*size2*size3, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_world, ierr)
  dat(:,:,:)=buf(:,:,:)
  deallocate(buf)
#endif
  !
END SUBROUTINE

SUBROUTINE para_merge_cmplx4(dat, size1, size2, size3, size4)
  !
  use constants, only: dp
  implicit none
  !
  integer :: size1, size2, size3, size4
  complex(dp) :: dat(:, :, :, :)
  integer ierr
  !
  complex(dp), allocatable :: buf(:, :, :, :)
  !
#if defined __MPI
  allocate(buf(1:size1, 1:size2, 1:size3, 1:size4))
  CALL mpi_allreduce(dat, buf, size1*size2*size3*size4, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_world, ierr)
  dat(:, :, :, :)=buf(:, :, :, :)
  deallocate(buf)
#endif
  !
END SUBROUTINE

SUBROUTINE finalize_para
  !
  implicit none
  !
#if defined __MPI
  integer ierr
  CALL mpi_finalize(ierr)
#endif
  !
END SUBROUTINE

END MODULE

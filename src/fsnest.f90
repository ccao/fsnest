PROGRAM wannchi
  !
  use constants,only : cmplx_0, stdout, dp, fout
  use para,     only : init_para, inode, distribute_k, para_merge, finalize_para
  use wanndata, only : read_ham, norb, finalize_wann
  use banddata, only : nkpt, nbnd, init_band, finalize_band
  use input,    only : read_input, seed, qvec, nqpt
  !
  implicit none
  !
  integer iq
  integer dnq
  real(dp), allocatable :: fsn(:)
  !
  CALL init_para
  CALL read_input
  CALL read_ham(seed)
  !
  nbnd=norb
  !
  CALL init_band
  !
  CALL distribute_k
  !
  CALL interpolate_bands
  !
  allocate(fsn(1:nqpt))
  fsn(:)=0.d0
  !
  open(unit=fout, file="fsnest.dat")
  !
  dnq=nqpt/100
  if (dnq<2) dnq=2
  !
  do iq=1, nqpt
    !
    CALL compute_fsnest(fsn(iq), qvec(:, iq))
    !
    if ((nqpt.ne.1).and.(mod(iq-1, dnq).eq.0)) then
      if (inode.eq.0) then
        write(stdout, *) " #... Percentage done: ", (iq-1)*100/nqpt, "%"
      endif
    endif
    !
    CALL output_fsnest(qvec(:, iq), fsn(iq))
    !
  enddo
  !
  deallocate(fsn)
  !
  close(unit=fout)
  !
  CALL finalize_wann
  !
  CALL finalize_band
  !
  CALL finalize_para
  !
END PROGRAM

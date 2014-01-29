SUBROUTINE compute_fsnest(fsn, qv) 
  !
  use constants, only : dp, eps9, cmplx_0, cmplx_i, stdout
  use para,      only : first_k, last_k, para_merge, inode
  use banddata,  only : kvec, nbnd, fs, nkpt
  use input,     only : omega, eps
  !
  implicit none
  !
  real(dp) qv(1:3)
  real(dp) fsn
  !
  real(dp) ikv(1:3), jkv(1:3)
  !
  integer ik, jk
  !
  integer kpt_index
  !
  do ik=first_k, last_k
    !
    if (fs(ik).ne.0) then
      !
      ikv(:)=kvec(:, ik)
      jkv(:)=ikv(:)+qv(:)
      jk=kpt_index(jkv)
      !
      fsn=fsn+fs(ik)*fs(jk)
      !
    endif
    !
  enddo ! ik
  !
  fsn=fsn/nkpt
  !
  CALL para_merge(fsn)
  !
END SUBROUTINE


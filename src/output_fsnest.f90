SUBROUTINE output_fsnest(qv, fsnest)
  !
  use constants, only : dp, fout
  use para,      only : inode
  !
  implicit none
  !
  real(dp) qv(1:3)
  real(dp) fsnest
  !
  if (inode.eq.0) then
    write(fout, '(3F12.8,1F22.16)') qv, fsnest
  endif
  !
END SUBROUTINE

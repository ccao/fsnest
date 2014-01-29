INTEGER FUNCTION kpt_index(kv)
  !
  use banddata, only : nkx, nky, nkz
  use constants, only : dp, eps6
  !
  implicit none
  !
  real(dp) kv(1:3)
  integer ik, ikx, iky, ikz
  !
  ikx=nint((kv(1)-floor(kv(1)))*nkx)
  iky=nint((kv(2)-floor(kv(2)))*nky)
  ikz=nint((kv(3)-floor(kv(3)))*nkz)
  !
  kpt_index=ikz*nkx*nky+iky*nkx+ikx+1
  return
  !
END FUNCTION

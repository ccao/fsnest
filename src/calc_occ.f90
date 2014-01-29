REAL(dp) FUNCTION calc_occ(en)
  !
  use constants, only : dp, eps6
  use input, only : temp
  !
  implicit none
  !
  real(dp) en
  !
  if ( temp < eps6 ) then ! Zero temperature
    if ( en > 0.d0 ) then
      calc_occ=0.d0
    else if ( en < 0.d0 ) then
      calc_occ=1.d0
    else
      calc_occ=0.5d0
    endif
  else  ! Fermi-Dirac distribution for non-zero temperature
    calc_occ=1.d0/(1.d0+exp(en/temp))
  endif
  !
END FUNCTION

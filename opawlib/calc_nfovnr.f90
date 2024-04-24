subroutine calc_nfovnr
  use opaw_mod, only : dx,dy,dz,p_fg
  use paw_mod, only : nfovnr
  implicit none

  nfovnr=max(ceiling(dx/p_fg),1)
  nfovnr=max(nfovnr,ceiling(dy/p_fg))
  nfovnr=max(nfovnr,ceiling(dz/p_fg))

end subroutine

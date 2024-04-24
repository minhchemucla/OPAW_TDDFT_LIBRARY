subroutine calc_eig(ham)
  use main_mod
  use mpi_lib_ours
  use ham_mod
  implicit none
  integer :: i
  real*8  :: denom, numer, tmp(nn)
  complex*16 :: p(nn), hp(nn)
  type(hamiltonian_obj) :: ham

  
  if(rank==0) then
    if(h_type .eq. 0) then
      write(*,*) 'eigs, <p|s^(-1)h|p>/<p|p>'
      open(unit=15,file='wf.txt')
    else if(h_type .eq. 1) then
      write(*,*) 'eigs, <p|s^(-1/2)hs^(-1/2)|p>/<p|p>'
      open(unit=15,file='wf_bar.txt')
    endif
    do i=1,11;read(15,*);enddo
    do i=1,nocc
      read(15,*);read(15,*) tmp
      p = tmp
      if(h_type .eq. 0) then
        call sh(ham,p,hp)
      else
        call shs(ham,p,hp)
      endif
      numer = sum(conjg(p)*hp)*dv
      denom = sum(conjg(p)*p)*dv
      write(*,*) eigs(i), numer/denom
    enddo
  endif


end subroutine calc_eig

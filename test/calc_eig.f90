subroutine calc_eig(ham)
  use main_mod
  use mpi_lib_ours
  use opaw_ham_mod
  implicit none
  integer :: i
  real*8  :: denom, numer, tmp(nn)
  complex*16 :: p(nn), hp(nn)
  type(opaw_ham_obj) :: ham

  
  if(rank==0) then
    write(*,*) 'eigs, <p|s^(-1/2)hs^(-1/2)|p>/<p|p>'
    open(unit=15,file='wf_bar.txt')
    do i=1,11;read(15,*);enddo
    do i=1,nocc
      read(15,*);read(15,*) tmp
      p = tmp
      call shs(ham,p,hp)
      numer = sum(conjg(p)*hp)*dv
      denom = sum(conjg(p)*p)*dv
      write(*,*) eigs(i), numer/denom
    enddo
  endif


end subroutine calc_eig

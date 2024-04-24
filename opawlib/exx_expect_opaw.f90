! Label the psuedo wavefunctions by |psi_i>
subroutine exx_expect_opaw(nn,nocc,nstates,psi_i,psi_j,psi_n,exx)
  use opaw_mod, only : scale_vh,dv
  use atom_mod
  use atom_mod, only : p=>pawinfo, at=> atominfo
  use libpaw_mod 
  use paw_mod
  use mpi_lib_ours
  implicit none
  integer :: nn,nocc,nstates
  integer :: i,st
  real*8  :: exx !expectation value of PAW exact exchange operator
  complex*16  :: psi_i(nn), psi_j(nn), psi_n(nn,nstates)
  real*8, allocatable :: rh(:), u(:), pot(:)

  if(rank==0) then
    !write(*,*) 'sum(abs(psi_j)), sum(abs(psi_i)), sum(abs(psi_n))'
    !write(*,*) sum(abs(psi_j)), sum(abs(psi_i)), sum(abs(psi_n))
    exx = 0d0
    allocate(u(nn),pot(nn),rh(nn),stat=st); if(st/=0) stop 'rh,u'
    do i=1,nocc
      rh = psi_j*psi_n(:,i)
      call vh_sub_opaw(rh,pot,scale_vh)
      rh = psi_i*psi_n(:,i)
      exx = exx - sum(rh*pot)*dv
    enddo
    write(6,*) 'i, exx: ', nocc, exx
  endif
end subroutine exx_expect_opaw

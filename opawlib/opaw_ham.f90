!apply S^-1/2 H S^-1/2 to a wavefunction
!S^-1/2 (K+Vloc) S^-1/2 + S^-1/2 D S^-1/2
subroutine opaw_ham_c16(ham,pin,pout)
    use mpi_lib_ours
    use opaw_mod, only : nx,ny,nz
    use paw_mod
    use atom_mod
    use atom_mod, only : p=> pawinfo, at=> atominfo
    use opaw_ham_mod
    implicit none

    type(opaw_ham_obj) :: ham
    complex*16,intent(in) :: pin(nx,ny,nz)
    complex*16 :: pout(nx,ny,nz)
    complex*16 :: tmp(nx,ny,nz)

    integer :: ib,jb,ia,is,js,it,ik
    integer :: ix,iy,iz,ix1,iy1,iz1,jx,jy,jz
    integer :: ixs,iys,izs,ixe,iye,ize

    call sn_phi(pin,tmp,-0.5d0)
    call paw_ham(ham,tmp,pout)
    tmp=pout
    call sn_phi(tmp,pout,-0.5d0)
end subroutine opaw_ham_c16

!I DO NOT RECOMMEND MAKING A _r8 version as the 
!input and out of sn_phi will have a complex part

!apply S^-1/2 H S^-1/2 to a wavefunction
!S^-1/2 (K+Vloc) S^-1/2 + S^-1/2 D S^-1/2
subroutine sh(ham,pin,pout)
!make sure the input is a complex and not a real because the
!implicit conversion from real to complex is buggy sometimes
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
    complex*16, allocatable :: ca(:)


    ik=1
    call h_phi(ham,pin,pout)
!pout= s^-1 h pin
    call time_print('after apply h')
    do ia=1,natom
        it=atom_map(ia)
        allocate(ca(p(it)%mstates),stat=stat)
        if(stat/=0) stop 'ca alloc problem in sphi'
        ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2)
        izs=at(ia)%ir_start(3)
        ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
        ize=izs+p(it)%nrough(3)-1
        call proj1(ia,pout,ca,p(it)%mstates,ik)
!        call proj(ia,pout,ca,p(it)%mstates)

!write(*,*) 'ia, ca', is, sum(abs(ca)), sum(abs(at(ia)%local_p3d1_c)), sum(at(ia)%sinv)
        if(.not.at(ia)%edge) then
            do is=1,p(it)%mstates
               pout(ixs:ixe,iys:iye,izs:ize)=&
                pout(ixs:ixe,iys:iye,izs:ize)+&
                 at(ia)%local_p3d1_c(:,:,:,is,ik)*at(ia)%sinv(is,ik)*&
!                at(ia)%local_p3d(:,:,:,is)*p(it)%ssqinvij(is)*&
                ca(is)
            enddo
        else
            do iz=izs,ize;do iy=iys,iye;do ix=ixs,ixe
                jx=mod(ix+nx-1,nx)+1
                jy=mod(iy+ny-1,ny)+1
                jz=mod(iz+nz-1,nz)+1
                pout(jx,jy,jz)=pout(jx,jy,jz)+&
                   sum(at(ia)%local_p3d1_c(ix-ixs+1,iy-iys+1,iz-izs+1,:,ik)*&
                   at(ia)%sinv(:,ik)*ca)
            enddo;enddo;enddo
        endif
        deallocate(ca)
    enddo
!    write(*,*) 'sum(pout) s', sum(pout)
end subroutine sh

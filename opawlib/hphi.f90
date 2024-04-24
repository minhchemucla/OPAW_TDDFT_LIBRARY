subroutine h_phi(ham,pin,hp)
!make sure the input is a complex and not a real because the
!implicit conversion from real to complex is buggy sometimes
    use opaw_mod, only : nx,ny,nz, nn, ek3d
    use atom_mod
    use atom_mod, only : p=> pawinfo, at_common=> atominfo
    use ham_mod
    use mpi_lib_ours
    implicit none

    integer :: ik
    complex*16 :: pin(nx,ny,nz)
    complex*16 :: hp(nx,ny,nz),tmp(nx,ny,nz)
    complex*16 :: cin(nx,ny,nz),cout(nx,ny,nz)
    type(hamiltonian_obj) :: ham
    !real*8     :: wtime

    ik = 1
    
    cin=pin
    call fft3d_forward(nx,ny,nz,cin,cout)
    cout=cout*ek3d(:,:,:)
    call fft3d_backward(nx,ny,nz,cout,tmp)
    hp=tmp
    tmp = reshape(ham%vks,(/nx,ny,nz/))
    hp=hp+pin*tmp
    call addvnl
contains
    subroutine addvnl
        implicit none

        integer :: ia,ig,ix,iy,iz,it,is,js,ix1,iy1,iz1
        integer :: ixs,iys,izs,ixe,iye,ize,jx,jy,jz
        complex*16, allocatable:: ca(:),arr(:)

        do ia=1,natom
            it=atom_map(ia)
            ixs=at_common(ia)%ir_start(1);iys=at_common(ia)%ir_start(2)
            izs=at_common(ia)%ir_start(3)
            ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
            ize=izs+p(it)%nrough(3)-1
            allocate(ca(p(it)%mstates),arr(p(it)%mstates),stat=stat)
            if(stat/=0) stop 'ca alloc problem in addvnl'
!            call time_print('before proj')
            call proj(ia,pin,ca,p(it)%mstates,ik)
!            call time_print('after proj')
            if(.not.at_common(ia)%edge) then
                do is=1,p(it)%mstates
                 hp(ixs:ixe,iys:iye,izs:ize)=&
                     hp(ixs:ixe,iys:iye,izs:ize)+&
                     at_common(ia)%local_p3d_c(:,:,:,is,ik)*&
                     sum(ham%at(ia)%dij(is,:)*ca(:))
                enddo
            else
                arr=matmul(ham%at(ia)%dij,ca)
                do iz=izs,ize;do iy=iys,iye;do ix=ixs,ixe
                    jx=mod(ix+nx-1,nx)+1
                    jy=mod(iy+ny-1,ny)+1
                    jz=mod(iz+nz-1,nz)+1
                      do is=1,p(it)%mstates
                        hp(jx,jy,jz)=hp(jx,jy,jz)+&
                        at_common(ia)%local_p3d_c(ix-ixs+1,iy-iys+1,iz-izs+1,is,ik)*&
                        arr(is) 
                      enddo
                enddo;enddo;enddo
            endif
!            call time_print('after add dij')
            deallocate(ca,arr)
        enddo
    end subroutine
end subroutine h_phi

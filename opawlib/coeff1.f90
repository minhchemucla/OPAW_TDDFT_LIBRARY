subroutine acc_coeff1(nx,ny,nz,nocc,wf,ham)
    use paw_mod
    use ham_mod
    use atom_mod, only : p => pawinfo, atom_map, ngrid, natom
    use mpi_lib_ours
    implicit none

    integer :: ia,ig,ix,iy,iz,is,ms,it,js,ib,iib
    integer :: nx,ny,nz,nocc
    complex*16  :: wf(nx,ny,nz,nocc)
    complex*16, allocatable :: ca(:),tmp(:,:)
    type(hamiltonian_obj) :: ham

    do ia=1,natom
      ham%at(ia)%rhoij=0d0
    enddo

    do ia=1,natom
      it=atom_map(ia)
      ms=p(it)%mstates
      allocate(ca(ms),tmp(ms,ms),stat=stat)
      if(stat/=0) stop 'ca alloc problem in coeff'
      do ib=1,nocc
        tmp=0d0
        call proj(ia,wf(:,:,:,ib),ca,ms,1) 
        do is=1,ms
          do js=1,ms
            tmp(is,js)=tmp(is,js)+conjg(ca(is))*ca(js)*2d0 !Note
          enddo
        enddo
        !write(*,*) 'ia,ib,ca',ia,ib,sum(ca)
        ham%at(ia)%rhoij=ham%at(ia)%rhoij+dble(tmp)
      enddo
      deallocate(ca,tmp)
    enddo
   ! enddo
    !do ia=1,natom
    !  write(*,*) 'ia, rhoij', ia, sum(ham%at(ia)%rhoij)
    !enddo
end subroutine acc_coeff1


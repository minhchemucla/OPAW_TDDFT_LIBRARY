subroutine get_dij(ham)
    use mpi_lib_ours
    use opaw_mod 
    use libpaw_mod
    use atom_mod
    use atom_mod, only : p => pawinfo, at => atominfo
    use paw_mod
    use opaw_ham_mod
    implicit none
    real*8  :: param_dij_max=500d0
    real*8, allocatable :: vxc2d(:,:)
    integer :: ia   
    type(opaw_ham_obj) :: ham

    allocate(vxc2d(nn,1))
    if(rank==0) then
      call set_rhoij(ham)

      vxc2d(:,1) = ham%vxc
      call pawdenpot(compch_sph,epaw,epawdc,0,ixc,&
          natom,natom,nspden,ntypat,nucdipmom,0,0,paw_an,paw_an,&
          paw_ij,pawang,pawprtvol,pawrad,pawrhoij,0,pawtab,xcdev,1.0d0,xclevel,1d-14,ucvol,znucl)

      call pawdij(cplex,enunit,gprimd,ipert,natom,natom,nfft,nfftot,&
        &     nspden,ntypat,paw_an,paw_ij,pawang,pawfgrtab,pawprtvol,&
        &     pawrad,pawrhoij,pawspnorb,pawtab,xcdev,qphon,spnorbscl,&
        &     ucvol,charge,ham%vks,vxc2d,xred)

      call transfer_dij
    endif


    do ia=1,natom
      call bcast_r8(ham%at(ia)%dij,size(at(ia)%dij),0)
 !     if(rank==0) write(*,*) 'ia, dij', ia, sum(at(ia)%dij)
    enddo
    deallocate(vxc2d)
contains
    subroutine transfer_dij
        implicit none

        integer :: i,j,ind,ndim,ndim1
        real*8 dd ! DN

        siz_dijall=0
        do ia=1,natom
            ndim=size(ham%at(ia)%rhoij,1)
            ndim1=size(paw_ij(ia)%dij,1)
            if (ndim*(ndim+1)/2/=ndim1) then
                write(*,*) 'dim of dij not match,ia,ndim,ndim1',ia,ndim,ndim1
                stop
            endif
            siz_dijall=siz_dijall+ndim1
        enddo

        if(allocated(dijall)) deallocate(dijall)
        if(.not.allocated(dijall)) then
            allocate(dijall(siz_dijall),stat=i)
            if(i/=0) stop 'dijall alloc problem'
        endif
       
        !write(*,*) 'siz_dijall',siz_dijall 

        siz_dijall=1
        do ia=1,natom
!            at(ia)%dij(i,j)=min(paw_ij(ia)%dij(ind,1),100d0)
!            at(ia)%dij(j,i)=min(paw_ij(ia)%dij(ind,1),100d0)
            ndim1=size(paw_ij(ia)%dij,1)
            dijall(siz_dijall:siz_dijall+ndim1-1)=paw_ij(ia)%dij(:,1)
            siz_dijall=siz_dijall+ndim1
        enddo

        if (siz_dijall/=size(dijall)+1) stop 'siz dij problem' 

        siz_dijall=1
        do ia=1,natom
            ndim=size(ham%at(ia)%rhoij,1)
            ndim1=size(paw_ij(ia)%dij,1)
            do j=1,ndim
                do i=j,ndim
                   ind=i*(i-1)/2+j
                   dd = dijall(siz_dijall+ind-1) !DN
                   if(abs(dd)>param_dij_max) dd = dd /abs(dd)*param_dij_max !DN 
                   ham%at(ia)%dij(i,j)=dd ! DN 
                   ham%at(ia)%dij(j,i)=dd ! DN
                   !if(i==ndim.and.j==ndim) then
                   !   write(6,*)' for atom ',ia,' dij_18_18 debug ',dd
                   !endif
                enddo
            enddo
            siz_dijall=siz_dijall+ndim1
        enddo
    end subroutine
end subroutine get_dij

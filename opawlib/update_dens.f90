subroutine update_dens_paw(n,nocc,wf,ham)
      use opaw_mod, only : nx,ny,nz, dv
      use atom_mod, only : natom, at => atominfo, atom_map
      use opaw_ham_mod
      use mpi_lib_ours
!                at_old => atominfo_old
      implicit none

      integer :: i,ia,it
      integer :: n,nocc
      complex*16 :: wf(n,nocc)
      real*8   :: dens(n),nhat(n)
      real*8  :: n1,n2,n3
      type(opaw_ham_obj) :: ham

      if(rank==0)call acc_coeff1(nx,ny,nz,nocc,wf,ham)
      ham%dens = 0d0
      do i=1,nocc
         ham%dens=ham%dens+2d0*abs(wf(:,i))**2d0
      enddo

      if(rank==0) then 
        !write(*,*) 'before allsum,rank,dens',rank,sum(ham%dens)*dv
        call get_nhat(nx,ny,nz,natom,ham%dens,ham%nhat,ham%at)
      endif
      call bcast_r8(ham%nhat,size(ham%nhat),0) 
!      !if(rank==0) write(*,*) 'max,min(dens)',rank,maxval(dens), minval(dens)
!      if(rank==0) write(*,*) 'max,min(nhat)',rank,maxval(ham%nhat), minval(ham%nhat)
end subroutine update_dens_paw      

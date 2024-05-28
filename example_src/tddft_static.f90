subroutine tddft_static
  use mpi_lib_ours
  use ham_obj, only : ham, ham_pert
  use main_mod, only : nn,dt,wfs,nt,nocc,nstates,sm
  use main_mod, only : nx,ny,nz,dx,dy,dz,ipol,dv
  use opaw_ham_mod
  implicit none
  integer :: it,st,i
  complex*16, allocatable :: wfs_tmp(:,:)
  
  allocate(wfs_tmp(nn,nstates),stat=st); if(st/=0) stop 'wfs_tmp'
  wfs_tmp = wfs
  do it=1,nt
    call rk4_static_prop_opaw_c16(nn,nstates,dt,wfs_tmp,ham)
    call plot_dip
  enddo

  contains
    subroutine plot_dip
      implicit none
      integer :: ix,iy,iz
      real*8 :: dip(3),dip_pert(3),d_dip(3)
      real*8 :: rr(3),a,b
      complex*16 :: tmp(nn), sntmp(nn)

      dip=0d0

      ham%dens_o = 0d0
      do i=1,nocc
        ham%dens_o(:) =  ham%dens_o(:) + 2d0*abs(wfs_tmp(:,i))**2d0
      enddo

      i=0
      do iz=1,nz
         do iy=1,ny
            do ix=1,nx
               i = i+1
               rr = (/ (ix-1-nx/2)*dx, (iy-1-ny/2)*dy, (iz-1-nz/2)*dz /)
               dip    = dip   + dv * rr(:) * ham%dens_o(i)
            enddo
         enddo
      enddo

      if(rank==0) write(*,*)it,(it-1)*dt, dip,'it, t,     dip_static'; call flush(6)
    end subroutine plot_dip
end subroutine tddft_static

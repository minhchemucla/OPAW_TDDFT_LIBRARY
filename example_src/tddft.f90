subroutine tddft
  use mpi_lib_ours
  use ham_obj, only : ham, ham_pert
  use main_mod, only : nn,dt,wfs,nt,nocc,nstates,sm
  use main_mod, only : nx,ny,nz,dx,dy,dz,ipol,dv
  use opaw_ham_mod
  implicit none
  integer :: it,st,i
  complex*16, allocatable :: wfs_pert(:,:)
  
  call perturb_wf
  do it=1,nt
    call opaw_make_ham_c16(nn,nocc,nstates,wfs,ham)
    call rk4_prop_opaw_c16(nn,nocc,nstates,dt,wfs,ham)

    call opaw_make_ham_c16(nn,nocc,nstates,wfs_pert,ham_pert)
    call rk4_prop_opaw_c16(nn,nocc,nstates,dt,wfs_pert,ham_pert)
    if(rank==0) call plot_dip
  enddo

  contains
    subroutine plot_dip
      implicit none
      integer :: ix,iy,iz
      real*8 :: dip(3),dip_pert(3),d_dip(3)
      real*8 :: rr(3),a,b
      complex*16 :: tmp(nn), sntmp(nn)

      dip=0d0
      dip_pert=0d0
      d_dip = 0d0

      ham%dens_o = 0d0
      ham_pert%dens_o = 0d0
      do i=1,nocc
        ham%dens_o(:) =  ham%dens_o(:) + 2d0*abs(wfs(:,i))**2d0
        ham_pert%dens_o(:) =  ham_pert%dens_o(:) + 2d0*abs(wfs_pert(:,i))**2d0
      enddo

      i=0
      do iz=1,nz
         do iy=1,ny
            do ix=1,nx
               i = i+1
               rr = (/ (ix-1-nx/2)*dx, (iy-1-ny/2)*dy, (iz-1-nz/2)*dz /)
               dip    = dip   + dv * rr(:) * ham%dens_o(i)
               dip_pert   = dip_pert  + dv * rr(:) * ham_pert%dens_o(i)
            enddo
         enddo
      enddo

      d_dip = (dip_pert - dip)/sm
      write(*,*)it,(it-1)*dt,  d_dip(:),'it, t,     d_dip '; call flush(6)
      write(*,*)it,(it-1)*dt, dip, dip_pert,'it, t,     dip,    dip_pert '; call flush(6)
    end subroutine plot_dip

    subroutine perturb_wf
      implicit none
      integer :: ix,iy,iz
      complex*16 :: tmp3d(nx,ny,nz), tmp3d_pert(nx,ny,nz)
      complex*16 :: ci=(0d0,1d0)
      real*8 :: rr(3)

      allocate(wfs_pert(nn,nstates),stat=st); if(st/=0) stop 'wfs_pert'
      do i=1,nstates
        tmp3d = reshape(wfs(:,i),(/nx,ny,nz/))
        tmp3d_pert = 0d0
        do iz=1,nz
          do iy=1,ny
            do ix=1,nx
              rr = (/ (ix-1-nx/2)*dx, (iy-1-ny/2)*dy, (iz-1-nz/2)*dz /)
              !wfs_pert(ix,iy,iz,:) = wfs(ix,iy,iz,:)* exp(-ci*sm*rr(ipol))
              tmp3d_pert(ix,iy,iz) = tmp3d(ix,iy,iz)* exp(-ci*sm*rr(ipol))
            enddo
          enddo
        enddo
        wfs_pert(:,i) = reshape(tmp3d_pert,(/nn/))
      enddo

    end subroutine perturb_wf

end subroutine tddft

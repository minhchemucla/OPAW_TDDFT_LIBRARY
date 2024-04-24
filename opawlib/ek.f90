subroutine ek3d_prep_opaw
  use opaw_mod
  !use tddft_mod, only : dt, tddft_flag
  use paw_mod
  use mpi_lib_ours
  implicit none
  integer i
  integer ix, iy, iz, ik, jk
  real*8  kx, ky, kz, kx1,ky1,kz1
  real*8  dkx, dky, dkz
  real*8  pi
  complex*16, parameter :: ci = (0d0,1d0)

  if(.not.allocated(ek3d)) allocate(ek3d(nx,ny,nz))
  pi = dacos(-1d0)
  
  dkx = 2d0*pi/(nx*dx)
  dky = 2d0*pi/(ny*dy)
  dkz = 2d0*pi/(nz*dz)

  i=0
  do iz=1,nz
     do iy=1,ny
        do ix=1,nx
           i=i+1

           kx = (ix-1)*dkx
           ky = (iy-1)*dky
           kz = (iz-1)*dkz
           
           if(kx>pi/dx) kx = kx - 2d0*pi/dx
           if(ky>pi/dy) ky = ky - 2d0*pi/dy
           if(kz>pi/dz) kz = kz - 2d0*pi/dz

           ek3d(ix,iy,iz) = min(ekcut,(kx**2+ky**2+kz**2)/2d0)
        enddo
     enddo
  enddo
  ek3d=ek3d/(dble(nx)*dble(ny)*dble(nz))
end subroutine ek3d_prep_opaw

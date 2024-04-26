subroutine read_wfs
  use main_mod
  use mpi_lib_ours, only : rank
  implicit none
  integer :: nsp_r
  integer :: st
  integer :: nx_r, ny_r, nz_r,i
  real*8  :: dx_r, dy_r, dz_r
  character*8 ch
  real*8  :: tmp(nn)

  close(441)
  open (441,file='wf_bar.txt',status='old');  
  rewind(441)
  
  read(441,*) ch, nx_r
  read(441,*) ch, ny_r
  read(441,*) ch, nz_r
  read(441,*) ch, dx_r
  read(441,*) ch, dy_r
  read(441,*) ch, dz_r
  read(441,*) ch, nsp_r
  read(441,*) ch, nstates

  call check(nx,nx_r,' nx, nx_r')
  call check(ny,ny_r,' ny, ny_r')
  call check(nz,nz_r,' nz, nz_r')
  call check(nsp_r,1,' 1, nsp_r')
  call check_r(dx,dx_r,' dx, dx_r')
  call check_r(dy,dy_r,' dy, dy_r')
  call check_r(dz,dz_r,' dz, dz_r')

  allocate(wfs(nn,nstates),eigs(nstates), stat=st);if(st/=0) stop 'wfs'

  read(441,*) ch !eigs
  read(441,*) eigs
  read(441,*) ch !skip orbitals

  !write(*,*) eigs
  
  dens = 0d0
  do i=1,nstates
    read(441,*) 
    read(441,*) tmp
    wfs(:,i) = tmp
    !if(rank==0) write(6,*) 'is, norm', i, sum(abs(wfs(:,i))**2d0)*dv, wfs(1,i)
    if(i <= nocc) dens = dens + 2d0*abs(wfs(:,i))**2d0
  enddo
  
  if(rank==0) write(*,*) 'sum(dens)*dv', sum(dens)*dv

  close(441)

end subroutine

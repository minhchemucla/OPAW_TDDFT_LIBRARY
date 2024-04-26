subroutine read_input
  use main_mod
  use mpi_lib_ours
  implicit none
  open(unit=13,file='dim_test.inp')

  if (rank==0) then
    ! number of grid points
    call fetch_i (nx,'nx',13)
    call fetch_i (ny,'ny',13)
    call fetch_i (nz,'nz',13)

    !length of the box that encloses system
    call fetch_r (dx,'dx',13)
    call fetch_r (dy,'dy',13)
    call fetch_r (dz,'dz',13)

    call fetch_l (periodic,'periodic',13)
    call fetch_i (funct, 'funct',13)
    call fetch_l (flg_bin   , 'flg_bin',13)
    call fetch_r (ekcut, 'ekcut',13)
    !call fetch_r (p_fg, 'finegrid',13)
    call fetch_r (dt, 'dt',13)
    call fetch_i (nt, 'nt',13)
    call fetch_r (sm, 'strength',13)
    call fetch_i(ipol, 'exct_pol',13)

  endif

  close(13)

  call bcast_scalar_i(nx)
  call bcast_scalar_i(ny)
  call bcast_scalar_i(nz)

  call bcast_scalar_r8(dx)
  call bcast_scalar_r8(dy)
  call bcast_scalar_r8(dz)

  call bcast_scalar_l(periodic)
  call bcast_scalar_i(funct)
  call bcast_scalar_l(flg_bin)

  call bcast_scalar_r8(ekcut)
  !call bcast_scalar_r8(p_fg)
  call bcast_scalar_r8(dt)
  call bcast_scalar_i(nt)
  call bcast_scalar_r8(sm)
  call bcast_scalar_i(ipol)

  nn = nx*ny*nz
  dv = dx*dy*dz
  if(periodic) scale_vh = 1
  if(.not.periodic) scale_vh = 2
  if(funct<0 .or. funct > 1) stop 'funct has to be 0 (LDA) or 1 (PBE)'
end subroutine

subroutine prepare_opaw
    use opaw_mod, only : nx,ny,nz,nhat
    use opaw_mod, only : xmax,ymax,zmax
    use opaw_mod, only : dx,dy,dz
    use mpi_lib_ours, only : rank, sync_mpi
    implicit none

    xmax=dble(nx)*dx/2d0
    ymax=dble(ny)*dy/2d0
    zmax=dble(nz)*dz/2d0

    call calc_nfovnr
    call prep_kpt
    call prep_paw
    call prep_libpaw
    call get_rnel_opaw
    !nhat = compensation charge
    allocate(nhat(nx,ny,nz))
    call vk_prep_opaw
    call ek3d_prep_opaw 
    call get_vloc_ncoret


    if(rank==0)then
      write(*,*) "==============================================="
      write(*,*) "           Finished preparing OPAW             "
      write(*,*) "==============================================="
    endif
    call sync_mpi
end subroutine prepare_opaw

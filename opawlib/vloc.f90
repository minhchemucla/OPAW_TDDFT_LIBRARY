subroutine get_vloc_ncoret
    use opaw_mod
    use libpaw_mod
    use paw_mod
    use atom_mod
    use atom_mod, only : p => pawinfo, at => atominfo
    use m_atm2fft
    use mpi_lib_ours
    implicit none

    real*8, allocatable :: ph1d(:,:),dens_tmp(:)
    real*8, allocatable :: vloc_tot3d(:,:,:), ncoret3d(:,:,:)
    integer :: mgrid

    if(allocated(vloc_tot)) deallocate(vloc_tot)
    allocate(vloc_tot(nn),stat=stat)
    if (stat/=0) stop 'vloc_tot alloc problem'
    allocate(vloc_tot3d(nx,ny,nz),ncoret3d(nx,ny,nz),stat=stat)
    if (stat/=0) stop 'vloc_tot3d alloc problem'
    allocate(ncoret(nn),dens_tmp(nn),stat=stat)
    if (stat/=0) stop 'ncoret, dens_tmp alloc problem'

    if(rank==0) then
        call get_mgrid
        allocate(ph1d(2,3*(2*mgrid+1)*natom),stat=stat)
        if(stat/=0) stop 'ph1d alloc problem'

        ph1d=0d0
        call getph(atindx,natom,nx,ny,nz,ph1d,xred,size(ph1d,1),size(ph1d,2))

        !Get the PS valence density (optn2=2)
        call atm2fft(atindx1,dens_tmp,vloc_tot3d,gmet,gprimd,gsqcut,mgrid,mqgrid_vl,&
            natom,nattyp,nn,ngfft,ntypat,pawtab,ph1d,&
            qgrid_vl,ucvol,vlspl,2)
        !Get the PS core density (optn2=1)
        call atm2fft(atindx1,ncoret3d,vloc_tot3d,gmet,gprimd,gsqcut,mgrid,mqgrid_vl,&
            natom,nattyp,nn,ngfft,ntypat,pawtab,ph1d,&
            qgrid_vl,ucvol,vlspl,1)

        deallocate(ph1d)

        vloc_tot = reshape(vloc_tot3d,(/nn/))
        ncoret   = reshape(ncoret3d,(/nn/))
    endif

    !Minh
    if(.not. periodic) then
        call vloc_tuma_prep_opaw
    endif
    call bcast_r8(ncoret,size(ncoret),0)
    call bcast_r8(vloc_tot,size(vloc_tot),0)
    
    if(rank==0) then
        write(*,*) 'finish preparing vloc,ncoret'
        write(*,*) '========================'
        write(*,*)
    endif
contains
    subroutine get_mgrid
        implicit none

        mgrid=nx
        if (ny>mgrid) then
            mgrid=ny
        endif
        if (nz>mgrid) then
            mgrid=nz
        endif

    end subroutine
end subroutine      

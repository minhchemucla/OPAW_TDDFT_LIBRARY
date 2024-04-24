subroutine prep_kpt
    use opaw_mod
    use mpi_lib_ours
    implicit none
    integer :: nds,nkpt

    nds = max(1,nodes)
    ! uses kpt = (0,0,0) in Brilluoin zone (gamma point)
    nkpt=nds 
    allocate(kpt(3,nkpt))
    kpt=0d0
end subroutine prep_kpt      

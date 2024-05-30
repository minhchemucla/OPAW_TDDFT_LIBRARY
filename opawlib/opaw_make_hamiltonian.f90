subroutine opaw_make_ham_c16(nn,nocc,nstates,wfs,ham)
    use opaw_ham_mod
    use mpi_lib_ours
    implicit none
    integer, intent(in) :: nn, nocc,nstates
    integer :: st, i, ns
    complex*16, intent(in) :: wfs(nn,nstates)!, dens(nn), vks(nn), nhat(nn)
    complex*16, allocatable :: tmp(:,:), sp(:,:), wfs_paw(:,:)!, wf_c(:,:)
    type(opaw_ham_obj) :: ham

    call init_ham(nn,ham) !will only alloc arrays/initilize if ham is fresh

    allocate(wfs_paw(nn,nstates),stat=st);if(st/=0) stop 'wfs_paw'
    call calc_paw_wf(nn,nstates,wfs,wfs_paw)
    call update_dens_paw(nn,nocc,wfs_paw(:,1:nocc),ham)
    call get_pot_opaw(ham)
    call get_dij(ham)
    deallocate(wfs_paw)
end subroutine opaw_make_ham_c16

subroutine opaw_make_ham_r8(nn,nocc,nstates,wfs,ham)
    use opaw_ham_mod
    use mpi_lib_ours
    implicit none
    integer, intent(in) :: nn, nocc,nstates
    integer :: st, i, ns
    real*8, intent(in) :: wfs(nn,nstates)!, dens(nn), vks(nn), nhat(nn)
    complex*16, allocatable :: wfs_c(:,:), tmp(:,:), sp(:,:), wfs_paw(:,:)!, wf_c(:,:)
    type(opaw_ham_obj) :: ham

    call init_ham(nn,ham) !will only alloc arrays/initilize if ham is fresh
    allocate(wfs_paw(nn,nstates),stat=st);if(st/=0) stop 'wfs_paw'
    allocate(wfs_c(nn,nstates),stat=st);if(st/=0) stop 'wfs_paw'
    wfs_c = wfs
    call calc_paw_wf(nn,nstates,wfs_c,wfs_paw)
    call update_dens_paw(nn,nocc,wfs_paw(:,1:nocc),ham)
    call get_pot_opaw(ham)
    call get_dij(ham)
    deallocate(wfs_paw,wfs_c)
end subroutine opaw_make_ham_r8

subroutine opaw_make_ham_tdh_c16(nn,nocc,nstates,wfs,ham)
    !Only call this after vxc is set
    use opaw_ham_mod
    use mpi_lib_ours
    implicit none
    integer, intent(in) :: nn, nocc,nstates
    integer :: st, i, ns
    complex*16, intent(in) :: wfs(nn,nstates)!, dens(nn), vks(nn), nhat(nn)
    complex*16, allocatable :: tmp(:,:), sp(:,:), wfs_paw(:,:)!, wf_c(:,:)
    type(opaw_ham_obj) :: ham

    if(.not.ham%alloc_flg) stop "need to have set vxc before hand for TDH c16"
    allocate(wfs_paw(nn,nstates),stat=st);if(st/=0) stop 'wfs_paw'
    call calc_paw_wf(nn,nstates,wfs,wfs_paw)
    call update_dens_paw(nn,nocc,wfs_paw(:,1:nocc),ham)
    call get_pot_static_xc_opaw(ham)
    call get_dij(ham)
    deallocate(wfs_paw)
end subroutine opaw_make_ham_tdh_c16

subroutine opaw_make_ham_tdh_r8(nn,nocc,nstates,wfs,ham)
    use opaw_ham_mod
    use mpi_lib_ours
    implicit none
    integer, intent(in) :: nn, nocc,nstates
    integer :: st, i, ns
    real*8, intent(in) :: wfs(nn,nstates)!, dens(nn), vks(nn), nhat(nn)
    complex*16, allocatable :: wfs_c(:,:), tmp(:,:), sp(:,:), wfs_paw(:,:)!, wf_c(:,:)
    type(opaw_ham_obj) :: ham

    if(.not.ham%alloc_flg) stop "need to have set vxc before hand for TDH r8"
    allocate(wfs_paw(nn,nstates),stat=st);if(st/=0) stop 'wfs_paw'
    allocate(wfs_c(nn,nstates),stat=st);if(st/=0) stop 'wfs_paw'
    wfs_c = wfs
    call calc_paw_wf(nn,nstates,wfs_c,wfs_paw)
    call update_dens_paw(nn,nocc,wfs_paw(:,1:nocc),ham)
    call get_pot_static_xc_opaw(ham)
    call get_dij(ham)
    deallocate(wfs_paw,wfs_c)
end subroutine opaw_make_ham_tdh_r8

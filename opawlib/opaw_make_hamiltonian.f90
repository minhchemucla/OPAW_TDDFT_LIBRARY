subroutine opaw_make_hamiltonian(nn,nocc,nstates,wfs,ham)
    use opaw_ham_mod
    use mpi_lib_ours
    implicit none
    integer, intent(in) :: nn, nocc,nstates
    integer :: st, i, ns, is, ie
    complex*16, intent(in) :: wfs(nn,nstates)!, dens(nn), vks(nn), nhat(nn)
    complex*16, allocatable :: tmp(:,:), sp(:,:), wfs_paw(:,:)!, wf_c(:,:)
    type(opaw_ham_obj) :: ham

    call calc_is_ie(is,ie,nstates)

    allocate(wfs_paw(nn,nstates),stat=st);if(st/=0) stop 'wfs_paw'
    call calc_soft_ortho_wf(nn,nstates,0,wfs,wfs_paw)
    call update_dens(nn,nocc,wfs_paw(:,1:nocc),ham)
    call get_pot_opaw(ham)
    call get_dij(ham)
    deallocate(wfs_paw)
end subroutine opaw_make_hamiltonian

subroutine opaw_make_hamiltonian_r8(nn,nocc,nstates,wfs,ham)
    use opaw_ham_mod
    use mpi_lib_ours
    implicit none
    integer, intent(in) :: nn, nocc,nstates
    integer :: st, i, ns, is, ie
    real*8, intent(in) :: wfs(nn,nstates)!, dens(nn), vks(nn), nhat(nn)
    complex*16, allocatable :: wfs_c(:,:), tmp(:,:), sp(:,:), wfs_paw(:,:)!, wf_c(:,:)
    type(opaw_ham_obj) :: ham

    call calc_is_ie(is,ie,nstates)

    allocate(wfs_paw(nn,nstates),stat=st);if(st/=0) stop 'wfs_paw'
    allocate(wfs_c(nn,nstates),stat=st);if(st/=0) stop 'wfs_paw'
    wfs_c = wfs
    call calc_soft_ortho_wf(nn,nstates,0,wfs_c,wfs_paw)
    call update_dens(nn,nocc,wfs_paw(:,1:nocc),ham)
    call get_pot_opaw(ham)
    call get_dij(ham)
    deallocate(wfs_paw,wfs_c)
end subroutine opaw_make_hamiltonian_r8

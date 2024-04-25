subroutine opaw_make_hamiltonian(nn,nocc,nstates,wfs,ham)
    !use opaw_mod, only : dens,nx, ny, nz, nn
    !use opaw_mod, only : nk_loc, h_type, nb, phit_tot, phit
    !use opaw_mod, only : is_opaw_nocc, ie_opaw_nocc, h_type
    use opaw_ham_mod
    use mpi_lib_ours
    implicit none
    integer, intent(in) :: nn, nocc,nstates
    integer :: st, i, ns, is, ie
    complex*16, intent(in) :: wfs(nn,nstates)!, dens(nn), vks(nn), nhat(nn)
    complex*16, allocatable :: tmp(:,:), sp(:,:), wf_soft(:,:)!, wf_c(:,:)
    type(opaw_ham_obj) :: ham

    !nn = ham%n
    !nocc = ham%nocc
    !do i=1,nocc
    !  if(rank==0) write(*,*) 'rank, i, phi_bar_pert', rank, i, sum(dble(wfs(:,i)))
    !enddo
    call calc_is_ie(is,ie,nstates)

    allocate(wf_soft(nn,nstates),stat=st);if(st/=0) stop 'wf_soft'
    !allocate(wf_c(nn,nocc),stat=st);if(st/=0) stop 'wf_c'
    !write(*,*) 'pre scatter', sum(wf(:,:))
    if(ham%h_type .eq. 1) then
      call calc_soft_ortho_wf(nn,nstates,0,wfs,wf_soft)
    else
      wf_soft = wfs
    endif
    !do i=1,nocc
    !  if(rank==0) write(*,*) 'rank, i, sn phi_bar_pert', rank, i, sum(dble(wf_soft(:,i)))
    !enddo
    !write(*,*) 'post scatter', sum(wf_soft(:,:))

    call update_dens(nn,nocc,wf_soft(:,1:nocc),ham)
    call get_pot_opaw(ham)
    call get_dij(ham)
    deallocate(wf_soft)
end subroutine opaw_make_hamiltonian

subroutine opaw_make_hamiltonian_r8(nn,nocc,nstates,wfs,ham)
    use opaw_mod, only : dv
    !use opaw_mod, only : nk_loc, h_type, nb, phit_tot, phit
    !use opaw_mod, only : is_opaw_nocc, ie_opaw_nocc, h_type
    use opaw_ham_mod
    use mpi_lib_ours
    implicit none
    integer, intent(in) :: nn, nocc,nstates
    integer :: st, i, ns, is, ie
    real*8, intent(in) :: wfs(nn,nstates)!, dens(nn), vks(nn), nhat(nn)
    complex*16, allocatable :: wfs_c(:,:), tmp(:,:), sp(:,:), wf_soft(:,:)!, wf_c(:,:)
    type(opaw_ham_obj) :: ham

    !nn = ham%n
    !nocc = ham%nocc
    !do i=1,nocc
    !  if(rank==0) write(*,*) 'rank, i, phi_bar_pert', rank, i, sum(dble(wfs(:,i)))
    !enddo
    call calc_is_ie(is,ie,nstates)

    allocate(wf_soft(nn,nstates),stat=st);if(st/=0) stop 'wf_soft'
    allocate(wfs_c(nn,nstates),stat=st);if(st/=0) stop 'wf_soft'
    !allocate(wf_c(nn,nocc),stat=st);if(st/=0) stop 'wf_c'
    !write(*,*) 'pre scatter', sum(wf(:,:))
    wfs_c = wfs
    if(ham%h_type .eq. 1) then
      call calc_soft_ortho_wf(nn,nstates,0,wfs_c,wf_soft)
    else
      wf_soft = wfs
    endif
    !do i=1,nocc
    !  if(rank==0) write(*,*) 'rank, i, sn phi_bar_pert', rank, i, sum(dble(wf_soft(:,i)))
    !enddo
    !write(*,*) 'post scatter', sum(wf_soft(:,:))

    call update_dens(nn,nocc,wf_soft(:,1:nocc),ham)
    call get_pot_opaw(ham)
    call get_dij(ham)
    deallocate(wf_soft,wfs_c)
end subroutine opaw_make_hamiltonian_r8

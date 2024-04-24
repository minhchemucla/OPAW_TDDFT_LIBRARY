program main
  use main_mod
  use opaw_mod, only : nhat
  use ham_obj
  use mpi_lib_ours
  implicit none
  integer :: st
  real*8  :: exx

  call prepare_mpi_lib
  call read_input
  call alloc_dens_pots
  call read_wfs

  !OPAW STUFF
  call opaw_libpaw_prepare
  call get_rnel_opaw
  call init_ham(nn,h_type,ham)
  call opaw_make_hamiltonian(nn,nocc,nstates,wfs,ham)

  call calc_eig(ham)
  call test_exx

  call tddft

  call finalize_mpi_lib
  contains
    subroutine alloc_dens_pots
      implicit none
      
      allocate(dens(nn),vxc(nn),vks(nn),stat=st);if(st/=0) stop 'alloc dens'
    end subroutine alloc_dens_pots

    subroutine test_exx
      implicit none
      complex*16 :: wf_o(nn,nstates)
      if(h_type.eq.0) then
        call calc_soft_ortho_wf(nn,nstates,1,wfs,wf_o)
      else if (h_type .eq. 1) then
        wf_o= wfs
      endif
      call exx_expect_opaw(nn,nocc,nstates,wf_o(:,nocc),wf_o(:,nocc),wf_o,exx)
      call sync_mpi
    end subroutine test_exx
end program main

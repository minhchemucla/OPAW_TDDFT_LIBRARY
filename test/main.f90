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
  call prepare_opaw
  call opaw_make_ham(nn,nocc,nstates,wfs,ham)

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
      call exx_expect_opaw(nn,nocc,nstates,wfs(:,nocc),wfs(:,nocc),wfs,exx)
      call sync_mpi
    end subroutine test_exx
end program main

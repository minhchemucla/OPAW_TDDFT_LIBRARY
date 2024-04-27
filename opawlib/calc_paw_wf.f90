subroutine calc_paw_wf(nn,nstates,opaw_wf,paw_wf)
  use mpi_lib_ours
  implicit none
  integer :: nn,nstates
  integer :: is,ie,i,st,ns
  real*8  :: n
  complex*16 :: opaw_wf(nn,nstates), paw_wf(nn,nstates)
  complex*16, allocatable :: tmp(:,:), sp(:,:)
  
  n=-0.5d0 
  call calc_is_ie(is,ie,nstates)
  !for hbar=shs, need to form phi_tilde=S^-1/2 phi_bar
  ns = ie-is+1
  !write(*,*) 'is,ie,ns', is_opaw_nocc, ie_opaw_nocc, ns
  allocate(sp(nn,ns), tmp(nn,ns), stat=st)
  if(st/=0) stop 'allocate sp, tmp make_ham'

  call scatterv_c16(opaw_wf,tmp,size(opaw_wf),size(tmp),0)
  !sp = tmp
  do i=1,ns
     call sn_phi(tmp(:,i),sp(:,i),n)
  enddo
  call gatherv_c16(sp,paw_wf,size(sp),size(paw_wf),0)

  deallocate(tmp,sp)
end subroutine calc_paw_wf

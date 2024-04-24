subroutine calc_soft_ortho_wf(nn,nstates,flag,ortho_wf,soft_wf)
  use mpi_lib_ours
  implicit none
  integer :: nn,nstates
  integer :: is,ie,i,st,ns,flag
  real*8  :: n
  complex*16 :: ortho_wf(nn,nstates), soft_wf(nn,nstates)
  complex*16, allocatable :: tmp(:,:), sp(:,:)
  
  if(flag.eq.0) then
    n=-0.5d0 
  else if(flag.eq.1) then
    n= 0.5d0 
  else
    stop 'flag in calc_soft ortho should be 0 (o->s) or 1 (s->o)'
  endif
  call calc_is_ie(is,ie,nstates)
  !for hbar=shs, need to form phi_tilde=S^-1/2 phi_bar
  ns = ie-is+1
  !write(*,*) 'is,ie,ns', is_opaw_nocc, ie_opaw_nocc, ns
  allocate(sp(nn,ns), tmp(nn,ns), stat=st)
  if(st/=0) stop 'allocate sp, tmp make_ham'

  call scatterv_c16(ortho_wf,tmp,size(ortho_wf),size(tmp),0)
  !sp = tmp
  do i=1,ns
     call sn_phi(tmp(:,i),sp(:,i),n)
  enddo
  call gatherv_c16(sp,soft_wf,size(sp),size(soft_wf),0)

  deallocate(tmp,sp)
end subroutine calc_soft_ortho_wf

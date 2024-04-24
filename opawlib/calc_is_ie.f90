subroutine calc_is_ie(is,ie,nstates)
!calculates index of slices of nstates over the nodes
!for example for 5 states using 2 nodes 
!   rank=0 is=1,ie=3
!   rank=1 is=2,ie=5
  use mpi_lib_ours, only : rank, nodes, sync_mpi
  implicit none
  integer, intent(out) :: is, ie
  integer, intent(in) :: nstates
  integer :: nds, ns, nextra


  nds = max(1,nodes)
  ns = floor(dble(nstates)/dble(nds)) 
  nextra = mod(nstates,nds)

  is = rank*ns + 1
  ie = is + ns - 1
  if(nextra > 0) then
    if(rank .eq. 0) then
      ie = ie + 1
    else if(rank < nextra .and. rank > 0) then
      is = is + rank
      ie = ie + rank + 1
    else
      is = is + nextra
      ie = ie + nextra
    endif
  endif

  !write(*,*) 'rank, is, ie', rank, is, ie
  !call flush(6)
  call sync_mpi
end subroutine

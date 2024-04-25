!This module contains a hamiltonian object which holds the minimum information
!necessary to describe a OPAW Hamiltonian for the system of interest
!Some stuff is assumed to be common through all Hamiltonians in this program namely:
!  ekcut, vloc_tot, the information in atom_mod besides rhoij and dij, the info in paw_mod
!wfs    = (o)paw wavefunctions
!rhoij  = sum_k<p_i|psi_k><psi_k|p_j> where psi_k is the soft paw wavefunction and
!              p_i are the untransformed paw projector functions
!dij    = D_ij non-local terms
!vks    = Kohn Sham potential = V_H+V_XC+V_loc
!vxc    = exchange-correlation potential
!dens   = density
!nhat   = charge compensation density
!h_type = h_type
module opaw_ham_mod
  use atom_mod
  implicit none
  type opaw_ham_obj
    integer :: h_type !ca*ca
    integer :: n, ns!, nocc
    real*8,  allocatable :: nhat(:)
    real*8,  allocatable :: dens(:) !density made from paw wfs
    real*8,  allocatable :: dens_o(:) !density made from opaw wfs
    real*8,  allocatable :: vks(:)
    real*8,  allocatable :: vxc(:)
    real*8,  allocatable :: vh(:)
    type(atom), allocatable :: at(:)  !just need for rhoij and dij
  end type
  contains
    subroutine init_ham(n,h_type,ham)
      !Only ran after prepare_paw
      implicit none
      integer :: st,it,ia,ms
      integer :: n,ns,h_type
      type(opaw_ham_obj) :: ham

      !ham%n = n
      !ham%ns = ns
      !ham%nocc = nocc
      ham%h_type = h_type
      if(h_type < 0 .or. h_type > 1) stop 'h_type should be 0 or 1'
      allocate(ham%nhat(n),stat=st); if(st/=0)stop 'ham%nhat alloc'
      allocate(ham%vks(n),stat=st); if(st/=0)stop 'ham%vks alloc'
      allocate(ham%vxc(n),stat=st); if(st/=0)stop 'ham%vxc alloc'
      allocate(ham%vh(n),stat=st); if(st/=0)stop 'ham%vh alloc'
      allocate(ham%dens(n),stat=st); if(st/=0)stop 'ham%dens alloc'
      allocate(ham%dens_o(n),stat=st); if(st/=0)stop 'ham%dens alloc'
      !allocate(ham%wfs(n,ns),stat=st); if(st/=0)stop 'ham%wfs alloc'
      allocate(ham%at(natom),stat=st); if(st/=0) stop 'ham%at alloc'

      do ia=1,natom
         it=atom_map(ia)
         ms=pawinfo(it)%mstates
         allocate(ham%at(ia)%rhoij(ms,ms),stat=stat)
         if(stat/=0) stop 'problem allocated rhoij, rank/=0' 
         allocate(ham%at(ia)%dij(ms,ms),stat=stat)
         if(stat/=0) stop 'problem allocated dij, rank/=0' 
         ham%at(ia)%rhoij=0d0
         ham%at(ia)%dij=0d0
      enddo

    end subroutine init_ham

    !subroutine set_ham_wfs(n,ns,ham,wfs)
    !  implicit none
    !  integer :: st,it,ia,ms
    !  integer :: n,ns
    !  real*8  :: wfs(n,ns)
    !  type(opaw_ham_obj) :: ham

    !  if(n.ne.ham%n) then
    !    write(*,*) 'n, ham%n wfs', n, ham%n
    !    stop
    !  endif
    !  if(ns.ne.ham%ns) then
    !    write(*,*) 'ns, ham%ns wfs', n, ham%ns
    !    stop
    !  endif
    !  ham%wfs = wfs
    !end subroutine set_ham_wfs

    !subroutine set_ham_var(n,ham,pot,flag)
    !  implicit none
    !  integer :: st,it,ia,ms
    !  integer :: n, flag
    !  real*8  :: pot(n)
    !  type(opaw_ham_obj) :: ham

    !  if(n.ne.ham%n) then
    !    write(*,*) 'n, ham%n ', n, ham%n
    !    stop
    !  endif
    !  select case(flag)
    !    case(1)
    !      ham%dens = pot
    !    case(2)
    !      ham%vks  = pot
    !    case(3)
    !      ham%vxc  = pot
    !    case(4)
    !      ham%vh  = pot
    !    case(5)
    !      ham%nhat  = pot
    !    case default
    !      write(*,*) 'invalid set_ham_var flag'
    !      stop
    !  end select
    !      ham%dens = pot
    !end subroutine set_ham_var

    !subroutine set_ham(n,ns,ham,wfs,dens,vks,vxc,vh,nhat)
    !  implicit none
    !  integer :: st,it,ia,ms
    !  integer :: n, ns
    !  real*8  :: wfs(n,ns),dens(n),vks(n),vxc(n),vh(n),nhat(n)
    !  type(opaw_ham_obj) :: ham

    !  if(n.ne.ham%n) then
    !    write(*,*) 'n, ham%n ham', n, ham%n
    !    stop
    !  endif
    !  if(ns.ne.ham%ns) then
    !    write(*,*) 'n, ham%ns ham', n, ham%ns
    !    stop
    !  endif
    !  call set_ham_wfs(n,ns,ham,wfs)
    !  call set_ham_var(n,ham,dens,1)
    !  call set_ham_var(n,ham,vks,2)
    !  call set_ham_var(n,ham,vxc,3)
    !  call set_ham_var(n,ham,vh,4)
    !  call set_ham_var(n,ham,nhat,5)
    !end subroutine set_ham
end module opaw_ham_mod

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
module opaw_ham_mod
  use atom_mod
  implicit none
  type opaw_ham_obj
    integer :: n, ns!, nocc
    real*8,  allocatable :: nhat(:)
    real*8,  allocatable :: dens(:) !density made from paw wfs
    real*8,  allocatable :: dens_o(:) !density made from opaw wfs
    real*8,  allocatable :: vks(:)
    real*8,  allocatable :: vxc(:)
    real*8,  allocatable :: vh(:)
    type(atom), allocatable :: at(:)  !just need for rhoij and dij
    logical              :: alloc_flg=.false.
  end type
  contains
    subroutine init_ham(n,ham)
      !Only ran after prepare_paw
      implicit none
      integer :: st,it,ia,ms
      integer :: n,ns
      type(opaw_ham_obj) :: ham

      if(ham%alloc_flg) return

      !ham%n = n
      !ham%ns = ns
      !ham%nocc = nocc
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

      ham%alloc_flg = .true.

    end subroutine init_ham
end module opaw_ham_mod

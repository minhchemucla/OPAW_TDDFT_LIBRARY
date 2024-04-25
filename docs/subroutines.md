 # OPAW Subroutines


The following are the main subroutines that you'll need to set up and use the OPAW library. There are a few subroutines that I expect to be redundant if you are trying to implement this library with any of the Neuhauser group's code so I renamed some subroutines such as vk_prep to vk_prep_opaw to avoid conflict. Feel free to change if using a specialized vk_prep routine. 

The  [opaw_libpaw_prepare](#opaw-libpaw-prepare) routine has to be called after the input file is read and the system variables (`nx,ny,nz,dx,...`) are set, before any of the following subroutines can be used.

## Subroutines
 1. [opaw_libpaw_prepare](#opaw-libpaw-prepare)
 3. [init_ham](#init_ham)
 4. [opaw_make_hamiltonian](#opaw_make_ham)
 5. [rk4_prop_opaw](#rk4_prop_opaw)
 6. [sn_phi](#sn_phi)
 7. [exx_expect_opaw](#exx_expect_opaw)


##  <a id="opaw-libpaw-prepare"></a> 1. opaw_libpaw_prepare
### input
none 

### output
none 

### description
Reads in the pawfile files, prepares the orthogonal projectors and matrix elements.
  It also prepares vk for applying the Coulomb potential and ek for applying the
  kinetic energy operator.


## <a id="init_ham"></a> 3. init_ham(nn,h_type,ham)

###   input 
    integer :: nn          ! number of grid points nn=nx*ny*nz
    integer :: h_type      ! hamiltonian type =0 (S^-1H) =1 (S^-1/2HS^-1/2)

### output
    hamiltonian_obj :: ham 
    
###  description
   This routine allocates the necessary information that would define a hamiltonian for a system. See the README.md for more information on the structure of the Hamiltonian.

## <a id="opaw_make_ham"></a> 4. opaw_make_hamiltonian(nn,nocc,nstates,wfs,ham)


###  input
    integer     :: nn           ! number of grid points nn=nx*ny*nz
    integer     :: nocc         ! number of occupied states
    integer     :: nstates      ! number of states
    complex*16  :: wfs(nn,nocc) ! non-orthogonal (h_type 0) or 
                                ! orthogonal (h_type 1) wavefunctions
                              
###  output
    hamiltonian_obj :: ham 

###  description
   This subroutine calculates the PAW density matrices and potentials.

## <a id="rk4_prop_opaw"></a> 5. rk4_prop_opaw(nn,nocc,nstates,dt,p,ham)


###  input
    integer     :: nn           ! number of grid points nn=nx*ny*nz
    integer     :: nocc         ! number of occupied states
    integer     :: nstates      ! number of states
    real*8      :: dt           ! time-step
    hamiltonian_obj :: ham      ! input hamiltonian

###  output
    complex*16  :: p(nn,nstates)   ! Propagated wavefunctions

###  description
   Given a hamiltonian and time step, use 4th order Runge-Kutta to 
    propagate all the wavefunctions a single time step. If 
    ham%h_type=0 the wavefunctions are assumed to be the
    non-orthogonal and propagated with the S^-1 H hamiltonian
    and for ham%h_type=1 orthogonal wavefunctions with S^-1/2HS^-1/2.

## <a id="sn_phi"></a> 6.  sn_phi(pin,sp,n)



### input
    integer     :: n               ! power in S^n
    complex*16  :: pin(nx,ny,nz)   ! input wavefunction

### output
    complex*16  :: sp(nx,ny,nz)    ! sp = S^n pin

### description
   Applies S^n to `pin` and returns the result as `sp`.

## <a id="exx_expect_opaw"></a> 7. exx_expect_opaw(nn,nocc,nstates,psi_i,psi_j,psi_n,exx)
 
 
### input
    integer     :: nn                ! number of grid points nn=nx*ny*nz
    integer     :: nocc              ! number of occupied states
    integer     :: nstates           ! number of states
    complex*16  :: psi_i(nn)         ! bra
    complex*16  :: psi_j(nn)         ! ket
    complex*16  :: psi_n(nn,nstates) ! orthogonal wavefunctions

### output
    real*8      :: exx               ! expectation value

### description
   Calculates the expectation value,  `exx`, of the exchange fock exchange operator
    using the orthogonal pseudowavefunctions=<psi_i|V^hat_x|psi_j>.

-----------------------------------------------------------------
              4.2 Overlapping subroutines
-----------------------------------------------------------------
I expect some subroutines to be redundant with other Neuhauser group
code. In particular:


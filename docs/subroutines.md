 # OPAW Subroutines
This document is split into two sections. There is a [main subroutines](#main) section that details the subroutines related to preparing, applying operators, and time propagation. Section [OPAW Hamiltonian](#ham) deals with the subroutines that go into calculating the Hamiltonian. All subroutines listed here are generally stored in a file with the subroutine name with the `.f90` extension.

## *Notes for Neuhauser Group Members Specifically*
There are a few subroutines related to calculating $V_{H}$ that I expect to be redundant with code from the Neuhauser group so I renamed some subroutines to avoid conflict. The following is a list of subroutines that were created and can be changed to the equivalent version in the main code in the library.

 1. `vh_sub_opaw`
 2. `vk_prep_opaw`
 3. `prep_vk_opaw`

The subroutine `vxc_libxc` also probably already exists somewhere in the main code, so remove it if needed.

# <a id="main"></a> $\color{Red}\rm{Main \ Subroutines}$

 1. [prepare_opaw](#opaw-libpaw-prepare): Reads OPAW files, processes projectors, and prepares for other OPAW routines.
 2. [opaw_make_hamiltonian](#opaw_make_ham) : calculates potentials and terms
 4. [rk4_prop_opaw](#rk4_prop_opaw) : 4th-order Runge-Kutta Time Propagation Step
 5. [sn_phi](#sn_phi) : Applies $S^n$
 6. [opaw_ham](#opaw_ham) : Applies $S^{-1/2}HS^{-1/2}$ 
 7.  [paw_ham](#paw_ham) : Applies $H$
 8. [exx_expect_opaw](#exx_expect_opaw): Calculates expectation value of exact Fock exchange operator
 9.  [proj_paw](#proj_paw) : calculates overlap of PAW projectors with wavefunction
 10. [proj_opaw](#proj_opaw): calculates overlap of OPAW projectors with wavefunction

##  <a id="opaw-libpaw-prepare"></a> $\color{blue}\rm{1.\  prepare\_-opaw}$

### usage
	call prepare_opaw
### input
none 

### output
none 

### description
Reads in the pawfile files and prepares the orthogonal projectors and matrix elements.
It also prepares vk for applying the Coulomb potential and ek for applying the
kinetic energy operator.

## <a id="opaw_make_ham"></a> $\color{blue}\rm{2.\ opaw\_-make\_- hamiltonian}$

### usage
	call opaw_make_hamiltonian(nn,nocc,nstates,wfs,ham)

###  input
    integer     :: nn           ! number of grid points nn=nx*ny*nz
    integer     :: nocc         ! number of occupied states
    integer     :: nstates      ! number of states
    complex*16  :: wfs(nn,nocc) ! OPAW Wavefunctions
                              
###  output
    opaw_ham_obj :: ham 

###  description
Takes in the OPAW wfs and from them makes the PAW functions by applying $S^{-1/2}$ then calculates the PAW density, density matrix, compensation charges, and potentials and stores this information in `ham`.


## <a id="rk4_prop_opaw"></a> $\color{blue}\rm{3.\ rk4\_-prop\_-opaw}$
### usage 
   	call rk4_prop_opaw(nn,nocc,nstates,dt,p,ham)

###  input
    integer      :: nn           ! number of grid points nn=nx*ny*nz
    integer      :: nocc         ! number of occupied states
    integer      :: nstates      ! number of states
    real*8       :: dt           ! time-step
    opaw_ham_obj :: ham         ! input hamiltonian

###  output
    complex*16   :: p(nn,nstates)   ! Propagated wavefunctions

###  description
   Given a Hamiltonian and time step, use 4th-order Runge-Kutta to propagate all the wavefunctions a single time step. 

## <a id="sn_phi"></a> $\color{blue}\rm{4.\ sn\_-phi}$

### usage 
   	call sn_phi(pin,sp,n)
   	
### input
    integer     :: n               ! power in S^n
    complex*16  :: pin(nx,ny,nz)   ! input wavefunction

### output
    complex*16  :: sp(nx,ny,nz)    ! sp = S^n pin

### description
   Applies $S^n=\sum_i\ket{\eta^a_i}\Big((1-o_i^a)^{n}-1\Big)\bra{\eta^a_i}$ to `pin` and returns the result as `sp`. The routine assumes the shape of `pin` and `sp` to be `(nx,ny,nz)`.

## <a id="opaw_ham"></a> $\color{blue}\rm{5.\ opaw\_-ham}$ 

### usage 
   	call opaw_ham(ham,pin,pout)
   	
### input
    opaw_ham_obj :: ham              ! PAW hamiltonian
    complex*16   :: pin(nx,ny,nz)    ! orbital for OPAW Ham to be applied onto

### output
    complex*16   :: pout(nx,ny,nz)    ! pout=S^-1/2HS^-1/2pin

### description
   Applies $S^{-1/2}HS^{-1/2}$ to `pin` and returns the result as `pout` under Hamiltonian `ham`. The routine assumes the shape of `pin` and `pout` to be `(nx,ny,nz)`.

## <a id="paw_ham"></a> $\color{blue}\rm{6.\ paw\_-ham}$ 

### usage 
   	 paw_ham(ham,pin,hp)
   	
### input
    opaw_ham_obj :: ham              ! PAW hamiltonian
    complex*16   :: pin(nx,ny,nz)    ! orbital for OPAW Ham to be applied onto

### output
    complex*16  :: hp(nx,ny,nz)      ! hp = H pin

### description
   Applies $H$ to `pin` and returns the result as `pout` under Hamiltonian `ham`. Starts with the kinetic energy, then the local terms, and finally the nonlocal $D^a_{ij}$ terms.

## <a id="exx_expect_opaw"></a>$\color{blue}\rm{7.\ exx\_-expect\_-opaw}$  
 
### usage 
	call exx_expect(nn,nocc,nstates,psi_i,psi_j,psi_n,exx)
	
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
   Calculates the expectation value of the exchange fock exchange operator with the orthogonal pseudowavefunctions, $\ket{\psi_i}$:  
   $$\rm{exx}=\braket{\psi_i|V_x|\psi_j}=\sum_n^{Nocc}\int drdr'\frac{\psi_i(r)\psi_n(r)\psi_j(r')\psi_n(r')}{|r-r'|}.$$
   *Note:* In the above equation the $\psi_n(r)$ and $\psi_n(r')$ should be complex conjugate but Github flavored Markdown LaTeX bugs out when I try.



## <a id="proj_paw"></a> $\color{blue}\rm{8.\ proj\_-paw}$  

### description
  Caculates the overlap of the original PAW projector functor onto `pin` and returns the results in `ca`, ca(i) = $\braket{p^a_i|p_{in}}$. For an example, refer to the [`paw_ham`](#paw_ham) subroutine and the subroutine `addvnl` where $V_{Nl}=\sum_{ij}\ket{p^a_i}D^a_{ij}\bra{p^a_j}$ is applied.
### usage 
   	call proj_paw(ia,pin,ca,ms,ik)
   	
### input
	integer     :: ia  !index of atom in system
	complex*16  :: pin(nx,ny,nz) !input wavefunction
	integer     :: ms  !number of projector states
	integer     :: ik  !legacy k-point related. It will be always be 1 in the code (gamma)
	
    
### output
	complex*16  :: ca(ms)     !overlap terms with PAW projectors and wavefunctions


## <a id="proj_opaw"></a> $\color{blue}\rm{9.\ proj\_-opaw}$  
### usage 
   	call proj_paw(ia,pin,ca,ms,ik)
   	
### input
	integer     :: ia  !atom index
	complex*16  :: pin(nx,ny,nz)
	integer     :: ms  !number of projector states
	integer     :: ik  !legacy k-point related. It will be always be 1 in the code (gamma)
	
    
### output
	complex*16  :: ca(ms)     !overlap terms with transformed projectors and wavefunctions

### description
  Calculates the overlap of the original PAW projector functor onto `pin` and returns the results in `ca`, ca(i) = $\braket{\eta^a_i|p_{in}}$. For an example, refer to how the projectors are applied in the [`sn_phi`](#sn_phi) subroutine code.

## <a id="ham"></a> $\color{red}\rm{OPAW \ Hamiltonian  \ Related \ Subroutines}$
To manipulate the individual parts of the Hamiltonian such as in Time-dependent Hartree Propagation where only the Hartree potential is updated, the related subroutines that are in [opaw_make_hamiltonian](#opaw_make_ham) are described here. 

1. [calc_paw_opaw_wf](#calc_paw_wf)
2. [update_dens_paw](#update_dens_paw)
3. [get_pot_opaw](#get_pot_opaw)
4. [get_dij](#get_dij)

## <a id="calc_paw_wf"></a> $\color{blue}\rm{1.\ calc\_-paw\_-wf}$  
### usage 
   	call calc_paw_wf(nn,nstates,opaw_wf,paw_wf)
   	
### input
    integer     :: nn                  ! number of grid points nn=nx*ny*nz
    integer     :: nstates             ! number of states
    complex*16  :: opaw_wf(nn,nstates) ! orthogonal PAW pseudowavefunctions
    

### output
    complex*16  :: paw_wf(nn,nstates)  ! PAW pseudowavefunctions

### description
   Applies $S^{-1/2}$ to `opaw_wf` to get `paw_wf` parallelizing over the number of nstates. 


## <a id="update_dens_paw"></a> $\color{blue}\rm{2.\ update\_-dens\_-paw}$ 
### usage 
   	call update_dens_paw(nn,nocc,wf,ham)
   	
### input
    integer     :: nn               ! number of grid points nn=nx*ny*nz
    integer     :: nocc             ! number of occupied states
    complex*16  :: wf(nn,nstates)   ! PAW pseudowavefunctions
    

### output
    opaw_ham_obj :: ham%dens    !PAW Density
    opaw_ham_obj :: ham%nhat    !PAW Compensation charge
    opaw_ham_obj :: ham%rhoij   !PAW Density matrix

### description
  Caculates the $n(r)$, $\rho^a_{ij}$, and $\hat{n}(r)$ from `wf` and stores the results in`ham%rhoij`, `ham%dens` and `ham%nhat`. 
 
 
## <a id="get_pot_opaw"></a>  $\color{blue}\rm{3.\ get \_- pot\_- opaw}$ 
### usage 
   	call get_pot_opaw(ham)
   	
### input
	opaw_ham_obj :: ham
    
### output
	opaw_ham_obj :: ham%vxc
	opaw_ham_obj :: ham%vks
	opaw_ham_obj :: ham%vh

### description
  Caculates the $V_{XC}$, $V_H$, and $V_{KS}$ from `ham%dens` and `ham%nhat` and stores the results in `ham%vxc`, `ham%vh`, `ham%vks`. The subroutine calls on `vxc_libxc` and `vh_sub_opaw` for the first two.
  
## <a id="get_dij"></a> $\color{blue}\rm{4.\ get\_-dij}$ 
### usage 
   	call get_dij(ham)
   	
### input
	opaw_ham_obj :: ham
    
### output
	opaw_ham_obj :: ham%dij     !PAW Dij nonlocal terms

### description
  Caculates  $D^a_{ij}$ from $\rho^a_{ij}$, $V_{KS}$, and $V_{XC}$ using ABINIT subroutines.


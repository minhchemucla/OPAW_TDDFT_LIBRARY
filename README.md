OPAW LIBRARY CODE DOCUMENTATION

This library code provides the necessary ingredients to
  1. prepare and apply the OPAW Hamiltonian,
  2. perform RK4 time propagation,
  3. calculculate the expectation value of the exact exchange operator.
The code only has gamma-point nspin=1 capabilities currently.

This library comes with a test program for you to look at to see
examples of how to use the library.


Chapters
  0. Prerequisites
  1. OPAW Files
  2. OPAW Interface
    2.1 Variables to import
    2.2 Internal OPAW Variables
  3. Structure of the code
  4. OPAW Subroutines
  5. Compiling

=================================================================
              0. Prerequisites
=================================================================
This library was mainly designed for members of the Neuhauser group
so there some particulars to sucessfully implementing this library.

  1. We have our own separate library that has various codes 
  to perform numerical functions and whatnot, but the most important
  is that we have a 0_library_mpi_module.f90 file that contains a
  module that simplifies the set up, execution, and deallocation of 
  relevant MPI stuff. You can find the 0_library_mpi_module.f90 file
  in the lib directory inside the example directory. The module 
  "mpi_libs_ours" is used extensively in the library for parallelization
  purposes. If you are not using the lib in the example directory,
  place the 0_library_mpi_module.f90 somewhere in your code and
  make sure to compile it before compiling the OPAW library.

  2. The library uses the libxc library to calculate the 
  exchange correlation potential. We have our own interface which is
  contained in the XCI directory in the example directory.

The below are numerical packages needed to make the library work.
You can find them on most suipercomputers.

  3. FFTW3

  4. BLAS

  5. LAPACK

=================================================================
              1. OPAW Files
=================================================================
The following are files and needed to make the library to work. 
  pawfiles:
    This file should contain a list of elements and the 
    corresponding PAW potential. For example:

    H     H.LDA_PW-JTH.xml
    C     C.LDA_PW-JTH.xml

  PAW potentials:
    The PAW potentials in the pawfiles should be in the same directory as the pawfiles itself.
  cnt.ini:
    Coordinate file of system in Bohr.

=================================================================
              2. OPAW Interface
=================================================================
-----------------------------------------------------------------
              2.1 Variables to import
-----------------------------------------------------------------
The 0_parameters.f90 file serves to store variables that are important in the OPAW library,
but also to be an interface of variables. Below is a list of variables that assumed to be imported
from the main code. If these variables are not in the main code, change them to be an internal variable.
Use pointers if necessary.
The term main_mod is just a generic name and should be changed as necessary.
The arrays from vloc_tot and below are also assumed to be allocated in the main program.

! 100% need to import
  use main_mod, only : nx,ny,nz,nn ! nn=nx*ny*nz, number of grid points !integer
  use main_mod, only : dx,dy,dz,dv ! grid volume elements !real*8
  use main_mod, only : ekcut       ! kinetic energy cutoff. If ekcut <= 0d0 then code will use  !real*8
  use main_mod, only : scale_vh    ! use Martyna-Tucker (scale_vh=2) or not (scale_vh=1) !integer
  use main_mod, only : nstates     ! total number of states !integer
  use main_mod, only : periodic    !periodic or not !logical
  use main_mod, only : funct       !0=lda, 1=pbe xc functional to be used !integer
  use main_mod, only : funct_c     !correlation functional for LIBXC !integer
  use main_mod, only : funct_x     !exchange functional for LIBXC !integer
                                   !   abinit subroutines to calculate ekcut.
! can be moved to be internal variables if needed
  use main_mod, only : rnel        ! # of electrons !real*8
  use main_mod, only : nocc        ! number of occupied states !integer
  use main_mod, only : vloc_tot  !shape (n) local part of the pseudopotential !real*8
  use main_mod, only : vk        !shape (n) related to applying the hartree potential !real*8
  use main_mod, only : ek        !shape (nx,ny,nz) kinetic energy in momentum space k^2/2 !real*8

-----------------------------------------------------------------
              2.2 Internal OPAW Variables
-----------------------------------------------------------------

Theres some internal parameters that generally do not need to be changed.
Feel free to play around with them if desired except for nk_loc

   real*8  :: p_fg=0.15d0   !parameter for the fine-grid (Ono-Hirose) feel free to change
                            !nfovnr=n_fine/n_rough = max(dx/p_fg,dy/p_fg,dz/p_fg)
   integer :: nk_loc=1      !do not change this. Required for gamma point calculation.
   real*8  :: rpad=1.0d0    !radial padding in calculating n_rough grid during ono-hirose
   real*8  :: rpad_r=2d0    !grid stuff
   real*8  :: tollsij=0.01  !see original OPAW DFT paper by Wenfei
   real*8  :: ek_factor=1d0  !for calculating ekcut from abinit PAW subroutines

=================================================================
              3. Structure of code
=================================================================

The library treats each Hamiltonian as an object. The file 4_ham.f90 has the ham_module
that contains details of the 'hamiltonian_obj' data type. Each instance of the
hamtilonian_obj contains the electronic density and potentials necessary to define a
Hamiltonian associated with a set of non-orthogonal pseudowavefunctions. The vloc_tot 
(the local part of the PAW potential) is a system constant and is not included in the obj.
The parts that define a differing hamiltonians is:
  dens     - non-orthogonal density
  dens_o   - orthogonal density
  nhat     - compensation charge
  v_h      - hartree potential
  v_xc     - exchange-correlation from libxc
  v_ks     = v_xc + v_h + vloc_tot
  rho^a_ij = sum_n <tilde \psi_n|p_i><p_j|tilde \psi_n>
  D^a_ij   - nonlocal terms that are made with rho^a_ij
  
Each instance of the hamiltonian_obj is initialized with the 'init_ham' subroutine and 
set with the 'opaw_make_hamiltonian' routines. See the 'OPAW Subroutines' chapter for more
details. 

Each hamiltonian_obj has an associated flag called 'h_type' that is declares what the 
hamiltonian is. For h_type=0 the associated hamiltonian is S^-1H, and for h_type=1 
the hamiltonian is S^(-1/2)HS^(-1/2).

All non-orthogonal and orthogonal wavefunctions in the code should be complex*16.

=================================================================
              4. OPAW Subroutines
=================================================================

The following are the main subroutines that you'll need to set up and use the OPAW library.
There are a few subroutines that I expect to be redundant if you are trying to implement this
library with any of the Neuhauser group's code so I renamed some subroutines such as vk_prep to
vk_prep_opaw to avoid conflict. Feel free to change if using a specialized vk_prep routine.
The routines 'opaw_libpaw_prepare' and 'get_rnel_opaw' can be called after the input file
is read and the system variables (nx,ny,nz,dx,...) are set.

Subroutines
 opaw_libpaw_prepare
 get_rnel_opaw
 init_ham
 opaw_make_hamiltonian
 rk4_prop_opaw
 sn_phi
 exx_expect_opaw

-----------------------------------------------------------------
              4.1 Main 
-----------------------------------------------------------------
Details

opaw_libpaw_prepare:
  Reads in the pawfile files, prepares the orthogonal projectors and matrix elements.
  It also prepares vk for applying the Coulomb potential and ek for applying the
  kinetic energy operator.


get_rnel_opaw
  description: calculates number of electrons (rnel) and occupied states (nocc) 

init_ham(nn,h_type,ham)
  description:
    This routine allocates nhat,dens,dens_o,vks,vxc,vh, and the atominfo%
   
  input: 
    integer :: nn          ! number of grid points nn=nx*ny*nz
    integer :: h_type      ! hamiltonian type =0 (S^-1H) =1 (S^-1/2HS^-1/2)

  output:
    hamiltonian_obj :: ham 

opaw_make_hamiltonian(nn,nocc,nstates,wfs,ham)
  description:
    This subroutine calculates the PAW density matrices and potentials.

  input: 
    integer     :: nn           ! number of grid points nn=nx*ny*nz
    integer     :: nocc         ! number of occupied states
    integer     :: nstates      ! number of states
    complex*16  :: wfs(nn,nocc) ! non-orthogonal (h_type 0) or 
                                ! orthogonal (h_type 1) wavefunctions
                              
  output:
    hamiltonian_obj :: ham 

rk4_prop_opaw(nn,nocc,nstates,dt,p,ham)
  description:
    Given a hamiltonian and time step, use 4th order Runge-Kutta to 
    propagate all the wavefunctions a single time step. If 
    ham%h_type=0 the wavefunctions are assumed to be the
    non-orthogonal and propagated with the S^-1 H hamiltonian
    and for ham%h_type=1 orthogonal wavefunctions with S^-1/2HS^-1/2.

  input: 
    integer     :: nn           ! number of grid points nn=nx*ny*nz
    integer     :: nocc         ! number of occupied states
    integer     :: nstates      ! number of states
    real*8      :: dt           ! time-step
    hamiltonian_obj :: ham      ! input hamiltonian

  output:
    complex*16  :: p(nn,nstates)   ! Propagated wavefunctions

sn_phi(pin,sp,n)
  description:
    Applies S^n to pin and returns the result as sp.

  input:
    integer     :: n               ! power in S^n
    complex*16  :: pin(nx,ny,nz)   ! input wavefunction

  output:
    complex*16  :: sp(nx,ny,nz)    ! sp = S^n pin


exx_expect_opaw(nn,nocc,nstates,psi_i,psi_j,psi_n,exx)
  description:
    Calculates the expectation value of the exchange fock exchange
    using the orthogonal pseudowavefunctions exx=<psi_i|V^hat_x|psi_j>.
 
  input: 
    integer     :: nn                ! number of grid points nn=nx*ny*nz
    integer     :: nocc              ! number of occupied states
    integer     :: nstates           ! number of states
    complex*16  :: psi_i(nn)         ! bra
    complex*16  :: psi_j(nn)         ! ket
    complex*16  :: psi_n(nn,nstates) ! orthogonal wavefunctions

  output:
    real*8      :: exx               ! expectation value

-----------------------------------------------------------------
              4.2 Overlapping subroutines
-----------------------------------------------------------------
I expect some subroutines to be redundant with other Neuhauser group
code. In particular:


=================================================================
              5. Compiling
=================================================================
Place this library in the main code you are trying to implement 
OPAW. I'm assuming you are using a Makefile to compile your program
and if there is not one, please make one. 

In the Makefile, define the make variables below

  opawlib = ./opawlib
  opaw_files = $(opawlib)/*.f90
  opawvloc_files = $(opawlib)/vloc/*.f90
  libpaw  = $(opawlib)/libpaw
  libpaw_  = $(libpaw)/0_m_libpaw_defs.o

Add the following rules below the all: rule

  $(libpaw_) :
    cd $(libpaw)  && $(CC) -c libpaw_libxc.c && $(FC) $(MPIFLG) -c *.F90 *.f90

  $(lib_) : #if 0_library_mpi_module.o does not exist run this rule
    cd $(libdir) && $(FC) $(MPIFLG) -c *.f90 *.f

  cleanlibpaw :
    touch $(libpaw)/a.o $(libpaw)/a.mod
    rm $(libpaw)/*.o $(libpaw)/*.mod

Then in the all rule compile $(libpaw_) first then $(lib_)

  all : $(lib_) $(libpaw_) main clean

See the examples Makefile to see an example of how to compile 
including how to compile with LIBXC.

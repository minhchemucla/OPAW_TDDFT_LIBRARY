

# OPAW LIBRARY 

This library is designed to implement our Orthogonal Projector Augmented Wave method into other electronic structure codes. Refer to Wenfei and Dannys' original [OPAW DFT paper](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.195118) and Minh, Tim, and Dannys' [OPAW TDDFT paper](https://pubs.aip.org/aip/jcp/article-abstract/160/14/144101/3281117/Time-dependent-density-functional-theory-with-the?redirectedFrom=fulltext) for theory. You can find both papers available [here](http://www.chem.ucla.edu/dept/Faculty/dxn/pages/publications.html) as well.

This library provides the necessary ingredients for 3 main functionalities:

  1. Apply the OPAW Hamiltonian.
  2. Perform RK4 time propagation.
  3. Calculate the expectation value of the exact exchange operator.

The code has the following limitations:

 

 - gamma-point (# k points = 1)
 - close-shelled (`nspin=1`)
 - LDA and PBE exchange-correlation functionals

 This library comes with a test program  with further explanation in the [test.md](docs/test.md) document. The library has its a parallelization scheme that will work for any number of cores. The parallelization scheme splits the occupied states over the number of cores, so I do not recommend running more cores than there are the number of occupied states. 


## Chapters

  1. [Prerequisite Libraries ](#prequisites) 
  2. [OPAW Input Files](#input)
  3. [How to use](#how_to)
  5. [Compiling](#compiling)  
  6. [Acknowledgements](#acknowledgements)

##      <a id="prequisites"></a>     1.   Prerequisite Libraries 
The following are scientific and mathematical libraries that are necessary. 

  1. **LIBXC**: Library for exchange-correlation functions. More information on how to install LIBXC can be found at this [link](https://libxc.gitlab.io/).

  2. **FFTW3** : Fastest Fourier Transform in the West version 3.x.x.

  3. **BLAS**: Basic Linear Algebra Subprograms.

  4. **LAPACK**: Linear Algebra Package.

The last three are typically available on supercomputers and Linux package tools such as `apt-get`.

##    <a id="input"></a>        2.     OPAW Input Files

The following input files are needed for the OPAW library to work. 

- **pawfiles:** This file should contain a list of elements and the corresponding PAW potential.   For example:

	    H     H.LDA_PW-JTH.xml
	    C     C.LDA_PW-JTH.xml

- **PAW potentials:** The PAW potentials in the pawfiles should be in the same directory as the pawfiles itself.
- **cnt.ini:** Coordinate file of the system in Bohr. An example of Naphthalene is provided in the `example` directory.


##  <a id="how_to"></a>   3. How to use

This section explains how to modify the code to use the library properly and the general procedures on how to use it. This section is an overview and the test program should be referred to as a sample of how to use the code.

***The OPAW wavefunctions in the code should be `complex*16` to avoid numerical issues.***

###  Variables interface

The `0_parameters.f90` file serves as an intermediate to store variables that are important in the OPAW library, but also to be an intermediate between variables in the library and main code that the library is being implemented into. Below is a list of variables that are assumed to be needed to be imported from the main code. If these variables are not in the main code, change them to be internal variables. Use pointers if necessary to make sure the names are spelled exactly as shown below. The module`main_mod` is just a generic name and should be changed as necessary. 

 
 #### Imported Variables
 
These variables have to be read in the main code and imported into the OPAW library:  

	use main_mod, only : nx,ny,nz,nn ! nn=nx*ny*nz, number of grid points !integer
	use main_mod, only : dx,dy,dz,dv ! grid volume elements !real*8
	use main_mod, only : ekcut       ! kinetic energy cutoff. If ekcut <= 0d0 then code will use  !real*8
	use main_mod, only : scale_vh    ! use Martyna-Tucker (scale_vh=2) or not (scale_vh=1) !integer
	use main_mod, only : nstates     ! total number of states !integer
	use main_mod, only : periodic    !periodic or not !logical
	use main_mod, only : funct       !0=lda, 1=pbe xc functional to be used !integer
	
	
The following variables can be made into internal variables if they are not needed in the main code.
	
	use main_mod, only : funct_c     !correlation functional for LIBXC !integer
	use main_mod, only : funct_x     !exchange functional for LIBXC !integer
	use main_mod, only : rnel        ! # of electrons !real*8
	use main_mod, only : nocc        ! number of occupied states !	integer
	use main_mod, only : vloc_tot  !shape (n) local part of the pseudopotential !real*8
	use main_mod, only : vk        !shape (n) related to applying the hartree potential !real*8

#### Internal OPAW Library Variables
Some internal parameters generally do not need to be changed. Feel free to play around with them if desired except for nk_loc

	real*8               :: p_fg=0.15d0   !parameter for the fine-grid (Ono-Hirose) feel free to change
                          !nfovnr=n_fine/n_rough = max(dx/p_fg,dy/p_fg,dz/p_fg) 
	integer              :: nk_loc=1      !do not change this. Required for gamma point calculation.
	real*8               :: rpad=1.0d0     !radial padding in calculating n_rough grid during ono-hirose
	real*8  	     :: rpad_r=2d0     !grid stuff
	real*8               :: tollsij=0.01   !see original OPAW DFT paper by Wenfei
	real*8               :: ek_factor=1d0  !for calculating ekcut from abinit PAW subroutines
	real*8, allocatable  :: ek3d(:,:,:)    !3d exp(-k^2)

#### OPAW Hamiltonian

The library treats the PAW Hamiltonian using a Fortran data_type. This is equivalent to using an object in python or a struct in C. The file `opawlib/4_ham.f90` has the `opaw_ham_mod` that contains details of the `opaw_ham_obj` data type. Each instance of the `opaw_ham_obj` contains the electronic density and potentials necessary to define a Hamiltonian associated with a set of non-orthogonal pseudowavefunctions. The local part of the PAW potential $V_{loc}^{tot}$ is assumed to be system constant. 

The parts that define a PAW Hamiltonian are:
 - $n(r)$         - The PAW density
 - $\hat{n}(r)$   -The compensation charge
 - $V_{H}$        - hartree potential
 - $V_{XC}$       - exchange-correlation from libxc
 - $V_{KS}$       = $V_{XC}$ + $V_{H}$ + $V_{loc}^{tot}$
 - $\rho^{a}_{ij}$ - PAW density matrix $=\Sigma_n \braket{\tilde{\psi}_n|p_i}\braket{ p_j|\tilde{\psi}_n}$: 
 - $D_{ij}^{a}$   - nonlocal terms that are made with $\rho^a_{ij}$


##   <a id="compiling"></a> 4. Compiling

The `opawlib` directory has the following subdirectories:
	
 - `vloc`
 - `libpaw`
 - `lib`
 - `XCI`

The latter two directories are designed for non-Neuhauser group users. The `vloc` directory contains routines to calculate the local part of the PAW potential and `libpaw` contains many Abinit Place this library in the main code you are trying to implement  OPAW. In the Makefile, define the make variables below

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

Then in the all rule compile `$(libpaw_)` first then `$(lib_)`

	all : $(lib_) $(libpaw_) main clean

See the `examples/makefile` to see an example of how to compile including how to compile with LIBXC.

##   <a id="acknowledgements"></a> 5.  Acknowledgements

This code is supported by the U.S. Department of Energy, Office of Science, Office of Advanced Scientific Computing Research and Office of Basic Energy Sciences, Scientific Discovery through Advanced Computing (SciDAC) program under Award Number DE-SC0022198.

Much of the PAW capabilities were adapted from the open source [ABINIT](https://www.abinit.org/) software version 8.0.8. The PAW algorithm in ABINIT was based on the work of [Torrent, Marc, et al. "_Implementation of the projector augmented-wave method in the ABINIT code: Application to the study of iron under pressure_" Computational Materials Science **42**, 337-351 (2008)](https://doi.org/10.1016/j.commatsci.2007.07.020).


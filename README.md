


# OPAW LIBRARY 

This FORTRAN library was designed to implement our Orthogonal Projector Augmented Wave method into other electronic structure codes. The library is built on representing wavefunctions in a 3D rectangular box. Refer to Wenfei and Dannys' original [OPAW DFT paper](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.195118) and Minh, Tim, and Dannys' [OPAW TDDFT paper](https://pubs.aip.org/aip/jcp/article-abstract/160/14/144101/3281117/Time-dependent-density-functional-theory-with-the?redirectedFrom=fulltext) for theory. You can find both papers available [here](http://www.chem.ucla.edu/dept/Faculty/dxn/pages/publications.html) as well.

This library provides the necessary ingredients for 3 main functionalities:

 - Apply the OPAW Hamiltonian.
 - Perform RK4 time propagation.
 - Calculate the expectation value of the exact exchange operator.

The package contains 4 directories.
 - [`docs`](docs) : Documentation files.
 - [`example_src`](docs/example.md) : A small program to highlight how to use the OPAW library.
 - [`opawlib`](opawlib): The library source files. 
 - [`example`](example): A Naphthalene example using the program compiled in `example_src`

The code has the following features:

 - gamma-point (# k points = 1).
 - close-shelled (`nspin=1`).
 - LDA and PBE exchange-correlation functionals.
 - Periodic and nonperiodic ([Martyna-Tuckermann](https://pubs.aip.org/aip/jcp/article-abstract/110/6/2810/474725/A-reciprocal-space-based-method-for-treating-long?redirectedFrom=fulltext) method)

# Chapters

  1. [Prerequisites](#prequisites) 
  2. [OPAW Input Files](#input)
  3. [How to use](#how_to)
  5. [Compiling](#compiling)  
  6. [Acknowledgements](#acknowledgements)

##      <a id="prequisites"></a>     1.   Prerequisites
This library is written in FORTRAN with MPI parallelization. Compile this code using the `mpif90` compiler wrapper or any other  MPI wrapper. There is a `libpaw_libxc.c` C file in the `opawlib/libpaw` directory so a C compiler like `gcc` will be need too.

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



This section explains how to modify the code to use the library properly and the general procedures on how to use it. This section is an overview and the example program should be referred to as a sample of how to use the code.

***The OPAW wavefunctions in the code should be `complex*16` to avoid numerical issues.***

### 3.1 Import variables and arrays into the library.

The `0_parameters.f90` file stores important variables in the OPAW library and functions as an interface to import variables and arrays from the main code. Use pointers if necessary so that the variables and arrays are spelled exactly as shown below. For example if `ekcut` is `ecut` in the main code, use `use main_mod, only : ekcut=>ecut`. The module `main_mod` is just a placeholder name. 
 
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

Some internal parameters generally do not need to be changed. Feel free to play around with them if desired except for nk_loc

	real*8               :: p_fg=0.15d0   !parameter for the fine-grid (Ono-Hirose) feel free to change
                          !nfovnr=n_fine/n_rough = max(dx/p_fg,dy/p_fg,dz/p_fg) 
	integer              :: nk_loc=1      !do not change this. Required for gamma point calculation.
	real*8               :: rpad=1.0d0     !radial padding in calculating n_rough grid during ono-hirose
	real*8  	     :: rpad_r=2d0     !grid stuff
	real*8               :: tollsij=0.01   !see original OPAW DFT paper by Wenfei
	real*8               :: ek_factor=1d0  !for calculating ekcut from abinit PAW subroutines
	real*8, allocatable  :: ek3d(:,:,:)    !3d exp(-k^2)

### 3.2 OPAW Hamiltonian

The library represents a PAW Hamiltonian using Fortran [derived types](https://fortran-lang.org/en/learn/quickstart/derived_types/). Fortran derived types are very similar to structs in C. The file `opawlib/4_ham.f90` contains the `opaw_ham_mod` that defines the `opaw_ham_obj` derived type. Each instance of the `opaw_ham_obj` type contains the electronic density and potentials necessary to define a PAW Hamiltonian. The local part of the PAW potential $V_{loc}^{tot}$ is assumed to be system constant and is not stored in the derived type.

The parts that define a unique PAW Hamiltonian for a system are:
 - $n(r)$         - The PAW density
 - $\hat{n}(r)$   -The compensation charge
 - $V_{H}$        - hartree potential
 - $V_{XC}$       - exchange-correlation from libxc
 - $V_{KS}$       = $V_{XC}$ + $V_{H}$ + $V_{loc}^{tot}$
 - $\rho^{a}_{ij}$ - PAW density matrix $=\Sigma_n \braket{\tilde{\psi}_n|p_i}\braket{ p_j|\tilde{\psi}_n}$: 
 - $D_{ij}^{a}$   - nonlocal terms that are made with $\rho^a_{ij}$

In the main code, you can create an instance of the `opaw_ham_obj` for example

    subroutine generic_name
      use opaw_ham_mod
      ...
      type(opaw_ham_obj) :: ham

Then given a set of OPAW wavefunctions you can use the `opaw_make_ham`routine to calculate the PAW potentials. Once the `ham` instance is prepared, you can use the `paw_ham` and `opaw_ham` subroutines to apply the PAW and OPAW Hamiltonian on a wavefunction. Additionally, to propagate it in time you can use the `rk4_prop_opaw` routine inputting the occupied states to perform a single time step under `ham`. Further details of these subroutines and how to use them can be found in [example.md](docs/example.md) and [subroutines.md](docs/subroutines.md).

***This library is parallelized over the occupied states so I do not recommend running more cores than the number of occupied states in the system.***

##   <a id="compiling"></a> 4. Compiling

The `opawlib` directory has the following subdirectories:
	
 - `vloc`
 - `libpaw`
 - `lib`
 - `XCI`

The latter two directories are for non-Neuhauser group users as the `opawlib/lib` directory has only the files necessary to make the OPAW library work while Neuhauser group members will have our own `lib` with more files. The `opawlib/vloc` directory contains routines to calculate the local part of the PAW potential and `libpaw` contains many ABINIT subroutines. The `XCI` directory contains subroutines to interface with LIBXC. Place this library in the source code you are trying to implement OPAW. 

A `makefile` in the `example_src` shows an example of how to compile and is recommended to be used as a template. In the compilation of the example program, there are 3 steps:

 1. Compile the ABINIT subroutines in `opawlib/libpaw`.
 2. Compile the Neuhauser library in `opawlib/lib`.
 3. Compile the main code with the library

In the Makefile, define the make variables below

	opawlib = ./opawlib
	opaw_files = $(opawlib)/*.f90
	opawvloc_files = $(opawlib)/vloc/*.f90
	libpaw  = $(opawlib)/libpaw
	libpaw_  = $(libpaw)/0_m_libpaw_defs.o
	libdir = ./lib
	lib_ =  $(libdir)/0_library_mpi_module.o

With these variables define the following rules

	$(libpaw_) :
	  cd $(libpaw)  && $(CC) -c libpaw_libxc.c && $(FC) $(MPIFLG) -c *.F90 *.f90

	$(lib_) : #if 0_library_mpi_module.o does not exist run this rule
	  cd $(libdir) && $(FC) $(MPIFLG) -c *.f90 *.f
	
	cleanlib :
	  touch $(libdir)/a.o $(libdir)/a.mod
	  rm lib/*.o lib/*.mod
	
	cleanlibpaw :
	  touch $(libpaw)/a.o $(libpaw)/a.mod
	  rm $(libpaw)/*.o $(libpaw)/*.mod

The `libpaw_` and `lib_` variables are there to tell GNU Make to not compile `libpaw` and `lib` if `0_m_libpaw_defs.o` and `0_library_mpi_module.o` are present i.e. if the object files in `libpaw` and `lib` have been compiled already. If you modify any files in `libpaw` or `lib`, you must run `make cleanlibpaw` and `make cleanlib` and then recompile for those changes to take place.

Shown below are other Make variables defining the MPI Fortran and C compilers, LIBXC, FFTW2, BLAS, and LAPACK  libraries, and some MPI flags.

	FC    = mpif90  
	CC    = gcc
	LIBXC = -I/home/minh/codes/LIBXC/include -L/home/minh/codes/LIBXC/lib -lxcf90 -lxc
	FFTFLG  = -lfftw3 -lfftw3f
	BLASFLG   = -lblas
	LAPACKFLG   = -llapack
	libs = $(FFTFLG) $(BLASFLG) $(LAPACKFLG) $(LIBXC)
	MPIFLG  = -DMPI -O3 # -g -fcheck=all -fbacktrace

The `-g`, ` -fcheck=all`, and `-fbacktrace` are debugging flags. 

To compile the library with the main program, base your `main` rule on the following example:

	main :
	  $(FC) -I $(libdir) $(libdir)/*.o -I $(libpaw) $(libpaw)/*.o \
	  $(param_file) $(opaw_files) $(opawvloc_files) $(src_files) $(XCI)\
	  -o $(output) $(libs)

Then define the following `all` rule to compile.

	all : $(lib_) $(libpaw_) main clean

See `example_src/makefile` to see a full example of a makefile.

##   <a id="acknowledgements"></a> 5.  Acknowledgements

This code is supported by the U.S. Department of Energy, Office of Science, Office of Advanced Scientific Computing Research and Office of Basic Energy Sciences, Scientific Discovery through Advanced Computing (SciDAC) program under Award Number DE-SC0022198.

Much of the PAW capabilities were adapted from the open source [ABINIT](https://www.abinit.org/) software version 8.0.8. The PAW algorithm in ABINIT was based on the work of [Torrent, Marc, et al. "_Implementation of the projector augmented-wave method in the ABINIT code: Application to the study of iron under pressure_" Computational Materials Science **42**, 337-351 (2008)](https://doi.org/10.1016/j.commatsci.2007.07.020).

This library was created and authored by Minh Nguyen based on the OPAW-DFT work done by Wenfei Li and Daniel Neuhauser and the OPAW-TDDFT work of Minh Nguyen, Tim Duong, and Daniel Neuhauser.

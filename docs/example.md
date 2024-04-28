# Example program

Included with this repository is an example directory where the program is run and an `example_src` directory where the source files are. This file will walk through both to highlight how the `opawlib` should be used. The example contains a `wf_bar.txt` file which is the output of a non-periodic LDA Naphthalene OPAW-DFT calculation. I am currently working on uploading the OPAW-TDDFT program on GitHub so you can generate your own OPAW-DFT output. To compile the program, refer to the [README.md](../README.md) document.


## example
This Naphthalene example contains:

 - `cnt.ini`: coordinate file of Naphthalene in Bohr.
 - `C.LDA_PW-JTH.xml` and `H.LDA_PW-JTH.xml`: PAW potentials for C and H.
 - `input`: Input file for the program. 
 - `pawfiles`: List of PAW potentials to be read in.
 - `wf_bar_lda_nonperiodic.bin`:  Nonperiodic LDA OPAW-DFT wavefunctions for Naphthalene in binary.
 - `wf_bar_lda_periodic.bin`: Periodic LDA OPAW-DFT wavefunctions for Naphthalene in binary.
 - `wf_bar_pbe_nonperiodic.bin`: Nonperiodic PBE OPAW-DFT wavefunctions for Naphthalene in binary.
 - `wf_bar_pbe_periodic.bin`: Periodic PBE OPAW-DFT wavefunctions for Naphthalene in binary.
 - `wf_bar.bin`: A symlink to OPAW-DFT wavefunctions for Naphthalene in binary.
 - `opaw_test.x`: A symlink to the program in `example_src`.

The example also includes 4 example output files: `output_lda_nonperiodic`, `output_pbe_nonperiodic`, `output_lda_periodic`, and `output_pbe_periodic` which were calculated using the corresponding wf_bar.bin files. These outputs are for you to compare to your compiled version of the example program. Change the symlink using `ln -s` to change what the `wf_bar.bin` file is pointing to.

Other PAW potentials are found on the [ABINIT](https://www.abinit.org/psp-tables) website. Once the `opaw_test.x` program is compiled in `example_src`, the symlink will be valid.

Assuming the program is compiled with mpif90, the program can be run with the following: 

	mpirun -n 1 ./opaw_test.x

## Algorithm in example_src
The program executes the following steps
```
 1. call prepare_mpi_lib           Prepare MPI using a library.
 2. call read_input                Read `input` for system parameters.
 3. call read_wfs                  Read in the OPAW wavefunctions from `wf_bar.bin`
 4. call prepare_opaw              Prepare OPAW-related terms.
 5. call opaw_make_ham             Calculates the PAW Hamiltonian from the OPAW wavefunctions.
 6. call exx_expect_opaw           Calculates $\braket{\psi_{nocc}|V_x|\psi_{nocc}}$.
 7. call tddft                     Performs OPAW-TDDFT using RK4.
```
In the OPAW-TDDFT algorithm, which is contained in the `tddft.f90` file, the OPAW wavefunctions are perturbed:

$$ \psi'_n(r)=e^{-i\gamma \hat{r}}\psi_n(r)$$

where $\hat{r}=x,y,z$ corresponding with `ipol=1,2,3` and $\gamma$ is the `strength` input in the `input` file. For each time-step, the PAW Hamiltonian is calculated then the RK4 propagation is performed on all the wavefunctions. In each time step the dipole-dipole correlation function, $\Delta n_i(r,t)$, is calculated from the OPAW wavefunctions and then printed out. 

$$\Delta n_i(r,t)=\frac{1}{\gamma}\Big(n^{\gamma}(r,t)-n^{\gamma=0}(r,t)\Big)$$

For more detail on $\Delta n_i(r,t)$ and how it is used to calculate absorption spectra see the [OPAW-TDDFT paper](https://pubs.aip.org/aip/jcp/article-abstract/160/14/144101/3281117/Time-dependent-density-functional-theory-with-the?redirectedFrom=fulltext). 

After running the code, some `fort.xxx` files will be generated that contain debugging info from the ABINIT subroutines.

## input file

Below is an input file highlight the parameters that should be imported into the OPAW library. Do not change any parameters from `dz` and above or else the example program will fail.

	nx       56             !number of grid points in x direction
	ny       52             !number of grid points in y direction
	nz       32             !number of grid points in z direction
	dx       0.5d0          !grid spacing in x direction
	dy       0.5d0          !grid spacing in y direction
	dz       0.5d0          !grid spacing in z direction
	periodic F              !Use periodic or non-periodic (Martyna-Tuckermann approach)
	funct    0              !0:pwlda, 1:pbe
	ekcut    -1d0           ! Kinetic energy cutt off. Negative values will use ABINIT PAW subroutines to calculate ekcut.
	dt       0.05d0         ! Time-step in atomic units
	nt       10             ! Number of time steps
	strength 1e-3           ! gamma
	exct_pol 1              ! Direction of perturbation



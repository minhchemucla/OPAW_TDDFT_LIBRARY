## compiler options
FC   	= mpif90 -O3 -g -fcheck=all -fbacktrace #-Wall
CC		= gcc 
LIBXC = -I/home/minh/codes/LIBXC/include -L/home/minh/codes/LIBXC/lib -lxcf90 -lxc
FFTFLG 	= -lfftw3 -lfftw3f
BLASFLG 	= -lblas
LAPACKFLG 	= -llapack
libs = $(FFTFLG) $(BLASFLG) $(LAPACKFLG) $(LIBXC)
MPIFLG  = -DMPI  

XCI  = XCI/*f90   # XCI includes a few interface files to LIBXC
opawlib = ./opawlib
opaw_files = $(opawlib)/*.f90
opawvloc_files = $(opawlib)/vloc/*.f90
libpaw  = $(opawlib)/libpaw
libpaw_  = $(libpaw)/0_m_libpaw_defs.o
libdir = ./lib
lib_ =  $(libdir)/0_library_mpi_module.o
output = opaw_test.x
src_files = $(filter-out 0_main.f90, $(sort $(wildcard *.f90)))
param_file = 0_main.f90

# compile
.PHONY: all clean cleanlib cleanx cleanlibpaw cleanall
all : $(lib_) $(libpaw_) main clean

$(libpaw_) : #if 0_m_libpaw_defs.o does not exist run this rule
	cd $(libpaw)  && $(CC) -c libpaw_libxc.c && $(FC) $(MPIFLG) -c *.F90 *.f90 

$(lib_) : #if 0_library_mpi_module.o does not exist run this rule
	cd $(libdir) && $(FC) $(MPIFLG) -c *.f90 *.f

main :
	$(FC) -I $(libdir) $(libdir)/*.o -I $(libpaw) $(libpaw)/*.o \
  $(param_file) $(opaw_files) $(opawvloc_files) $(src_files) $(XCI)\
  -o $(output) $(libs) 

#	$(FC) -c -I $(LIB_OPAW) 1_libpaw_mod.f90
#	$(FC) -I ./lib ./lib/*.o -I $(LIB_OPAW) $(LIB_OPAW)/*.o   *.f90 vloc/*.f90 $(OPAW)/*.f90 -o opaw.x \
#        $(XCID) $(FFTFLG) -lblas -llapack $(LIBXC) 

cleanall : clean cleanx cleanlib cleanlibpaw
      
cleanx :
	touch a.x
	rm *.x

clean :
	touch empty.o empty.mod
	rm *.o *mod 

cleanlib :
	touch $(libdir)/a.o $(libdir)/a.mod
	rm lib/*.o lib/*.mod

cleanlibpaw :
	touch $(libpaw)/a.o $(libpaw)/a.mod
	rm $(libpaw)/*.o $(libpaw)/*.mod

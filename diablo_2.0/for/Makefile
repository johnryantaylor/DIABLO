# This is the makefile for diablo.
# To compile the code, just type make.  Such an approach makes
# recompilation of the code easy, recompiling only as necessary
# to account for recent changes to the code.
#
# As the user, set the following definitions:

#**********   User Defined Variables Below *********
# Fortran 90 complier to use:
COMPILER = ifort

# Any desired options for the compiler (e.g. -O2, -g, etc.)
USEROPTS =-fpp -u
#USEROPTS = -O3 -fpp -mcmodel=medium -u

# Location where fftw and netcdf (optional) libraries are installed
LINKDIR = /mhome/damtp/atmos/jrt51/local/lib

# Location where the optional netcdf include file (netcdf.inc) is installed
INCLUDEDIR = /usr/local/include

# Option to compile with MPI libraries
PARALLEL = TRUE

# Option to enable the LES model (loads the required variables into memory)
LES = TRUE

# Option to complile with the HDF5 libarary
HDF5 = TRUE
# **********    END of user definitions ************

ifeq ($(ARCH),ifort)
COMPILER = ifort
USEROPTS = -O3 -fpp
else ifeq ($(ARCH),gfortran)
COMPILER = gfortran
USEROPTS = -O3 -cpp
endif

ifeq ($(LES),TRUE)
LES_o = les.o
else
LES_o = no_les.o
endif

ALL2ALL=1

ifeq ($(PARALLEL),TRUE)
COMPILER = mpif90 
MPI = mpi.o
MPIF = mpi.f
ifeq ($(HDF5),TRUE)
HDF5_o = hdf5.o 
COMPILER = h5pfc
HDF5OPTS=-DHDF5 
endif
else
ifeq ($(HDF5),TRUE)
HDF5_o = hdf5s.o 
COMPILER = h5fc
HDF5OPTS=-DHDF5
endif
MPI = mpi_serial.o
endif

MAIN = diablo.f
HEADER = header

# Use the parameters to set flags
ifeq ($(NETCDF),TRUE)
COMPOPTS = $(USEROPTS) $(HDF5OPTS) -I$(INCLUDEDIR)
LINKOPTS = -L$(LINKDIR) -lrfftw -lfftw -lnetcdf
NETCDF_o = netcdf.o
else
COMPOPTS = $(USEROPTS) $(HDF5OPTS)
LINKOPTS = -L$(LINKDIR) -lrfftw -lfftw
NETCDF_o = no_netcdf.o
endif

diablo: $(MAIN) diablo_io.o periodic.o channel.o courant.o $(LES_o) \
	duct.o cavity.o fft.o rand.o set_ics.o save_stats.o filter.o user_rhs.o $(MPI) \
	$(HEADER) grid_def $(HDF5_o)
	$(COMPILER) $(COMPOPTS) $(MAIN) -o diablo \
	diablo_io.o periodic.o channel.o courant.o $(LES_o) \
	duct.o cavity.o fft.o rand.o set_ics.o save_stats.o filter.o user_rhs.o \
	$(MPI) $(LINKOPTS) $(HDF5_o)

diablo_io.o: diablo_io.f header grid_def
	$(COMPILER) $(COMPOPTS) -c diablo_io.f

periodic.o: periodic.f fft.o header grid_def
	$(COMPILER) $(COMPOPTS) -c periodic.f

channel.o: channel.f fft.o $(MPI) header grid_def
	$(COMPILER) $(COMPOPTS) -c channel.f

set_ics.o: set_ics.f channel.o fft.o $(MPI) header grid_def
	$(COMPILER) $(COMPOPTS) -c set_ics.f

save_stats.o: save_stats.f channel.o fft.o $(MPI) header grid_def
	$(COMPILER) $(COMPOPTS) -c save_stats.f

filter.o: filter.f channel.o fft.o $(MPI) header grid_def
	$(COMPILER) $(COMPOPTS) -c filter.f

user_rhs.o: user_rhs.f channel.o fft.o $(MPI) header grid_def
	$(COMPILER) $(COMPOPTS) -c user_rhs.f

courant.o: courant.f channel.o fft.o $(MPI) header grid_def
	$(COMPILER) $(COMPOPTS) -c courant.f

ifeq ($(LES),TRUE) 
les.o: les.f fft.o header header_les grid_def
	$(COMPILER) $(COMPOPTS) -c les.f
else
no_les.o: dummy_code/no_les.f
	$(COMPILER) $(COMPOPTS) -c dummy_code/no_les.f
endif

ifeq ($(PARALLEL),TRUE)
mpi.o: $(MPIF) header 
	$(COMPILER) $(COMPOPTS) -c -o mpi.o $(MPIF) 
else
mpi_serial.o: dummy_code/mpi_serial.f header grid_def
	$(COMPILER) $(COMPOPTS) -c dummy_code/mpi_serial.f
endif

header : header_mpi grid_def

header_mpi : grid_mpi

hdf5.o : hdf5.f
	$(COMPILER) $(COMPOPTS) -c hdf5.f

hdf5s.o : hdf5s.f
	$(COMPILER) $(COMPOPTS) -c hdf5s.f

duct.o: duct.f header grid_def
	$(COMPILER) $(COMPOPTS) -c duct.f

cavity.o: cavity.f header grid_def
	$(COMPILER) $(COMPOPTS) -c cavity.f

fft.o:  fft.f header grid_def
	$(COMPILER) $(COMPOPTS) -c fft.f

rand.o:  rand.f
	$(COMPILER) $(COMPOPTS) -c rand.f

clean:
	rm -f *.o fort.* *~ diablo core

# Compiler specific notes:
#
# Compilation with Absoft Linux Fortran 77 appears to be impossible, as it
# cannot handle the INTEGER*8 option required by FFTW.  If someone finds
# a way around this, please let me know.
# 
# Compilation with Absoft Linux Fortran 90 is possible, but the option
# -YEXT_NAMES=LCS must be used as one of the link options so the compiler
# can find the lowercase external library function names.
#
# Compilation with Lahey Fortran 95 (lf95) is possible, but there is an
# underscore incompatability with the FFTW libraries, which are compiled
# with g77.  To get around this, you need to go into fft.f and add 
# trailing underscores to the name of every fftw function where they
# appear throughout the code.


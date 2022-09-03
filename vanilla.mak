#
#  Default build configuration file. 
#
#  See README.txt for further instructions.
#
BUILD_ID :="Plain vanilla, built on $(shell hostname) at $(shell date)"
  ACT = sed -e 's/^!\*nq/    /'                     # Disable quad-math statements
# ACT = sed -e 's/^!\*qd/    /' -e 's/^!\*lq/    /' # Enable quad-math statements everywhere
 ACT2 =                                             # Disable FFTW3
#ACT2 = -e 's/^!\*ft/    /'                         # Enable FFTW3 calls
#  Fortran compiler. These optimization flags are most likely grossly inadequate;
#  see examples in the configs/ subdirectory for suggestions.
F90 = gfortran -O -fopenmp -cpp -I. -D__BUILD_ID__='$(BUILD_ID)'
#  Fortran linker
F90L = $(F90) 
#  BLAS and LAPACK libraries
#  If you enable OpenMP parallel execution, the libraries you supply MUST support
#  multi-threaded execution as well!
 LAPACK = -llapack -lblas
#LAPACK = -llapack -lblas ~/lib64/libquadlapack_gfortran.a
#  Any additional support libraries
 LIBEXTRA = 
#LIBEXTRA = -lfftw3_omp -lfftw3f_omp -lfftw3 -lfftw3f

BUILD_ID :="Optimized gfortran, built on $(shell hostname) at $(shell date)"
  ACT = sed -e 's/^!\*nq/    /'                     # Disable quad-math statements
# ACT = sed -e 's/^!\*qd/    /' -e 's/^!\*lq/    /' # Enable quad-math statements everywhere
#ACT2 =                                             # Disable FFTW3
 ACT2 = -e 's/^!\*ft/    /'                         # Enable FFTW3 calls
# Add m64 -mavx -mavx2 etc non-native instruction sets are wanted.
F90 = gfortran-13 -I. \
      -O3 -flto -fprotect-parens -march=native -mtune=native -fopenmp \
      -ffast-math -fcx-fortran-rules -mrecip \
      -fno-realloc-lhs -fbacktrace -g \
      -static-libgfortran \
      -cpp -D__BUILD_ID__='$(BUILD_ID)' -Wall
#
F90L = $(F90) 
#  BLAS and LAPACK libraries
#  If you enable OpenMP parallel execution, the libraries you supply MUST support
#  multi-threaded execution as well!
 LAPACK = -llapack -lblas
#LAPACK = -llapack -lblas ~/lib64/libquadlapack_gfortran.a
#  Any additional support libraries
 LIBEXTRA = -lfftw3_omp -lfftw3f_omp -lfftw3 -lfftw3f

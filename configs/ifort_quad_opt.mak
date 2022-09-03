BUILD_ID :="Optimized ifort (quad), built on $(shell hostname) at $(shell date)"
#
# ACT = sed -e 's/^!\*nq/    /'                     # Disable quad-math statements
  ACT = sed -e 's/^!\*qd/    /' -e 's/^!\*lq/    /' # Enable quad-math statements everywhere
#ACT2 =                                             # Disable FFTW3
 ACT2 = -e 's/^!\*ft/    /'                         # Enable FFTW3 calls
# -xHost blows on a Threadripper
F90 = ifort \
            -qopenmp -qmkl=sequential -warn -assume buffered_io -assume protect_parens -heap-arrays 32 \
            -Ofast -ipo8 -xAVX2 -complex_limited_range -fp-model fast=1 -ftz \
            -debug full -debug extended -traceback -static-intel \
            -cpp -D__BUILD_ID__='$(BUILD_ID)'
#
F90L = $(F90) 
#  BLAS and LAPACK libraries - use MKL where possible
 LAPACK = ~/lib64/libquadlapack_intel.a
#  Any additional support libraries
 LIBEXTRA = -lfftw3_omp -lfftw3f_omp -lfftw3 -lfftw3f

.PHONY: goal clean
goal: makefile.dep
	  make hydrogen_tunnel_v2.x
	  make general_tunnel.x

MAKEFLAGS = -r

.SUFFIXES: .f90 .o .x .c .dep

#
#  This is the default; a config file may override it.
#
  ACT = sed -e 's/^!\*nq/    /'                     # Disable quad-math statements
# ACT = sed -e 's/^!\*qd/    /' -e 's/^!\*lq/    /' # Enable quad-math statements everywhere
 ACT2 =                                             # Disable FFTW3
#ACT2 = -e 's/^!\*ft/    /'                         # Enable FFTW3 calls
#
# System-specific overrides
#
  include vanilla.mak
# include configs/ifort_opt.mak
# include configs/ifort_quad_opt.mak
# include configs/gfortran_opt.mak
# include configs/gfortran_quad_opt.mak

#
# Finish the set-up
#
LIBS = $(LAPACK) $(LAPACK) $(LIBEXTRA)

#
# Compiling and archiving rules
#
.f90.o:
	$(ACT) $(ACT2) $< >preprocess/$<
	$(F90) -c preprocess/$<

dgefa.o:        dgefa.f
	$(F90) -c dgefa.f

dgedi.o:        dgedi.f
	$(F90) -c dgedi.f

clean:
	-/bin/rm -f *.{o,mod,x,il,a} *__genmod.f90 checkpoint_{field,main}.* makefile.dep *.optrpt ./preprocess/*.f90

makefile.dep: $(shell echo *.f90)
	./make-depend.sh $^ > $@

#
# Explicit dependencies
#

LIBHYDROGEN += accuracy.o
LIBHYDROGEN += constants.o
LIBHYDROGEN += derivative_tools.o
LIBHYDROGEN += dgedi.o
LIBHYDROGEN += dgefa.o
LIBHYDROGEN += find_minimum.o
LIBHYDROGEN += fftw.o
LIBHYDROGEN += find_root.o
LIBHYDROGEN += lapack.o
LIBHYDROGEN += math.o
LIBHYDROGEN += poly_tools.o
LIBHYDROGEN += sort_tools.o
LIBHYDROGEN += timer.o
LIBHYDROGEN += tridiagonal_tools.o
LIBHYDROGEN += tridiagonal_cmtql1.o
LIBHYDROGEN += versions.o
LIBHYDROGEN += general_tunnel_bound.o
LIBHYDROGEN += general_tunnel_continuum.o
LIBHYDROGEN += general_tunnel_asymptotic.o
LIBHYDROGEN += general_tunnel_data.o
LIBHYDROGEN += general_tunnel_dump.o
LIBHYDROGEN += general_tunnel_nonadiabatic.o
LIBHYDROGEN += general_tunnel_potential.o

#
# Building the binaries
#
hydrogen_tunnel_v2.x: hydrogen_tunnel_v2.o $(LIBHYDROGEN)
	$(F90) -o hydrogen_tunnel_v2.x hydrogen_tunnel_v2.o $(LIBHYDROGEN) $(LIBS)

general_tunnel.x: general_tunnel.o $(LIBHYDROGEN)
	$(F90) -o general_tunnel.x general_tunnel.o $(LIBHYDROGEN) $(LIBS)

#
# Automatically-generated dependencies
#
include makefile.dep

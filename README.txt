exit 1;

Last updated: 2022 September 3rd
------------

Contents
========

0. Preamble and references
1. Installation pre-requisites
2. Installation
3. Advanced installation
4. Input
4.1. hydrogen_tunnel_v2.x
4.2. general_tunnel.x

Static-field tunneling in a central potential
=============================================

The canonical location for the source code is at:

https://github.com/MBI-Theory/hydrogen-tunnel.git

Please report issues, ask questions, and make feature requests through the
github issue interface. Alternatively, you can reach the main author at:

Serguei Patchkovskii, Serguei.Patchkovskii@mbi-berlin.de

0. Preamble and references

TO BE ADDED

1. Installation pre-requisites

Building the codes requires a reasonably modern Fortran compiler (Fortran 95,
with some Fortran 2003 extensions), and a working LAPACK library. GNU Fortran 7
through 12 and Intel Fortran 2021.5.0 are known to work. LAPACK 3.5.0 is known
to work. Most later versions should also work.

Building for parallel execution requires OpenMP (2.0 or better) support, and
LAPACK/BLAS built with re-entrancy support.

WARNING: Some vendor-supplied LAPACK produce incorrect results or hang 
WARNING: when used inside an OpenMP parallel region.

Quadruple-precision compilation requires a special, non-standard version
of the LAPACK library (see "Advanced installation" below).

Some optional analysis features depend on the FFTW3 library being available.

2. Installation

By default, all system-specific configuration variables controlling build
process are assembled in the file "vanilla.mak", found in the top-level
installation directory. Additional configuration sets are found in configs/
subdirectory; see "Advanced installation" below.

A build configuration file contains the following variables:

 BUILD_ID - Character string identifying a specific build; it will be printed
            each time the resulting binary is executed. See subroutine "start"
            in "hydrogen_tunnel_v2.f90" and "general_tunnel.f90" for more 
            details.

 ACT and
 ACT2     - Defines pre-processing commands for the specific variant of the
            code. The main possibilities are building for double or quadruple 
            precision, or enabling FFTW calls.

            Activated sources are kept in the ./preprocess/ subdirectory.

 F90      - A command line used to invoke Fortran compiler, including all 
            optimization options. Among other things, this command line should
            enable Fortran source pre-processor, and define the variable
            __BUILD_ID__.

 F90L     - A command used to invoke Fortran linker. Normally, this is the same
            command line used to invoke the compiler:
              F90L = $(F90)

 LAPACK   - Link instructions for LAPACK and BLAS libraries. If quad-precision
            support is enabled, this line should also include the special
            quad-precision LAPACK library (see "Prerequisites" above). If OpenMP
            parallelization is enabled, the libraries MUST be re-entrant;  they
            will be called from inside the parallel region.

 LIBEXTRA - Any other libraries needed to build the executable. For example,
            FFTW3 libraries should appear here.

Please note that the default build configuration file emphasizes portability
over the performance and features. In particular, it disables quadruple-precision 
support, FFTW3 support, as well as all system-specific optimizations. 

You will almost certainly get a faster code by using system-specific optimization 
options.

Finally, say "make" and wait for the dust to settle. If nothing goes wrong,
make will produce two executable files, "./hydrogen_tunnel_v2.x" and
"./general_tunnel.x".

It is highly recommended to run at least some test cases, supplied in the
"examples/" subdirectory. Just stepping into the examples directory and running
make should do the trick.

3. Advanced installation

It is possible to build the codes to use quadruple-precision arithmetics on
systems which support it. Quadruple-precision option is especially useful for
multi-channel general_tunnex.x runs in weaker fields, where double precision
is often insufficient.

Doing so requires a quadruple-precision versions of LAPACK and BLAS. Building
these libraries from source is described in the file "doc/README-quad.txt" of
the SCID-TDSE repository at https://github.com/MBI-Theory/scid-tdse.git

If you are feeling adventurous, a few examples of advanced build configuration 
files can be found in ./configs/ subdirectory, including:

 gfortran_opt.mak        - GNU Fortran (version 7 or later)
 gfortran_quad_opt.mak   - GNU Fortran (version 7 or later), quadruple precision
 ifort_opt.mak           - Intel Fortran (2021 or later)
 ifort_quad_opt.mak      - Intel Fortran (2021 or later), quadruple precision

Please note that the advanced build configurations will likely require further 
debugging/testing.

4. Input

Two different codes are available: the purely-Coulombic code, essentially
following the prescriptions of Kolosov, and the generalized code, which
includes more general central potential and offers additional analysis options.

4.1. hydrogen_tunnel_v2.x

Evaluation of scattering solutions in the static field for a hydrogenic atom
using Kolosov's method. All input for this code is through the Fortran
namelist "ht_v2", in the format (the default value is listed in the example;
definition of the keywords is given below):

 &ht_v2
   verbose          = 1
   omp_num_threads  = 0
   comment          = ''
   task             = 'energy'
   mval             = 0
   efield           = 0.0
   znuc             = 1.0
   energy_single    = (-0.5,0.0)
   energy_guess     = (-0.5,0.0)
   zeta_guess       = (2.0,0.0)
   zeta_tol         = -1.0
   zeta_stab_cos    =  0.5
   energy_tol       = -1.0
   xi2_max          = 50.0
   xi_npts          =100
   xi_maxorder      = 60
   eta2_max         = 100.0
   eta_npts         = 200
   eta_maxorder     = 60
   asymp_order      = 200
   file_bound       = 'bound.table'
   file_continuum   = 'continuum.table'         
 /

All functionality of this code is also available through the more powerful
"general_tunnel.x" interface, which may however be harder to prepare the input
for.
   
verbose: Verbosity level of the output. Values above 1 are potentially useful
for debugging, but are not enlightening without a close look at the source
code.

omp_num_threads: Force the maximum number of threads. If this parameter is not
given, or set to zero, we will use as many threads we could make use of. The
code supports a very limited threading, with at most three parallel threads
being useful in some parts. Setting this variable should be equivalent to
setting the OMP_THREAD_LIMIT environment variable.

comment: A descriptive string, to be copied to the ouput.

task: Either 'energy' or 'outgoing'. For task='energy', a scattering solution
at the (complex) energy given by the energy_single parameter will be
calculated.  For task='outgoing', a purely-outgoing solution will be located
iteratively, starting from the complex energy given by energy_guess parameter.

mval: Magnetic quantum number of the desired solution. This value must be
non-negative.

efield: Static electric field, in atomic units. This value must be positive.

znuc: Nuclear charge, in units of the proton charge. Hydrogen atom is +1.

energy_single: Desired complex energy of the scattering solution
(task='energy'). For weak fields, resonances are expected close to
znuc**2/(2*n**2), where (b) is the principal quantum number.

energy_guess: Initial guess for an outgoing-solution search. Adding a (snall)
negative imaginary component is usually helpful.

zeta_guess: Initial guess for the zeta1 (bound component quasi-energy)
separation parameter. For weak fields, znuc*(2+k), k=0, 1, ... is a good
choice. Using larger zeta_guess values will produce solutions with more nodes
along the bound coordinate; these solutions may need larger xi2_max values to
be represented properly.

zeta_tol: Convergence tolerance for the zeta1 separation parameter.
Non-positive values request "machine precision"; the actual convergence
criterion in this case is 2 decimal digits less than the machine epsilon.
Tigher criteria will be capped at this value. Relevant for all task choices.

zeta_stab_cos: Critrion for the stability analysis of the bound-component
solution. Once the the direction cosine between the solution and its gradient
for the two perturbed solution is smaller than zeta_stab_cos, the stability is
lost.

energy_tol: Convergence tolerance for the outgoing-solution energy (similar to
zeta_tol). The parameter is only relevant for task='outgoing'. If optimization
appears to be stuck, it may be a good idea to use a less-tight convergence
criterion.

xi2_max: Maximum extent of the numerical grid for the bound part of the
solution. The grid is uniform in the squared coordinate. The first point of the
grid is at the origin; the last point is at sqrt(xi2_max). Higher resonances
and higher underlying numerical accuracy (quad-precision build) may require
larger radial grids; excessively large grids may however lead to numerical
problems due to the accuracy loss in the radial ODE solver or floating-point
overflows; it is advisable to visually examine the bound part of the
wavefunction for the final solutions. The solver will try to estimate the
largest extent of the grid where the solution is still well-behaved, but it
is not perfect.

xi_npts: Number of grid points for the bound part of the solution.  Please note
that due to the use of the high-order integration scheme, the Nyqvist limit
does not apply to the grid spacing!

xi_maxorder: Maximum derivative order in integration along the bound
coordinate. Using higher xi_maxorder allows coarser grids to be used, without
accuracy loss.

eta2_max, eta_npts, eta_maxorder: Similar to xi* variables above, but for the
unbound coordinate. The unbound-coordinate typically needs to extend further
away from the origin; at eta=sqrt(eta2_max) it is matched to the asymptotic
expansion (see asymp_order below). Note that *lower* electric fields may
require *larger* eta2_max values; make sure to examine the final continuum
solution.

asymp_order: Order of the asymptotic expation of the wavefunction along the
unbound coordinate (eta). Higher orders allow the matching point, determined by
eta2_max, to be closer to the origin. 

file_bound: File containing the final solution along the bound coordinate eta.

file_continuum: File containing the final solution along the continuum
coordinate xi.

4.2. general_tunnel.x

Due to a considerably more flexible potential support and additional solution
and analysis options, the input for this code is rather more complicated. The
main input is through the GT namelist, described below. Additional input may
follow the namelist, see below. The values listed in the namelist are the 
defaults; the description of the keywords follows below. 

 &gt
   verbose             = 1
   comment             = ""
   omp_num_threads     = 0
   task                = "energy"
   boundary            = "asymptotic"
   boundary_solution   = "incoming"
   boundary_index      = 1
   boundary_phase      = (1.0,0.0)
   outgoing_solver     = "min"
   mval                =   0
   efield              =   0.0
   znuc                =   1.0
   potential           = "hydrogenic"
   pot_param_real      =  6*0.0
   pot_param_int       =  3*0
   energy_single       = (-0.5,0.0)
   energy_guess        = (-0.5,0.0)
   energy_step         =   0.1
   zeta_tol            =  -1.0
   zeta_stab_cos       =   0.9
   energy_tol          =  -1.0
   wavefunction_tol    =  -1e4
   boundary_tol        =  -1.0
   xi_max              =   9.0
   xi_npts             =  90
   xi_maxorder         =  60
   xi_npts_guess       = 400
   eta_max             =  12.0
   eta_npts            = 360
   eta_maxorder        =  60
   eta_maxorder2       =   6
   eta_points2         =   9
   eta_order2          =   8
   nonad_points        =   7
   nonad_order         =   6
   asymp_order         =  60
   base_channel        =   1
   n_channels          =   1
   main_channel        =   1
   file_bound          = "('bound_chan', i0, '.table')"
   file_continuum      = "('continuum_chan', i0, '.table')"
   file_coupling       = "('coupling_', a, '.table')"
   file_total_mode     = "wavefunction"
   file_total          = "solution.table"
   fourier_centre      =  -1.0
   fourier_width       =  -1.0
   file_cartesian      = ""
   cartesian_dx        =   0.2
   cartesian_phi       =   0.0
   cartesian_npts      = -100, 100, -1, 1, -200,  100
   cartesian_ref       =  2*0.0, 0.5
   cart_interp_points  =   9
   cart_laplace_order  =   3
   file_bohm           = ""
   file_husimi         = ""
   husimi_ndim         =   1
   husimi_coord        =   3, 1
   husimi_width        =   2.5066282746310002
   husimi_detail       = .T.
  /

The keywords are:

verbose: Verbosity level of the output. Values above 1 are potentially useful
for debugging, but are not enlightening without a close look at the source
code.

comment: A descriptive string, to be copied to the ouput.

omp_num_threads: Limit the maximum number of threads to this value. If this
value is not given, or is zero, we will create as many threads as we appear
to be able to make use of - which, depending on the input parameters, may be
a lot more than the number of cores/CPUs available. For hydrogenic potentials,
we could make use of up to 3 threads; for general potentials, we might be able
to benefit from up to (3*eta_npts) threads. Setting this variable should be
equivalent to setting the OMP_THREAD_LIMIT environment variable.

task: Can be 'energy', 'outgoing', 'minflux real', or 'minflux imag'. For the
task='energy', determine the solution at a fixed energy, given by
energy_single. For the remaining three otions, the complex energy ('outgoing')
or the real/imaginary parts of the energu ('minflux real'/'minflux imag') are
optimized to minimize the outgoing flux. The initial guess for the energy is
given by energy_guess.

boundary: Can be 'origin all', 'origin single', 'asymptotic', or 'flux'. This
option determines the handling of the individual channel contributions in the
total solution. 

The simplest (but generally not physically relevant) case is boundary='origin
all'. In this case, (n_channels) complex values are expected to follow the
namelist on the standard input. These values are used as is, to integrate the
continuum part of the solution outwards, until the matching point with the
asymptotic solution. The boundary is taken to refer the mval-th derivative of
the wavefunction at the origin (all lower derivatives must vanish).

boundary='origin single' is generally similar to the 'origin all', except that
the inner boundary for the channel (main_channel) is set to one, and all others
to zero. This solution is physically relevant for purely Coulombic potentials,
where channels do not mix.

boundary='asymptotic' requires that the solution remains factorizable into the
bound and continuum part on the asymptotic region, for large values of the eta
coordinate. The channel selected by (main_channel) will have real, positive
amplitude at the origin. This choice is usually the desired one, and is the
default.

boundary='flux' requests that either the incoming or outgoing flux is
stationary with respect to the boundary amplitudes. The resulting solutions are
generally not factorizable at infinity. This choice is numerically delicate,
and is not guaranteed to converge.

The choices 'aymptotic' and 'flux' are further precised by boundary_solution
and boundary_index variables, see below.

boundary_solution: This parameter is relevant when (boundary) is set to
'asymptotic' or 'flux'. The possible choices are:

boundary_solution='magnitude', which is only applicable for
boundary='asymptotic', chooses the solution with the (boundary_index)-th lowest
magnitude. This is a good choice for task='outgoing'. For real energies, where
all boundary solutions have the same magnitude, this choice is unstable.

boundary_solution='phase', which is only applicable for boundary='asymptotic',
chooses the solution with the scattering phase closest to the value of
(boundary_phase). This choice is most useful for real-energy solutions.

boundary_solution='incoming' or 'outgoing' chooses the solution with
(boundary_index)-th smallest total incoming (or outgoing) amplitude. This is
generally the most robust choice. boundary_solution='incoming' is the default.

boundary_index: Solution index for boundary_solution='magnitude', 'incoming', or
'outgoing'. Using boundary_index/=1 for any task other than 'energy' will likely
not converge, or produce unexpected results.

boundary_phase: Desired scattering phase for boundary_solution='phase'

outgoing_solver: Solver used for all tasks other than 'energy'. Either 'root'
or 'min' is possible. outgoing_solver='root' is a good choice for
task='outgoing'.  For all other tasks outgoing_solver='min' is preferable. The
root finder is potentially more accurate, while the minimum finder tends to be
more robust, and is the default.

mval: Magnetic quantum number of the desired solution. This value must be
non-negative.

efield: Static electric field, in atomic units. This value must be positive.

znuc: Nuclear charge, in units of the proton charge. Hydrogen atom is +1.

potential: Central potential. The only long-range part of the potential must be
Coulombic, in the form -znuc/r. The following choices are available for the
short-range part (see subroutine potential_u_r in
general_tunnel_potential.f90):

potential='hydrogenic': No short-range part. Some special-case optimizations
will be applied where appropriate.

potential='[Tong05] He', '[Tong05] Ne', '[Tong05] Ar': Effective one-electron
potentials for the noble gases from Tong and Lin, J Phys B 38, 2593 (2005).

potential='[Tong05] Xe': Effective Xenon potential from Zhang, Lan, and Lu,
Phys Rev A 90, 043410 (2014).

potential='[SFG99] Li', '[SFG99] Na', '[SFG99] K', '[SFG99] Rb', '[SFG99] Cs':
Effective one-electron potentials from Schweizer, Fassbinder, and
Gonzalez-Ferez, At. Data Nucl. Data Tables 72, 33-55 (1999).

potential='yukawa': Generalized Yukawa potential, in the form:

  v(r) = Sum a0*r^n1*exp(-a1*r)

where quantities a0 and a1 are specified in the pot_param_real array (up to 6
values) and n1 in pot_param_int (up to 3 values).

pot_param_real: 
pot_param_int: Additional parameters for potential='yukawa'.

energy_single: Desired complex solution energy for task='energy'

energy_guess: Initial complex energy for all task values other than 'energy'

energy_step: For outgoing_solver='root', the largest step permitted during root
optimization. For outgoing_solver='min', characteristic initial step for the
linear search (see code in find_minimum.f90). Smaller values will improve
numerical stability, but likely make convergence slower.

zeta_tol: Convergence tolerance for the zeta1 separation parameter.
Non-positive values request "machine precision"; the actual convergence
criterion in this case is 2 decimal digits less than the machine epsilon.
Tigher criteria will be capped at this value. Relevant for all task choices.

zeta_stab_cos: Critrion for the stability analysis of the bound-component
solution. Once the the direction cosine between the solution and its gradient
for the two perturbed solution is smaller than zeta_stab_cos, the stability is
lost.

energy_tol: Convergence tolerance for the outgoing-solution energy (similar to
zeta_tol). The parameter is only relevant for task='outgoing'. If optimization
appears to be stuck, it may be a good idea to use a less-tight convergence
criterion.

wavefunction_tol: Desired tolerance on the tail of wavefunction. Positive
values are used as is. For negative values, tolerance is taken as:
abs(wavefunction_tol)*spacing(1._rk). This value is only used to estimate 
the radial size of the grid needed to contain the solution. The default
should be adequate for most cases.

boundary_tol: Desired convergence in iterative refinement of the boundary
conditions when boundary_solution='magnitude', 'incoming', or 'outgoing'.
Negative values imply machine accuracy.

xi_max: Maximum grid extent for the bound coordinate. Using a large number
of channels would typically require an increase in xi_max.

xi_npts: Number of grid points for the bound coordinate. This code only
supports uniform grids, so that xi_max and xi_npts fully specify the grid.
Please note that due to the use of the high-order integration scheme, the
Nyqvist limit does not apply to the grid spacing!

xi_maxorder: Maximum derivative order to use in integrating along the bound
coordinate.

xi_npts_guess: Number of grid points used for a separation-parameter guess. The
guess requires diagonalization of a tri-diagonal, complex-symmetrix matrix, so
that excessively large values will make calculation expensive. Values on the
order of 5*xi_npts should usually produce a sufficiently-accurate guess to use
in the (very accurate) shooting diagonalizer.

eta_max: Extent of the grid for the continuum coordinate. This value must be
large enough for the short-range part of the potential to vanish; only the
Coulombic potentials can be handled in the asymptotic region.

eta_npts: Number of grid points for the continuum coordinate. Because some of
non-adiabatic terms along the continuum coordinate are handled numerically,
denser grids are necessary for adequate convergence, compared to the
bound-coordinate grid. Please note that due to the use of the high-order
integration scheme, the Nyqvist limit does not apply to the grid spacing!

eta_maxorder: Maximum derivative order along the continuum coordinate, for the
terms known analytically.

eta_maxorder2: Maximum derivative order along the continuum coordinate, for the
terms which are evaluated by numerical differentiation.

eta_points2: Number of grid points used in numerical differentiation of the
terms known only numerically.

eta_order2: Order of the fit used for the numerical differentiation.

nonad_points: Number of grid points used for evaluating non-adiabatic matrix
elements. Odd values are recommended.

nonad_order: Fit order in evaluating non-adiabatic matrix elements.

asymp_order: Order used for calculating the asymptotic solutions. If warnings
on the Wronskian deviating from the exact solution start appearing in the
output, increasing asymp_order _may_ help.

base_channel: The first channel to include in the expansion. Channels are
sorted in the increasing order of the zeta separation parameter. Using
base_channel/=1 is probably not very meaningful for potentials other than
hydrogenic.

n_channels: Number of channels to include in the expansion. This number should
be much smaller than xi_npts_guess; using more than a handful of channels will
probably fail due to numerical issues. The cost of the calculation increases at
least quadratically with the number of channels.

main_channel: The channel we consider to be asymptotically dominant. This value
affects the phase and normalization of the solutions; see code in
general_tunnel_continuum.f90.

file_bound: Template for the wavefunctions in the bound channel. The value will
be used a Fortran format string, with the (integer) channel index given as the
argument. There is one wavefunction per grid coordinate in the continuum
channel. Blank argument suppresses the output.

file_continuum: Template for the wavefunctions in the continuum channel. A
Fortran format string, which will be used with an integer channel number.
Blank argument suppresses the output.

file_coupling: Template for the channel-coupling parameters. A Fortran format
string, which will be used with a character argument ("chan_alp", "chan_w1", or
"chan_w2"). Blank suppresses the output.

file_total_mode: The data to be included in the total solution in file_total.
Can be one of the following:

  'wavefunction'        = Total wavefunction, solving the Hamiltonian
  'asymptotic outgoing' = Asymptotic outgoing wavefunction, continued
                          inwards. This wavefunction is hydrogenic in
                          nature, and will diverge as we approach the 
                          origin.
  'asymptotic incoming' = ditto, incoming component
  'fourier outgoing'    = Total wavefunction, filtered to preserve only
                          the outgoing Fouriers components of the unbound
                          solution on the grid. This quantity is NOT an
                          eigenfunction of the Hamiltonian!
  'fourier incoming'    = ditto, incoming component
  'reverse outgoing'    = Outgoing part of the wavefunction, integrated
                          backwards from the matching point.
  'reverse incoming'    = ditto, incoming component
  'reverse total'       = ditto, total wavefunction. Please note that
                          reverse inegration leads to inherently less
                          accurate solutions compared to integration
                          from the origin.

The "fourier" options are only available if the code was built with the FFT
support; see above.

file_total: Filename for the final scattering solution on the squared-parabolic
grid. Blank suppresses the output.

fourier_centre: Position of the Blackman-Harris Fourier filtering window for
file_total_mode="fourier". The filter is applied to the unbound (eta)
coordinate.  This parameter works together with fourier_width.  The default is
to cover the entire unbound coordinate.

fourier_width: Width of the Fourier window for file_total_mode="fourier".
Negative values will use the entire range of the unbound coordinate. Please
make sure to read and understand the code in build_fft in
general_tunnel_dump.f90 before adjusting the defaults.

file_cartesian: Filename for the final wavefunction interpolated on a uniform
Cartesian grid. The grid parameters are defined by (cartesian_dx),
(cartesian_phi), and (cartesian_npts). Additional processing will occur if
file_bohm and/or file_husimi are also specified. Blank supresses the output.

cartesian_dx: Uniform Cartesian grid spacing. See (cartesian_npts) and
(cartesian_ref) below.

cartesian_phi: This angle will be added to the parabolic phi angle determined
from the Cartesian coordiantes, effectively rotating the atom by cartesian_phi
around the Z axis (the electric-field direction). This value has no effect when
mval=0, see code in "interpolate_wavefunction" in general_tunnel_dump.f90.

cartesian_npts: Six integers, giving the lowest and highest indices along the
X, Y, and Z axes, in this order. See also (cartesian_dx) and (cartesian_ref).
The grid positions are are (cartesian_ref(:)+(/ix,iy,iz/))*cartesian_dx(:).

cartesian_ref: Central position of the Cartesian intepolation grid, in units of
cartesian_dx(). It is not advisable to place an interpolation grid point at the
coordinate origin.

cart_interp_points: Number of the points along the unbound (eta) coordinate to
use for interpolating the bound part of the solution. If this sounds confusing,
check the code in interpolate_wavefunction in general_tunnel_dump.f90. Excessively
high orders will likely spoil the solution close to the origin!

cart_laplace_order: Order of the finite-difference expression for the Laplacian
on the interpolated Cartesian grid; can be 1, 2, or 3. For sensible results,
(cart_interp_points) must be at lease (2*cart_laplace_order+1).

file_bohm: If (file_bohm) is not blank, evaluate Bohmian velocity and Bohmian
quantum potential at each grid Cartesian grid point, and store the result in
(file_bohm). Note that for non-Coulombic potentials, the quantum potential
appears to be much more sensitive to the number of channels than the velocity.
The grid is determined by the (cart_interp_points) etc, even when
(file_cartesian) is blank. Additional processing, specified by
(file_total_mode) will be applied _before_ Bohmian trajectories are evaluated.
See subroutine build_bohm_potential_and_velocity in general_tunnel_dump.f90.

file_husimi: If (file_husimi) is not blank, evaluate Husimi distribution for
the wavefunction. In the 3D space, the Husimi distribution is 6-dimensional.
We have no hope to either store or sensibly visualize that much data.  Instead,
we'll apply the Husimi transformation along 1- or 2-dimensional slices only,
and leave the wavefunction untransformed along the remaining direction.

husimi_ndim: Dimensionality of the Husimi transform. Only the default choice
(1) is implemented.

husimi_coord: Coordinates, along which the Husimi transform is applied. Only
the default choice (3, a.k.a. Z) is implemented.

husimi_width: Width of the Husimi coordinate filter, in Bohr. See subroutine
husimi_1D in general_tunnel_dump.f90.

husimi_detail: .false. = print only the maximizing momentum and momentum
expectation.  .true.  = print the Husimi distribution as well. Can get -very-
large.

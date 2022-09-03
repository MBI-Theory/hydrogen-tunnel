!
!   hydrogen-tunnel: Static-field tunneling in a central potential
!   Copyright (C) 2018-2022 Serguei Patchkovskii, Serguei.Patchkovskii@mbi-berlin.de
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!  Data structures use by general_tunnel.f90 and its siblings
!
module general_tunnel_data
  use accuracy
  use constants
  !
  implicit none
  public
  public rcsid_general_tunnel_data
  !
  character(len=clen), save :: rcsid_general_tunnel_data = "$Id: $"
  !
  integer(ik)              :: verbose               = 1_ik              ! How verbose do we need to be?
  integer(ik)              :: omp_num_threads       = 0_ik              ! Non-zero value will cause number of OpenMP threads
                                                                        ! to be set explicitly. It looks like some libraries
                                                                        ! mess with the environment variables, and stop
                                                                        ! OMP_NUM_THREADS from working.
  character(len=clen)      :: comment               = ' '               ! Descriptive string, to be copied to the output
  character(len=clen)      :: task                  = 'energy'          ! What to do. Can be one of:
                                                                        ! 'energy'       - Find solution for a given (complex) energy
                                                                        ! 'outgoing'     - Solve for complex outgoing solution at a given field
                                                                        !                  The outgoing solution is sought for in the channel
                                                                        !                  specified by "main_channel" below. What happens to
                                                                        !                  other channels depends on the value of "boundary"
                                                                        ! 'minflux real' - Minimize incoming flux while keeping the imaginary
                                                                        !                  part of the energy fixed.
                                                                        ! 'minflux imag' - Same, but keep the real part of the energy fixed
  character(len=clen)      :: outgoing_solver       = 'min'             ! Either 'root' (use root finder, best for task='outgoing')
                                                                        ! or 'min' (use minimum finder, slower but potentially more robust)
                                                                        ! Despite the name, also applies to 'minflux real' and 'minflux imaginary'
  character(len=clen)      :: boundary              = 'asymptotic'      ! Boundary conditions to apply to the solution. Can be one of:
                                                                        ! 'origin single' - Only the "main_channel" has amplitude at the 
                                                                        !                   boundary. It is fixed at unity in the lowest
                                                                        !                   allowed derivate term.
                                                                        ! 'origin all'    - (Complex) amplitudes of the channels at the 
                                                                        !                   origin follow the namelist on input, in free 
                                                                        !                   format.
                                                                        ! 'asymptotic'    - The solution is asymptotically factorizable.
                                                                        !                   Channel selected by "main_channel" has positive
                                                                        !                   real amplitude at the origin.
                                                                        ! 'flux'          - Make either incoming or outgoing flux (controlled
                                                                        !                   by boundary_solution and boundary_index) stationary.
                                                                        !                   Please keep in mind that solutions for the 
                                                                        !                   boundary_index>1 will be less accurate; see the code.
                                                                        !                   Solving for outgoing/incoming solutions with 
                                                                        !                   boundary='flux' is generally less accurate than
                                                                        !                   boundary='asymptotic', and should be avoided unless
                                                                        !                   really essential.
                                                                        ! For pure Coulombic potential, 'origin single', 'asymptotic', and 'flux'
                                                                        ! are equivalent.
  character(len=clen)      :: boundary_solution     = 'incoming'        ! The choice of the boundary condition in a multi-channel calculation
                                                                        ! when boundary=='asymptotic' or 'flux'. The choices are:
                                                                        ! 'magnitude' - choose the (boundary_index)-th smallest solution.
                                                                        !               This is a good choice for task='outgoing' when energy
                                                                        !               is away from the real axis. On the real axis, all 
                                                                        !               solutions necessarily have the same magnitude (unity),
                                                                        !               and this choice becomes unstable. Not valid for 'flux'
                                                                        ! 'phase'     - choose the eigenvalue which has the phase most similar
                                                                        !               to the phase of (boundary_phase). This is a useful choice
                                                                        !               choice on the real axis. It may however yield a solution
                                                                        !               with a higher flux than the lowest-flux (most stable) 
                                                                        !               solution. Note that close to resonances, the phase of 
                                                                        !               of the eigenvalue may change rapidly! Not valid for 'flux'
                                                                        ! 'incoming'  - choose the (boundary_index)-th smallest total incoming 
                                                                        !               amplitude across all channels. This is generally the
                                                                        !               most robust option, especially close to the real axis.
                                                                        ! 'outgoing'  - ditto, outgoing
                                                                        ! Only 'incoming' or 'outgoing' may be used for boundary='flux'
  integer(ik)              :: boundary_index        = 1                 ! Parameter for boundary_solution=='magnitude'
  complex(rk)              :: boundary_phase        = (1._rk,0._rk)     ! Parameter for boundary_solution=='phase'
  integer(ik)              :: mval                  = 0                 ! Magnetic quantum number. Must be non-negative.
                                                                        ! Physically, mval<0 is allowed, but is the same as -mval>0;
  real(rk)                 :: efield                = 0._rk             ! Electric field magnitude. MUST be non-zero for the code to work!
  real(rk)                 :: znuc                  = 1._rk             ! Nuclear charge for the long-range part
  character(len=clen)      :: potential             = 'hydrogenic'      ! Central potential. Possible values are:
                                                                        ! 'hydrogenic'  - Purely Coulombic potential, defined by znuc
                                                                        !                 v(r) = -znuc/r
                                                                        ! 'yukawa'      - Yukawa-like potential, in the form:
                                                                        !                 v(r) = -znuc/r + a0*r**n1*exp(-a1*r) + ...
                                                                        !                 n1>=-1
                                                                        ! '[Tong05] He' - One-electron helium-like potential
                                                                        !                 Requires znuc=1 to produce meaningful results
                                                                        ! See potential_u_r() in general_tunnel_potential.f90 for the
                                                                        ! comprehensive list.
                                                                        ! The additional numerical parameters are in pot_param
  real(rk)                 :: pot_param_real(6)     = (/0,0,0,0,0,0/)   ! Meaning deprends in the value of (potential)
  integer(ik)              :: pot_param_int (3)     = (/0,0,0/)         ! Meaning deprends in the value of (potential)
  complex(rk)              :: energy_single         = (-0.5_rk,0.0_rk)  ! Energy to solve for (only for task='energy')
  complex(rk)              :: energy_guess          = (-0.5_rk,0.0_rk)  ! Energy to solve for (only for task='outgoing')
  real(rk)                 :: wavefunction_tol      = -1e4_rk           ! Desired tolerance on the tail of wavefunction.
                                                                        ! Positive values are used as is. For negative values,
                                                                        ! tolerance is taken as: abs(wavefunction_tol)*spacing(1._rk)
  real(rk)                 :: zeta_tol              = -1._rk            ! Desired convergence for the separation parameter
  real(rk)                 :: zeta_stab_cos         =  0.9_rk           ! Smallest cosine between the directions of the two
                                                                        ! perturbed solutions for zeta and their derivatives 
                                                                        ! before declaring stability loss. Setting zeta_stab_cos
                                                                        ! to -1 will effectively disable solution stability analysis.
  real(rk)                 :: energy_tol            = -1._rk            ! Desired convergence for the energy
  real(rk)                 :: boundary_tol          = -1._rk            ! Desired convergence for boundary='asymptotic'
                                                                        ! Negative values imply machine accuracy
  real(rk)                 :: energy_step           = 0.1_rk            ! For outgoing_solver='root': Limit on the size of a single step
                                                                        !                             in root optimization.
                                                                        ! For outgoing_solver='min: Characteristic size of the _first_
                                                                        !                           step during energy minimization.
  !                                                                     
  real(rk)                 :: xi_max                = 9._rk             ! Extent of the grid for the bound coordinate
                                                                        ! xi_max must be large enough to accommodate the highest
                                                                        ! channel. Using too large xi_max may cause numerical difficulties
                                                                        ! for low channels, so this is a compromise.
  integer(ik)              :: xi_npts               = 90_ik             ! Number of grid points for the bound coordinate
  integer(ik)              :: xi_maxorder           = 60_ik             ! Maximum derivative order to include in integration
  integer(ik)              :: xi_npts_guess         = 400_ik            ! Number of grid points to use for the initial guess
  !                                                                     
  real(rk)                 :: eta_max               = 12._rk            ! Extent of the grid for the continuum coordinate
                                                                        ! eta_max must be larger than the radius of the short-range
                                                                        ! part of the potential, or the results will be incorrect.
  integer(ik)              :: eta_npts              = 360_ik            ! Number of grid points for the continuum coordinate
  integer(ik)              :: eta_maxorder          = 60_ik             ! Maximum derivative order to include in integration, terms 
                                                                        ! known analytically
  integer(ik)              :: eta_maxorder2         = 6_ik              ! Maximum derivative order, terms known only numerically
  !
  integer(ik)              :: eta_points2           = 9_ik              ! Number of grid points used to evaluate numerical derivatives
                                                                        ! of the parameters known only numerically (chan_alp, chan_w1,
                                                                        ! and chan_w2).
  integer(ik)              :: eta_order2            = 8_ik              ! Order of the fit in derivative evaluation
  !
  integer(ik)              :: nonad_points          = 7_ik              ! Number of grid points used to evaluate wavefunction 
                                                                        ! derivatives for non-adiabatic matrix elements. Odd number
                                                                        ! of points is probably a good idea.
  integer(ik)              :: nonad_order           = 6_ik              ! Order of the fit in derivative evaluation
                                                                        ! Make sure to read through dt_savgol and understand how it works
                                                                        ! before adjusting nonad_points and nonad_order!
  !
  integer(ik)              :: asymp_order           = 60_ik             ! Order of the asymptotic expansion for the continuum solution
  !
  integer(ik)              :: base_channel          = 1_ik              ! First channel to include in the expansion. Channels are sorted by
                                                                        ! the real part of the separation parameter.
  !
  integer(ik)              :: n_channels            = 1_ik              ! Number of separable solutions ("channels") included in the expansion
                                                                        ! This number should not exceed half of (xi_npts_guess), or the initial
                                                                        ! guess for the channel energy will likely fail.
  integer(ik)              :: main_channel          = 1_ik              ! The dominant channel of the solution. The channel indices are within
                                                                        ! the (n_channels) window.  We may need a more sophisticated starting 
                                                                        ! guess in some situations ...
  character(len=clen)      :: file_bound            = '("bound_chan",i0,".table")'
  character(len=clen)      :: file_continuum        = '("continuum_chan",i0,".table")'
                                                                        ! Filename templates for the bound and countinuum channel
                                                                        ! wavefunctions output. Set to blank to suppress output.
  character(len=clen)      :: file_coupling         = '("coupling_",a,".table")'
                                                                        ! Filename template for the chan_alp, chan_w1, and chan_w2
                                                                        ! report. Set to blank to suppress output.
  character(len=clen)      :: file_total_mode       = 'wavefunction'    ! The information to include in the wavefunction (file_total)
                                                                        !   'wavefunction'        = Total wavefunction, solving the Hamiltonian
                                                                        !   'asymptotic outgoing' = Asymptotic outgoing wavefunction, continued
                                                                        !                           inwards. This wavefunction is hydrogenic in
                                                                        !                           nature, and will diverge as we approach the 
                                                                        !                           origin.
                                                                        !   'asymptotic incoming' = ditto, incoming component
                                                                        !   'fourier outgoing'    = Total wavefunction, filtered to preserve only
                                                                        !                           the outgoing Fouriers components of the unbound
                                                                        !                           solution on the grid. This quantity is NOT an
                                                                        !                           eigenfunction of the Hamiltonian!
                                                                        !   'fourier incoming'    = ditto, incoming component
                                                                        !   'reverse outgoing'    = Outgoing part of the wavefunction, integrated
                                                                        !                           backwards from the matching point.
                                                                        !   'reverse incoming'    = ditto, incoming component
                                                                        !   'reverse total'       = ditto, total wavefunction. Please note that
                                                                        !                           reverse inegration leads to inherently less
                                                                        !                           accurate solutions compared to integration
                                                                        !                           from the origin.
  character(len=clen)      :: file_total            = 'solution.table'  ! Filename for the final wavefunction on the squared-parabolic grid. 
                                                                        ! Blank will suppress the output.
  character(len=clen)      :: file_cartesian        = ' '               ! Filename for the final wavefunction interpolated on a uniform Cartesian
                                                                        ! grid, defined by cartesian_dx, cartesian_phi, and cartesian_npts() below
                                                                        ! Blank will suppress the output. Cartesian output can be further processed
                                                                        ! if file_bohm or file_husimi are not blank.
  real(rk)                 :: cartesian_dx          = 0.2_rk            ! Uniform Cartesian grid spacing
  real(rk)                 :: cartesian_phi         = 0.0_rk            ! Rotation angle along the Z axis, defining the local Cartesian coordinates
  integer(ik)              :: cartesian_npts(2,3)   = reshape( (/ -100_ik, 100_ik, -1_ik, 1_ik, -200_ik, 100_ik /), (/ 2, 3 /) )
                                                                        ! Lowest and highest point indices along the local Cartesian X, Y, and Z axes
  real(rk)                 :: cartesian_ref(3)      = (/ 0._rk, 0._rk, 0.5_rk /)
                                                                        ! The interpolation grid is centered on the point (cartesian_dx*cartesian_ref(:))
                                                                        ! Placing interpolation grid point at the origin potentially leads to difficulties
                                                                        ! with the singular potential and delta-function term in the Laplacian.
  integer(ik)              :: cart_interp_points    = 9_ik              ! Number of points to use for interpolating bound solutions in the 
                                                                        ! (file_cartesian) output. Using very high interpolation orders will
                                                                        ! likely spoil the solution close to the origin.
  integer(ik)              :: cart_laplace_order    = 3_ik              ! Order of the finite-differences operator for the Laplacian; either 1, 2, or 3
                                                                        ! Note that for a sensible result, we need (cart_interp_points>2*cart_laplace_order+1)
  real(rk)                 :: fourier_centre        = -1                ! Position of the Fourier window for file_total_mode="fourier ..."
  real(rk)                 :: fourier_width         = -1                ! Width of the Fourier window. Negative values will use the entire
                                                                        ! range of the unbound coordinate. Please make sure to read and understand 
                                                                        ! the code in build_fft in general_tunnel_dump.f90 before adjusting the
                                                                        ! defaults.
  character(len=clen)      :: file_bohm             = ' '               ! If (file_bohm) is not blank, evaluate Bohmian velocity and Bohmian
                                                                        ! quantum potential at each grid Cartesian grid point, and store the
                                                                        ! result in (file_bohm). Note that for non-Coulombic potentials,
                                                                        ! the quantum potential appears to be much more sensitive to the 
                                                                        ! number of channels than the velocity.
  character(len=clen)      :: file_husimi           = ' '               ! If (file_husimi) is not blank, evaluate Husimi distribution for the
                                                                        ! wavefunction. In the 3D space, the Husimi distribution is 6-dimensional.
                                                                        ! We have no hope to either store or sensibly visualize that much data.
                                                                        ! Instead, we'll apply the Husimi transformation along 1- or 2-dimensional
                                                                        ! slices only, and leave the wavefunction untransformed along the
                                                                        ! remaining direction.
  integer(ik)              :: husimi_ndim           = 1_ik              ! Dimensionality of the Husimi transform. 
  integer(ik)              :: husimi_coord(2)       = (/ 3_ik, 1_ik /)  ! Coordinates along which Husimi transform is applied.
  real(rk)                 :: husimi_width          = sqrt(2._rk*pi)    ! Width of the coordinate filter (FWHM)
  logical                  :: husimi_detail         = .true.            ! .false. = print only the maximizing momentum and momentum expectation
                                                                        ! .true.  = print the Husimi distribution as well. Can get -very- large
  !
  real(rk), allocatable, target  :: xi_tab(:)                           ! Radial grid for the bound coordinate (xi)
  real(rk), allocatable, target  :: eta_tab(:)                          ! Radial grid for the continnum coordinate (eta)
                                                                        ! The first point is always at zero. The grids must be
                                                                        ! uniform, so that xi_tab(2) and eta_tab(2) are grid
                                                                        ! spacings.
  real(rk), allocatable    :: tab_u_halfx2(:,:,:)                       ! Cached values of potential_u_halfx2() at all (xi,eta)
                                                                        ! grid points. This function is quite expensive to
                                                                        ! evaluate repeatedly!
  !
  !  We have multiple bound and countinuum solutions in works at any given time;
  !  We therefore have to differentiate between common data fields (prepared once during global
  !  initialization time) and data fields specific to individual solutions
  !
  !  We use a uniform description of the wavefunction for both bound and
  !  continuum coordinates. As the result, not all fields are valid under
  !  all circumstances.
  !
  type wavefunction
    complex(rk)              :: energy      ! Energy of the solution
    complex(rk)              :: zeta        ! Separation parameter for the solution
    real(rk)                 :: efield      ! Field value used; f>0 gives bound solutions; f<0 is continuum
    character(len=20)        :: norm        ! Current normalization of the solution; one of:
                                            ! 'invalid' - No solution present
                                            ! 'natural' - Lowest-order non-zero solution coefficient is
                                            !             unity at the origin
                                            ! 'bound'   - Normalized to unit total probability
                                            ! 'flux'    - Normalized to 1/(2 pi) outgoing flux
                                            ! 'square'  - Square of the function (not square modulus!) is unity
    complex(rk)              :: c_in        ! Amplitude of the incoming solution
    complex(rk)              :: c_out       ! Amplitude of the outgoing solution
                                            ! c_in and c_out are only valid for the continuum coordinate!
    integer(ik)              :: ipt_max     ! For bound solutions, the effective end of the grid. The
                                            ! solutions are guaranteed to fall below the significance
                                            ! threshold at this point. This value is determined from the
                                            ! outer turning point of the bound potential.
    integer(ik)              :: ipt_stop    ! Grid index at which radial integration stopped to be reliable. 
                                            ! Although the rest of the array is still calculated, it should
                                            ! not be treated as zero. Only relevant for the bound solutions. 
                                            ! ipt_stop is determined by the analysis of the perturbed bound solutions
                                            ! (see st_up/st_dn below).
    real(rk), pointer        :: r(:)        ! Grid positions. Points to either (xi_tab) or (eta_tab)
    integer(ik), allocatable :: iord(:)     ! Maximum power-series order at each grid point
    complex(rk), allocatable :: wf(:,:)     ! Scaled wavefunction and its derivatives at grid points
                                            ! First index:
                                            !   0 = wavefunction at the grid point
                                            !   i = i-th derivative of the wavefunction, times factorial of i (i!), 
                                            !       times i-th power of the grid spacing
                                            ! Second index: grid position, matches r() above
  end type wavefunction
  !
  !  Solution for a single channel
  !
  type channel
    integer(ik)                     :: ichan      ! Channel index within the bound eigenspectrum
    type(wavefunction), allocatable :: bound(:)   ! Bound solutions for this channel, one per grid point 
                                                  ! of the continuum solution. We have to use "pointer" in
                                                  ! order to be able to use pointer to this structure element later on
    type(wavefunction), allocatable :: st_up(:)   ! Bound solutions, zeta separation parameter perturbed "up"
    type(wavefunction), allocatable :: st_dn(:)   ! Bound solutions, zeta separation parameter perturbed "down"
  end type channel
  !
  !  Channel table. The wavefunction is a sum over all channels
  !
  type(channel), allocatable, target      :: channels(:)       ! There are n_channels of these bad boys.
  type(wavefunction), allocatable, target :: continuum(:)      ! Continuum solutions for each channel
  type(wavefunction), allocatable, target :: pseudo(:)         ! Pseudo-wavefunction for each channel, obtained by
                                                               ! integrating the incoming/outbound part of the
                                                               ! solution backwards from the matching point. Is used
                                                               ! by file_total_mode=='reverse ...' and is not calculated
                                                               ! otherwise. For low fields, accuracy of these solutions
                                                               ! may be low. Only the function and its derivative are
                                                               ! evaluated at the origin, regardless of the (mval) value
                                                               ! - inverse integration does not honour boundary conditions
                                                               ! at the origin.
  !
  !  Per-channel amplitudes at the origin
  !
  complex(rk), allocatable    :: inner_boundary(:) 
  !
  !  Additional parameters we need for solving the continuum equations.
  !  These are derived from the bound solutions in channels().
  !
  complex(rk), allocatable    :: chan_alp(:,:,  :)   ! Adiabatic separation parameters and their derivatives
                                                     ! chan_alp(0,:,:) is taken from channels(:)%bound(:)%zeta
  complex(rk), allocatable    :: chan_w1 (:,:,:,:)   ! Gradient non-adiabatic terms and their derivatives
  complex(rk), allocatable    :: chan_w2 (:,:,:,:)   ! 2-nd order non-adiabatic terms and derivatives
                                                     ! In alp, w1, and w2:
                                                     ! The first index is the order of the derivative; 0 to eta_maxorder2
                                                     !   0 = the value at the eta grid point
                                                     !   i = i-th derivative, times factorial of i (i!)
                                                     !
                                                     ! The second (alp) and second+third (w1, w2) indices are the channels,
                                                     ! from 1 to n_channel. The wavefunction derivatives in w1 and w2 are 
                                                     ! taken in the channel specified by the third index.
                                                     !
                                                     ! The last index is the eta grid point index.
  !
  !  Simulation parameters; we are collecting variables from many modules.
  !
  namelist /gt/ &
             ! Parameters defined locally
             verbose, comment, &
             omp_num_threads, &
             task, boundary, boundary_solution, boundary_index, boundary_phase, &
             outgoing_solver, &
             mval, efield, znuc, &
             potential, pot_param_real, pot_param_int, &
             energy_single, energy_guess, energy_step, &
             zeta_tol, zeta_stab_cos, energy_tol, wavefunction_tol, boundary_tol, &
             xi_max, xi_npts, xi_maxorder, xi_npts_guess, &
             eta_max, eta_npts, eta_maxorder, &
             eta_maxorder2, eta_points2, eta_order2, &
             nonad_points, nonad_order, &
             asymp_order, &
             base_channel, n_channels, main_channel, &
             file_bound, file_continuum, file_coupling, &
             file_total_mode, file_total, &
             fourier_centre, fourier_width, &
             file_cartesian, cartesian_dx, cartesian_phi, cartesian_npts, &
             cartesian_ref, cart_interp_points, cart_laplace_order, &
             file_bohm, file_husimi, husimi_ndim, husimi_coord, husimi_width, &
             husimi_detail
  !
  contains
  !
  subroutine init_wavefunction_1D(state,rtab,maxorder)
    type(wavefunction), intent(inout) :: state      ! Wavefunction structure to initialize
    real(rk), target                  :: rtab(:)    ! Radial grid
    integer(ik), intent(in)           :: maxorder   ! Order of Taylor expansion
    !
    integer(ik) :: alloc, npts
    !
    npts = size(rtab)
    !
    allocate (state%wf(0:maxorder,npts),state%iord(npts),stat=alloc)
    if (alloc/=0) then
      write (out,"('Allocation failed in init_wavefunction_1D. code = ',i0)") alloc
      stop 'general_tunnel%init_wavefunction_1D - allocation failed' 
    end if
    !
    state%r      => rtab
    !
    !  Initialize scalar parameters
    !
    state%ipt_max  = npts
    state%ipt_stop = npts
    state%energy   = huge(1._rk)
    state%zeta     = huge(1._rk)
    state%efield   = huge(1._rk)
    state%norm     = 'invalid'
  end subroutine init_wavefunction_1D
  !
  subroutine destroy_wavefunction_1D(state)
    type(wavefunction), intent(inout) :: state      ! Wavefunction structure to release
    !
    deallocate (state%wf,state%iord)
    nullify(state%r)
  end subroutine destroy_wavefunction_1D
  !
  subroutine init_grids
    integer(ik) :: alloc ! Allocation status
    
    allocate (xi_tab(xi_npts),eta_tab(eta_npts),stat=alloc)
    if (alloc/=0) stop 'general_tunnel%init_grids - error allocating grids'
    !
    call fill_grid(xi_tab, xi_max, xi_npts)
    call fill_grid(eta_tab,eta_max,eta_npts)
    !
    contains
    subroutine fill_grid(r,rmax,npts)
      real(rk), intent(out)   :: r(:)
      real(rk), intent(in)    :: rmax
      integer(ik), intent(in) :: npts
      real(rk)                :: h
      integer(ik)             :: ipt
      !
      if (size(r)/=npts) stop 'general_tunnel%init_grids - assumption failure'
      !
      h = rmax/(npts-1)
      fill_radial: do ipt=1,npts
        r(ipt) = h * (ipt-1)
      end do fill_radial
    end subroutine fill_grid
  end subroutine init_grids
  !
  subroutine init_global_arrays
    integer(ik) :: ic    ! Channel
    integer(ik) :: ipt   ! Point within a channel
    integer(ik) :: alloc ! Allocation status
    !
    call init_grids
    allocate (channels(n_channels), &
              continuum(n_channels), &
              pseudo(n_channels), &
              inner_boundary(n_channels), &
              chan_alp(0:eta_maxorder2,n_channels,           eta_npts), &
              chan_w1 (0:eta_maxorder2,n_channels,n_channels,eta_npts), &
              chan_w2 (0:eta_maxorder2,n_channels,n_channels,eta_npts), &
              tab_u_halfx2(0:xi_maxorder,xi_npts,eta_npts), &
              stat=alloc)
    if (alloc/=0) stop 'general_tunnel%init_global_arrays - error allocating channels (1)'
    init_basis: do ic=1,n_channels
      channels(ic)%ichan = base_channel + ic - 1
      allocate (channels(ic)%bound(eta_npts), &
                channels(ic)%st_up(eta_npts), &
                channels(ic)%st_dn(eta_npts),stat=alloc)
      if (alloc/=0) stop 'general_tunnel%init_global_arrays - error allocating channels (2)'
      !
      call init_wavefunction_1D(continuum(ic),eta_tab,eta_maxorder) 
      !
      call init_wavefunction_1D(pseudo(ic),eta_tab,eta_maxorder) 
      !
      init_bound_channels: do ipt=1,eta_npts
        call init_wavefunction_1D(channels(ic)%bound(ipt),xi_tab,xi_maxorder)
        call init_wavefunction_1D(channels(ic)%st_up(ipt),xi_tab,xi_maxorder)
        call init_wavefunction_1D(channels(ic)%st_dn(ipt),xi_tab,xi_maxorder)
      end do init_bound_channels
    end do init_basis
    !
  end subroutine init_global_arrays
end module general_tunnel_data

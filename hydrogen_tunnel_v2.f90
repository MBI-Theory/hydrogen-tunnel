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
!  Ionization of a hydrogen atom in an intense static field, squared-parabolic coordinates
!  The code follows Kolosov 1987 paper [J Phys B 20, 2359 (1987)].
!  This is an version solves directly for the position of the complex pole. 
!  IEEE double precision should yield accurate results for all electric fields where
!  ionization rate is measurable.
!
module hydrogen_tunnel_v2
  use accuracy
  use constants
  use find_root
  use poly_tools
  use timer
  implicit none
  private
  public start
  public rcsid_hydrogen_tunnel_v2
  !
  character(len=clen), save :: rcsid_hydrogen_tunnel_v2 = "$Id: $"
  !
  integer(ik)              :: verbose               = 1_ik              ! How verbose do we need to be?
  integer(ik)              :: omp_num_threads       = 0_ik              ! Non-zero value will cause number of OpenMP threads
                                                                        ! to be set explicitly. It looks like some libraries
                                                                        ! mess with the environment variables, and stop
                                                                        ! OMP_NUM_THREADS from working.
  character(len=clen)      :: comment               = ' '               ! Descriptive string, to be copied to the output
  character(len=clen)      :: task                  = 'energy'          ! What to do. Can be one of:
                                                                        ! 'energy'   - Find solution for a given (complex) energy
                                                                        ! 'outgoing' - Solve for complex outgoing solution at a given fielx
  character(len=clen)      :: file_bound            = 'bound.table'     ! Wavefunction for the bound coordinate
  character(len=clen)      :: file_continuum        = 'continuum.table' ! Wavefunction for the continuum coordinate
  integer(ik)              :: mval                  = 0                 ! Magnetic quantum number. Must be non-negative.
                                                                        ! Physically, mval<0 is allowed, but is the same as -mval>0;
  real(rk)                 :: efield                = 0._rk             ! Electric field magnitude
  real(rk)                 :: znuc                  = 1._rk             ! Nuclear charge
  complex(rk)              :: energy_single         = (-0.5_rk,0.0_rk)  ! Energy to solve for (only for task='energy')
  complex(rk)              :: energy_guess          = (-0.5_rk,0.0_rk)  ! Initial energy for the solution search (only for task='outgoing')
  complex(rk)              :: zeta_guess            = (2.0_rk,0.0_rk)   ! Initial guess for the bound zeta1 separation parameter
                                                                        ! For hydrogenic states, good initial values are close to
                                                                        ! znuc*(2+n), for n=0...
  real(rk)                 :: zeta_tol              = -1._rk            ! Desired convergence for the separation parameter
                                                                        ! Zero and negative values correspond to machine accuracy
  real(rk)                 :: zeta_stab_cos         =  0.5_rk           ! Smallest cosine between the directions of the two
                                                                        ! perturbed solutions for zeta and their derivatives 
                                                                        ! before declaring stability loss. Setting zeta_stab_cos
                                                                        ! to -1 will effectively disable solution stability analysis.
  real(rk)                 :: energy_tol            = -1._rk            ! Desired convergence for the energy
                                                                        ! Negative values imply machine accuracy
  !                                                                     
  real(rk)                 :: xi2_max               = 50._rk            ! Extent of the grid for the bound coordinate
  integer(ik)              :: xi_npts               =100_ik             ! Number of grid points for the bound coordinate
  integer(ik)              :: xi_maxorder           = 60_ik             ! Maximum derivative order to include in integration
  !                                                                     
  real(rk)                 :: eta2_max              = 100._rk           ! Extent of the grid for the continuum coordinate
  integer(ik)              :: eta_npts              = 200_ik            ! Number of grid points for the continuum coordinate
  integer(ik)              :: eta_maxorder          = 60_ik             ! Maximum derivative order to include in integration
  !
  integer(ik)              :: asymp_order           = 200_ik            ! Order of the asymptotic expansion for the continuum solution
  !
  !  We may have multiple bound and countinuum solutions in works at any given time;
  !  We therefor have to differentiate between common data fields (prepared once during global
  !  initialization time) and data fields specific to individual solutions
  !
  !  We use a uniform description of the wavefunction for both bound and
  !  continuum states. As the results, not all fields are valid under
  !  all circumstances.
  !
  type wavefunction
    complex(rk)          :: energy      ! Energy of the solution
    complex(rk)          :: zeta        ! Separation parameter for the solution
    real(rk)             :: efield      ! Field value used; f>0 gives bound solutions; f<0 is continuum
    character(len=20)    :: norm        ! Current normalization of the solution; one of:
                                        ! 'invalid' - No solution present
                                        ! 'natural' - Lowest-order non-zero solution coefficient is
                                        !             unity at the origin
                                        ! 'bound'   - Normalized to unit total probability
                                        ! 'flux'    - Normalized to 1/(2 pi) outgoing flux
    complex(rk)          :: c_in        ! Amplitude of the incoming solution
    complex(rk)          :: c_out       ! Amplitude of the outgoing solution
                                        ! c_in and c_out are only valid for the continuum coordinate!
    integer(ik)          :: ipt_stop    ! Grid index at which radial integration stopped to be reliable. 
                                        ! Although the rest of the array is still calculated, it should
                                        ! not be treated as zero. Only relevant for the bound solutions. 
    real(rk), pointer    :: r(:)        ! Grid positions. Grids start at zero.
    complex(rk), pointer :: wf(:,:)     ! Scaled wavefunction and its derivatives at grid points
                                        ! First index:
                                        !   0 = wavefunction at the grid point
                                        !   i = i-th derivative of the wavefunction, times i!
                                        ! Second index: grid position, matches r() above
    integer(ik), pointer :: wf_scale(:) ! Wavefunction can have a very large dynamic range;
                                        ! rather than using multiple-precision arithmetics,
                                        ! it is easier to carry additional exponent parameter.
                                        ! The actual wavefunction is (assumung infinite dynamic range):
                                        !   wf(0,:) * radix(1._rk)**wf_scale(:)
  end type wavefunction
  !
  !  Simulation parameters; we are collecting variables from many modules.
  !
  namelist /ht_v2/ &
             ! Parameters defined locally
             verbose, omp_num_threads, comment, &
             task, &
             mval, efield, znuc, &
             energy_single, energy_guess, &
             zeta_guess, &
             zeta_tol, zeta_stab_cos, energy_tol, &
             xi2_max, xi_npts, xi_maxorder, &
             eta2_max, eta_npts, eta_maxorder, &
             asymp_order, &
             file_bound, file_continuum
  !
  contains
  !
  subroutine init_wavefunction(state,r2max,npts,maxorder)
    type(wavefunction), intent(inout) :: state
    real(rk), intent(in)              :: r2max      ! Since we use square-parabolic coordinates,
                                                    ! it is easier to specify the *square* of the
                                                    ! maximum grid extent
    integer(ik), intent(in)           :: npts
    integer(ik), intent(in)           :: maxorder
    !
    integer(ik) :: alloc, ipt
    real(rk)    :: h
    !
    if (npts<2 .or. maxorder<1 .or. r2max<=0._rk) stop 'hydrogen_tunnel_v2%init_wavefunction - nonsense call'
    if (mval<0) stop 'hydrogen_tunnel_v2%init_wavefunction - negative m values are not supported'
    !
    allocate (state%r(npts),state%wf(0:maxorder,npts),state%wf_scale(npts),stat=alloc)
    if (alloc/=0) stop 'hydrogen_tunnel_v2%init_wavefunction - allocation failed' 
    !
    !  Initialize scalar parameters
    !
    state%energy   = huge(1._rk)
    state%zeta     = huge(1._rk)
    state%efield   = huge(1._rk)
    state%norm     = 'invalid'
    state%ipt_stop = npts
    !
    !  Fill the radial grid. For the moment, we choose grid which becomes uniform
    !  after it is squared.
    !
    h = r2max/(npts-1)
    fill_radial: do ipt=1,npts
      state%r(ipt) = sqrt(h * (ipt-1))
    end do fill_radial
  end subroutine init_wavefunction
  !
  subroutine destroy_wavefunction(state)
    type(wavefunction), intent(inout) :: state
    !
    deallocate (state%r,state%wf,state%wf_scale)
    nullify (state%r,state%wf,state%wf_scale)
  end subroutine destroy_wavefunction
  !
  subroutine integrate_radial(state,energy,zeta,efield)
    type (wavefunction), intent(inout) :: state
    complex(rk), intent(in)            :: energy ! Energy of the solution
    complex(rk), intent(in)            :: zeta   ! Separation parameter
    real(rk), intent(in)               :: efield ! Electric field
    !
    integer(ik) :: npts, maxorder
    integer(ik) :: ipt, iord
    complex(rk) :: temp
    real(rk)    :: xi, step
    complex(rk) :: c2, c3, c4, c5, c6, c7  ! Fixed part of the recursion coefficients
    !
    !  Store solution parameters
    !
    state%energy = energy
    state%zeta   = zeta
    state%efield = efield
    state%norm   = 'natural'
    !
    npts         = size(state%wf,dim=2)
    maxorder     = ubound(state%wf,dim=1)
    !
    !  The origin of the grid requires special handling due to the
    !  presence of singularities in the Hamiltonian.
    !
    state%wf_scale(1) = 0
    state%wf (   :,1) = 0
    state%wf (mval,1) = 1
    fill_origin: do iord=mval+1,maxorder
      temp = 0
      if (iord>=2) temp = temp + state%wf(iord-2,1) * state%zeta
      if (iord>=4) temp = temp + state%wf(iord-4,1) * state%energy * 2
      if (iord>=6) temp = temp + state%wf(iord-6,1) * state%efield
      state%wf(iord,1) = -temp / (iord**2 - mval**2)
    end do fill_origin
    call poly_adjust_scale(state%wf(:,1),state%wf_scale(1))
    !
    !  We can use general code for the rest of the radial grid
    !
    propagate: do ipt=2,npts
      !
      ! Use power-series expansion at the previous point to get wavefunction
      ! and its gradient
      !
      step = state%r(ipt) - state%r(ipt-1)
      call poly_power_series(state%wf(:,ipt-1),step,state%wf(0:1,ipt))
      !
      ! Calculate power-series coefficients at the new point
      !
      xi = state%r(ipt)
      c2 = -(mval**2) + xi**2*state%zeta + 2*xi**4*state%energy + xi**6*state%efield
      c3 = 2*xi*state%zeta + 8*xi**3*state%energy + 6*xi**5*state%efield
      c4 = state%zeta + 12*xi**2*state%energy + 15*xi**4*state%efield
      c5 = 8*xi*state%energy + 20*xi**3*state%efield
      c6 = 2*state%energy + 15*xi**2*state%efield
      c7 = 6*xi*state%efield
      !
      !  From the second term on, we can use the general expression
      !
      propagate_power: do iord=2,maxorder
                     temp =        state%wf(iord-1,ipt) * xi*(iord-1)*(2*iord-3)
        if (iord>=2) temp = temp + state%wf(iord-2,ipt) * ((iord-2)**2 + c2)
        if (iord>=3) temp = temp + state%wf(iord-3,ipt) * c3
        if (iord>=4) temp = temp + state%wf(iord-4,ipt) * c4
        if (iord>=5) temp = temp + state%wf(iord-5,ipt) * c5
        if (iord>=6) temp = temp + state%wf(iord-6,ipt) * c6
        if (iord>=7) temp = temp + state%wf(iord-7,ipt) * c7
        if (iord>=8) temp = temp + state%wf(iord-8,ipt) * state%efield
        !
        state%wf(iord,ipt) = - temp / (iord*(iord-1)*xi**2)
      end do propagate_power
      !
      ! Carry over scale, and adjust the scale of the result
      !
      state%wf_scale(ipt) = state%wf_scale(ipt-1)
      call poly_adjust_scale(state%wf(:,ipt),state%wf_scale(ipt))
    end do propagate
  end subroutine integrate_radial
  ! 
  !  Find the value of the separation parameter which yields a bound solution at
  !  a given energy. solve_bound() is a little tricky: the outer part of the
  !  solution is decaying exponentially, and we potentially need to stop the
  !  outward integration before reaching the outer boundary - otherwise, it
  !  starts growing exponentially again. However, until we are very close to
  !  the solution for the zeta, the wavefunction keeps growing exponentially,
  !  making our life quite miserable.
  !
  subroutine solve_bound(verbose,state,st_up,st_dn,energy)
    integer(ik), intent(in)           :: verbose
    type(wavefunction), intent(inout) :: state   ! Wavefunction
    type(wavefunction), intent(inout) :: st_up   ! Perturbed wavefunctions for stability analysis
    type(wavefunction), intent(inout) :: st_dn   ! 
    complex(rk), intent(in)           :: energy
    !
    type (fr_state) :: root
    real(rk)        :: tail_re, tail_im
    integer(ik)     :: npts
    real(rk)        :: vars(2)
    !
    call TimerStart('Bound direction')
    vars(1) =  real(zeta_guess,kind=rk)
    vars(2) = aimag(zeta_guess)
    call fr_initialize(verbose,root,neqs=2_ik,vars=vars,epsstep=zeta_tol,diffstep=zeta_tol)
    !
    root_iterations: do
      call fr_step(root)
      select case (root%action)
        case default
          write (out,"('hydrogen_tunnel_v2%solve_bound: unexpected request from equation solver: ',i0)") root%action
          stop 'hydrogen_tunnel_v2%solve_bound - state machine corrupted'
        case (fr_act_evaluate)
          call evaluate_function
        case (fr_act_done)
          if (verbose>=1) write (out,"('Optimization of separation parameter zeta1 converged.')")
          exit root_iterations
        case (fr_act_failed)
          write (out,"('WARNING: Optimization of separation parameter zeta1 failed to converge.')")
          write (out,"('WARNING: Continuing with possibly unconverged zeta1 value.')")
          flush (out)
          exit root_iterations
      end select
    end do root_iterations
    call evaluate_function
    !
    if (verbose>=0) then
      npts    = state%ipt_stop
      tail_re = scale( real(state%wf(0,npts),kind=rk),state%wf_scale(npts))
      tail_im = scale(aimag(state%wf(0,npts)),        state%wf_scale(npts))
      write (out,"(/'Final bound-state solution')")
      write (out,"('     Stop at point = ',i0,' = ',g44.34e3)") npts, state%r(npts)
      write (out,"('             Zeta1 = ',2g44.34e3)") state%zeta
      write (out,"('            Energy = ',2g44.34e3)") state%energy
      write (out,"(' Wavefunction tail = ',2g44.34e3)") tail_re, tail_im
      flush (out)
    end if
    !
    call fr_finalize(root)
    call TimerStop('Bound direction')
    !
    contains
    !
    !  Our optimization goal is the wavefunction at the last valid grid point.
    !  We locate the last valid point by running the radial integral three
    !  times: once with the desired zeta, then with zeta a bit higher and
    !  a bit lower. We stop at the point where the argument of the two perturbed
    !  solutions no longer point in the same direction.
    !
    subroutine evaluate_function
      complex(rk) :: zeta, zeta_up, zeta_dn
      integer(ik) :: npts, ipt
      real(rk)    :: f1, f2, cos_f, cos_g
      !
      zeta    = cmplx(root%vars(1),root%vars(2),kind=rk)
      zeta_up = zeta * (1._rk+spacing(5._rk))
      zeta_dn = zeta * (1._rk-spacing(5._rk))
      !$omp parallel sections num_threads(3)
      !$omp section  ! The solution we want
      call integrate_radial(state,energy,zeta,   -efield)
      !$omp section  ! The perturbed-zeta solutions
      call integrate_radial(st_up,energy,zeta_up,-efield)
      !$omp section  
      call integrate_radial(st_dn,energy,zeta_dn,-efield)
      !$omp end parallel sections
      npts = ubound(state%wf,dim=2)
      state%ipt_stop = npts  ! Start by assuming that the entire grid is OK
                             ! We'll skip the first point, which may be zero
      find_the_end: do ipt=2,npts 
        cos_f = poly_direction_angle(st_up%wf(0,ipt),st_dn%wf(0,ipt))
        cos_g = poly_direction_angle(st_up%wf(1,ipt),st_dn%wf(1,ipt))
        if (cos_f<zeta_stab_cos .and. cos_g<zeta_stab_cos) then
          state%ipt_stop = ipt
          exit find_the_end
        end if
      end do find_the_end
      !
      npts = state%ipt_stop
      f1   = scale( real(state%wf(0,npts),kind=rk),state%wf_scale(npts))
      f2   = scale(aimag(state%wf(0,npts)),        state%wf_scale(npts))
      root%eqs(:) = (/ f1, f2 /)
      if (verbose>=2) then
        write (out,"('Z1 = ',2g32.20e3,' WF @ ',i0,' = ',2g32.20e3)") zeta, state%ipt_stop, root%eqs
      end if
    end subroutine evaluate_function
  end subroutine solve_bound
  !
  !  Constructs an asymptotic solution as a power-series expansion:
  !
  !   sqrt(x)*psi = Sum cn(k) x**(-k) exp( pm I [(sqrt(F)/3)*x**3 + (E/sqrt(F))*x ] )
  !
  !  where pm = +1 for the outgoing solution and -1 for the incoming wave.
  !
  !  The coefficients cn(k) are determined using recursively.
  !
  subroutine asymptotic_solution(efield,energy,zeta2,pm,cn)
    real(rk), intent(in)     :: efield             ! Electric field strength
    complex(rk), intent(in)  :: energy             ! Energy of the solution
    complex(rk), intent(in)  :: zeta2              ! Separation parameter 
    real(rk), intent(in)     :: pm                 ! +1 for outgoing solution; -1 for an incoming
    complex(rk), intent(out) :: cn(1:asymp_order)  ! Power series coefficients for the solution
                                                   ! cn(1) is arbitrarily chosen as 1.
    !
    integer(ik) :: l
    complex(rk) :: temp
    !
    if (efield<=0) stop 'hydrogen_tunnel_v2%asymptotic_solution - efield must be positive'
    !
    cn(1) = 1._rk
    recurrence: do l=1,asymp_order-1
               temp =        pm*(0,1)*((energy**2)/efield - zeta2)      * cn(l)
      if (l>1) temp = temp - 2*(energy/sqrt(efield))*(l-1)              * cn(l-1)
      if (l>2) temp = temp - pm*(0,1)*((l-2)*(l-1) + (0.25_rk-mval**2)) * cn(l-2)
      cn(l+1) = temp/(2*l*sqrt(efield))
    end do recurrence
    !
    if (verbose>=3) then
      write (out,"(/'Asymptotic solution for pm=',f3.0)") pm
      write (out,"(2(1x,g34.22e3,1x,g34.22e3,1x))") cn
      write (out,"()")
    end if
  end subroutine asymptotic_solution
  !
  !  Evaluates asymptotic series for the gradient of the incoming/outgoing solution
  !
  subroutine asymptotic_gradient(efield,energy,pm,cn,gr)
    real(rk), intent(in)     :: efield                ! Electric field strength
    complex(rk), intent(in)  :: energy                ! Energy of the solution
    real(rk), intent(in)     :: pm                    ! +1 for outgoing solution; -1 for an incoming
    complex(rk), intent(in)  :: cn( 1:asymp_order)    ! Solution coefficients
    complex(rk), intent(out) :: gr(-1:asymp_order+1)  ! Gradient coefficients
    !
    integer(ik) :: l
    complex(rk) :: temp
    !
    if (efield<=0) stop 'hydrogen_tunnel_v2%asymptotic_gradient - efield must be positive'
    !
    gradient: do l=-1,asymp_order+1
      temp = 0
      if (l>= 1 .and. l<=asymp_order  ) temp = temp + pm*(0,1)*(energy/sqrt(efield)) * cn(l)
      if (l>= 2 .and. l<=asymp_order+1) temp = temp - (l-1)                          * cn(l-1)
      if (l>=-1 .and. l<=asymp_order-2) temp = temp + pm*(0,1)*sqrt(efield)          * cn(l+2)
      gr(l) = temp
    end do gradient
    !
    if (verbose>=3) then
      write (out,"(/'Asymptotic gradient for pm=',f3.0)") pm
      write (out,"(2(1x,g34.22e3,1x,g34.22e3,1x))") gr
      write (out,"()")
    end if
  end subroutine asymptotic_gradient
  !
  function evaluate_asymptote(efield,energy,pm,lead,cn,x) result (v)
    real(rk), intent(in)     :: efield ! Electric field strength
    complex(rk), intent(in)  :: energy ! Energy of the solution
    real(rk), intent(in)     :: pm     ! +1 for outgoing solution; -1 for an incoming
    integer(ik), intent(in)  :: lead   ! Leading power of x
    complex(rk), intent(in)  :: cn(:)  ! Solution coefficients
    real(rk), intent(in)     :: x      ! Coordinate where we need the solution
    complex(rk)              :: v
    !
    real(rk)    :: xm1
    integer(ik) :: nterms, i
    complex(rk) :: phase
    !
    nterms = size(cn)
    xm1    = 1._rk/x
    v      = cn(nterms)
    sum_terms: do i=nterms-1,1,-1
      v = cn(i) + xm1*v
    end do sum_terms
    v     = v * xm1
    if (lead/=0) v = v * x**lead
    phase = (sqrt(efield)/3._rk)*x**3 + (energy/sqrt(efield))*x
    v     = v * exp(pm*(0,1)*phase)
  end function evaluate_asymptote
  !
  subroutine solve_continuum(verbose,state,energy,zeta)
    integer(ik), intent(in)           :: verbose
    type(wavefunction), intent(inout) :: state
    complex(rk), intent(in)           :: energy
    complex(rk), intent(in)           :: zeta
    !
    complex(rk) :: cna( 1:asymp_order)   ! Asymptotic outgoing solution
    complex(rk) :: cnb( 1:asymp_order)   ! Asymptotic incoming solution
    complex(rk) :: gra(-1:asymp_order+1) ! Gradient of the asymptotic outgoing solution
    complex(rk) :: grb(-1:asymp_order+1) ! Gradient of the asymptotic incoming solution
    integer(ik) :: npts
    real(rk)    :: eta                   ! Matching point
    complex(rk) :: va, ga, vb, gb        ! Values of the asymptotic solutions and gradients at the matching point
    complex(rk) :: wrn                   ! Wronskian of the asymptotic solution
    complex(rk) :: ewrn                  ! Error in the numerical Wronskian at the matching point
    complex(rk) :: vs, gs                ! Function and gradient in scaled coordinates
    real(rk)    :: scl                   ! Scale factor for the incoming/outgoing amplitudes
    !
    call TimerStart('Continuum direction')
    !
    !  Integrate Schroedinger equation from the origin
    !
    call integrate_radial(state,energy,zeta,efield)
    !
    !  Form asymptotic solutions. These soltuions contain an additional factor eta**(1/2),
    !  which is needed to eliminate gradient term from the Hamiltonian.
    !
    call asymptotic_solution(efield,energy,zeta, 1._rk,cna)
    call asymptotic_solution(efield,energy,zeta,-1._rk,cnb)
    call asymptotic_gradient(efield,energy, 1._rk,cna,gra)
    call asymptotic_gradient(efield,energy,-1._rk,cnb,grb)
    npts = size(state%wf,dim=2)
    eta  = state%r(npts)
    va   = evaluate_asymptote(efield,energy, 1._rk,0,cna,eta)
    ga   = evaluate_asymptote(efield,energy, 1._rk,2,gra,eta)
    vb   = evaluate_asymptote(efield,energy,-1._rk,0,cnb,eta)
    gb   = evaluate_asymptote(efield,energy,-1._rk,2,grb,eta)
    wrn  = ga * vb - gb * va
    ewrn = wrn - (0,2)*sqrt(efield)
    !
    !  Convert origin-based solution to the scaled coordinates. 
    !  Don't forget the overall scale factor.
    !
    vs  = sqrt(eta)*state%wf(0,npts) 
    gs  = sqrt(eta)*state%wf(1,npts) + 0.5_rk*state%wf(0,npts)/sqrt(eta)
    scl = scale(1._rk,state%wf_scale(npts))
    !
    !  Match the solutions, and we are done
    !
    state%c_out = scl*(gs*vb - vs*gb)/wrn
    state%c_in  = scl*(gs*va - vs*ga)/wrn
    !
    if (verbose>=2) then
      write (out,"()")
      write (out,"('            Matching at eta = ',g44.34e3)") eta
      write (out,"('          Outgoing solution = ',2g44.34e3)") va
      write (out,"(' Outgoing solution gradient = ',2g44.34e3)") ga
      write (out,"('          Incoming solution = ',2g44.34e3)") vb
      write (out,"(' Incoming solution gradient = ',2g44.34e3)") gb
      write (out,"('                  Wronskian = ',2g44.34e3)") wrn
      write (out,"('     Error in the Wronskian = ',2g44.34e3)") ewrn
    end if
    !
    if (verbose>=0) then
      write (out,"()")
      write (out,"('Outgoing solution amplitude = ',2g44.34e3)") state%c_out
      write (out,"('Incoming solution amplitude = ',2g44.34e3)") state%c_in
      write (out,"()")
      flush (out)
    end if
    !
    if (abs(ewrn)>=20*spacing(sqrt(efield))) then
      write (out,"('WARNING: Wronskian deviates from the exact solution [2*I*sqrt(efield)] by ',g14.6e3)") abs(ewrn)
      write (out,"('WARNING: Incoming and outgoing amplitudes are inaccurate')")
      flush (out)
    end if
    !
    call TimerStop('Continuum direction')
  end subroutine solve_continuum
  !
  subroutine dump_wavefunction(tag,file,state)
    character(len=*), intent(in)   :: tag      ! Descriptive tag
    character(len=*), intent(in)   :: file     ! Name of the file
    type(wavefunction), intent(in) :: state    ! Wavefunction to report
    !
    integer(ik) :: ipt, max_order
    integer(ik) :: iu
    !
    if (file==' ') return
    !
    open (newunit=iu,form='formatted',recl=1024,action='write',position='rewind',status='replace',file=trim(file))
    write (iu,"('# ',a)") trim(comment)
    write (iu,"('# ',a)") trim(tag)
    write (iu,"('#              E = ',2g36.26e3)") state%energy
    write (iu,"('#              Z = ',2g36.26e3)") state%zeta
    write (iu,"('#              F = ',g36.26e3)") state%efield
    write (iu,"('#  Normalization = ',a)") state%norm
    write (iu,"('#       ipt_stop = ',i0)") state%ipt_stop
    write (iu,"('#    R(ipt_stop) = ',g36.26e3)") state%r(state%ipt_stop)
    write (iu,"('#       FP radix = ',i0)") radix(1._rk)
    write (iu,"('#')")
    write (iu,"('#',a8,1x,a26,1x,a10,5(1x,a26))") &
           ' I ', ' R ', ' Exp ', ' Re[wf] ', ' Im[wf] ', ' Re[d wf/d R] ', ' Im[d wf/d R] ', '...'
    max_order = min(4,ubound(state%wf,dim=1))
    dump_points: do ipt=1,size(state%r)
      write (iu,"(1x,i8,1x,g26.16e3,1x,i10,5(2x,g26.16e3,1x,g26.16e3))") &
             ipt, state%r(ipt), state%wf_scale(ipt), state%wf(0:max_order,ipt)
    end do dump_points
    close (iu)
  end subroutine dump_wavefunction
  !
  subroutine solve_single_energy
    type(wavefunction) :: bound
    type(wavefunction) :: st_up, st_dn
    type(wavefunction) :: continuum
    complex(rk)        :: zeta1, zeta2
    !
    call init_wavefunction(bound,     xi2_max, xi_npts, xi_maxorder)
    call init_wavefunction(st_up,     xi2_max, xi_npts, xi_maxorder)
    call init_wavefunction(st_dn,     xi2_max, xi_npts, xi_maxorder)
    call init_wavefunction(continuum,eta2_max,eta_npts,eta_maxorder)
    !
    call solve_bound(verbose,bound,st_up,st_dn,energy_single)
    call dump_wavefunction('bound',file_bound,bound)
    zeta1 = bound%zeta
    zeta2 = 4*znuc - zeta1
    call solve_continuum(verbose,continuum,energy_single,zeta2)
    call dump_wavefunction('continuum',file_continuum,continuum)
    !
    call destroy_wavefunction(continuum)
    call destroy_wavefunction(st_dn)
    call destroy_wavefunction(st_up)
    call destroy_wavefunction(bound)
  end subroutine solve_single_energy
  !
  subroutine solve_outgoing
    type(wavefunction) :: bound
    type(wavefunction) :: st_up, st_dn
    type(wavefunction) :: continuum
    type (fr_state)    :: root
    real(rk)           :: guess(2)
    !
    call TimerStart('Outgoing solution')
    !
    call init_wavefunction(bound,     xi2_max, xi_npts, xi_maxorder)
    call init_wavefunction(st_up,     xi2_max, xi_npts, xi_maxorder)
    call init_wavefunction(st_dn,     xi2_max, xi_npts, xi_maxorder)
    call init_wavefunction(continuum,eta2_max,eta_npts,eta_maxorder)
    !
    guess(1) =  real(energy_guess,kind=rk)
    guess(2) = aimag(energy_guess)
    call fr_initialize(verbose,root,neqs=2_ik,vars=guess,epsstep=energy_tol,diffstep=energy_tol)
    !
    root_iterations: do
      call fr_step(root)
      select case (root%action)
        case default
          write (out,"('hydrogen_tunnel_v2%solve_outgoing: unexpected request from equation solver: ',i0)") root%action
          stop 'hydrogen_tunnel_v2%solve_outgoing - state machine corrupted'
        case (fr_act_evaluate)
          call evaluate_function
        case (fr_act_done)
          if (verbose>=1) write (out,"('Optimization of outgoing solution converged.')")
          exit root_iterations
        case (fr_act_failed)
          write (out,"('WARNING: Optimization of outgoing solution failed to converge.')")
          write (out,"('WARNING: Reporting unconverged values below.')")
          flush (out)
          exit root_iterations
      end select
    end do root_iterations
    call evaluate_function
    !
    call dump_wavefunction('bound',file_bound,bound)
    call dump_wavefunction('continuum',file_continuum,continuum)
    !
    write (out,"()")
    write (out,"('Outgoing solution is at the energy = ',2g44.34e3)") bound%energy
    write (out,"('                Outgoing amplitude = ',2g44.34e3)") continuum%c_out
    write (out,"('                Incoming amplitude = ',2g44.34e3)") continuum%c_in
    write (out,"()")
    flush (out)
    !
    call fr_finalize(root)
    !
    call destroy_wavefunction(continuum)
    call destroy_wavefunction(st_dn)
    call destroy_wavefunction(st_up)
    call destroy_wavefunction(bound)
    call TimerStop('Outgoing solution')
    !
    contains
    !
    subroutine evaluate_function
      complex(rk) :: energy
      complex(rk) :: zeta2
      real(rk)    :: f1, f2
      !
      energy = cmplx(root%vars(1),root%vars(2),kind=rk)
      call solve_bound(verbose-1,bound,st_up,st_dn,energy)
      zeta2 = 4*znuc - bound%zeta
      call solve_continuum(verbose-1,continuum,energy,zeta2)
      f1   =  real(continuum%c_in,kind=rk)
      f2   = aimag(continuum%c_in)
      root%eqs(:) = (/ f1, f2 /)
    end subroutine evaluate_function
  end subroutine solve_outgoing
  !
  subroutine start
    use ISO_FORTRAN_ENV
    !$ use OMP_LIB
    external :: versions
    !
    write (out,"('Version: ',a/)") __BUILD_ID__
    write (out,"('Compiler: ',a)") compiler_version()
    write (out,"('Build flags: ',a/)") compiler_options()
    !
    write (out,"('  Integer kind = ',i0,' (',i0,' decimals)')") kind(1_ik), int(log10(huge(1_ik)+1._rk))
    write (out,"('     Real kind = ',i0,' (',i0,' decimals)')") kind(1._rk), precision(1._rk)
    write (out,"('Aux. real kind = ',i0,' (',i0,' decimals)')") kind(1._xk), precision(1._xk)
    write (out,"()")
    !
    write (out,"(/t5,a)") trim(rcsid_hydrogen_tunnel_v2)
    call versions
    !
    call TimerStart('start')
    !
    read (input,nml=ht_v2)
    write (out,"(' ===== begin simulation parameters ===== ')")
    write (out,nml=ht_v2)
    write (out,"(' ====== end simulation parameters ====== ')")
    !
    write (out,"(/a/)") trim(comment)
    !
    !$ if (omp_num_threads/=0) then
    !$   write (out,"('WARNING: Please set OMP_THREAD_LIMIT=',i0,' environment variable')") omp_num_threads
! omp target is a very recent addition to the OpenMP spec; it is not supported before gfortran 12
!   !$   write (out,"('Forcing maximum number of OpenMP threads to ',i0)") omp_num_threads
!   !$omp target thread_limit(omp_num_threads)
!   !$omp end target 
    !$ end if
    !$ call omp_set_nested(.true.)
    !
    select case (task)
      case default
        stop 'hydrogen_tunnel_v2%start - bad task'
      case ('energy')
        call solve_single_energy
      case ('outgoing')
        call solve_outgoing
    end select
    !
    call TimerStop('start')
    call TimerReport
  end subroutine start
end module hydrogen_tunnel_v2
!
program main
  use hydrogen_tunnel_v2
  call start
end program main

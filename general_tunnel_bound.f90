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
!  Routines dealing with the bound coordinate in general_tunnel.f90
!
module general_tunnel_bound
  use accuracy
  use constants
  use derivative_tools
  use find_root
  use general_tunnel_data
  use general_tunnel_dump
  use general_tunnel_potential
  use lapack
  use math
  use poly_tools
  use tridiagonal_cmtql1
  use timer
  !$ use OMP_LIB
  implicit none
  private
  public bound_overlap, guess_zeta1, channel_prepare_xi
  public check_channel_overlaps
  public rcsid_general_tunnel_bound
  !
  character(len=clen), save :: rcsid_general_tunnel_bound = "$Id: $"
  !
  contains
  !
  real(rk) function wavefunction_eps ()
    if (wavefunction_tol>0) then
      wavefunction_eps = wavefunction_tol
    else
      wavefunction_eps = abs(wavefunction_tol)*spacing(1._rk)
    end if
  end function wavefunction_eps
  !
  !  Unlike Kolosov's hydrogen code, it is no longer sensible to keep radial integrators
  !  for the bound and continuum solutions unified.
  !
  !  This routine accounts for the bulk of our runtime for general potentials (and almost
  !  half for hydrogenic solutions), so we care about the speed. In order to simplify
  !  cut-offs for the power series and prevent avoidable overflows, we'll fold the grid 
  !  step into the power series coefficients. To make this process efficient, we have to 
  !  limit ourselves to uniform grids.
  !
  subroutine integrate_radial_bound(state,eta_pt,energy,zeta,efield)
    type (wavefunction), intent(inout) :: state
    integer(ik), intent(in)            :: eta_pt ! The continuum grid point
    complex(rk), intent(in)            :: energy ! Energy of the solution
    complex(rk), intent(in)            :: zeta   ! Separation parameter
    real(rk), intent(in)               :: efield ! Electric field. Negative values correspond to discrete-spectrum solutions
    !
    integer(ik) :: npts, maxorder
    integer(ik) :: ipt, iord, idm
    complex(rk) :: temp
    real(rk)    :: xi
    complex(rk) :: c2, c3, c4, c5, c6, c7        ! Fixed part of the recursion coefficients
    real(rk)    :: udm(0:ubound(state%wf,dim=1)) ! Power series of 2 u(r) 
    real(rk)    :: dx(0:ubound(state%wf,dim=1))  ! Powers of the grid spacing; needed to avoid repeated
                                                 ! evaluation in calculating the recursions.
    real(rk)    :: abs_c                         ! Absolute value of the power series coefficient
    real(rk)    :: max_c                         ! Largest power series coefficient seen so far
    real(rk)    :: max_wf                        ! Largest wavefunction value so far
    integer(ik) :: zero_c_count                  ! Number of consequitive negligible power series coefficients
  ! real(rk)    :: eps                           ! Tolerance on the wavefunction and derivative
    !
    !  Store solution parameters
    !
    state%energy = energy
    state%zeta   = zeta
    state%efield = efield
    state%norm   = 'natural'
    !
    npts         = state%ipt_max
    maxorder     = ubound(state%wf,dim=1)
  ! eps          = wavefunction_eps()
    !
    dx(0) = 1._rk
    dx(1) = state%r(2) - state%r(1)
    !
    fill_step_powers: do iord=2,maxorder
      dx(iord) = dx(iord/2) * dx(iord-iord/2) ! Minimize error propagation by using binary-tree formula
    end do fill_step_powers
    !
    state%iord(:) = maxorder
    !
    !  The origin of the grid requires special handling due to the
    !  presence of singularities in the Hamiltonian.
    !
    state%wf (   :,1) = 0
    state%wf (mval,1) = dx(mval)
    ! call potential_u_halfx2(eta,0._rk,maxorder,udm)
    ! udm(:) = 4.0_rk * udm(:)
    udm(:) = 4.0_rk * tab_u_halfx2(:,1,eta_pt)
    !
    max_c = abs(state%wf(mval,1))
    zero_c_count = 0
    fill_origin: do iord=mval+1,maxorder
      temp = 0
      if (iord>=2) temp = temp + dx(2) * state%wf(iord-2,1) * state%zeta
      if (iord>=4) temp = temp + dx(4) * state%wf(iord-4,1) * state%energy * 2
      if (iord>=6) temp = temp + dx(6) * state%wf(iord-6,1) * state%efield
      potential_at_origin: do idm=0,iord-2
        temp = temp + dx(2+idm) * state%wf(iord-2-idm,1) * udm(idm)
      end do potential_at_origin
      state%wf(iord,1) = -temp / (iord**2 - mval**2)
      !
      !  Terminate power series early if possible
      !
      abs_c = abs(state%wf(iord,1))
      if (abs_c>max_c) max_c = abs_c
      if (abs_c>=spacing(max_c)) then
        zero_c_count = 0
      else
        zero_c_count = zero_c_count + 1
        if (zero_c_count>=6) then
          state%iord(1) = iord
          state%wf(iord+1:,1) = 0
          exit fill_origin
        end if
      end if
    end do fill_origin
    !
    !  We can use general code for the rest of the radial grid
    !
    max_wf = 0
    state%ipt_stop = state%ipt_max
    propagate: do ipt=2,npts
      !
      ! Use power-series expansion at the previous point to get wavefunction
      ! and its gradient. Since we've scaled our derivatives, the effective
      ! step size is unity.
      !
      iord = state%iord(ipt-1)
      call poly_power_series(state%wf(:iord,ipt-1),1._rk,state%wf(0:1,ipt))
      !
      max_wf = max(max_wf,abs(state%wf(0,ipt))**2)
    ! 27.08.2022: Do not stop the integration; we'll use perturbation analysis
    !             to determine the effective end of the grid.
    ! !
    ! ! Check for the possible early exit if wavefunction and gradient vanished
    ! !
    ! if (maxval(abs(state%wf(0:1,ipt))**2)<eps*max_wf) then
    !    state%ipt_stop = ipt
    !    state%wf(2:,ipt) = 0
    !    state%iord(ipt) = 1
    !    exit propagate
    ! end if
      !
      ! Calculate power-series coefficients at the new point
      !
      xi = state%r(ipt)
      ! call potential_u_halfx2(eta,xi,maxorder,udm)
      ! udm(:) = 4.0_rk * udm(:)
      udm(:) = 4.0_rk * tab_u_halfx2(:,ipt,eta_pt)
      c2 = -(mval**2) + xi**2*state%zeta + 2*xi**4*state%energy + xi**6*state%efield
      c3 = 2*xi*state%zeta + 8*xi**3*state%energy + 6*xi**5*state%efield
      c4 = state%zeta + 12*xi**2*state%energy + 15*xi**4*state%efield
      c5 = 8*xi*state%energy + 20*xi**3*state%efield
      c6 = 2*state%energy + 15*xi**2*state%efield
      c7 = 6*xi*state%efield
      !
      !  From the second term on, we can use the general expression
      !
      max_c = maxval(abs(state%wf(0:1,ipt)))
      zero_c_count = 0
      propagate_power: do iord=2,maxorder
                     temp =        state%wf(iord-1,ipt) * dx(1) * xi*(iord-1)*(2*iord-3)
        if (iord>=2) temp = temp + state%wf(iord-2,ipt) * dx(2) * ((iord-2)**2 + c2)
        if (iord>=3) temp = temp + state%wf(iord-3,ipt) * dx(3) * c3
        if (iord>=4) temp = temp + state%wf(iord-4,ipt) * dx(4) * c4
        if (iord>=5) temp = temp + state%wf(iord-5,ipt) * dx(5) * c5
        if (iord>=6) temp = temp + state%wf(iord-6,ipt) * dx(6) * c6
        if (iord>=7) temp = temp + state%wf(iord-7,ipt) * dx(7) * c7
        if (iord>=8) temp = temp + state%wf(iord-8,ipt) * dx(8) * state%efield
        potential_1: do idm=0,iord-2
          temp = temp + xi**2 * dx(2+idm) * state%wf(iord-2-idm,ipt) * udm(idm)
        end do potential_1
        potential_2: do idm=0,iord-3
          temp = temp + 2._rk * xi * dx(3+idm) * state%wf(iord-3-idm,ipt) * udm(idm)
        end do potential_2
        potential_3: do idm=0,iord-4
          temp = temp + dx(4+idm) * state%wf(iord-4-idm,ipt) * udm(idm)
        end do potential_3
        !
        state%wf(iord,ipt) = - temp / (iord*(iord-1)*xi**2)
        !
        !  Terminate power series early if possible
        !
        abs_c = abs(state%wf(iord,ipt))
        if (abs_c>max_c) max_c = abs_c
        if (abs_c>=spacing(max_c)) then
          zero_c_count = 0
        else
          zero_c_count = zero_c_count + 1
          if (zero_c_count>=8) then
            state%iord(ipt) = iord
            state%wf(iord+1:,ipt) = 0
            exit propagate_power
          end if
        end if
      end do propagate_power
    end do propagate
    !
    !  Clear all points beyond the final solution point; they are zero for our purposes.
    !
    state%wf(:,state%ipt_stop+1:) = 0
  end subroutine integrate_radial_bound
  !
  !  Find the value of the separation parameter which yields a bound solution at
  !  a given energy. 
  !
  subroutine solve_bound(verbose,state,st_up,st_dn,eta_pt,energy)
    integer(ik), intent(in)           :: verbose
    type(wavefunction), intent(inout) :: state     ! state%zeta must already contain a reasonable 
                                                   ! guess for the separation parameter
    type(wavefunction), intent(inout) :: st_up     ! Perturbed wavefunctions for stability analysis
    type(wavefunction), intent(inout) :: st_dn     ! 
    integer(ik), intent(in)           :: eta_pt    ! Continuum grid point we are solving at
    complex(rk), intent(in)           :: energy    ! Energy of the solution.
    !
    type (fr_state) :: root
    real(rk)        :: tail_re, tail_im
    real(rk)        :: vars_init(2)
    integer(ik)     :: npts
    !
    vars_init(1) =  real(state%zeta,kind=rk)
    vars_init(2) = aimag(state%zeta)
    call fr_initialize(verbose-10,root,neqs=2_ik,vars=vars_init,epsstep=zeta_tol,diffstep=zeta_tol)
    root_iterations: do
      call fr_step(root)
      select case (root%action)
        case default
          write (out,"('general_tunnel%solve_bound: unexpected request from equation solver: ',i0)") root%action
          stop 'general_tunnel%solve_bound - state machine corrupted'
        case (fr_act_evaluate)
          call evaluate_function
        case (fr_act_done)
          if (verbose>=3) write (out,"('Optimization of separation parameter zeta1 converged.')")
          exit root_iterations
        case (fr_act_failed)
          write (out,"('WARNING: Optimization of separation parameter zeta1 failed to converge.')")
          write (out,"('WARNING: Continuing with possibly unconverged zeta1 value.')")
          call flush_wrapper (out)
          exit root_iterations
      end select
    end do root_iterations
    call evaluate_function
    call fr_finalize(root)
    !
    if (verbose>=3) then
      npts    = state%ipt_stop ! May be less than the full size of the array!
      tail_re =  real(state%wf(0,npts),kind=rk)
      tail_im = aimag(state%wf(0,npts))
      write (out,"(/'Final bound-state solution for eta point = ',i0,' = ',g44.34e3)") eta_pt, state%r(npts)
      write (out,"('       Guess zeta1 = ',2g44.34e3)") vars_init
      write (out,"('   Optimized zeta1 = ',2g44.34e3)") state%zeta
      write (out,"('            Energy = ',2g44.34e3)") state%energy
      write (out,"(' Wavefunction tail = ',2g44.34e3)") tail_re, tail_im
      call flush_wrapper (out)
    end if
    !
    contains
    !
    !  Our optimization goal is the wavefunction at the last grid point
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
      call integrate_radial_bound(state,eta_pt,energy,zeta,   -efield)
      !$omp section  ! The perturbed-zeta solutions
      call integrate_radial_bound(st_up,eta_pt,energy,zeta_up,-efield)
      !$omp section  
      call integrate_radial_bound(st_dn,eta_pt,energy,zeta_dn,-efield)
      !$omp end parallel sections
      npts = state%ipt_max
      state%ipt_stop = npts  ! Start by assuming that the entire usable grid is OK
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
      f1   =  real(state%wf(0,npts),kind=rk)
      f2   = aimag(state%wf(0,npts))
      root%eqs(:) = (/ f1, f2 /)
      if (verbose>=4) then
        write (out,"('Z1 = ',2g32.20e3,' WF = ',2g32.20e3)") zeta, root%eqs
      end if
    end subroutine evaluate_function
  end subroutine solve_bound
  !
  !  Calculate overlap between two bound solutions. Because we work with a complex
  !  symmetric Hamiltonian (rather than the usual Hermitian), no conjugation is
  !  involved, and bound_overlap(f1,f2) = bound_overlap(f2,f1).
  !
  function bound_overlap(f1,f2) result(ov)
    type(wavefunction), intent(in) :: f1, f2 ! The two wavefunctions
    complex(rk)                    :: ov
    !
    integer(ik) :: ipt, npts, nder
    complex(rk) :: poly1(2*ubound(f1%wf,dim=1)+3) ! Buffer for the product polynomial
    complex(rk) :: poly2(2*ubound(f1%wf,dim=1)+3) ! Buffer for the product polynomial
    complex(rk) :: term(1), ve(0:1)
    !
    if (any(lbound(f1%wf)/=(/0,1/)) .or. any(lbound(f2%wf)/=(/0,1/)) .or. any(ubound(f1%wf)/=ubound(f2%wf))) then
      stop 'general_tunnel%bound_overlap - bad call'
    end if
    !
    ! npts = size(f1%wf,dim=2)
    ! We do not want to integrate past the point where the solution lost significance!
    npts = min(f1%ipt_stop,f2%ipt_stop)
    nder = ubound(f1%wf,dim=1)
    !
    !  Our wavefunctions are piece-wise polynomials. 
    !
    ov = 0._rk
    process_pieces: do ipt=1,npts-1
      ! Multiply polynomials representing the wavefunctions
      call poly_multiply(f1%wf(:,ipt),f2%wf(:,ipt),poly1(:2*nder+1))
      ! Multiply by the volume element
      ve(0) = f1%r(ipt)
      ve(1) = f1%r(2)   ! That's the grid spacing
      call poly_multiply(poly1(:2*nder+1),ve,poly2(:2*nder+2))
      ! Integrate the polynomial
      call poly_integrate(poly2(:2*nder+2),poly1(:))
      ! Accumulate the integral. The effective grid spacing is always unity
      call poly_power_series(poly1,1._rk,term)
      ov  = ov + term(1)
    end do process_pieces
    ov = ov * f1%r(2) ! Multiply by the uniform grid spacing, and we are done.
  end function bound_overlap
  !
  !  Guess initial values of the separation parameter zeta1.
  !  Zeta1 are the discrete eigenvalues of one-dimensional eigenproblem along the bound coordinate
  !  (Eq. 12 of [1]). Because the guess does not need to be terribly accurate (it will be refined
  !  later), we will discretize the eigenproblem using naive 3-point formulae, and solve it using
  !  brute-force diagonalization.
  !
  !  We need to find the eigenvalues of the operator:
  !
  !  G = (d^2/d x^2) + (1/x)(d/d x) - m^2/x^2 - Fx^4 + 2Ex^2 + 2u(y^2/2+x^2/2)
  !
  !  Due to the presence of the (d/d x) term, this operator is not very friendly to
  !  tri-diagonal solvers. Because we don't care about the eigenvectors at all, we'll
  !  instead look for the eigenvalues of the operator:
  !
  !  G' = (d^2/d x^2) - (m^2-1/4)/x^2 - Fx^4 + 2Ex^2 + 2u(y^2/2+x^2/2)
  !
  !  G' and G have identical eigenvalues, since:
  !
  !  x^(1/2) G [ x^(-1/2) g(x) ] = G' g(x)
  !
  !  However, G' is (complex!) symmetric tri-diagonal, and hence is easier to diagonalize.
  !
  subroutine guess_zeta1(verbose,en,eta,z1_tab)
    integer(ik), intent(in)  :: verbose
    complex(rk), intent(in)  :: en          ! Energy of the solution
    real(rk), intent(in)     :: eta         ! Value of the continuum coordinate
    complex(rk), intent(out) :: z1_tab(:)   ! Guess at the separation parameter, n_channel values.
                                            ! Note that we skip over the lower base_channel-1 channels.
    !
    real(rk)                 :: d                    ! Uniform grid spacing. Grid points are at d*i, i=1 .. xi_npts_guess
    complex(rk)              :: gp_d(xi_npts_guess)  ! Diagonal part of the discretized G'
    complex(rk)              :: gp_e(xi_npts_guess)  ! Sub-diagonal part of the discretized G'
    integer(ik)              :: i
    real(rk)                 :: xi
    integer(ik)              :: ierr
    real(rk)                 :: ptab(0:0)
    !
    if (xi_npts_guess<n_channels+base_channel-1) then
      stop 'general_tunnel%guess_zeta1 - xi_npts_guess must be greater than (n_channels+base_channel-1)'
    end if
    !
    d = xi_max / xi_npts_guess
    !
    !  Prepare linear system matrix. The sub-/super-diagonal elements are real. The diagonal
    !  is generally complex. The kinetic energy operator at the first point requires special
    !  treatment: we need to satisfy the cusp condition at the origin:
    !
    !    g(x) = C x**(m+1/2), when x->0
    !
    !  Additionally, we want to keep the Hamiltonian matrix symmetric (this simplifies 
    !  diagonalization). Taken together, these two conditions fix the discretized
    !  Laplacian coefficient to the somewhat ungainly expression below.
    !
    !  The accuracy of this Laplacian is not great in the m=0 channel. However, it is
    !  sufficient to give us the starting point for optimizing the separation
    !  parameter using the shooting method.
    !
    fill_diagonal: do i=1,xi_npts_guess
      xi = d * i
      if (i==1) then
        gp_d(i) = (mval**2-0.25_rk-2._rk**(mval+0.5_rk))/d**2
      else
        gp_d(i) = -2.0_rk/d**2 
      end if
      ! call potential_u_r(0.5_rk*(eta**2+xi**2),0_ik,ptab)
      call potential_u_halfx2(eta,xi,0_ik,ptab)
      gp_d(i) = gp_d(i) - (mval**2-0.25_rk)/xi**2 - efield*xi**4 + 2*en*xi**2 + 4._rk * ptab(0)
    end do fill_diagonal
    gp_e(:) = 1.0_rk/d**2
    !
    call cmtql1(gp_d,gp_e,ierr)
    if (ierr/=0) then
      write (out,"('guess_zeta1: Diagonalization of the trial Hamiltonian failed. Error core = ',i0)") ierr
      stop 'general_tunnel%guess_zeta1 - diagonalization failed'
    end if
    !
    !  Our version of cmtql1 returns eigenvalues in the order of increasing real part.
    !  Take the last (n_channels) values with the sign changed, and we are done.
    !
    pick_zeta1_guesses: do i=1,n_channels
      z1_tab(i) = -gp_d(xi_npts_guess-(base_channel-1)-(i-1))
    end do pick_zeta1_guesses
    !
    if (verbose>=3) then
      write (out,"(/t5,'Initial guess for channel separation parameters zeta1 at eta = ',g20.12/)") eta
      write (out,"((1x,a8,1x,2(1x,a20)))") 'Channel', 'Re[zeta1]', 'Im[zeta1]', &
                                           '-------', '---------', '---------'
      print_zeta1_guesses: do i=1,n_channels
        write (out,"(1x,i8,1x,2(1x,g20.12))") i, z1_tab(i)
      end do print_zeta1_guesses
      write (out,"()") 
    end if
  end subroutine guess_zeta1
  !
  function guess_integration_boundary(verbose,eta,energy,efield,zeta) result(ind)
    integer(ik), intent(in) :: verbose
    real(rk), intent(in)    :: eta    ! Value of the unbound coourdinate
    complex(rk), intent(in) :: energy ! Energy of the solution
    real(rk), intent(in)    :: efield ! Static field
    complex(rk), intent(in) :: zeta   ! Separation parameter
    integer(ik)             :: ind    ! Index within xi_tab() where we should stop integration
    !
    real(rk), parameter    :: eps = 1e-3_rk
    integer(ik), parameter :: max_iter = 100_ik
    type(fr_state)         :: root
    real(rk)               :: x0, x
    !
    !  To find the outer turning point, let's solve the equation:
    !
    !    Re[-(m**2-1/4)/x**2 - F*x**4 + 2E*x**2 + 4*u((x**2+eta**2)/2) + zeta] == 0
    !
    !  (this is the potential we use in guess_zeta1)
    !
    !  We can use the outer grid boundary as the starting point. We do not need terribly
    !  high accuracy here.
    !
    call fr_initialize(verbose-10_ik,root,neqs=1_ik,vars=(/xi_max/),epsstep=eps,diffstep=eps,maxiter=max_iter)
    turn_iterations: do
      call fr_step(root)
      select case (root%action)
        case default
          stop 'general_tunnel%guess_integration_boundary - state machine corrupted (1)'
        case (fr_act_evaluate)
          root%eqs(1) = evaluate_potential(root%vars(1))
        case (fr_act_done,fr_act_stopped)
          x0 = root%vars(1)
          exit turn_iterations
        case (fr_act_failed)
          write (out,"('WARNING: Optimization of the outer turning point failed to converge. Using maximum.')")
          x0 = xi_max
          exit turn_iterations
      end select
    end do turn_iterations
    call fr_finalize(root)
    !
    !  From the outer turning point, the wavefunction decays asymptotically as
    !
    !   f(x) = (1/x) * exp[ -(sqrt(F)/3) * x**3 + (E/sqrt(F)) * x ]
    !
    !  We need to solve the equation:
    !
    !    f(x) <= eps * f(x0)
    !
    !  The easiest way is to take the logarithm and again to use the root finder. The starting 
    !  point for the search requires a bit of thinking: When E>0, the linear term will initially 
    !  dominate, which we may need to take into account.
    !
    x = x0
    if (real(energy,kind=rk)>0) then
      x = max(x0,sqrt(3._rk*real(energy,kind=rk)/efield))
    end if
    call fr_initialize(verbose-10_ik,root,neqs=1_ik,vars=(/x/),epsstep=eps,diffstep=eps,maxiter=max_iter)
    tail_iterations: do
      call fr_step(root)
      select case (root%action)
        case default
          stop 'general_tunnel%guess_integration_boundary - state machine corrupted (2)'
        case (fr_act_evaluate)
          root%eqs(1) = evaluate_tail(root%vars(1))
        case (fr_act_done,fr_act_stopped)
          x = root%vars(1)
          exit tail_iterations
        case (fr_act_failed)
          write (out,"('WARNING: Optimization of the asymptotic tail failed to converge. Using maximum.')")
          x = xi_max
          exit tail_iterations
      end select
    end do tail_iterations
    call fr_finalize(root)
    !
    !  Find the index. We know that our grid is linear and starts at zero ...
    !
    if (x>xi_tab(xi_npts)) then
      ind = xi_npts
    else
      ind = min(nint(x/xi_tab(2)),xi_npts)
    end if
    !
    if (verbose>=3) then
      write (out,"('Boundary E = ',2g16.6,' F = ',g16.6,' Z = ',2g16.6,' x0 = ',g16.6,' x = ',g16.6)") &
             energy, efield, zeta, x0, x
    end if
    !
    contains 
    function evaluate_potential(xi) result(v)
      real(rk), intent(in) :: xi
      real(rk)             :: v
      !
      real(rk)             :: ptab(0:0)
      !
      if (xi<=0) then
        v = 1e3_rk
      else
        ! call potential_u_r(0.5_rk*(eta**2+xi**2),0_ik,ptab)
        call potential_u_halfx2(eta,xi,0_ik,ptab)
        v = real(zeta - (mval**2-0.25_rk)/xi**2 - efield*xi**4 + 2*energy*xi**2 + 4._rk*ptab(0),kind=rk)
      end if
    end function evaluate_potential
    !
    function evaluate_tail(xi) result(v)
      real(rk), intent(in) :: xi
      real(rk)             :: v
      !
      real(rk) :: logfx0, logfx, logeps
      !
      logfx0 = -(sqrt(efield)/3._rk)*x0**3 + (real(energy,kind=rk)/sqrt(efield))*x0 - log(x0)
      logfx  = -(sqrt(efield)/3._rk)*xi**3 + (real(energy,kind=rk)/sqrt(efield))*xi - log(xi)
      logeps = log(wavefunction_eps())
      v = logfx - (logfx0 + logeps)
    end function evaluate_tail
    !
  end function guess_integration_boundary
  !
  subroutine prepare_bound_solution(verbose,bound,st_up,st_dn,eta_pt,energy)
    integer(ik), intent(in)           :: verbose
    type(wavefunction), intent(inout) :: bound       ! Desired bound solution
    type(wavefunction), intent(inout) :: st_up       ! Perturbed solutions for stability analysis
    type(wavefunction), intent(inout) :: st_dn       ! Perturbed solutions for stability analysis
    integer(ik), intent(in)           :: eta_pt      ! Continuum coordinate grid point
    complex(rk), intent(in)           :: energy      ! Solution energy
    !
    complex(rk) :: overlap
    integer(ik) :: ipass
    !
    !  Bound solutions may run out of precision in the separation parameter before reaching 
    !  the end of the grid. This causes the solutions beoynd this point to grow exponentially,
    !  spoiling the overall accuracy. The simple remedy is to stop integration at a point
    !  where the bound solution is known to be sufficiently small.
    !
    !  Because our initial guess for the separation parameter is not terribly good,
    !  especially in the m=0 channel, we'll adjust the boundary once the solution
    !  is found, and rerun the solver. The new guess is supposed to be much more
    !  accurate, so the extra time should be negligible.
    !
    solution_passes: do ipass=1,2
      bound%ipt_max = guess_integration_boundary(verbose,eta_tab(eta_pt),energy,efield,bound%zeta)
      call solve_bound(verbose,bound,st_up,st_dn,eta_pt,energy)
    end do solution_passes
    !
    !  Normalize the solution. We require square (not square modulus!) normalization.
    !
    overlap       = bound_overlap(bound,bound)
    bound%wf(:,:) = (1._rk/sqrt(overlap)) * bound%wf(:,:)
    bound%norm    = 'square'
  end subroutine prepare_bound_solution
  !
  subroutine copy_bound_solution(src,dst)
    type(wavefunction), intent(in)    :: src
    type(wavefunction), intent(inout) :: dst
    !
    if (size(src%r)/=size(dst%r) .or. any(ubound(src%wf)/=ubound(dst%wf))) then
      stop 'general_tunnel_bound%copy_bound_solution - bad call'
    end if
    !
    dst%energy   = src%energy
    dst%zeta     = src%zeta
    dst%efield   = src%efield
    dst%norm     = src%norm
    dst%c_in     = src%c_in
    dst%c_out    = src%c_out
    dst%ipt_max  = src%ipt_max
    dst%ipt_stop = src%ipt_stop
    ! We deliberately do not copy r, since it is a pointer
    dst%wf       = src%wf
    dst%iord     = src%iord
  end subroutine copy_bound_solution
  !
  !  Build solutions for the bound variable for each grid point of the
  !  dissociative channel. This may involve a fair bit of work - but
  !  all points are independent, so that we can use as many threads
  !  are we could scrounge. Having limited nested parallelism also
  !  helps, since we solve each bound channel 3 times to evaluate the
  !  stability limits.
  !
  subroutine channel_prepare_xi(verbose,ichan,chan,energy)
    integer(ik), intent(in)      :: verbose
    integer(ik), intent(in)      :: ichan  ! Channel index for chan below
    type(channel), intent(inout) :: chan   ! Channel we are working on
    complex(rk), intent(in)      :: energy ! Desired solution energu
    !
    integer(ik) :: ipt, ipt_stop
    !
    call TimerStart('Bound solutions')
    if (potential=='hydrogenic') then
      !
      !  For potential='hydrogenic', bound solutions at all eta points
      !  are identical. Save some work by only running the solver once,
      !  then copy the results. This means we run serially.
      !
      call prepare_bound_solution(verbose,chan%bound(1),chan%st_up(1),chan%st_dn(1),1_ik,energy)
      !$omp parallel do default(none) num_threads(eta_npts-1) shared(eta_npts,chan) private(ipt)
      eta_points_copy: do ipt=2,eta_npts
        call copy_bound_solution(chan%bound(1),chan%bound(ipt))
        call copy_bound_solution(chan%st_up(1),chan%st_up(ipt))
        call copy_bound_solution(chan%st_dn(1),chan%st_dn(ipt))
      end do eta_points_copy
      !$omp end parallel do
    else
      !
      !  General case - solutions are different at each eta point
      !
      !$omp parallel do default(none) num_threads(eta_npts) shared(eta_npts,verbose,chan,energy) private(ipt)
      eta_points: do ipt=1,eta_npts
        call prepare_bound_solution(verbose,chan%bound(ipt),chan%st_up(ipt),chan%st_dn(ipt),ipt,energy)
      end do eta_points
      !$omp end parallel do
    end if
    ! 
    if (verbose>=1) then
      !
      !  Report to standard output
      !
      call report_channel_bound_solutions(ichan,channels(ichan))
      !
      !  ... and to external files
      !
      call dump_bound_wavefunction(ichan,channels(ichan)) ! general_tunnel_dump.f90
    end if
    !
    !  Once the bound wavefunctions have been reported for debugging, we should zero 
    !  out the solutions past ipt_stop, because later parts of the code rely on it.
    !
    if (verbose>1) then
      write (out,"(/'Forcing bound wavefunctions past ipt_stop to zero.')")
    end if
    !$omp parallel do default(none) shared(eta_npts,channels,ichan) private(ipt,ipt_stop)
    zero_tail: do ipt=1,eta_npts
      ipt_stop = channels(ichan)%bound(ipt)%ipt_stop
      channels(ichan)%bound(ipt)%wf(:,ipt_stop+1:) = 0
    end do zero_tail 
    !$omp end parallel do
    !
    call TimerStop('Bound solutions')
  end subroutine channel_prepare_xi
  !
  subroutine report_channel_bound_solutions(ichan,chan)
    integer(ik), intent(in)   :: ichan ! Channel index
    type(channel), intent(in) :: chan  ! Channel data
    !
    integer(ik) :: eta_ipt, xi_end, xi_stop
    complex(rk) :: zeta, tail
    real(rk)    :: eta
    !
    write (out,"(/t5,'Bound solutions for channel ',i0,' (',i0,')'/)") ichan, chan%ichan
    write (out,"((1x,a5,1x,a26,1x,a5,1x,a5,1x,4(1x,a26)))") &
           ' I ', ' eta ', 'maxpt', 'stop', ' Re[zeta1] ', ' Im[zeta1] ', ' Re[Psi-N] ', ' Im[Psi-N] ', &
           '---', '-----', '-----', '----', '-----------', '-----------', '-----------', '-----------'
    scan_bound: do eta_ipt=1,eta_npts
      eta     = eta_tab(eta_ipt)
      zeta    = chan%bound(eta_ipt)%zeta
      xi_end  = chan%bound(eta_ipt)%ipt_max
      xi_stop = chan%bound(eta_ipt)%ipt_stop
      tail    = chan%bound(eta_ipt)%wf(0,xi_stop)
      write (out,"(1x,i5,1x,g26.16,1x,i5,1x,i5,1x,4(1x,g26.16e3))") eta_ipt, eta, xi_end, xi_stop, zeta, tail
    end do scan_bound
    write (out,"()")
  end subroutine report_channel_bound_solutions
  !
  subroutine check_channel_overlaps(verbose)
    integer(ik), intent(in)  :: verbose
    integer(ik)              :: ichan, jchan, ipt
    integer(ik)              :: max_pt, alloc
    complex(rk)              :: overlap, max_error
    complex(rk), allocatable :: ovtab(:,:,:)
    !
    !  A bit of sanity checking: at least, the overlaps must be diagonal
    !
    if (verbose>=1) then
       write (out,"(/t5,'Overlaps of bound channel solutions'/)")
       write (out,"((1x,a5,1x,a16,1x,a3,1x,a3,1x,2(1x,a26)))") &
              ' IPT ', ' ETA ', ' I ', ' J ', ' Re[(I,J)]-delta(I,J) ', ' Im[(I,J)] ', &
              '-----', '-----', '---', '---', '----------------------', '-----------'
    end if
    !
    max_pt = eta_npts
    if (potential=='hydrogenic') max_pt = 1
    !
    allocate (ovtab(n_channels,n_channels,max_pt),stat=alloc)
    if (alloc/=0)  then
      stop 'general_tunnel_bound%check_channel_overlaps - allocation failed'
    end if
    !
    !  Calculate overlaps, running in parallel as much as possible
    !
    !$omp parallel do collapse(3) default(none) shared(max_pt,n_channels,channels,ovtab) private(ipt,ichan,jchan)
    build_overlaps: do ipt=1,max_pt
      build_chan1: do ichan=1,n_channels
        build_chan2: do jchan=1,n_channels
          ovtab(ichan,jchan,ipt) = bound_overlap(channels(ichan)%bound(ipt),channels(jchan)%bound(ipt))
        end do build_chan2
      end do build_chan1
    end do build_overlaps
    !$omp end parallel do
    !
    !  Determine maximum error first, to know what to print
    ! 
    max_error = 0
    error_check_overlaps: do ipt=1,max_pt
      error_ro_chan1: do ichan=1,n_channels
        error_ro_chan2: do jchan=1,n_channels
          overlap = ovtab(ichan,jchan,ipt)
          if (ichan==jchan) then
            overlap = overlap - 1.0_rk
          end if
          if (abs(overlap)>abs(max_error)) max_error = overlap
        end do error_ro_chan2
      end do error_ro_chan1
    end do error_check_overlaps
    !
    !  Now we are good to print the stuff :)
    !
    check_overlaps: do ipt=1,max_pt
      ro_chan1: do ichan=1,n_channels
        ro_chan2: do jchan=1,n_channels
          overlap = ovtab(ichan,jchan,ipt)
          if (ichan==jchan) then
            overlap = overlap - 1.0_rk
          end if
          if (verbose>=2 .or. (verbose>=1 .and. abs(overlap)>=0.1_rk*abs(max_error)) ) then
            write (out,"(1x,i5,1x,f16.8,1x,i3,1x,i3,1x,2(1x,g26.16e3))") ipt, eta_tab(ipt), ichan, jchan, overlap
          end if
        end do ro_chan2
      end do ro_chan1
    end do check_overlaps
    if (abs(max_error)>spacing(100._rk)) then
      write (out,"(/'WARNING: Adiabatic bound solutions are not orthonormal. Max. error = ',2g16.6e3/)") max_error
      if (abs(max_error)>1e-5_rk) then
        write (out,"(/'ERROR: Unacceptable errors in adiabatic bound solutions. Calculation is meaningless.'/)")
        stop 'general_tunnel%solve_single_energy - accuracy lost'
      end if
    end if
    deallocate (ovtab)
  end subroutine check_channel_overlaps
end module general_tunnel_bound

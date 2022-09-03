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
!  Ionization of 1-electron system in combination of a general central potential
!  and uniform static field. The code mostly follows [Kolosov, J Phys B 20, 2359 (1987)], 
!  combined with the Ansatz from [Batishev, Tolstikhin, and Morishita, PRA 82, 023416 (2010)]
!
!  We are working in the squared-parabolic coordinates. In these coordinates, we have:
!
!   x = Xi*Eta*cos(phi)
!   y = Xi*Eta*sin(phi)
!   z = (Xi**2 - Eta**2)/2
!   r = (Xi**2 + Eta**2)/2
!
!  In squared-parabolic coordinates, pure Coulombic problems are separable. General central
!  potentials are not, and need to be expanded over multiple separable solutions ("channels")
!
!  For simplicity, we will use a uniform grid in Xi and Eta (this is not the same
!  choice used in "hydrogen_tunnel_v2.f90".
!
!  Please keep in mind that normalization of the solutions is not conventional - neither
!  the total density (Hermitian normalization) nor the square of the wavefunction
!  (Siegert/complex-symmetric normalization) integrate to one.
!
module general_tunnel
  use accuracy
  use constants
  use derivative_tools
  use find_root
  use find_minimum
  use general_tunnel_bound
  use general_tunnel_continuum
  use general_tunnel_data
  use general_tunnel_dump
  use general_tunnel_potential
  use general_tunnel_nonadiabatic
  use lapack
  use math
  use timer
  !$ use OMP_LIB
  implicit none
  private
  public start
  public rcsid_general_tunnel
  !
  character(len=clen), save :: rcsid_general_tunnel = "$Id: $"
  !
  contains
  !
  subroutine solve_single_energy(energy,verbose)
    complex(rk), intent(in) :: energy   ! Desired energy of the solution
    integer(ik), intent(in) :: verbose
    !
    integer(ik)  :: ichan, ipt
    complex(rk)  :: z1_guess(n_channels)
    !
    call TimerStart('Solve single energy')
    !
    !  Calculation of the initial guess for the separation parameter is shared
    !  between the channels, so let's do it once. The guess is obtained by
    !  diagonalizing a discretized, tri-diagonal Hamiltonian for the bound 
    !  coordinate. It is not terribly accurate, but is sufficient to bootstrap
    !  the (much more accurate) shooting solver.
    !
    call TimerStart('Guess zeta1')
    !
    !  For hydrogenic potential, the solution is separable. Save some time
    !  by guessing only once. In this case, we can just as well run serially.
    !
    if (potential=='hydrogenic') then
      call guess_zeta1(verbose,energy,eta_tab(1),z1_guess)  ! general_tunnel_bound.f90
      prepare_guess_h: do ipt=1,eta_npts
        fill_guess_h: do ichan=1,n_channels
          channels(ichan)%bound(ipt)%zeta = z1_guess(ichan)
        end do fill_guess_h
      end do prepare_guess_h
    else
      !$omp parallel do default(none) num_threads(eta_npts) &
      !$omp& shared(eta_npts,potential,verbose,energy,eta_tab,n_channels,channels) &
      !$omp& private(ipt,z1_guess,ichan)
      prepare_guess: do ipt=1,eta_npts
        call guess_zeta1(verbose,energy,eta_tab(ipt),z1_guess)
        fill_guess: do ichan=1,n_channels
          channels(ichan)%bound(ipt)%zeta = z1_guess(ichan)
        end do fill_guess
      end do prepare_guess
      !$omp end parallel do
    end if
    call TimerStop('Guess zeta1')
    !
    !  Loop below contains enough parallelism inside channel_prepare_xi
    !
    prepare_channel: do ichan=1,n_channels
      call channel_prepare_xi(verbose,ichan,channels(ichan),energy)  ! general_tunnel_bound.f90
    end do prepare_channel
    !
    !  check_channel_overlaps may abort if bound solutions are not orthonormal enough
    !
    call check_channel_overlaps(verbose)   ! general_tunnel_bound.f90
    !
    !  Evaluate adiabatic separation parameters, non-adiabatic terms, and their derivatives
    !
    call prepare_channel_coupling(verbose) ! general_tunnel_nonadiabatic.f90
    !
    call solve_continuum(verbose,energy)   ! general_tunnel_continuum.f90
    !
    if (verbose>=1) then
      write (out,"(/'Using solutions to generate ',a,' for plotting'/)") trim(file_total_mode)
      if (file_total_mode(1:8)=='reverse ') then
        !
        !  Reverse integration from the boundary requested; do it
        !
        call invert_continuum(verbose,energy)        ! general_tunnel_continuum.f90
      end if
      call dump_total_wavefunction(energy)           ! general_tunnel_dump.f90
      call dump_total_cartesian_wavefunction(energy) ! general_tunnel_dump.f90
    end if
    !
    call TimerStop('Solve single energy')
  end subroutine solve_single_energy
  !
  subroutine solve_outgoing_root
    type (fr_state)    :: root
    real(rk)           :: guess(2)
    complex(rk)        :: energy
    integer(ik)        :: nvars
    !
    call TimerStart('Outgoing solution (root)')
    !
    select case (task)
      case ('outgoing')
        guess(1) =  real(energy_guess,kind=rk)
        guess(2) = aimag(energy_guess)
        nvars    = 2
      case ('minflux real')
        guess(1) =  real(energy_guess,kind=rk)
        nvars    = 1
      case ('minflux imag')
        guess(1) = aimag(energy_guess)
        nvars    = 1
    end select
    call fr_initialize(verbose,root,neqs=2_ik*n_channels,vars=guess(:nvars),epsstep=energy_tol, &
                       diffstep=energy_tol,maxstep=energy_step)
    !
    root_iterations: do
      call fr_step(root)
      select case (root%action)
        case default
          write (out,"('general_tunnel%solve_outgoing: unexpected request from equation solver: ',i0)") root%action
          stop 'general_tunnel%solve_outgoing - state machine corrupted'
        case (fr_act_evaluate)
          call evaluate_function
        case (fr_act_done)
          if (verbose>=1) write (out,"('Optimization of outgoing solution converged.')")
          exit root_iterations
        case (fr_act_stopped)
          if (verbose>=1) write (out,"('Optimization of outgoing solution stopped due to zero search direction.')")
          exit root_iterations
        case (fr_act_failed)
          write (out,"('WARNING: Optimization of outgoing solution failed to converge.')")
          write (out,"('WARNING: Reporting unconverged values below.')")
          call flush_wrapper (out)
          exit root_iterations
      end select
    end do root_iterations
    !
    ! For the final evaluation, provide more detail (same as we would have for
    ! a single-energy calculation).
    !
    verbose = verbose + 1
    call evaluate_function
    call fr_finalize(root)
    !
    write (out,"()")
    write (out,"('    Outgoing solution is at the energy = ',2g44.34e3)") energy
    write (out,"('Outgoing amplitude in the main channel = ',2g44.34e3)") continuum(main_channel)%c_out
    write (out,"('Incoming amplitude in the main channel = ',2g44.34e3)") continuum(main_channel)%c_in
    write (out,"()")
    call flush_wrapper (out)
    !
    call TimerStop('Outgoing solution (root)')
    !
    contains
    !
    subroutine evaluate_function
      real(rk)    :: scl
      !
      select case (task)
        case ('outgoing')
          energy = cmplx(root%vars(1),root%vars(2),kind=rk)
        case ('minflux real')
          energy = cmplx(root%vars(1),aimag(energy_guess),kind=rk)
        case ('minflux imag')
          energy = cmplx(real(energy_guess,kind=rk),root%vars(1),kind=rk)
      end select
      !
      if (verbose>=1) then
        write (out,"('Solving for energy = ',2g34.25e3)") energy
      end if
      call solve_single_energy(energy,verbose-1)
      !
      !  The goal for optimization is a bit tricky: even when our main_channel dominates
      !  at the origin, most of the incoming flux may be in some other channel, leading to
      !  very slow optimization progress. 
      !
      if (verbose>=1) then
        scl = sqrt(sum(abs(continuum(:)%c_in)**2))
        write (out,"(/'At energy = ',2(g0.25,1x),'incoming = ',2(g0.20,1x),'(main) ',g0.20,' (total)')") &
                    energy, continuum(main_channel)%c_in, scl
        call flush_wrapper(out)
      end if
      root%eqs(  :n_channels) =  real(continuum(:)%c_in,kind=rk)
      root%eqs(n_channels+1:) = aimag(continuum(:)%c_in)
    end subroutine evaluate_function
  end subroutine solve_outgoing_root
  !
  subroutine solve_outgoing_min
    type (fm_state)    :: root
    real(rk)           :: guess(2)
    complex(rk)        :: energy
    integer(ik)        :: nvars
    !
    call TimerStart('Outgoing solution (min)')
    !
    select case (task)
      case ('outgoing')
        guess(1) =  real(energy_guess,kind=rk)
        guess(2) = aimag(energy_guess)
        nvars    = 2
      case ('minflux real')
        guess(1) =  real(energy_guess,kind=rk)
        nvars    = 1
      case ('minflux imag')
        guess(1) = aimag(energy_guess)
        nvars    = 1
    end select
    call fm_initialize(verbose,root,vars=guess(:nvars),epsstep=energy_tol,diffstep=energy_tol,maxstep=energy_step)
    !
    root_iterations: do
      call fm_step(root)
      select case (root%action)
        case default
          write (out,"('general_tunnel%solve_outgoing: unexpected request from equation solver: ',i0)") root%action
          stop 'general_tunnel%solve_outgoing - state machine corrupted'
        case (fm_act_evaluate)
          call evaluate_function
        case (fm_act_done)
          if (verbose>=1) write (out,"('Optimization of outgoing solution converged.')")
          exit root_iterations
        case (fm_act_failed)
          write (out,"('WARNING: Optimization of outgoing solution failed to converge.')")
          write (out,"('WARNING: Reporting unconverged values below.')")
          call flush_wrapper (out)
          exit root_iterations
      end select
    end do root_iterations
    !
    ! For the final evaluation, provide more detail (same as we would have for
    ! a single-energy calculation).
    !
    verbose = verbose + 1
    call evaluate_function
    call fm_finalize(root)
    !
    write (out,"()")
    write (out,"('    Outgoing solution is at the energy = ',2g44.34e3)") energy
    write (out,"('Outgoing amplitude in the main channel = ',2g44.34e3)") continuum(main_channel)%c_out
    write (out,"('Incoming amplitude in the main channel = ',2g44.34e3)") continuum(main_channel)%c_in
    write (out,"()")
    call flush_wrapper (out)
    !
    call TimerStop('Outgoing solution (min)')
    !
    contains
    !
    subroutine evaluate_function
      real(rk)    :: goal
      !
      select case (task)
        case ('outgoing')
          energy = cmplx(root%vars(1),root%vars(2),kind=rk)
        case ('minflux real')
          energy = cmplx(root%vars(1),aimag(energy_guess),kind=rk)
        case ('minflux imag')
          energy = cmplx(real(energy_guess,kind=rk),root%vars(1),kind=rk)
      end select
      !
      if (verbose>=1) then
        write (out,"('Solving for energy = ',2g34.25e3)") energy
      end if
      call solve_single_energy(energy,verbose-1)
      !
      goal = sum(abs(continuum(:)%c_in)**2)
      if (verbose>=1) then
        write (out,"(/'At energy = ',2(g0.25,1x),'incoming = ',2(g0.20,1x),'(main) ',g0.20,' (total)')") &
                    energy, continuum(main_channel)%c_in, goal
        call flush_wrapper(out)
      end if
      root%func = goal
    end subroutine evaluate_function
  end subroutine solve_outgoing_min
  !
  subroutine math_init
    real(rk) :: dummy
    !
    if (log10(huge(1._rk))>=4930._rk) then
      ! Quad precision; don't forget "factorial_slack"
      dummy = MathFactorial(1750_ik-5_ik)
      dummy = MathLogFactorial(20000_ik)
    else if (log10(huge(1._rk))>=308._rk) then
      ! Double precision
      dummy = MathFactorial(170_ik-5_ik)
      dummy = MathLogFactorial(10000_ik)
    else
      ! Single precision
      dummy = MathFactorial(34_ik-5_ik)
      dummy = MathLogFactorial(500_ik)
    end if
  end subroutine math_init
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
    write (out,"(/t5,a)") trim(rcsid_general_tunnel)
    call versions  ! versions.f90
    !
    call TimerStart('start')
    !
    read (input,nml=gt)
    write (out,"(' ===== begin simulation parameters ===== ')")
    write (out,nml=gt)
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
    !  Do common initialization
    !
    call math_init    
    !
    call init_global_arrays          ! general_tunnel_data.f90
    !
    call check_asymptotic_potential  ! general_tunnel_potential.f90
    !
    call cache_potential             ! general_tunnel_potential.f90
    !
    select case (boundary) 
      case default
        stop 'general_tunnel%start - bad boundary'
      case ('origin single','asymptotic','flux')
        if (main_channel<1 .or. main_channel>n_channels) stop 'general_tunnel%start - bad main_channel'
        inner_boundary(:)            = 0
        inner_boundary(main_channel) = 1
      case ('origin all')
        read(input,*) inner_boundary(:)
    end select
    !
    select case (task)
      case default
        stop 'general_tunnel%start - bad task'
      case ('energy')
        call solve_single_energy(energy_single,verbose)
      case ('outgoing','minflux real','minflux imag')
        select case (outgoing_solver)
          case default
            stop 'general_tunnel%start - bad outgoing_solver'
          case ('root')
            call solve_outgoing_root
          case ('min')
            call solve_outgoing_min
        end select
    end select
    !
    call TimerStop('start')
    call TimerReport
  end subroutine start
end module general_tunnel
!
program main
  use general_tunnel
  call start
end program main

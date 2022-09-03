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
!  Simple Newton root finder with line searches and backtracking
!  The implementation follows discussion in Numerical Recipies
!
!  Fortran does not offer a natural way of passing around an opaque
!  object, so that the implementation below has to use an explicit 
!  state machine. This looks weird, but makes the code general and
!  thread-safe :P
!
!  After fr_initialize() was called, all interaction with the calling
!  routine occurs through the fr_state data structure.
!
module find_root
  use accuracy
  use constants
  use lapack
  use timer
  implicit none
  private
  public fr_state 
  public fr_initialize, fr_step, fr_finalize
  public fr_act_evaluate, fr_act_none, fr_act_done, fr_act_stopped, fr_act_failed
  public rcsid_find_root
  !
  character(len=clen), save :: rcsid_find_root = "$Id: find_root.f90,v 1.15 2021/08/10 13:17:38 ps Exp ps $"
  !
  !  User activity codes (see fr_state%action below)
  !
  integer(ik), parameter :: fr_act_none     = 2300_ik ! No action required; call again
  integer(ik), parameter :: fr_act_evaluate = 2302_ik ! Evaluate equations for the given value of vars()
  integer(ik), parameter :: fr_act_done     = 2303_ik ! Roots located successfully
  integer(ik), parameter :: fr_act_stopped  = 2304_ik ! Search direction vanished. Roots may have been located successfully
  integer(ik), parameter :: fr_act_failed   = 2305_ik ! Roots optimization failed to converge
  !
  !  State-machine state codes (see fr_state%mstate below)
  !
  integer(ik), parameter :: fr_mc_waiting   = 3310_ik ! Waiting to perform an optimization step
  integer(ik), parameter :: fr_mc_jacobian  = 3320_ik ! Evaluating the Jacobian
  integer(ik), parameter :: fr_mc_linear    = 3330_ik ! Performing linear search
  integer(ik), parameter :: fr_mc_invalid   = 3930_ik ! State is not initialized
  !
  !  Additional linear-search state machine code (see linear_search_step() below)
  !
  integer(ik), parameter :: fr_ls_direction = 6725_ik ! Have search direction; haven't stepped yet
  integer(ik), parameter :: fr_ls_2point    = 6735_ik ! Done Newton step; ready for 2-point intepolation
  integer(ik), parameter :: fr_ls_3point    = 6745_ik ! Backtracked already; ready for 3-point intepolation
  integer(ik), parameter :: fr_ls_invalid   = 6755_ik ! Linear search is not initialized
  !
  type fr_state
    !
    !  Values which are of interest for the caller
    !
    integer(ik)       :: action          ! Action required from the caller. Can be one of:
                                         ! fr_act_none, fr_act_evaluate, fr_act_done, fr_act_stopped, or fr_act_failed.
    real(rk), pointer :: vars(:)         ! fr_act_evaluate: Values of the parameters
                                         ! fr_act_done: Final roots of the equations
                                         ! fr_act_stopped: Search direction vanished identically. The position
                                         !                 may or may not correspond to a true root.
                                         ! fr_act_failed: The best guess at the roots so far
    real(rk), pointer :: eqs(:)          ! fr_act_evaluate: Function values at vars(:)
    !
    !  Internal state; should not be accessed from outside
    !
    integer(ik)       :: verbose         ! Verbosity level
    integer(ik)       :: nvars           ! Number of variables
    integer(ik)       :: neqs            ! Number of equations to solve
    integer(ik)       :: maxiter         ! Maximum number of iterations allowed
    real(rk)          :: epsstep         ! Tolerance on the displacement
    real(rk)          :: diffstep        ! Displacement to use for numerical differentiation
    real(rk)          :: maxstep         ! Largest allowed step
    !                                    
    integer(ik)       :: mstate          ! Current state of the state machine
    integer(ik)       :: nr_iter         ! Current Newton iteration
    real(rk), pointer :: vars_last(:)    ! Variables from the last optimization step
    integer(ik)       :: jcb_index       ! Current energy evaluation step for Jacobian
    real(rk), pointer :: jcb_values(:,:) ! (neqs,0:2*nvars)
                                         ! First index: function values
                                         ! Second index:
                                         !   0     = function values at vars_last(:)
                                         !   2*i-1 = variable (i) is incremented by diffstep
                                         !   2*i   = variable (i) is decremented by diffstep
    real(rk), pointer :: search(:)       ! Current search direction
    integer(ik)       :: lsstate         ! Mini-state-machine state variable for the linear search;
                                         ! see linear_search_step() below for the specific values
    !
    ! Variables driving linear search; see linear_search_step() for details
    !
    real(rk)          :: g0              ! (half of) 2-norm at zero displacement
    real(rk)          :: gp0             ! Gradient of the (half of) 2-norm at zero displacement
    real(rk)          :: l1              ! Latest displacement
    real(rk)          :: g1              ! (half of) 2-norm at l1 displacement
    real(rk)          :: l2              ! Older displacement
    real(rk)          :: g2              ! (half of) 2-norm at l2 displacement
  end type fr_state
  !
  contains
  !
  subroutine fr_initialize(verbose,state,neqs,vars,maxiter,epsstep,diffstep,maxstep)
    type(fr_state), intent(inout)      :: state     ! Solver state
    integer(ik), intent(in)            :: verbose   ! Desired verbosity level
    integer(ik), intent(in)            :: neqs      ! Number of equations to solve
    real(rk), intent(in)               :: vars(:)   ! Initial values of the parameters
    integer(ik), intent(in), optional  :: maxiter   ! Maximum number of iterations allowed
    real(rk), intent(in), optional     :: epsstep   ! Tolerance on the displacement
    real(rk), intent(in), optional     :: diffstep  ! Displacement to use for numerical differentiation
    real(rk), intent(in), optional     :: maxstep   ! Largest allowed displacement in a single step
    !
    integer(ik) :: alloc
    !
    state%verbose  = verbose
    state%nvars    = size(vars)
    state%neqs     = neqs 
    !
    !  Some defaults
    !
    state%maxiter  = 100
    state%epsstep  = spacing(100._rk)
    state%diffstep = spacing(100._rk)
    state%maxstep  = 0.1_rk * maxval(abs(vars))
    !
    if (present(maxiter )) state%maxiter  = maxiter
    if (present(epsstep )) state%epsstep  = max(epsstep,state%epsstep)
    if (present(diffstep)) state%diffstep = max(diffstep,state%diffstep)
    if (present(maxstep )) state%maxstep  = maxstep
    ! 
    allocate (state%vars(state%nvars),state%eqs(state%neqs), &
              state%vars_last(state%nvars),state%jcb_values(state%neqs,0:2*state%nvars), &
              state%search(state%nvars), &
              stat=alloc)
    if (alloc/=0) then
      write (out,"('find_root%fr_initialize: allocation failed with code ',i0)") alloc
      stop 'find_root%fr_initialize - no memory'
    end if
    !
    state%mstate       = fr_mc_waiting
    state%nr_iter = 1
    state%action       = fr_act_none
    state%vars         = vars
    state%lsstate      = fr_ls_invalid
  end subroutine fr_initialize
  !
  subroutine fr_finalize(state)
    type(fr_state), intent(inout) :: state

    deallocate (state%vars,state%eqs,state%vars_last,state%jcb_values,state%search)
    nullify (state%vars,state%eqs,state%vars_last,state%jcb_values,state%search)
  end subroutine fr_finalize
  !
  subroutine fr_step(state)
    type(fr_state), intent(inout) :: state
    !
    integer(ik) :: ivar
    !
    !  Run the state machine ...
    !
    select case (state%mstate)
      case default
        write (out,"('find_root%fr_step: Unexpected state ',i0)") state%mstate
        stop 'find_root%fr_step - state machine corrupted'
      case (fr_mc_waiting)
        !
        !  Start each step by evaluating the Jacobian
        ! 
        state%mstate    = fr_mc_jacobian
        state%jcb_index = 0
        state%action    = fr_act_evaluate
        state%vars_last = state%vars
        if (state%verbose>2) then 
          write (out,"(/'Evaluating the Jacobian')")
        end if
      case (fr_mc_jacobian)
        !
        !  We received another function value back from the caller; remember it
        !
        state%jcb_values(:,state%jcb_index) = state%eqs(:)
        !
        !  See whether we need another value
        !
        state%jcb_index = state%jcb_index + 1
        if (state%jcb_index<=2*state%nvars) then
          !
          !  Yes, more values are needed. Calculate the displacement, and
          !  ask user for more.
          !
          state%vars = state%vars_last
          ivar = (state%jcb_index+1)/2
          if (state%jcb_index/=2*ivar) then
            state%vars(ivar) = state%vars(ivar) + state%diffstep
          else
            state%vars(ivar) = state%vars(ivar) - state%diffstep
          end if
        else
          !
          !  No, everything we need to build the Jacobian is already here
          !
          call build_search_direction(state)
          !
          !  Occationally, it may happen that the search direction is exactly
          !  zero, so that no firther improvement is possible. When this
          !  happens, we have to declare convergence.
          !
          if (all(state%search==0._rk)) then
            if (state%verbose>=1) then
              write (out,"('Search direction vanished after ',i0,' iterations. Optimization complete.')") state%nr_iter
              state%mstate = fr_mc_invalid
              state%action = fr_act_stopped
              state%vars   = state%vars_last
            end if
            state%mstate = fr_mc_invalid
            state%action = fr_act_done
          else
            state%mstate   = fr_mc_linear
            state%lsstate = fr_ls_direction
            call linear_search_step(state)
          end if
        end if
      case (fr_mc_linear)
        call linear_search_step(state)
        if (state%mstate==fr_mc_waiting) then
          !
          !  Linear search exited. We need to decide whether we reached convergence,
          !  failed, or need to do more iterations.
          !
          if (all(abs(state%vars-state%vars_last)<=state%epsstep)) then
            if (state%verbose>=1) then
              write (out,"('Optimization converged after ',i0,' iterations. Final residue = ',g20.8e3)") state%nr_iter, state%g1
            end if
            state%mstate = fr_mc_invalid
            state%action = fr_act_done
          else
            state%nr_iter = state%nr_iter + 1
            if (state%nr_iter>state%maxiter) then
              if (state%verbose>=1) then
              write (out,"('Optimization failed after ',i0,' iterations. Final residue = ',g20.8e3)") state%nr_iter, state%g1
              end if
              state%mstate = fr_mc_invalid
              state%action = fr_act_failed
            end if
          end if
        end if
    end select
  end subroutine fr_step
  !
  subroutine build_search_direction(state)
    type(fr_state), intent(inout) :: state
    !
    integer(ik) :: ivar
    real(rk)    :: jcb(state%neqs,state%nvars)
    real(rk)    :: gradf(state%nvars)
    real(rk)    :: rhs(max(state%neqs,state%nvars),1)
    !
    build_jacobian: do ivar=1,state%nvars
      jcb(:,ivar) = (state%jcb_values(:,2*ivar-1)-state%jcb_values(:,2*ivar))/(2*state%diffstep)
    end do build_jacobian
    !
    !  Half of the squared norm; we'll need to to perform linear search later
    !
    state%g0 = 0.5_rk * sum(state%jcb_values(:,0)**2)
    !
    !  Gradient of the half of squared norm; we'll need it to calculate gradient
    !  along the search direction
    !
    gradf = matmul(state%jcb_values(:,0),jcb)
    !
    rhs(:state%neqs,1) = -state%jcb_values(:,0)
    call lapack_gelss(jcb,rhs)
    state%search(:) = rhs(:state%nvars,1)
    state%gp0       = sum(gradf*state%search)
    if (state%verbose>2) then
      write (out,"('Gradient of (1/2) squared norm:')")
      write (out,"(4(1x,g32.20e3))") gradf
      write (out,"('Search direction vector:')")
      write (out,"(4(1x,g32.20e3))") state%search
      write (out,"('Gradient of (1/2) squared norm along search direction = ',g32.20e3)") state%gp0
    end if
  end subroutine build_search_direction
  !
  subroutine linear_search_step(state)
    type(fr_state), intent(inout) :: state
    !
    linsearch: select case (state%lsstate)
      case default
        write (out,"('find_root%linear_search_step: state ',i0,' is not recognized')") state%lsstate
        stop 'find_root%linear_search_step - state machine corrupted'
      case (fr_ls_direction)
        state%l1   = 1._rk
        call limit_step_size(state%l1)
        if (state%verbose>2) then
          write (out,"('Initial step: ',g24.12,' x search direction')") state%l1
        end if
        state%lsstate = fr_ls_2point 
        state%vars    = state%vars_last + state%l1 * state%search
      case (fr_ls_2point)
        if (.not.step_accepted()) then
          !
          !  Step rejected; perform 2-point backtrack
          !
          call backtrack_2point
          if (state%verbose>2) then
            write (out,"('2-point backtrack: ',g24.12,' x search direction')") state%l1
          end if
          state%lsstate = fr_ls_3point 
          state%vars    = state%vars_last + state%l1 * state%search
          call limit_backtrack
        end if
      case (fr_ls_3point)
        if (.not.step_accepted()) then
          !
          !  Step rejected; perform 3-point backtrack
          !
          call backtrack_3point
          if (state%verbose>2) then
            write (out,"('3-point backtrack: ',g24.12,' x search direction')") state%l1
          end if
          state%vars = state%vars_last + state%l1 * state%search
          call limit_backtrack
        end if
    end select linsearch
    !
    contains
    !
    !  Makes sure step size does not exceed the allowable maximum
    !
    subroutine limit_step_size(step)
      real(rk), intent(inout) :: step
      !
      real(rk) :: cs ! Cartesian step size
      !
      cs = sqrt(sum(state%search**2))*abs(step)
      if (cs>state%maxstep) step = step * (state%maxstep/cs)
    end subroutine limit_step_size
    !
    !  If the proposed step is too small, we may need to terminate the search
    !  and revert to our starting position. In order to abort the linear search,
    !  we'll use a tigher convergence criterion here.
    ! 
    subroutine limit_backtrack
      if (all(abs(state%l1 * state%search)<=1e-4_rk*state%epsstep)) then
        if (state%verbose>=1) then
          write (out,"('Unable to improve the solution along linear search direction.')")
        end if
        state%vars    = state%vars_last
        state%g1      = state%g0
        state%mstate  = fr_mc_waiting
        state%lsstate = fr_ls_invalid
      end if
    end subroutine limit_backtrack
    !
    !  Verify whether current step brings us sufficiently closer
    !  to the solution.
    !
    logical function step_accepted ()
      real(rk) :: target, extrap
      !
      state%g1 = 0.5_rk * sum(state%eqs**2)
      extrap   = state%g0 +  0.5_rk * state%gp0 * state%l1
      target   = state%g0 + 1e-4_rk * state%gp0 * state%l1
      if (state%g1<target) then
        if (state%verbose>2) then
          write (out,"('Accepted step ',g24.12,' residue ',g24.12e3,' extrap. ',g24.12e3)") &
                 state%l1, state%g1, extrap
        end if
        !
        !  Terminate linear search; the final displacement is already in state%vars
        !
        state%mstate  = fr_mc_waiting
        state%lsstate = fr_ls_invalid
        step_accepted = .true.
      else
        if (state%verbose>2) then
          write (out,"('Rejected step ',g24.12,' residue ',g24.12e3,' extrap. ',g24.12e3,' threshold ',g24.12e3)") &
                 state%l1, state%g1, extrap, target
        end if
        step_accepted = .false.
      end if
    end function step_accepted
    !
    subroutine backtrack_2point
      real(rk) :: num, den
      !
      state%l2 = state%l1 
      state%g2 = state%g1
      ! Try to prevent an overflow at all costs
      num      = -0.5_rk * state%l2 * state%gp0
      den      = state%g2 - state%g0 - state%gp0
      if (abs(den)>abs(num/sqrt(huge(num)))) state%l1 = num/den
      !
      !  Limit the backtrack to between 30% and 90% of the initial step
      !
      state%l1 = min(max(0.3_rk*state%l2,state%l1),0.9_rk*state%l2)
    end subroutine backtrack_2point
    !
    !  3-point backtrack is probably an overkill, but what the heck
    !  See Numerical Recipes eqs. 9.7.13 and 9.7.14
    !
    subroutine backtrack_3point
      real(rk) :: r1, r2
      real(rk) :: a, b, lambda, det
      !
      r1  = state%g1 - state%gp0*state%l1 - state%g0
      r2  = state%g2 - state%gp0*state%l2 - state%g0
      ! At this point, we know that l1/=l2, so we are safe
      a   = ( r1/state%l1**2 - r2/state%l2**2                   ) / (state%l1-state%l2)
      b   = (-r1*state%l2/state%l1**2 + r2*state%l1/state%l2**2 ) / (state%l1-state%l2)
      det = b**2-3*a*state%gp0
      det = -b + sqrt(abs(det))
      if (abs(det)<1.5_rk*state%l1*abs(a)) then
        lambda = det/(3*a)
      else
        lambda = 1._rk
      end if
      lambda = min(0.5_rk*state%l1,lambda)
      lambda = max(0.1_rk*state%l1,lambda)
      !
      state%l2 = state%l1 
      state%g2 = state%g1
      state%l1 = lambda
    end subroutine backtrack_3point
  end subroutine linear_search_step
  !
end module find_root

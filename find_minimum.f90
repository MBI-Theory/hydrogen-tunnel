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
!  Simple conjugate-gradient minimization routine. The implementation follows 
!  discussion in section 10.6 of Numerical Recipes. Gradients are calculated
!  by numerical differentiation of the target function.
!
!  After fm_initialize() was called, all interaction with the calling
!  routine occurs through the fm_state data structure.
!
!  All routines here are thread-safe.
!
module find_minimum
  use accuracy
  implicit none
  private
  public fm_state 
  public fm_initialize, fm_step, fm_finalize
  public fm_act_evaluate, fm_act_gradient, fm_act_done, fm_act_failed
  public rcsid_find_minimum
  !
  character(len=clen), save :: rcsid_find_minimum = "$Id: $"
  !
  !  User activity codes (see fm_state%action below)
  !
  integer(ik), parameter :: fm_act_none     = 4300_ik ! No action required; call again
  integer(ik), parameter :: fm_act_evaluate = 4302_ik ! Evaluate function for the given value of vars()
  integer(ik), parameter :: fm_act_gradient = 4303_ik ! Evaluate function gradient for the given value of vars()
                                                      ! Will only be called if fm_initialize() indicates that
                                                      ! gradients are supported.
  integer(ik), parameter :: fm_act_done     = 4304_ik ! Minimum located successfully
  integer(ik), parameter :: fm_act_failed   = 4305_ik ! Optimization failed to converge
  !
  !  State-machine state codes (see fm_state%mstate below)
  !
  integer(ik), parameter :: fm_mc_waiting   = 5310_ik ! Waiting to start optimization
  integer(ik), parameter :: fm_mc_p1        = 5320_ik ! Function evaluation; initialization
  integer(ik), parameter :: fm_mc_p2        = 5330_ik ! Gradient evaluation; initialization
  integer(ik), parameter :: fm_mc_p3        = 5340_ik ! Linear search
  integer(ik), parameter :: fm_mc_p4        = 5350_ik ! Gradient evaluation; in-function

  integer(ik), parameter :: fm_mc_invalid   = 5500_ik ! Optimization exited
  !
  !  Linear-search state codes (see fm_state%lstate below)
  !
  integer(ik), parameter :: ls_mc_waiting   = 6310_ik ! Waiting to start linear search
  integer(ik), parameter :: ls_mc_brak_p1   = 6320_ik ! Braketing, first function evaluation
  integer(ik), parameter :: ls_mc_brak_p2   = 6330_ik ! Braketing, second function evaluation
  integer(ik), parameter :: ls_mc_brak_p3   = 6340_ik ! Braketing, third function evaluation
  integer(ik), parameter :: ls_mc_min_p1    = 6410_ik ! Minimum search on the braketed interval
  integer(ik), parameter :: ls_mc_done      = 6500_ik ! Linear search complete
  !
  type fm_state
    !
    !  Values which are of interest for the caller
    !
    integer(ik)           :: action          ! Action required from the caller. Can be one of:
                                             ! fm_act_none, fm_act_evaluate, fm_act_gradient 
                                             ! (only if havegrad is .true.), fm_act_done, or 
                                             ! fm_act_failed.
    real(rk), pointer     :: vars(:)         ! fm_act_evaluate: Values of the parameters
                                             ! fm_act_done: Final minimum point 
                                             ! fm_act_failed: The best guess at the roots so far
    real(rk)              :: func            ! fm_act_evaluate: Function value at vars(:)
                                             ! fm_act_done: Function at the final minimum point
    real(rk), allocatable :: grad(:)         ! fm_act_gradient: Function gradient at vars(:)
    !
    !  Internal state; should not be accessed from outside
    !
    integer(ik)           :: verbose         ! Verbosity level
    integer(ik)           :: nvars           ! Number of variables
    integer(ik)           :: maxiter         ! Maximum number of iterations allowed
    real(rk)              :: epsstep         ! Tolerance on the step.
    real(rk)              :: diffstep        ! Displacement to use for numerical differentiation
    real(rk)              :: maxstep         ! Largest allowed step
    logical               :: havegrad        ! True if function gradients are supported
    !                                        
    integer(ik)           :: mstate          ! Current state of the state machine
    integer(ik)           :: iter            ! Current congugate gradient step
    integer(ik)           :: grad_index      ! Current gradient evaluation step.
                                             ! 2*i-1 = increment variable (i)
                                             ! 2*i   = decrement variable (i)
    real(rk)              :: grad_temp       ! Temporary used by gradient evaluation
    !
    real(rk), allocatable :: vars_save(:)    ! Temporary, used by gradient and linear search
    !
    integer(ik)           :: converged_steps ! Number of consequitive linear searches satisfying 
                                             ! the displacement tolerance
    !
    !  Variables used in the linear search
    !
    integer(ik)           :: lstate          ! Current state of the linear search
    real(rk)              :: step_sz         ! Initial step size used for braketing
                                             ! We'll begin with the maximum step size, and then
                                             ! update after each linear minimization step
    !
    !  Variables used to braket the minimum and do the linear search
    !
    real(rk)              :: ax, bx, cx, ux  ! Points braketing the minimum. We requires ax<bx<cx
    real(rk)              :: fa, fb, fc, fu  ! Function values at the braketing points
    real(rk)              :: sx              ! The sign of the displacement; either +1 or -1
    !
    !  Variables used in the CG algorithm
    !
    real(rk)              :: fret
    real(rk)              :: dgg             ! 
    real(rk)              :: fp              !
    real(rk)              :: gam             !
    real(rk)              :: gg              !
    real(rk), allocatable :: g (:)           ! 
    real(rk), allocatable :: xi(:)           ! Used to evaluate function derivative
    real(rk), allocatable :: h (:)           ! 
  end type fm_state
  !
  contains
  !
  subroutine fm_initialize(verbose,state,vars,maxiter,epsstep,diffstep,maxstep,havegrad)
    type(fm_state), intent(inout)      :: state     ! Solver state
    integer(ik), intent(in)            :: verbose   ! Desired verbosity level
    real(rk), intent(in)               :: vars(:)   ! Initial values of the parameters
    integer(ik), intent(in), optional  :: maxiter   ! Maximum number of iterations allowed
    real(rk), intent(in), optional     :: epsstep   ! Tolerance on the displacement
    real(rk), intent(in), optional     :: diffstep  ! Displacement to use for numerical differentiation
    real(rk), intent(in), optional     :: maxstep   ! Largest allowed displacement in a single step
    logical, intent(in), optional      :: havegrad ! Caller supports function gradients
                                                   ! The default is .false.
    !
    integer(ik) :: alloc, nv
    !
    state%verbose  = verbose
    state%nvars    = size(vars)
    !
    !  Some defaults
    !
    state%maxiter  = 100
    state%epsstep  = spacing(100._rk)
    state%diffstep = spacing(100._rk)
    state%maxstep  = 0.1_rk * maxval(abs(vars))
    state%havegrad = .false.
    !
    if (present(maxiter )) state%maxiter  = maxiter
    if (present(epsstep )) state%epsstep  = max(epsstep,state%epsstep)
    if (present(diffstep)) state%diffstep = max(diffstep,state%diffstep)
    if (present(maxstep )) state%maxstep  = maxstep
    if (present(havegrad)) state%havegrad = havegrad
    ! 
    nv = state%nvars
    allocate (state%vars(nv),state%grad(nv),state%vars_save(nv),state%g(nv),state%xi(nv),state%h(nv),stat=alloc)
    if (alloc/=0) then
      write (out,"('find_minimum%fm_initialize: allocation failed with code ',i0)") alloc
      stop 'find_minimum%fm_initialize - no memory'
    end if
    !
    state%converged_steps = 0
    state%mstate  = fm_mc_waiting
    state%action  = fm_act_none
    state%vars    = vars
  end subroutine fm_initialize
  !
  subroutine fm_finalize(state)
    type(fm_state), intent(inout) :: state

    deallocate (state%vars,state%grad,state%vars_save,state%g,state%h,state%xi)
  end subroutine fm_finalize
  !
  subroutine fm_step(state)
    type(fm_state), intent(inout) :: state
    !
    real(rk) :: disp
    !
    !  The state machine below is a literal translation of NR frprmn() routine
    !  It looks a little ugly since we need to simulate function calls 
    !
    state_machine: select case (state%mstate)
      case default
        write (out,"('find_minimum%fm_step: Unexpected state ',i0)") state%mstate
        stop 'find_minimum%fm_step - state machine corrupted'
      case (fm_mc_waiting)
        state%action = fm_act_evaluate
        state%mstate = fm_mc_p1         ! Evaluate fp - function value at the start of the optimization
      case (fm_mc_p1)
        state%fp         = state%func
        state%grad_index = 0
        call update_gradient(state)
        state%mstate = fm_mc_p2         ! Evaluate gradient; this requires a _lot_ of callbacks
                                        ! if analytical gradients are not available
      case (fm_mc_p2)
        call update_gradient(state)
        if (state%action/=fm_act_none) exit state_machine
        ! Gradient is now in xi
        if (sum(state%xi**2)==0._rk) then 
          ! Function gradient vanishes at the initial point. Exit optimization
          state%func   = state%fp
          state%action = fm_act_done
          state%mstate = fm_mc_invalid
          exit state_machine
        end if
        state%g      = -state%xi
        state%h      =  state%g
        state%xi     =  state%h
        ! Ready to perform the first linear search
        state%iter    = 1
        state%lstate  = ls_mc_waiting
        state%step_sz = state%maxstep
        call linear_search(state)
        state%mstate = fm_mc_p3         ! Perform linear search; this requires a _lot_ of callbacks
      case (fm_mc_p3)
        call linear_search(state)
        if (state%action/=fm_act_none) exit state_machine
        ! Linear search complete; new variables are in vars(:). Corresponding function value is in fret
        ! linear_search() also leaves the values of variables before the search in vars_save()
        disp = sqrt(sum((state%vars-state%vars_save)**2))
        if (state%verbose>=1) then
          write (out,"('On the linear search step ',i0,' displacement is ',g0.12,' function = ',g0.24)") &
                 state%iter, disp, state%fret
        end if
        if (disp<state%epsstep) then
          state%converged_steps = state%converged_steps + 1
          if (state%converged_steps>=min(3_ik,state%nvars)) then
            state%func   = state%fret
            state%action = fm_act_done
            state%mstate = fm_mc_invalid
            exit state_machine
          end if
        else
          state%converged_steps = 0
        end if
        state%fp = state%fret
        state%grad_index = 0
        call update_gradient(state)
        state%mstate = fm_mc_p4     ! Evaluate gradient; this potentially requires a _lot_ of callbacks
      case (fm_mc_p4)
        call update_gradient(state)
        if (state%action/=fm_act_none) exit state_machine
        state%gg  = sum(state%g**2)
        state%dgg = sum((state%xi+state%g)*state%xi)
        if (state%gg==0._rk) then 
          ! Function gradient vanishes. Exit optimization
          state%func   = state%fret
          state%action = fm_act_done
          state%mstate = fm_mc_invalid
          exit state_machine
        end if
        ! Calculate the new search direction
        state%gam = state%dgg / state%gg
        state%g   = -state%xi
        state%h   = state%g + state%gam*state%h
        state%xi  = state%h
        !
        if (state%verbose>=2) then
          write (out,"(/'On the optimization step ',i0,' gradient is:')") state%iter
          write (out,"((5(1x,g24.12)))") -state%g
          write (out,"('Search direction:')")
          write (out,"((5(1x,g24.12)))")  state%xi
          write (out,"()")
        end if
        !
        if (sum(state%xi**2)==0._rk) then 
          ! Search direction vanished. We are done
          state%func   = state%fret
          state%action = fm_act_done
          state%mstate = fm_mc_invalid
          exit state_machine
        end if
        !
        state%iter = state%iter + 1
        if (state%iter>state%maxiter) then
          ! Too many iterations. Optimization failed
          state%func   = state%fret
          state%action = fm_act_failed
          state%mstate = fm_mc_invalid
          exit state_machine
        end if
        !
        state%lstate  = ls_mc_waiting
        call linear_search(state)
        state%mstate = fm_mc_p3         ! Perform linear search; this requires a _lot_ of callbacks
    end select state_machine
  end subroutine fm_step
  !
  !  update_gradient() evaluates target function gradient by symmetric numerical
  !  differentiation. On the first call, state%grad_index must be set to zero.
  !  Once gradient evaluation is complete, state%action will be set to fm_act_none,
  !  and the gradient will be in state%xi
  !
  subroutine update_gradient(state)
    type(fm_state), intent(inout) :: state
    !
    integer(ik) :: ivar
    !
    if (state%havegrad) then
      !
      !  Analytical gradients are available; evaluate them
      !
      if (state%grad_index==0) then
        ! Ask for the gradient to be placed in state%grad(:)
        state%action = fm_act_gradient
        state%grad_index = 1
      else
        ! We are done
        state%action  = fm_act_none
        state%xi(:)   = state%grad(:)
      end if
    else ! state%havegrad==.false.
      !
      !  Analytical gradients are not available; evalate gradient using finite differences
      !
     if (state%grad_index>0) then
       !
       !  We come in with a function value available
       !
       ivar = (state%grad_index+1)/2
       if (mod(state%grad_index,2)==0) then
         !
         !  Second half-step. Update the derivative (in state%xi)
         !
         state%xi(ivar) = (state%grad_temp - state%func) / (2._rk*state%diffstep)
       else
         !
         !  First half-step. Remember function value.
         !
         state%grad_temp = state%func
       end if
     else
       !
       !  We are about to state the displacements - remember the variables
       !
       state%vars_save(:) = state%vars(:)
     end if
     !
     !  Advance to the next function value
     !
     state%grad_index = state%grad_index + 1
     if (state%grad_index>2*state%nvars) then
       !
       !  Gradient evaluation complete, restore variables and return
       !
       state%vars(:) = state%vars_save(:)
       state%action  = fm_act_none
     else
       ivar = (state%grad_index+1)/2
       if (mod(state%grad_index,2)==0) then
         !
         !  Decrement variable
         !
         state%vars(:)    = state%vars_save(:)
         state%vars(ivar) = state%vars(ivar) - state%diffstep
       else
         !
         !  Increment variable
         !
         state%vars(:)    = state%vars_save(:)
         state%vars(ivar) = state%vars(ivar) + state%diffstep
       end if
       state%action = fm_act_evaluate
     end if
    end if
  end subroutine update_gradient
  !
  !  Linear search along direction in state%xi, starting from state%vars
  !  At the first call, state%lstate must be ls_mc_waiting
  !  At the end of the linear search, state%vars will be updated to the new
  !  position of the minimum. state%fret will be function value at these
  !  coordinates. 
  !  The end of the linear search is reached when state%action is fm_act_none
  !
  subroutine linear_search(state)
    type(fm_state), intent(inout) :: state
    !
    real(rk), parameter :: golden_ratio = 0.5_rk*(sqrt(5._rk)+1._rk)
    real(rk)            :: ca ! Length of the braketing interval
    !
    state_machine: select case (state%lstate)
      case default
        write (out,"('find_minimum%linear_search: Unexpected state ',i0)") state%lstate
        stop 'find_minimum%find_minimum - state machine corrupted'
      case (ls_mc_waiting)
        !
        if (state%verbose>=3) then 
          write (out,"(/'Linear search direction:')")
          write (out,"((5(1x,g24.14)))") state%xi
        end if
        !
        state%action    = fm_act_evaluate
        !
        !  We'll begin by bracketing the minimum, using the golden-ratio interval expansion
        !
        state%sx        = 1._rk
        state%ax        = 0._rk
        state%bx        = state%step_sz/(1._rk+golden_ratio) + spacing(1._rk)
        state%vars_save = state%vars
        call displace(state%ax)
        state%lstate = ls_mc_brak_p1
      case (ls_mc_brak_p1)
        state%fa = state%func
        call displace(state%bx)
        state%lstate = ls_mc_brak_p2
      case (ls_mc_brak_p2)
        state%fb = state%func
        if (state%fb>state%fa) then
          call swap(state%fa,state%fb)
          call swap(state%ax,state%bx)
          state%ax = -state%ax
          state%bx = -state%bx
          state%sx = -1._rk
        end if
        ! We would like to maintain ax<bx<cx; this makes the code less messy
        state%cx = state%bx + golden_ratio*(state%bx-state%ax)
        call displace(state%cx)
        state%lstate = ls_mc_brak_p3
      case (ls_mc_brak_p3)
        ! At this point, we know that fa>=fb
        state%fc = state%func
        !
        if (state%verbose>=3) then
          write (out,"(/'Bracketing the minimum')")
          write (out,"('s, a, b, c = ',f4.0,1x,3g24.12)") state%sx, state%ax, state%bx, state%cx
          write (out,"('  fa,fb,fc = ',4x,  1x,3g24.12)")           state%fa, state%fb, state%fc
        end if
        !
        if (state%fb<state%fc) then
          ! We are done. The minimum is between ax and cx
          ! Form the guess for the minimum based on quadratic interpolation between ax, bx, and cx
          ! The guess position is in ux
          state%ux = quadratic_minimum(state%ax,state%bx,state%cx,state%fa,state%fb,state%fc,state%epsstep)
          call displace(state%ux)
          state%lstate = ls_mc_min_p1
          exit state_machine
        end if
        ! No braket. Discard point a, and continue searching
        state%ax = state%bx ; state%fa = state%fb
        state%bx = state%cx ; state%fb = state%fc
        state%cx = state%bx + golden_ratio*(state%bx-state%ax)
        call displace(state%cx)
        state%lstate = ls_mc_brak_p3
      case (ls_mc_min_p1)
        state%fu = state%func
        !
        if (state%verbose>=3) then
          write (out,"(/'Locating the minimum')")
          write (out,"('s, a, b, c, u = ',f3.0,1x,4g24.12)") state%sx, state%ax, state%bx, state%cx, state%ux
          write (out,"('  fa,fb,fc,fu = ',3x,  1x,4g24.12)")           state%fa, state%fb, state%fc, state%fu
        end if
        ! We are guaranteed to have the points in the order: a .. b .. c
        ! u is guaranteed to be between a and c, and does not coincide with any of a, b, or c
        if (state%ux<state%bx) then
          ! a .. u .. b .. c
          if (state%fu<state%fb) then
            ! The new interval: a .. u .. b
            state%cx = state%bx ; state%fc = state%fb
            state%bx = state%ux ; state%fb = state%fu
          else
            ! The new interval: u .. b .. c
            state%ax = state%ux ; state%fa = state%fu
          end if
        else
          ! a .. b .. u .. c
          if (state%fu<state%fb) then
            ! The new interval: b .. u .. c
            state%ax = state%bx ; state%fa = state%fb
            state%bx = state%ux ; state%fb = state%fu
          else
            ! The new interval: a .. b .. u
            state%cx = state%ux ; state%fc = state%fu
          end if
        end if
        !
        !  Are we converged?
        !
        ca = state%cx - state%ax
        if (ca <= 0.5_rk * state%epsstep) then
          call displace(state%bx)
          state%fret    = state%fb
          state%step_sz = abs(state%bx)
          state%action  = fm_act_none
          state%lstate  = ls_mc_done
        end if
        !
        !  Not yet converged; continue linear search
        !
        state%ux = quadratic_minimum(state%ax,state%bx,state%cx,state%fa,state%fb,state%fc,state%epsstep)
        call displace(state%ux)
    end select state_machine
    !
    contains
    subroutine displace(dx)
      real(rk), intent(in) :: dx
      real(rk)             :: lxi
      !
      lxi = sqrt(sum(state%xi**2))
      if (lxi<=0._rk) stop 'find_minimum%displace - Search direction vanished'
      state%vars = state%vars_save + state%xi * state%sx * dx / lxi
    end subroutine displace
    subroutine swap(a,b)
      real(rk), intent(inout) :: a, b
      real(rk)                :: t
      t = a ; a = b ; b = t
    end subroutine swap
  end subroutine linear_search
  !
  function quadratic_minimum(ax,bx,cx,fa,fb,fc,eps) result(ux)
    real(rk), intent(in) :: ax, bx, cx ! Points to fit to; bx should be between ax and cx. 
                                       ! The minimum must be braketed by ax and cx.
    real(rk), intent(in) :: fa, fb, fc ! Function values at ax, bx, and cx
    real(rk), intent(in) :: eps
    real(rk)             :: ux         ! Position where quadratic minimum should be.
    !
    real(rk) :: da, dc, dfa, dfc       ! Displacements relative to (bx,fb)
    real(rk) :: a, b                   ! Quadratic and linear coefficients around (b)
                                       ! We'll skip the (immaterial, common) denominator
    !
    !  Some sanity checking could not hurt!
    !
    if (ax>=bx .or. cx<=bx) then
      write (out,"('Braketing points are not ordered: ',3(g34.24,1x))") ax, bx, cx
      stop 'find_minimum%quadratic_minimum - ordering failure'
    end if
    if (fa<fb .or. fc<fb) then
      write (out,"('Braketing points do not contain a minimum: ',3(g34.24,1x))") fa, fb, fc
      stop 'find_minimum%quadratic_minimum - braketing failure'
    end if
    da  = ax - bx
    dc  = cx - bx
    dfa = fa - fb
    dfc = fc - fb
    a   = dfc*da    - dfa*dc    ! We drop the common denominator: da*dc*(da-dc)
    b   = dfc*da**2 - dfa*dc**2
    !
    ! (a) should not normally vanish, since da and dc have opposite signs, while dfc and dfa have the same sign
    !     However, if the function is very flat, we could have a situation where both dfa and dfc vanish.
    !     In this case, we must pick ux at the middle of the interval.
    !
    if (a==0._rk) then
      ux = bx
    else
      ux = bx + 0.5_rk*b/a
    end if
    !
    !  If ux coincides with any of the points, or falls outside the interval, 
    !  arbitrarily displace it.
    !
    if (ux<=ax .or. ux>=cx) then
      if (cx-bx>bx-ax) then
        ux = 0.5_rk*(bx+cx)
      else
        ux = 0.5_rk*(ax+bx)
      end if
    end if
    if (abs(ux-bx)<=eps) then
      if (cx-bx>bx-ax) then
        ux = bx + 0.1_rk*(cx-bx)
      else
        ux = bx - 0.1_rk*(bx-ax)
      end if
    end if
  end function quadratic_minimum
  !
end module find_minimum
!
! program test
!   use accuracy
!   use find_minimum
!   !
!   type (fm_state) :: minimum
!   real(rk)        :: mat(3,3)
!   !
!   mat(1,:) = (/ 1.0_rk, 0.3_rk, 0.4_rk /)
!   mat(2,:) = (/ 0.6_rk, 2.0_rk, 0.5_rk /)
!   mat(3,:) = (/ 0.7_rk, 0.8_rk, 3.0_rk /)
!   !
!   call fm_initialize(4_ik,minimum,vars=(/1._rk,2._rk,3._rk/),havegrad=.true.,diffstep=spacing(1e5_rk),epsstep=spacing(1._rk))
!   !
!   min_iterations: do
!     call fm_step(minimum)
!     select case (minimum%action)
!       case default
!         stop 'test - bad action'
!       case (fm_act_evaluate)
!         call evaluate_function
!       case (fm_act_gradient)
!         call evaluate_gradient
!       case (fm_act_done)
!         write (out,"('Minimization converged.')")
!         exit min_iterations
!       case (fm_act_failed)
!         write (out,"('Minimization failed to converge.')")
!         exit min_iterations
!     end select
!   end do min_iterations
!   call evaluate_function
!   !
!   contains
!   subroutine evaluate_function
!     minimum%func = sum(minimum%vars * matmul(mat,minimum%vars))
!     write (out,"('At vars = ',3g24.16,' func = ',g24.16)") minimum%vars, minimum%func
!   end subroutine evaluate_function
!   subroutine evaluate_gradient
!     minimum%grad = matmul(mat,minimum%vars) + matmul(minimum%vars,mat)
!     write (out,"('At vars = ',3g24.16,' grad = ',3g24.16)") minimum%vars, minimum%grad
!   end subroutine evaluate_gradient
! end program test

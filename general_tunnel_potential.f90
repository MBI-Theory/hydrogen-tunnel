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
!  Central potentials for use with general_tunnel.f90 and friemds
!
module general_tunnel_potential
  use accuracy
  use constants
  use general_tunnel_data
  use math
  use poly_tools
  use timer
  !$ use OMP_LIB
  implicit none
  private
! public potential_u_r ! We should not use u(r) - this is just one special case of the general
!                      ! u(eta,xi) that we can handle in the code!
  public potential_u_halfx2, check_asymptotic_potential
  public cache_potential
  public rcsid_general_tunnel_potential
  !
  character(len=clen), save :: rcsid_general_tunnel_potential = "$Id: $"
  !
  !  Parameter tables from Tong and Lin, J Phys B 38, 2593 (2005), except for Xe, which is
  !  from Zhang, Lan, and Lu, Phys Rev A 90, 043410 (2014).
  !
  !  The potential is in the form: -(1/r)*(Z + a1*exp(-a2*r) + a3*r*exp(-a4*r) + a5*exp(-a6*r) + a7*r**2*exp(-a8*r))
  !
  !  The last term is not actually in Tong&Lin.
  !
  !                                              Zc         a1        a2          a3        a4         a5        a6     a7     a8
  real(rk), parameter :: tong05_he(0:8) = (/ 1.0_rk,  1.231_rk, 0.662_rk,  -1.325_rk, 1.236_rk, -0.231_rk, 0.480_rk, 0._rk, 1._rk /)
  real(rk), parameter :: tong05_ne(0:8) = (/ 1.0_rk,  8.069_rk, 2.148_rk,  -3.570_rk, 1.986_rk,  0.931_rk, 0.602_rk, 0._rk, 1._rk /)
  real(rk), parameter :: tong05_ar(0:8) = (/ 1.0_rk, 16.039_rk, 2.007_rk, -25.543_rk, 4.525_rk,  0.961_rk, 0.443_rk, 0._rk, 1._rk /)
  real(rk), parameter :: tong05_xe(0:8) = (/ 1.0_rk, 51.356_rk, 2.112_rk, -99.927_rk, 3.737_rk,  1.644_rk, 0.431_rk, 0._rk, 0._rk /)
  real(rk), parameter :: test_he  (0:8) = (/ 1.0_rk,  1.000_rk, 2.1325_rk,  0.000_rk, 0.000_rk,  0.000_rk, 0.000_rk, 0._rk, 0._rk /)
  real(rk), parameter :: test_3he (0:8) = (/ 1.0_rk, 1.2407_rk, 1.6527_rk,  0.000_rk, 0.000_rk,  0.000_rk, 0.000_rk, 0._rk, 0._rk /)
  !
  real(rk), parameter :: sfg99_li(0:8)  = (/ 1.0,  2.0, 3.395,  3.212, 3.207, 0., 0., 0., 0. /)
  real(rk), parameter :: sfg99_na(0:8)  = (/ 1.0, 10.0, 7.902, 23.51 , 2.688, 0., 0., 0., 0. /)
  real(rk), parameter :: sfg99_k (0:8)  = (/ 1.0, 18.0, 3.491, 10.591, 1.730, 0., 0., 0., 0. /)
  real(rk), parameter :: sfg99_rb(0:8)  = (/ 1.0, 36.0, 3.431, 10.098, 1.611, 0., 0., 0., 0. /)
  real(rk), parameter :: sfg99_cs(0:8)  = (/ 1.0, 54.0, 3.294, 11.005, 1.509, 0., 0., 0., 0. /)
  !
  contains
  !
  !  Evaluate the quantity:
  !
  !   (1/m!) (d^m/d r^m) r^n exp(-alpha r)
  !
  !  We'll use recursions to get there :)
  !
  subroutine rexptab(nmax,mmax,r,alpha,rtab)
    integer(ik), intent(in) :: nmax  ! Desired order of the polynomial term (must be >= 0)
    integer(ik), intent(in) :: mmax  ! Maximum desired order of the derivative (min is 0)
    real(rk), intent(in)    :: r     ! Point at which we evaluate the derivatives
    real(rk), intent(in)    :: alpha ! The exponent
    real(rk), intent(out)   :: rtab(0:mmax) 
    !
    integer(ik) :: n, m
    real(rk)    :: tab(0:nmax,0:mmax)
    !
    if (nmax<0 .or. mmax<0) stop 'general_tunnel%rexptab - bad arguments'
    !
    if (alpha*r>=log(sqrt(huge(1._rk)))) then
      tab(0,0) = 0._rk
    else
      tab(0,0) = exp(-alpha*r) 
    end if
    !
    power_only: do n=1,nmax
      tab(n,0) = r * tab(n-1,0)
    end do power_only
    !
    order_only: do m=1,mmax
      tab(0,m) = (-alpha/real(m,kind=rk)) * tab(0,m-1)
    end do order_only
    !
    order: do m=1,mmax
      power: do n=1,nmax
        tab(n,m) = (real(n,kind=rk)/real(m,kind=rk)) * tab(n-1,m-1) &
                           - (alpha/real(m,kind=rk)) * tab(n,m-1)
      end do power
    end do order
    rtab(:) = tab(nmax,:)
  end subroutine rexptab
  !
  subroutine tong05(a,order,r0,ptab)
    real(rk), intent(in)    :: a(0:8) ! Parameters of the potential; see comment above
    integer(ik), intent(in) :: order 
    real(rk), intent(in)    :: r0
    real(rk), intent(out)   :: ptab(0:order)
    real(rk)                :: term(0:order)
    !
    !  We can synthesize Tong05 potentials from the generalized Yukawa terms
    !  provided by rexptab above. 
    !
    if (znuc/=a(0)) then
      write (out,"('ERROR: znuc (',g0.12,') does not match long-range part of the potential (',g0.12,')')") znuc, a(0)
      stop 'general_tunnel_potential%tong05 - bad znuc'
    end if
    call rexptab(0_ik,order,r0,a(2),term) ; ptab(:) =           a(1) * term(:)
    call rexptab(1_ik,order,r0,a(4),term) ; ptab(:) = ptab(:) + a(3) * term(:)
    call rexptab(0_ik,order,r0,a(6),term) ; ptab(:) = ptab(:) + a(5) * term(:)
    call rexptab(2_ik,order,r0,a(8),term) ; ptab(:) = ptab(:) + a(7) * term(:)
  end subroutine tong05
  !
  !  Our internal representation of the potential is in terms of (u), defined as:
  !
  !   u(r) = -r*v(r) - znuc
  !
  !  where v(r) is the potential (attractive=negative), and znuc is the long-range
  !  Coulombic part. Thus, a purely Coulombic potential gives u(r)=0. All other
  !  physical potentials we are concerned with give short-range u(r).
  !
  subroutine potential_u_r(r0,order,ptab)
    real(rk), intent(in)    :: r0
    integer(ik), intent(in) :: order         ! Desired derivative order
    real(rk), intent(out)   :: ptab(0:order) ! ptab(n) = n! d^n u(r)/d r^n, so that
                                             ! u(r-r0) = Sum ptab(n) (r-r0)^n
    real(rk) :: lp(0:order)
    !
    select case (potential)
      case default
        write (out,"('Potential ',a,' is not implemented')") trim(potential)
        stop 'general_tunnel%potential_u - unknown potential'
      case ('hydrogenic')
        ptab = 0.0
      case ('yukawa')
        !
        !  v(r) = a0*r**n1*exp(-a1*r)
        !
        ptab(:) = 0
        if (pot_param_real(1)/=0._rk) then
          call rexptab(pot_param_int(1)+1_ik,order,r0,pot_param_real(2),lp)
          ptab(:) = ptab(:) -pot_param_real(1) * lp(:)
        end if
        if (pot_param_real(3)/=0._rk) then
          call rexptab(pot_param_int(2)+1_ik,order,r0,pot_param_real(4),lp)
          ptab(:) = ptab(:) -pot_param_real(3) * lp(:)
        end if
        if (pot_param_real(5)/=0._rk) then
          call rexptab(pot_param_int(3)+1_ik,order,r0,pot_param_real(6),lp)
          ptab(:) = ptab(:) -pot_param_real(5) * lp(:)
        end if
      case ('[Tong05] He')
        call tong05(tong05_he,order,r0,ptab)
      case ('[Tong05] Ne')
        call tong05(tong05_ne,order,r0,ptab)
      case ('[Tong05] Ar')
        call tong05(tong05_ar,order,r0,ptab)
      case ('[Tong05] Xe')
        call tong05(tong05_xe,order,r0,ptab)
      case ('[SFG99] Li')
        call tong05(sfg99_li,order,r0,ptab)
      case ('[SFG99] Na')
        call tong05(sfg99_na,order,r0,ptab)
      case ('[SFG99] K' )
        call tong05(sfg99_k, order,r0,ptab)
      case ('[SFG99] Rb')
        call tong05(sfg99_rb,order,r0,ptab)
      case ('[SFG99] Cs')
        call tong05(sfg99_cs,order,r0,ptab)
      case ('Test He')
        call tong05(test_he,order,r0,ptab)
      case ('Test 3He')
        call tong05(test_3he,order,r0,ptab)
    end select
  end subroutine potential_u_r
  !
  !  Radial integration from the origin requires power series with respect
  !  to a different parameter:
  !
  !   r+dr = 0.5*eta**2 + 0.5*(x+dx)**2
  !
  !  Evaluate the corresponding coefficients using chain rule.
  !
  subroutine potential_u_halfx2(eta,x0,order,ptab)
    real(rk), intent(in)    :: eta           ! The value of the orthogonal parameter
    real(rk), intent(in)    :: x0            ! Expansion point
    integer(ik), intent(in) :: order         ! Desired derivative order
    real(rk), intent(out)   :: ptab(0:order) ! ptab(n) = n! d^n u(0.5*x**2)/d x^n, so that
                                             ! u(x-x0) = Sum ptab(n) (x-x0)^n
    !
    integer(ik) :: m, k
    real(rk)    :: prtab(0:order)
    real(rk)    :: logx0, log2
    real(rk)    :: fact
    !
    call potential_u_r(0.5_rk*(eta**2+x0**2),order,prtab)
    if (x0==0._rk) then
      !
      !  x0 = 0 requires special treatment in the chain rule
      !
      ptab(:) = 0
      ptab_at_zero: do m=0,order,2
        ptab(m) = prtab(m/2) * 0.5_rk**(m/2)
      end do ptab_at_zero
    else
      !
      !  General case
      !  There is a real potential for overflow here; let's use logarithmic form
      !
      logx0 = log(x0)
      log2  = log(2._rk)
      order_x: do m=0,order
        ptab(m) = 0
        order_r: do k=(m+1)/2,m
          if (prtab(k)==0._rk) cycle order_r
          fact    = MathLogFactorial(k) + real(2*k-m,kind=rk)*logx0 &
                  - MathLogFactorial(m-k) - MathLogFactorial(2*k-m) &
                  - real(m-k,kind=rk)*log2
          ptab(m) = ptab(m) + prtab(k)*exp(fact)
        end do order_r
      end do order_x
    end if
  end subroutine potential_u_halfx2
  !
  !  Our potential must vanish at the matching point: asymptotic solutions
  !  assume purely hydrogenic potential. Although we could include the long-
  !  range (inverse-power of r) terms in the asymptotic solution, we can't
  !  handle then in the coupling equations - so there is no point in doing
  !  this.
  !
  subroutine check_asymptotic_potential
    real(rk)     :: ptab(0:asymp_order)
    real(rk)     :: max_ptab
    integer(ik)  :: order
    !
    call potential_u_halfx2(0._rk,eta_max,asymp_order,ptab)
    max_ptab = maxval(abs(ptab))
    if (verbose>=0) then
      write (out,"('Asymptotic continuum solutions are used from point eta = ',g26.16e3)") eta_max
      write (out,"('Effective potential u = -r*v(r)-Q and its derivatives at the matching point:')")
      write (out,"(3x,a5,1x,a26)") 'Order', 'o! (d^o/d x^o) u(x^2/2)', &
                                   '-----', '-----------------------'
      report_potential: do order=0,asymp_order
        write (out,"(3x,i5,1x,g28.16e4)") order, ptab(order)
      end do report_potential
      write (out,"(/'Max = ',g28.16e4/)") max_ptab
    end if
    if (max_ptab>spacing(1._rk)) then
      write (out,"(/'WARNING: Short-range potential is not negligible at the asymptotic matching point.')")
      write (out,"( 'WARNING: Diagnostic = ',g28.16e4/)") max_ptab
    end if
  end subroutine check_asymptotic_potential
  !
  subroutine cache_potential
    integer(ik) :: xi_pt, eta_pt
    real(rk)    :: eta, xi
    !
    call TimerStart('Cache potential')
    !
    !$omp parallel do collapse(2) default(none) &
    !$omp& shared(eta_npts,xi_npts,xi_maxorder,eta_tab,xi_tab,tab_u_halfx2) &
    !$omp& private(eta_pt,eta,xi_pt,xi)
    eta_points: do eta_pt=1,eta_npts
      xi_points: do xi_pt=1,xi_npts
        eta = eta_tab(eta_pt)
        xi  = xi_tab(xi_pt)
        call potential_u_halfx2(eta,xi,xi_maxorder,tab_u_halfx2(:,xi_pt,eta_pt))
      end do xi_points
    end do eta_points
    !$omp end parallel do
    !
    call TimerStop('Cache potential')
  end subroutine cache_potential
end module general_tunnel_potential

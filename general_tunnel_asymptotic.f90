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
!  Routines dealing with the continuum-direction integration for general_tunnel.f90
!
module general_tunnel_asymptotic
  use accuracy
  use constants
  use general_tunnel_data
  use timer
  !$ use OMP_LIB
  implicit none
  public
  public rcsid_general_tunnel_asymptotic
  !
  character(len=clen), save :: rcsid_general_tunnel_asymptotic = "$Id: $"
  !
  contains
  !
  !  Constructs an asymptotic solution as a power-series expansion:
  !
  !   sqrt(x)*psi = Sum cn(k) x**(-k) exp( pm I [(sqrt(F)/3)*x**3 + (E/sqrt(F))*x ] )
  !
  !  where pm = +1 for the outgoing solution and -1 for the incoming wave.
  !
  !  The coefficients cn(k) are determined using recursively.
  !
  !  The coefficients in the power series can grow rather large. In order to avoid
  !  exceeding the dynamic range of the underlying floating-point type, we'll
  !  scale coefficient cn(k) by x**-k, where x is the coordinate where we plan
  !  to evaluate the series. This way, all we need to do at the end is to sum
  !  up the coefficients.
  !
  subroutine asymptotic_solution(efield,energy,zeta2,pm,x,cn)
    real(rk), intent(in)     :: efield             ! Electric field strength
    complex(rk), intent(in)  :: energy             ! Energy of the solution
    complex(rk), intent(in)  :: zeta2              ! Separation parameter 
    real(rk), intent(in)     :: pm                 ! +1 for outgoing solution; -1 for an incoming
    real(rk), intent(in)     :: x                  ! Coordinate where we need to evaluate the series
    complex(rk), intent(out) :: cn(1:asymp_order)  ! Power series coefficients for the solution
                                                   ! cn(1) is arbitrarily chosen as 1.
    !
    integer(ik) :: l
    complex(rk) :: temp
    real(rk)    :: xm    ! 1/x
    !
    if (efield<=0) stop 'general_tunnel%asymptotic_solution - efield must be positive'
    !
    xm    = 1._rk/x
    cn(1) = xm
    recurrence: do l=1,asymp_order-1
               temp =        pm*(0,1)*((energy**2)/efield - zeta2)      * cn(l)   * xm
      if (l>1) temp = temp - 2*(energy/sqrt(efield))*(l-1)              * cn(l-1) * xm**2
      if (l>2) temp = temp - pm*(0,1)*((l-2)*(l-1) + (0.25_rk-mval**2)) * cn(l-2) * xm**3
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
  subroutine asymptotic_gradient(efield,energy,pm,x,cn,gr)
    real(rk), intent(in)     :: efield                ! Electric field strength
    complex(rk), intent(in)  :: energy                ! Energy of the solution
    real(rk), intent(in)     :: pm                    ! +1 for outgoing solution; -1 for an incoming
    real(rk), intent(in)     :: x                     ! Coordinate where we need to evaluate the series
    complex(rk), intent(in)  :: cn( 1:asymp_order)    ! Solution coefficients
    complex(rk), intent(out) :: gr(-1:asymp_order+1)  ! Gradient coefficients
    !
    integer(ik) :: l
    complex(rk) :: temp
    real(rk)    :: xm     ! 1/x
    !
    if (efield<=0) stop 'general_tunnel%asymptotic_gradient - efield must be positive'
    !
    xm = 1._rk / x
    gradient: do l=-1,asymp_order+1
      temp = 0
      if (l>= 1 .and. l<=asymp_order  ) temp = temp + pm*(0,1)*(energy/sqrt(efield)) * cn(l)
      if (l>= 2 .and. l<=asymp_order+1) temp = temp - (l-1)                          * cn(l-1) * xm
      if (l>=-1 .and. l<=asymp_order-2) temp = temp + pm*(0,1)*sqrt(efield)          * cn(l+2) * x**2
      gr(l) = temp
    end do gradient
    !
    if (verbose>=4) then
      write (out,"(/'Asymptotic gradient for pm=',f3.0)") pm
      write (out,"(2(1x,g34.22e3,1x,g34.22e3,1x))") gr
      write (out,"()")
    end if
  end subroutine asymptotic_gradient
  !
  function evaluate_asymptote(efield,energy,pm,cn,x) result (v)
    real(rk), intent(in)     :: efield ! Electric field strength
    complex(rk), intent(in)  :: energy ! Energy of the solution
    real(rk), intent(in)     :: pm     ! +1 for outgoing solution; -1 for an incoming
    complex(rk), intent(in)  :: cn(:)  ! Solution coefficients
    real(rk), intent(in)     :: x      ! Coordinate where we need the solution
    complex(rk)              :: v
    !
    complex(rk) :: phase
    !
    v     = sum(cn)
    phase = (sqrt(efield)/3._rk)*x**3 + (energy/sqrt(efield))*x
    v     = v * exp(pm*(0,1)*phase)
  end function evaluate_asymptote
  !
end module general_tunnel_asymptotic

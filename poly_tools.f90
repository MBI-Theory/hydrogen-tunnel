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
!  Simple tools for dealing with polynomials
!
module poly_tools
  use accuracy
  use constants
  implicit none
  private
  public poly_adjust_scale, poly_phase_unwrap, poly_power_series, poly_multiply, poly_integrate
  public poly_direction_angle
  public rcsid_poly_tools
  !
  character(len=clen), save :: rcsid_poly_tools = "$Id: $"
  !
  !  Adjust dynamic range of a block of data
  !
  interface poly_adjust_scale
    module procedure poly_adjust_scale_complex_array
    module procedure poly_adjust_scale_complex
  end interface poly_adjust_scale
  !
  contains
  !
  !  Adjust dynamic range of the data so that the largest component
  !  (real or imaginary) of the first element of the array has zero 
  !  exponent.
  !
  subroutine poly_adjust_scale_complex(value,scl)
    complex(rk), intent(inout) :: value  ! Singe value to adjust
    integer(ik), intent(inout) :: scl    ! Overall exponent
    !
    real(rk)    :: re, im, big
    integer(ik) :: extra
    !
    re    = real (value,kind=rk)
    im    = aimag(value)
    big   = max(abs(re),abs(im))
    extra = exponent(big)
    scl   = scl + extra
    re    = scale(re,-extra)
    im    = scale(im,-extra)
    value = cmplx(re,im,kind=rk)
  end subroutine poly_adjust_scale_complex
  !
  subroutine poly_adjust_scale_complex_array(values,scl)
    complex(rk), intent(inout) :: values(:)  ! Values to adjust
    integer(ik), intent(inout) :: scl        ! Overall exponent
    !
    integer(ik) :: first_scale, i
    real(rk)    :: re, im
    !
    if (size(values)<1) return ! Nothing to do
    !
    first_scale = 0
    call poly_adjust_scale(values(1),first_scale)
    scl = scl + first_scale
    !
    scale_rest: do i=2,size(values)
      re        = real (values(i),kind=rk)
      im        = aimag(values(i))
      re        = scale(re,-first_scale)
      im        = scale(im,-first_scale)
      values(i) = cmplx(re,im,kind=rk)
    end do scale_rest
  end subroutine poly_adjust_scale_complex_array
  !
  subroutine poly_phase_unwrap(phases)
    real(rk), intent(inout) :: phases(:)  ! Phases to unwrap
    !
    integer(ik) :: ipt, npt
    real(rk)    :: d1, d2 ! Possible changes in the phase
    !
    npt = size(phases)
    unwrap_phases: do ipt=2,npt
      d1 = modulo(phases(ipt) - phases(ipt-1),pi)
      d2 = d1 - pi
      if (abs(d2)<d1) then ! d1 is guaranteed to be non-negative
        d1 = d2
      end if
      phases(ipt) = phases(ipt-1) + d1
    end do unwrap_phases
  end subroutine poly_phase_unwrap
  !
  !  Calculate: Sum cn(i) * dx**(i-1) and Sum cn(i)*(i-1)*x**(i-2)
  !
  subroutine poly_power_series(cn,dx,fg)
    complex(rk), intent(in)  :: cn(:) ! Power series coefficients
    real(rk),    intent(in)  :: dx    ! Step
    complex(rk), intent(out) :: fg(:) ! Function value and gradient
    !
    integer(ik) :: it
    complex(rk) :: t
    !
    if (size(fg)>=1) then
      t = cn(size(cn))
      descend_terms_function: do it=size(cn)-1,1,-1
        t = cn(it) + dx*t
      end do descend_terms_function
      fg(1) = t
      !
      if (size(fg)>=2) then
        t = cn(size(cn))*(size(cn)-1)
        descend_terms_gradient: do it=size(cn)-1,2,-1
          t = cn(it)*(it-1) + dx*t
        end do descend_terms_gradient
        fg(2) = t
      end if
    end if
  end subroutine poly_power_series
  !
  subroutine poly_multiply(p1,p2,pp)
    complex(rk), intent(in)  :: p1(:), p2(:) ! Polynomial coefficients to multiply. First coefficient must be for x^0
    complex(rk), intent(out) :: pp(:)        ! Product polynomial
    !
    integer(ik) :: sz1, sz2, szp ! Sizes of the polynomials
    integer(ik) :: i2            ! Index within the second polynomial
    !
    sz1 = size(p1)
    sz2 = size(p2)
    szp = size(pp)
    if (szp<sz1+sz2-1) stop 'poly_tools%poly_multiply - output array too small'
    !
    pp = 0
    scan_terms2: do i2=1,sz2
      pp(i2:i2+sz1-1) = pp(i2:i2+sz1-1) + p1(:)*p2(i2)
    end do scan_terms2
  end subroutine poly_multiply
  !
  subroutine poly_integrate(p1,pint)
    complex(rk), intent(in)  :: p1(:)    ! Polynomial to integrate. First coefficient must be for x^0.
    complex(rk), intent(out) :: pint(:)  ! The integral. Must have one element more than p1(:)
    !
    integer(ik) :: i1
    !
    if (size(pint)/=size(p1)+1) stop 'poly_tools%poly_integrate - bad array sizes'
    !
    pint(1) = 0
    scan_terms: do i1=1,size(p1)
      pint(i1+1) = p1(i1)/real(i1,kind=rk)
    end do scan_terms
  end subroutine poly_integrate
  !
  !  This is not really a polynomial-handling function, but ...
  !
  function poly_direction_angle(v1,v2) result(cos_ang)
    complex(rk), intent(in) :: v1, v2   ! Two complex numbers, treated as vectors
    real(rk)                :: cos_ang  ! cosine of the angle between v1 and v2 directions
    complex(rk)             :: a, b, c  ! a = v1/|v1|; b = v2/|v2|; c = b-a
    real(rk)                :: r1, r2, r3
    !
    r1 = abs(v1)
    r2 = abs(v2)
    if (r1==0._rk .or. r2==0._rk) then
      cos_ang = 1._rk ! Degenerate case; zero-length vectors are always collinear with all other vectors
      return
    end if
    a       = v1/r1   ! |a| = 1 by construction
    b       = v2/r2   ! |b| = 1
    c       = b-a
    r3      = abs(c)  ! r3^2 = 2 - 2*cos(alpha)
    cos_ang = 1._rk - 0.5_rk*r3**2
    if (cos_ang<-1._rk) cos_ang = -1._rk
  end function poly_direction_angle
end module poly_tools

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
! function MathInterpolateComplex(x,xtab,vtab) result (v)
!   real(xk), intent(in)    :: x       ! Point where the interpolant is desired
!   real(xk), intent(in)    :: xtab(:) ! Positions to interpolate
!   complex(xk), intent(in) :: vtab(:) ! Values to interpolate
!   complex(xk)             :: v       ! Interpolant
!   !
!   complex(xk) :: c(size(xtab)), d(size(xtab)), scl
    real(xk)    :: ho, hp, den
    integer(ik) :: ipt, icol, npts, nelem, i
    !
    npts = size(xtab)
    if (size(vtab)/=npts) stop 'math%MathInterpolate - bad input array sizes'
    !
    c    = vtab
    d    = vtab
    ipt  = minloc(abs(x-xtab),dim=1)
    v    = vtab(ipt)
    ipt  = ipt - 1
    tableau_columns: do icol=1,npts-1
      nelem = npts-icol
      update_tableau: do i=1,nelem
        ho   = xtab(i)-x
        hp   = xtab(i+icol)-x
        den  = xtab(i) - xtab(i+icol)
        if (den==0._xk) stop 'math%MathInterpolate - division by zero'
        scl  = (c(i+1) - d(i)) / den
        d(i) = hp*scl
        c(i) = ho*scl
      end do update_tableau
      if (2*ipt<nelem) then
        v   = v + c(ipt+1)
      else
        v   = v + d(ipt)
        ipt = ipt - 1
      end if
    end do tableau_columns
! end function MathInterpolateComplex

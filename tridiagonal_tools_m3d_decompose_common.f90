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
! subroutine m3d_decompose_r(m,mf,fail)
!   real(rk), intent(in)           :: m (:,:) ! Tridiagonal matrix, assumed to be diagonally dominant
!                                             ! As the result, we do not need to bother with pivoting
!   real(rk), intent(out)          :: mf(:,:) ! Factorized tridiagonal matrix
!   logical, intent(out), optional :: fail    ! Set to .true. if decomposition fails; 
!                                             ! if fail is absent, abort on decomposition failure.
!   !
!   real(rk)    :: denom
    integer(ik) :: i, sz
    !
    sz = size(m,dim=2)
    if (size(m,dim=1)<3 .or. size(mf,dim=1)<3 .or. size(mf,dim=2)/=sz) then
      stop 'tridiagonal_tools%m3d_decompose_common - bad input sizes'
    end if
    if (present(fail)) fail = .false.
    !
    mf(1,1) = 1._rk/m(1,1)
    mf(2,1) = 0._rk
    mf(3,1) = m(3,1)*mf(1,1)
    factor_m3d: do i=2,sz-1
      denom   = m(1,i)-m(2,i-1)*mf(3,i-1)
      if (too_small(abs(denom))) return
      mf(1,i) = 1._rk/denom
      mf(2,i) = -m(2,i-1)*mf(1,i)
      mf(3,i) = m(3,i)*mf(1,i)
    end do factor_m3d
    if (sz<2) return
    denom    = m(1,sz)-m(2,sz-1)*mf(3,sz-1)
    if (too_small(abs(denom))) return
    mf(1,sz) = 1._rk/denom
    mf(2,sz) = -m(2,sz-1)*mf(1,sz)
    mf(3,sz) = 0._rk
    !
    contains
    logical function too_small(x)
      real(kind(m)), intent(in) :: x
      !
      too_small = .false.
      if (x>100*tiny(x)) return
      too_small = .true.
      if (present(fail)) then
        fail = .true.
      else
        write (out,"('Fatal error in m3d_decompose_common: denominator ',g34.16e3,' is too small.')") x
        stop 'tridiagonal_tools%m3d_decompose_common - decomposition failed'
      end if
    end function too_small
! end subroutine m3d_decompose_r

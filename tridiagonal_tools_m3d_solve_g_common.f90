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
! subroutine m3d_solve_rr_g(mf,r,x)
!   real(rk), intent(in)  :: mf(:,:) ! Factorized tridiagonal matrix
!   real(rk), intent(in)  :: r (:,:) ! Right-hand sizes
!   real(rk), intent(out) :: x (:,:) ! Solution vectors
    !
    integer(ik) :: i, sz, ir
    !
    sz = size(mf,dim=2)
    if (size(mf,dim=1)<3 .or. size(r,dim=1)/=sz .or. size(x,dim=1)/=sz .or. size(x,dim=2)/=size(x,dim=2)) then
      stop 'tridiagonal_tools%m3d_solve_g_common - bad input sizes'
    end if
    ! A decent compiler will rearrange this loop nest as needed ....
    scan_rhs: do ir=1,size(r,dim=2)
      x(1,ir) = mf(1,1)*r(1,ir)
      transform_rhs: do i=2,sz
        x(i,ir) = mf(1,i)*r(i,ir) + mf(2,i)*x(i-1,ir)
      end do transform_rhs
      !
      backsubstitute: do i=sz-1,1,-1
        x(i,ir) = x(i,ir) - mf(3,i)*x(i+1,ir)
      end do backsubstitute
    end do scan_rhs
! end subroutine m3d_solve_rr_g

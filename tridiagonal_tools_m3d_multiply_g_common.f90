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
! subroutine m3d_multiply_rr_g(m,v,mv)
!   real(rk), intent(in)  :: m(:,:)  ! Tri-diagonal matrix
!   real(rk), intent(in)  :: v(:,:)  ! Vectors
!   real(rk), intent(out) :: mv(:,:) ! Vectors
    !
    integer(ik) :: sz, iv
    !
    sz = size(v,dim=1)
    if (size(m,dim=1)<3 .or. size(m,dim=2)/=sz .or. size(mv,dim=1)/=sz .or. size(v,dim=2)/=size(mv,dim=2)) then
      stop 'tridiagonal_tools%m3d_multiply_g_common - bad input sizes'
    end if
    scan_vectors: do iv=1,size(v,dim=2)
      mv(1:sz  ,iv) =                 m(1,:    )*v(1:sz  ,iv) 
      mv(2:sz  ,iv) = mv(2:sz  ,iv) + m(2,:sz-1)*v(1:sz-1,iv)
      mv(1:sz-1,iv) = mv(1:sz-1,iv) + m(3,:sz-1)*v(2:sz  ,iv)
    end do scan_vectors
! end subroutine m3d_multiply_rr_g

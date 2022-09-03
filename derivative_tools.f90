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
!  Implementation of Savitzky-Golay smoothing filter (NR 14.8).
!  Here, we use it as a way of calculating arbitrary-order derivatives 
!  for data on a uniform grid.
!
!  Additionally, a few hard-coded formulae for differentiation on uniform 3d grids.
!
module derivative_tools
  use accuracy
  use lapack
  implicit none
  private
  public dt_savgol
  public dt_lap3, dt_lap5, dt_lap7
  public dt_grad3, dt_grad5, dt_grad7
  public rcsid_derivative_tools
  !
  interface dt_lap3
    module procedure lap3_real
    module procedure lap3_complex
  end interface dt_lap3
  interface dt_lap5
    module procedure lap5_real
    module procedure lap5_complex
  end interface dt_lap5
  interface dt_lap7
    module procedure lap7_real
    module procedure lap7_complex
  end interface dt_lap7
  !
  interface dt_grad3
    module procedure grad3_real
  end interface dt_grad3
  interface dt_grad5
    module procedure grad5_real
  end interface dt_grad5
  interface dt_grad7
    module procedure grad7_real
  end interface dt_grad7
  !
  character(len=clen), save :: rcsid_derivative_tools = "$Id: $"
  !
  contains
  !
  !  3-point finite difference
  !
  function lap3_real(dx,f) result(lap)
    real(rk), intent(in) :: dx                ! Grid spacing
    real(rk), intent(in) :: f(-1:1,-1:1,-1:1) ! Function values on a grid
    real(rk)             :: lap               ! Laplacian at the central point
    !
    real(rk), parameter :: w3c = -6._rk
    real(rk), parameter :: w3o =  1._rk
    !
    lap = w3c * f( 0, 0, 0) &
        + w3o * f(-1, 0, 0) + w3o * f(+1, 0, 0) &
        + w3o * f( 0,-1, 0) + w3o * f( 0,+1, 0) &
        + w3o * f( 0, 0,-1) + w3o * f( 0, 0,+1)
    lap = lap / dx**2
  end function lap3_real
  !
  function lap3_complex(dx,f) result(lap)
    real(rk), intent(in)    :: dx                ! Grid spacing
    complex(rk), intent(in) :: f(-1:1,-1:1,-1:1) ! Function values on a grid
    complex(rk)             :: lap               ! Laplacian at the central point
    !
    real(rk), parameter :: w3c = -6._rk
    real(rk), parameter :: w3o =  1._rk
    !
    lap = w3c * f( 0, 0, 0) &
        + w3o * f(-1, 0, 0) + w3o * f(+1, 0, 0) &
        + w3o * f( 0,-1, 0) + w3o * f( 0,+1, 0) &
        + w3o * f( 0, 0,-1) + w3o * f( 0, 0,+1)
    lap = lap / dx**2
  end function lap3_complex
  !
  !  5-point finite difference
  !
  function lap5_real(dx,f) result(lap)
    real(rk), intent(in) :: dx                ! Grid spacing
    real(rk), intent(in) :: f(-2:2,-2:2,-2:2) ! Function values on a grid
    real(rk)             :: lap               ! Laplacian at the central point
    !
    real(rk), parameter :: w5c = -(90._rk/12._rk)
    real(rk), parameter :: w5s =  (16._rk/12._rk)
    real(rk), parameter :: w5d = -( 1._rk/12._rk)
    !
    lap = w5c * f( 0, 0, 0) &
        + w5s * f(-1, 0, 0) + w5s * f(+1, 0, 0) &
        + w5s * f( 0,-1, 0) + w5s * f( 0,+1, 0) &
        + w5s * f( 0, 0,-1) + w5s * f( 0, 0,+1) &
        + w5d * f(-2, 0, 0) + w5d * f(+2, 0, 0) &
        + w5d * f( 0,-2, 0) + w5d * f( 0,+2, 0) &
        + w5d * f( 0, 0,-2) + w5d * f( 0, 0,+2)
    lap = lap / dx**2
  end function lap5_real
  !
  function lap5_complex(dx,f) result(lap)
    real(rk), intent(in)    :: dx                ! Grid spacing
    complex(rk), intent(in) :: f(-2:2,-2:2,-2:2) ! Function values on a grid
    complex(rk)             :: lap               ! Laplacian at the central point
    !
    real(rk), parameter :: w5c = -(90._rk/12._rk)
    real(rk), parameter :: w5s =  (16._rk/12._rk)
    real(rk), parameter :: w5d = -( 1._rk/12._rk)
    !
    lap = w5c * f( 0, 0, 0) &
        + w5s * f(-1, 0, 0) + w5s * f(+1, 0, 0) &
        + w5s * f( 0,-1, 0) + w5s * f( 0,+1, 0) &
        + w5s * f( 0, 0,-1) + w5s * f( 0, 0,+1) &
        + w5d * f(-2, 0, 0) + w5d * f(+2, 0, 0) &
        + w5d * f( 0,-2, 0) + w5d * f( 0,+2, 0) &
        + w5d * f( 0, 0,-2) + w5d * f( 0, 0,+2)
    lap = lap / dx**2
  end function lap5_complex
  !
  !  7-point finite difference
  !
  function lap7_real(dx,f) result(lap)
    real(rk), intent(in) :: dx                ! Grid spacing
    real(rk), intent(in) :: f(-3:3,-3:3,-3:3) ! Function values on a grid
    real(rk)             :: lap               ! Laplacian at the central point
    !
    real(rk), parameter :: w7c = -(1470._rk/180._rk)
    real(rk), parameter :: w7s =  ( 270._rk/180._rk)
    real(rk), parameter :: w7d = -(  27._rk/180._rk)
    real(rk), parameter :: w7t =  (   2._rk/180._rk)
    !
    lap = w7c * f( 0, 0, 0) &
        + w7s * f(-1, 0, 0) + w7s * f(+1, 0, 0) &
        + w7s * f( 0,-1, 0) + w7s * f( 0,+1, 0) &
        + w7s * f( 0, 0,-1) + w7s * f( 0, 0,+1) &
        + w7d * f(-2, 0, 0) + w7d * f(+2, 0, 0) &
        + w7d * f( 0,-2, 0) + w7d * f( 0,+2, 0) &
        + w7d * f( 0, 0,-2) + w7d * f( 0, 0,+2) &
        + w7t * f(-3, 0, 0) + w7t * f(+3, 0, 0) &
        + w7t * f( 0,-3, 0) + w7t * f( 0,+3, 0) &
        + w7t * f( 0, 0,-3) + w7t * f( 0, 0,+3)
    lap = lap / dx**2
  end function lap7_real
  !
  function lap7_complex(dx,f) result(lap)
    real(rk), intent(in)    :: dx                ! Grid spacing
    complex(rk), intent(in) :: f(-3:3,-3:3,-3:3) ! Function values on a grid
    complex(rk)             :: lap               ! Laplacian at the central point
    !
    real(rk), parameter :: w7c = -(1470._rk/180._rk)
    real(rk), parameter :: w7s =  ( 270._rk/180._rk)
    real(rk), parameter :: w7d = -(  27._rk/180._rk)
    real(rk), parameter :: w7t =  (   2._rk/180._rk)
    !
    lap = w7c * f( 0, 0, 0) &
        + w7s * f(-1, 0, 0) + w7s * f(+1, 0, 0) &
        + w7s * f( 0,-1, 0) + w7s * f( 0,+1, 0) &
        + w7s * f( 0, 0,-1) + w7s * f( 0, 0,+1) &
        + w7d * f(-2, 0, 0) + w7d * f(+2, 0, 0) &
        + w7d * f( 0,-2, 0) + w7d * f( 0,+2, 0) &
        + w7d * f( 0, 0,-2) + w7d * f( 0, 0,+2) &
        + w7t * f(-3, 0, 0) + w7t * f(+3, 0, 0) &
        + w7t * f( 0,-3, 0) + w7t * f( 0,+3, 0) &
        + w7t * f( 0, 0,-3) + w7t * f( 0, 0,+3)
    lap = lap / dx**2
  end function lap7_complex
  !
  !  1D gradient on a uniform grid
  !
  function grad3_real(dx,f) result(grad)
    real(rk), intent(in) :: dx      ! Grid spacing
    real(rk), intent(in) :: f(-1:1) ! Function values on a grid
    real(rk)             :: grad    ! Gradient at the central point
    !
    real(rk), parameter :: w3 = 0.5_rk
    !
    grad = w3*f(+1) - w3*f(-1)
    grad = grad / dx
  end function grad3_real
  !
  function grad5_real(dx,f) result(grad)
    real(rk), intent(in) :: dx      ! Grid spacing
    real(rk), intent(in) :: f(-2:2) ! Function values on a grid
    real(rk)             :: grad    ! Gradient at the central point
    !
    real(rk), parameter :: w5s =  2._rk/3._rk
    real(rk), parameter :: w5d = -1._rk/12._rk
    !
    grad = w5d*f(+2) + w5s*f(+1) - w5s*f(-1) - w5d*f(-2)
    grad = grad / dx
  end function grad5_real
  !
  function grad7_real(dx,f) result(grad)
    real(rk), intent(in) :: dx      ! Grid spacing
    real(rk), intent(in) :: f(-3:3) ! Function values on a grid
    real(rk)             :: grad    ! Gradient at the central point
    !
    real(rk), parameter :: w7s =  3._rk/4._rk
    real(rk), parameter :: w7d = -3._rk/20._rk
    real(rk), parameter :: w7t =  1._rk/60._rk
    !
    grad = w7t*f(+3) + w7d*f(+2) + w7s*f(+1) - w7s*f(-1) - w7d*f(-2) - w7t*f(-3)
    grad = grad / dx
  end function grad7_real
  !
  !  Evaluate coefficients in the polynomial expansion:
  !
  !     f(x) = Sum D[j] x^j 
  !
  !  from values of f(x) on a uniform grid.
  !
  !  WARNING: This routine gets rather unstable for high orders - beware.
  !
  subroutine dt_savgol(nord,nl,d,min_der,coeff)
                                          ! npts is the first dimension of coeff(:,:)
    integer(ik), intent(in) :: nord       ! Maximum order of the filter. Must be at least 1.
    integer(ik), intent(in) :: nl         ! Number of grid points to the left of the point where derivative(s) are needed
    real(rk), intent(in)    :: d          ! Uniform grid spacing
    integer(ik), intent(in) :: min_der    ! Lowest order of the derivative needed; min_der=0 is the function itself
    real(rk), intent(out)   :: coeff(:,:) ! First index: point index. 
                                          ! 1 = the left-most point (nl points to the left of the reference point), etc
                                          ! Second index: derivative order. 
                                          ! 1 = min_der, etc
    !
    real(rk), allocatable :: amat(:,:)    ! The design matrix
    real(rk), allocatable :: bmat(:,:)    ! (A^T A) and its inverse
    real(rk), allocatable :: cmat(:,:)    ! The unscaled fit matrix
    integer(ik)           :: npts         ! Total width of the filter, in points
    integer(ik)           :: nder         ! Total number of derivative orders needed
    integer(ik)           :: ipt, icol, iord, alloc
    real(rk)              :: scl
    !
    !  A bit of sanity checking
    !
    npts = size(coeff,dim=1)
    nder = size(coeff,dim=2)
    if (npts<=0 .or. nord<=0 .or. nord>npts-1 .or. nl>npts-1 .or. nl<0 .or. min_der<0 .or. min_der+nder-1>nord .or. d<=0) then
      write (out,"('ERROR: Inconsistent arguments in a call to derivative_tools%dt_savgol')")
      write (out,"('npts = ',i0)") npts
      write (out,"('nder = ',i0)") nder
      write (out,"('nl   = ',i0)") nl
      write (out,"('iord = ',i0)") iord
      write (out,"('d    = ',g16.6)") d
      stop 'derivative_tools%dt_savgol - bad arguments'
    end if
    !
    allocate (amat(npts,nord+1),bmat(nord+1,nord+1),cmat(nord+1,npts),stat=alloc)
    if (alloc/=0) then
      stop 'derivative_tools%dt_savgol - allocation failed'
    end if
    !
    !  Fill the design matrix
    !
    amat(:,1) = 1._rk ! Zeroth order of anything is one
    fill_displacements: do ipt=1,npts
      amat(ipt,2) = (ipt-1_ik) - nl
    end do fill_displacements
    fill_powers: do icol=3,nord+1
      amat(:,icol) = amat(:,icol-1) * amat(:,2)
    end do fill_powers
    !
    !  Solve for the fit matrix: (A^T A)^{-1} A^T
    !
    bmat = matmul(transpose(amat),amat)
    call lapack_ginverse(bmat,eps_=0.0_rk)
    cmat = matmul(bmat,transpose(amat))
    !
    copy_weights: do iord=min_der,min_der+nder-1
      scl = 1._rk
      if (iord>0) scl = (1._rk/d)**iord
      coeff(:,iord-min_der+1) = scl * cmat(iord+1,:)
    end do copy_weights
    !
    deallocate (amat,bmat,cmat)
  end subroutine dt_savgol
end module derivative_tools
!
!program test
!  use accuracy
!  use derivative_tools
!  !
!  integer(ik) :: npts, nord, nl, min_der, nder, ipt, ider
!  real(rk)    :: d
!  real(rk), allocatable :: coeff(:,:)
!  !
!  read(input,*) npts, nord, nl, d, min_der, nder
!  !
!  allocate (coeff(npts,nder))
!  !
!  call dt_savgol(nord,nl,d,min_der,coeff)
!  !
!  do ider=1,nder
!    write (out,*) ' derivative order = ',min_der+ider-1
!    write (out,*) coeff(:,ider)
!  end do
!end program test

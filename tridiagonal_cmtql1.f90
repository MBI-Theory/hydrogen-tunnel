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
!  Implementation of CMTQL1 routine from Cullum and Willoughby, 
!  "Lanczos Algorithms for Large Symmetric Eigenvalue Computations",
!  vol. 2
!
!  The code below closely follows the routine on pages 504-507 of the
!  book available at www.netlib.org/lanczos/
!
module tridiagonal_cmtql1
  use accuracy
  implicit none
  private
  public cmtql1
  public rcsid_tridiagonal_cmtql1
  !
  character(len=clen), save :: rcsid_tridiagonal_cmtql1 = "$Id: $"
  !
  integer(ik), parameter :: max_iterations = 100_ik
  !
  contains

  subroutine cmtql1(d,e,ierr)
    complex(rk), intent(inout) :: d(:)  ! In:  Diagonal of the matrix
                                        ! Out: Eigenvalues
    complex(rk), intent(inout) :: e(:)  ! In:  e(1:size(d)-1) - subdiagonal
                                        !      The total size must be the same as d(:)
                                        ! Out: Destroyed
    integer(ik), intent(out)   :: ierr  ! Error code. 
                                        ! 0  = successful completion.
                                        ! +n = eigenvalue (n) failed to converge after (max_iterations) iterations
                                        ! -n = QL decomposition failed for eigenvalue (n)
    !
    integer(ik) :: ntot                ! Dimension of the matrix
    integer(ik) :: iev                 ! Eigenvalue we are working on
                                       ! d(1:iev-1) contains converged eigenvalues
    integer(ik) :: itop                ! Start of the first diagonal subblock above the current eigenvalue
    integer(ik) :: iter                ! Iteration count for the current eigenvalue
    integer(ik) :: ic                  ! Currenty active column
    complex(rk) :: shift               ! value of the shift
    complex(rk) :: g, s, c, p, f, b, r !
    !
    ntot = size(d)
    if (size(e)/=ntot) then
      stop 'tridiagonal_cmtql1%cmtql1 - Bad array sizes'
    end if
    ierr = 0
    if (ntot<=1) return
    !
    e(ntot) = 0
    !
    eigenvalues: do iev=1,ntot-1
      iter = 0
      iterations: do 
        !
        !  If the 2x2 sub-block starting at iev is diagonal, current eigenvalue converged
        !
        if (is_diagonal2x2(iev)) then
          call store_eigenvalue(iev)
          cycle eigenvalues
        end if
        iter = iter + 1
        if (iter>max_iterations) then
          ierr = iev
          return
        end if
        !
        !  Try to find the first diagonal 2x2 sub-block. This is our starting point
        !  for the transformation. If no diagonal sub-block exists, we'll start at 
        !  the last column of the matrix.
        !
        find_break: do itop=iev+1,ntot-1
          if (is_diagonal2x2(itop)) exit find_break
        end do find_break
        call closest_2x2eigenvalue(d(iev),d(iev+1),e(iev),shift)
        !
        !  We now perform a QL iteration.
        !
        g = d(itop) - shift
        s =  1._rk
        c = -1._rk
        p =  0._rk
        transform: do ic=itop-1,iev,-1
          f =  s*e(ic)
          b = -c*e(ic)
          if (.not.eliminate(f,g,c,s)) then
            ierr = -iev 
            return
          end if
          e(ic+1) = s*f + c*g
          g       = d(ic+1) - p
          r       = (d(ic)-g)*s + 2._rk*c*b
          p       = s*r
          d(ic+1) = g + p
          g       = b - c*r
        end do transform
        d(iev) = d(iev) - p
        e(iev) = g
        e(itop) = 0._rk
      end do iterations
    end do eigenvalues
    !
    !  Once we've reached the last element on the diagonal, it must be an eigenvalue
    !
    call store_eigenvalue(ntot)
    !
    contains 
      !
      !  Return true if the (2x2) sub-block starting at column (ic) is diagonal
      !
      logical function is_diagonal2x2(ic)
        integer(ik), intent(in) :: ic   ! First column of the sub-block
        real(rk)                :: diag ! Sum of the absolute values on the diagonal
        !
        diag           = abs(d(ic)) + abs(d(ic+1))
        is_diagonal2x2 = abs(e(ic)) <= spacing(diag)
        return 
      end function is_diagonal2x2
      !
      !  Diagonalize the 2x2 symmetric matrix, and pick the eigenvalue
      !  closest to the value in the upper left corner.
      !
      subroutine closest_2x2eigenvalue(d1,d2,e,ev)
        complex(rk), intent(in)  :: d1, d2  ! Diagonal of the matrix
        complex(rk), intent(in)  :: e       ! Off-diagonal matrix element
        complex(rk), intent(out) :: ev      ! Eigenvalue of the 2x2 sub-block closest in magnitude to d1
        complex(rk)              :: delta   ! Difference between the diagonal elements
        complex(rk)              :: e2      ! Square of the off-diagonal element
        complex(rk)              :: det     ! Square root of the determinant/4
        complex(rk)              :: lambda1 ! Smaller eigenvalue
        complex(rk)              :: lambda2 ! Larger eigenvalue
        real(rk)                 :: scale   ! Overall scaling factor
        !
        !  For a 2x2 matrix in the form:
        !
        !   / 0 e     \
        !   \ e delta /
        !
        !  the eigenvalues are:
        !
        !    lambda1 = delta + sqrt(delta**2+e**2)
        !    lambda2 = delta - sqrt(delta**2+e**2)
        !
        !  We'll compute the larger eigenvalue directly, then determine the smaller eigenvalue
        !  (which is the one we require) from the usual relation:
        !
        !    lambda1*lambda2 = -e**2
        !
        !  We'll also rescale the 2x2 matrix such that the largest of (e) and (delta) has
        !  unit magnitude. This minimizes chances of overflows.
        delta = d2 - d1
        scale = max(abs(delta),abs(e))
        delta = delta/scale
        e2    = (e/scale)**2
        det   = sqrt(delta**2+e2)
        if (abs(delta+det)>=abs(delta-det)) then
          lambda2 = delta+det
        else
          lambda2 = delta-det
        end if
        lambda1 = -scale*e2/lambda2
        ev = d1 + lambda1
      end subroutine closest_2x2eigenvalue
      !
      !  Solve system of equations:
      !
      !     c*f - s*g = 0
      !     c**2 + s**2 = 1
      !
      !  This system of equations is equivalent to a rotation eliminating
      !  off-diagonal matrix element (f):
      !
      !    / c -s \  / ? f \  =  / ? 0 \
      !    \ s  c /  \ ? g /     \ ? ? /
      !
      !  The formal solution is:
      !
      !    c = g / sqrt(f**2 + g**2)
      !    s = f / sqrt(g**2 + f**2)
      !
      !  Unfortunately, we can't use Fortran 2008 hypot intrinsic, which
      !  is limited to real arguments. Shame.
      !
      !  Unlike the real case, these equations may have no finite solutions.
      !  This happens when f = +/-I*g, and the denominator of the formal
      !  solution vanishes. 
      !
      function eliminate(f,g,c,s) result(ok)
        complex(rk), intent(in)  :: f, g 
        complex(rk), intent(out) :: c, s  ! Values are not defined unless the return value is ok
        logical                  :: ok    ! .true. if sensible solutions exist
        complex(rk)              :: hypot
        !
        if (abs(f)>=abs(g)) then
          hypot = f*sqrt(1._rk+(g/f)**2)
        else
          hypot = g*sqrt(1._rk+(f/g)**2)
        end if
        if (abs(hypot)<=max(abs(f),abs(g))/sqrt(huge(1._rk))) then
          ok = .false.
        else
          ok = .true.
          c  = g / hypot
          s  = f / hypot
        end if
      end function eliminate
      !
      !  Store the new eigenvalue in the desired sorting order.
      !  All prior eigenvalues are already in the desired order.
      !  We'll sort in the order of increasing real part of the eigenvalue.
      !
      subroutine store_eigenvalue(itop)
        integer(ik),  intent(in) :: itop  ! initial location of the new eigenvalue. 
                                          ! eigenvalues 1..itop-1 are already sorted
        integer(ik)              :: ipos  ! insertion position
        complex(rk)              :: ev
        !
        ev = d(itop)
        find_spot: do ipos=1,itop-1
          if (real(ev,kind=rk)<real(d(ipos),kind=rk)) exit find_spot
        end do find_spot
        d(ipos+1:itop) = d(ipos:itop-1)  ! Shift eigenvalues above the insertion point
        d(ipos)        = ev
      end subroutine store_eigenvalue
  end subroutine cmtql1
end module tridiagonal_cmtql1

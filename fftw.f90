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
!  FFTW 3 interface for 1DD complex to complex transforms
!  This is a subset of the fftw_fftw3.f90 module from multigrid
!
!  On input, we take uniformly spaced grids. On output, we
!  spit out frequency components (as you would expect from FFT ;-).
!  For each index, the frequencies increase monotonically from
!  -Pi to Pi. The zero-frequency component is found at element
!  1+(N-1)/2, where N is the number of points along a given dimension
!
 module fftw
   use accuracy
   implicit none

   integer, parameter :: hik         = selected_int_kind(15)      ! "Pointer" integers - sufficient to store
                                                                  ! memory address

   private
   public fftw_1d
   public rcsid_fftw

   interface fftw_1d
     module procedure cfftw_1d_inplace
     module procedure zfftw_1d_inplace
   end interface ! fftw_1d

   integer(ik), parameter :: FFTW_FORWARD        = -1
   integer(ik), parameter :: FFTW_BACKWARD       = +1
   integer(ik), parameter :: FFTW_MEASURE        =  0
   integer(ik), parameter :: FFTW_ESTIMATE       = 64
   integer(ik), parameter :: FFTW_DESTROY_INPUT  =  1
   integer(ik), parameter :: FFTW_PRESERVE_INPUT = 16

   character(len=clen), save :: rcsid_fftw = "$Id: $"

   contains
   !
   subroutine cfftw_1d_inplace(n1,n2,a,invert)
     integer(ik), intent(in)      :: n1, n2
     complex(srk), intent(inout)  :: a(n1,n2) ! Due to the idiotic way FFTW 3 defines plans ...
     logical, intent(in)          :: invert
     !
!*ft integer(hik) :: plan         ! Large enough to store a pointer
!*ft integer(ik)  :: direction, i !
!*ft complex(srk) :: buf(n1)      ! Local buffer for 1D transforms
!*ft !
!*ft direction = FFTW_FORWARD
!*ft if (invert) direction = FFTW_BACKWARD
!*ft !
!*ft call sfftw_plan_dft_1d(plan, n1, buf, buf, direction, FFTW_ESTIMATE+FFTW_DESTROY_INPUT) 
!*ft column_loop: do i=1,n2
!*ft   if (invert) then
!*ft     buf = cshift(a(:,i),+(n1-1)/2)
!*ft   else
!*ft     buf = a(:,i)
!*ft   end if
!*ft   !
!*ft   call sfftw_execute(plan) 
!*ft   !
!*ft   if (invert) then
!*ft     a(:,i) = buf
!*ft   else
!*ft     a(:,i) = cshift(buf,-(n1-1)/2)
!*ft   end if
!*ft end do column_loop
!*ft call sfftw_destroy_plan(plan)
!*ft return
     stop 'fftw%cfftw_1d_inplace - FFT is not enabled at compile time'
   end subroutine cfftw_1d_inplace

   subroutine zfftw_1d_inplace(n1,n2,a,invert)
     integer(ik), intent(in)     :: n1, n2
     complex(drk), intent(inout) :: a(n1,n2) ! Due to the idiotic way FFTW 3 defines plans ...
     logical, intent(in)         :: invert
     !
!*ft integer(hik) :: plan         ! Large enough to store a pointer
!*ft integer(ik)  :: direction, i !
!*ft complex(drk) :: buf(n1)      ! Local buffer for 1D transforms
!*ft !
!*ft direction = FFTW_FORWARD
!*ft if (invert) direction = FFTW_BACKWARD
!*ft !
!*ft call dfftw_plan_dft_1d(plan, n1, buf, buf, direction, FFTW_ESTIMATE+FFTW_DESTROY_INPUT) 
!*ft column_loop: do i=1,n2
!*ft   if (invert) then
!*ft     buf = cshift(a(:,i),+(n1-1)/2)  ! Circular shift to the left
!*ft   else
!*ft     buf = a(:,i)
!*ft   end if
!*ft   !
!*ft   call dfftw_execute(plan) 
!*ft   !
!*ft   if (invert) then
!*ft     a(:,i) = buf
!*ft   else
!*ft     a(:,i) = cshift(buf,-(n1-1)/2)  ! Circular shift to the right
!*ft   end if
!*ft end do column_loop
!*ft call dfftw_destroy_plan(plan)
!*ft return
     stop 'fftw%zfftw_1d_inplace - FFT is not enabled at compile time'
   end subroutine zfftw_1d_inplace

 end module fftw

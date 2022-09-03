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
!  Report versions of all modules in this build
!
 subroutine versions
   use accuracy
   use constants
   use derivative_tools
   use find_minimum
   use fftw
   use find_root
   use lapack
   use math
   use poly_tools
   use sort_tools
   use timer
   use tridiagonal_tools
   use tridiagonal_cmtql1
   use general_tunnel_bound
   use general_tunnel_continuum
   use general_tunnel_asymptotic
   use general_tunnel_data
   use general_tunnel_dump
   use general_tunnel_nonadiabatic
   use general_tunnel_potential
   !
   write (out,"(t5,a)") trim(rcsid_accuracy)
   write (out,"(t5,a)") trim(rcsid_constants)
   write (out,"(t5,a)") trim(rcsid_derivative_tools)
   write (out,"(t5,a)") trim(rcsid_find_minimum)
   write (out,"(t5,a)") trim(rcsid_fftw)
   write (out,"(t5,a)") trim(rcsid_find_root)
   write (out,"(t5,a)") trim(rcsid_lapack)
   write (out,"(t5,a)") trim(rcsid_math)
   write (out,"(t5,a)") trim(rcsid_poly_tools)
   write (out,"(t5,a)") trim(rcsid_sort_tools)
   write (out,"(t5,a)") trim(rcsid_timer)
   write (out,"(t5,a)") trim(rcsid_tridiagonal_tools)
   write (out,"(t5,a)") trim(rcsid_tridiagonal_cmtql1)
   write (out,"(t5,a)") trim(rcsid_general_tunnel_bound)
   write (out,"(t5,a)") trim(rcsid_general_tunnel_continuum)
   write (out,"(t5,a)") trim(rcsid_general_tunnel_asymptotic)
   write (out,"(t5,a)") trim(rcsid_general_tunnel_data)
   write (out,"(t5,a)") trim(rcsid_general_tunnel_dump)
   write (out,"(t5,a)") trim(rcsid_general_tunnel_nonadiabatic)
   write (out,"(t5,a)") trim(rcsid_general_tunnel_potential)
 end subroutine versions

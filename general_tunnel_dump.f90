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
!  Saving intermediate quantities, mostly for debugging of general_tunnel.f90
!
module general_tunnel_dump
  use accuracy
  use constants
  use derivative_tools
  use find_root
  use general_tunnel_data
  use general_tunnel_potential
  use general_tunnel_asymptotic
  use lapack
  use math
  use derivative_tools
  use poly_tools
  use fftw
  use timer
  !$ use OMP_LIB
  implicit none
  public
  public rcsid_general_tunnel_dump
  !
  character(len=clen), save :: rcsid_general_tunnel_dump = "$Id: $"
  !
  contains
  !
  subroutine dump_bound_wavefunction(ichan,chan)
    integer(ik), intent(in)           :: ichan    ! Channel index, needed to create the file name
    type(channel), intent(in), target :: chan     ! Channel data
    !
    integer(ik)                 :: iu, eta_pt, xi_pt, max_order
    real(rk)                    :: eta
    character(len=clen)         :: file
    type(wavefunction), pointer :: st   ! Syntax sugar
    !
    if (file_bound==' ') return
    !
    write (file,fmt=trim(file_bound)) ichan
    !
    open (newunit=iu,form='formatted',recl=1024,action='write',position='rewind',status='replace',file=trim(file))
    write (iu,"('# ',a)") trim(comment)
    write (iu,"('#        Channel = ',i0)") ichan
    write (iu,"('#              M = ',i0)") mval
    write (iu,"('#           ZNUC = ',g36.26e3)") znuc
    write (iu,"('#      Potential = ',a)") trim(potential)
    write (iu,"('#             a1 = ',g36.26e3)") pot_param_real(1)
    write (iu,"('#             a2 = ',g36.26e3)") pot_param_real(2)
    write (iu,"('#             n1 = ',i0)") pot_param_int(1)
    dump_eta_points: do eta_pt=1,size(chan%bound)
      eta =  eta_tab(eta_pt)
      st  => chan%bound(eta_pt)
      write (iu,"('#         ETA PT = ',i0)") eta_pt
      write (iu,"('#            ETA = ', g36.26e3)") eta
      write (iu,"('#              E = ',2g36.26e3)") st%energy
      write (iu,"('#              Z = ',2g36.26e3)") st%zeta
      write (iu,"('#              F = ', g36.26e3)") st%efield
      write (iu,"('#  Normalization = ',a)") trim(st%norm)
      write (iu,"('#        ipt_max = ',i0)") st%ipt_max
      write (iu,"('#     R(ipt_max) = ',g36.26e3)") st%r(st%ipt_max)
      write (iu,"('#       ipt_stop = ',i0)") st%ipt_stop
      write (iu,"('#    R(ipt_stop) = ',g36.26e3)") st%r(st%ipt_stop)
      write (iu,"('# WARNING: Solution past ipt_stop is reset to zero after this data dump')")
      write (iu,"('#',a8,1x,a26,a8,1x,a26,5(1x,a26))") &
             ' I_ETA ', ' ETA ', ' I_XI ', ' XI ', ' Re[wf] ', ' Im[wf] ', ' Re[d wf/d R] ', ' Im[d wf/d R] ', '...'
      max_order = min(4,ubound(st%wf,dim=1))
      dump_xi_points: do xi_pt=1,size(st%r)
        write (iu,"(1x,i8,1x,g26.16e3,1x,i8,1x,g26.16e3,5(2x,g26.16e3,1x,g26.16e3))") &
               eta_pt, eta, xi_pt, st%r(xi_pt), st%wf(0:max_order,xi_pt)
      end do dump_xi_points
      write (iu,"()")
    end do dump_eta_points
    close (iu)
  end subroutine dump_bound_wavefunction
  !
  subroutine dump_continuum_wavefunction(ichan,st)
    integer(ik), intent(in)        :: ichan    ! Channel index, needed to create the file name
    type(wavefunction), intent(in) :: st       ! Channel data
    !
    integer(ik)                 :: iu, eta_pt, max_order
    real(rk)                    :: eta
    character(len=clen)         :: file
    !
    if (file_continuum==' ') return
    !
    write (file,fmt=trim(file_continuum)) ichan
    !
    open (newunit=iu,form='formatted',recl=1024,action='write',position='rewind',status='replace',file=trim(file))
    write (iu,"('# ',a)") trim(comment)
    write (iu,"('#        Channel = ',i0)") ichan
    write (iu,"('#              M = ',i0)") mval
    write (iu,"('#           ZNUC = ',g36.26e3)") znuc
    write (iu,"('#      Potential = ',a)") trim(potential)
    write (iu,"('#             a1 = ',g36.26e3)") pot_param_real(1)
    write (iu,"('#             a2 = ',g36.26e3)") pot_param_real(2)
    write (iu,"('#             n1 = ',i0)") pot_param_int(1)
    write (iu,"('#              E = ',2g36.26e3)") st%energy
    write (iu,"('#              F = ', g36.26e3)") st%efield
    write (iu,"('#  Normalization = ',a)") trim(st%norm)
    write (iu,"('#        ipt_max = ',i0)") st%ipt_max
    write (iu,"('#     R(ipt_max) = ',g36.26e3)") st%r(st%ipt_max)
    write (iu,"('#       ipt_stop = ',i0)") st%ipt_stop
    write (iu,"('#    R(ipt_stop) = ',g36.26e3)") st%r(st%ipt_stop)
    write (iu,"('#   Asymptotic Z = ',2g36.26e3)") st%zeta
    write (iu,"('#',a8,1x,a26,5(1x,a26))") &
           ' I_ETA ', ' ETA ', ' Re[wf] ', ' Im[wf] ', ' Re[d wf/d R] ', ' Im[d wf/d R] ', '...'
    dump_eta_points: do eta_pt=1,eta_npts
      eta =  eta_tab(eta_pt)
      max_order = min(4,ubound(st%wf,dim=1))
      write (iu,"(1x,i8,1x,g26.16e3,5(2x,g26.16e3,1x,g26.16e3))") &
             eta_pt, eta, st%wf(0:max_order,eta_pt)
    end do dump_eta_points
    close (iu)
  end subroutine dump_continuum_wavefunction
  !
  !  This routine will dump the total wavefunction.
  !
  subroutine dump_total_wavefunction(energy)
    complex(rk), intent(in)  :: energy
    !
    integer(ik)              :: xi_pt, eta_pt         ! Grid indices 
    integer(ik)              :: iu, alloc                 
    real(rk)                 :: eta, xi               ! Squared-parabolic coordinates
    real(rk)                 :: dv                    ! Volume element
    real(rk)                 :: x, z                  ! Cartesian coordinates
    complex(rk), allocatable :: wf(:,:,:)             ! Total wavefunction, and derivatives, 2D grid
                                                      ! Last index: 
                                                      !    0 = wavefunction
                                                      !    1 = xi gradient
                                                      !    2 = xi second derivative
                                                      !    3 = eta gradient
                                                      !    4 = eta second derivative
    !
    if (file_total==' ') return
    !
    call TimerStart('Dump total')
    allocate (wf(xi_npts,eta_npts,0:4),stat=alloc)
    !
    !  Build the quantity we want to visualize
    !
    select case (file_total_mode)
      case default
        write (out,"('general_tunnel_dump%dump_total_wavefunction: file_total_mode=',a,' is not recognized')") &
               trim(file_total_mode)
        stop 'general_tunnel_dump%dump_total_wavefunction - bad file_total_mode'
      case ('wavefunction')
        call build_wavefunction(continuum)
      case ('asymptotic outgoing')
        call build_asymptotic(outgoing=.true.)
      case ('asymptotic incoming')
        call build_asymptotic(outgoing=.false.)
      case ('fourier outgoing')
        call build_fft(outgoing=.true.)
      case ('fourier incoming')
        call build_fft(outgoing=.false.)
      case ('reverse total','reverse outgoing','reverse incoming')
        call build_wavefunction(pseudo)
    end select
    !
    open (newunit=iu,form='formatted',recl=1024,action='write',position='rewind',status='replace',file=trim(file_total))
    write (iu,"('# ',a)") trim(comment)
    write (iu,"('# Total wavefunction')")
    write (iu,"('#              M = ',i0)") mval
    write (iu,"('#           ZNUC = ',g36.26e3)") znuc
    write (iu,"('#      Potential = ',a)") trim(potential)
    write (iu,"('#             a1 = ',g36.26e3)") pot_param_real(1)
    write (iu,"('#             a2 = ',g36.26e3)") pot_param_real(2)
    write (iu,"('#             n1 = ',i0)") pot_param_int(1)
    write (iu,"('#              E = ',2g36.26e3)") energy
    write (iu,"('#              F = ', g36.26e3)") efield
    write (iu,"('#           Mode = ',a)") trim(file_total_mode)
    !
    write (iu,"(('#',2(1x,a5,1x,a26),5(1x,a26)))") &
           'I_XI', 'XI', 'I_ETA', 'ETA', ' dV ', ' X ', ' Z ', ' Re[psi] ', ' Im[psi] ', &
           '----', '--', '-----', '---', '----', '---', '---', '---------', '---------'
    !
    dump_eta: do eta_pt=1,eta_npts
      dump_xi: do xi_pt=1,xi_npts
        xi      = xi_tab(xi_pt)
        eta     = eta_tab(eta_pt)
        dv      = eta * xi * (eta**2 + xi**2)
        x       = xi * eta
        z       = 0.5_rk * (xi**2 - eta**2)
        write (iu,"(2(1x,i5,1x,g26.16e3),7(1x,g26.16e3))") xi_pt, xi, eta_pt, eta, dv, x, z, wf(xi_pt,eta_pt,0)
      end do dump_xi
      write (iu,"()")
    end do dump_eta
    !
    close (iu)
    !
    deallocate (wf)
    !
    call TimerStop('Dump total')
    !
    contains
    !
    subroutine build_wavefunction(continuum)
      type(wavefunction), intent(in) :: continuum(:)
      integer(ik) :: eta_pt, xi_pt, ichan
      !
      !  Calculate total wavefunction
      !
      wf(:,:,0) = 0
      assemble_channels: do ichan=1,n_channels
        !$omp parallel do collapse(2) default(none) &
        !$omp& shared(eta_npts,xi_npts,wf,continuum,channels,ichan) private(eta_pt,xi_pt)
        assemble_continuum: do eta_pt=1,eta_npts
          assemble_bound: do xi_pt=1,xi_npts
            wf(xi_pt,eta_pt,0) = wf(xi_pt,eta_pt,0) + continuum(ichan)%wf(0,eta_pt)*channels(ichan)%bound(eta_pt)%wf(0,xi_pt)
          end do assemble_bound
        end do assemble_continuum
        !$omp end parallel do
      end do assemble_channels
    end subroutine build_wavefunction
    !
    subroutine build_asymptotic(outgoing)
      logical, intent(in) :: outgoing
      integer(ik)         :: eta_pt, xi_pt, ichan
      real(rk)            :: eta, dir
      complex(rk)         :: wgt
      complex(rk)         :: cn(1:asymp_order)
      complex(rk)         :: as_wf(eta_npts,n_channels)
      !
      !  Calculate asymptotic wavefunctions in each channel. The code is 
      !  cribbed from general_tunnel_continuum.f90
      !
      if (outgoing) then
        dir = 1._rk
      else
        dir = -1._rk
      end if
      as_wf(1,:) = 0   ! The first point is at zero; the asymptotic solution is ill-defined there
      !$omp parallel do collapse(2) default(none) &
      !$omp& shared(n_channels,eta_npts,eta_tab,outgoing,continuum,dir,energy,efield,as_wf) &
      !$omp& private(ichan,eta_pt,eta,cn,wgt)
      asymptotic_channels: do ichan=1,n_channels
        evaluate_as_wf: do eta_pt=2,eta_npts
          eta = eta_tab(eta_pt)
          if (     outgoing) wgt = continuum(ichan)%c_out
          if (.not.outgoing) wgt = continuum(ichan)%c_in
          call asymptotic_solution(efield,energy,continuum(ichan)%zeta,dir,eta,cn)
          as_wf(eta_pt,ichan) = wgt * evaluate_asymptote(efield,energy,dir,cn,eta) / sqrt(eta)
        end do evaluate_as_wf
      end do asymptotic_channels
      !$omp end parallel do
      !
      wf(:,:,0) = 0
      assemble_channels: do ichan=1,n_channels
        !$omp parallel do collapse(2) default(none) &
        !$omp&  shared(eta_npts,xi_npts,wf,as_wf,channels,ichan) private(eta_pt,xi_pt)
        assemble_continuum: do eta_pt=1,eta_npts
          assemble_bound: do xi_pt=1,xi_npts
            wf(xi_pt,eta_pt,0) = wf(xi_pt,eta_pt,0) + as_wf(eta_pt,ichan)*channels(ichan)%bound(eta_pt)%wf(0,xi_pt)
          end do assemble_bound
        end do assemble_continuum
        !$omp end parallel do
      end do assemble_channels
    end subroutine build_asymptotic
    !
    subroutine build_fft(outgoing)
      logical, intent(in) :: outgoing
      integer(ik)         :: eta_pt, xi_pt, ichan, izero
      complex(drk)        :: proc_wf(eta_npts,n_channels)
      real(rk)            :: pos, width, eta
      !
      !  Prepare for applying Fourier filter to the data
      !
      pos   = fourier_centre
      width = fourier_width
      if (pos<0) pos = eta_max / 2._rk
      if (width<0) width = 2._rk * min(pos,eta_max-pos)
      !
      !  Process the continuum wavefunction, picking out only the "outgoing"/"incoming"
      !  component according to FFT. Please keep in mind that these are NOT true incoming/outgoing solutions!
      !  We'll pad the data with zeros on the left and the right.
      !
      fill_outgoing: do ichan=1,n_channels
        ! Type conversion is necessary: FFTW is not available in quad precision, so we'll have to
        ! use doubles here even if we are doing everything else in quad precision.
        proc_wf(:,ichan) = cmplx(continuum(ichan)%wf(0,:),kind=drk)
        apply_fft_filter: do eta_pt=1,eta_npts
          eta = eta_tab(eta_pt)
          proc_wf(eta_pt,ichan) = cmplx(proc_wf(eta_pt,ichan) * blackman_harris(eta,pos,width),kind=kind(proc_wf))
        end do apply_fft_filter
      end do fill_outgoing
      !
      call fftw_1d(eta_npts,n_channels,proc_wf,invert=.false.)
      izero = 1 + (eta_npts-1)/2  ! Position of the zero-frequency component.
      proc_wf(izero,:) = cmplx(0.5_rk * proc_wf(izero,:),kind=kind(proc_wf))
      if (outgoing) then
        ! Keep only the positive-frequency components
        proc_wf(:izero-1,:) = 0
      else
        ! Keep only the negative-frequency components
        proc_wf(izero+1:,:) = 0
      end if
      call fftw_1d(eta_npts,n_channels,proc_wf,invert=.true.)
      ! Restore wavefunction normalization
      proc_wf(:,:) = proc_wf(:,:)/eta_npts
      !
      wf(:,:,0) = 0
      assemble_channels: do ichan=1,n_channels
        !$omp parallel do collapse(2) default(none) &
        !$omp&  shared(eta_npts,xi_npts,wf,proc_wf,channels,ichan) private(eta_pt,xi_pt)
        assemble_continuum: do eta_pt=1,eta_npts
          assemble_bound: do xi_pt=1,xi_npts
            wf(xi_pt,eta_pt,0) = wf(xi_pt,eta_pt,0) + proc_wf(eta_pt,ichan)*channels(ichan)%bound(eta_pt)%wf(0,xi_pt)
          end do assemble_bound
        end do assemble_continuum
        !$omp end parallel do
      end do assemble_channels
    end subroutine build_fft
    !
    !  Blackman-Harris Fourier window function
    !
    function blackman_harris(x,pos,width) result(wgt)
      real(rk), intent(in) :: x, pos, width 
      real(rk)             :: wgt
      real(rk)             :: a, dx
      !
      a  = 2._rk * pi / width
      dx = x-pos
      if (abs(dx)>0.5_rk*width) then
        wgt = 0._rk
      else
        wgt = 0.35875_rk + 0.48829_rk*cos(a*dx) + 0.14128_rk*cos(2*a*dx) + 0.01168_rk*cos(3*a*dx)
      end if
    end function blackman_harris
    !
  end subroutine dump_total_wavefunction
  !
  subroutine dump_channel_parameter(tag,val)
    character(len=*), intent(in) :: tag      ! Channel tag; needed to create the filename
    complex(rk), intent(in)      :: val(:,:) ! Data field to dump. First index is the value and its derivatives (x n!)
                                             !                     Second index is the eta point
    !
    integer(ik)                 :: iu, eta_pt, max_order
    character(len=clen)         :: file
    !
    if (file_coupling==' ') return
    !
    write (file,fmt=trim(file_coupling)) trim(tag)
    !
    open (newunit=iu,form='formatted',recl=2048,action='write',position='rewind',status='replace',file=trim(file))
    write (iu,"('# ',a)") trim(comment)
    write (iu,"('# ',a)") trim(tag)
    write (iu,"('#',a8,1x,a26,5(2x,a26,1x,a26))") &
           ' I_ETA ', ' ETA ', ' Re[F] ', ' Im[F] ', ' Re[d F/d R] ', ' Im[d F/d R] ', '...'
    max_order = min(20,size(val,dim=1))
    dump_eta_points: do eta_pt=1,eta_npts
      write (iu,"(1x,i8,1x,g26.16e3,20(2x,g26.16e3,1x,g26.16e3))") &
             eta_pt, eta_tab(eta_pt), val(:max_order,eta_pt)
    end do dump_eta_points
    close (iu)
  end subroutine dump_channel_parameter
  !
  !  This routine will dump the total wavefunction, interpolated on a uniform Cartesian grid.
  !
  subroutine dump_total_cartesian_wavefunction(energy)
    complex(rk), intent(in)  :: energy
    !
    real(rk), allocatable    :: xyz  (:,:,:,:)  ! Local cartesian coordinates of the grid points: X, Y, Z
    real(rk), allocatable    :: pa2  (:,:,:,:)  ! Squared-parabolic coodrinates for the grid points: ETA, XI, PHI
    real(rk), allocatable    :: pot    (:,:,:)  ! Potential at the grid points
    complex(rk), allocatable :: wfn    (:,:,:)  ! Wavefunction at the grid points
    complex(rk), allocatable :: lap    (:,:,:)  ! Laplacian of the wavefunction, calculated using the five-point numerical
                                                ! differences. 
    complex(rk), allocatable :: res    (:,:,:)  ! Local residue of the SE solution
    real(rk), allocatable    :: bhm_r  (:,:,:)  ! Modulus of the wavefunction, needed for the Bohmian trajectory
    real(rk), allocatable    :: bhm_s  (:,:,:)  ! Phase of the wavefunction, ditto
    real(rk), allocatable    :: bhm_u  (:,:,:)  ! Bohmian quantum potential
    real(rk), allocatable    :: bhm_v(:,:,:,:)  ! Bohmian velocity field
    integer(ik)              :: x0, y0, z0      ! Minimal indices of the local Cartesian grid
    integer(ik)              :: xn, yn, zn      ! Maximal indices of the local Cartesian grid
    integer(ik)              :: ix, iy, iz      ! Running Cartesian indices
    integer(ik)              :: alloc                 
    !
    if (file_cartesian==' ' .and. file_bohm==' ' .and. file_husimi==' ') return
    !
    call TimerStart('Dump Cartesian')
    !
    x0 = cartesian_npts(1,1) ; y0 = cartesian_npts(1,2) ; z0 = cartesian_npts(1,3) ;
    xn = cartesian_npts(2,1) ; yn = cartesian_npts(2,2) ; zn = cartesian_npts(2,3) ;
    !
    !  We'll evaluate one extra point on each side, to make sure we have the Laplacian
    !
    allocate (xyz  (3,x0-3:xn+3,y0-3:yn+3,z0-3:zn+3), &
              pa2  (3,x0-3:xn+3,y0-3:yn+3,z0-3:zn+3), &
              pot  (  x0-3:xn+3,y0-3:yn+3,z0-3:zn+3), &
              wfn  (  x0-3:xn+3,y0-3:yn+3,z0-3:zn+3), &
              lap  (  x0-3:xn+3,y0-3:yn+3,z0-3:zn+3), &
              res  (  x0-3:xn+3,y0-3:yn+3,z0-3:zn+3), &
              bhm_r(  x0-3:xn+3,y0-3:yn+3,z0-3:zn+3), &
              bhm_s(  x0-3:xn+3,y0-3:yn+3,z0-3:zn+3), &
              bhm_u(  x0-3:xn+3,y0-3:yn+3,z0-3:zn+3), &
              bhm_v(3,x0-3:xn+3,y0-3:yn+3,z0-3:zn+3), &
              stat=alloc)
    if (alloc/=0) then
      write (out,"('dump_total_cartesian_wavefunction: Error ',i0,' allocating temporaries')") alloc
      stop 'dump_total_cartesian_wavefunction - alloc'
    end if
    !
    !  Fill coordinates of all grid points and potential at these points
    !
    call build_coordinates
    !
    !  Build the quantity we want to visualize
    !
    select case (file_total_mode)
      case default
        write (out,"('general_tunnel_dump%dump_total_cartesian_wavefunction: file_total_mode=',a,' is not recognized')") &
               trim(file_total_mode)
        stop 'general_tunnel_dump%dump_total_wavefunction - bad file_total_mode'
      case ('wavefunction')
        call build_wavefunction(continuum,use_eta_origin=.true.)
      case ('asymptotic outgoing','asymptotic incoming','fourier outgoing','fourier incoming')
        write (out,"('WARNING: file_total_mode=',a,' is not implemented in dump_total_cartesian_wavefunction')") &
               trim(file_total_mode)
        write (out,"('WARNING: Continuing with the total wavefunction instead.')")
        call build_wavefunction(continuum,use_eta_origin=.true.)
      case ('reverse total','reverse outgoing','reverse incoming')
        call build_wavefunction(pseudo,use_eta_origin=.false.)
    end select
    !
    if (file_cartesian/=' ') then
      call build_laplacian_and_residue
      call cartesian_output
    end if
    !
    if (file_bohm/=' ') then
      call build_bohm_potential_and_velocity
      call bohmian_output
    end if
    !
    if (file_husimi/=' ') then
      call husimi_compute_and_output ! Husimi distribution is potentially too large to store completely!
    end if
    !
    deallocate (xyz,pa2,pot,wfn,lap,res,bhm_r,bhm_s,bhm_u,bhm_v)
    !
    call TimerStop('Dump Cartesian')
    !
    contains
    !
    subroutine build_coordinates
      !$omp parallel do collapse(3) default(none) private(ix,iy,iz) &
      !$omp& shared(x0,xn,y0,yn,z0,zn,xyz,pa2,pot,cartesian_dx,cartesian_ref)
      fill_coord_z: do iz=z0-3,zn+3
        fill_coord_y: do iy=y0-3,yn+3
          fill_coord_x: do ix=x0-3,xn+3
            xyz(:,ix,iy,iz) = cartesian_dx * (cartesian_ref(:) + real( (/ix,iy,iz/), kind=rk))
            pa2(:,ix,iy,iz) = squared_parabolic(xyz(:,ix,iy,iz))
            pot(  ix,iy,iz) = total_potential(xyz(:,ix,iy,iz))
          end do fill_coord_x
        end do fill_coord_y
      end do fill_coord_z
      !$omp end parallel do
    end subroutine build_coordinates
    !
    !  Transformation from local Cartesian to squared-parabolic coordinates.
    !
    function squared_parabolic(xyz) result(pa2)
      real(rk), intent(in) :: xyz(:)     ! X, Y, Z
      real(rk)             :: pa2(3)     ! Eta, Xi, Phi
      real(rk)             :: r
      !
      r = sqrt(sum(xyz**2))
      !
           if (xyz(3)<-r) then   ! Could just happen due to rounding errors
        pa2(1) = sqrt(2._rk*r)   ! Eta
        pa2(2) = 0._rk           ! Xi
      else if (xyz(3)> r) then   ! Ditto
        pa2(1) = 0._rk           ! Eta
        pa2(2) = sqrt(2._rk*r)   ! Xi
      else                       ! General case
        pa2(1) = sqrt(r-xyz(3))  ! Eta
        pa2(2) = sqrt(r+xyz(3))  ! Xi
      end if
      pa2(3) = atan2(xyz(2),xyz(1)) + cartesian_phi ! Phi
    end function squared_parabolic
    !
    !  Total potential, including the long-range part and static external field
    !
    function total_potential(xyz) result(pot)
      real(rk), intent(in) :: xyz(:)     ! X, Y, Z
      real(rk)             :: pot
      !
      real(rk)             :: r
      real(rk)             :: utab(1)
      real(rk)             :: pa2(3)
      !
      r = sqrt(sum(xyz**2))
      if (abs(r)<=spacing(100._rk)) then
        pot = 0._rk  ! Total potential is singular at the origin
      else
        ! call potential_u_r(r,0_ik,utab)
        pa2 = squared_parabolic(xyz)
        call potential_u_halfx2(pa2(1),pa2(2),0_ik,utab)
        pot = -(znuc+utab(1))/r + efield*xyz(3)
      end if
    end function total_potential
    !
    subroutine build_wavefunction(continuum,use_eta_origin)
      type(wavefunction), intent(in) :: continuum(:)
      logical, intent(in)            :: use_eta_origin
      !
      !$omp parallel do collapse(3) default(none) private(ix,iy,iz) &
      !$omp& shared(x0,y0,z0,xn,yn,zn,continuum,pa2,wfn,bhm_r,bhm_s) &
      !$omp& shared(use_eta_origin,cart_laplace_order)
      build_z: do iz=z0-cart_laplace_order,zn+cart_laplace_order
        build_y: do iy=y0-cart_laplace_order,yn+cart_laplace_order
          build_x: do ix=x0-cart_laplace_order,xn+cart_laplace_order
            wfn  (ix,iy,iz) = interpolate_wavefunction(continuum,pa2(:,ix,iy,iz),use_eta_origin)
            bhm_r(ix,iy,iz) = abs(wfn(ix,iy,iz))
            bhm_s(ix,iy,iz) = atan2(aimag(wfn(ix,iy,iz)),real(wfn(ix,iy,iz),kind=rk))
          end do build_x
        end do build_y
      end do build_z
      !$omp end parallel do
    end subroutine build_wavefunction
    !
    subroutine build_laplacian_and_residue
      !$omp parallel do collapse(3) default(none) private(ix,iy,iz) &
      !$omp& shared(x0,xn,y0,yn,z0,zn,cartesian_dx,cart_laplace_order,energy,wfn,lap,pot,res)
      build_laplacian_z: do iz=z0,zn
        build_laplacian_y: do iy=y0,yn
          build_laplacian_x: do ix=x0,xn
            select case (cart_laplace_order)
              case default ; stop 'general_tunnel_dump%build_laplacian_and_residue - bad cart_laplace_order'
              case (1) ; lap(ix,iy,iz) = dt_lap3(cartesian_dx,wfn(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1))
              case (2) ; lap(ix,iy,iz) = dt_lap5(cartesian_dx,wfn(ix-2:ix+2,iy-2:iy+2,iz-2:iz+2))
              case (3) ; lap(ix,iy,iz) = dt_lap7(cartesian_dx,wfn(ix-3:ix+3,iy-3:iy+3,iz-3:iz+3))
            end select
            res(ix,iy,iz) = -0.5_rk * lap(ix,iy,iz) + (pot(ix,iy,iz)-energy) * wfn(ix,iy,iz)
          end do build_laplacian_x
        end do build_laplacian_y
      end do build_laplacian_z
      !$omp end parallel do
    end subroutine build_laplacian_and_residue
    !
    !  Bohm's quantum potential is defined as (-1/2) (\Delta R) / R, where R is the density
    !  Bohm's velocity is defined as \Grad S, where S is the phase of the wavefunction
    !  There is a slight complication with the velocity: it needs to be unwrapped before
    !  we differentiate it.
    !
    subroutine build_bohm_potential_and_velocity
      real(rk) :: unwrap(-3:3,3)
      !$omp parallel do collapse(3) default(none) private(ix,iy,iz,unwrap) &
      !$omp& shared(x0,xn,y0,yn,z0,zn,cartesian_dx,cart_laplace_order,bhm_r,bhm_s,bhm_u,bhm_v)
      build_laplacian_z: do iz=z0,zn
        build_laplacian_y: do iy=y0,yn
          build_laplacian_x: do ix=x0,xn
            unwrap(-3:3,1) = bhm_s(ix-3:ix+3,    iy   ,    iz   )
            unwrap(-3:3,2) = bhm_s(    ix   ,iy-3:iy+3,    iz   )
            unwrap(-3:3,3) = bhm_s(    ix   ,    iy   ,iz-3:iz+3)
            call poly_phase_unwrap(unwrap(:,1))
            call poly_phase_unwrap(unwrap(:,2))
            call poly_phase_unwrap(unwrap(:,3))
            select case (cart_laplace_order)
              case default ; stop 'general_tunnel_dump%build_bohm_potential_and_velocity - bad cart_laplace_order'
              case (1) 
                bhm_u(  ix,iy,iz) = -0.5_rk * dt_lap3 (cartesian_dx,bhm_r(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1))
                bhm_v(1,ix,iy,iz) = dt_grad3(cartesian_dx,unwrap(-1:1,1))
                bhm_v(2,ix,iy,iz) = dt_grad3(cartesian_dx,unwrap(-1:1,2))
                bhm_v(3,ix,iy,iz) = dt_grad3(cartesian_dx,unwrap(-1:1,3))
              case (2) 
                bhm_u(  ix,iy,iz) = -0.5_rk * dt_lap5 (cartesian_dx,bhm_r(ix-2:ix+2,iy-2:iy+2,iz-2:iz+2))
                bhm_v(1,ix,iy,iz) = dt_grad5(cartesian_dx,unwrap(-2:2,1))
                bhm_v(2,ix,iy,iz) = dt_grad5(cartesian_dx,unwrap(-2:2,2))
                bhm_v(3,ix,iy,iz) = dt_grad5(cartesian_dx,unwrap(-2:2,3))
              case (3) 
                bhm_u(  ix,iy,iz) = -0.5_rk * dt_lap7 (cartesian_dx,bhm_r(ix-3:ix+3,iy-3:iy+3,iz-3:iz+3))
                bhm_v(1,ix,iy,iz) = dt_grad7(cartesian_dx,unwrap(-3:3,1))
                bhm_v(2,ix,iy,iz) = dt_grad7(cartesian_dx,unwrap(-3:3,2))
                bhm_v(3,ix,iy,iz) = dt_grad7(cartesian_dx,unwrap(-3:3,3))
            end select
            ! bhm_r is guaranteed to be non-negative, but may be zero. Watch out!
            if (bhm_r(ix,iy,iz)>spacing(bhm_u(ix,iy,iz))) then
              bhm_u(ix,iy,iz) = bhm_u(ix,iy,iz) / bhm_r(ix,iy,iz)
            else
              bhm_u(ix,iy,iz) = 0._rk
            end if
          end do build_laplacian_x
        end do build_laplacian_y
      end do build_laplacian_z
      !$omp end parallel do
    end subroutine build_bohm_potential_and_velocity
    !
    function interpolate_wavefunction(continuum,pa2,use_eta_origin) result(wfn)
      type(wavefunction), intent(in) :: continuum(:)
      real(rk), intent(in)           :: pa2(:)          ! Eta, Xi, Phi
      logical, intent(in)            :: use_eta_origin  ! Series expansion of the continuum (eta) solution at the origin is valid.
      complex(rk)                    :: wfn
      !
      integer(ik) :: ichan              ! Channel number
      integer(ik) :: iord               ! Polynomial order of the piece-wise function
      integer(ik) :: eta_min            ! Permissible lower bondary of the eta (continuum) grid
      integer(ik) :: eta_pt             ! Reference grid point in eta
      integer(ik) :: xi_pt              ! Reference grid point in xi
      integer(ik) :: eta_iter           ! Eta-position iterator
      integer(ik) :: eta_low            ! Lowest eta point for the interpolation window
      real(rk)    :: eta, xi, phi       ! Desired eta, xi, and phi values
      real(rk)    :: d_eta              ! Displacement from the eta reference point, fraction of the grid spacing
      real(rk)    :: d_xi               ! ditto, xi grid
      complex(rk) :: wfn_xi (2)         ! Bound part of the wavefunction and its gradient
      complex(rk) :: wfn_eta(2)         ! Continuum part of the wavefunction and its gradient
      complex(rk) :: wfn_xi_tab(cart_interp_points) 
                                        ! Bound (xi) solutions for interpolation
      !
      eta_min = 2
      if (use_eta_origin) eta_min = 1
      !
      !  Our first task is to find the closest grid position in both eta and xi.
      !  Luckily, both grids are uniform, so this is only a question of scaling.
      !
      eta    = pa2(1)
      xi     = pa2(2)
      phi    = pa2(3)
      eta_pt = min(max(eta_min,1_ik+nint(eta/eta_tab(2),kind=ik)),eta_npts)
      xi_pt  = min(max(   1_ik,1_ik+nint( xi/ xi_tab(2),kind=ik)), xi_npts)
      d_eta  = (eta-eta_tab(eta_pt))/eta_tab(2)
      d_xi   = ( xi- xi_tab( xi_pt))/ xi_tab(2)
      !
      wfn = 0
      if (abs(d_eta)<=1._rk .and. abs(d_xi)<=1._rk) then
        !
        !  We won't even try for the points where extrapolation is requested -
        !  our piece-wise solution is not well-defined there!
        !
        assemble_channels: do ichan=1,n_channels
          !
          !  Bound part first (wavefunctuion along xi). This is quite a bit of work,
          !  since we do not have a solution at the desired eta point. Instead, we will
          !  need to extrapolate xi solutions at neighbouring eta grid points, then 
          !  interpolate them at the desired eta point.
          !
          !  For bound interpolation, we are always allowed to use the solution
          !  at the eta origin, regardless of the value of use_eta_origin.
          !
          eta_low = min(max(1_ik,eta_pt-cart_interp_points/2),eta_npts-cart_interp_points+1)
          if (eta_low<1_ik) stop 'dump_total_cartesian_wavefunction%interpolate_wavefunction - can''t interpolate!'
          !
          prepare_xi_tab: do eta_iter=1,cart_interp_points
            iord =                 channels(ichan)%bound(eta_low+eta_iter-1)%iord(    xi_pt)
            call poly_power_series(channels(ichan)%bound(eta_low+eta_iter-1)%wf(:iord,xi_pt),d_xi,wfn_xi)
            wfn_xi_tab(eta_iter) = wfn_xi(1)
          end do prepare_xi_tab
          wfn_xi(1) = MathInterpolate(eta,eta_tab(eta_low:eta_low+cart_interp_points-1),wfn_xi_tab)
          !
          !  Continuum part (wavefunction along eta)
          !
          iord = continuum(ichan)%iord(eta_pt)
          call poly_power_series(continuum(ichan)%wf(:iord,eta_pt),d_eta,wfn_eta)
          !
          !  Colect the result
          !
          wfn = wfn + wfn_xi(1)*wfn_eta(1)
        end do assemble_channels
        wfn = wfn * exp(cmplx(0_ik,mval,kind=rk)*phi)
      end if
    end function interpolate_wavefunction
    !
    subroutine cartesian_output
      integer(ik) :: iu
      !
      open (newunit=iu,form='formatted',recl=1024,action='write',position='rewind',status='replace',file=trim(file_cartesian))
      call write_common_header(iu)
      !
      write (iu,"(('#',3(1x,a8),1x,13(1x,a26)))") &
             '  1 ', '  2 ', '  3 ', ' 4 ', ' 5 ', ' 6 ', ' 7 ', '    8    ', '    9    ', &
             '   10    ', '   11    ', '   12    ', '   13    ', '  10 ', ' 11 ', '  12 ', &
             ' IX ', ' IY ', ' IZ ', ' X ', ' Y ', ' Z ', ' V ', ' Re[Psi] ', ' Im[Psi] ', &
             ' Re[Res] ', ' Im[Res] ', ' Re[Lap] ', ' Im[Lap] ', ' ETA ', ' XI ', ' PHI ', &
             '----', '----', '----', '---', '---', '---', '---', '---------', '---------', &
             '---------', '---------', '---------', '---------', '-----', '----', '-----'
      !
      dump_z: do iz=z0,zn
        dump_x: do ix=x0,xn
          dump_y: do iy=y0,yn
            write (iu,"(1x,3(1x,i8),1x,13(1x,g26.16e3))") &
                   ix, iy, iz, xyz(1:3,ix,iy,iz), pot(ix,iy,iz), wfn(ix,iy,iz), res(ix,iy,iz), lap(ix,iy,iz), pa2(1:3,ix,iy,iz)
          end do dump_y
          if (yn-y0>1) write (iu,"()")
        end do dump_x
        if (xn-x0>1) write (iu,"()")
      end do dump_z
      !
      close (iu)
    end subroutine cartesian_output
    !
    subroutine bohmian_output
      integer(ik) :: iu
      !
      open (newunit=iu,form='formatted',recl=1024,action='write',position='rewind',status='replace',file=trim(file_bohm))
      call write_common_header(iu)
      !
      write (iu,"(('#',3(1x,a8),1x,10(1x,a26)))") &
             '  1 ', '  2 ', '  3 ', ' 4 ', ' 5 ', ' 6 ', ' 7 ', '    8    ', '    9    ', & 
             '   10    ', '   11    ', '   12    ', '   13    ', &
             ' IX ', ' IY ', ' IZ ', ' X ', ' Y ', ' Z ', ' V ', ' Re[Psi] ', ' Im[Psi] ', &
             ' UQuant. ', '   VX    ', '   VY    ', '   VZ    ', &
             '----', '----', '----', '---', '---', '---', '---', '---------', '---------', &
             '---------', '---------', '---------', '---------'
      !
      dump_z: do iz=z0,zn
        dump_x: do ix=x0,xn
          dump_y: do iy=y0,yn
            write (iu,"(1x,3(1x,i8),1x,10(1x,g26.16e3))") &
                   ix, iy, iz, xyz(1:3,ix,iy,iz), pot(ix,iy,iz), wfn(ix,iy,iz), bhm_u(ix,iy,iz), bhm_v(1:3,ix,iy,iz)
          end do dump_y
          if (yn-y0>1) write (iu,"()")
        end do dump_x
        if (xn-x0>1) write (iu,"()")
      end do dump_z
      !
      close (iu)
    end subroutine bohmian_output
    !
    subroutine husimi_compute_and_output
      integer(ik) :: iu
      !
      open (newunit=iu,form='formatted',recl=1024,action='write',position='rewind',status='replace',file=trim(file_husimi))
      call write_common_header(iu)
      write (iu,"('# Husimi transform dimensionality = ',i24)") husimi_ndim
      write (iu,"('#    Husimi transform coordinates = ',3i8)") husimi_coord
      write (iu,"('#              Husimi filter FWHM = ',g24.12)") husimi_width
      select case (husimi_ndim)
        case default
          stop 'general_tunnel_dump%husimi_compute_and_output - bad husimi_ndim value'
        case (1)
          select case (husimi_coord(1))
            case default ; stop 'general_tunnel_dump%husimi_compute_and_output - bad husimi_coord(1)'
            case (3) ; call husimi_1D_Z(iu)
          end select
      end select
      close (iu)
    end subroutine husimi_compute_and_output
    !
    !  1D Husimi transform along the Z direction. First, we need a wrapper, to handle
    !  various permutations of the index positions.
    !
    subroutine husimi_1D_Z(iu)
      integer(ik), intent(in)  :: iu
      real(rk), allocatable    :: ptab(:)  ! Momentum grid corresponding to the coordinate grid being FFTed
      real(rk), allocatable    :: pmax(:)  ! Position of the maximum of hus(:,iz), interpolated from grid values
      real(rk), allocatable    :: pexp(:)  ! Expectation of P at each IZ
      complex(rk), allocatable :: sec(:)   ! Section of the wavefunction along the Z direction
      real(rk), allocatable    :: hus(:,:) ! Husimi distribution derived from the section. First index is P, second index is Z
      integer(ik)              :: ip       ! Index for the P coordinate
      !
      allocate (sec(z0:zn),ptab(z0:zn),pmax(z0:zn),pexp(z0:zn),hus(z0:zn,z0:zn),stat=alloc)
      if (alloc/=0) stop 'general_tunnel_dump%husimi_1D_Z - no memory'
      !
      write (iu,"(('#@#',3(1x,a8),10x,5(1x,a26)))") &
           '  1 ', '  2 ', '  3 ', ' 4 ', ' 5 ', ' 6 ', '   7  ', '   8  ', &
           ' IX ', ' IY ', ' IZ ', ' X ', ' Y ', ' Z ', ' Pmax ', ' Pexp ', &
           '----', '----', '----', '---', '---', '---', '------', '------'
      if (husimi_detail) &
      write (iu,"(('#  ',4(1x,a8),1x,5(1x,a26)))") &
           '  1 ', '  2 ', '  3 ', '  4 ', ' 5 ', ' 6 ', ' 7 ', ' 8 ', '      9      ', &
           ' IX ', ' IY ', ' IZ ', ' IP ', ' X ', ' Y ', ' Z ', ' P ', ' Pr(X,Y,Z,P) ', &
           '----', '----', '----', '----', '---', '---', '---', '---', '-------------'
      hc1z_x: do ix=x0,xn
        hc1z_y: do iy=y0,yn
          sec(z0:zn) = wfn(ix,iy,z0:zn)
          ! Formally, the call below was initially wrong, with xyz(3,1,1,:) as the third
          ! argument. However, because the Cartesian grid is uniform and only coordinate
          ! differences matter in husimi_1D, the result was nonetheless correct.
          call husimi_1D(z0,zn,xyz(3,1,1,z0:zn),sec,ptab,hus,pmax,pexp)
          report_z: do iz=z0,zn
            write (iu,"('#@ ',3(1x,i8),10x,5(1x,g26.16e3))") ix, iy, iz, xyz(1:3,ix,iy,iz), pmax(iz), pexp(iz)
            if (.not.husimi_detail) cycle report_z
            report_pz: do ip=z0,zn
              write (iu,"(3x,4(1x,i8),1x,5(1x,g26.16e3))") ix, iy, iz, ip, xyz(1:3,ix,iy,iz), ptab(ip), hus(ip,iz)
            end do report_pz
            write (iu,"()")
          end do report_z
        end do hc1z_y
      end do hc1z_x
      deallocate (sec,ptab,pmax,pexp,hus)
    end subroutine husimi_1D_Z
    !
    subroutine husimi_1D(c0,cn,cval,psi,ptab,hus,pmax,pexp)
      integer(ik), intent(in)   :: c0                ! Lowest coordinate index
      integer(ik), intent(in)   :: cn                ! Highest coordinate index
      real(rk), intent(in)      :: cval(c0:cn)       ! Spatial coordinates
      complex(rk), intent(in)   :: psi (c0:cn)       ! Wavefunction to analyze
      real(rk), intent(out)     :: ptab(c0:cn)       ! Output momentum grid
      real(rk), intent(out)     :: hus (c0:cn,c0:cn) ! Output Husimi distribution (probability of finding electron with
                                                     ! momentum P at coordinate X). The first index is momentum; the
                                                     ! second index is coordinate
      real(rk), intent(out)     :: pmax(c0:cn)       ! For each spatial coordinate, the momentum P where probability is
                                                     ! at maximum.
      real(rk), intent(out)     :: pexp(c0:cn)       ! For each spatial coordinate, expectation of the momentum P
      !                        
      integer(ik)               :: csz               ! Size of the coordinate buffer
      integer(ik)               :: ic, ip            ! Iterators for the coordinate and momentum
      integer(ik)               :: ip0               ! Position of the zero frequence
      real(rk)                  :: dp                ! Resolution of the P grid
      real(rk)                  :: alp               ! Exponent of the Gaussian filter
      real(rk), allocatable     :: gauss(:)          ! Coordinate mask used in the Husimi distribution (shared)
      complex(drk), allocatable :: psim (:,:)        ! Masked wavefunction (per-thread)
      !
      csz = cn - c0 + 1
      allocate (gauss(-csz:csz),psim(c0:cn,c0:cn),stat=alloc)
      if (alloc/=0) stop 'general_tunnel_dump%husimi_1D - allocate 1'
      !
      dp  = 2._rk*pi / (cval(cn)-cval(c0))
      ip0 = c0 + (csz-1_ik)/2_ik
      fill_pgrid: do ip=c0,cn
        ptab(ip) = dp*(ip-ip0)
      end do fill_pgrid
      !
      !  Our filter is a Gaussian, with FWHM of husimi_width, normalized to 1:
      !    sqrt(alp/pi) Exp(-alp*x**2), with alp=4ln(2)/w0**2
      !  We'll fold the normalization factor for the FFT into the filter as well,
      !  adding an extra factor of
      !    cartesian_dx/sqrt(2*pi)
      !
      alp = 4._rk*log(2._rk)/husimi_width**2
      fill_gaussian: do ic=-csz,csz
        gauss(ic) = (cartesian_dx/sqrt(2._rk*pi)) * sqrt(alp/pi) * exp(-alp*(ic*cartesian_dx)**2)
      end do fill_gaussian
      !
      !$omp parallel default(none) num_threads(cn-c0+1) &
      !$omp& private(ic) &
      !$omp& shared(c0,cn,psi,gauss,csz,hus,pmax,pexp,ptab,psim)
        !$omp do
        fill_filter_position: do ic=c0,cn
          psim(c0:cn,ic) = cmplx(psi(c0:cn) * gauss(c0-ic:cn-ic),kind=kind(psim))
        end do fill_filter_position
        !$omp end do
        !
        !  For some reason, the code below leads to a segmentation violation in FFTW library.
        !  As the result, we are forced to perform FFT in a single thread. Oops.
        !
        !!!!!$omp do
        !!!!fft_filtered_position: do ic=c0,cn
        !!!!  call fftw_1d(csz,csz,psim(c0:cn,ic:ic),invert=.false.)
        !!!!end do fft_filtered_position
        !!!!!$omp end do
        !
        !$omp single
          call fftw_1d(csz,csz,psim(c0:cn,c0:cn),invert=.false.)
        !$omp end single
        !$omp barrier
        !
        !$omp do
        save_filter_position: do ic=c0,cn
          hus(:,ic) = abs(psim(:,ic))**2
          hus(:,ic) = hus(:,ic) / maxval(hus(:,ic),dim=1)
          pmax(ic)  = interpolate_max_position(ptab,hus(:,ic))
          pexp(ic)  = sum(ptab * hus(:,ic)**2) / sum(hus(:,ic)**2)
        end do save_filter_position
        !$omp end do
      !$omp end parallel
      deallocate (gauss,psim)
    end subroutine husimi_1D
    !
    function interpolate_max_position(x,y) result(xmax)
      real(rk), intent(in) :: x(:)   ! Grid positions
      real(rk), intent(in) :: y(:)   ! Function values
      real(rk)             :: xmax   ! Position where (interpolated) function reaches the maximum
      !                              
      integer(ik)          :: loc    ! Index of central point we use for 3-point interpolation
      real(rk)             :: yp, ym ! Increase in the function at the "right" and "left" points
      real(rk)             :: xp, xm ! Displacement of the "right" and "left" points
      real(rk)             :: dx     ! Position where the derivative of the quadratic interpolant vanishes
      real(rk)             :: a, b   ! Numerator and denominator
      !
      if (size(x)/=size(y)) stop 'eneral_tunnel_dump%interpolate_max_position - bad array'
      loc = maxloc(y,dim=1)
      if (loc<=1) loc = 2
      if (loc>=size(y)) loc = size(y)-1
      !
      xp = x(loc+1)-x(loc)
      xm = x(loc-1)-x(loc)
      yp = y(loc+1)-y(loc)
      ym = y(loc-1)-y(loc)
      !
      a  = yp * xm**2 - ym * xp**2
      b  = yp * xm    - ym * xp
      if (abs(b)>spacing(a)) then
        dx = 0.5_rk * a / b
      else
        dx = 0
      end if
      if (dx<xm) dx = xm
      if (dx>xp) dx = xp
      xmax = x(loc) + dx
    end function interpolate_max_position
    !
    subroutine write_common_header(iu)
      integer(ik), intent(in) :: iu
      !
      write (iu,"('# ',a)") trim(comment)
      write (iu,"('# Total wavefunction')")
      write (iu,"('#              M = ',i0)") mval
      write (iu,"('#           ZNUC = ',g36.26e3)") znuc
      write (iu,"('#      Potential = ',a)") trim(potential)
      write (iu,"('#             a1 = ',g36.26e3)") pot_param_real(1)
      write (iu,"('#             a2 = ',g36.26e3)") pot_param_real(2)
      write (iu,"('#             n1 = ',i0)") pot_param_int(1)
      write (iu,"('#              E = ',2g36.26e3)") energy
      write (iu,"('#              F = ', g36.26e3)") efield
      write (iu,"('#           Mode = ',a)") trim(file_total_mode)
    end subroutine write_common_header
  end subroutine dump_total_cartesian_wavefunction
end module general_tunnel_dump

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
!  Evaluation of channel coupling parameters for general_tunnel.f90
!
module general_tunnel_nonadiabatic
  use accuracy
  use constants
  use derivative_tools
  use general_tunnel_data
  use general_tunnel_bound
  use general_tunnel_dump
  use timer
  !$ use OMP_LIB
  implicit none
  private
  public prepare_channel_coupling
  public rcsid_general_tunnel_nonadiabatic
  !
  character(len=clen), save :: rcsid_general_tunnel_nonadiabatic = "$Id: $"
  !
  contains
  !
  !  Compute non-adiabatic couplings for a single right-hand-side channel
  !  This subroutine is safe to call from a parallel region.
  !
  subroutine channel_nonadiabatic_single(jchan,eta_start,eta_end)
    integer(ik), intent(in)    :: jchan                 ! Channel in which we apply wavefunction derivative
    integer(ik), intent(in)    :: eta_start, eta_end    ! Block of eta points we are interested in. Intended
                                                        ! for running multiple threads in parallel
    !
    integer(ik)                :: eta_ipt               ! Current target point
    integer(ik)                :: eta_1st               ! First point within the differentiation window
    integer(ik)                :: save_nl               ! Number of "left" points in the current wder
    integer(ik)                :: ichan
    integer(ik)                :: alloc
    real(rk)                   :: grid_d                ! Uniform grid spacing. We can't handle non-uniform grids here!
    real(rk)                   :: wder(nonad_points,2)  ! Weights for derivative evaluation
    complex(rk), allocatable   :: wf_buf(:,:,:)         ! Wavefunction buffer for differentiation
    type(wavefunction)         :: wf_der                ! Derivative of the bound wavefunction (1st or second) with respect to eta.
    !
    !  In order to evaluate the non-adiabatic matrix elements, we need to calculate
    !  the derivatives of the bound solutions within each channel
    !
    if (eta_npts<nonad_points) then
      stop 'general_tunnel%channel_nonadiabatic_single - not enough points to differentiate'
    end if
    if (eta_start<1 .or. eta_end>eta_npts) then
      stop 'general_tunnel%channel_nonadiabatic_single - bad eta range'
    end if
    !
    allocate (wf_buf(0:xi_maxorder,nonad_points,xi_npts),stat=alloc)
    if (alloc/=0) stop 'general_tunnel%prepare_channel_coupling - allocation failed (1)'
    call init_wavefunction_1D(wf_der,xi_tab,xi_maxorder)  ! general_tunnel_data
    !
    grid_d  = eta_tab(2) ! The first point is guaranteed to be at zero, so no worries here
    save_nl = -1_ik      ! An impossible value to force evaluation first-time around
    target_eta: do eta_ipt=eta_start,eta_end
      !
      !  Position the differentiation window within the grid and prepare the
      !  weigts. For the interior points, we can skip evaluation of the weights
      !  most of the time.
      !
      eta_1st = min(max(1_ik,eta_ipt-(nonad_points/2_ik)),eta_npts-nonad_points+1_ik)
      if (save_nl /= eta_ipt-eta_1st) then
        save_nl = eta_ipt-eta_1st
        call dt_savgol(nonad_order,save_nl,grid_d,1_ik,wder) ! 1 is the gradient; we also get the laplacian
        wder(:,2) = 2.0_rk * wder(:,2)  ! Laplacian: multiply by 2! = 2
      end if
      call fetch_wavefunction(channels(jchan),eta_1st,nonad_points,wf_buf)
      !
      !  1st derivative
      ! 
      call differentiate(wder(:,1),wf_buf,wf_der%wf)
      gradient_elements: do ichan=1,n_channels
        chan_w1(0,ichan,jchan,eta_ipt) = bound_overlap(channels(ichan)%bound(eta_ipt),wf_der)
      end do gradient_elements
      !
      !  2nd derivative
      !
      call differentiate(wder(:,2),wf_buf,wf_der%wf)
      laplacian_elements: do ichan=1,n_channels
        chan_w2(0,ichan,jchan,eta_ipt) = bound_overlap(channels(ichan)%bound(eta_ipt),wf_der)
      end do laplacian_elements
    end do target_eta
    call destroy_wavefunction_1D(wf_der)
    deallocate (wf_buf)
    !
    contains
    !
    subroutine fetch_wavefunction(chan,e1,ne,wf)
      type(channel), intent(in) :: chan
      integer(ik), intent(in)   :: e1                           ! The first eta index we need
      integer(ik), intent(in)   :: ne                           ! Number of eta indices
      complex(rk), intent(out)  :: wf(0:xi_maxorder,ne,xi_npts) ! Wavefunction; 
                                                                !    first index = derivative order
                                                                !    second index = eta points
                                                                !    last index = xi (bound) points; 
      !
      integer(ik) :: xi_ipt, eta_ind
      !
      xi_points: do xi_ipt=1,xi_npts
        !
        !  Fetch wavefunction; wavefunction past ipt_stop is zero, so we don't need to worry about it.
        !
        eta_points: do eta_ind=1,ne
          wf(:,eta_ind,xi_ipt) = chan%bound(e1+eta_ind-1)%wf(:,xi_ipt)
        end do eta_points
      end do xi_points
    end subroutine fetch_wavefunction
    !
    subroutine differentiate(wgt,wf_buf,der)
      real(rk), intent(in)     :: wgt   (  :  )   ! Weights of eta points in the derivative
      complex(rk), intent(in)  :: wf_buf(:,:,:)   ! [1] = order in xi; [2] = eta point; [3] = xi point
      complex(rk), intent(out) :: der   (:,  :)   ! [1] = order in xi; [2] = xi point
      !
      integer(ik) :: iord, xi_ipt
      !
      if (size(wgt)/=size(wf_buf,dim=2) .or. size(wf_buf,dim=1)/=size(der,dim=1) .or. size(wf_buf,dim=3)/=size(der,dim=2)) then
        stop 'general_tunnel%channel_nonadiabatic_single%differentiate - incompatible arrays'
      end if
      !
      xi_points: do xi_ipt=1,size(der,dim=2)
        orders: do iord=1,size(der,dim=1)
          der(iord,xi_ipt) = sum(wgt(:)*wf_buf(iord,:,xi_ipt))
        end do orders
      end do xi_points
    end subroutine differentiate
  end subroutine channel_nonadiabatic_single
  !
  !  Differentiate chan_alp, chan_w1, and chan_w2 for a block of points.
  !  Thread-safe.
  !
  subroutine channel_nonadiabatic_differentiate(eta_start,eta_end)
    integer(ik), intent(in)    :: eta_start, eta_end              ! Block of eta points we are interested in. Intended
                                                                  ! for running multiple threads in parallel
    !                                                             
    integer(ik)                :: eta_ipt                         ! Current target point
    integer(ik)                :: eta_1st                         ! First point within the differentiation window
    integer(ik)                :: eta_x                           ! Last point in same
    integer(ik)                :: save_nl                         ! Number of "left" points in the current wder
    integer(ik)                :: ichan, jchan                    ! Channel indices
    real(rk)                   :: grid_d                          ! Uniform grid spacing. We can't handle non-uniform grids here!
    real(rk)                   :: wder(eta_points2,eta_maxorder2) ! Weights for derivative evaluation
    !
    !  We need at least eta_points2 in the eta grid - otherwise, there is not enough data to fit to!
    !
    if (eta_npts<eta_points2) then
      stop 'general_tunnel%channel_nonadiabatic_differentiate - not enough points to differentiate'
    end if
    if (eta_start<1 .or. eta_end>eta_npts) then
      stop 'general_tunnel%channel_nonadiabatic_differentiate - bad eta range'
    end if
    !
    grid_d  = eta_tab(2) ! The first point is guaranteed to be at zero, so no worries here
    save_nl = -1_ik      ! An impossible value to force evaluation first-time around
    target_eta: do eta_ipt=eta_start,eta_end
      !
      !  Position the differentiation window within the grid and prepare the
      !  weigts. For the interior points, we can skip evaluation of the weights
      !  most of the time.
      !
      eta_1st = min(max(1_ik,eta_ipt-(eta_points2/2_ik)),eta_npts-eta_points2+1_ik)
      eta_x   = eta_1st+eta_points2-1
      if (save_nl /= eta_ipt-eta_1st) then
        save_nl = eta_ipt-eta_1st
        call dt_savgol(eta_order2,save_nl,grid_d,1_ik,wder) ! 1 is the gradient; we get eta_maxorder2 orders altogether
      end if
      !
      channel_rhs: do jchan=1,n_channels
          chan_alp(1:,jchan,      eta_ipt) = matmul(chan_alp(0,      jchan,eta_1st:eta_x),wder)
        channel_lhs: do ichan=1,n_channels
          chan_w1 (1:,ichan,jchan,eta_ipt) = matmul(chan_w1 (0,ichan,jchan,eta_1st:eta_x),wder)
          chan_w2 (1:,ichan,jchan,eta_ipt) = matmul(chan_w2 (0,ichan,jchan,eta_1st:eta_x),wder)
        end do channel_lhs
      end do channel_rhs
    end do target_eta
  end subroutine channel_nonadiabatic_differentiate
  !
  !  Initialize channel separation parameters, coupling terms and their derivatives
  !
  subroutine prepare_channel_coupling(verbose)
    integer(ik), intent(in) :: verbose
    integer(ik)             :: ichan, jchan, threads
    integer(ik)             :: eta_start, eta_end, eta_block
    character(len=clen)     :: tag
    !
    call TimerStart('Channel coupling')
    !
    !  Copy the separation parameters
    !
    fill_separation: do jchan=1,n_channels
      chan_alp(0,jchan,:) = channels(jchan)%bound(:)%zeta
    end do fill_separation
    !
    !  Calculate non-adiabatic matrix elements. This is only needed for non-hydrogenic
    !  potentials; otherwise we know these matrix elements are all zero.
    !
    if (potential=='hydrogenic') then
      chan_alp(1:,:,:)  = 0._rk
      chan_w1 (:,:,:,:) = 0._rk
      chan_w2 (:,:,:,:) = 0._rk
    else
      threads   = 1
      !$ threads = omp_get_max_threads()
      eta_block = (eta_npts + 2*threads - 1) / (2*threads)
      if (verbose>=2) then
        write (out,"('Using eta_block = ',i0)") eta_block
      end if
      if (eta_block<1) stop 'general_tunnel%prepare_channel_coupling - oops'
      !$omp parallel do collapse(2) default(private) shared(n_channels,eta_npts,eta_block)
      rhs_channels: do jchan=1,n_channels
        eta_blocks: do eta_start=1,eta_npts,eta_block
          eta_end = min(eta_npts,eta_start+eta_block-1)
          call channel_nonadiabatic_single(jchan,eta_start,eta_end)
        end do eta_blocks
      end do rhs_channels
      !$omp end parallel do
      !
      !  Numerically differentiate non-adiabatic matrix elements
      !
      !$omp parallel do default(private) shared(eta_npts,eta_block)
      eta_grad_blocks: do eta_start=1,eta_npts,eta_block
        eta_end = min(eta_npts,eta_start+eta_block-1)
        call channel_nonadiabatic_differentiate(eta_start,eta_end)
      end do eta_grad_blocks
      !$omp end parallel do
      !
      !  Debugging output
      !
      if (verbose>=1) then
        report_rhs_channel: do jchan=1,n_channels
          write (tag,"('alp_',i0)") jchan
          call dump_channel_parameter(trim(tag),chan_alp(:,jchan,:))
          report_lhs_channel: do ichan=1,n_channels
            write (tag,"('w1_',i0,'_',i0)") ichan, jchan
            call dump_channel_parameter(trim(tag),chan_w1(:,ichan,jchan,:))
            write (tag,"('w2_',i0,'_',i0)") ichan, jchan
            call dump_channel_parameter(trim(tag),chan_w2(:,ichan,jchan,:))
          end do report_lhs_channel
        end do report_rhs_channel
      end if
      !
      !  Gradient coupling within the same channel must vanish (but may be
      !  non-zero due to numerical noise). We'll report it for diagnostic
      !  purposes, but zero it out for the calculations of the wavefunction.
      !
      if (verbose>=1) then
        write (out,"(/'Setting same-channel gradient coupling to zero'/)") 
      end if
      clear_same_channel: do jchan=1,n_channels
        chan_w1(:,jchan,jchan,:) = 0._rk
      end do clear_same_channel
    endif
    !
    call TimerStop('Channel coupling')
    !
  end subroutine prepare_channel_coupling
end module general_tunnel_nonadiabatic

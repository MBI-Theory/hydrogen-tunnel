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
!  Routines dealing with the continuum-direction integration for general_tunnel.f90
!
module general_tunnel_continuum
  use accuracy
  use constants
  use derivative_tools
  use find_root
  use find_minimum
  use general_tunnel_data
  use general_tunnel_asymptotic
  use general_tunnel_dump
  use lapack
  use math
  use poly_tools
  use sort_tools
  use timer
  !$ use OMP_LIB
  implicit none
  private
  public solve_continuum, invert_continuum
  public rcsid_general_tunnel_continuum
  !
  character(len=clen), save :: rcsid_general_tunnel_continuum = "$Id: $"
  !
  contains
  !
  !  Evaluation of the power-series is quite messy, so it is cleaner to keep it separate.
  !  We need two separate routines: one for the origin (where we have do deal with a singularity)
  !  and one for the general interior point.
  !
  !  This routine evaluates a single order of the power series at the origin:
  !
  subroutine radial_continuum_origin_step(chn,energy,efield,alp,w1,w2,iord,dx,max_c)
    type (wavefunction), intent(inout) :: chn     (n_channels)   ! Per-channel wavefunctions; see global channels for the description
    complex(rk), intent(in)            :: energy                 ! Energy of the solution
    real(rk), intent(in)               :: efield                 ! Electric field. Expected to be > 0
    complex(rk), intent(in)            :: alp(0:eta_maxorder2,n_channels,           eta_npts)
    complex(rk), intent(in)            :: w1 (0:eta_maxorder2,n_channels,n_channels,eta_npts)
    complex(rk), intent(in)            :: w2 (0:eta_maxorder2,n_channels,n_channels,eta_npts)
    integer(ik), intent(in)            :: iord                   ! Coefficient order we are currently working on. Must be at least mval+1!
    real(rk), intent(in)               :: dx(0:eta_maxorder)     ! Powers of grid spacing
    real(rk), intent(inout)            :: max_c                  ! Largest power series coefficient seen so far
    !
    integer(ik) :: ic, jc                        ! Channel numbers
    integer(ik) :: jord, jord_max                ! Order index
    complex(rk) :: temp
    !
    if (iord<=mval) stop 'general_tunnel_continuum%radial_continuum_origin_step - called for an invalid iord'
    !
    !$omp parallel do default(none) num_threads(n_channels) &
    !$omp&  shared(n_channels,chn,energy,efield,alp,w1,w2,iord,dx) &
    !$omp&  shared(znuc,eta_maxorder2,mval) &
    !$omp&  private(temp,ic,jc,jord,jord_max) reduction(max:max_c)
    origin_channels: do ic=1,n_channels
      temp = 0
      !
      !  Intra-channel contributions
      !
      if (iord>=2) temp = temp + dx(2) * chn(ic)%wf(iord-2,1) * 4._rk * znuc
      if (iord>=4) temp = temp + dx(4) * chn(ic)%wf(iord-4,1) * 2._rk * energy
      if (iord>=6) temp = temp + dx(6) * chn(ic)%wf(iord-6,1) * efield
      !
      !  Adiabatic separation parameter, same channel.
      !
      jord_max = min(iord-2,eta_maxorder2)
      origin_adiabatic: do jord=0,jord_max
        temp = temp - dx(2+jord) * chn(ic)%wf(iord-jord-2,1) * alp(jord,ic,1)
      end do origin_adiabatic
      !
      !  Non-adiabatic terms, cross-channel.
      !
      origin_cross_channels: do jc=1,n_channels
        jord_max = min(iord-1,eta_maxorder2)
        origin_gradient: do jord=0,jord_max
          temp = temp + dx(1+jord) * chn(jc)%wf(iord-jord-1,1) * (2._rk*(iord-jord)-1._rk) * w1(jord,ic,jc,1)
        end do origin_gradient
        !
        jord_max = min(iord-2,eta_maxorder2)
        origin_laplacian: do jord=0,jord_max
          temp = temp + dx(2+jord) * chn(jc)%wf(iord-jord-2,1) * w2(jord,ic,jc,1)
        end do origin_laplacian
      end do origin_cross_channels
      !
      chn(ic)%wf(iord,1) = -temp / (iord**2 - mval**2)
      max_c = max(max_c,abs(chn(ic)%wf(iord,1)))
    end do origin_channels
    !$omp end parallel do
  end subroutine radial_continuum_origin_step
  !
  !  This routine evaluates a single order in the power series at an interior grid point:
  !
  subroutine radial_continuum_interior_step(chn,energy,efield,alp,w1,w2,ipt,iord,dx,max_c)
    type (wavefunction), intent(inout) :: chn     (n_channels)   ! Per-channel wavefunctions; see global channels for the description
    complex(rk), intent(in)            :: energy                 ! Energy of the solution
    real(rk), intent(in)               :: efield                 ! Electric field. Expected to be > 0
    complex(rk), intent(in)            :: alp(0:eta_maxorder2,n_channels,           eta_npts)
    complex(rk), intent(in)            :: w1 (0:eta_maxorder2,n_channels,n_channels,eta_npts)
    complex(rk), intent(in)            :: w2 (0:eta_maxorder2,n_channels,n_channels,eta_npts)
    integer(ik), intent(in)            :: ipt                    ! Grid index where we are currently working on. Can't be 1!
    integer(ik), intent(in)            :: iord                   ! Coefficient order we are currently working on. Must be at least 2!
    real(rk), intent(in)               :: dx(0:eta_maxorder)     ! Powers of grid spacing
    real(rk), intent(inout)            :: max_c                  ! Largest power series coefficient seen so far
    !
    integer(ik) :: ic, jc                        ! Channel numbers
    integer(ik) :: jord, jord_max                ! Order index
    complex(rk) :: temp, temp2
    real(rk)    :: eta
    real(rk)    :: c7                            ! Fixed part of the recursion coefficients, real
    complex(rk) :: c2, c3, c4, c5, c6            ! ditto, complex
    !
    if (ipt <=1) stop 'general_tunnel_continuum%radial_continuum_interior_step - called for an invalid ipt'
    if (iord<=1) stop 'general_tunnel_continuum%radial_continuum_interior_step - called for an invalid iord'
    !
    eta = eta_tab(ipt)
    c2  = -(mval**2) + 4._rk*znuc*eta**2 + 2._rk*energy*eta**4 + efield*eta**6
    c3  = 8._rk*znuc*eta + 8._rk*energy*eta**3 + 6._rk*efield*eta**5
    c4  = 4._rk*znuc + 12._rk*energy*eta**2 + 15._rk*efield*eta**4
    c5  = 8._rk*energy*eta + 20._rk*efield*eta**3
    c6  = 2._rk*energy + 15._rk*efield*eta**2
    c7  = 6._rk*efield*eta
    !
    !$omp parallel do default(none) num_threads(n_channels) &
    !$omp&  shared(n_channels,chn,energy,efield,alp,w1,w2,ipt,iord,dx,eta,c2,c3,c4,c5,c6,c7) &
    !$omp&  shared(eta_maxorder2) &
    !$omp&  private(temp,temp2,ic,jc,jord,jord_max) reduction(max:max_c)
    propagate_channels: do ic=1,n_channels
      temp = 0
      ! 
      !  Intra-channel contributions
      !
                   temp = temp + dx(1) * chn(ic)%wf(iord-1,ipt) * eta*(iord-1._rk)*(2._rk*iord-3._rk)
                   temp = temp + dx(2) * chn(ic)%wf(iord-2,ipt) * ((iord-2)**2 + c2)
      if (iord>=3) temp = temp + dx(3) * chn(ic)%wf(iord-3,ipt) * c3
      if (iord>=4) temp = temp + dx(4) * chn(ic)%wf(iord-4,ipt) * c4
      if (iord>=5) temp = temp + dx(5) * chn(ic)%wf(iord-5,ipt) * c5
      if (iord>=6) temp = temp + dx(6) * chn(ic)%wf(iord-6,ipt) * c6
      if (iord>=7) temp = temp + dx(7) * chn(ic)%wf(iord-7,ipt) * c7
      if (iord>=8) temp = temp + dx(8) * chn(ic)%wf(iord-8,ipt) * efield
      !
      !  Adiabatic separation parameter, same channel.
      !  We have three distinct terms appearing here.
      !
      temp2    = 0
      jord_max = min(iord-2,eta_maxorder2)
      propagate_adiabatic_1: do jord=0,jord_max
        temp2 = temp2 + dx(2+jord) * chn(ic)%wf(iord-jord-2,ipt) * alp(jord,ic,ipt)
      end do propagate_adiabatic_1
      temp     = temp - eta**2 * temp2
      !
      temp2    = 0
      jord_max = min(iord-3,eta_maxorder2)
      propagate_adiabatic_2: do jord=0,jord_max
        temp2 = temp2 + dx(3+jord) * chn(ic)%wf(iord-jord-3,ipt) * alp(jord,ic,ipt)
      end do propagate_adiabatic_2
      temp     = temp - 2._rk*eta * temp2
      !
      jord_max = min(iord-4,eta_maxorder2)
      propagate_adiabatic_3: do jord=0,jord_max
        temp = temp - dx(4+jord) * chn(ic)%wf(iord-jord-4,ipt) * alp(jord,ic,ipt)
      end do propagate_adiabatic_3
      !
      !  Non-adiabatic terms, cross-channel.
      !  We have four distinct terms here.
      !
      propagate_cross_channels: do jc=1,n_channels
        jord_max = min(iord-1,eta_maxorder2)
        propagate_nonadiabatic_1: do jord=0,jord_max
          temp2 = (2._rk*(iord-jord)-2._rk)*w1(jord,ic,jc,ipt)*eta**2
          temp  = temp + dx(1+jord) * chn(jc)%wf(iord-jord-1,ipt) * temp2
        end do propagate_nonadiabatic_1
        !
        jord_max = min(iord-2,eta_maxorder2)
        propagate_nonadiabatic_2: do jord=0,jord_max
          temp2 = (4._rk*(iord-jord)-7._rk)*w1(jord,ic,jc,ipt)*eta + w2(jord,ic,jc,ipt)*eta**2
          temp  = temp + dx(2+jord) * chn(jc)%wf(iord-jord-2,ipt) * temp2
        end do propagate_nonadiabatic_2
        !
        jord_max = min(iord-3,eta_maxorder2)
        propagate_nonadiabatic_3: do jord=0,jord_max
          temp2 = (2._rk*(iord-jord)-5._rk)*w1(jord,ic,jc,ipt) + 2._rk*w2(jord,ic,jc,ipt)*eta
          temp  = temp + dx(3+jord) * chn(jc)%wf(iord-jord-3,ipt) * temp2
        end do propagate_nonadiabatic_3
        !
        jord_max = min(iord-4,eta_maxorder2)
        propagate_nonadiabatic_4: do jord=0,jord_max
          temp  = temp + dx(4+jord) * chn(jc)%wf(iord-jord-4,ipt) * w2(jord,ic,jc,ipt)
        end do propagate_nonadiabatic_4
        !
      end do propagate_cross_channels
      chn(ic)%wf(iord,ipt) = - temp / (iord*(iord-1)*eta**2)
      max_c = max(max_c,abs(chn(ic)%wf(iord,ipt)))
    end do propagate_channels
    !$omp end parallel do
  end subroutine radial_continuum_interior_step
  !
  !  Ourward integration for the continuum solution. This routine is a bit messier than 
  !  its bound counterpart, since we now have multiple, coupled channels.
  !
  !  This routine requires eta grid spacing to be uniform.
  !
  subroutine integrate_radial_continuum_outbound(chn,boundary,energy,efield,alp,w1,w2)
    type (wavefunction), intent(inout) :: chn     (n_channels)   ! Per-channel wavefunctions; see global channels for the description
    complex(rk), intent(in)            :: boundary(n_channels)   ! Per-channel boundary conditions at the origin. 
    complex(rk), intent(in)            :: energy                 ! Energy of the solution
    real(rk), intent(in)               :: efield                 ! Electric field. Expected to be > 0
    complex(rk), intent(in)            :: alp(0:eta_maxorder2,n_channels,           eta_npts)  
    complex(rk), intent(in)            :: w1 (0:eta_maxorder2,n_channels,n_channels,eta_npts)  
    complex(rk), intent(in)            :: w2 (0:eta_maxorder2,n_channels,n_channels,eta_npts)  
                                                                 ! alp: Adiabatic separation parameters for the bound soulutions
                                                                 ! w1:  Gradient non-adiabatic coupling
                                                                 ! w2:  Laplacian non-adiabatic coupling
                                                                 ! See global chan_alp, chan_w1, and chan_w2 for the 
                                                                 ! description of the indices.
    !
    integer(ik) :: ic                            ! Channel number
    integer(ik) :: ipt                           ! Grid index
    integer(ik) :: iord                          ! Order index
    real(rk)    :: dx(0:eta_maxorder)            ! Powers of grid spacing
    real(rk)    :: max_c                         ! Largest power series coefficient seen so far
    integer(ik) :: zero_c_count                  ! Number of consequitive negligible power series coefficients
    !
    !  Sanity checking
    !
    if (efield<=0._rk) then
      stop 'general_tunnel%integrate_radial_continuum_outbound - efield must be positive'
    end if
    !
    dx(0) = 1._rk
    dx(1) = eta_tab(2) - eta_tab(1)
    !
    fill_step_powers: do iord=2,eta_maxorder
      dx(iord) = dx(iord/2) * dx(iord-iord/2) ! Minimize error propagation by using binary-tree formula
    end do fill_step_powers
    !
    !  Prepare for integration
    !
    init_channels: do ic=1,n_channels
      chn(ic)%energy      = energy
      chn(ic)%efield      = efield
      chn(ic)%norm        = 'natural'
      chn(ic)%wf (   :,1) = 0
      chn(ic)%wf (mval,1) = boundary(ic) * dx(mval)
      chn(ic)%iord(:)     = eta_maxorder
    end do init_channels
    !
    !  The origin of the grid requires special handling due to the
    !  presence of singularities in the Hamiltonian.
    !
    max_c = maxval(abs(boundary(:))) * dx(mval)
    zero_c_count = 0
    origin_orders: do iord=mval+1,eta_maxorder
      call radial_continuum_origin_step(chn,energy,efield,alp,w1,w2,iord,dx,max_c)
      !
      !  If enough coefficients become negligible, we can drop out early
      !
      origin_early_exit: do ic=1,n_channels
        if (abs(chn(ic)%wf(iord,1))>spacing(max_c)) exit origin_early_exit
      end do origin_early_exit
      if (ic<=n_channels) then
        ! One of the coefficients is larger than the numerical epsilon
        zero_c_count = 0
      else  
        ! All coefficienta are less than numerical epsilon
        zero_c_count = zero_c_count + 1
        if (zero_c_count>=6) then
          origin_stop: do ic=1,n_channels
            chn(ic)%iord(1) = iord
            chn(ic)%wf(iord+1:,1) = 0
          end do origin_stop
          exit origin_orders
        end if
      end if
    end do origin_orders
    !
    !  We can use general code for the rest of the radial grid
    !
    propagate: do ipt=2,eta_npts
      !
      !  Use power-series expansion at the previous point to get wavefunction and its gradient
      !
      max_c = 0._rk
      propagate_channel_step: do ic=1,n_channels
        iord = chn(ic)%iord(ipt-1)
        call poly_power_series(chn(ic)%wf(:iord,ipt-1),1._rk,chn(ic)%wf(0:1,ipt))
        max_c = max(max_c,maxval(abs(chn(ic)%wf(0:1,ipt))))
      end do propagate_channel_step
      !
      !  Calculate fixed power-series coefficients at the new point
      !  From the second term on, we can use the general expression
      !
      zero_c_count = 0
      propagate_orders: do iord=2,eta_maxorder
        call radial_continuum_interior_step(chn,energy,efield,alp,w1,w2,ipt,iord,dx,max_c)
        !
        !  If enough coefficients become negligible, we can drop out early
        !
        early_exit: do ic=1,n_channels
          if (abs(chn(ic)%wf(iord,ipt))>spacing(max_c)) exit early_exit
        end do early_exit
        if (ic<=n_channels) then
          ! One of the coefficients is larger than the numerical epsilon
          zero_c_count = 0
        else  
          ! All coefficienta are less than numerical epsilon
          zero_c_count = zero_c_count + 1
          if (zero_c_count>=8) then
            general_stop: do ic=1,n_channels
              chn(ic)%iord(ipt) = iord
              chn(ic)%wf(iord+1:,ipt) = 0
            end do general_stop
            exit propagate_orders
          end if
        end if
      end do propagate_orders
    end do propagate
  end subroutine integrate_radial_continuum_outbound
  ! 
  !  Inbound integration of the radial continuum from the outer matching point,
  !  This routine is only used once, in the reporting part. There is not much
  !  sense in optimizing it too much!
  !
  subroutine integrate_radial_continuum_inbound(chn,energy,efield,alp,w1,w2)
    type (wavefunction), intent(inout) :: chn     (n_channels)   ! Per-channel wavefunctions; see global channels for the description
    complex(rk), intent(in)            :: energy                 ! Energy of the solution
    real(rk), intent(in)               :: efield                 ! Electric field. Expected to be > 0
    complex(rk), intent(in)            :: alp(0:eta_maxorder2,n_channels,           eta_npts)  
    complex(rk), intent(in)            :: w1 (0:eta_maxorder2,n_channels,n_channels,eta_npts)  
    complex(rk), intent(in)            :: w2 (0:eta_maxorder2,n_channels,n_channels,eta_npts)  
                                                                 ! alp: Adiabatic separation parameters for the bound soulutions
                                                                 ! w1:  Gradient non-adiabatic coupling
                                                                 ! w2:  Laplacian non-adiabatic coupling
                                                                 ! See global chan_alp, chan_w1, and chan_w2 for the 
                                                                 ! description of the indices.
    !
    integer(ik) :: ic                    ! Channel number
    integer(ik) :: ipt                   ! Grid index
    integer(ik) :: iord                  ! Order index
    real(rk)    :: dx(0:eta_maxorder)    ! Powers of grid spacing
    real(rk)    :: max_c                 ! Largest power series coefficient seen so far. Not used, but required by worker routines
    !
    !  Sanity checking
    !
    if (efield<=0._rk) then
      stop 'general_tunnel%integrate_radial_continuum_inbound - efield must be positive'
    end if
    !
    dx(0) = 1._rk
    dx(1) = eta_tab(2) - eta_tab(1)
    !
    fill_step_powers: do iord=2,eta_maxorder
      dx(iord) = dx(iord/2) * dx(iord-iord/2) ! Minimize error propagation by using binary-tree formula
    end do fill_step_powers
    !
    !  Prepare for integration
    !
    init_channels: do ic=1,n_channels
      chn(ic)%energy      = energy
      chn(ic)%efield      = efield
      chn(ic)%norm        = 'natural'
      chn(ic)%iord(:)     = eta_maxorder
    end do init_channels
    !
    !  Start at the outrmost point, and integrate inward
    !
    propagate: do ipt=eta_npts,2,-1
      !
      !  We already have function and it gradient. Use recursions to evaluate all higher orders
      !  Do not bother trying to bail out early: this routine is not on the critical path!
      !
      max_c = 0._rk
      propagate_orders: do iord=2,eta_maxorder
        call radial_continuum_interior_step(chn,energy,efield,alp,w1,w2,ipt,iord,dx,max_c)
      end do propagate_orders
      ! 
      !  Step to the previous point
      !
      propagate_channel_step: do ic=1,n_channels
        call poly_power_series(chn(ic)%wf(:,ipt),-1._rk,chn(ic)%wf(0:1,ipt-1))
      end do propagate_channel_step
    end do propagate
    !
    !  Note that there is no point in trying to evaluate higher-order terms at the origin:
    !  we have almost certainly already violated the boundary conditions at the origin,
    !  making our recursions (which rely on these boundary conditions to be satisfied!)
    !  meaningless.
    !
  end subroutine integrate_radial_continuum_inbound
  ! 
  subroutine channel_asymptotic_solutions(ichan,efield,energy,zeta,va,ga,vb,gb,wrn)
    integer(ik), intent(in)  :: ichan           ! Channel number (for reporting purposes; not needed otherwise)
    real(rk), intent(in)     :: efield          ! Electric field (must be positive)
    complex(rk), intent(in)  :: energy          ! Energy of the solution
    complex(rk), intent(in)  :: zeta            ! Hydrogenic separation parameter
    complex(rk), intent(out) :: va, ga          ! Outgoing solution and gradient
    complex(rk), intent(out) :: vb, gb          ! Incoming solution and gradient
    complex(rk), intent(out) :: wrn             ! Wronskian for the solution
    !
    complex(rk) :: cna( 1:asymp_order)   ! Asymptotic outgoing solution
    complex(rk) :: cnb( 1:asymp_order)   ! Asymptotic incoming solution
    complex(rk) :: gra(-1:asymp_order+1) ! Gradient of the asymptotic outgoing solution
    complex(rk) :: grb(-1:asymp_order+1) ! Gradient of the asymptotic incoming solution
    complex(rk) :: ewrn                  ! Error in the Wronskian
    !
    !  Form asymptotic solutions. These solutions contain an additional factor eta**(1/2),
    !  which is needed to eliminate gradient term from the Hamiltonian.
    !
    call asymptotic_solution(efield,energy,zeta, 1._rk,eta_max,cna)
    call asymptotic_solution(efield,energy,zeta,-1._rk,eta_max,cnb)
    call asymptotic_gradient(efield,energy, 1._rk,eta_max,cna,gra)
    call asymptotic_gradient(efield,energy,-1._rk,eta_max,cnb,grb)
    va   = evaluate_asymptote(efield,energy, 1._rk,cna,eta_max)
    ga   = evaluate_asymptote(efield,energy, 1._rk,gra,eta_max)
    vb   = evaluate_asymptote(efield,energy,-1._rk,cnb,eta_max)
    gb   = evaluate_asymptote(efield,energy,-1._rk,grb,eta_max)
    wrn  = ga * vb - gb * va
    ewrn = wrn - (0,2)*sqrt(efield)
    !
    if (abs(ewrn)>=20*spacing(sqrt(efield))) then
      write (out,"(/'WARNING: Wronskian deviates from the exact solution [2*I*sqrt(efield)] by ',g14.6e3)") abs(ewrn)
      write (out,"( 'WARNING: Calculated Wronskian: ',2g45.32e3)") wrn
      write (out,"( 'WARNING:      Exact Wronskian: ',2g45.32e3)") (0,2)*sqrt(efield)
      write (out,"( 'WARNING: Incoming and outgoing amplitudes in channel ',i0,' may be inaccurate'/)") ichan
      call flush_wrapper (out)
    end if
  end subroutine channel_asymptotic_solutions
  !
  !  All routines with the names starting with "solve_continuum_" are meant to be used
  !  only by solve_sontinuum. They could be converted to contained subroutines, but
  !  this makes the overall routine too unwieldy for my taste.
  !
  subroutine solve_continuum_asymptotic_solutions(verbose,energy,va,ga,vb,gb,wrn)
    integer(ik), intent(in)  :: verbose
    complex(rk), intent(in)  :: energy
    complex(rk), intent(out) :: va(n_channels)        ! Outgoing solutions at the matching point
    complex(rk), intent(out) :: ga(n_channels)        ! Gradients of the outgoing solutions
    complex(rk), intent(out) :: vb(n_channels)        ! Incoming solutions at the matching point
    complex(rk), intent(out) :: gb(n_channels)        ! Gradients of the incoming solutions
    complex(rk), intent(out) :: wrn(n_channels)       ! The Wronskian of the incoming/outgoing solution pair
    !
    integer(ik) :: ic
    complex(rk) :: zeta   ! Hydrogenic separation parameter; taken from the last radial point
    !
    call TimerStart('Asymptotic solutions')
    !$omp parallel do default(shared) private(ic,zeta) num_threads(n_channels)
    asymptotic_channels: do ic=1,n_channels
      zeta = 4._rk * znuc - channels(ic)%bound(eta_npts)%zeta
      continuum(ic)%zeta = zeta
      call channel_asymptotic_solutions(ic,efield,energy,zeta,va(ic),ga(ic),vb(ic),gb(ic),wrn(ic))
    end do asymptotic_channels
    !$omp end parallel do
    if (verbose>=1) then
      write (out,"(/t5,'Asymptotic solutions at ETA = ',g16.8/)") eta_max
      write (out,"((1x,a5,5(2x,a16,1x,a16)))") &
             'Chan', 'Outgoing', '', 'Grad Outgoing', '', 'Incoming', '', 'Grad Incoming', '', 'Wronskian', '', &
             '----', '--------', '', '-------------', '', '--------', '', '-------------', '', '---------', ''
      report_asymptotics: do ic=1,n_channels
        write (out,"(1x,i5,5(2x,g16.8,1x,g16.8))") ic, va(ic), ga(ic), vb(ic), gb(ic), wrn(ic)
      end do report_asymptotics
      write (out,"()")
      call flush_wrapper (out)
    end if
    call TimerStop('Asymptotic solutions')
  end subroutine solve_continuum_asymptotic_solutions
  !
  subroutine solve_continuum_integrate_outbound(verbose,energy,va,ga,vb,gb,wrn,boundary,continuum)
    integer(ik), intent(in)           :: verbose
    complex(rk), intent(in)           :: energy
    complex(rk), intent(in)           :: va(n_channels)     ! Outgoing solutions at the matching point
    complex(rk), intent(in)           :: ga(n_channels)     ! Gradients of the outgoing solutions
    complex(rk), intent(in)           :: vb(n_channels)     ! Incoming solutions at the matching point
    complex(rk), intent(in)           :: gb(n_channels)     ! Gradients of the incoming solutions
    complex(rk), intent(in)           :: wrn(n_channels)    ! The Wronskian of the incoming/outgoing solution pair
    complex(rk), intent(in)           :: boundary(n_channels)
    type(wavefunction), intent(inout) :: continuum(n_channels)
    !
    integer(ik) :: ic
    complex(rk) :: vs, gs    ! Function and gradient in scaled coordinates
    !
    call TimerStart('Continuum outbound')
    call integrate_radial_continuum_outbound(continuum,boundary,energy,efield,chan_alp,chan_w1,chan_w2)
    !
    !  Match the two solutions at the exterior edge of the grid
    !
    match_channels: do ic=1,n_channels
      !
      !  The asymptotic solutions include an extra eta**(1/2) prefactor. 
      !  Include the same for the outward-integrated solution.
      !  Don't forget to correct the derivative - we scaled it by the grid spacing!
      !
      vs  = sqrt(eta_max)*continuum(ic)%wf(0,eta_npts) 
      gs  = sqrt(eta_max)*continuum(ic)%wf(1,eta_npts)/eta_tab(2) + 0.5_rk*continuum(ic)%wf(0,eta_npts)/sqrt(eta_max)
      !
      !  Match the solutions
      !
      continuum(ic)%c_out = (gs*vb(ic) - vs*gb(ic))/wrn(ic)
      continuum(ic)%c_in  = (gs*va(ic) - vs*ga(ic))/wrn(ic)
    end do match_channels
    !
    if (verbose>=1) then
      write (out,"(/t5,'Asymptotic amplitudes'/)")
      write (out,"((1x,a5,3(2x,a24,1x,a24)))") &
            'Chan', 'Re[in/out]', 'Im[in/out]', 'Re[outgoing]', 'Im[outgoing]', 'Re[incoming]', 'Im[incoming]', &
            '----', '----------', '----------', '------------', '------------', '------------', '------------'
      report_composition: do ic=1,n_channels
        write (out,"(1x,i5,3(2x,g24.14e3,1x,g24.14e3))") &
               ic, continuum(ic)%c_in/continuum(ic)%c_out, continuum(ic)%c_out, continuum(ic)%c_in
      end do report_composition
      call flush_wrapper (out)
    end if
    call TimerStop('Continuum outbound')
  end subroutine solve_continuum_integrate_outbound
  !
  !  Fixing the inner boundary amounts to solving a generalized eigenequation:
  !
  !   Ain C = lambda Aout C
  !
  !  where Ain (Aout) is the matrix giving the incoming (outgoing) amplitudes as 
  !  a function of the value at the inner boundary. We'll build the Ain/Aout
  !  matrices column-by-column, by repreatedly integrating the SE with unit boundary 
  !  in each of the channels.
  !
  subroutine solve_continuum_guess_inner_boundary_asymptotic(verbose,energy,va,ga,vb,gb,wrn)
    integer(ik), intent(in)           :: verbose
    complex(rk), intent(in)           :: energy
    complex(rk), intent(in)           :: va(n_channels)     ! Outgoing solutions at the matching point
    complex(rk), intent(in)           :: ga(n_channels)     ! Gradients of the outgoing solutions
    complex(rk), intent(in)           :: vb(n_channels)     ! Incoming solutions at the matching point
    complex(rk), intent(in)           :: gb(n_channels)     ! Gradients of the incoming solutions
    complex(rk), intent(in)           :: wrn(n_channels)    ! The Wronskian of the incoming/outgoing solution pair
    !
    complex(rk), allocatable :: a_in(:,:), a_out(:,:), a_in_save(:,:), a_out_save(:,:), alpha(:), beta(:)
    complex(rk), allocatable :: amp_in(:,:), amp_out(:,:) ! Our guess at the incoming/outgoing amplitudes for the solutions
    real(rk), allocatable    :: keys(:)                   ! Keys for sorting
    complex(rk)              :: scl
    integer(ik)              :: ichan, jchan, ib
    integer(ik)              :: alloc
    integer(ik)              :: order(n_channels)
    real(rk)                 :: delta_phase(n_channels), ref_phase
    !
    allocate (a_in     (n_channels,n_channels),a_out     (n_channels,n_channels), &
              a_in_save(n_channels,n_channels),a_out_save(n_channels,n_channels), &
              amp_in   (n_channels,n_channels),amp_out   (n_channels,n_channels), &
              alpha(n_channels),beta(n_channels),keys(n_channels),stat=alloc)
    if (alloc/=0) then
      stop 'general_tunnel%solve_continuum_fix_inner_boundary_asymptotic - allocation failed'
    end if
    !
    build_eigenproblem: do jchan=1,n_channels
      inner_boundary(:)     = 0
      inner_boundary(jchan) = 1
      call solve_continuum_integrate_outbound(verbose-1,energy,va,ga,vb,gb,wrn,inner_boundary,continuum)
      fill_columns: do ichan=1,n_channels
        a_in (ichan,jchan) = continuum(ichan)%c_in
        a_out(ichan,jchan) = continuum(ichan)%c_out
      end do fill_columns
    end do build_eigenproblem
    !
    !  ggev() overwrites the matrices; we need them to get estimations of the fluxes for each
    !         of the solutions
    !
    a_in_save  = a_in
    a_out_save = a_out
    call lapack_ggev(a_in,a_out,alpha,beta)
    !
    !  Normalize eigenvalues and adjust the phase of the dominant channel to make its
    !  coefficient real.
    !
    normalize: do jchan=1,n_channels
      scl = 1.0_rk / sqrt(sum(abs(a_out(:,jchan))**2))
      scl = scl * (abs(a_out(main_channel,jchan))/a_out(main_channel,jchan))
      a_out(:,jchan) = scl * a_out(:,jchan)
    end do normalize
    !
    !  Get an estimation of the incoming and outgoing amplitudes
    !
    amp_in  = matmul(a_in_save, a_out)
    amp_out = matmul(a_out_save,a_out)
    !
    !  Choosing the desired solution may be a little tricky.
    !
    select case (boundary_solution)
      case default
        write (out,"('Unrecognized boundary_solution = ''',a,'''')") trim(boundary_solution)
        stop 'general_tunnel%solve_continuum_fix_inner_boundary_asymptotic - unknonw boundary_solution'
      case ('magnitude')
        if (boundary_index<1 .or. boundary_index>n_channels) then
          write (out,"('boundary_index = ',i0,' is not between 1 and ',i0)") boundary_index, n_channels
          stop 'general_tunnel%solve_continuum_fix_inner_boundary_asymptotic - bad boundary index'
        end if
        keys(:) = abs(alpha(:)/beta(:))
        call order_keys(keys,order)
        ichan = order(boundary_index)
      case ('incoming')
        if (boundary_index<1 .or. boundary_index>n_channels) then
          write (out,"('boundary_index = ',i0,' is not between 1 and ',i0)") boundary_index, n_channels
          stop 'general_tunnel%solve_continuum_fix_inner_boundary_asymptotic - bad boundary index'
        end if
        keys(:) = sum(abs(amp_in(:,:))**2,dim=1)
        call order_keys(keys,order)
        ichan = order(boundary_index)
      case ('outgoing')
        if (boundary_index<1 .or. boundary_index>n_channels) then
          write (out,"('boundary_index = ',i0,' is not between 1 and ',i0)") boundary_index, n_channels
          stop 'general_tunnel%solve_continuum_fix_inner_boundary_asymptotic - bad boundary index'
        end if
        keys(:) = sum(abs(amp_out(:,:))**2,dim=1)
        call order_keys(keys,order)
        ichan = order(boundary_index)
      case ('phase')
        ref_phase      = atan2(aimag(boundary_phase),  real(boundary_phase,kind=rk))
        delta_phase(:) = atan2(aimag(alpha(:)/beta(:)),real(alpha(:)/beta(:),kind=rk)) - ref_phase
        where (delta_phase> pi) ; delta_phase = delta_phase - 2._rk * pi ; end where
        where (delta_phase<-pi) ; delta_phase = delta_phase + 2._rk * pi ; end where
        ichan = minloc(abs(delta_phase),dim=1)
    end select
    inner_boundary(:) = a_out(:,ichan)
    !
    if (verbose>=1) then
      write (out,"(/'Choosing solution # ',i0,' for the boundary'/)") ichan
    end if
    !
    report_boundaries: do ib=1,n_channels
      if (verbose>=2 .or. (verbose>=1 .and. ib==ichan)) then
         write (out,"(/t2,'Boundary choice ',i5,' eigenvalue = ',2(1x,g24.14e3),' abs = ',g24.14e3/)") &
                ib, alpha(ib)/beta(ib), abs(alpha(ib)/beta(ib))
         write (out,"( t2,'Guess eigenvector and amplitudes:'/)")
         write (out,"((1x,a5,3(2x,a38,1x,a38)))") &
                    'Channel', 'Re[Ci]', 'Im[Ci]', 'Re[In]', 'Im[In]', 'Re[out]', 'Im[out]', &
                    '-------', '------', '------', '------', '------', '-------', '-------'
         print_amplitudes: do ichan=1,n_channels
           write (out,"(1x,i5,3(2x,g38.24e3,1x,g38.24e3))") ichan, a_out(ichan,ib), amp_in(ichan,ib), amp_out(ichan,ib)
         end do print_amplitudes
         write (out,"()")
         call flush_wrapper(out)
      end if
    end do report_boundaries
    !
    deallocate (a_in,a_out,a_in_save,a_out_save,amp_in,amp_out,alpha,beta,keys)
  end subroutine solve_continuum_guess_inner_boundary_asymptotic
  !
  !  Once we have a decent initial guess for the eigenvector, we can improve it
  !  using general equation solver together with the shooting method (same as
  !  we did for the bound solutions).
  !
  subroutine solve_continuum_refine_inner_boundary_asymptotic(verbose,energy,va,ga,vb,gb,wrn)
    integer(ik), intent(in)           :: verbose
    complex(rk), intent(in)           :: energy
    complex(rk), intent(in)           :: va(n_channels)     ! Outgoing solutions at the matching point
    complex(rk), intent(in)           :: ga(n_channels)     ! Gradients of the outgoing solutions
    complex(rk), intent(in)           :: vb(n_channels)     ! Incoming solutions at the matching point
    complex(rk), intent(in)           :: gb(n_channels)     ! Gradients of the incoming solutions
    complex(rk), intent(in)           :: wrn(n_channels)    ! The Wronskian of the incoming/outgoing solution pair
    !
    integer(ik)    :: ichan
    type(fr_state) :: root
    real(rk)       :: vars_init(2*n_channels-2)
    !
    call extract_guess ! Fills vars_init from inner_boundary
    call fr_initialize(verbose-10,root,neqs=2*n_channels-2,vars=vars_init,epsstep=boundary_tol,diffstep=boundary_tol)
    root_iterations: do
      call fr_step(root)
      select case (root%action)
        case default
          write (out,"('general_tunnel%solve_continuum_refine_inner_boundary_asymptotic:"// &
                     " unexpected request from equation solver: ',i0)") &
                 root%action
          stop 'general_tunnel%solve_continuum_refine_inner_boundary_asymptotic - state machine corrupted'
        case (fr_act_evaluate)
          call evaluate_function
        case (fr_act_done)
          if (verbose>=1) then
            write (out,"('Optimization of boundary amplitudes converged after ',i0,' iterations')") root%nr_iter
          end if
          exit root_iterations
        case (fr_act_failed)
          write (out,"('WARNING: Optimization of boundary amplitudes failed to converge.')")
          write (out,"('WARNING: Continuing with possibly unconverged values.')")
          call flush_wrapper (out)
          exit root_iterations
      end select
    end do root_iterations
    call evaluate_function  ! Updates inner_boundary
    call fr_finalize(root)
    !
    if (verbose>=0) then
       write (out,"(/t2,'Refined boundary conditions:'/)")
       write (out,"((1x,a5,1x,2(1x,a38)))") 'Channel', 'Re[Ci]', 'Im[Ci]', '-------', '------', '------'
       print_amplitudes: do ichan=1,n_channels
         write (out,"(1x,i5,1x,2(1x,g38.24e3))") ichan, inner_boundary(ichan)
       end do print_amplitudes
       write (out,"()")
       call flush_wrapper(out)
    end if
    !
    contains
    !
    subroutine extract_guess
      integer(ik)    :: ichan, ivar
      complex(rk)    :: scl
      !
      scl  = 1._rk / inner_boundary(main_channel)
      ivar = 1
      copy_guess: do ichan=1,n_channels
        if (ichan==main_channel) cycle copy_guess
        vars_init(ivar+0) =  real(scl*inner_boundary(ichan),kind=rk)
        vars_init(ivar+1) = aimag(scl*inner_boundary(ichan))
        ivar = ivar + 2
      end do copy_guess
    end subroutine extract_guess
    !
    subroutine evaluate_function
      integer(ik) :: ichan, ivar
      real(rk)    :: scl
      complex(rk) :: a_in(n_channels), a_out(n_channels) ! Incoming/outgoing amplitudes in each channel
      complex(rk) :: ref_ratio, ratio, error             ! incomong/outgoing ratios
      real(rk)    :: ref_weight, weight
      !
      ivar = 1
      copy_boundary: do ichan=1,n_channels
        if (ichan==main_channel) then
          inner_boundary(ichan) = 1._rk
        else
          inner_boundary(ichan) = cmplx(root%vars(ivar+0),root%vars(ivar+1),kind=rk)
          ivar = ivar + 2
        end if
      end do copy_boundary
      scl = 1._rk / sqrt(sum(abs(inner_boundary(:))**2))
      inner_boundary(:) = scl * inner_boundary(:)
      !
      call solve_continuum_integrate_outbound(verbose-1,energy,va,ga,vb,gb,wrn,inner_boundary,continuum)
      !
      extract_amplitudes: do ichan=1,n_channels
        a_in (ichan) = continuum(ichan)%c_in
        a_out(ichan) = continuum(ichan)%c_out
      end do extract_amplitudes
      !
      ref_ratio  = a_in(main_channel)/a_out(main_channel)
      ref_weight = sqrt(abs(a_in(main_channel))**2 + abs(a_out(main_channel))**2)
      !
      ivar = 1
      set_equations: do ichan=1,n_channels
        if (ichan==main_channel) cycle set_equations
        !
        ratio  = a_in(ichan)/a_out(ichan)
        weight = sqrt(abs(a_in(ichan))**2 + abs(a_out(ichan))**2)
        error  = (weight/ref_weight) * (ratio - ref_ratio)
        root%eqs(ivar+0) =  real(error,kind=rk)
        root%eqs(ivar+1) = aimag(error)
        !
        ivar = ivar + 2
      end do set_equations
    end subroutine evaluate_function
  end subroutine solve_continuum_refine_inner_boundary_asymptotic
  !
  !  Making the flux stationary requires solving a Hermitian eigenproblem
  !
  !   A^H A C = lambda C
  !
  !  where A is the matrix giving the incoming (or outgoing) amplitudes as 
  !  a function of the value at the inner boundary. 
  !
  !  The subroutine below is a transplant of solve_continuum_guess_inner_boundary_asymptotic()
  !
  subroutine solve_continuum_guess_inner_boundary_flux(verbose,energy,va,ga,vb,gb,wrn)
    integer(ik), intent(in)           :: verbose
    complex(rk), intent(in)           :: energy
    complex(rk), intent(in)           :: va(n_channels)     ! Outgoing solutions at the matching point
    complex(rk), intent(in)           :: ga(n_channels)     ! Gradients of the outgoing solutions
    complex(rk), intent(in)           :: vb(n_channels)     ! Incoming solutions at the matching point
    complex(rk), intent(in)           :: gb(n_channels)     ! Gradients of the incoming solutions
    complex(rk), intent(in)           :: wrn(n_channels)    ! The Wronskian of the incoming/outgoing solution pair
    !
    complex(rk), allocatable :: a_in(:,:), a_out(:,:), aha(:,:), amp_in(:,:), amp_out(:,:)
    real(rk), allocatable    :: lambda(:)
    complex(rk)              :: scl
    integer(ik)              :: ichan, jchan, ib
    integer(ik)              :: alloc
    !
    allocate (a_in   (n_channels,n_channels), &
              a_out  (n_channels,n_channels), &
              aha    (n_channels,n_channels), &
              amp_in (n_channels,n_channels), &
              amp_out(n_channels,n_channels), &
              lambda (n_channels),stat=alloc)
    if (alloc/=0) then
      stop 'general_tunnel%solve_continuum_fix_inner_boundary_flux - allocation failed'
    end if
    !
    build_eigenproblem: do jchan=1,n_channels
      inner_boundary(:)     = 0
      inner_boundary(jchan) = 1
      call solve_continuum_integrate_outbound(verbose-1,energy,va,ga,vb,gb,wrn,inner_boundary,continuum)
      fill_columns: do ichan=1,n_channels
        a_in (ichan,jchan) = continuum(ichan)%c_in
        a_out(ichan,jchan) = continuum(ichan)%c_out
      end do fill_columns
    end do build_eigenproblem
    !
    select case (boundary_solution)
      case default
        write (out,"('boundary_solution=''',a,''' is not valid for boundary=''flux''')") trim(boundary_solution)
        stop 'general_tunnel%solve_continuum_fix_inner_boundary_flux - bad boundary_solution'
      case ('incoming')
        aha = matmul(conjg(transpose(a_in )),a_in )
      case ('outgoing')
        aha = matmul(conjg(transpose(a_out)),a_out)
    end select
    aha = 0.5_rk * (aha + conjg(transpose(aha)))  ! Enforce hermiticity
    call lapack_heev(aha,lambda)
    !
    !  Normalize eigenvalues and adjust the phase of the dominant channel to make its
    !  coefficient real.
    !
    normalize: do jchan=1,n_channels
      scl = 1.0_rk / sqrt(sum(abs(aha(:,jchan))**2))
      scl = scl * (abs(aha(main_channel,jchan))/aha(main_channel,jchan))
      aha(:,jchan) = scl * aha(:,jchan)
    end do normalize
    !
    !  Get an estimation of the incoming and outgoing amplitudes
    !
    amp_in  = matmul(a_in, aha)
    amp_out = matmul(a_out,aha)
    !
    if (boundary_index<1 .or. boundary_index>n_channels) then
      write (out,"('boundary_index = ',i0,' is not between 1 and ',i0)") boundary_index, n_channels
      stop 'general_tunnel%solve_continuum_fix_inner_boundary_flux - bad boundary index'
    end if
    inner_boundary(:) = aha(:,boundary_index)
    !
    report_boundaries: do ib=1,n_channels
      if (verbose>=2 .or. (verbose>=1 .and. ib==boundary_index)) then
         write (out,"(/t2,'Boundary choice ',i5,' eigenvalue = ',(1x,g24.14e3)/)") ib, lambda(ib)
         write (out,"( t2,'Guess eigenvector and amplitudes:'/)")
         write (out,"((1x,a5,3(2x,a38,1x,a38)))") &
                    'Channel', 'Re[Ci]', 'Im[Ci]', 'Re[In]', 'Im[In]', 'Re[out]', 'Im[out]', &
                    '-------', '------', '------', '------', '------', '-------', '-------'
         print_amplitudes: do ichan=1,n_channels
           write (out,"(1x,i5,3(2x,g38.24e3,1x,g38.24e3))") ichan, aha(ichan,ib), amp_in(ichan,ib), amp_out(ichan,ib)
         end do print_amplitudes
         write (out,"(1x,a5,3(2x,g38.24e3,1x,38x))") &
                'Norm', sum(abs(aha(:,ib))**2), sum(abs(amp_in(:,ib))**2), sum(abs(amp_out(:,ib))**2)
         write (out,"()")
         call flush_wrapper(out)
      end if
    end do report_boundaries
    !
    deallocate (a_in,a_out,aha,amp_in,amp_out,lambda)
  end subroutine solve_continuum_guess_inner_boundary_flux
  !
  !  Refine the boundary conditions, using direct minimization of the flux
  !  magnitude. This approach only works when boundary_index=1 (which, luckily,
  !  is where the loss of accuracy in the direct solution is most severe)
  !
  subroutine solve_continuum_refine_inner_boundary_flux(verbose,energy,va,ga,vb,gb,wrn)
    integer(ik), intent(in)           :: verbose
    complex(rk), intent(in)           :: energy
    complex(rk), intent(in)           :: va(n_channels)     ! Outgoing solutions at the matching point
    complex(rk), intent(in)           :: ga(n_channels)     ! Gradients of the outgoing solutions
    complex(rk), intent(in)           :: vb(n_channels)     ! Incoming solutions at the matching point
    complex(rk), intent(in)           :: gb(n_channels)     ! Gradients of the incoming solutions
    complex(rk), intent(in)           :: wrn(n_channels)    ! The Wronskian of the incoming/outgoing solution pair
    !
    integer(ik)    :: ichan
    type(fm_state) :: root
    real(rk)       :: vars_init(2*n_channels-2)
    !
    call extract_guess ! Fills vars_init from inner_boundary
    call fm_initialize(verbose-10,root,vars=vars_init,epsstep=boundary_tol,diffstep=boundary_tol)
    root_iterations: do
      call fm_step(root)
      select case (root%action)
        case default
          write (out,"('general_tunnel%solve_continuum_refine_inner_boundary_flux: unexpected request from equation solver: ',i0)")&
                 root%action
          stop 'general_tunnel%solve_continuum_refine_inner_boundary_flux - state machine corrupted'
        case (fm_act_evaluate)
          call evaluate_function
        case (fm_act_done)
          if (verbose>=1) then
            write (out,"('Optimization of boundary amplitudes converged.')")
          end if
          exit root_iterations
        case (fm_act_failed)
          write (out,"('WARNING: Optimization of boundary amplitudes failed to converge.')")
          write (out,"('WARNING: Continuing with possibly unconverged values.')")
          call flush_wrapper (out)
          exit root_iterations
      end select
    end do root_iterations
    call evaluate_function  ! Updates inner_boundary
    call fm_finalize(root)
    !
    if (verbose>=0) then
       write (out,"(/t2,'Refined boundary conditions:'/)")
       write (out,"((1x,a5,1x,2(1x,a38)))") 'Channel', 'Re[Ci]', 'Im[Ci]', '-------', '------', '------'
       print_amplitudes: do ichan=1,n_channels
         write (out,"(1x,i5,1x,2(1x,g38.24e3))") ichan, inner_boundary(ichan)
       end do print_amplitudes
       write (out,"()")
       call flush_wrapper(out)
    end if
    !
    contains
    !
    subroutine extract_guess
      integer(ik)    :: ichan, ivar
      complex(rk)    :: scl
      !
      scl  = 1._rk / inner_boundary(main_channel)
      ivar = 1
      copy_guess: do ichan=1,n_channels
        if (ichan==main_channel) cycle copy_guess
        vars_init(ivar+0) =  real(scl*inner_boundary(ichan),kind=rk)
        vars_init(ivar+1) = aimag(scl*inner_boundary(ichan))
        ivar = ivar + 2
      end do copy_guess
    end subroutine extract_guess
    !
    subroutine evaluate_function
      integer(ik) :: ichan, ivar
      real(rk)    :: scl
      complex(rk) :: a_in(n_channels), a_out(n_channels) ! Incoming/outgoing amplitudes in each channel
      !
      ivar = 1
      copy_boundary: do ichan=1,n_channels
        if (ichan==main_channel) then
          inner_boundary(ichan) = 1._rk
        else
          inner_boundary(ichan) = cmplx(root%vars(ivar+0),root%vars(ivar+1),kind=rk)
          ivar = ivar + 2
        end if
      end do copy_boundary
      scl = 1._rk / sqrt(sum(abs(inner_boundary(:))**2))
      inner_boundary(:) = scl * inner_boundary(:)
      !
      call solve_continuum_integrate_outbound(verbose-1,energy,va,ga,vb,gb,wrn,inner_boundary,continuum)
      !
      extract_amplitudes: do ichan=1,n_channels
        a_in (ichan) = continuum(ichan)%c_in
        a_out(ichan) = continuum(ichan)%c_out
      end do extract_amplitudes
      !
      select case (boundary_solution)
        case default
          stop 'general_tunnel%solve_continuum_refine_inner_boundary_flux - bad boundary_solution'
        case ('incoming')
          root%func = sum(abs(a_in)**2)
        case ('outgoing')
          root%func = sum(abs(a_out)**2)
      end select
      if (verbose>=2) then
        write (out,"(/'Boundary refinement step: ',g48.36e3)") root%func
        call flush_wrapper(out)
      end if
    end subroutine evaluate_function
  end subroutine solve_continuum_refine_inner_boundary_flux
  !
  subroutine solve_continuum(verbose,energy)
    integer(ik), intent(in)           :: verbose
    complex(rk), intent(in)           :: energy
    !
    integer(ik) :: ichan
    complex(rk) :: va(n_channels)        ! Outgoing solutions at the matching point
    complex(rk) :: ga(n_channels)        ! Gradients of the outgoing solutions
    complex(rk) :: vb(n_channels)        ! Incoming solutions at the matching point
    complex(rk) :: gb(n_channels)        ! Gradients of the incoming solutions
    complex(rk) :: wrn(n_channels)       ! The Wronskian of the incoming/outgoing solution pair
    !
    call TimerStart('Continuum direction')
    !
    !  Prepare the asymptotic solutions. We'll keep reusing them, so there is no
    !  sense in recomputing them all the time.
    !
    call solve_continuum_asymptotic_solutions(verbose,energy,va,ga,vb,gb,wrn)
    !
    boundaries: select case (boundary)
      case default
        ! Boundary has been fixed at input, nothing to do
      case ('asymptotic')
        if (n_channels==1) then
          inner_boundary(1) = 1
          exit boundaries
        end if
        !
        !  We want to make the solution factorizable in the asymtotic region.
        !  Finding the necessary boundary conditions at the origin requires some extra work.
        !
        call solve_continuum_guess_inner_boundary_asymptotic(verbose,energy,va,ga,vb,gb,wrn)
        !
        !  Unfortunately, LAPACK routines (which we use above) have very poor
        !  accuracy for the eigenvectors of general complex eigenproblems, even
        !  when the problem is well-conditioned. We will still need to refine the 
        !  eigenvector we just determined if we want full numerical accuracy.
        !
        call solve_continuum_refine_inner_boundary_asymptotic(verbose,energy,va,ga,vb,gb,wrn)
      case ('flux')
        if (n_channels==1) then
          inner_boundary(1) = 1
          exit boundaries
        end if
        !
        !  We want to make the total flux stationary
        !
        call solve_continuum_guess_inner_boundary_flux(verbose,energy,va,ga,vb,gb,wrn)
        !
        !  Unfortunately, we had to solve an eigenproblem for a square dependency matrix.
        !  This eigenproblem has an awfully bad condition number, and our boundary is
        !  not very accurate. We'll need to refine it before calculating the final
        !  solution by directly minimizing it. This can only be done for boundary_index=1
        !
        if (boundary_index/=1) then
          if (verbose>=0) then
            write (out,"('WARNING: Solution refinement for boundary=''flux'' and boundary_index>1 is not implenmented')")
            write (out,"('WARNING: Expect decreased accuracy in the final solutions')")
          end if
          exit boundaries
        end if
        call solve_continuum_refine_inner_boundary_flux(verbose,energy,va,ga,vb,gb,wrn)
    end select boundaries
    !
    !  The inner boundary is fixed; calculate the continuum part of the wavefunction
    !
    call solve_continuum_integrate_outbound(verbose,energy,va,ga,vb,gb,wrn,inner_boundary,continuum)
    !
    if (verbose>=1) then
      report_continuum_wavefunctions: do ichan=1,n_channels
        call dump_continuum_wavefunction(ichan,continuum(ichan))
      end do report_continuum_wavefunctions
    end if
    !
    call TimerStop('Continuum direction')
  end subroutine solve_continuum
  !
  !  (re-)Calculate continuum wavefunction by integrating backwards from the matching point.
  !
  subroutine invert_continuum(verbose,energy)
    integer(ik), intent(in)           :: verbose
    complex(rk), intent(in)           :: energy
    !
    integer(ik) :: ichan
    complex(rk) :: va(n_channels)        ! Outgoing solutions at the matching point
    complex(rk) :: ga(n_channels)        ! Gradients of the outgoing solutions
    complex(rk) :: vb(n_channels)        ! Incoming solutions at the matching point
    complex(rk) :: gb(n_channels)        ! Gradients of the incoming solutions
    complex(rk) :: wrn(n_channels)       ! The Wronskian of the incoming/outgoing solution pair
    !
    call TimerStart('Reverse continuum')
    !
    !  Prepare the asymptotic solutions. We need these soltions to construct the
    !  outer boundary conditions. Note that our asymptotic solutions need to be
    !  divided by sqrt(eta_max). Furthermore, the gradient must be scaled by the
    !  grid spacing.
    !
    call solve_continuum_asymptotic_solutions(verbose,energy,va,ga,vb,gb,wrn)
    ! 
    !  eta_tab(2) is the grid spacing
    !
    ga = eta_tab(2) * ( ga / sqrt(eta_max) - 0.5_rk * va / eta_max**1.5_rk )
    gb = eta_tab(2) * ( gb / sqrt(eta_max) - 0.5_rk * vb / eta_max**1.5_rk )
    va = va / sqrt(eta_max)
    vb = vb / sqrt(eta_max)
    !
    select case (file_total_mode)
      case default
        write (out,"('invert_continuum: Can''t deal with file_total_mode = ',a)") trim(file_total_mode)
        stop 'general_tunnel_continuum%invert_continuum - bad file_total_mode'
      case ('reverse incoming')
        copy_incoming: do ichan=1,n_channels
          pseudo(ichan)%wf(0,eta_npts) = continuum(ichan)%c_in * vb(ichan)
          pseudo(ichan)%wf(1,eta_npts) = continuum(ichan)%c_in * gb(ichan)
        end do copy_incoming
      case ('reverse outgoing')
        copy_outgoing: do ichan=1,n_channels
          pseudo(ichan)%wf(0,eta_npts) = continuum(ichan)%c_out * va(ichan)
          pseudo(ichan)%wf(1,eta_npts) = continuum(ichan)%c_out * ga(ichan)
        end do copy_outgoing
      case ('reverse total')
        copy_total: do ichan=1,n_channels
          pseudo(ichan)%wf(0:1,eta_npts) = continuum(ichan)%wf(0:1,eta_npts)
        end do copy_total
    end select
    !
    !  The outer boundary is fixed; calculate the continuum part of the wavefunction
    !
    call integrate_radial_continuum_inbound(pseudo,energy,efield,chan_alp,chan_w1,chan_w2)
    !
    if (verbose>=1) then
      report_continuum_wavefunctions: do ichan=1,n_channels
        call dump_continuum_wavefunction(ichan,pseudo(ichan))
      end do report_continuum_wavefunctions
    end if
    !
    call TimerStop('Reverse continuum')
  end subroutine invert_continuum
end module general_tunnel_continuum

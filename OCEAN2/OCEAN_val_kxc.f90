! Copyright (C) 2020 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!


module OCEAN_val_kxc
  use AI_kinds, only : DP


  implicit none
  private
  save

  public :: OCEAN_val_kxc_init

  contains


  subroutine OCEAN_val_kxc_init( sys, VR, ierr )
    use OCEAN_system, only : o_system
    !
    type(o_system), intent( in ) :: sys
    real(DP), intent( out ) :: VR(:)
    integer, intent( inout ) :: ierr
    !
    real(DP), allocatable :: rho( : )


    allocate( rho( sys%nxpts ), STAT=ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_get_rho( sys%xmesh, sys%celvol, rho, ierr )
    if( ierr .ne. 0 ) return

    call alda( sys, rho, VR, ierr )

  end subroutine OCEAN_val_kxc_init



  
  subroutine alda( sys, rho, VR, ierr )
    use OCEAN_system, only : o_system
    use OCEAN_val_states, only : nxpts, startx
    use OCEAN_mpi, only : myid, root, comm, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
    type(o_system), intent( in ) :: sys
    real(DP), intent( in ) :: rho(:)
    real(DP), intent( out ) :: VR(:)
    integer, intent( inout ) :: ierr


    real(DP) :: prefac, localDen, n, dn, ex, ec, val(-2:2), totalDen, totalDen2, ux1, ux2, uc1, uc2 
    integer :: i, j

    prefac = 1.0_DP / ( sys%nkpts * sys%celvol )

    totalDen = 0.0_DP
    totalDen2 = sum( rho(:) ) / real(sys%nxpts, DP ) * sys%celvol
    write(6,*) nxpts, sys%nxpts

    if( myid .eq. root ) then
      open(unit=99,file='kxc.txt' )
    endif

    do i = 1, nxpts
      localDen = rho( i + startx - 1 )
      totalDen = totalDen + localDen
      
      dn = 0.01_DP * localDen
      do j = -2, 2
        n = localDen + real( j, DP ) * dn

!        call mod_cacorr( n, ex, ec )
        call cacorr( n, ex, ec, ux1, ux2, uc1, uc2 )
        val( j ) = n * ( ex + ec )
      enddo

      VR( i ) = ( 16.0d0 * ( val(-1) + val(1) ) - ( val(-2) + val(2) ) - 30.0d0 * val(0) ) & 
              / ( 12.0d0 * dn ** 2 )
      if( myid .eq. root ) write(99,*) localDen, VR(i)
    enddo
    VR(:) = VR(:) * prefac
    call MPI_ALLREDUCE( MPI_IN_PLACE, totalDen, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
    totalDen = totalDen  / real(sys%nxpts, DP ) * sys%celvol
    if( myid .eq. root ) write( 6,*) 'Total density: ', totalDen, totalDen2

    if( myid .eq. root ) then
      close(99)
    endif

  end subroutine alda

!  exchange correlation routine, by Ceperley-Alder, as parametrized by
!  Perdew and Zunger, Phys. Rev. B 23, 5048.  we use their interpolation
!  between the unpolarized and polarized gas for the correlation part.
!
  subroutine mod_cacorr( xn, ex, ec )
    use OCEAN_constants, only : PI_DP
    real(DP), intent( in ) :: xn
    real(DP), intent( out ) :: ex, ec

    real(DP) :: ux1,ux2,uc1,uc2
    real(DP) :: xn1, xn2, rs, zeta, exchfactor, gamma, beta1, beta2, denom, ecu, ucu, ecp, ucp, & 
                xlr, rlr, au, bu, cu, du, ap, bp, cp, f, dfdz, ddp, rootr
    real(DP) :: fe1, fu1, befactor, fp, b2, eta, xl, ex1, ex2, fe2, fu2, beta

    real(DP), parameter :: trd=1.0_DP/3.0_dp
    real(DP), parameter :: ft=4.0_dp/3.0_dp
!    real(DP), parameter :: rel = 0.

!
!  get the n's, and the rs.
!
      fp=4.d0*pi_dp
      xn1=xn/2
      xn2=xn/2

!  effect cutoff, to avoid overflow

      if ( xn .lt. 0.00000001d0 ) then

        ex=0.d0
        ec=0.d0
        ux1=0.d0
        ux2=0.d0
        uc1=0.d0      
        uc2=0.d0

      else

        rs=(3.d0/(fp*xn))**trd
        zeta=(xn1-xn2)/xn
!       exchfactor=-0.930525736d0
        exchfactor=-1.5d0*(0.75d0/pi_dp)**trd

        befactor=(9.d0*pi_dp/4.d0)**trd/137.03599976d0
        if (xn1.eq.0.d0) then
          fe1=1.d0
          fu1=1.d0
          ex1=0.d0
          ux1=0.d0
        else
          beta=befactor/rs
          b2=beta*beta
          eta=dsqrt(1.d0+b2)
          xl=dlog(beta+eta)
          fe1=1.d0-1.5d0*((beta*eta-xl)/b2)**2.d0
          fu1=-0.5d0+1.5d0*xl/beta/eta
          ex1=exchfactor*xn1**trd
          ux1=4.d0*ex1/3.d0
        endif
        if (xn2.eq.0.d0) then
          fe2=1.d0
          fu2=1.d0
          ex2=0.d0
          ux2=0.d0
        else
          beta=befactor/rs
          b2=beta*beta
          eta=dsqrt(1.d0+b2)
          xl=dlog(beta+eta)
          fe2=1.d0-1.5d0*((beta*eta-xl)/b2)**2.d0
          fu2=-0.5d0+1.5d0*xl/beta/eta
          ex2=exchfactor*xn2**trd
          ux2=4.d0*ex2/3.d0
        endif
!  these next lines do the Ceperley-Alder correlation
        if (rs.ge.1.d0) then

          rootr=dsqrt(rs)

          gamma=-0.1423d0
          beta1=1.0529d0
          beta2=0.3334d0
          denom=(1.d0+beta1*rootr+beta2*rs)
          ecu=gamma/denom
          ucu=ecu*(1.d0+7.d0/6.d0*beta1*rootr+ft*beta2*rs)/denom

          gamma=-0.0843d0
          beta1=1.3981d0
          beta2=0.2611d0
          denom=(1.d0+beta1*rootr+beta2*rs)
          ecp=gamma/denom
          ucp=ecp*(1.d0+7.d0/6.d0*beta1*rootr+ft*beta2*rs)/denom

        else

          xlr=dlog(rs)
          rlr=rs*xlr

          au= 0.0311d0
          bu=-0.048d0
          cu= 0.002d0
          du=-0.0116d0
          ecu=au*xlr+bu+cu*rlr+du*rs
          ucu=au*xlr+(bu-au/3.d0)+2.d0/3.d0*cu*rlr+(2.d0*du-cu)*rs/3.d0

          ap= 0.01555d0
          bp=-0.0269d0
          cp= 0.0007d0
          ddp=-0.0048d0
          ecp=ap*xlr+bp+cp*rlr+ddp*rs
          ucp=ap*xlr+(bp-ap/3.d0)+2.d0/3.d0*cp*rlr+(2.d0*ddp-cp)*rs/3.d0

        endif

!  if we are nonrel, turn off the MacDonald-Vosko correction.
!
!        if (rel.eq.0.d0) then
          fe1=1.d0
          fu1=1.d0
          fe2=1.d0
          fu2=1.d0
!        endif

!  interpolate the correlation energies.

        denom=2.d0**ft-2.d0
        f=((1.d0+zeta)**ft+(1.d0-zeta)**ft-2.d0)/denom
        dfdz=ft/denom*((1.d0+zeta)**trd-(1.d0-zeta)**trd)
        ec=ecu+f*(ecp-ecu)
        uc1=ucu+f*(ucp-ucu)+(ecp-ecu)*(1.d0-zeta)*dfdz
        uc2=ucu+f*(ucp-ucu)-(ecp-ecu)*(1.d0+zeta)*dfdz        
!
!  get the final functional and potential.
!
        ex=(xn1*fe1*ex1+xn2*fe2*ex2)/xn
        ux1=fu1*ux1
        ux2=fu2*ux2
        uc1=uc1
        uc2=uc2
      endif
!
      return
      end subroutine mod_cacorr
    

end module


!
subroutine radial_exciton( iunit, rsDelta, rsMax, specpnt, targ_loc, inverseA, xmesh, Rstart, Rmesh, nspn, Rspace_exciton )
  use ocean_interpolate
  implicit none
  
  integer, parameter :: DP = kind(1.0d0 )
  integer, intent( in ) :: iunit
  real(DP), intent( in ) :: rsDelta
  real(DP), intent( in ) :: rsMax
  integer, intent( in ) :: specpnt
  real(DP), intent( in ) :: targ_loc(3)
  real(DP), intent( in ) :: inverseA(3,3)
  integer, intent( in ) :: xmesh(3), Rstart(3), Rmesh(3), nspn
  complex(DP), intent( in ) :: Rspace_exciton( xmesh(3)*Rmesh(3), xmesh(2)*Rmesh(2), xmesh(1)*Rmesh(1), nspn )

  real(DP), allocatable :: wgt(:), ang(:,:)
  integer :: nang, iang, nr, ir, xtarg(3), ixx, ixxx, iyy, iyyy, ispn
  real(DP) :: su, rs, rvec(3), temp1(3), distance(3)
  complex(DP) :: Pgrid(4), P(4,4), Qgrid(4), Q(4), R
  character(len=10) :: spectFileName

  write( spectFileName, '(a8,i0)' ) 'specpnt.', specpnt
  open(unit=100, file=spectFileName, form='formatted', status='old')
  read(100,*) nang
  allocate( ang(3,nang), wgt(nang) )
  su = 0.0_DP
  do iang = 1, nang
    read(100,*) ang(:,iang), wgt(iang)
    su = su + wgt(iang)
  enddo
  close(100)

  su = 1.0_DP/ su
  wgt(:) = wgt(:) * su

  write(6,*) 'su ', su, Rspace_exciton( 1,1,1,1)

  nr = ( (rsMax - rsDelta) / rsDelta ) + 1

  do ir = 1, nr
    rs = real( ir, DP ) * rsDelta
    su = 0.0_DP
    do iang = 1, nang
      rvec(:) = targ_loc(:) + rs * ang(:,iang )

      temp1 = matmul( inverseA, rvec )
      temp1(:) = temp1(:) * dble(xmesh(:)) - Rstart(:)*xmesh(:) + 1
      xtarg(:) = floor( temp1(:) )
      distance(:) = temp1(:) - xtarg(:)
!      write(6,*) xtarg

      do ispn = 1, nspn
        do ixx = xtarg(1)-1, xtarg(1)+2
          ixxx = ixx
          if( ixx .lt. 1 ) ixxx = 1
          if( ixx .gt. xmesh(1)*Rmesh(1) ) ixxx = xmesh(1)*Rmesh(1)
          do iyy = xtarg(2)-1, xtarg(2)+2
            iyyy = iyy
            if( iyy .lt. 1 ) iyyy = 1
            if( iyy .gt. xmesh(2)*Rmesh(2) ) iyyy = xmesh(2)*Rmesh(2)
            call MakeLagrange( 4, xtarg(3), iyyy, ixxx, Rspace_exciton(:,:,:,1), &
                               Pgrid(:) )
            P( iyy -xtarg(2) + 2, ixx-xtarg(1)+2 ) = evalLagrange( 4, distance(3), 1.0_DP, Pgrid(:) )
          enddo
        enddo

        do ixx = 1, 4
          call MakeLagrange( 4, P(:,ixx), Qgrid(:) )
          Q(ixx) = evalLagrange( 4, distance(2), 1.0_DP, Qgrid(:) )
        enddo
        call makeLagrange( 4, Q, Qgrid )
        R = evalLagrange( 4, distance(1), 1.0_DP, Qgrid )

        su = su + (R * conjg(R) ) * wgt( iang )
      enddo
    enddo
    write(99,*) rs, su
  enddo

  
end subroutine 
    

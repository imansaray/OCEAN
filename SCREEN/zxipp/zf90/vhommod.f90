! Copyright (C) 2015, 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program vhommod
  implicit none
  integer, parameter :: stdin = 5, stdout = 6
  integer, parameter :: dp = kind( 1.0d0 )
  !
  integer :: i, j, k, nr, nq, idum, nshells, nsites, isite, ierr, xmesh( 3 ), ids
  integer, allocatable :: siteIndex( : )
  real( kind=dp ) :: epsq0, qf, wpsqd, wf, lam, pi, eps
  real( kind=dp ) :: qq, q, dq, qmax, r1, r2, dr, rmax
  real( kind=dp ) :: jqr1, jqr2, tj1, tj2, th1, th2, j1, s
  real( kind=dp ) :: nav, ds, ss, dn, n, avec( 3, 3 ), omega, nel
  real( kind=dp ) :: d2, d3, rd3, q2, q3, q4, q5, delta, sinqr, sinqd, cosqd
  !
  real( kind=dp ), allocatable :: v( :, : ), vh( :, : ), dv( : )
  real( kind=dp ), allocatable :: tab1( :, : ), bq( : ), shells( : ), tabjqr2( :, : )
  !
  character(len=2), allocatable :: siteSymbol( : )
  character(len=4) :: smode
  character(len=80) :: avgName, reoptName
  !
  real( kind=dp ), external :: levlou, sphj0, sphj1
  !
  logical :: havenav, valenceGrid, thickShell, ex, doubleShell, solidShell
  character(len=80) :: dummy
  character(len=1), parameter :: delim_ = '/'
  !
  pi = 4.d0 * datan( 1.d0 )
  !
!  read ( stdin, * ) r2 !, dr, nr
!  read ( stdin, * ) dq, qmax   ! dq and qmax in units of qfermi
  !
  inquire(file='screen.doubleshell', exist=doubleShell )
  if( doubleShell ) then
    open(unit=99, file='screen.doubleshell', form='formatted', status='old' )
    read( 99, * ) doubleShell
    close(99)
  endif
  if( doubleShell ) then
    ids = 2
  else
    ids = 1
  endif

  open( unit=99, file='shells', form='formatted', status='old' )
  rewind 99
  read( 99, * ) nshells
  allocate( shells( nshells * ids ) )
  do i = 1, nshells
    read( 99, * ) shells( i )
  enddo
  close( 99 )

  if( doubleShell ) then
    do i = 1, nshells
      shells( i + nshells ) = 4.0_DP * shells( i ) / 5.0_dp
    enddo
  endif

  inquire(file='screen.mode', exist=valenceGrid )
  if( valenceGrid ) then
    open( unit=99, file='screen.mode', form='formatted', status='old' )
    read( 99, * ) smode
    close( 99 )
    if( smode .eq. 'grid' ) then
      valenceGrid = .true.
    else
      valenceGrid = .false.
    endif
  endif

  delta = 0.0_DP
  inquire( file='screen.shellwidth', exist=thickShell )
  if( thickShell ) then
    open( unit=99, file='screen.shellwidth', form='formatted', status='old' )
    read( 99, * ) delta
    close( 99 )
    if( delta .lt. 0.000001_dp .or. delta .gt. maxval( shells ) ) then
      delta = 0.0_dp
      thickShell = .false.
    endif
  endif

  if( delta .lt. 0.000001_dp ) then
    inquire( file='screen.solidshell', exist=solidShell )
    if( solidShell ) then
      open(99,file='screen.solidshell', form='formatted', status='old' )
      read( 99, * ) solidShell
      close(99)
    endif
  endif



  if( valenceGrid ) then
    open( unit=99, file='xmesh.ipt', form='formatted', status='old' )
    read( 99, * ) xmesh(:)
    close( 99 )
    nsites = product( xmesh( : ) )
    allocate( siteSymbol( 0 ), siteIndex( 0 ) )
  else  !core-level
    open( unit=99, file='sitelist', form='formatted', status='old' )
    read( 99, * ) nsites
    allocate( siteSymbol( nsites ), siteIndex( nsites ) )
    do i = 1, nsites
      read( 99, * ) siteSymbol( i ), idum, siteIndex( i )
    enddo
    close( 99 )
  endif

  open( unit=99, file='screen.final.rmax', form='formatted', status='old' )
  rewind 99
  read( 99, * ) rmax
  close( 99 )
  open( unit=99, file='screen.final.dr', form='formatted', status='old' )
  read( 99, * ) dr
  close( 99 )
  nr = rmax / dr

  open( unit=99, file='screen.model.dq', form='formatted', status='old' )
  rewind 99
  read( 99, * ) dq
  close( 99 )
  open( unit=99, file='screen.model.qmax', form='formatted', status='old' )
  rewind 99
  read( 99, * ) qmax
  close( 99 )
  open( unit=99, file='epsilon', form='formatted', status='old' )
  rewind 99
  read ( 99, * ) epsq0
  close( 99 )
  open( unit=99, file='avecsinbohr.ipt', form='formatted', status='old' )
  rewind 99
  read ( 99, * ) avec( :, : )
  close( unit=99 )
  call getomega( avec, omega )
  open ( unit=99, file='rhoofg', form='formatted', status='old' )
  rewind 99
  read ( 99, * ) idum
  read ( 99, * ) idum, idum, idum, nel
  close( unit=99 )
  inquire( file='fake_nav.ipt', exist=havenav)
  if( havenav ) then
    open( unit=99, file='fake_nav.ipt', form='formatted', status='old' )
    rewind 99
    read( 99, * ) nav
  else
    nav = dble( nel ) / omega
  endif
  !
!  read ( stdin, * ) epsq0 
  allocate( v( nr, nshells * ids ), vh( nr, nshells * ids ), dv( nr ) )
  !
  qf = ( 3 * pi ** 2 * nav ) ** ( 1.d0 / 3.d0 )
  wpsqd = 4 * pi * nav
  wf = qf ** 2 / 2
  lam = dsqrt( wpsqd / ( wf ** 2 * ( epsq0 - 1 ) ) )
  !
  nq = qmax / dq
  dq = qf * dq
  allocate( tab1( nr, nq ), bq( nq ), tabjqr2( nq, nshells * ids ) )
  !
  if( thickShell ) then
    d2 = delta * delta
    d3 = delta * d2
    do k = 1, nshells * ids
      r2 = shells( k )
      rd3 = ( r2 + delta ) ** 3
      if( solidShell ) then
        delta = shells( k ) / 2.0_DP
        r2 = shells( k ) / 2.0_DP
        d2 = delta * delta
        d3 = d2 * delta
        rd3 = shells( k ) ** 3
      endif
      do i = 1, nq
        q = dq * ( i - 0.5d0 )
        cosqd = cos( q * delta )
        sinqr = sin( q * r2 )
        sinqd = sin( q * delta )
        q2 = q * q
        q3 = q2 * q
        q4 = q2 * q2
        q5 = q3 * q2

        tabjqr2( i, k ) = -3.0_dp * cosqd * sinqr / ( q3 * rd3 ) &
                        - 90.0_dp * cosqd * sinqr / ( q5 * d2 * rd3 ) &
                        - 3.0_dp * r2**2 * cosqd * sinqr / ( q3 * d2 * rd3 ) &
                        - 270.0_dp * cosqd * sinqr / ( q5 * r2 * delta * rd3 ) &
                        - 9.0_DP * r2 * cosqd * sinqr / (q3 * delta * rd3 ) &
                        + 15.0_DP * delta * sinqr * cosqd / ( q3 * r2 * rd3 ) &
                        - 105.0_DP * sinqr * sinqd / ( q4 * r2 * rd3 ) &
                        + 90.0_dp * sinqr * sinqd / ( q3 * q3 * d3 * rd3 ) &
                        + 3.0_DP * r2**2 * sinqr * sinqd / ( q4 * d3 * rd3 ) &
                        + 270.0_DP * sinqr * sinqd / ( q3 * q3 * r2 * d2 * rd3 ) &
                        + 9.0_DP * r2 * sinqr * sinqd / ( q4 * d2 * rd3 ) &
                        - 27.0_DP * sinqr * sinqd / ( q4 * delta * rd3 )
        write(6,*) tabjqr2( i, k ), sphj0( q * r2 )
      enddo
    enddo
  else
    do k = 1, nshells * ids
      r2 = shells( k )
      do i = 1, nq
        q = dq * ( i - 0.5d0 )
        tabjqr2( i, k ) = sphj0( q * r2 )
      enddo
    enddo
  endif
  !
  v( :, : ) = 0
  do k = 1, nshells * ids
    r2 = shells( k )
    do i = 1, nq
       q = dq * ( i - 0.5d0 )
!       jqr2 = sphj0( q * r2 )
       qq = q / qf
       eps = 1 / levlou( qq, qf, lam )
       bq( i ) = ( 1 - 1 / eps ) / ( 4 * pi )
       r1 = dr / 2
       do j = 1, nr
          jqr1 = sphj0( q * r1 )
!          v( j, k ) = v( j, k ) + 8 * dq * bq( i ) * jqr1 * jqr2
          v( j, k ) = v( j, k ) + 8 * dq * bq( i ) * jqr1 * tabjqr2( i, k )
          r1 = r1 + dr
       end do
    end do
  end do
  vh( :, : ) = v( :, : )
  !
  do i = 1, nq 
     q = dq * ( i - 0.5d0 )
     r1 = dr / 2
     do j = 1, nr
        tab1( j, i ) = sphj0( q * r1 )
        r1 = r1 + dr
     end do
  end do

  do isite = 1, nsites
    s = 0.00001d0
    ds = 0.10d0

    if( isite .gt. 1 ) then
      v( :, : ) = vh( :, : )
    endif

    if( valenceGrid ) then
      write( avgName, '(A3,I6.6)' ) 'avg', isite
      inquire( file=avgName, exist=ex )
      if( .not. ex ) then
        write( avgName, '(A1,I6.6,A1,A3)' ) 'x', isite, delim_, 'avg' 
      endif
    else
      write( avgName, '(A3,A2,I4.4)' ) 'avg', siteSymbol( isite ), siteIndex( isite )
      inquire( file=avgName, exist=ex )
      if( .not. ex ) then
        write( avgName, '(A1,A2,I4.4,A1,A3)' ) 'z', siteSymbol(isite), siteIndex( isite ), delim_, 'avg'
      endif
    endif
    write( 6, * ) trim( avgName )

    open( unit=99, file=trim(avgName), form='formatted', status='old' )
    rewind 99

    do 
      read( 99, *, iostat=ierr ) dummy, dummy, ss, n
      if( ierr .lt. 0 ) then
        goto 100
      endif
      dn = n - nav
      if ( abs( s - ss ) .gt. 0.000001d0 ) then
        write(6,*) s, ss
        stop 'jive!'
      endif

      do k = 1, nshells * ids
        r2 = shells( k ) + delta
        if( solidShell ) then
          delta = shells( k ) / 2.0_DP
          r2 = shells( k ) 
        endif 
        th2 = 0
        if( s .gt. r2 ) then
          th2 = 1.0_DP
!!!! FIX HERE FOR SOLID
        elseif( solidShell ) then
          r2 = shells( k ) / 2.0_DP
          th2 = 3.0_DP * r2**2 + 12.0_DP * r2 * delta + 5.0_DP * delta**2
          th2 = 3.0_DP * ( r2 + delta ) * th2
          th2 = th2 - 2.0_DP * s * ( r2 + 3.0_DP * delta ) * ( 7.0_DP * r2 + 5.0_DP * delta )
          th2 = th2 + 5.0_DP * s**2 * ( r2 + 3.0_DP * delta )
          th2 = th2 * s**2 * ( s + delta - r2 )**2
          th2 = th2 / ( 16.0_DP * r2 * delta**3 * (r2 + delta )**3 )
        elseif( thickShell .and. ( s .gt. (shells(k) - delta ) ) ) then
          th2 = 3.0_DP * r2**2 + 12.0_DP * r2 * delta + 5.0_DP * delta**2 
          th2 = 3.0_DP * ( r2 + delta ) * th2
          th2 = th2 - 2.0_DP * s * ( r2 + 3.0_DP * delta ) * ( 7.0_DP * r2 + 5.0_DP * delta )
          th2 = th2 + 5.0_DP * s**2 * ( r2 + 3.0_DP * delta )
          th2 = th2 * s**2 * ( s + delta - r2 )**2
          th2 = th2 / ( 16.0_DP * r2 * delta**3 * (r2 + delta )**3 )
        endif
        if( th2 .lt. 0.0_DP .or. th2 .gt. 1.0_DP ) then
          write(6,*) 'bad th2', th2, r2, s
          stop
        endif
        dv = 0
        do i = 1, nq
          q = dq * ( i - 0.5d0 )
          j1 = sphj1( q * s )
          tj2 = j1 * th2
!          jqr2 = sphj0( q * r2 )
          r1 = dr / 2
          do j = 1, nr
            th1 = 0
            if ( s .gt. r1 ) th1 = 1
            tj1 = j1 * th1
            jqr1 = tab1( j, i )
            dv( j ) = dv( j ) + q * bq( i ) * ( tj1 * tabjqr2( i, k ) + tj2 * jqr1 )
!            dv( j ) = dv( j ) + q * bq( i ) * ( tj1 * jqr2 + tj2 * jqr1 )
            r1 = r1 + dr
          end do
        end do
        v( :, k ) = v( :, k ) + dv( : ) * ds * dn * 4 * dq / nav
      enddo

      s = s + ds
    end do
100 continue
    close( 99 )

    do k = 1, nshells
      if( valenceGrid ) then
        write( reoptName, '(A5,I6.6,A2,F4.2)' ) 'reopt', isite, &
                                                   '.R', shells( k )
      else
        if( NINT( shells(k) * 100 ) < 1000 ) then
          write( reoptName, '(A5,A2,I4.4,A2,F4.2)' ) 'reopt', siteSymbol( isite ), siteIndex( isite ), &
                                                     '.R', shells( k )
        else
          write( reoptName, '(A5,A2,I4.4,A2,F5.2)' ) 'reopt', siteSymbol( isite ), siteIndex( isite ), &
                                                     '.R', shells( k )
        endif
      endif

      open( unit=99, file=reoptName, form='formatted', status='unknown' )
      rewind( 99 )
      r1 = dr / 2
      do j = 1, nr
        if( doubleShell ) then
          write ( 99, '(7(1x,1e15.8))' ) r1, 5.0_dp * vh( j, k ) - 4.0_dp*vh( j, k + nshells ), &
                                             5.0_dp*v( j, k )- 4.0_dp*v( j, k + nshells ), &
              vh( j, k ), vh( j, k + nshells ), v( j, k ), v( j, k + nshells )
        else
          write ( 99, '(3(1x,1e15.8))' ) r1, vh( j, k ), v( j, k )
        endif
        r1 = r1 + dr
      end do
      close( unit=99 )
    enddo

  enddo

#if 0
!  open( unit=99, file='avden', form='formatted', status='unknown' )
!  rewind 99
  do while ( s .lt. 40 )
     read ( 99, * ) dummy, dummy, ss, n
     dn = n - nav
     if ( abs( s - ss ) .gt. 0.000001d0 ) stop 'jive!'
     th2 = 0
     if ( s .gt. r2 ) th2 = 1
     dv = 0
     do i = 1, nq
        q = dq * ( i - 0.5d0 )
        j1 = sphj1( q * s )
        tj2 = j1 * th2
        jqr2 = sphj0( q * r2 )
        r1 = dr / 2
        do j = 1, nr
           th1 = 0
           if ( s .gt. r1 ) th1 = 1
           tj1 = j1 * th1
           jqr1 = tab1( j, i )
           dv( j ) = dv( j ) + q * bq( i ) * ( tj1 * jqr2 + tj2 * jqr1 )
           r1 = r1 + dr
        end do
     end do
     v( : ) = v( : ) + dv( : ) * ds * dn * 4 * dq / nav
     s = s + ds
  end do
  close( unit=99 )
  !
  open( unit=99, file='reopt', form='formatted', status='unknown' )
  rewind 99
  r1 = dr / 2
  do j = 1, nr
     write ( 99, '(3(1x,1e15.8))' ) r1, vh( j ), v( j )
     r1 = r1 + dr
  end do
  close( unit=99 )
#endif
  !
  deallocate( v, vh, tab1, shells, siteSymbol, siteIndex )
  !
end program vhommod

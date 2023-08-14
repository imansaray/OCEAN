! Copyright (C) 2017 - 2019 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
module screen_centralPotential
  use ai_kinds, only : DP

  implicit none
  private
  save

  type potential
    real(DP), allocatable :: pot(:)
    real(DP), allocatable :: rad(:)
    integer :: z, n, l
!    real(DP) :: deltaR  ! this is the spacing on the written out induced potential
  end type potential



  type( potential ), allocatable, target :: znlPotentials( : )

  public :: potential
  public :: screen_centralPotential_newScreenShell
  public :: screen_centralPotential_prepAll, screen_centralPotential_freePot
  public :: screen_centralPotential_loadInternal, screen_centralPotential_freeInternal
  public :: screen_centralPotential_findNextByZ, screen_centralPotential_countByZ
  public :: screen_centralPotential_dumplog
!  public :: screen_centralPotential_NOutGrid, screen_centralPotential_MakeOutGrid

  contains

#if 0
  pure function screen_centralPotential_NOutGrid( Pot, rmax ) result( n )
    type( potential ), intent( in ) :: Pot
    real( DP ), intent( in ) :: rmax
    integer :: n

    n = max( floor( rmax / Pot%deltaR ) + 1, 1 )
  end function screen_centralPotential_NOutGrid

  subroutine screen_centralPotential_MakeOutGrid( Pot, rmax, rlist, ierr )
    type( potential ), intent( in ) :: Pot
    real( DP ), intent( in ) :: rmax
    real( DP ), intent( out ) :: rlist( : )
    integer, intent( inout ) :: ierr

    integer :: i, n  

    n = screen_centralPotential_NOutGrid( Pot, rmax )

    if( size( rlist, 1 ) .lt. n ) then
      ierr = 1
      return
    endif

    do i = 1, n
      rlist( i ) = real( i - 1, DP ) * Pot%deltaR
    enddo

  end subroutine screen_centralPotential_MakeOutGrid
#endif

  pure function screen_centralPotential_countByZ( Z ) result( n )
    integer, intent( in ) :: Z
    integer :: n
    !
    integer :: PotSize, i
    !
    n = 0
    PotSize = size( znlPotentials )
    do i = 1, PotSize
      if( znlPotentials( i )%z .eq. Z ) n = n + 1
    enddo
    
  end function screen_centralPotential_countByZ

  subroutine screen_centralPotential_findNextByZ( Z, i, PotPointers, ierr )
    integer, intent( in ) :: Z
    integer, intent( inout ) :: i
    type( potential ), pointer, intent( out ) :: PotPointers
    integer, intent( inout ) :: ierr
    !
    integer :: PotSize, j

!    if( size( PotPointers ) .lt. screen_centralPotential_countByZ( Z ) ) then
!      write(6,*) 'PotPointers too small!'
!      ierr = 1
!      return
!    endif

    PotSize = size( znlPotentials )
    
    do j = i+1, PotSize
      if( znlPotentials( j )%z .eq. Z ) then
        PotPointers => znlPotentials( j )
        return
      endif
    enddo

    ierr = 1

  end subroutine screen_centralPotential_findNextByZ
    

  subroutine screen_centralPotential_loadInternal( ierr )
    use ocean_mpi, only : myid, root
    integer, intent( inout ) :: ierr
    !
    integer, allocatable :: tmp_znl( :, : )
    integer :: znlLength, i
    !
    ! load up potentials to be screened later
    ! 2000 should be more than enough for the entire periodic table at the moment
    if( myid .eq. root ) then
      allocate( tmp_znl( 3, 2000 ), stat=ierr )
    else
      allocate( tmp_znl( 3, 1 ), stat=ierr )
    endif
    call screen_centralPotential_prepAll( tmp_znl, znlLength, ierr )
    if( ierr .ne. 0 ) return
    if( myid .eq. root ) then
      do i = 1, znlLength
        write(6,*) tmp_znl( :, i )
      enddo
    endif
!    write(6,*) "znlPotentials", znlLength
    allocate( znlPotentials( znlLength ), stat=ierr )
    if( ierr .ne. 0 ) return
!    call screen_centralPotential_loadAll( znlLength, znlPotentials, tmp_znl, ierr )
!    if( ierr .ne. 0 ) return
    do i = 1, znlLength
      call screen_centralPotential_loadSingle( znlPotentials( i ), tmp_znl( :, i ), ierr )
      if( ierr .ne. 0 ) return
    enddo
    deallocate( tmp_znl )

  end subroutine screen_centralPotential_loadInternal

  subroutine screen_centralPotential_freeInternal()
    integer :: i

    do i = 1, size( znlPotentials ) 
      call screen_centralPotential_freePot( znlPotentials( i ) )
    enddo

    if( allocated( znlPotentials ) ) deallocate( znlPotentials )
  end subroutine screen_centralPotential_freeInternal

  subroutine screen_centralPotential_freePot( pot )
    type( potential ), intent( inout ) :: pot


    if( allocated( pot%pot ) ) deallocate( pot%pot )
    if( allocated( pot%rad ) ) deallocate( pot%rad )

  end subroutine screen_centralPotential_freePot

  subroutine screen_centralPotential_loadSingle( Pot, znl, ierr )
    use ocean_mpi, only : myid, root, comm, MPI_INTEGER
    type( potential ), intent( inout ) :: Pot
    integer, intent( in ) :: znl(:)
    integer, intent( inout ) :: ierr

    integer :: ierr_

    ierr_ = 0

    if( myid .eq. root ) then
      call screen_centralPotential_load( znl(1), znl(2), znl(3), Pot, ierr )
      write(6,*) znl(:)
      write(6,*) Pot%z, Pot%n, Pot%l
    endif
#ifdef MPI
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
#endif
    if( ierr .ne. 0 ) return
    if( ierr_ .ne. 0 ) then
      ierr = ierr_
      return
    endif
    call screen_centralPotential_share( Pot, myid, root, comm, ierr )
    if( ierr .ne. 0 ) return

!      write(3000+myid, * ) pot%z, pot%n, pot%l
!      flush(3000+myid)


  end subroutine screen_centralPotential_loadSingle


  subroutine screen_centralPotential_loadAll( Nznl, znlPot, znl, ierr )
    use ocean_mpi, only : myid, root, comm, MPI_INTEGER
    integer, intent( in ) :: Nznl
    type( potential ), intent( inout ) :: znlPot( : )
    integer, intent( in ) :: znl(:,:)
    integer, intent( inout ) :: ierr

    integer :: i, ierr_

!    Nznl = size( znlPot, 1 )
    ierr_ = 0
    write(6,*) myid, Nznl

    do i = 1, Nznl
      if( myid .eq. root ) then
        call screen_centralPotential_load( znl(1,i), znl(2,i), znl(3,i), znlPot(i), ierr )
      endif
#ifdef MPI
      call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
#endif
      if( ierr .ne. 0 ) return
      if( ierr_ .ne. 0 ) then
        ierr = ierr_
        return
      endif
      call screen_centralPotential_share( znlPot(i), myid, root, comm, ierr )
      if( ierr .ne. 0 ) return
    enddo
    
!    call screen_centralPotential_dumpLog( znlPot, myid )
#ifdef PRINTLOG
    write(1000+myid, * ) 'Central Potential: ', Nznl

    do i = 1, Nznl
      write(1000+myid, * ) znlpot( i )%z, znlpot( i )%n, znlpot( i )%l
    enddo
#endif

  end subroutine screen_centralPotential_loadAll

  subroutine screen_centralPotential_dumpLog( pot, myid )
    type( potential ), intent( in ) :: pot
    integer, intent( in ) :: myid
    integer :: i

#ifdef PRINTLOG
    write(1000+myid, * ) 'Central Potential: ', size( pot%pot, 1 )
    write(1000+myid, * ) pot%z, pot%n, pot%l
    do i = 1, size( pot%pot, 1 )
      write(1000+myid, * ) pot%rad(i), pot%pot(i)
    enddo
#endif
  end subroutine screen_centralPotential_dumpLog


  subroutine screen_centralPotential_share( pot, myid_, root_, comm_, ierr )
    use ocean_mpi, only : MPI_INTEGER, MPI_DOUBLE_PRECISION
    type( potential ), intent( inout ) :: pot
    integer, intent( in ) :: myid_, root_, comm_
    integer, intent( inout ) :: ierr
    !
    integer :: intArray( 4 )

    if( myid_ .eq. root_ ) then
      intArray(4) = size( pot%pot, 1 )
      intArray(1) = pot%z
      intArray(2) = pot%n
      intArray(3) = pot%l
    endif

    call MPI_BCAST( intArray, 4, MPI_INTEGER, root_, comm_, ierr )
    if( ierr .ne. 0 ) return
    if( myid_ .ne. root_ ) then
!      write(6,*) 'pot_share'
!      write(6,*) intArray(:)
      allocate( pot%pot( intArray(4) ), pot%rad( intArray(4) ) )
      pot%z = intArray(1)
      pot%n = intArray(2)
      pot%l = intArray(3)
    endif
    call MPI_BCAST( pot%pot, intArray(4), MPI_DOUBLE_PRECISION, root_, comm_, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( pot%rad, intArray(4), MPI_DOUBLE_PRECISION, root_, comm_, ierr )
    if( ierr .ne. 0 ) return

  end subroutine screen_centralPotential_share

      

  subroutine screen_centralPotential_prepAll( znl, n, ierr )
    use ocean_mpi, only : myid, root, comm, MPI_INTEGER
    use screen_system, only : screen_system_mode
    integer, intent( out ) :: znl(:,:)
    integer, intent( out ) :: n
    integer, intent( inout ) :: ierr

    integer :: i, maxSize, j

    n = 0
    if( screen_system_mode() .eq. 'grid' ) then
      n = 1
      znl(1,1) = 0
      znl(2,1) = 0
      znl(3,1) = 0
      return
    endif

    if( myid .ne. root ) then
#ifdef MPI
      call MPI_BCAST( n, 1, MPI_INTEGER, root, comm, ierr )
#endif
      return
    endif

    maxSize = size( znl, 2 )
    open( unit=99, file='edgelist', form='formatted', status='old' )

    do i = 1, maxSize
      read(99,*,iostat=j) znl(1:3, i )
      if( j .lt. 0 ) then
        n = i - 1
        close( 99 )
#ifdef MPI
        call MPI_BCAST( n, 1, MPI_INTEGER, root, comm, ierr )
#endif
        return
      elseif( j .gt. 0 ) then
        ierr = j
!        return
      endif
    enddo
    close( 99 )

    ! if this isn't here the other mpi will hang above ...
#ifdef MPI
        call MPI_BCAST( n, 1, MPI_INTEGER, root, comm, ierr )
#endif

    write(6,'(A,I8)') 'Too many unique edges! Greater than: ', maxSize
    ierr = 1
  end subroutine screen_centralPotential_prepAll

! Need to add options to do other things with the screening
! 1. plain shell
! 2. shell with width
! 3. double shell method to give 0 potential inside R/2
  subroutine screen_centralPotential_newScreenShell( pot, newPot, rad, ierr, width, double_shell )
    type( potential ), intent( in ) :: pot
    type( potential ), intent( out ) :: newPot
    real(DP), intent( in ) :: rad
    integer, intent( inout ) :: ierr
    real(DP), intent( in ), optional :: width
    logical, intent( in ), optional :: double_shell
    logical, intent( in ), optional :: sawtooth

    integer :: i
    real(DP) :: invRad, ww, denom, f3, f4, f5, x
    logical :: do_double_shell, do_sawtooth
!#ifdef DEBUG
    character(len=11) :: filename
!#endif
  
    if( present( width ) ) then
      ww = width
      if( ww .lt. 0.0 ) ww = 0.0_DP
    else
      ww = 0.0_DP
    endif

    if( present( double_shell ) ) then
      do_double_shell = double_shell
    else
      do_double_shell = .false.
    endif
    if( do_double_shell ) then
      ww = 0.0_DP
    endif

    if( present( sawtooth ) ) then
      do_sawtooth = sawtooth
    else 
      do_sawtooth = .false.
    endif
    if( do_sawtooth ) then
      do_double_shell = .false.
      ww = 0.0_DP
    endif

    if( ( .not. allocated( pot%pot ) ) .or. ( .not. allocated( pot%rad ) ) ) then
      ierr = 5
      return
    endif

    allocate( newPot%pot( size( pot%pot ) ), newPot%rad( size( pot%rad ) ) )
    newPot%rad(:) = pot%rad(:)


    invRad = 1.0_DP / rad

    if( ww .gt. 0.0_DP ) then
      denom = 4.0_DP * rad * ww**2 
      denom = denom * ( rad + ww )**3
      denom = 1.0_DP / denom
      f3 = ( rad**2 + 4.0_DP * rad*ww + 5.0_DP*ww**2 ) * denom

      denom = denom / ( 4.0_DP * ww )
      f4 = ( rad + 3.0_DP * ww ) * ( rad + 5.0_DP * ww ) * denom
      f5 = ( rad + 3.0_DP * ww ) * denom
    else
      f3 = 0.0_DP
      f4 = 0.0_DP
      f5 = 0.0_DP
    endif
    
    if( do_double_shell ) then
      do i = 1, size( pot%pot )
        if( newPot%rad( i ) .lt. ( 4.0_dp * rad / 5.0_dp ) ) then
          newPot%pot( i ) = pot%pot( i )
        elseif( newPot%rad( i ) .le. rad ) then
          newPot%pot( i ) = pot%pot( i ) + 5.0_dp * invrad - 4.0_dp / newPot%rad( i )
        else
          newPot%pot( i : size( newPot%pot ) ) = 0.0_DP
          exit
        endif
      enddo
    elseif( do_sawtooth ) then
      if( rad .le. 2.0_DP ) then
        alpha = 0.0_DP
        beta = rad / 3.0_DP
        gamm = 2.0_DP * rad / 3.0_DP
        delta = rad
      else
        alpha = rad - 2.0_DP
        beta = alpha + 2.0_DP/3.0_DP
        gamm = alpha + 4.0_DP/3.0_DP
        delta = rad
      endif
          aaa = -(1.0_DP/6.0_DP)*(-beta**3 + 3*beta*gamm**2 - (-3 + gamm)*gamm**2 - 3*(1 + gamm)*delta**2 + &
          2*delta**3)/ &
        (((-alpha**4 + 4*alpha*beta**3 - 2*(-2 + beta)*beta**3 - 4*(1 + beta)*gamm**3 + 3*gamm**4)* &
             (beta**3 - 3*beta*gamm**2 + (-3 + gamm)*gamm**2 + 3*(1 + gamm)*delta**2 - 2*delta**3))/72. - &
          ((-alpha**3 + 3*alpha*beta**2 - (-3 + beta)*beta**2 - 3*(1 + beta)*gamm**2 + 2*gamm**3)* &
             (beta**4 - 4*beta*gamm**3 + 2*(-2 + gamm)*gamm**3 + 4*(1 + gamm)*delta**3 - 3*delta**4))/72.)
    bbb = (-12.0_DP*(alpha**3 - 3.0_DP*beta**2 - 3.0_DP*alpha*beta**2 + beta**3 &
                            + 3.0_DP*gamm**2 + 3.0_DP*beta*gamm**2 - 2.0_DP*gamm**3))/ &
               (alpha**4*beta**3 - alpha**3*beta**4 - beta**6 - alpha*beta**6 + beta**7 &
                - 3.0_DP*alpha**4*gamm**2 - 3.0_DP*alpha**4*beta*gamm**2  &
                + 12.0_DP*beta**3*gamm**2 + 12.0_DP*alpha*beta**3*gamm**2 + 3.0_DP*beta**4*gamm**2 &
                + 12.0_DP*alpha*beta**4*gamm**2 - 9.0_DP*beta**5*gamm**2 + 4.0_DP*alpha**3*gamm**3 &
                + alpha**4*gamm**3 + 4.0_DP*alpha**3*beta*gamm**3 - 12.0_DP*beta**2*gamm**3 &
                - 12.0_DP*alpha*beta**2*gamm**3 - 8.0_DP*beta**3*gamm**3 &
                - 16.0_DP*alpha*beta**3*gamm**3 + 12.0_DP*beta**4*gamm**3 - 2.0_DP*alpha**3*gamm**4 &
                + 6.0_DP*beta**2*gamm**4 + 6.0_DP*alpha*beta**2*gamm**4 - 5.0_DP*beta**3*gamm**4 &
                - gamm**6 - beta*gamm**6 + gamm**7 + 3.0_DP*alpha**4*delta**2 &
                - 12.0_DP*beta**3*delta**2 - 12.0_DP*alpha*beta**3*delta**2 &
                + 6.0_DP*beta**4*delta**2 + 3.0_DP*alpha**4*gamm*delta**2 &
                - 12.0_DP*beta**3*gamm*delta**2 - 12.0_DP*alpha*beta**3*gamm*delta**2 &
                + 6.0_DP*beta**4*gamm*delta**2 + 12.0_DP*gamm**3*delta**2 &
                + 12.0_DP*beta*gamm**3*delta**2 + 3.0_DP*gamm**4*delta**2 &
                + 12.0_DP*beta*gamm**4*delta**2 - 9.0_DP*gamm**5*delta**2 - 4.0_DP*alpha**3*delta**3 &
                - 2.0_DP*alpha**4*delta**3 + 12.0_DP*beta**2*delta**3 &
                + 12.0_DP*alpha*beta**2*delta**3 + 4.0_DP*beta**3*delta**3 &
                + 8.0_DP*alpha*beta**3*delta**3 - 4.0_DP*beta**4*delta**3 &
                - 4.0_DP*alpha**3*gamm*delta**3 + 12.0_DP*beta**2*gamm*delta**3 &
                + 12.0_DP*alpha*beta**2*gamm*delta**3 - 4.0_DP*beta**3*gamm*delta**3 &
                - 12.0_DP*gamm**2*delta**3 - 12.0_DP*beta*gamm**2*delta**3 &
                - 12_DP*gamm**3*delta**3 - 20.0_DP*beta*gamm**3*delta**3 + 14.0_DP*gamm**4*delta**3 &
                + 3.0_DP*alpha**3*delta**4 - 9.0_DP*beta**2*delta**4 - 9.0_DP*alpha*beta**2*delta**4 &
                + 3.0_DP*beta**3*delta**4 + 9.0_DP*gamm**2*delta**4 + 9.0_DP*beta*gamm**2*delta**4 &
                - 6.0_DP*gamm**3*delta**4)

      do i = 1, size( pot%pot )
        if( newPot%rad( i ) .le. alpha ) then
          newPot%pot( i ) = pot%pot( i )
        elseif( newPot%rad( i ) .le. beta ) then
          newPot%pot( i ) = pot%pot( i ) &
                          - (1.0_D_*/(12.0_DP*newPot%rad( i ))) * &
                            ( aaa*( newPot%rad( i )**4 - 2.0_DP*newPot%rad( i )**3*alpha &
                                  - 2.0_DP*newPot%rad( i )*beta**2*(-3.0_DP-3.0_dp*alpha*beta) &
                                  + 6.0_DP*newPot%rad( i )*(1.0_DP+beta)*gamm**2 &
                                  + 4.0_DP * newPot%rad( i ) * gamm**3 )
                            + 2.0_DP*bbb*newPot%rad( i )* ( beta**3 - 3.0_DP*beta*gamm**2 &
                                  +(gamm-3.0_DP)*gamm**2 + 3.0_DP*(1.0_DP*gamm)*delta**2  &
                                  - 2.0_dp * delta**3 ) )
        elseif( newPot%rad( i ) .le. gamm ) then
          newPot%pot( i ) = pot%pot( i ) &
                    - (-(aaa*(x**4 + alpha**4 - 2.0_DP*newPot%rad( i )**3*(1.0_DP + beta) &
                        + 2.0_DP*beta**3*(-2.0_dp - 2.0_DP*alpha + beta) +  &
                        2.0_DP*newPot%rad( i )*(3.0_DP + 3.0_DP*beta - 2.0_DP*gamm)*gamm**2)) + &
                   bbb*(-newPot%rad( i )**4 + 2.0*newPot%rad( i )**3*beta + beta**4 +  &
                      2.0_DP*newPot%rad( i )*(gamm**2*(-3.0_DP - 3.0_DP*beta + gamm) + &
                      3.0_DP*(1.0_DP + gamm)*delta**2 - 2.0D_*delta**3)))/(12.0_DP*newPot%rad( i ))
        elseif( newPot%rad( i ) .lt. delta ) then
          newPot%pot( i ) = pot%pot( i ) & 
            -(-(aaa*(alpha**4 - 4.0_DP*alpha*beta**3 + 2.0_DP*(-2.0_DP + beta)*beta**3 &
                    + 4.0_DP*(1.0_DP + beta)*gamm**3 - 3.0_DP*gamm**4)) + &
                     bbb*(newPot%rad( i )**4 + beta**4 - 2.0_DP*newPot%rad( i )**3*(1.0_DP + gamm) &
                          + 2.0_DP*gamm**3*(-2.0_DP - 2.0_DP*beta + gamm) + &
                        2.0_DP*newPot%rad( i )*(3.0_DP + 3.0_DP*gamm - 2.0_DP*delta)*delta**2)) &
                /(12.0_DP*newPot%rad( i ))
                                    
        else                        
          newPot%pot( i : size( newPot%pot ) ) = 0.0_DP 
          exit
        endif
      enddo        

    else
    do i = 1, size( pot%pot ) 
      if( newPot%rad( i ) .lt. rad - ww) then
        newPot%pot( i ) = pot%pot( i ) + invRad
      elseif( newPot%rad( i ) .le. rad + ww ) then
        x = newPot%rad( i ) - ( rad - ww )
        newPot%pot( i ) = x*f5
        newPot%pot( i ) = x*( f4 - newPot%pot( i ) )
        newPot%pot( i ) = x**3 * ( f3 - newPot%pot( i ) )
        newPot%pot( i ) = pot%pot( i ) + invRad - newPot%pot( i )
        ! 
      else
        newPot%pot( i : size( newPot%pot ) ) = 0.0_DP
        exit
      endif
    enddo
    endif

#ifdef DEBUG
    write(6,*) 'newScreenShell', rad, newPot%rad( size( newPot%rad ) )
    if( rad .le. 9.99) then
      write(filename,'(A,F4.2)') 'vpert.', rad
    else
      write(filename,'(A,F5.2)') 'vpert.', rad
    endif
    open(unit=99,file=filename, form='formatted' )
    do i = 1, size( pot%pot ) 
      write(99,* ) newPot%rad( i ), newPot%pot( i ), pot%pot( i )
    enddo
    close( 99 )
#endif

    newPot%z = pot%z
    newPot%n = pot%n
    newPot%l = pot%l
!    newPot%deltaR = pot%deltaR

  end subroutine screen_centralPotential_newScreenShell

  subroutine screen_centralPotential_load( z, n, l, pot, ierr )
    integer, intent( in ) :: z, n, l
    type( potential ), intent( out ) :: pot
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: tmpPot(:), tmpRad(:), tmp(:)
    integer :: maxLength, curLength, fh

    character(len=26) :: fileName

    pot%z = z
    pot%n = n
    pot%l = l
!    pot%deltaR = 0.02_DP
    
    if( z .lt. 0 ) then
      call screen_centralPotential_makeBare( z, pot, ierr )
      return
    elseif( z .eq. 0 ) then
      call screen_centralPotential_makeValence( z, pot, ierr )
      return
#ifdef FAKEVAL
    else
      write(6,*) 'Using fake-valence!' 
      call screen_centralPotential_makeValence( z, pot, ierr )
      return
#endif 
    endif

    write(fileName,'(A17,I3.3,A1,I2.2,A1,I2.2)') 'zpawinfo/vc_barez', z, 'n', n, 'l', l
    fh = 99
    open( unit=fh, file=fileName, form='formatted', status='old', iostat=ierr )
    if( ierr .ne. 0 ) then
      write(6,*) 'Failed to open ', fileName, ' with error ', ierr
      return
    endif
    rewind( 99 )

    maxLength = 4000
    curLength = 1

    allocate( tmpPot( maxLength ), tmpRad( maxLength ), STAT=ierr )
    if( ierr .ne. 0 ) return

    do while( maxLength .le. 512000 )

      call doRead( fh, curLength, MaxLength, tmpRad, tmpPot, ierr )

      if( curLength .lt. 0 ) goto 11

      allocate( tmp( maxLength ), stat=ierr )
      if( ierr .ne. 0 ) return
      tmp(:) = tmpPot(:)
      deallocate( tmpPot )
      allocate( tmpPot( maxLength*2 ), stat=ierr )
      if( ierr .ne. 0 ) return
      tmpPot(1:maxLength) = tmp(:)
      tmp(:) = tmpRad(:)
      deallocate( tmpRad )
      allocate( tmpRad( maxLength*2 ), stat=ierr )
      tmpRad(1:maxLength) = tmp(:)

      curLength = maxLength + 1
      maxLength = 2*maxLength
      deallocate( tmp )
    enddo
    !    If we exit the loop then we've given up not reached the end of the file
    ierr = -1
    return

11  continue
    curLength = abs( curLength )
    allocate( pot%rad( curLength ), pot%pot( curLength ), stat=ierr )
    if( ierr .ne. 0 ) return
    pot%rad(:) = tmpRad(1:curLength)
    pot%pot(:) = tmpPot(1:curLength)
    deallocate( tmpPot, tmpRad )

    close( fh )

    write(6,* ) z, n, l
    write(6,*) curLength
    write(6,*) pot%rad(1), pot%pot(1)
    write(6,*) pot%rad(curLength), pot%pot(curLength)

  end subroutine screen_centralPotential_load

  subroutine doRead( fh, curLength, maxLength, rad, pot, ierr )
    integer, intent( in ) :: fh, maxLength
    integer, intent( inout ) :: curLength
    real(DP), intent( inout ) :: rad(:), pot(:)
    integer, intent( inout ) :: ierr
    !
    integer :: i, j, start

    start = curLength

    do i = start, maxLength
      read(fh,*,iostat=j) rad(i), pot(i)
      if( j .lt. 0 ) then
        curLength = -( i - 1 )
        return
      elseif( j .gt. 0 ) then
        ierr = j
        return
      endif
      ! if j .eq. 0 then keep going!
    enddo

  end subroutine doRead

  subroutine screen_centralPotential_makeBare( z, pot, ierr )
    integer, intent( in ) :: z
    type( potential ), intent( out ) :: pot
    integer, intent( inout ) :: ierr

    real(DP) :: rat, xrat, xr1, dl
    integer :: i
    real(DP), parameter :: rmin = 0.00000001d0
    real(DP), parameter :: rmax = 800.d0
    integer, parameter :: nr = 4096

    
    allocate( pot%rad( nr ), pot%pot( nr ), stat=ierr )
    if( ierr .ne. 0 ) return

    rat=rmax/rmin
    dl=dlog(rat)/dble(nr)
    xrat=dexp(dl)
    xr1=dsqrt(xrat)-dsqrt(1.d0/xrat)
    do i=1,nr
      pot%rad(i)=rmin*xrat**dble(i)
      pot%pot(i)=-1.0_dp / pot%rad(i)
    end do

  end subroutine screen_centralPotential_makeBare

  subroutine screen_centralPotential_makeValence( z, pot, ierr )
    use screen_system, only : screen_system_Volume, screen_system_xmesh
    use ocean_constants, only : PI_DP
    integer, intent( in ) :: z
    type( potential ), intent( out ) :: pot
    integer, intent( inout ) :: ierr

    real(DP) :: rat, xrat, dl, valRadius, invR, invR2
    integer :: i, xmesh(3), temp
    real(DP), parameter :: rmin = 0.00000001d0
    real(DP), parameter :: rmax = 800.d0
    integer, parameter :: nr = 4096


    allocate( pot%rad( nr ), pot%pot( nr ), stat=ierr )
    if( ierr .ne. 0 ) return

    xmesh(:) = screen_system_xmesh()
    i = product( xmesh(:) )
    if( i .lt. 1 ) then
      ierr = -5
      return
    endif
    ! V = 4 pi / 3 R^3
    ! R = ( 3 V / 4 pi ) ^(1/3)
    valRadius = ( 0.75_DP * screen_system_Volume() ) / ( PI_DP * real( i, DP ) )
    valRadius = valRadius ** (1.0_DP/3.0_DP)

    write(6,*) 'val radius: ', valRadius
    

    rat=rmax/rmin
    dl=dlog(rat)/dble(nr)
    xrat=dexp(dl)
!    xr1=dsqrt(xrat)-dsqrt(1.d0/xrat)

    invR = 0.5_DP / valRadius
    invR2 = 1.0_DP / ( valRadius**2)
    temp = floor( ( log( valRadius ) - log( rmin ) ) / dl )
    temp = min( temp, nr )
    ! V = - 1 / (2 R ) *( 3 - r^2/R^2)
    do i=1,temp
      pot%rad(i)=rmin*xrat**dble(i)
      pot%pot(i)= -invR * ( 3.0_DP - invR2 * pot%rad(i)**2 )
    end do
    do i=temp,nr
      pot%rad(i)=rmin*xrat**dble(i)
      pot%pot(i)=-1.0_dp / pot%rad(i)
    end do

  end subroutine screen_centralPotential_makeValence

end module screen_centralPotential

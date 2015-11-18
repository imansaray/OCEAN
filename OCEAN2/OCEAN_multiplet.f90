! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
module OCEAN_multiplet
  use AI_kinds
!  use OCEAN_constants, only : eVtoHartree

  private
  save

  real( DP ), allocatable :: mpcr( :, :, :, :, :, : )
  real( DP ), allocatable :: mpci( :, :, :, :, :, : )
  real( DP ), allocatable :: mpm( :, :, : )
  real( DP ), allocatable :: mhr( : )
  real( DP ), allocatable :: mhi( : )
  real( DP ), allocatable :: cms( : ), cml( : ), vms( : )
  real( DP ), allocatable :: hcml( : ), hvml( : ), hcms( : ), hvms( : )
  real( DP ), allocatable :: somelr( :, : )
  real( DP ), allocatable :: someli( :, : )

#ifdef __INTEL_COMPILER
!DIR$ attributes align: 64 :: mpcr, mpci, mpm, mhr, mhi, somelr, someli
#endif


  integer, allocatable :: nproj( : ), hvnu(:)
  integer, allocatable :: ibeg( : )
  integer, allocatable :: jbeg( : )
  integer, allocatable :: mham( : )

  real(DP ) :: xi

  integer :: npmax
  integer :: nptot
  integer :: lmin
  integer :: lmax
  integer :: lvl
  integer :: lvh
  integer :: jtot
  integer :: itot

  logical :: is_init = .false.
  logical :: do_staggered_sum 



  integer :: push_psi_sum
  integer :: pull_psi_sum
  integer :: share_psi_sum

  integer :: ct_n
  integer :: fg_n
  integer :: so_n

  integer, allocatable :: ct_list( :, : )
  integer, allocatable :: fg_list( : )
  integer, allocatable :: so_list( :, : )

  public OCEAN_create_central, OCEAN_mult_act, OCEAN_mult_single, OCEAN_mult_slice

  contains

  subroutine OCEAN_mult_create_par( sys, ierr )
    use OCEAN_mpi
    use OCEAN_system
    implicit none

    type(O_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    integer(8), allocatable :: proc_load( : )
    integer, allocatable :: proc_list(:)
    integer, allocatable :: ct_list_temp( :, : ), fg_list_temp( : ), so_list_temp( :, : )


    integer( 8 ) :: total, temp_load
    integer :: ia, ja, il, im, dumi, ip, jp, least_proc

    ct_n = 0
    fg_n = 0
    so_n = 0


    total = 0
    allocate( proc_load( -1:nproc-1 ) )
    proc_load( : ) = 0

    dumi = 0
    do il = lvl, lvh
      dumi = dumi + ( 2 * il + 1 )
    enddo

!    allocate( ct_list_temp( 3, sys%nalpha * ( 2 * ( lvh - lvl + 1 ) + 1 ) ), &
    allocate( ct_list_temp( 3, sys%nalpha * dumi ), &
              fg_list_temp( lvh - lvl + 1 ), so_list_temp( 2, sys%nalpha**2 ) )
!fg
    do il = lvl, lvh

      least_proc = 0
      do ip = 1, nproc-1
        if( proc_load( ip ) < proc_load( least_proc ) ) least_proc = ip
      enddo
      
      dumi = 2 * sys%nalpha * ( 2 * il + 1 ) * nproj( il ) * sys%nkpts * sys%num_bands &
           + ( 4 * ( 2 * sys%cur_run%ZNL(3) + 1 ) * ( 2 * il + 1 ) * nproj( il ) )**2
      proc_load( least_proc ) = proc_load( least_proc ) + dumi
      total = total + dumi

      if( least_proc .eq. myid ) then
        fg_n = fg_n + 1
        fg_list_temp( fg_n ) = il
      endif

    enddo
! ct

    do ia = 1, sys%nalpha
      do il = lvl, lvh
        do im = -il, il
          
          least_proc = 0
          do ip = 1, nproc-1
            if( proc_load( ip ) < proc_load( least_proc ) ) least_proc = ip
          enddo
          
          dumi = 2 * nproj( il ) * sys%nkpts * sys%num_bands + nproj( il ) + nproj( il )**2
          
          proc_load( least_proc ) = proc_load( least_proc ) + dumi
          total = total + dumi

          if( least_proc .eq. myid ) then
            ct_n = ct_n + 1
            ct_list_temp( :, ct_n ) = (/ ia, il, im /)
          endif

        enddo
      enddo
    enddo


! s-o
    if( sys%cur_run%ZNL(3) > 0 ) then
      do ia = 1, sys%nalpha
        do ja = 1, sys%nalpha

          least_proc = 0
          do ip = 1, nproc-1
            if( proc_load( ip ) < proc_load( least_proc ) ) least_proc = ip
          enddo

          dumi = sys%nkpts * sys%num_bands
          proc_load( least_proc ) = proc_load( least_proc ) + dumi
          total = total + dumi

          if( least_proc .eq. myid ) then
            so_n = so_n + 1
            so_list_temp( :, so_n ) = (/ ia, ja /)
          endif

        enddo
      enddo
    endif

    allocate( proc_list( 0:nproc-1 ) )
    do ip = 0, nproc - 1
      proc_list(ip) = ip
    enddo

    do ip = 1, nproc-1
      jp = ip
      do while( ( jp > 0 ) .and. ( proc_load( jp-1 ) > proc_load( jp ) ) )
        temp_load = proc_load( jp )
        proc_load( jp ) = proc_load( jp-1 )
        proc_load( jp-1 ) = temp_load

        dumi = proc_list( jp )
        proc_list( jp ) = proc_list( jp - 1 )
        proc_list( jp-1 ) = dumi

        jp = jp - 1
      enddo
    enddo
    
    if( myid .eq. 0 ) then
      do ip = 0, nproc-1
        write(6,*) proc_list( ip ), proc_load( ip ), dble(proc_load( ip ))/dble(total)
      enddo
    endif

    dumi = -1
    do ip = 0, nproc-1
      if( proc_list( ip ) .eq. myid ) then
        dumi = ip
      endif
    enddo

    if( dumi .gt. 0 ) then
      pull_psi_sum = proc_list( dumi - 1 )
    else
      pull_psi_sum = -1
    endif

    if( dumi .lt. nproc - 1 ) then
      push_psi_sum = proc_list( dumi + 1 )
    else
      push_psi_sum = -1
    endif

    share_psi_sum = proc_list( nproc-1 )


    if( ct_n .gt. 0 ) then
      allocate( ct_list( 3, ct_n ) )
      ct_list( :, 1:ct_n ) = ct_list_temp( :, 1:ct_n )
    else
      allocate( ct_list(1,1) )
    endif

    if( fg_n .gt. 0 ) then
      allocate( fg_list( fg_n ) )
      fg_list( 1:fg_n ) = fg_list_temp( 1:fg_n )
    else
      allocate( fg_list( 1 ) )
    endif

    if( so_n .gt. 0 ) then
      allocate( so_list( 2, so_n ) )
      so_list( :, 1:so_n ) = so_list_temp( :, 1:so_n )
    else
      allocate( so_list( 1, 1 ) )
    endif
    
    deallocate( so_list_temp, fg_list_temp, ct_list_temp, proc_load, proc_list )    

!    do_staggered_sum = .true.
    do_staggered_sum = .false.

  end subroutine OCEAN_mult_create_par


  subroutine OCEAN_create_central( sys, ierr )
    use OCEAN_mpi
    use OCEAN_system
    implicit none
    
    type(O_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    integer :: lv, ic, icms, icml, ivms, ii, ivml, jj, nu, i, iband, ikpt, ip, ispn
    integer :: l, m, tmp_n
    real( DP ), allocatable :: pcr(:,:,:,:), pci(:,:,:,:), list(:)
    
    character(len=4) :: add04
    character(len=10) :: add10
    character(len=12) :: add12
    character(len=11) :: s11
    character(len=5) :: s5

    character(len=11) :: cks_filename
    character(len=5) :: cks_prefix

    ierr = 0
    if( is_init .and. associated( sys%cur_run%prev_run ) ) then
      if( sys%cur_run%ZNL(1) .ne. sys%cur_run%prev_run%ZNL(1) ) is_init = .false.
      if( sys%cur_run%ZNL(1) .ne. sys%cur_run%prev_run%ZNL(1) ) is_init = .false. 
      if( sys%cur_run%ZNL(1) .ne. sys%cur_run%prev_run%ZNL(1) ) is_init = .false. 
      if( sys%cur_run%indx .ne. sys%cur_run%prev_run%indx ) is_init = .false.
!      if( ( myid .eq. root ) .and. ( .not. is_init ) ) then
      if( .not. is_init )  then
        if( myid .eq. root ) write(6,*) 'Multiplets are reloading'
        deallocate( mpcr, mpci, mpm, mhr, mhi, cms, cml, vms, hcml, hvml, hcms, hvms, somelr, someli )
        deallocate( nproj, hvnu, ibeg, jbeg, mham )
        deallocate( ct_list, fg_list, so_list )
      endif
    endif

    if( is_init ) return

    write(add12 , '(A2,I2.2,A1,A1,I2.2,A1,I2.2)' ) sys%cur_run%elname, sys%cur_run%indx, &
            '_', 'n', sys%cur_run%ZNL(2), 'l', sys%cur_run%ZNL(3)

    if( myid .eq. root ) then
!      write(deflinz,'(a6,a1,i2.2,a1,i2.2,a1,i2.2)') 'deflin', 'z', sys%cur_run%ZNL(1), 'n', sys%cur_run%ZNL(2), 'l', sys%cur_run%ZNL(3)
!      open(unit=99,file=deflinz,form='formatted',status='old')
!      open(unit=99,file='derp_control',form='formatted', status='old' )
      open(unit=99,file='bse.in',form='formatted', status='old')
      rewind 99
      read( 99, * ) tmp_n, lvl, lvh
      read(99,*) xi
      close( 99 )
      xi = xi / 27.2114d0

      write( add04, '(1a1,1i3.3)' ) 'z', sys%cur_run%ZNL(1)
      write( add10, '(1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'z', sys%cur_run%ZNL(1), &
                  'n', sys%cur_run%ZNL(2), 'l', sys%cur_run%ZNL(3)
      write ( s11, '(1a7,1a4)' ) 'prjfile', add04
      write(6,*) s11
      open( unit=99, file=s11, form='formatted', status='old' )
      rewind 99
      read ( 99, * ) lmin, lmax
      allocate( nproj( lmin : lmax ) )
      read ( 99, * ) nproj
      close( unit=99 )
      if ( ( lvl .lt. lmin ) .or. ( lvh .gt. lmax ) ) stop 'bad lvl lvh'

      allocate( ibeg( lvl : lvh ), jbeg( lvl : lvh ), mham( lvl : lvh ) )
      itot = 0
      jtot = 0
      do lv = lvl, lvh
         ibeg( lv ) = itot + 1
         jbeg( lv ) = jtot + 1
         mham( lv ) = 4 * ( 2 * sys%cur_run%ZNL(3) + 1 ) * ( 2 * lv + 1 ) * nproj( lv )
         itot = itot + mham( lv )
         jtot = jtot + mham( lv ) ** 2
      end do

      allocate( cml( sys%cur_run%nalpha ), cms( sys%cur_run%nalpha ), vms( sys%cur_run%nalpha ) )
      allocate( hcml( itot ), hcms( itot ) )
      allocate( hvml( itot ), hvms( itot ) )
      allocate( hvnu( itot ) )
      allocate( mhr( jtot ), mhi( jtot ) )
!      allocate( pwr( itot ), pwi( itot ) )
!      allocate( hpwr( itot ), hpwi( itot ) )

      open( unit=99, file='channelmap', form='formatted', status='unknown' )
      rewind 99
      ic = 0
      do icms = -1, 1, 2
         do icml = -sys%cur_run%ZNL(3), sys%cur_run%ZNL(3)
            do ivms = -1, 1, 2
               ic = ic + 1
               cms( ic ) = 0.5d0 * icms
               cml( ic ) = icml
               vms( ic ) = 0.5d0 * ivms
               write ( 99, '(3f8.2,2i5,1a11)' ) cms( ic ), cml( ic ), vms( ic ), ic, sys%cur_run%nalpha, 'cms cml vms'
            end do
         end do
      end do
      !
      ii = 0
      do lv = lvl, lvh
         do ic = 1, 4 * ( 2 * sys%cur_run%ZNL(3) + 1 )
            do ivml = -lv, lv
               do nu = 1, nproj( lv )
                  ii = ii + 1
                  hcms( ii ) = cms( ic )
                  hcml( ii ) = cml( ic )
                  hvms( ii ) = vms( ic )
                  hvml( ii ) = ivml
                  hvnu( ii ) = nu
               end do
            end do
         end do
      end do
      !
      call nbsemkcmel( add04, add12 )
      do lv = lvl, lvh
         ii = ibeg( lv )
         jj = jbeg( lv )
!         call nbsemhsetup( sys%cur_run%ZNL(3), lv, nproj( lv ), mham( lv ), hcms( ii:ii+mham( lv )-1 ) , &
!              hcml( ii:ii+mham( lv )-1 ), hvms( ii:ii+mham( lv )-1 ), hvml( ii:ii+mham( lv )-1 ), hvnu( ii:ii+mham( lv )-1 ), &
!              mhr( jj:jj+mham( lv )*mham( lv )-1 ), mhi( jj:jj+mham( lv )*mham( lv )-1 ), add10 )
          call nbsemhsetup2( sys%cur_run%ZNL(3), lv, nproj( lv ), mham( lv ), & 
              hcms( ii:ii+mham(lv)-1 ) , hcml( ii:ii+mham(lv)-1 ), hvms( ii:ii+mham(lv)-1 ), &
              hvml( ii:ii+mham(lv)-1 ), hvnu( ii:ii+mham(lv)-1 ), &
              mhr( jj:jj-1+mham(lv)*mham(lv) ), mhi( jj:jj-1+mham(lv)*mham(lv) ), add10 )
      end do
      mhr = mhr / 27.2114d0
      mhi = mhi / 27.2114d0
      write ( 6, * ) 'multiplet hamiltonian set up'
      write(6,*) 'Number of spins = ', sys%nspn
      write ( 6, * ) 'n, nc, nspn', sys%num_bands*sys%nkpts, sys%cur_run%nalpha, sys%nspn


      npmax = maxval( nproj( lmin : lmax ) )
      allocate( mpcr( sys%num_bands, sys%nkpts, npmax, -lmax : lmax, lmin : lmax, sys%nspn ) )
      allocate( mpci( sys%num_bands, sys%nkpts, npmax, -lmax : lmax, lmin : lmax, sys%nspn ) )
!      allocate( ampr( npmax ), ampi( npmax ), hampr( npmax ), hampi( npmax ) )
      ! 
      ip = 0
      do l = lmin, lmax
         do m = -l, l
            do nu = 1, nproj( l )
               ip = ip + 1
            end do 
         end do
      end do
      nptot = ip
    ! Spin needed here too!
      allocate( pcr( nptot, sys%num_bands, sys%nkpts, sys%nspn ), pci( nptot, sys%num_bands, sys%nkpts, sys%nspn ) )
!JTV
      select case ( sys%cur_run%calc_type)
      case( 'XES' )
        cks_prefix = 'cksv.'
      case( 'XAS' )
        cks_prefix = 'cksc.'
      case default
        cks_prefix = 'cksc.'
      end select
      write(cks_filename,'(A5,A2,I4.4)' ) cks_prefix, sys%cur_run%elname, sys%cur_run%indx
      open(unit=99,file=cks_filename,form='unformatted', status='old' )
!      open( unit=99, file='ufmi', form='unformatted', status='old' )
      rewind 99
      read ( 99 )
      read ( 99 )
      read ( 99 ) pcr 
      read ( 99 ) pci
      do ispn = 1, sys%nspn
        do ikpt = 1, sys%nkpts
        do iband = 1, sys%num_bands 
          ip = 0
          do l = lmin, lmax
            do m = -l, l
               do nu = 1, nproj( l )
                  ip = ip + 1
                  mpcr( iband, ikpt, nu, m, l, ispn ) = pcr( ip, iband, ikpt, ispn )
                  mpci( iband, ikpt, nu, m, l, ispn ) = pci( ip, iband, ikpt, ispn )
               end do 
            end do
          end do
        end do
        enddo
      enddo
      close( unit=99 )
      deallocate( pcr, pci )
      write ( 6, * ) 'projector coefficients have been read in'
      !
      allocate( list( npmax * npmax ), mpm( npmax, npmax, lmin : lmax ) )
      mpm( :, :, : ) = 0.0_DP
      do l = lmin, lmax
         if ( l .gt. 9 ) stop 'l exceeds 9 in multip'
         write ( s5, '(1a4,1i1)' ) 'cmel', l
         open( unit=99, file=s5, form='formatted', status='old' )
         rewind 99
         read( 99, * ) list( 1: nproj( l ) * nproj( l ) )
         i = 0
         do nu = 1, nproj( l ) 
            mpm( 1:nproj(l), nu, l ) = list( i + 1 : i + nproj( l ) )
            i = i + nproj( l )
         end do
         close( unit=99 )
      end do
      mpm( :, :, : ) = mpm( :, :, : ) / 27.2114d0
      deallocate( list )
      write ( 6, * ) 'central projector matrix elements have been read in'
!      write(100,*) mpcr(:,1,1,0,0,1)
!      write(101,*) mpci(1,1,1,0,0,1)
!      write(102,*) mpm(:,1,0)

    endif

#ifdef MPI
    call MPI_BARRIER( comm, ierr )
    call MPI_BCAST( itot, 1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( jtot, 1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( npmax, 1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( lmin, 1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( lmax, 1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( lvl, 1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( lvh, 1, MPI_INTEGER, root, comm, ierr )


    if( myid .ne. root ) then
      allocate( cml( sys%cur_run%nalpha ), cms( sys%cur_run%nalpha ), vms( sys%cur_run%nalpha ) )
      allocate( hcml( itot ), hcms( itot ) )
      allocate( hvml( itot ), hvms( itot ) )
      allocate( hvnu( itot ) )
      allocate( mhr( jtot ), mhi( jtot ) )
      allocate( nproj( lmin : lmax ) )
      allocate( ibeg( lvl : lvh ), jbeg( lvl : lvh ), mham( lvl : lvh ) )
      allocate( mpcr( sys%num_bands, sys%nkpts, npmax, -lmax : lmax, lmin : lmax, sys%nspn ) )
      allocate( mpci( sys%num_bands, sys%nkpts, npmax, -lmax : lmax, lmin : lmax, sys%nspn ) )
      allocate( mpm( npmax, npmax, lmin : lmax ) )
    endif
    call MPI_BARRIER( comm, ierr )
    if( myid .eq. root ) write(6,*) 'HERE'
    call MPI_BARRIER( comm, ierr )
    call MPI_BCAST( cml, sys%cur_run%nalpha, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( cms, sys%cur_run%nalpha, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( vms, sys%cur_run%nalpha, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( hcml, itot, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( hcms, itot, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( hvml, itot, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( hvms, itot, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( hvnu, itot, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( mhr, jtot, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( mhi, jtot, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( nproj( lmin ), lmax-lmin+1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( ibeg( lvl ), lvh-lvl+1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( jbeg( lvl ), lvh-lvl+1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( mham( lvl ), lvh-lvl+1, MPI_INTEGER, root, comm, ierr )

    call MPI_BCAST( mpcr, sys%num_bands * sys%nkpts * npmax * (2*lmax+1) * (lmax-lmin+1) *sys%nspn, &
                    MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( mpci, sys%num_bands * sys%nkpts * npmax * (2*lmax+1) * (lmax-lmin+1) *sys%nspn, &
                    MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( mpm, npmax*npmax*(lmax-lmin+1), MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( myid .eq. root ) write( 6, * ) 'multiplet hamilotian shared'
#endif
  
    call OCEAN_soprep( sys, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_mult_create_par( sys, ierr )
    if( ierr .ne. 0 ) return

    is_init = .true.

  end subroutine OCEAN_create_central

  subroutine OCEAN_soprep( sys, ierr )
    use OCEAN_mpi
    use OCEAN_system
    implicit none
    
    type( O_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    integer :: ic, jc, i
    complex( kind = kind( 1.0d0 ) ) :: ctmp, rm1
    complex( kind = kind( 1.0d0 ) ), external :: jimel
    real( kind = kind( 1.d0 ) )  :: life_time( 2 ), l_alpha, l_beta, delta_so( 2 )
    logical :: broaden_exist
    !
    integer :: nsphpt, isphpt
    real( kind = kind( 1.0d0 ) ) :: sphsu
    real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: xsph, ysph, zsph, wsph

    ierr = 0
    allocate( somelr( sys%cur_run%nalpha, sys%cur_run%nalpha), someli( sys%cur_run%nalpha, sys%cur_run%nalpha) )

    if( myid .eq. root ) then

      write(6,*) 'SO PREP', xi

      open( unit=99, file='sphpts', form='formatted', status='old' )
      rewind 99
      read ( 99, * ) nsphpt
      allocate( xsph( nsphpt ), ysph( nsphpt ), zsph( nsphpt ), wsph( nsphpt ) )
      do isphpt = 1, nsphpt
         read ( 99, * ) xsph( isphpt ), ysph( isphpt ), zsph( isphpt ), wsph( isphpt )
      end do
      close( unit=99 )
      sphsu = sum( wsph( : ) )
      wsph( : ) = wsph( : ) * ( 4.0d0 * 4.0d0 * atan( 1.0d0 ) / sphsu )
      write ( 6, * ) nsphpt, ' points with weights summing to four pi '

      inquire( file='core_broaden.ipt', exist= broaden_exist )
      if( broaden_exist ) then
        write(6,*) 'spin-orbit dependent broadening'
        open( unit=99, file='core_broaden.ipt', form='formatted', status='old' )
        read( 99, * ) life_time( : )
        close( 99 )
        life_time( : ) = life_time( : ) / 27.2114d0

        delta_so( 1 ) = 0.5d0 * ( ( dble( sys%cur_run%ZNL(3) ) + 0.5 ) * ( dble( sys%cur_run%ZNL(3) ) + 1.5 ) - & 
                                    dble( sys%cur_run%ZNL(3) ) * dble( sys%cur_run%ZNL(3) + 1 ) - 0.75d0 )
        delta_so( 2 ) = 0.5d0 * ( ( dble( sys%cur_run%ZNL(3) ) - 0.5 ) * ( dble( sys%cur_run%ZNL(3) ) + 0.5 ) - &
                                    dble( sys%cur_run%ZNL(3) ) * dble( sys%cur_run%ZNL(3) + 1 ) - 0.75d0 )
        !
        l_alpha = ( life_time( 2 ) - life_time( 1 ) ) / ( - xi * ( delta_so( 2 ) - delta_so( 1 ) ) )
        l_beta = ( life_time( 2 ) * delta_so( 1 ) - life_time( 1 ) * delta_so( 2 ) ) / ( delta_so( 1 ) - delta_so( 2 ) )
        write(6,*) l_alpha, l_beta
      else
        l_alpha = 0.d0
        l_beta  = 0.d0
      endif
      somelr( :, : ) = 0.0_DP
      someli( :, : ) = 0.0_DP
      rm1 = -1; rm1 = sqrt( rm1 )
      do ic = 1, sys%cur_run%nalpha
         do jc = 1, sys%cur_run%nalpha
            if ( vms( ic ) .eq. vms( jc ) ) then
!????
!               call limel( sys%cur_run%ZNL(3), nint( cml( jc ) ), nint( cml( ic ) ), & 
!                           vrslt, nsphpt, xsph, ysph, zsph, wsph, prefs )
!????
               ctmp = 0
               do i = 1, 3
!????                  ctmp = ctmp + vrslt( i ) * jimel( 0.5d0, cms( jc ), cms( ic ), i )
                 ctmp = ctmp + jimel( dble( sys%cur_run%ZNL(3) ), cml( jc ), cml( ic ), i ) & 
                             * jimel( 0.5d0, cms( jc ), cms( ic ), i )
               end do
               write(20,*) ic, jc, real( ctmp ) !, aimag( ctmp )
               write(21,*) jimel( dble( sys%cur_run%ZNL(3) ), cml( jc ), cml( ic ), 1 ), &
                  sys%cur_run%ZNL(3), nint( cml( jc ) ), nint( cml( ic ) )
!               flush(20)
!               write(21,*) vrslt( : ), jimel( 0.5d0, cms( jc ), cms( ic ), 1 )
!               ctmp = -xi * ctmp
               somelr( ic, jc ) = -xi * real( ctmp ) !- aimag( ctmp ) * l_alpha
               someli( ic, jc ) = xi * aimag( ctmp ) !- real( ctmp ) * l_alpha
!               if( ic .eq. jc ) then
!                 someli( ic, jc ) = someli( ic, jc ) - l_beta !* real( ctmp )
!               endif
             end if
         end do
      end do
    endif

#ifdef MPI
    if( nproc .gt. 1 ) then
      call MPI_BCAST( someli, sys%cur_run%nalpha*sys%cur_run%nalpha, MPI_DOUBLE_PRECISION, root, comm, ierr )
      call MPI_BCAST( somelr, sys%cur_run%nalpha*sys%cur_run%nalpha, MPI_DOUBLE_PRECISION, root, comm, ierr )
    endif
#endif

  end subroutine OCEAN_soprep
    
  
  subroutine OCEAN_mult_act( sys, inter, in_vec, out_vec )
    use OCEAN_system
    use OCEAN_psi
    use OCEAN_mpi
    implicit none
    !
    type(O_system), intent( in ) :: sys
    real( DP ), intent( in ) :: inter
    type( ocean_vector ), intent( in ) :: in_vec
    type( ocean_vector ), intent( inout ) :: out_vec

!    integer :: take_longer

!    do take_longer = 1, 10

    out_vec%r(:,:,:) = 0.0_DP
    out_vec%i(:,:,:) = 0.0_DP

    if( do_staggered_sum ) then
      call OCEAN_ctact_dist( sys, inter, in_vec, out_vec )
!      call OCEAN_ctact_dist( sys, inter, in_vec%r, in_vec%i, out_vec%r, out_vec%i )
      call fgact_dist( sys, inter, in_vec, out_vec )
      if( sys%ZNL(2) .gt. 0 ) then
        call OCEAN_soact_dist( sys, in_vec, out_vec )
      endif
    elseif (myid .eq. 0 ) then
      call OCEAN_ctact( sys, inter, in_vec, out_vec )
      call fgact( sys, inter, in_vec, out_vec )
      if( sys%ZNL(2) .gt. 0 ) then
        call OCEAN_soact( sys, in_vec, out_vec )
      endif
    endif
!    call OCEAN_soact( sys, in_vec, out_vec )

!  enddo

  end subroutine OCEAN_mult_act


  subroutine OCEAN_ctact( sys, inter, in_vec, out_vec )
    use OCEAN_system
    use OCEAN_psi
    implicit none
    !
    type(O_system), intent( in ) :: sys
    real( DP ), intent( in ) :: inter
    type( ocean_vector ), intent( in ) :: in_vec
    type( ocean_vector ), intent( inout ) :: out_vec
    !
    !
    integer :: ialpha, l, m, nu, ispn, ikpt, ibnd
    real( DP ) :: mul
    real( DP ), dimension( npmax ) :: ampr, ampi, hampr, hampi
  !

    
    mul = inter / ( dble( sys%nkpts ) * sys%celvol )
!$OMP PARALLEL DO &
!$OMP DEFAULT( NONE ) &
!$OMP PRIVATE( ialpha, l, nu, m, ampr, ampi, hampr, hampi, ispn ) &
!$OMP SHARED( mul, sys, lmin, lmax, nproj ) &
!$OMP SHARED( mpcr, mpci, mpm, in_vec, out_vec )
    do ialpha = 1, sys%cur_run%nalpha
      ispn = 2 - mod( ialpha, 2 )
      if( sys%nspn .eq. 1 ) then
        ispn = 1
      endif
      do l = lmin, lmax
        do m = -l, l
          ampr( : ) = 0
          ampi( : ) = 0
          do nu = 1, nproj( l )
            do ikpt = 1, sys%nkpts
              do ibnd = 1, sys%num_bands
                ampr( nu ) = ampr( nu ) &
                           + in_vec%r( ibnd, ikpt, ialpha ) * mpcr( ibnd, ikpt, nu, m, l, ispn ) &
                           - in_vec%i( ibnd, ikpt, ialpha ) * mpci( ibnd, ikpt, nu, m, l, ispn )
                ampi( nu ) = ampi( nu ) &
                           + in_vec%r( ibnd, ikpt, ialpha ) * mpci( ibnd, ikpt, nu, m, l, ispn ) &
                           + in_vec%i( ibnd, ikpt, ialpha ) * mpcr( ibnd, ikpt, nu, m, l, ispn )
               enddo
            enddo
          enddo


          hampr( : ) = 0
          hampi( : ) = 0
          do nu = 1, nproj( l )
            hampr( : ) = hampr( : ) - mpm( :, nu, l ) * ampr( nu ) 
            hampi( : ) = hampi( : ) - mpm( :, nu, l ) * ampi( nu ) 
          enddo
          hampr( : ) = hampr( : ) * mul
          hampi( : ) = hampi( : ) * mul
          do nu = 1, nproj( l )
            do ikpt = 1, sys%nkpts 
              do ibnd = 1, sys%num_bands
                out_vec%r( ibnd, ikpt, ialpha ) = out_vec%r( ibnd, ikpt, ialpha ) &
                                                + mpcr( ibnd, ikpt, nu, m, l, ispn ) * hampr( nu )  &
                                                + mpci( ibnd, ikpt, nu, m, l, ispn ) * hampi( nu )
                out_vec%i( ibnd, ikpt, ialpha ) = out_vec%i( ibnd, ikpt, ialpha ) &
                                                + mpcr( ibnd, ikpt, nu, m, l, ispn ) * hampi( nu ) &
                                                - mpci( ibnd, ikpt, nu, m, l, ispn ) * hampr( nu )
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine OCEAN_ctact

  subroutine OCEAN_ctact_dist( sys, inter, in_vec, out_vec )
    use OCEAN_system
    use OCEAN_psi
    implicit none
    !
    type(O_system), intent( in ) :: sys
    real( DP ), intent( in ) :: inter
    type( ocean_vector ), intent( in ) :: in_vec
    type( ocean_vector ), intent( inout ) :: out_vec
!    real( DP ), dimension( sys%num_bands, sys%nkpts, sys%nalpha ), intent( in ) :: &
!      in_vec_r, in_vec_i
!    real( DP ), dimension( sys%num_bands, sys%nkpts, sys%nalpha ), intent( inout ) :: &
!      out_vec_r, out_vec_i

    !
    !
    integer :: ialpha, l, m, nu, ispn, ikpt, ibnd, i
    real( DP ) :: mul
    real( DP ), dimension( npmax ) :: ampr, ampi, hampr, hampi
    !
#ifdef __INTEL_COMPILER
!DIR$ attributes align: 64 :: ampr, ampi, hampr, hampi
#endif
    !
    !
    if( ct_n .lt. 1 ) return

    mul = inter / ( dble( sys%nkpts ) * sys%celvol )


!$OMP PARALLEL &
!$OMP DEFAULT( NONE ) &
!$OMP PRIVATE( i, ialpha, l, m, ispn, nu, ikpt, ibnd ) &
!$OMP SHARED( ampr, ampi, in_vec, mpcr, mpci, sys, out_vec ) &
!$OMP SHARED( ct_n, ct_list, mul, nproj, mpm, hampi, hampr ) 


      do i = 1, ct_n

        ialpha = ct_list( 1, i )
        l = ct_list( 2, i )
        m = ct_list( 3, i )
        ispn = 2 - mod( ialpha, 2 )
        if( sys%nspn .eq. 1 ) then
          ispn = 1
        endif

        ampr( : ) = 0
        ampi( : ) = 0
        do nu = 1, nproj( l )
          call ctact_dist_in( sys%nkpts, sys%num_bands, ampr( nu ), ampi( nu ), &
                              in_vec%r( :, :, ialpha ), in_vec%i( :, :, ialpha ), &
                              mpcr( :, :, nu, m, l, ispn ), mpci( :, :, nu, m, l, ispn ) )
        enddo

!$OMP SINGLE
        hampr( : ) = 0
        hampi( : ) = 0
        do nu = 1, nproj( l )
          hampr( : ) = hampr( : ) - mpm( :, nu, l ) * ampr( nu )
          hampi( : ) = hampi( : ) - mpm( :, nu, l ) * ampi( nu )
        enddo
        hampr( : ) = hampr( : ) * mul
        hampi( : ) = hampi( : ) * mul
!$OMP END SINGLE


        do nu = 1, nproj( l )
          call ctact_dist_out( sys%nkpts, sys%num_bands, hampr(nu), hampi(nu), &
                               out_vec%r( :, :, ialpha ), out_vec%i( :, :, ialpha ), &
                               mpcr( :, :, nu, m, l, ispn ), mpci( :, :, nu, m, l, ispn ) )
        enddo

    enddo

!$OMP END PARALLEL 

  end subroutine OCEAN_ctact_dist

  pure subroutine ctact_dist_in( nkpts, num_bands, ampr, ampi, in_vec_r, in_vec_i, &
                                 sm_mpcr, sm_mpci )
    implicit none
    !
    integer, intent( in ) :: nkpts, num_bands
    real(DP), intent( out ) :: ampr, ampi
    real(DP), dimension(nkpts*num_bands), intent( in ) :: in_vec_r, in_vec_i, sm_mpcr, sm_mpci
    !
    integer :: iter
    !
!$OMP DO SIMD &
!$OMP SCHEDULE( STATIC ) &
!$OMP ALIGNED(in_vec_r,in_vec_i,sm_mpcr,sm_mpci: 32 ) &
!$OMP REDUCTION(+:ampr,ampi)
    do iter = 1, num_bands * nkpts
      ampr = ampr + in_vec_r( iter ) * sm_mpcr( iter )  &
           - in_vec_i( iter ) * sm_mpci( iter )
      ampi = ampi + in_vec_r( iter ) * sm_mpci( iter ) &
           + in_vec_i( iter ) * sm_mpcr( iter )
    enddo
!$OMP END DO SIMD

  end subroutine ctact_dist_in

  pure subroutine ctact_dist_out( nkpts, num_bands, hampr, hampi, out_vec_r, out_vec_i, &
                                  sm_mpcr, sm_mpci )
    implicit none
    !
    integer, intent( in ) :: nkpts, num_bands
    real(DP), intent( in ) :: hampr, hampi
    real(DP), dimension( nkpts * num_bands ), intent( in ) :: sm_mpcr, sm_mpci
    real(DP), dimension( nkpts * num_bands ), intent(inout) :: out_vec_r, out_vec_i
    !
    integer :: iter
    
!$OMP DO SIMD &
!$OMP SCHEDULE( STATIC ) &
!$OMP ALIGNED(out_vec_r,out_vec_i,sm_mpcr,sm_mpci: 32 ) 
      do iter = 1, nkpts * num_bands
        out_vec_r( iter ) = out_vec_r( iter ) + sm_mpcr( iter ) * hampr &
                          + sm_mpci( iter ) * hampi
        out_vec_i( iter ) = out_vec_i( iter ) + sm_mpcr( iter ) * hampi &
                          - sm_mpci( iter ) * hampr
      enddo
!$OMP END DO SIMD

  
  end subroutine ctact_dist_out


  subroutine fgact( sys, inter, in_vec, out_vec )
    use OCEAN_system
    use OCEAN_psi
    implicit none
    !
    type(O_system), intent( in ) :: sys
    real( DP ), intent( in ) :: inter
    type( ocean_vector ), intent( in ) :: in_vec
    type( ocean_vector ), intent( inout ) :: out_vec
    !
!    real( DP ), dimension( n, nc, 2 ) :: v, hv
!    real( DP ), dimension( n, npmax, -lmax : lmax, lmin : lmax, nspn ) :: mpcr, mpci
!    real( DP ), dimension( jtot ) :: mhr, mhi
    !
    integer :: lv, ii, ic, ivml, nu, j1, jj, ispn, ikpt, ibnd
    real( DP ) :: mul
    real( DP ), allocatable, dimension( : ) :: pwr, pwi, hpwr, hpwi
    !
    allocate( pwr(itot), pwi( itot), hpwr(itot), hpwi(itot) )
    ! lvl, lvh is not enough for omp
    ! should probably pull thread spawning out of the do loop though

    mul = inter / ( dble( sys%nkpts ) * sys%celvol )
!    write(6,*) mul
    do lv = lvl, lvh
       pwr( : ) = 0
       pwi( : ) = 0
       hpwr( : ) = 0
       hpwi( : ) = 0
!$OMP PARALLEL &
!$OMP PRIVATE( ic, ivml, nu, ii, jj, j1, ispn ) &
!$OMP FIRSTPRIVATE( lv, jbeg, mham, nproj, mul ) &
!$OMP SHARED( pwr, pwi, mhr, mhi, mpcr, mpci, hpwr, hpwi, in_vec, out_vec, sys ) &
!$OMP DEFAULT( NONE )

!$OMP DO 
      do ic = 1, sys%cur_run%nalpha
        ispn = 2 - mod( ic, 2 )
  !       ispn = 1 + mod( ic, 2 )
        if( sys%nspn .eq. 1 ) then
          ispn = 1
        endif
        do ivml = -lv, lv
          do nu = 1, nproj( lv )
            ii = nu + ( ivml + lv ) * nproj( lv ) + ( ic - 1 ) * ( 2 * lv + 1 ) * nproj( lv )
            do ikpt = 1, sys%nkpts
              do ibnd = 1, sys%num_bands
                pwr( ii ) = pwr( ii ) &
                          + in_vec%r( ibnd, ikpt, ic ) * mpcr( ibnd, ikpt, nu, ivml, lv, ispn ) &
                          - in_vec%i( ibnd, ikpt, ic ) * mpci( ibnd, ikpt, nu, ivml, lv, ispn ) 
                pwi( ii ) = pwi( ii ) &
                          + in_vec%r( ibnd, ikpt, ic ) * mpci( ibnd, ikpt, nu, ivml, lv, ispn ) &
                          + in_vec%i( ibnd, ikpt, ic ) * mpcr( ibnd, ikpt, nu, ivml, lv, ispn ) 
              enddo
             enddo
           enddo
         enddo
       enddo
!$OMP END DO


!$OMP DO
      do ii = 1, mham( lv )
        do jj = 1, mham( lv )
          j1 = jbeg( lv ) + ( jj -1 ) + ( ii - 1 ) * mham( lv )
          hpwr( ii ) = hpwr( ii ) + mhr( j1 ) * pwr( jj ) - mhi( j1 ) * pwi( jj )
          hpwi( ii ) = hpwi( ii ) + mhr( j1 ) * pwi( jj ) + mhi( j1 ) * pwr( jj )
        end do
        hpwr( ii ) = hpwr( ii ) * mul
        hpwi( ii ) = hpwi( ii ) * mul
      end do
!$OMP END DO


!$OMP DO
      do ic = 1, sys%cur_run%nalpha
        ispn = 2 - mod( ic, 2 )
  !       ispn = 1 + mod( ic, 2 )
        if( sys%nspn .eq. 1 ) then
          ispn = 1
        endif
        do ivml = -lv, lv
          do nu = 1, nproj( lv )
            ii = nu + ( ivml + lv ) * nproj( lv ) + ( ic - 1 ) * ( 2 * lv + 1 ) * nproj( lv )
            do ikpt = 1, sys%nkpts
              do ibnd = 1, sys%num_bands
                out_vec%r( ibnd, ikpt, ic ) = out_vec%r( ibnd, ikpt, ic )  &
                                            + hpwr( ii ) * mpcr( ibnd, ikpt, nu, ivml, lv, ispn ) &
                                            + hpwi( ii ) * mpci( ibnd, ikpt, nu, ivml, lv, ispn )
                out_vec%i( ibnd, ikpt, ic ) = out_vec%i( ibnd, ikpt, ic )  &
                                            + hpwi( ii ) * mpcr( ibnd, ikpt, nu, ivml, lv, ispn ) &
                                            - hpwr( ii ) * mpci( ibnd, ikpt, nu, ivml, lv, ispn )
              enddo
            enddo
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL 
    end do
    !

  end subroutine fgact


  subroutine fgact_dist( sys, inter, in_vec, out_vec )
    use OCEAN_system
    use OCEAN_psi
    implicit none
    !
    type(O_system), intent( in ) :: sys
    real( DP ), intent( in ) :: inter
    type( ocean_vector ), intent( in ) :: in_vec
    type( ocean_vector ), intent( inout ) :: out_vec
    !
    !
    integer :: lv, ii, ic, ivml, nu, j1, jj, ispn, i
    real( DP ) :: mul
    real( DP ), allocatable, dimension( : ) :: pwr, pwi, hpwr, hpwi
    !
#ifdef __INTEL_COMPILER
!DIR$ attributes align: 64 :: pwr, pwi, hpwr, hpwi
#endif
    !
    !


    if( fg_n .lt. 1 ) return
    
    allocate( pwr(itot), pwi( itot), hpwr(itot), hpwi(itot) )
    ! lvl, lvh is not enough for omp
    ! should probably pull thread spawning out of the do loop though

    mul = inter / ( dble( sys%nkpts ) * sys%celvol )


!$OMP PARALLEL &
!$OMP PRIVATE( ic, ivml, nu, ii, jj, j1, ispn, i, lv ) &
!$OMP SHARED( mul, nproj, mham, jbeg, pwr, pwi, mhr, mhi, mpcr, mpci, hpwr, hpwi ) &
!$OMP SHARED( in_vec, out_vec, sys, fg_n, fg_list ) &
!$OMP DEFAULT( NONE )
  
    do i = 1, fg_n
      lv = fg_list( i )
      pwr( : ) = 0
      pwi( : ) = 0
      hpwr( : ) = 0
      hpwi( : ) = 0

      do ic = 1, sys%cur_run%nalpha
        ispn = 2 - mod( ic, 2 )
        if( sys%nspn .eq. 1 ) then
          ispn = 1
        endif
        do ivml = -lv, lv
          do nu = 1, nproj( lv )
            ii = nu + ( ivml + lv ) * nproj( lv ) + ( ic - 1 ) * ( 2 * lv + 1 ) * nproj( lv )
    
            call ctact_dist_in( sys%nkpts, sys%num_bands, pwr( ii ), pwi( ii ), &
                                in_vec%r( :, :, ic ), in_vec%i( :, :, ic ), &
                              mpcr( :, :, nu, ivml, lv, ispn ), mpci( :, :, nu, ivml, lv, ispn ) )
           enddo
         enddo
       enddo


!$OMP DO
      do ii = 1, mham( lv )
        do jj = 1, mham( lv )
          j1 = jbeg( lv ) + ( jj -1 ) + ( ii - 1 ) * mham( lv )
          hpwr( ii ) = hpwr( ii ) + mhr( j1 ) * pwr( jj ) - mhi( j1 ) * pwi( jj )
          hpwi( ii ) = hpwi( ii ) + mhr( j1 ) * pwi( jj ) + mhi( j1 ) * pwr( jj )
        end do
        hpwr( ii ) = hpwr( ii ) * mul
        hpwi( ii ) = hpwi( ii ) * mul
      end do
!$OMP END DO


      do ic = 1, sys%cur_run%nalpha
        ispn = 2 - mod( ic, 2 )
        if( sys%nspn .eq. 1 ) then
          ispn = 1
        endif
        do ivml = -lv, lv
          do nu = 1, nproj( lv )
            ii = nu + ( ivml + lv ) * nproj( lv ) + ( ic - 1 ) * ( 2 * lv + 1 ) * nproj( lv )
            call ctact_dist_out( sys%nkpts, sys%num_bands, hpwr( ii ), hpwi( ii ), &
                                 out_vec%r( :, :, ic ), out_vec%i( :, :, ic ), &
                               mpcr( :, :, nu, ivml, lv, ispn ), mpci( :, :, nu, ivml, lv, ispn ) )
          enddo
        enddo
      enddo
    end do
!$OMP END PARALLEL 
    !

  end subroutine fgact_dist

  subroutine OCEAN_soact_dist( sys, in_vec, out_vec )
    use OCEAN_system
    use OCEAN_psi
    implicit none

    type( O_system ), intent( in ) :: sys
    type( OCEAN_vector ), intent( in ) :: in_vec
    type( OCEAN_vector ), intent( inout ) :: out_vec

    integer :: ic, jc, ikpt, iter
    real( DP ) :: melr, meli

    if( so_n .lt. 1 ) return

!JTV Maybe want to consider blocking this routine for cache reuse

! !$OMP PARALLEL &
! !$OMP DEFAULT( NONE ) &
! !$OMP PRIVATE( iter, ic, jc, ikpt, melr, meli ) &
! !$OMP SHARED( sys, in_vec, out_vec, somelr, someli, so_n, so_list, vms )
    do iter = 1, so_n
      ic = so_list( 1, iter )
      jc = so_list( 2, iter )
      melr = somelr( ic, jc )
      meli = someli( ic, jc )
      if ( vms( ic ) .eq. vms( jc ) ) then
! !$OMP DO 
        do ikpt = 1, sys%nkpts
          out_vec%r( 1:sys%num_bands, ikpt, ic ) = out_vec%r( 1:sys%num_bands, ikpt, ic )  &
                                   + melr * in_vec%r( 1:sys%num_bands, ikpt, jc ) &
                                   - meli * in_vec%i( 1:sys%num_bands, ikpt, jc )
!                                   + somelr( ic, jc ) * in_vec%r( 1:sys%num_bands, ikpt, jc ) &
!                                   - someli( ic, jc ) * in_vec%i( 1:sys%num_bands, ikpt, jc )
          out_vec%i( 1:sys%num_bands, ikpt, ic ) = out_vec%i( 1:sys%num_bands, ikpt, ic )  &
                                   + melr * in_vec%i( 1:sys%num_bands, ikpt, jc ) &
                                   + meli * in_vec%r( 1:sys%num_bands, ikpt, jc )
!                                   + somelr( ic, jc ) * in_vec%i( 1:sys%num_bands, ikpt, jc ) &
!                                   + someli( ic, jc ) * in_vec%r( 1:sys%num_bands, ikpt, jc )
        enddo
! !$OMP END DO
      endif
    enddo
! !$OMP END PARALLEL


  end subroutine OCEAN_soact_dist




  subroutine OCEAN_soact( sys, in_vec, out_vec )
    use OCEAN_system
    use OCEAN_psi
    implicit none

    type( O_system ), intent( in ) :: sys
    type( OCEAN_vector ), intent( in ) :: in_vec
    type( OCEAN_vector ), intent( inout ) :: out_vec

    integer :: ic, jc, ikpt
!    complex( DP ) :: ctmp
!    real( DP ) :: melr, meli

! !$OMP PARALLEL DO &
! !$OMP DEFAULT( NONE ) &
! !$OMP PRIVATE( ic, jc, ikpt ) &
! !$OMP SHARED( sys, in_vec, out_vec, somelr, someli )
    do ic = 1, sys%cur_run%nalpha
      do jc = 1, sys%cur_run%nalpha
        do ikpt = 1, sys%nkpts
          out_vec%r( 1:sys%num_bands, ikpt, ic ) = out_vec%r( 1:sys%num_bands, ikpt, ic )  & 
                                   + somelr( ic, jc ) * in_vec%r( 1:sys%num_bands, ikpt, jc ) &
                                   - someli( ic, jc ) * in_vec%i( 1:sys%num_bands, ikpt, jc )
          out_vec%i( 1:sys%num_bands, ikpt, ic ) = out_vec%i( 1:sys%num_bands, ikpt, ic )  &
                                   + somelr( ic, jc ) * in_vec%i( 1:sys%num_bands, ikpt, jc ) &
                                   + someli( ic, jc ) * in_vec%r( 1:sys%num_bands, ikpt, jc )
        enddo
      enddo
    enddo
! !$OMP END PARALLEL DO
!    do ic = 1, sys%cur_run%nalpha
!      do jc = 1, sys%cur_run%nalpha
!        if( vms( ic ) .eq. vms( jc ) ) then
!          ctmp = 0.0_DP
!          do i = 1, 3
!            ctmp = ctmp + jimel( dble( lc ), cml( jc ), cml( ic ), i ) &
!                        * jimel( 0.5_DP, cms( jc ), cms( ic ), i )
!          enddo
!          melr = real(ctmp * -xi)
!          meli = aimag( ctmp * -xi )
!          do ikpt = 1, sys%nkpts
!            out_vec%r( 1:sys%num_bands, ikpt, ic ) = out_vec%r( 1:sys%num_bands, ikpt, ic ) &
!                                                   + melr( 



  end subroutine OCEAN_soact


  subroutine nbsemhsetup2( lc, lv, np, mham, cms, cml, vms, vml, vnu, mhr, mhi, add10 )
    implicit none
    !
    integer :: lc, lv, np, mham
    real( kind = kind( 1.0d0 ) ) :: cms( mham ), cml( mham )
    real( kind = kind( 1.0d0 ) ) :: vms( mham ), vml( mham )
    real( kind = kind( 1.0d0 ) ) :: mhr( mham, mham ), mhi( mham, mham )
    integer :: vnu( mham )
    !
    integer :: npt
    real( kind = kind( 1.0d0 ) ) :: pi, su, yp( 0 : 1000 )
    complex( kind = kind( 1.0d0 ) ) :: rm1
    real( kind = kind( 1.0d0 ) ), allocatable :: x( : ), w( : )
    !
    integer :: kfl, kfh, kgl, kgh
    real( kind = kind( 1.0d0 ) ), allocatable :: fk( :, :, : ), scfk( : )
    real( kind = kind( 1.0d0 ) ), allocatable :: gk( :, :, : ), scgk( : )
    !
    integer :: i, i1, i2, nu1, nu2
    integer :: l1, m1, s1, l2, m2, s2, l3, m3, s3, l4, m4, s4, k, mk
    real( kind = kind( 1.0d0 ) ) :: ggk, ffk
    complex( kind = kind( 1.0d0 ) ) :: f1, f2, ctmp
    logical, parameter :: no = .false., yes = .true.
    logical :: tdlda
    !
    character * 10 :: add10
    character * 15 :: filnam
    !
    include 'sphsetnx.h.f90'
    !
    include 'sphsetx.h.f90'
    !
     write(6,*) 'nbsemhsetup', lv
    call newgetprefs( yp, max( lc, lv ), nsphpt, wsph, xsph, ysph, zsph )
    rm1 = -1
    rm1 = sqrt( rm1 )
    pi = 4.0d0 * atan( 1.0d0 )
    !
    write ( 6, * ) ' add10 = ', add10 
    open( unit=99, file='Pquadrature', form='formatted', status='old' )
    rewind 99
    read ( 99, * ) npt
    allocate( x( npt ), w( npt ) )
    su = 0
    do i = 1, npt
       read ( 99, * ) x( i ), w( i )
       su = su + w( i )
    end do
    close( unit=99 )
    w = w * 2 / su
    !
    if ( lv .gt. 9 ) stop 'lv must be single digit'
    if ( lc .gt. 9 ) stop 'lc must be single digit'
    !
    ! TDLDA
    inquire( file='tdlda', exist=tdlda )
    if( tdlda ) then
      open(unit=99, file='tdlda', form='formatted', status='old' )
      read( 99, * ) tdlda
      close( 99 )
    endif
    kfh = min( 2 * lc, 2 * lv )
    kfl = 2
    if( tdlda ) then
       kfl = 0
       kfh = 0
       allocate( fk( np, np, kfl : kfh ), scfk( kfl : kfh ) )
       fk = 0; scfk = 0
       do k = kfl, kfh, 2
          write ( filnam, '(1a2,3i1,1a10)' ) 'kk', lc, lv, k, add10
          write( 6, * ) filnam
          open( unit=99, file=filnam, form='formatted', status='old' )
          rewind 99
          read ( 99, * ) fk( :, :, k ), scfk( k )
          close( unit=99 )
       end do
    else
    if ( kfh .ge. kfl ) then
       allocate( fk( np, np, kfl : kfh ), scfk( kfl : kfh ) )
       fk = 0; scfk = 0
       do k = kfl, kfh, 2
          write ( filnam, '(1a2,3i1,1a10)' ) 'fk', lc, lv, k, add10
          open( unit=99, file=filnam, form='formatted', status='old' )
          rewind 99
          read ( 99, * ) fk( :, :, k ), scfk( k )
          close( unit=99 )
       end do
    end if
    end if
    !
    kgh = lc + lv
    kgl = abs( lc - lv )
    if ( kgh .ge. kgl ) then
       allocate( gk( np, np, kgl : kgh ), scgk( kgl : kgh ) )
       gk = 0; scgk = 0
       do k = kgl, kgh, 2
          write ( filnam, '(1a2,3i1,1a10)' ) 'gk', lc, lv, k, add10
          open( unit=99, file=filnam, form='formatted', status='old' )
          rewind 99
          read ( 99, * ) gk( :, :, k ), scgk( k )
        close( unit=99 )
       end do
    end if
    !
    mhr = 0
    mhi = 0
    do i1 = 1, mham
       nu1 = vnu( i1 )
       do i2 = 1, mham
          nu2 = vnu( i2 )
          !
          if( kfh .ge. kfl ) then
             l1 = lc; m1 = nint( cml( i2 ) ); s1 = nint( 2 * cms( i2 ) )
             l2 = lv; m2 = nint( vml( i1 ) ); s2 = nint( 2 * vms( i1 ) )
             l3 = lc; m3 = nint( cml( i1 ) ); s3 = nint( 2 * cms( i1 ) )
             l4 = lv; m4 = nint( vml( i2 ) ); s4 = nint( 2 * vms( i2 ) )
             if ( ( s1 .eq. s3 ) .and. ( s2 .eq. s4 ) ) then
                mk = m1 - m3
                if ( m1 + m2 .eq. m3 + m4 ) then
                   do k = kfl, kfh, 2
                      if ( abs( mk ) .le. k ) then
                         ffk = scfk( k ) * fk( nu1, nu2, k ) 
                         call threey( l1, m1, k, mk, l3, m3, no, npt, x, w, yp, f1 )
                         call threey( l2, m2, k, mk, l4, m4, yes, npt, x, w, yp, f2 )
                         ctmp = - ffk * f1 * f2 * ( 4 * pi / ( 2 * k + 1 ) )
                         mhr( i1, i2 ) = mhr( i1, i2 ) + real(ctmp)
                         mhi( i1, i2 ) = mhi( i1, i2 ) + aimag(ctmp) !- ctmp * rm1
                      end if
                   end do
                end if
             end if
          end if
          !
          if( kgh .ge. kgl ) then
             l1 = lc; m1 = nint( cml( i2 ) ); s1 = nint( 2 * cms( i2 ) )
             l2 = lv; m2 = nint( vml( i1 ) ); s2 = nint( 2 * vms( i1 ) )
             l3 = lv; m3 = nint( vml( i2 ) ); s3 = nint( 2 * vms( i2 ) )
             l4 = lc; m4 = nint( cml( i1 ) ); s4 = nint( 2 * cms( i1 ) )
             if ( ( s1 .eq. s3 ) .and. ( s2 .eq. s4 ) ) then
                mk = m1 - m3
                if ( m1 + m2 .eq. m3 + m4 ) then
                   do k = kgl, kgh, 2
                      if ( abs( mk ) .le. k ) then
                         ggk = scgk( k ) * gk( nu1, nu2, k ) 
                         call threey( l1, m1, k, mk, l3, m3, no, npt, x, w, yp, f1 )
                         call threey( l2, m2, k, mk, l4, m4, yes, npt, x, w, yp, f2 )
                         ctmp = ggk * f1 * f2 * ( 4 * pi / ( 2 * k + 1 ) )
                         mhr( i1, i2 ) = mhr( i1, i2 ) + real(ctmp)
                         mhi( i1, i2 ) = mhi( i1, i2 ) + aimag(ctmp)! - ctmp * rm1
                      end if
                   end do
                end if
             end if
          end if
          !
       end do
    end do
    !
    return
  end subroutine nbsemhsetup2


  subroutine OCEAN_mult_single( sys, bse_me, inter, ib, ik, ia, jb, jk, ja )
    use OCEAN_system
    implicit none
    !
    type(O_system), intent( in ) :: sys
    real(DP), intent( in ) :: inter
    integer, intent( in ) :: ib, ik, ia, jb, jk, ja
    complex(DP), intent( inout ) :: bse_me

    call OCEAN_ct_single( sys, bse_me, inter, ib, ik, ia, jb, jk, ja )
    call OCEAN_fg_single( sys, bse_me, inter, ib, ik, ia, jb, jk, ja )

!    call OCEAN_soact( sys, in_vec, out_vec )


  end subroutine OCEAN_mult_single

  subroutine OCEAN_mult_slice( sys, out_vec, inter, jb, jk, ja )
    use OCEAN_system
    use OCEAN_psi
    implicit none
    !
    type(O_system), intent( in ) :: sys
    type( ocean_vector ), intent( inout ) :: out_vec
    real(DP), intent( in ) :: inter
    integer, intent( in ) :: jb, jk, ja


    call OCEAN_ct_slice( sys, out_vec, inter, jb, jk, ja )
    call OCEAN_fg_slice( sys, out_vec, inter, jb, jk, ja )

  end subroutine OCEAN_mult_slice

  subroutine OCEAN_fg_slice( sys, out_vec, inter, jb, jk, ja )
    use OCEAN_system
    use OCEAN_psi
    implicit none
    !
    type(O_system), intent( in ) :: sys
    type( ocean_vector ), intent( inout ) :: out_vec
    real(DP), intent( in ) :: inter
    integer, intent( in ) :: jb, jk, ja
    !
    !
    integer :: lv, ii, ivml, nu, j1, jj, ispn, ia, ibnd, ikpt
    real( DP ) :: mul
    real( DP ), allocatable, dimension( : ) :: pwr, pwi, hpwr, hpwi
    !
    allocate( pwr(itot), pwi( itot), hpwr(itot), hpwi(itot) )
    ! lvl, lvh is not enough for omp
    ! should probably pull thread spawning out of the do loop though

!JTV clean this up too
    mul = inter / ( dble( sys%nkpts ) * sys%celvol )

    ispn = 2 - mod( ja, 2 )
    if( sys%nspn .eq. 1 ) then
      ispn = 1
    endif

    do lv = lvl, lvh
      pwr( : ) = 0.0_DP
      pwi( : ) = 0.0_DP
      hpwr( : ) = 0.0_DP
      hpwi( : ) = 0.0_DP


      do ivml = -lv, lv
        do nu = 1, nproj( lv )
          ii = nu + ( ivml + lv ) * nproj( lv ) + ( ja - 1 ) * ( 2 * lv + 1 ) * nproj( lv )

          pwr( ii ) = pwr( ii ) + mpcr( jb, jk, nu, ivml, lv, ispn )
          pwi( ii ) = pwi( ii ) + mpci( jb, jk, nu, ivml, lv, ispn )

         enddo
       enddo

      do ii = 1, mham( lv )
        do jj = 1, mham( lv )
          j1 = jbeg( lv ) + ( jj -1 ) + ( ii - 1 ) * mham( lv )
          hpwr( ii ) = hpwr( ii ) + mhr( j1 ) * pwr( jj ) - mhi( j1 ) * pwi( jj )
          hpwi( ii ) = hpwi( ii ) + mhr( j1 ) * pwi( jj ) + mhi( j1 ) * pwr( jj )
        end do
        hpwr( ii ) = hpwr( ii ) * mul
        hpwi( ii ) = hpwi( ii ) * mul
      end do

      do ia = 1, sys%nalpha !ja-1
        ispn = 2 - mod( ia, 2 )
        if( sys%nspn .eq. 1 ) then
          ispn = 1
        endif
        do ivml = -lv, lv
          do nu = 1, nproj( lv )
          ii = nu + ( ivml + lv ) * nproj( lv ) + ( ia - 1 ) * ( 2 * lv + 1 ) * nproj( lv )
            do ikpt = 1, sys%nkpts
              do ibnd = 1, sys%num_bands
                out_vec%r( ibnd, ikpt, ia ) = out_vec%r( ibnd, ikpt, ia )  &
                                            + hpwr( ii ) * mpcr( ibnd, ikpt, nu, ivml, lv, ispn ) &
                                            + hpwi( ii ) * mpci( ibnd, ikpt, nu, ivml, lv, ispn )
                out_vec%i( ibnd, ikpt, ia ) = out_vec%i( ibnd, ikpt, ia )  &
                                            + hpwi( ii ) * mpcr( ibnd, ikpt, nu, ivml, lv, ispn ) &
                                            - hpwr( ii ) * mpci( ibnd, ikpt, nu, ivml, lv, ispn )
              enddo
            enddo
          enddo
        enddo
      enddo

!      ia = ja
!      ispn = 2 - mod( ia, 2 )
!      if( sys%nspn .eq. 1 ) then
!        ispn = 1
!      endif
!      do ivml = -lv, lv
!        do nu = 1, nproj( lv )
!        ii = nu + ( ivml + lv ) * nproj( lv ) + ( ia - 1 ) * ( 2 * lv + 1 ) * nproj( lv )
!          do ikpt = 1, jk-1
!            do ibnd = 1, sys%num_bands
!              out_vec%r( ibnd, ikpt, ia ) = out_vec%r( ibnd, ikpt, ia )  &
!                                          + hpwr( ii ) * mpcr( ibnd, ikpt, nu, ivml, lv, ispn ) &
!                                          + hpwi( ii ) * mpci( ibnd, ikpt, nu, ivml, lv, ispn )
!              out_vec%i( ibnd, ikpt, ia ) = out_vec%i( ibnd, ikpt, ia )  &
!                                          + hpwi( ii ) * mpcr( ibnd, ikpt, nu, ivml, lv, ispn ) &
!                                          - hpwr( ii ) * mpci( ibnd, ikpt, nu, ivml, lv, ispn )
!            enddo
!          enddo
!
!          ikpt = jk
!            do ibnd = 1, jb
!              out_vec%r( ibnd, ikpt, ia ) = out_vec%r( ibnd, ikpt, ia )  & 
!                                          + hpwr( ii ) * mpcr( ibnd, ikpt, nu, ivml, lv, ispn ) & 
!                                          + hpwi( ii ) * mpci( ibnd, ikpt, nu, ivml, lv, ispn ) 
!              out_vec%i( ibnd, ikpt, ia ) = out_vec%i( ibnd, ikpt, ia )  & 
!                                          + hpwi( ii ) * mpcr( ibnd, ikpt, nu, ivml, lv, ispn ) & 
!                                          - hpwr( ii ) * mpci( ibnd, ikpt, nu, ivml, lv, ispn ) 
!            enddo
!
!        enddo
!      enddo


    enddo
    !

  end subroutine OCEAN_fg_slice


  subroutine OCEAN_ct_slice( sys, out_vec, inter, jb, jk, ja )
    use OCEAN_system
    use OCEAN_psi
    implicit none
    !
    type(O_system), intent( in ) :: sys
    type( ocean_vector ), intent( inout ) :: out_vec
    real(DP), intent( in ) :: inter
    integer, intent( in ) :: jb, jk, ja
    !
    !
    integer :: l, m, nu, ispn, ikpt, ibnd
    real( DP ) :: mul
    real( DP ), dimension( npmax ) :: ampr, ampi, hampr, hampi
  !

!   CT is diagonal in alpha

    mul = inter / ( dble( sys%nkpts ) * sys%celvol )

    ispn = 2 - mod( ja, 2 )
    if( sys%nspn .eq. 1 ) then
      ispn = 1
    endif
!JTV clean this up
    do l = lmin, lmax
      do m = -l, l
!        ampr( : ) = 0
!        ampi( : ) = 0
        hampr( : ) = 0
        hampi( : ) = 0
!             do i = 1, n
        do nu = 1, nproj( l )
          ampr( nu ) = mpcr( jb, jk, nu, m, l, ispn )
          ampi( nu ) = mpci( jb, jk, nu, m, l, ispn )
        enddo

        do nu = 1, nproj( l )
          hampr( : ) = hampr( : ) - mpm( :, nu, l ) * ampr( nu )
          hampi( : ) = hampi( : ) - mpm( :, nu, l ) * ampi( nu )
        enddo
        hampr( : ) = hampr( : ) * mul
        hampi( : ) = hampi( : ) * mul

        do nu = 1, nproj( l )
          ! only need upper triangle and ialpha = jalpha
          do ikpt = 1, sys%nkpts !jk-1
            do ibnd = 1, sys%num_bands
              out_vec%r( ibnd, ikpt, ja ) = out_vec%r( ibnd, ikpt, ja ) &
                                          + mpcr( ibnd, ikpt, nu, m, l, ispn ) * hampr( nu )  &
                                          + mpci( ibnd, ikpt, nu, m, l, ispn ) * hampi( nu )
              out_vec%i( ibnd, ikpt, ja ) = out_vec%i( ibnd, ikpt, ja ) &
                                          + mpcr( ibnd, ikpt, nu, m, l, ispn ) * hampi( nu ) &
                                          - mpci( ibnd, ikpt, nu, m, l, ispn ) * hampr( nu )
            enddo
          enddo
!          ikpt = jk
!          do ibnd = 1, jb
!            out_vec%r( ibnd, ikpt, ja ) = out_vec%r( ibnd, ikpt, ja ) &
!                                        + mpcr( ibnd, ikpt, nu, m, l, ispn ) * hampr( nu )  & 
!                                        + mpci( ibnd, ikpt, nu, m, l, ispn ) * hampi( nu ) 
!            out_vec%i( ibnd, ikpt, ja ) = out_vec%i( ibnd, ikpt, ja ) &
!                                        + mpcr( ibnd, ikpt, nu, m, l, ispn ) * hampi( nu ) &
!                                        - mpci( ibnd, ikpt, nu, m, l, ispn ) * hampr( nu )
!          enddo
        enddo
      enddo
    enddo


  end subroutine OCEAN_ct_slice



  subroutine OCEAN_ct_single( sys, bse_me, inter, ib, ik, ia, jb, jk, ja )
    use OCEAN_system
    implicit none
    !
    type(O_system), intent( in ) :: sys
    real(DP), intent( in ) :: inter
    integer, intent( in ) :: ib, ik, ia, jb, jk, ja
    complex(DP), intent( inout ) :: bse_me
    !
    !
    integer :: l, m, nu, ispn
    real( DP ) :: mul, outr, outi
    real( DP ), dimension( npmax ) :: ampr, ampi, hampr, hampi
  !

!   CT is diagonal in alpha
    if( ia .ne. ja ) return

    outr = 0.0_DP
    outi = 0.0_DP

    mul = inter / ( dble( sys%nkpts ) * sys%celvol )

    ispn = 2 - mod( ia, 2 )
    if( sys%nspn .eq. 1 ) then
      ispn = 1
    endif
!JTV clean this up
    do l = lmin, lmax
      do m = -l, l
!        ampr( : ) = 0
!        ampi( : ) = 0
        hampr( : ) = 0
        hampi( : ) = 0
!             do i = 1, n
        do nu = 1, nproj( l )
          ampr( nu ) = mpcr( ib, ik, nu, m, l, ispn ) 
          ampi( nu ) = mpci( ib, ik, nu, m, l, ispn )
        enddo

        do nu = 1, nproj( l )
          hampr( : ) = hampr( : ) - mpm( :, nu, l ) * ampr( nu )
          hampi( : ) = hampi( : ) - mpm( :, nu, l ) * ampi( nu )
        enddo
        hampr( : ) = hampr( : ) * mul
        hampi( : ) = hampi( : ) * mul

        do nu = 1, nproj( l )
          outr = outr + mpcr( jb, jk, nu, m, l, ispn ) * hampr( nu ) &
                      + mpci( jb, jk, nu, m, l, ispn ) * hampi( nu )
          outi = outi + mpcr( jb, jk, nu, m, l, ispn ) * hampi( nu ) &
                      - mpci( jb, jk, nu, m, l, ispn ) * hampr( nu )

        enddo
      enddo
    enddo

    bse_me = bse_me + CMPLX( outr, -outi, DP )

  end subroutine OCEAN_ct_single

  subroutine OCEAN_fg_single( sys, bse_me, inter, ib, ik, ia, jb, jk, ja )
    use OCEAN_system
    implicit none
    !
    type(O_system), intent( in ) :: sys
    real(DP), intent( in ) :: inter
    integer, intent( in ) :: ib, ik, ia, jb, jk, ja
    complex(DP), intent( inout ) :: bse_me
    !
    !
    integer :: lv, ii, ivml, nu, j1, jj, ispn
    real( DP ) :: mul, outr, outi
    real( DP ), allocatable, dimension( : ) :: pwr, pwi, hpwr, hpwi
    !
    allocate( pwr(itot), pwi( itot), hpwr(itot), hpwi(itot) )
    ! lvl, lvh is not enough for omp
    ! should probably pull thread spawning out of the do loop though

!JTV clean this up too
    outr = 0.0_DP
    outi = 0.0_DP
    mul = inter / ( dble( sys%nkpts ) * sys%celvol )

    do lv = lvl, lvh
      pwr( : ) = 0.0_DP
      pwi( : ) = 0.0_DP
      hpwr( : ) = 0.0_DP
      hpwi( : ) = 0.0_DP

      ispn = 2 - mod( ia, 2 )
      if( sys%nspn .eq. 1 ) then
        ispn = 1
      endif
      
      do ivml = -lv, lv
        do nu = 1, nproj( lv )
          ii = nu + ( ivml + lv ) * nproj( lv ) + ( ia - 1 ) * ( 2 * lv + 1 ) * nproj( lv )

          pwr( ii ) = pwr( ii ) + mpcr( ib, ik, nu, ivml, lv, ispn )
          pwi( ii ) = pwi( ii ) + mpci( ib, ik, nu, ivml, lv, ispn )

         enddo
       enddo

      do ii = 1, mham( lv )
        do jj = 1, mham( lv )
          j1 = jbeg( lv ) + ( jj -1 ) + ( ii - 1 ) * mham( lv )
          hpwr( ii ) = hpwr( ii ) + mhr( j1 ) * pwr( jj ) - mhi( j1 ) * pwi( jj )
          hpwi( ii ) = hpwi( ii ) + mhr( j1 ) * pwi( jj ) + mhi( j1 ) * pwr( jj )
        end do
        hpwr( ii ) = hpwr( ii ) * mul
        hpwi( ii ) = hpwi( ii ) * mul
      end do

      ispn = 2 - mod( ja, 2 )
      if( sys%nspn .eq. 1 ) then
        ispn = 1
      endif
      do ivml = -lv, lv
        do nu = 1, nproj( lv )
          ii = nu + ( ivml + lv ) * nproj( lv ) + ( ja - 1 ) * ( 2 * lv + 1 ) * nproj( lv )
          outr = outr + hpwr( ii ) * mpcr( jb, jk, nu, ivml, lv, ispn ) &
                      + hpwi( ii ) * mpci( jb, jk, nu, ivml, lv, ispn )
          outi = outi + hpwi( ii ) * mpcr( jb, jk, nu, ivml, lv, ispn ) &
                      - hpwr( ii ) * mpci( jb, jk, nu, ivml, lv, ispn )
        enddo
      enddo

    enddo
    !
    bse_me = bse_me + CMPLX( outr,-outi, DP ) 

  end subroutine OCEAN_fg_single


subroutine limel( l, m, mp, vrslt, nsphpt, xsph, ysph, zsph, wsph, prefs )
  implicit none
  !
  integer :: nsphpt, l, m, mp
  real( kind = kind( 1.0d0 ) ) :: prefs( 0 : 1000 )
  real( kind = kind( 1.0d0 ) ), dimension( nsphpt ) :: xsph, ysph, zsph, wsph
  complex( kind = kind( 1.0d0 ) ) :: vrslt( 3 )
  !
  integer :: i, j, mrun
  real( kind = kind( 1.0d0 ) ) :: xx, yy, zz, xarg, yarg, zarg
  complex( kind = kind( 1.0d0 ) ) :: su1( -l : l, 3 ), su2( -l : l, 3 ), ylm, ylmp, yrun
  !
  su1( :, : ) = 0.0d0
  su2( :, : ) = 0.0d0
  do j = 1, nsphpt
     xx = xsph( j ); yy = ysph( j ); zz = zsph( j )
     call newgetylm( l, m, xx, yy, zz, ylm, prefs )
     call newgetylm( l, mp, xx, yy, zz, ylmp, prefs )
     do i = 1, 3
        select case( i )
        case( 1 )
           zarg = xx; xarg = yy; yarg = zz
        case( 2 )
           zarg = yy; xarg = zz; yarg = xx
        case( 3 )
           zarg = zz; xarg = xx; yarg = yy
        end select
        do mrun = -l, l
           call newgetylm( l, mrun, xarg, yarg, zarg, yrun, prefs )
           su1( mrun, i ) = su1( mrun, i ) + wsph( j ) * conjg( yrun ) * ylm
           su2( mrun, i ) = su2( mrun, i ) + wsph( j ) * conjg( yrun ) * ylmp
        end do
     end do
  end do
  vrslt( : ) = 0.0d0
  do i = 1, 3
     do mrun = -l, l
        vrslt( i ) = vrslt( i ) + mrun * conjg( su1( mrun, i ) ) * su2( mrun, i )
     end do
  end do
  !
  return
end subroutine limel



end module OCEAN_multiplet
program buildAtomCenter
  use ai_kinds
  use periodic, only : get_atom_number
  use screen_opf, only : screen_opf_loadAll, screen_opf_getNrad, screen_opf_getRad
  implicit none


  integer :: natoms
  character(len=2), allocatable :: elnames(:)
  integer, allocatable :: elnumbers(:), zeelist(:)
  real(DP) :: kappa(3), avecs(3,3), xbohr(3), coord(3)
  real(DP), allocatable :: xred(:,:), angular(:,:)
  complex(DP), allocatable ::  full_grid(:,:)
  
  character(len=10) :: spcpnt

  integer : i, n, iatom, zee, j, k, nUniqueZ, ierr, nang, nr

  ierr = 0

  ! Load list of atoms and locations
  open( unit=99, file='xyz.wyck', form='formatted', status='old')
  read(99,*) natoms
  allocate(elnames(natoms), xred(3,natoms), elnumbers(natoms), zeelist(natoms) )
  do i = 1, natoms
    read( elnames(i), xred(:,i) )
    get_Atom_number( elnumbers(i), elnames(i) )
  enddo
  close(99)

  open( unit=99, file='avecsinbohr.ipt', form='formatted', status='old')
  read(99,*) avecs(:,:)
  close(99)

  nUniqueZ = 0
  do i = 1, natoms 
    do j = 1, nUniqueZ
      if( elnumbers(i) .eq. zeelist(j) ) goto 11
    enddo
    nUniqueZ = nUniqueZ + 1
    zeelist(nUniqueZ) = elnumbers(i)
11  continue
  enddo

  open( unit=99, file='ang.inp', form='formatted', status='old' )
  read( 99, * ) i
  close(99)
  write( spcpnt, '(A,I0') 'specpnt.', i
  open(unit=99, file=spcpnt, form='formatted', status='old' )
  read( 99, * ) nang
  allocate( angular(4,nang) )
  do i= 1, nang
    read(99,*) angular(:,i)
  enddo
  close(99)

  ! direction and magnitude of kappa, out-going electron
  open( unit=99, file='kappa', form='formatted', status='old' )
  read( 99, * ) kappa(:)
  close( 99 )

  call screen_opf_loadAll( ierr, zeelist(1:nUniqueZ) )

  !TODO, only deallocate/re-allocate if grid different sizes
  itarg = 1
  do iatom = 1, natoms
    call screen_opf_getNrad( elnumbers(iatom), nr, ierr, itarg )
    if( ierr .ne. 0 ) goto 111
    allocate( rad( nr ) )

    call screen_opf_getRad( elnumbers(iatom), rad, ierr, itarg )
    if( ierr .ne. 0 ) goto 111

    ! allocate angular x radial grid
    allocate( full_grid( nang, nr ) )

    ! step through and fill in X = R + Tau
    ! Then for each can calculate exp i Kappa (r+Tau)
    xbohr(:) = matmul( avecs, xred(:,iatom) )
    !
    do ir = 1, nr
      do iang = 1, nang
        coord( : ) = xbohr(:) + rad(ir) * angular(1:3, iang )
        full_grid(iang,ir) = exp( cmplx( 0.0_DP, dot_product( kappa, coord ), DP ) )
      enddo
    enddo

    
    ! Now for each n, l, m of the projectors
    
    

  enddo

111 continue

end program

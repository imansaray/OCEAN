program buildAtomCenter
  use ai_kinds
  use periodic, only : get_atom_number
  use screen_opf, only : screen_opf_loadAll, screen_opf_getNrad, screen_opf_getRad
  implicit none


  integer :: natoms
  character(len=2), allocatable :: elnames(:)
  integer, allocatable :: elnumbers(:), zeelist(:)
  real(DP) :: kappa(4)
  real(DP), allocatable :: xred(:,:), angular(:,:)
  
  character(len=10) :: spcpnt

  integer : i, n, iatom, zee, j, k, nUniqueZ, ierr

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
  read( 99, * ) n
  allocate( angular(4,n) )
  do i= 1, n
    read(99,*) angular(:,i)
  enddo
  close(99)

  ! direction and magnitude of kappa, out-going electron
  open( unit=99, file='kappa', form='formatted', status='old' )
  read( 99, * ) kappa(:)
  close( 99 )

  call screen_opf_loadAll( ierr, zeelist(1:nUniqueZ) )

  itarg = 1
  do iatom = 1, natoms
    call screen_opf_getNrad( elnumbers(iatom), nr, ierr, itarg )
    if( ierr .ne. 0 ) goto 111
    allocate( rad( nr ) )

    call screen_opf_getRad( elnumbers(iatom), rad, ierr, itarg )
    if( ierr .ne. 0 ) goto 111

    ! allocate angular x radial grid
    ! step through and fill in X = R + Tau
    ! Then for each can calculate exp i Kappa (r+Tau)
    ! just calculate directly, 

  enddo

111 continue

end program

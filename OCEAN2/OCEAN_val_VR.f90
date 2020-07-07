! Copyright (C) 2020 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!


module OCEAN_val_VR
  use AI_kinds, only : DP


  implicit none
  save
  private

  logical :: is_init = .false.

  real(DP), allocatable :: VR(:) !> A contact interation in real-space

  public :: OCEAN_val_VR_act, OCEAN_val_VR_init, OCEAN_val_VR_clean

  contains

  subroutine OCEAN_val_VR_init( sys, ierr )
    use OCEAN_system, only : o_system
    use OCEAN_val_states, only : nxpts
    use OCEAN_kxc
    !
    type(o_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    !

    allocate( VR(nxpts), STAT=ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_kxc( sys, VR, ierr )
    if( ierr .ne. 0 ) return

    is_init = .true.
  end subroutine OCEAN_val_VR_init


  subroutine OCEAN_val_VR_clean()
    if( allocated( VR ) ) deallocate( VR )
    is_init = .falst.
  end subroutine OCEAN_val_VR_clean


  subroutine OCEAN_val_VR_act( sys, psi, psi_out, ierr )
    use OCEAN_psi, only : OCEAN_vector
    use OCEAN_system, only : o_system
    implicit none
    !
    type(o_system), intent( in ) :: sys
    type(OCEAN_vector), intent( in ) :: psi
    type(OCEAN_vector), intent( inout ) :: psi_out
    integer, intent( inout ) :: ierr
    !
    integer :: i, j, ibeta, cspn, vspn

    ibeta = 0
    do i = 1, sys%valence_ham_spin
      vspn = min( i, sys%nspn )
      do j = 1, sys%valence_ham_spin
        cspn = min( j, sys%nspn )
        ibeta = ibeta + 1

        call OCEAN_ladder_act_single( sys, psi, psi_out, ibeta, cspn, vspn, ierr )
        if( ierr .ne. 0 ) return

      enddo
    enddo

  end subroutine OCEAN_ladder_act

  subroutine act_single( sys, psi, psi_out, psi_spn, cspn, vspn, ierr )
    use OCEAN_psi
    use OCEAN_mpi
    use OCEAN_val_states, only : nxpts, nxpts_pad, re_val, im_val, re_con, im_con, &
                                 nbv, nbc, max_nxpts, nxpts_by_mpiID, startx_by_mpiID
    use OCEAN_system

    type(o_system), intent( in ) :: sys
    type(OCEAN_vector), intent( in ) :: psi
    type(OCEAN_vector), intent( inout ) :: psi_out
    integer, intent( in ) :: psi_spn, cspn, vspn
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: re_a_mat(:,:), im_a_mat(:,:), re_b_mat(:), im_b_mat(:)
    real(DP) :: spin_prefac, minus_spin_prefac
    integer :: ik, ib, psi_con_pad


    spin_prefac = 1.0_DP
    minus_spin_prefac = - spin_prefac
    
    call OCEAN_psi_returnBandPad( psi_con_pad, ierr )
    if( ierr .ne. 0 ) return


    allocate( re_a_mat( nxpts, nbv ), im_a_mat( nxpts, nbv ) )
    allocate( re_b_mat( nxpts ), im_b_mat( nxpts ) )
    re_b_mat(:) = 0.0_DP
    im_b_mat(:) = 0.0_DP

    do ik = 1, nkpts
      call DGEMM( 'N', 'N', nxpts, nbv, nbc, one, re_con( 1, 1, ik, cspn ), nxpts_pad, &
                  psi%valr( 1, ib, ik, psi_spn ), psi_con_pad, zero, re_a_mat, nxpts )
      call DGEMM( 'N', 'N', nxpts, nbv, nbc, minusone, im_con( 1, 1, ik, cspn ), nxpts_pad, &
                  psi%vali( 1, ib, ik, psi_spn ), psi_con_pad, one, re_a_mat, nxpts )

      call DGEMM( 'N', 'N', nxpts, nbv, nbc, one, re_con( 1, 1, ik, cspn ), nxpts_pad, &
                        psi%vali( 1, ib, ik, psi_spn ), psi_con_pad, zero, im_a_mat, nxpts )
      call DGEMM( 'N', 'N', nxpts, nbv, nbc, one, im_con( 1, 1, ik, cspn ), nxpts_pad, &
                        psi%valr( 1, ib, ik, psi_spn ), psi_con_pad, one, im_a_mat, nxpts )

      do ib = 1, nbv
        ! valence bands are starred
        re_b_mat( : ) = re_b_mat( : ) + re_a_mat( :, ib ) * re_val( 1:nxpts, ib, ik, vspn ) &
                      + im_a_mat( :, ib ) * im_val( :, ib, ik, vspn )

        im_b_mat( : ) = im_b_mat( : ) + im_a_mat( :, ib ) * re_val( 1:nxpts, ib, ik, vspn ) &
                      - re_a_mat( :, ib ) * im_val( :, ib, ik, vspn )

    enddo
      

    re_b_mat(:) = re_b_mat(:) * VR(:)
    im_b_mat(:) = im_b_mat(:) * VR(:)


    do ik = 1, nkpt
      do ib = 1, nbv
        re_a_mat( :, ib ) = re_b_mat( : ) * re_val( 1:nxpts, ib, ik, vspn ) &
                          - im_b_mat( : ) * im_val( 1:nxpts, ib, ik, vspn )
        im_a_mat( :, ib ) = re_b_mat( : ) * im_val( 1:nxpts, ib, ik, vspn ) &
                          + im_b_mat( : ) * re_val( 1:nxpts, ib, ik, vspn )
      enddo

      ! conduction bands are starred
      call DGEMM( 'T', 'N', nbc, nbv, nxpts, spin_prefac, re_con( 1, 1, ik, cspn ), nxpts_pad, &
                  re_a_mat, nxpts, one, psi_out%valr( 1, 1, ik, psi_spn ), psi_con_pad )
      call DGEMM( 'T', 'N', nbc, nbv, nxpts, spin_prefac, im_con( 1, 1, ik, cspn ), nxpts_pad, &
                  im_a_mat, nxpts, one, psi_out%valr( 1, 1, ik, psi_spn ), psi_con_pad )

      call DGEMM( 'T', 'N', nbc, nbv, nxpts, spin_prefac, re_con( 1, 1, ik, cspn ), nxpts_pad, &
                        im_a_mat, nxpts, one, psi_out%vali( 1, 1, ik, psi_spn ), psi_con_pad )
      call DGEMM( 'T', 'N', nbc, nbv, nxpts, minus_spin_prefac, im_con( 1, 1, ik, cspn ), nxpts_pad, &
                        re_a_mat, nxpts, one, psi_out%vali( 1, 1, ik, psi_spn ), psi_con_pad )

    enddo

    deallocate( re_b_mat, im_b_mat, re_a_mat, im_a_mat )
                  
  end subroutine act_single

end module OCEAN_val_VR


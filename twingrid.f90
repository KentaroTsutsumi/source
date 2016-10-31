!
! Copyright (C) 2001-2014 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------
!
! This program set q points grid in irr. BZ with a harf shift and
! twin k grid in whole BZ (with and without a harf shift).
!
!------------------------------------------------------------------------------
module twingrid_vals
  !----------------------------------------------------------------------------
  !
  ! This module contains global variables for twingrid.x
  !
  implicit none
  !
  integer,save :: &
  & nsym, & ! # of symmetries
  & nk,   & ! Total # of k
  & ng(3)   ! Grid for k-points ng(1) * ng(2) * ng(3)
  !
  integer,allocatable,save :: &
  & sym(:,:,:)   ! (3,3,nsym) Symmetry operators
  !
  real(8),save :: &
  & bvec(3,3) ! Reciplocal lattice vector [alat]
  !
end module twingrid_vals
!
!-------------------------------------------------------------------------------
module twingrid_routines
  !-----------------------------------------------------------------------------
  !
  ! This module contains subroutines for twingrid.x
  !
  implicit none
  !
contains
!
!------------------------------------------------------------------
subroutine read_dat()
  !----------------------------------------------------------------
  !
  ! This routine reads symmetry and bvec from data-file.xml
  !  ng(1:3) from standard input
  !
  use twingrid_vals, only : nsym, sym, bvec, nk, ng
  use iotk_module
  !
  integer :: isym, fi = 10
  logical :: invsym
  character(iotk_namlenx) :: attr, csym
  !
  ! Open datafile.xml
  !
  write(*,'(a,a)') "   open data-file.xml"
  call iotk_open_read(fi,"data-file.xml")
  !
  call iotk_scan_begin(fi,"CELL")
  !
  ! Read reciprocal lattice vecor
  !
  call iotk_scan_begin(fi,"RECIPROCAL_LATTICE_VECTORS")
  call iotk_scan_dat(fi,"b1",bvec(1:3,1))
  call iotk_scan_dat(fi,"b2",bvec(1:3,2))
  call iotk_scan_dat(fi,"b3",bvec(1:3,3))
  call iotk_scan_end(fi,"RECIPROCAL_LATTICE_VECTORS")
  !
  write(*,*) "  Reciprocal lattice vector[2pi/a] : "
  write(*,'(3f10.5)') bvec(1:3,1:3)
  !
  call iotk_scan_end(fi,"CELL")
  !
  ! Read # of point symmetry
  !
  call iotk_scan_begin(fi,"SYMMETRIES")
  call iotk_scan_dat(fi,"NUMBER_OF_SYMMETRIES",nsym)
  call iotk_scan_dat(fi,"INVERSION_SYMMETRY",invsym)
  if(invsym) then
     allocate(sym(3,3,nsym))
     write(*,*) "  # of BZ symmetry : ", nsym
  else
     allocate(sym(3,3,nsym * 2))
     write(*,*) "  Inversion symmetry is added."
     write(*,*) "  # of BZ symmetry : ", nsym * 2
  end if
  !
  ! Read symmmetry operators
  !
  do isym = 1, nsym
     !
     write(csym,*) isym
     write(attr,'(a,a)') "SYMM.", trim(adjustl(csym))
     !
     call iotk_scan_begin(fi,trim(attr))
     call iotk_scan_dat(fi,"ROTATION",sym(1:3,1:3,isym))
     call iotk_scan_end(fi,trim(attr))
     !
  end do
  !
  if(.not. invsym) then
     sym(1:3,1:3,nsym + 1:nsym + nsym) = - sym(1:3,1:3,1:nsym)
     nsym = nsym * 2
  end if
  !
  call iotk_scan_end(fi,"SYMMETRIES")
  !
  call iotk_close_read(fi)
  !
  write(*,'(a)', advance = "no") " q grid ? "
  read(*,*) ng(1:3)
  !
  nk = product(ng(1:3))
  !
end subroutine read_dat
!
!-----------------------------------------------------------------------------
subroutine irr_bz()
  !---------------------------------------------------------------------------
  !
  ! This routine compute the irreducible BZ and output to a file (q.out).
  !
  use twingrid_vals, only : nk, ng, nsym, sym, bvec
  !
  integer :: i1, i2, i3, ik, isym, nk0
  real(8) :: kv0(3,nk), kv1(3), kv2(3)
  !
  nk0 = 0
  do i1 = 1, ng(1)
     do i2 = 1, ng(2)
        do i3 = 1, ng(3)
           !
           if(ng(3) /= 1) then
             kv1(1:3) = (dble((/i1, i2, i3/)) - 0.5d0) / dble(ng(1:3))
           else
             kv1(1) = (dble(i1) - 0.5d0) / dble(ng(1))
             kv1(2) = (dble(i2) - 0.5d0) / dble(ng(2))
             kv1(3) = (dble(i3) - 1.0d0) / dble(ng(3))  ! modified to calc. in case of FeSe
           end if
           !
           kv1(1:3) = kv1(1:3) - dble(floor(kv1(1:3) + 0.5d0 + 1d-4))
           !
           do isym = 1, nsym
              !
              kv2(1:3) = matmul(dble(sym(1:3,1:3,isym)), kv1(1:3))
              kv2(1:3) = kv2(1:3) - dble(floor(kv2(1:3) + 0.5d0 + 1d-4))
              !
              do ik = 1, nk0
                 if(all(abs(kv2(1:3) - kv0(1:3,ik)) < 1d-8)) goto 10
              end do
              !
           end do
           !
           nk0 = nk0 + 1
           kv0(1:3,nk0) = kv1(1:3)
           !
10         continue
           !
        end do
     end do
  end do
  !
  write(*,*) "  # of q : ", nk0
  !
  ! Write input q for ph.x
  !
  open(20, file = "q.out")
  !
  do ik = 1, nk0
     kv1(1:3) = matmul(bvec(1:3,1:3), kv0(1:3,ik))
     write(20,'(i5,10f20.15)') ik, kv1(1:3)
  end do
  !
  close(20)
  !
end subroutine irr_bz
!
!------------------------------------------------------------------------------------
subroutine write_kvec()
  !----------------------------------------------------------------------------------
  !
  ! This routine outputs the twin k grid in whole BZ (with and without a harf shift)
  !
  use twingrid_vals, only : nk, bvec, ng
  !
  integer :: i1, i2, i3
  real(8) :: kv1(3)
  !
  ! Write input k for pw.x
  !
  open(20, file = "k.out")
  !
  write(20,*) nk * 2
  !
  do i3 = 1, ng(3)
     do i2 = 1, ng(2)
        do i1 = 1, ng(1)
           !
           kv1(1:3) = (dble((/i1, i2, i3/)) - 1d0) / dble(ng(1:3))
           kv1(1:3) = kv1(1:3) - dble(floor(kv1(1:3) + 0.5d0 + 1d-4))
           kv1(1:3) = matmul(bvec(1:3,1:3), kv1(1:3))
           write(20,'(10f20.15)') kv1(1:3), 1d0
           !
        end do
     end do
  end do
  !
  do i3 = 1, ng(3)
     do i2 = 1, ng(2)
        do i1 = 1, ng(1)
           !
           if(ng(3) /= 1) then
             kv1(1:3) = (dble((/i1, i2, i3/)) - 0.5d0) / dble(ng(1:3))
           else
             kv1(1) = (dble(i1) - 0.5d0) / dble(ng(1))
             kv1(2) = (dble(i2) - 0.5d0) / dble(ng(2))
             kv1(3) = (dble(i3) - 1.0d0) / dble(ng(3))  ! modified to calc. in case of FeSe
           end if
           !
           kv1(1:3) = kv1(1:3) - dble(floor(kv1(1:3) + 0.5d0 + 1d-4))
           kv1(1:3) = matmul(bvec(1:3,1:3), kv1(1:3))
           write(20,'(10f20.15)') kv1(1:3), 1d0
           !
        end do
     end do
  end do
  !
  close(20)
  !
end subroutine write_kvec
!
end module twingrid_routines
!
!--------------------------------------------------------------------------------
program twingrid
  !------------------------------------------------------------------------------
  !
  ! This program set q points grid in irr. BZ with a harf shift and
  ! twin k grid in whole BZ (with and without a harf shift).
  !
  use twingrid_routines, only : read_dat, irr_bz, write_kvec 
  implicit none
  !
  call read_dat()
  !
  call irr_bz()
  !
  call write_kvec()
  !
end program twingrid

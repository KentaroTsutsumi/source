module get_eigval_vals
  !
  implicit none
  !
  integer,save :: &
  & nk, &
  & nb, &
  & nsym, &
  ng(3)
  !
  integer,allocatable,save :: &
  & sym(:,:,:)
  !
  real(8),allocatable,save :: &
  & kvec(:,:), &
  & eig(:,:), &
  & eig2(:,:,:,:)
  !
end module get_eigval_vals
!
!
!
module get_eigval_routines
  !
  implicit none
  !
contains
!
! Read dat
!
subroutine read_dat()
  !
  use get_eigval_vals, only : nk, nb, nsym, sym, kvec, eig, ng
  use iotk_module
  !
  integer :: isym, fi = 10, ik
  real(8) :: a0, ef, avec(3,3), bvec(3,3)
  logical :: invsym
  character(iotk_namlenx) :: attr, csym
  !
  ! Open datafile.xml
  !
  write(*,'(a,a)') "   open data-file.xml"
  call iotk_open_read(fi,"data-file.xml")
  !
  ! Read CELL PARAMETER
  !
  call iotk_scan_begin(fi,"CELL")
  call iotk_scan_dat(fi,"LATTICE_PARAMETER",a0)
  write(*,*) "  Lattice parameter[a.u.] : ", a0
  !
  ! Read Lattice vector
  !
  call iotk_scan_begin(fi,"DIRECT_LATTICE_VECTORS")
  call iotk_scan_dat(fi,"a1",avec(:,1))
  call iotk_scan_dat(fi,"a2",avec(:,2))
  call iotk_scan_dat(fi,"a3",avec(:,3))
  avec(:,:) = avec(:,:) / a0
  call iotk_scan_end(fi,"DIRECT_LATTICE_VECTORS")
  write(*,*) "  Direct lattice vector[a] : "
  write(*,'(3f10.5)') avec(:,:)
  !
  ! Read reciprocal lattice vecor
  !
  call iotk_scan_begin(fi,"RECIPROCAL_LATTICE_VECTORS")
  call iotk_scan_dat(fi,"b1",bvec(:,1))
  call iotk_scan_dat(fi,"b2",bvec(:,2))
  call iotk_scan_dat(fi,"b3",bvec(:,3))
  call iotk_scan_end(fi,"RECIPROCAL_LATTICE_VECTORS")
  !
  write(*,*) "  Reciprocal lattice vector[2pi/a] : "
  write(*,'(3f10.5)') bvec(:,:)
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
     call iotk_scan_dat(fi,"ROTATION",sym(:,:,isym))
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
  ! Read Monkhorst-Pack grid
  !
  call iotk_scan_begin(fi,"BRILLOUIN_ZONE")
  !
  attr=""
  call iotk_scan_empty(fi,"MONKHORST_PACK_GRID",attr)
  call iotk_scan_attr(attr,"nk1",ng(1))
  call iotk_scan_attr(attr,"nk2",ng(2))
  call iotk_scan_attr(attr,"nk3",ng(3))
  !
  write(*,*) "  MONKHORST_PACK_GRID : ", ng(:)
  call iotk_scan_end(fi,"BRILLOUIN_ZONE")
  !
  ! Read # of k, # of band , Ef
  !
  call iotk_scan_begin(fi,"BAND_STRUCTURE_INFO")
  !
  call iotk_scan_dat(fi,"NUMBER_OF_K-POINTS",nk)
  write(*,*) "  # of k-vector : ", nk
  !
  call iotk_scan_dat(fi,"NUMBER_OF_BANDS",nb)
  write(*,*) "  # of band : ", nb
  !
  call iotk_scan_dat(fi,"FERMI_ENERGY",Ef)
  write(*,*) "  Fermi energy[Ryd] : ", Ef * 2
  !
  call iotk_scan_end(fi,"BAND_STRUCTURE_INFO")
  !
  allocate(kvec(3,nk),eig(nb,nk))
  !
  ! Read k-vector
  !
  call iotk_scan_begin(fi,"EIGENVALUES")
  !
  do ik = 1, nk
     !
     write(csym,*) ik
     write(attr,'(a,a)') "K-POINT.", trim(adjustl(csym))
     !
     call iotk_scan_begin(fi,trim(attr))
     call iotk_scan_dat(fi,"K-POINT_COORDS",kvec(1:3,ik))
     call iotk_scan_end(fi,trim(attr))
     !
     kvec(1:3,ik) = matmul(kvec(1:3,ik), avec(1:3,1:3))
     !
  end do
  !
  call iotk_scan_end(fi,"EIGENVALUES")
  !
  call iotk_close_read(fi)
  !
  ! Read eigenvalue
  !
  do ik = 1, nk
     !
     write(attr,'(a,i5.5,a)') "K", ik, "/eigenval.xml"
     call iotk_open_read(fi, trim(attr))
     call iotk_scan_dat(fi,"EIGENVALUES",eig(1:nb,ik) )
     call iotk_close_read(fi)
     !
  enddo
  !
  ! Htr -> Ryd . Mesured from Fermi energy.
  !
  eig(:,:) = ( eig(:,:) - Ef ) * 2d0
  !
end subroutine read_dat
!
!  Write Symm opt
!
subroutine write_symm()
  !
  use get_eigval_vals, only : nsym, sym
  !
  integer :: fo = 10
  !
  write(*,*) "    open symm.dat"
  open(fo, file = 'symm.dat')
  !
  write(fo,*) nsym
  write(fo,*) sym(1:3,1:3,1:nsym)
  !
  close(fo)
  !
end subroutine write_symm
!
! 
!
subroutine rotate_k()
  !
  use get_eigval_vals, only : nk, nb, nsym, ng, sym, kvec, eig, eig2
  !
  integer :: isym, ik, ikv(3)
  real(8) :: kv(3)
  logical :: ldone(ng(1),ng(2),ng(3))
  !
  ldone(:,:,:) = .false.
  !
  do ik = 1, nk
     !
     do isym = 1, nsym
        !
        kv(1:3) = matmul(dble(sym(1:3,1:3,isym)),kvec(1:3,ik)) * dble(ng(1:3))
        ikv(1:3) = nint(kv(1:3))
        !
        if(any(abs(kv(1:3) - dble(ikv(1:3))) > 1d-8)) cycle
        !
        ikv(1:3) = modulo(ikv(1:3), ng(1:3)) + 1
        eig2(1:nb,ikv(1),ikv(2),ikv(3)) = eig(1:nb,ik)
        ldone(ikv(1),ikv(2),ikv(3)) = .true.
        !
     end do ! End isym
     !
  end do ! End ik
  !
  ! Check
  !
  if(count( .not. ldone) /= 0) &
  &     write(*,*)  "  # of elements that are not done : ", count( .not. ldone)
  !
end subroutine rotate_k
!
end module get_eigval_routines
!
!main module
!
program get_eigval
  !
  use get_eigval_vals, only : nb, ng, eig2
  use get_eigval_routines, only : read_dat, rotate_k, write_symm
  !
  implicit none
  !
  ! read symmmetry
  !
  write(*,*) ""
  write(*,*) "#####  Read data  #####"
  write(*,*) ""
  call read_dat()
  !
  allocate(eig2(nb, ng(1), ng(2), ng(3)) )
  !
  ! Rotate k-points
  !
  write(*,*) ""
  write(*,*) "##### Rotate & write k-poinnts #####"
  write(*,*) ""
  call rotate_k()
  !
  !  Write eigenvalue
  !
  open(11, file = 'eigval.dat')
  write(*,*) "  write to eigval.dat"
  !
  write(11,*) eig2(1:nb,1:ng(1),1:ng(2),1:ng(3))
  !
  close(11)
  !
  ! Write Symm. opt.
  !
  call write_symm()
  !
  write(*,*) ""
  write(*,*) "##### Done #####"
  write(*,*) ""
  !
end program get_eigval

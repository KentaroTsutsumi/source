module ultrasonic_vals
  !
  implicit none
  !
  real(8),parameter :: &
  & pi = acos(-1d0) ! Pi = 3.14159265358979...
  !
  integer,save :: &
  & ltetra,       & ! Type of tetra
  & ivvec(3,20,6),& ! points for tetrahedron method
  & ncf,          & ! # of points for Chebyshev interpolation
  & nbee,         & ! # of Bands for Vc
  & fbee,         & ! First band for Vc
  & lbee,         & ! Last band for Vc
  & nbep,         & ! # of bands for elph
  & fbep,         & ! First band for elph
  & lbep,         & ! Last band for elph
  & fstk,         & ! First k of this PE
  & lstk,         & ! Last k of this PE
  & my_rank,      & ! Index of PE
  & petot,        & ! # of PEs
  & nmf,          & ! # of Matsubara frequencoes
  & nt,           & ! Total # of Delta, Xi, ...
  & nb,           & ! # of bands
  & nbf,          & ! # of bands contais FS
  & fstb,         & ! first band contais FS
  & lstb,         & ! last band contais FS
  & ng(3),        & ! k grid for matrix elements
  & nk,           & ! ng(1) * ng(2) * ng(3). # of total k points for matrix elements
  & nk0,          & ! # of irreducible k points
  & ngd(3),       & ! k grid for DOS
  & nkd,          & ! ngd(1) * ngd(2) * ngd(3). # of total k points for DOS
  & nm,           & ! # of modes
  & nsym            ! # of symmetries
  !
  real(8),save :: &
  & avec(3,3), & ! Lattice vector
  & wlsm(4,20),& ! Weights for tetrahedron method
  & xic,       & !
  & bvec(3,3), & ! Reciplocal lattice vector
  & beta         ! inversed temperature [Ry]
  !
  integer,allocatable,save :: &
  & dgrid(:,:),  & ! (3,nkd)
  & kindx(:),    & ! (nt) k point for gap equation
  & bindx(:),    & ! (nt) band index for gap equation
  & sym(:,:,:)     ! (3,3,nsym) Symmetrical operators
  !
  real(8),allocatable,save :: &
  & rkv(:,:),       & ! (3,nk0) irreducible k points
  & grid(:,:),      & ! (3,nk) k points
  & mf(:),          & ! (nmf) Matsubara frequencies
  & wmf(:),         & ! (nmf) Weight for frequency integration
  & Z(:,:),         & ! (nbf,nkd) Renormalization fuctor
  & omg0(:,:),      & ! (nm,nk0) Phonon frequencies [Ry]
  & omgf(:,:,:),    & ! (nm,nk,nkd) Phonon frequencies [Ry]
  & gg0(:,:,:,:,:), & ! (nm,nb,nb,nk,nk0) El-Ph matrix element [Ry]
  & Vc0(:,:,:,:,:), & ! (0:nmf,nb,nb,nk,nk0) Screened Coulomb matrix element [Ry]
  & ggf(:,:,:,:,:), & ! (nm,nb,nk,nbf,nkd) El-Ph matrix element [Ry]
  & Vcf(:,:,:,:,:), & ! (0:nmf,nb,nk,nbf,nkd) Screened Coulomb matrix element [Ry]
  & delta(:),       & ! (nt) Kohn-Sham gap functions [Ry]
  & dltf(:,:),      & ! (nbf,nkd) Kohn-Sham gap functions [Ry]
  & xi(:),          & ! (nt) Kohn-Sham energy [Ry]
  & dk(:),          & ! (nt) Weight of k
  & vf(:,:,:),      & ! (nbf,nkd)
  & eig(:,:,:,:)      ! (nb,ngd(1),ngd(2),ngd(3)) Kohn-Sham energy [Ry]
  !
end module ultrasonic_vals
!
module ultrasonic_routines
  !
  implicit none
  !
contains
  !
 ! Standard input
!
subroutine read_stdin()
  !
  use mpi
  use ultrasonic_vals, only : my_rank, beta, xic, fbee, lbee, nbee, xic, nmf, ltetra
  !
  integer :: ierr
  !
  namelist /input/ beta, fbee, lbee, xic, nmf, ltetra
  !
  ltetra = 1
  xic = -1d0
  !
  if(my_rank == 0) then
     !
     read(*,input,err=100)
     write(*,*) '              Temparature[K] : ', beta
     beta = 157887d0 / beta
     write(*,*) '    Inverse temparature[/Ry] : ', beta
     write(*,*) '                   Xi cutoff : ', xic
     write(*,*) '                  First band : ', fbee
     write(*,*) '                   Last band : ', lbee
     write(*,*) "  # of Matsubara frequencies : ", nmf     
     write(*,*) '                      ltetra : ', ltetra
     !
  end if
  !
  call MPI_BCAST(beta, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(xic,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(fbee, 1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(lbee, 1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  nbee = lbee - fbee + 1
  call MPI_BCAST(nmf,  1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ltetra, 1, MPI_INTEGER,        0, MPI_COMM_WORLD, ierr)
  !
  return
  !
100 write(*,*) "Stop in read_stdin. reading namelist file"
  write(*,*) '              Temparature[K] : ', beta
  beta = 157887d0 / beta
  write(*,*) '    Inverse temparature[/Ry] : ', beta
  write(*,*) '                   Xi cutoff : ', xic
  write(*,*) '                  First band : ', fbee
  write(*,*) '                   Last band : ', lbee
  write(*,*) "  # of Matsubara frequencies : ", nmf     
  write(*,*) '                      ltetra : ', ltetra
  call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  call MPI_FINALIZE(ierr)
  stop
  !
end subroutine read_stdin
!
! Read data-file.xml
!
subroutine read_dat()
  !
  use mpi
  use iotk_module
  use ultrasonic_vals, only : my_rank, nsym, sym, nb, ngd, nkd, eig, avec, bvec, nbf, &
  &                  fstb, lstb, fstk, lstk
  !use ultrasonic_routines, only : cnt_and_dsp
  !
  integer :: isym, fi = 10, ib
  real(8) :: alat
  logical :: invsym
  character(iotk_namlenx) :: attr
  !
  if(my_rank == 0) then
     !  
     ! Open datafile.xml
     !
     write(*,'(a,a)') "   open data-file.xml"
     call iotk_open_read(fi,"data-file.xml")
     !
     ! Read reciprocal lattice vecor
     !
     call iotk_scan_begin(fi,"CELL")
     !
     call iotk_scan_dat(fi,"LATTICE_PARAMETER",alat)
     write(*,*) "  Lattice constant [a.u.] : ", alat
     !
     call iotk_scan_begin(fi,"DIRECT_LATTICE_VECTORS")
     call iotk_scan_dat(fi,"a1", avec(1:3,1))
     call iotk_scan_dat(fi,"a2", avec(1:3,2))
     call iotk_scan_dat(fi,"a3", avec(1:3,3))
     call iotk_scan_end(fi,"DIRECT_LATTICE_VECTORS")
     !avec(1:3,1:3) = avec(1:3,1:3) / alat
     write(*,*) "  Lattice vector[a] : "
     write(*,'(3f10.5)') avec(1:3,1:3) / alat
     !
     call iotk_scan_begin(fi,"RECIPROCAL_LATTICE_VECTORS")
     call iotk_scan_dat(fi,"b1",bvec(1:3,1))
     call iotk_scan_dat(fi,"b2",bvec(1:3,2))
     call iotk_scan_dat(fi,"b3",bvec(1:3,3))
     call iotk_scan_end(fi,"RECIPROCAL_LATTICE_VECTORS")
     call iotk_scan_end(fi,"CELL")
     !
     write(*,*) "  Reciprocal lattice vector[2pi/a] : "
     write(*,'(3f10.5)') bvec(1:3,1:3)
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
        write(attr,*) isym
        write(attr,'(a,a)') "SYMM.", trim(adjustl(attr))
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
     attr=""
     call iotk_scan_empty(fi,"MONKHORST_PACK_GRID",attr)
     call iotk_scan_attr(attr,"nk1",ngd(1))
     call iotk_scan_attr(attr,"nk2",ngd(2))
     call iotk_scan_attr(attr,"nk3",ngd(3))
     call iotk_scan_end(fi,"BRILLOUIN_ZONE")
     !
     nkd = product(ngd(1:3))
     !
     write(*,*) "   Dense grid : ", ngd(1:3)
     write(*,*) " # of Dense k : ", nkd
     !
     ! Read # of band
     !
     call iotk_scan_begin(fi,"BAND_STRUCTURE_INFO")
     call iotk_scan_dat(fi,"NUMBER_OF_BANDS",nb)
     call iotk_scan_end(fi,"BAND_STRUCTURE_INFO")
     !
     write(*,*) "  # of bands : ", nb
     !
     call iotk_close_read(fi)
     !
     ! Read eigenvalue
     !
     allocate(eig(nb, ngd(1), ngd(2), ngd(3)))
     open(fi, file = 'eigval.dat')
     write(*,*) "  read from eigval.dat"
     read(fi,*) eig(1:nb,1:ngd(1),1:ngd(2),1:ngd(3))
     close(fi)
     !
  end if ! if(my_rank == 0)
  !
  call MPI_BCAST(nsym, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, isym)
  if(my_rank /= 0) allocate(sym(3,3,nsym))
  call MPI_BCAST(sym,  9 * nsym, MPI_INTEGER, 0, MPI_COMM_WORLD, isym)
  call MPI_BCAST(nb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, isym)
  call MPI_BCAST(ngd, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, isym)
  call MPI_BCAST(nkd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, isym)
  call cnt_and_dsp(nkd,lstk,fstk)
  lstk = fstk + lstk
  fstk = fstk + 1
  if(my_rank /= 0) allocate(eig(nb, ngd(1), ngd(2), ngd(3)))
  call MPI_BCAST(eig, nb * nkd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, isym)
  call MPI_BCAST(avec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, isym)
  call MPI_BCAST(bvec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, isym)
  !
  ! Which band contains Fermi level ?
  !
  do ib = 1, nb
     !
     if(minval(eig(ib,1:ngd(1),1:ngd(2),1:ngd(3))) < 0d0 .and. &
     &  maxval(eig(ib,1:ngd(1),1:ngd(2),1:ngd(3))) > 0d0) then
        fstb = ib
        exit
     end if
     !
  end do
  !
  do ib = fstb, nb
     !
     if(minval(eig(ib,1:ngd(1),1:ngd(2),1:ngd(3))) < 0d0 .and. &
     &  maxval(eig(ib,1:ngd(1),1:ngd(2),1:ngd(3))) > 0d0) then
        lstb = ib
     end if
     !
  end do
  !
  nbf = lstb - fstb + 1
  !
  if(my_rank == 0) then
     write(*,*) "   Lowest band which contains FS : ", fstb
     write(*,*) "  Highest band which contains FS : ", lstb
     write(*,*) "    # of bands which contains FS : ", nbf
  end if
  !
end subroutine read_dat
!
! Read elph*.dat
!
subroutine read_elph()
  !
  use mpi
  use ultrasonic_vals, only : my_rank, ng, nk, nk0, nm, gg0, omg0, rkv, nbep, fbep, lbep
  !use ultrasonic_routines, only : cnt_and_dsp
  use omp_lib
  !
  integer :: nb0, ik, ierr, fi = 10, iqv(3), cnt, dsp
  real(8) :: qvec(3)
  character(100) :: fname
  !
  if(my_rank == 0) then
     !
     ! Read # of k, bands, modes
     !
     open(fi, file = "elph1.dat" )
     !
     read(fi,*) ng(1:3)
     write(*,*) "  k grid for MEs ", ng(1:3)
     !
     read(fi,*) nbep, fbep
     !
     read(fi,*) qvec(1:3)
     !
     read(fi,*) nm
     write(*,*) "  # of modes : ", nm
     close(fi)
     !
  end if
  !
  call MPI_BCAST(ng,   3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  nk = product(ng(1:3))
  call MPI_BCAST(nm,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nbep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(fbep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  lbep = fbep + nbep - 1
  !
  ! Search irreducible k points
  !
  call irr_bz()
  !
  call cnt_and_dsp(nk0,cnt,dsp)
  !
  allocate(gg0(nm,fbep:lbep,fbep:lbep,nk,cnt), omg0(nm,cnt))
  !
  ! Read omega & gg
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(cnt,dsp,nk0,nm,nbep,fbep,lbep,nk,ng,rkv,omg0,gg0) &
  !$OMP & PRIVATE(fi,ik,fname,nb0,qvec,iqv,ierr)
  !
  fi = OMP_GET_THREAD_NUM() + 10
  !
  !$OMP DO
  do ik = 1, cnt
     !
     write(fname,*) dsp + ik
     write(fname,'(3a)') "elph", trim(adjustl(fname)), ".dat"
     open(fi, file = trim(fname))
     !
     read(fi,*) iqv(1:3)
     if(.not. all(iqv(1:3) == ng(1:3))) then
        write(*,*) "Stop in read_elph. k grid. ik = ", dsp + ik
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        call MPI_FINALIZE(ierr)
        stop
     end if
     !
     read(fi,*) iqv(1:2)
     if(iqv(1) /= nbep) then
        write(*,*) "Stop in read_elph. # of bands. ik = ", dsp + ik
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        call MPI_FINALIZE(ierr)
        stop
     end if
     if(iqv(2) /= fbep) then
        write(*,*) "Stop in read_elph. # of bands. ik = ", dsp + ik
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        call MPI_FINALIZE(ierr)
        stop
     end if
     !
     read(fi,*) qvec(1:3)
     !
     if(any(abs(qvec(1:3) - rkv(1:3,ik + dsp)) > 1d-5)) then
        write(*,*) "Stop in read_elph.  rkv = ", &
        &          rkv(1:3,ik + dsp), " ik = ", dsp + ik
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        call MPI_FINALIZE(ierr)
        stop
     end if
     !
     read(fi,*) nb0
     if(nb0 /= nm) then
        write(*,*) "Stop in read_elph. # of modes. ik = ", dsp + ik
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        call MPI_FINALIZE(ierr)
        stop
     end if
     !     
     read(fi,*) omg0(1:nm,ik)
     read(fi,*) gg0(1:nm,fbep:lbep,fbep:lbep,1:nk,ik)
     !
     close(fi)
     !
  end do
  !$OMP END DO
  !
  !$OMP END PARALLEL
  !
end subroutine read_elph
!
! Compute cnt and dsp
!
subroutine cnt_and_dsp(n,cnt1,dsp1)
  !
  use ultrasonic_vals, only : petot, my_rank
  !
  integer,intent(in) :: n
  integer,intent(out) :: cnt1, dsp1
  !
  integer :: ii
  integer :: cnt(0:petot-1), dsp(0:petot-1)
  !
  cnt(0:petot-1)        = n / petot
  cnt(0:mod(n,petot)-1) = n / petot + 1
  dsp(0) = 0
  do ii = 1, petot - 1
     dsp(ii) = dsp(ii - 1) + cnt(ii - 1)
  end do
  !
  cnt1 = cnt(my_rank)
  dsp1 = dsp(my_rank)
  !
end subroutine cnt_and_dsp
!
! Irreducible Brillouin zone
!
subroutine irr_bz()
  !
  use ultrasonic_vals, only : my_rank, nk, nk0, ng, rkv, grid, nsym, sym
  !
  integer :: i1, i2, i3, ik, isym, nkk
  real(8) :: kv0(3,nk), kv1(3), kv2(3)
  !
  allocate(grid(3,nk))
  !
  nk0 = 0
  nkk = 0
  !
  do i3 = 1, ng(3)
     do i2 = 1, ng(2)
        do i1 = 1, ng(1)
           !
           nkk = nkk + 1
           grid(1:3,nkk) = dble((/i1, i2, i3/) - 1) / dble(ng(1:3))
           !
           kv1(1:3) = (dble((/i1, i2, i3/)) - 0.5d0) / dble(ng(1:3))
           kv1(1:3) = kv1(1:3) - dble(floor(kv1(1:3) + 0.5d0 + 1d-4))
           !
           do isym = 1, nsym
              !
              kv2(1:3) = matmul(dble(sym(1:3,1:3,isym)), kv1(1:3))
              kv2(1:3) = kv2(1:3) - dble(floor(kv2(1:3) + 0.5d0 + 1d-4))
              !
              do ik = 1, nk0
                 if(all(abs(kv2(1:3) - kv0(1:3,ik)) < 1d-8)) goto 10
              end do ! ik
              !
           end do ! isym
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
  if(my_rank == 0) write(*,*) "  # of irreducible k points : ", nk0
  !
  allocate(rkv(3,nk0))
  !
  rkv(1:3,1:nk0) = kv0(1:3,1:nk0)
  !
end subroutine irr_bz
!
! Read Screened Coulomb matrix elements
!
subroutine read_Coulomb()
  !
  use mpi
  use ultrasonic_vals, only : my_rank, nk0, nk, nb, Vc0, rkv, ng, ncf, nmf
  !use ultrasonic_routines, only : cnt_and_dsp
  use omp_lib
  !
  integer :: fi, ik, iqv(3), nb0, ierr, cnt, dsp
  real(8) :: qvec(3)
  character(100) :: fname
  !
  if(my_rank == 0) then
     !
     open(10, file = "vel1.dat", form = 'unformatted')
     !
     read(10) iqv(1:3)
     !
     read(10) nb0
     !
     read(10) qvec(1:3)
     !
     read(10) ncf
     write(*,*) "  Dimention of Chebyshev interpolation : ", ncf
     !
     close(10)
     !
  end if
  !
  call MPI_BCAST(ncf, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  !
  call cnt_and_dsp(nk0,cnt,dsp)
  allocate(vc0(ncf,nb,nb,nk,cnt))
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(cnt,dsp,nk0,nk,nb,vc0,ng,rkv,ncf) &
  !$OMP & PRIVATE(ierr,fi,ik,fname,nb0,iqv,qvec)
  !
  fi = OMP_GET_THREAD_NUM() + 10
  !
  !$OMP DO
  do ik = 1, cnt
     !
     write(fname,*) dsp + ik
     write(fname,'(3a)') "vel", trim(adjustl(fname)), ".dat"
     open(fi, file = trim(fname), form = 'unformatted')
     !
     read(fi) iqv(1:3)
     if(.not. all(iqv(1:3) == ng(1:3))) then
        write(*,*) "Stop in rea_coulomb. kgrid. ik = ", dsp + ik
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        call MPI_FINALIZE(ierr)
        stop
     end if
     !
     read(fi) nb0
     if(nb /= nb0) then
        write(*,*) "Stop in read_coulomb. # of bands. ik = ", dsp + ik
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        call MPI_FINALIZE(ierr)
        stop
     end if
     !
     read(fi) qvec(1:3)
     !  
     if(any(abs(qvec(1:3) - rkv(1:3,ik + dsp)) > 1d-5)) then
        write(*,*) "Stop in raad_coulomb. q point. ik = ", dsp + ik        
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        call MPI_FINALIZE(ierr)
        stop
     end if
     !
     read(fi) nb0
     if(nb0 /= ncf) then
        write(*,*) "Stop in read_coulomb. # of Chevyshev interpolation. ik = ", dsp + ik        
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        call MPI_FINALIZE(ierr)
        stop
     end if
     !
     read(fi) vc0(1:ncf,1:nb,1:nb,1:nk,ik)
     !
     close(fi)
     !
  end do
  !$OMP END DO
  !
  !$OMP END PARALLEL
  !
  if(nmf < 0) then
     vc0(1:ncf,1:nb,1:nb,1:nk,1:cnt) = 0d0
     nmf = 0
  end if
  !
end subroutine read_Coulomb
!
! Average matrix in grid
!
subroutine average_matrix()
  !
  use mpi
  use omp_lib
  use ultrasonic_vals, only : petot, my_rank, nk0, nk, nkd, ngd, fstk, lstk, &
  &                  fstb, lstb, nbf, nbep, fbep, lbep, nm, &
  &                  ncf, Vc0, gg0, omg0, Vcf, ggf, omgf, nbee, fbee, lbee
  !use ultrasonic_routines, only : interpol_indx_weight, &
  !&                  average_matrix_average, average_matrix_symmetrize
  !
  integer :: iq, ik, jk, ik2, ii, cnt, dsp, cntmax, ierr, ipe, nindmax, &
  &          nqs(nk0), nind(nk,nk), i1, i2, i3
  real(8) :: dgrid(3,nkd)
  integer,allocatable :: qindx(:), iks(:,:,:), ind(:,:,:,:)
  real(8),allocatable :: gg1(:,:,:,:,:), Vc1(:,:,:,:,:), omg1(:,:), wght(:,:), wght0(:,:)
  !
  ! Dense grid
  !
  ik = 0
  do i3 = 1, ngd(3)
     do i2 = 1, ngd(2)
        do i1 = 1, ngd(1)
           ik = ik + 1         
           dgrid(1:3,ik) = ((/i1, i2, i3/) - 1) / dble(ngd(1:3))
        end do
     end do
  end do
  !
  call cnt_and_dsp(nk0,  cnt,  dsp)
  !
  ! #####  Expand with Symmetries  #####
  !
  call average_matrix_average(nindmax,nind,ind)
  !
  ! #####  Communicate  #####
  !
  cntmax = cnt
  call MPI_allREDUCE(MPI_IN_PLACE, cntmax, 1, MPI_INTEGER, &
  &                  MPI_MAX, MPI_COMM_WORLD,ierr)
  !
  allocate(Vc1(ncf,fbee:lbee,nbf,nk,cntmax), &
  &        gg1( nm,fbep:lbep,nbf,nk,cntmax), &
  &       omg1( nm,                 cntmax), &
  &      qindx(                   cntmax+1))
  !
  allocate(Vcf(ncf,fbee:lbee,nk,nbf,fstk:lstk), &
  &        ggf( nm,fbep:lbep,nk,nbf,fstk:lstk), &
  &       omgf( nm,          nk,    fstk:lstk), &
  &      wght0(              nk,    fstk:lstk)  )
  !
  Vcf(1:ncf,fbee:lbee,1:nk,1:nbf,fstk:lstk) = 0d0
  ggf( 1:nm,fbep:lbep,1:nk,1:nbf,fstk:lstk) = 0d0
  omgf(1:nm,          1:nk,      fstk:lstk) = 0d0
  wght0(              1:nk,      fstk:lstk) = 0d0
  do ipe = 0, petot - 1
     !
     Vc1(1:ncf,fbee:lbee,1:nbf,1:nk,1:cntmax) = 0d0
     gg1( 1:nm,fbep:lbep,1:nbf,1:nk,1:cntmax) = 0d0
     omg1(1:nm,                     1:cntmax) = 0d0
     qindx(                         1:cntmax) = 1
     !
     if(ipe == my_rank) then
        !
        Vc1(1:ncf,fbee:lbee,1:nbf,1:nk,1:cnt) = Vc0(1:ncf,fbee:lbee,fstb:lstb,1:nk,1:cnt)
        gg1( 1:nm,fbep:lbep,1:nbf,1:nk,1:cnt) = gg0( 1:nm,fbep:lbep,fstb:lstb,1:nk,1:cnt)
        omg1(1:nm,                     1:cnt) = omg0(1:nm,                         1:cnt)
        !
        qindx(cntmax + 1) = cnt
        do iq = 1, cnt
           qindx(iq) = dsp + iq
        end do
        !
     end if
     !
     call MPI_BCAST(Vc1, ncf * nbee * nbf * nk * cntmax, MPI_DOUBLE_PRECISION, &
     &              ipe, MPI_COMM_WORLD, ierr)
     !
     call MPI_BCAST(gg1, nm * nbep * nbf * nk * cntmax, MPI_DOUBLE_PRECISION, &
     &              ipe, MPI_COMM_WORLD, ierr)
     !
     call MPI_BCAST(omg1, nm * cntmax, MPI_DOUBLE_PRECISION, &
     &              ipe, MPI_COMM_WORLD, ierr)
     !
     call MPI_BCAST(qindx, cntmax + 1, MPI_INTEGER, &
     &              ipe, MPI_COMM_WORLD, ierr)
     !
     !$OMP PARALLEL DEFAULT(NONE) &
     !$OMP & SHARED(fstk,lstk,cntmax,qindx,ncf,fbee,lbee,nm,nbf,fbep,lbep, &
     !$OMP &        Vc1,Vcf,omg1,omgf,gg1,ggf,dgrid,nind,nindmax,ind,wght0) &
     !$OMP & PRIVATE(ik,jk,iq,ik2,ii,iks,wght,nqs)
     !
     !$OMP DO
     do ik = fstk, lstk
        !
        ! #####  Interpol & Symmetrize  #####
        !
        call average_matrix_symmetrize(dgrid(1:3,ik),nind,nindmax,ind,nqs,iks,wght)
        !
        do iq = 1, qindx(cntmax + 1)
           !
           do ii = 1, nqs(qindx(iq))
              !
              jk  = iks(1, ii, qindx(iq))
              ik2 = iks(2, ii, qindx(iq))
              !
              wght0(jk,ik) = wght0(jk,ik) + wght(ii, qindx(iq))
              !
              Vcf(1:ncf,fbee:lbee,jk,1:nbf,ik) = Vcf(1:ncf,fbee:lbee,jk,1:nbf,ik) &
              &    + Vc1(1:ncf,fbee:lbee,1:nbf,ik2,iq) * wght(ii, qindx(iq))
              !
              ggf(1:nm,fbep:lbep,jk,1:nbf,ik) = ggf(1:nm,fbep:lbep,jk,1:nbf,ik) &
              &        + gg1(1:nm,fbep:lbep,1:nbf,ik2,iq) * wght(ii, qindx(iq))
              !
              omgf(1:nm,jk,ik) = omgf(1:nm,jk,ik) &
              &        + omg1(1:nm,iq) * wght(ii, qindx(iq))
              !
           end do ! ii
           !
        end do ! iq
        !
        deallocate(iks,wght)
        !
     end do ! ik
     !$OMP END DO
     !$OMP END PARALLEL
     !
  end do ! ipe
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fstk,lstk,ncf,nm,fbee,lbee,nbf,nk,Vcf,ggf,omgf,wght0,fbep,lbep) &
  !$OMP & PRIVATE(ik,jk)
  !
  !$OMP DO
  do ik = fstk, lstk
     do jk = 1, nk
        Vcf(1:ncf,fbee:lbee,jk,1:nbf,ik) = Vcf(1:ncf,fbee:lbee,jk,1:nbf,ik) / wght0(jk,ik)
        ggf( 1:nm,fbep:lbep,jk,1:nbf,ik) = ggf( 1:nm,fbep:lbep,jk,1:nbf,ik) / wght0(jk,ik)
        omgf(1:nm,          jk,      ik) = omgf(1:nm,          jk,      ik) / wght0(jk,ik)
     end do ! jk
  end do ! ik
  !$OMP END DO
  !$OMP END PARALLEL
  !
  deallocate(gg0, Vc0, omg0)
  deallocate(gg1, Vc1, omg1, ind, wght0)
  !
end subroutine average_matrix
!
! Average
!
subroutine average_matrix_average(nindmax,nind,ind)
  !
  use mpi
  use ultrasonic_vals, only : nk, nk0, nsym, sym, rkv, grid, ng, my_rank
  !
  integer,intent(out) :: nind(nk,nk), nindmax
  integer,intent(out),allocatable :: ind(:,:,:,:)
  !
  integer :: ik, iq, isym, ierr, jk1, ik1, ikv1(3)
  real(8) :: kv1(3)
  !
  ! Query of memory size
  !
  nind(1:nk,1:nk) = 0
  !
  do ik = 1, nk
     ! 
     do iq = 1, nk0
        !
        do isym = 1, nsym
           !
           ! rotate k
           !
           kv1(1:3) = matmul(dble(sym(1:3,1:3,isym)), &
           &                 grid(1:3,ik)             ) * dble(ng(1:3))
           ikv1(1:3) = nint(kv1(1:3))
           !
           if(any(abs(kv1(1:3) - dble(ikv1(1:3))) > 1d-5)) cycle
           !
           ikv1(1:3) = modulo(ikv1(1:3), ng(1:3))
           ik1 = 1 + ikv1(1) + ng(1) * ikv1(2) + ng(1) * ng(2) * ikv1(3)
           !
           ! Rotate k + q
           !
           kv1(1:3) = matmul(dble(sym(1:3,1:3,isym)),  &
           &                 grid(1:3,ik) + rkv(1:3,iq)) * dble(ng(1:3))
           kv1(1:3) = kv1(1:3) - 0.5d0
           ikv1(1:3) = nint(kv1(1:3))
           !
           if(any(abs(kv1(1:3) - dble(ikv1(1:3))) > 1d-5)) cycle
           !
           ikv1(1:3) = modulo(ikv1(1:3), ng(1:3))
           jk1 = 1 + ikv1(1) + ng(1) * ikv1(2) + ng(1) * ng(2) * ikv1(3)
           !
           nind(jk1,ik1) = nind(jk1,ik1) + 1
           !
        end do ! isym
        !
     end do ! iq
     !
  end do ! ik
  !
  nindmax = maxval(nind(1:nk,1:nk))
  allocate(ind(2,nindmax,nk,nk))
  !
  ind(1:2,1:nindmax,1:nk,1:nk) = 0
  nind(1:nk,1:nk) = 0
  !
  call MPI_allREDUCE(nindmax, ik, 1, MPI_INTEGER, &
  &                  MPI_MAX, MPI_COMM_WORLD,ierr)
  if(my_rank == 0) write(*,*) "  nindmax : ", nindmax
  !
  do ik = 1, nk
     ! 
     do iq = 1, nk0
        !
        do isym = 1, nsym
           !
           ! rotate k
           !
           kv1(1:3) = matmul(dble(sym(1:3,1:3,isym)), &
           &                 grid(1:3,ik)             ) * dble(ng(1:3))
           ikv1(1:3) = nint(kv1(1:3))
           !
           if(any(abs(kv1(1:3) - dble(ikv1(1:3))) > 1d-5)) cycle
           !
           ikv1(1:3) = modulo(ikv1(1:3), ng(1:3))
           ik1 = 1 + ikv1(1) + ng(1) * ikv1(2) + ng(1) * ng(2) * ikv1(3)
           !
           ! Rotate k + q
           !
           kv1(1:3) = matmul(dble(sym(1:3,1:3,isym)),  &
           &                 grid(1:3,ik) + rkv(1:3,iq)) * dble(ng(1:3))
           kv1(1:3) = kv1(1:3) - 0.5d0
           ikv1(1:3) = nint(kv1(1:3))
           !
           if(any(abs(kv1(1:3) - dble(ikv1(1:3))) > 1d-5)) cycle
           !
           ikv1(1:3) = modulo(ikv1(1:3), ng(1:3))
           jk1 = 1 + ikv1(1) + ng(1) * ikv1(2) + ng(1) * ng(2) * ikv1(3)
           !
           nind(jk1,ik1) = nind(jk1,ik1) + 1
           ind(1:2,nind(jk1,ik1),jk1,ik1) = (/ik, iq/)
           !
        end do ! isym
        !
     end do ! iq
     !
  end do ! ik
  !
end subroutine average_matrix_average
!
! Symmetrize
!
subroutine average_matrix_symmetrize(dgrid,nind,nindmax,ind,nqs,iks,wght)
  !
  use mpi
  use omp_lib
  use ultrasonic_vals, only : nk, nk0, fstk, lstk, nsym, grid, ng, sym, nsym, my_rank,  nkd, ngd
  !use ultrasonic_routines, only : interpol_indx_weight
  !
  integer,intent(in) :: nind(nk,nk), nindmax, ind(2,nindmax,nk,nk)
  real(8),intent(in) :: dgrid(3)
  integer,intent(out) :: nqs(nk0)
  integer,intent(out),allocatable :: iks(:,:,:)
  real(8),intent(out),allocatable :: wght(:,:)
  !
  integer :: ii, jj, jk, isym, ik2, jk1, iq, jkv1(3), iki(8), nqsmax
  real(8) :: wi(8), kv1(3)
  !
  nqs(1:nk0) = 0
  !
  ! Query of memory size
  !
  do isym = 1, nsym
     !
     kv1(1:3) = matmul(dble(sym(1:3,1:3,isym)), dgrid(1:3))
     call interpol_indx_weight(ng,kv1,iki,wi)
     !
     do jk = 1, nk
        !
        kv1(1:3) = matmul(dble(sym(1:3,1:3,isym)), grid(1:3,jk) + 0.5d0 / dble(ng(1:3)))
        kv1(1:3) = kv1(1:3) * dble(ng(1:3)) - 0.5d0
        jkv1(1:3) = nint(kv1(1:3))
        !
        if(any(abs(kv1(1:3) - dble(jkv1(1:3))) > 1d-5)) cycle
        !
        jkv1(1:3) = modulo(jkv1(1:3), ng(1:3))
        jk1 = 1 + jkv1(1) + ng(1) * jkv1(2) + ng(1) * ng(2) * jkv1(3)
        !
        do ii = 1, 8
           do jj = 1, nind(jk1,iki(ii))
              !
              iq = ind(2,jj,jk1,iki(ii))
              !
              nqs(iq) = nqs(iq) + 1
              !
           end do ! jj
        end do ! ii
        !
     end do ! jk
     !
  end do ! isym
  !
  nqsmax = maxval(nqs(1:nk0))
  !
  allocate(wght(nqsmax,nk0), iks(2,nqsmax,nk0))
  !
  iks(1:2,1:nqsmax,1:nk0) = 0
  wght(   1:nqsmax,1:nk0) = 0d0
  nqs(             1:nk0) = 0
  !
  do isym = 1, nsym
     !
     kv1(1:3) = matmul(dble(sym(1:3,1:3,isym)), dgrid(1:3))
     call interpol_indx_weight(ng,kv1,iki,wi)
     !
     do jk = 1, nk
        !
        kv1(1:3) = matmul(dble(sym(1:3,1:3,isym)), grid(1:3,jk) + 0.5d0 / dble(ng(1:3)))
        kv1(1:3) = kv1(1:3) * dble(ng(1:3)) - 0.5d0
        jkv1(1:3) = nint(kv1(1:3))
        !
        if(any(abs(kv1(1:3) - dble(jkv1(1:3))) > 1d-5)) cycle
        !
        jkv1(1:3) = modulo(jkv1(1:3), ng(1:3))
        jk1 = 1 + jkv1(1) + ng(1) * jkv1(2) + ng(1) * ng(2) * jkv1(3)
        !
        do ii = 1, 8
           do jj = 1, nind(jk1,iki(ii))
              !
              ik2 = ind(1,jj,jk1,iki(ii))
              iq  = ind(2,jj,jk1,iki(ii))
              !
              nqs(iq) = nqs(iq) + 1
              iks( 1:2, nqs(iq), iq) = (/jk, ik2/)
              wght(     nqs(iq), iq) = wi(ii)
              !
           end do ! jj
        end do ! ii
        !
     end do ! jk
     !
  end do ! isym
  !
end subroutine average_matrix_symmetrize
!
! Bi-linear interpolation
!
subroutine interpol_indx_weight(ng,ko,iki,wi)
  !
  !
  integer,intent(in)  :: ng(3)
  real(8),intent(in)  :: ko(3)
  integer,intent(out) :: iki(8)
  real(8),intent(out) :: wi(8)
  !
  integer :: ikv0(3), ikv1(3), i1, i2, i3, ii
  real(8) :: x(3)
  !
  ! Search nearest neighbor grid points.
  !
  x(1:3) = ko(1:3) * dble(ng(1:3))
  ikv0(1:3) = floor(x(1:3))
  x(1:3) = x(1:3) - dble(ikv0(1:3))
  !
  ! Interpolation
  !
  ii = 0
  do i1 = 0, 1
     do i2 = 0, 1
        do i3 = 0, 1
           !
           ii = ii + 1
           !
           ikv1(1:3) = ikv0(1:3) + (/i1, i2, i3/)
           ikv1(1:3) = modulo(ikv1(1:3), ng(1:3))
           iki(ii) = 1 + ikv1(1) + ng(1) * ikv1(2) + ng(1) * ng(2) * ikv1(3)
           !
           wi(ii) = x(1)**i1 * (1d0 - x(1))**(1 - i1) &
           &      * x(2)**i2 * (1d0 - x(2))**(1 - i2) &
           &      * x(3)**i3 * (1d0 - x(3))**(1 - i3) 
           !
        end do
     end do
  end do
  !
end subroutine interpol_indx_weight
!
! Define initial delta
!
subroutine read_delta()
  !
  use mpi
  use ultrasonic_vals, only : my_rank, nt, xi, delta, dk, kindx, bindx
  !
  integer :: it, fi = 10, is, ierr 
  real(8) :: Z0
  character(1) :: tmp
  !
  if(my_rank == 0) then
     !
     open(fi, file = 'delta.dat',status="old", action = 'read',iostat = is)
     !
     read(fi,*) tmp, nt
     !
     write(*,*) "  # of total points for gap equation : ", nt
     allocate(xi(nt), delta(nt), dk(nt), kindx(nt), bindx(nt))
     !
     do it = 1, nt
        !
        read(fi,*) xi(it), delta(it), Z0, dk(it), kindx(it), bindx(it)
        !
     enddo
     !
     close(fi)
     !
  end if
  !
  call MPI_BCAST(nt,   1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  if(my_rank /= 0) allocate(xi(nt), delta(nt), dk(nt), kindx(nt), bindx(nt))
  call MPI_BCAST(xi,    nt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(delta, nt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(dk,    nt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(kindx, nt, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(bindx, nt, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  !
end subroutine read_delta
!
! Calc Z_{n k}
!
subroutine make_Z()
  !
  use mpi
  use ultrasonic_vals, only : nkd, nbf, ngd, nm, nt, xi, dk, kindx, bindx, omgf, ggf, Z, &
  &                  fstk, lstk, fbep, lbep, beta
  !use ultrasonic_routines, only : Zweight
  !
  integer :: ik, ib, jt, jk, jb, im, ierr
  real(8) :: xp, om, z0, txp, tom
  !
  allocate(Z(nbf,nkd))
  Z(1:nbf,1:nkd) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fstk,lstk,nbf,nm,nt,xi,dk,kindx,bindx,Z,ggf,omgf,fbep,lbep,beta) &
  !$OMP & PRIVATE(ik, ib, jt, jk, jb, im, xp, om, z0, txp, tom)
  !
  !$OMP DO
  do ik = fstk, lstk
     !
     do ib = 1, nbf
        !
        do jt = 1, nt
           !
           xp = abs(xi(jt) * beta * 0.5d0)
           txp = tanh(xp)
           jk = kindx(jt)
           jb = bindx(jt)
           !
           if(jb < fbep .or. lbep < jb) cycle
           !
           do im = 1, nm
              !
              om = abs(omgf(im,jk,ik) * beta * 0.5d0)
              tom = tanh(om)
              z0 = beta**2 * Zweight(xp, om, txp, tom)
              !    
              Z(ib,ik) = Z(ib,ik) &
              &        + dk(jt) * ggf(im,jb,jk,ib,ik) * Z0
              !
           end do ! im
           !
        enddo ! jt
        !
        Z(ib,ik) = - Z(ib,ik)
        !       
     end do ! ib
     !
  end do ! ik
  !$OMP END DO
  !
  !$OMP END PARALLEL
  !
  call MPI_allREDUCE(MPI_IN_PLACE, Z, nbf * nkd, &
  &                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  !
end subroutine make_Z
!
! Weight for Z
!
function Zweight(y,z,ty,tz) result(Wz)
  !
  !
  real(8),intent(in) :: y, z, ty, tz
  real(8) :: Wz
  !
  real(8) :: thr = 1d-8, mp, pp
  !
  pp = 1d0 / (  y + z)
  mp = 1d0 / (- y + z)
  !
  if(abs(y + z) < thr) then
     !
     Wz = 0.125d0 * ( &
     &      (2d0 - 3d0 * ty**2) * (1d0 - ty**2) / (3d0 * ty) &
     &    + ((tz**2 - 1d0) + (-2d0 + (tz + ty) * tz + 2d0 * (tz - ty) * mp) * mp / tz) * mp &
     &              )
     !
  else if(abs(- y + z) < thr) then
     !
     Wz = 0.125d0 * ( &
     &      ((tz**2 - 1d0) + (-2d0 + (tz - ty) * tz + 2d0 * (tz + ty) * pp) * pp / tz) * pp &
     &    + (2d0 - 3d0 * ty**2) * (1d0 - ty**2) / (- 3d0 * ty) &
     &              )
     !
  else
     !
     Wz = 0.125d0 * ( &
     &      ((tz**2 - 1d0) + (-2d0 + (tz - ty) * tz + 2d0 * (tz + ty) * pp) * pp / tz) * pp &
     &    + ((tz**2 - 1d0) + (-2d0 + (tz + ty) * tz + 2d0 * (tz - ty) * mp) * mp / tz) * mp &
     &              )
     !
  end if
  !
end function Zweight
!
! Make Kel, Kph & a part of Jacobian
!
subroutine gapeq()
  !
  use mpi
  use ultrasonic_vals, only : pi, my_rank, nkd, nt, ngd, nb, nbf, nm, ncf, omgf, ggf, vcf, xic, &
  &                  xi, dk, delta, kindx, bindx, dltf, Z, beta, fbep, lbep, nmf, mf, wmf, &
  &                  fstk, lstk
  !use ultrasonic_routines, only : Kweight, calc_Kel, weightspoints_gl
  !
  integer :: ik, ib, jt, jk, jb, im, ierr, nh
  real(8) :: xp, om, eqp, Kph, Kel, Kave, dlth, dosh, xmax, bxp, txp, tom
  !
  allocate(dltf(nbf,nkd))
  dltf(1:nbf,1:nkd) = 0d0
  !
  ! Set Matsubara sum
  !
  allocate(mf(nmf), wmf(nmf))
  !
  call weightspoints_gl(nmf,mf,wmf)
  !
  wmf(1:nmf) = 2d0 / (pi * (mf(1:nmf)**2 + 1d0)) * wmf(1:nmf)
  mf(1:nmf) = (1d0 + mf(1:nmf)) / (1d0 - mf(1:nmf))
  !
  ! Extrapolation of high energy region
  !
  nh = count(xi(1:nt) > xic)
  dlth = sum(delta(1:nt) / xi(1:nt), xi(1:nt) > xic) &
  &    / sum(1d0 / xi(1:nt)**2,      xi(1:nt) > xic)
  !
  xmax = maxval(xi(1:nt))
  dosh = sum(dk(1:nt), xi(1:nt) > xic) / (xmax - xic)
  !
  if(my_rank == 0) write(*,*) dlth / xic * 13.6d3
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fstk, lstk, nt, nb, nbf, nm, ncf, omgf, ggf, vcf, &
  !$OMP &        xi, dk, delta, kindx, bindx, dltf, Z, beta, fbep, lbep, &
  !$OMP &        nh, dlth, dosh, xic, xmax) &
  !$OMP & PRIVATE(ik, ib, jt, jk, jb, im, bxp, txp, tom, &
  !$OMP &         xp, om, eqp, Kph, Kel, Kave)
  !
  !$OMP DO
  do ik = fstk, lstk
     !
     do ib = 1, nbf
        !
        Kave = 0d0
        !
        do jt = 1, nt
           !
           xp = abs(xi(jt))
           jk = kindx(jt)
           jb = bindx(jt)
           bxp = 0.5d0 * beta * xp
           txp = tanh(bxp)
           !
           Kph = 0d0
           !
           if(fbep <= jb .and. jb <= lbep) then
              !
              do im = 1, nm
                 !
                 om = abs(0.5d0 * beta * omgf(im,jk,ik))
                 tom = tanh(om)
                 !
                 Kph = Kph + ggf(im,jb,jk,ib,ik) &
                 &   * beta * Kweight(bxp, om, txp, tom)
                 !
              end do ! im
              Kph = Kph * 2d0
              !
           end if
           !
           Kel = calc_Kel(xp,Vcf(1:ncf,jb,jk,ib,ik))
           !
           eqp = sqrt(delta(jt)**2 + xp**2)
           !
           dltf(ib,ik) = dltf(ib,ik) &
           &                 - 0.5d0 * dk(jt) * (Kph + Kel) / (1d0 + Z(ib,ik)) &
           &                   * delta(jt) / eqp * tanh(0.5d0 * beta * eqp)
           !
           if(xi(jt) > xic) Kave = Kave + (Kph + Kel)
           !
        enddo ! jt
        !
        if(xic > 0) then
           !
           Kave = Kave / dble(nh)
           !
           dltf(ib,ik) = dltf(ib,ik) - 0.5d0 * Kave / (1d0 + Z(ib,ik)) &
           &                          * dosh * dlth / xmax
           !
        end if ! (xic > 0)
        !
     end do ! ib
     !
  end do ! ik
  !$OMP END DO
  !
  !$OMP END PARALLEL
  !
  call MPI_allREDUCE(MPI_IN_PLACE, dltf, nbf * nkd, &
  &                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  !
  dltf(1:nbf,1:nkd) = dltf(1:nbf,1:nkd) * beta

!eqp = maxval( dltf(1:nbf,1:nkd))
!dltf(1:nbf,1:nkd) = dltf(1:nbf,1:nkd) * (dltf(1:nbf,1:nkd) / eqp)**3

  !
  deallocate(ggf, vcf, omgf)
  !
end subroutine gapeq
!
! Weights & Points for Gauss-Legendre method
!
subroutine weightspoints_gl(n,x,w)
  !
  use mpi
  use ultrasonic_vals, only: pi
  !use ultrasonic_routines, only : legendre
  !
  integer,intent(in) :: n
  real(8),intent(out) :: x(n), w(n)
  !
  integer :: i, itr, ierr
  real(8) :: Pn, Pnm1, PP
  !
  do i = 1, n
     !
     x(i) = cos(pi * (dble(n - i + 1) - 0.25d0) / (dble(n) + 0.5d0))
     !
     do itr = 1, 100
        !
        call legendre(n    , x(i), Pn  )
        call legendre(n - 1, x(i), Pnm1)
        !
        PP = dble(n) * (Pnm1 - x(i) * Pn) / (1d0 - x(i)**2)
        !
        x(i) = x(i) - Pn / PP
        !
        if(abs(Pn/PP) < 1d-10) exit
        !
     end do
     !
     if(itr >= 100) then
        write(*,*) "Stop in weightspoints_gl. newton."
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        call MPI_FINALIZE(itr)
        stop
     end if
     !
     call legendre(n - 1, x(i), Pnm1)     
     w(i) = 2d0 * (1d0 - x(i)**2) / (dble(n) * Pnm1)**2
     !
  end do
  !
end subroutine weightspoints_gl
!
! Legendre polynomial
!
subroutine legendre(n,x,P)
  !
  !
  integer,intent(in) :: n
  real(8),intent(in) :: x
  real(8),intent(out) :: P
  !
  integer :: i
  real(8) :: Pold, Pnew
  !
  if(n == 0) then
     P = 1
     return
  else if(n == 1) then
     P = x
     return
  end if
  !
  Pold = 1d0
  Pnew = x
  !
  do i = 1, n - 1
     !
     P = (dble(2 * i + 1) * x * Pnew - dble(i) * Pold) / dble(i + 1)
     Pold = Pnew
     Pnew = P
     !
  end do
  !
end subroutine legendre
!
! Weight for K
!
function Kweight(y,z,ty,tz) result(Wk)
  !
  !
  real(8),intent(in) :: y, z, ty, tz
  real(8) :: Wk
  !
  real(8) :: thr = 1d-4, zp, zm
  !
  if(abs(y) < thr) then
     !
     Wk = ((-1d0 / z + 1d0 / tz) / z - 0.5d0) / z
     !
  else
     !
     zp = 1d0 / (  y + z)
     zm = 1d0 / (- y + z)
     !
     if (abs(y + z) < thr) then
        !
        Wk = ( &
        &        (- 1d0 + ty * tz + (- ty + tz) * zm) * (- zm) &
        &    ) / (4d0 * ty * tz)
        !
     else if (abs(- y + z) < thr) then
        !
        Wk = ( &
        &        (- 1d0 - ty * tz + (  ty + tz) * zp) * (  zp) &
        &    ) / (4d0 * ty * tz)
        !
     else
        !
        Wk = ( &
        &        (-1d0 - ty * tz + (  ty + tz) * zp) * (  zp) &
        &      + (-1d0 + ty * tz + (- ty + tz) * zm) * (- zm) &
        &    ) / (4d0 * ty * tz)
        !
     end if
     !
  end if
  !
end function Kweight
!
! Compute Kel
!
function calc_Kel(xp,Vci) result(Kel)
  !
  use ultrasonic_vals, only : pi, ncf, nmf, mf, wmf
  !
  real(8),intent(in) :: xp, Vci(ncf)
  real(8) :: Kel
  !
  real(8),parameter :: &
  &  ints(1:20) = (/ 0d0,  1d0,  2d0,  3d0,  4d0, &
  &                  5d0,  6d0,  7d0,  8d0,  9d0, &
  &                 10d0, 11d0, 12d0, 13d0, 14d0, &
  &                 15d0, 16d0, 17d0, 18d0, 19d0  /)
  !
  integer :: imf
  real(8) :: xpxp, Vel0, Vel, x0, xcf
  !
  x0 = cos(pi / dble(2 * ncf))
  !
  xpxp = abs(xp)
  !
  xcf = acos(-x0)
  Vel0 = dot_product(Vci(1:ncf), cos(ints(1:ncf) * xcf))
  Kel = Vel0
  !
  do imf = 1, nmf
     !
     xcf = mf(imf) * xpxp
     xcf = (xcf - 1d0) / (xcf + 1d0) * x0
     xcf = acos(xcf)
     !
     Vel = dot_product(Vci(1:ncf), cos(ints(1:ncf) * xcf))
     !
     Kel = Kel + wmf(imf) * (Vel - Vel0)
     !
  end do
  !
end function calc_Kel
!
! Fermi verocity
!
subroutine calc_fvel()
  !
  use mpi
  use ultrasonic_vals, only : nkd, ngd, nbf, avec, eig, vf, dgrid, fstb, my_rank, fstk, lstk
  !
  !
  integer :: ib, ii, ik, ikp(3), ikm(3), ierr, i1, i2, i3
  real(8) :: de(3)
  !
  allocate(vf(3,nbf, nkd), dgrid(3,nkd))
  vf(1:3,1:nbf,1:nkd) = 0d0
  !
  ik = 0
  do i3 = 1, ngd(3)
     do i2 = 1, ngd(2)
        do i1 = 1, ngd(1)
           ik = ik + 1         
           dgrid(1:3,ik) = (/i1, i2, i3/) - 1
        end do
     end do
  end do
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fstk,lstk,ngd,nbf,avec,eig,vf,dgrid,fstb) &
  !$OMP & PRIVATE(ik,ib,ii,ikp,ikm,de)
  !
  !$OMP DO
  do ik = fstk, lstk
     do ib = 1, nbf
        !
        do ii = 1, 3
           !
           ikp(1:3) = dgrid(1:3,ik)
           ikp(ii) = ikp(ii) + 1
           ikp(1:3) = modulo(ikp(1:3), ngd(1:3)) + 1
           !                  
           ikm(1:3) = dgrid(1:3,ik)
           ikm(ii) = ikm(ii) - 1
           ikm(1:3) = modulo(ikm(1:3), ngd(1:3)) + 1
           !                 
           de(ii) = eig(fstb - 1 + ib,ikp(1),ikp(2),ikp(3)) &
           &      - eig(fstb - 1 + ib,ikm(1),ikm(2),ikm(3))
           de(ii) = 0.5d0 * de(ii) * dble(ngd(ii))
           !
        end do
        !
        de(1:3) = matmul(avec(1:3,1:3), de(1:3))
        vf(1,ib,ik) = de(1)
        vf(2,ib,ik) = de(1) + de(2)
        vf(3,ib,ik) = de(3)
        !
     end do
     !
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  !
  call MPI_allREDUCE(MPI_IN_PLACE, vf, 3 * nbf * nkd, &
  &                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
if(my_rank == 0 ) write(60,'(3e15.5)') vf(1:3,1:nbf,1:nkd)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !
end subroutine calc_fvel
!
! define shortest diagonal line & define type of tetragonal
!
subroutine tetra_type()
  !
  use ultrasonic_vals, only : my_rank, ngd, bvec, wlsm, ivvec, ltetra
  !
  integer,parameter :: DP = 8
  integer :: itype, i1, i2, i3, it, divvec(4,4), ivvec0(4), ii
  real(DP) :: l(4), bvec2(3,3), bvec3(3,4)
  !
  do ii = 1, 3
     bvec2(1:3,ii) = bvec(1:3,ii) / real(ngd(ii), DP)
  end do
  !
  bvec3(1:3,1) = -bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
  bvec3(1:3,2) =  bvec2(1:3,1) - bvec2(1:3,2) + bvec2(1:3,3)
  bvec3(1:3,3) =  bvec2(1:3,1) + bvec2(1:3,2) - bvec2(1:3,3)
  bvec3(1:3,4) =  bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
  !
  ! length of delta bvec
  !
  do ii = 1, 4
     l(ii) = dot_product(bvec3(1:3,ii),bvec3(1:3,ii))
  end do
  !
  itype = minloc(l(1:4),1)
  !
  ! start & last
  !
  ivvec0(1:4) = (/ 0, 0, 0, 0 /)
  !
  divvec(1:4,1) = (/ 1, 0, 0, 0 /)
  divvec(1:4,2) = (/ 0, 1, 0, 0 /)
  divvec(1:4,3) = (/ 0, 0, 1, 0 /)
  divvec(1:4,4) = (/ 0, 0, 0, 1 /)
  !
  ivvec0(itype) = 1
  divvec(itype, itype) = - 1
  !
  it = 0
  do i1 = 1, 3
     do i2 = 1, 3
        if(i2 == i1) cycle
        do i3 = 1, 3
           if(i3 == i1 .or. i3 == i2) cycle
           !
           it = it + 1
           !
           ivvec(1:3,1,it) = ivvec0(1:3)
           ivvec(1:3,2,it) = ivvec(1:3,1,it) + divvec(1:3,i1)
           ivvec(1:3,3,it) = ivvec(1:3,2,it) + divvec(1:3,i2)
           ivvec(1:3,4,it) = ivvec(1:3,3,it) + divvec(1:3,i3)
           !
        end do
     end do
  end do
  !
  ivvec(1:3, 5,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6)
  ivvec(1:3, 6,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6)
  ivvec(1:3, 7,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6)
  ivvec(1:3, 8,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6)
  !
  ivvec(1:3, 9,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,3,1:6)
  ivvec(1:3,10,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,4,1:6)
  ivvec(1:3,11,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,1,1:6)
  ivvec(1:3,12,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,2,1:6)
  !
  ivvec(1:3,13,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,4,1:6)
  ivvec(1:3,14,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,1,1:6)
  ivvec(1:3,15,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,2,1:6)
  ivvec(1:3,16,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,3,1:6)
  !
  ivvec(1:3,17,1:6) =  ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6) + ivvec(1:3,2,1:6)
  ivvec(1:3,18,1:6) =  ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6) + ivvec(1:3,3,1:6)
  ivvec(1:3,19,1:6) =  ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6) + ivvec(1:3,4,1:6)
  ivvec(1:3,20,1:6) =  ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6) + ivvec(1:3,1,1:6)
  !
  if(ltetra == 1) then
     !
     if(my_rank == 0) write(*,*) "Linear tetrahedron method is used."
     !
     wlsm(1:4,1:20) = 0.0_dp
     wlsm(1,1:4) = real((/1, 0, 0, 0/), DP)
     wlsm(2,1:4) = real((/0, 1, 0, 0/), DP)
     wlsm(3,1:4) = real((/0, 0, 1, 0/), DP)
     wlsm(4,1:4) = real((/0, 0, 0, 1/), DP)
     !
  else if(ltetra == 2) then
     !
     if(my_rank == 0) write(*,*) "Improved (leveled) tetrahedron method is used."
     !
     wlsm(1, 1: 4) = dble((/1440,    0,   30,    0/))
     wlsm(2, 1: 4) = dble((/   0, 1440,    0,   30/))
     wlsm(3, 1: 4) = dble((/  30,    0, 1440,    0/))
     wlsm(4, 1: 4) = dble((/   0,   30,    0, 1440/))
     !
     wlsm(1, 5: 8) = dble((/ -38,    7,   17,  -28/))
     wlsm(2, 5: 8) = dble((/ -28,  -38,    7,   17/))
     wlsm(3, 5: 8) = dble((/  17,  -28,  -38,    7/))
     wlsm(4, 5: 8) = dble((/   7,   17,  -28,  -38/))
     !
     wlsm(1, 9:12) = dble((/ -56,    9,  -46,    9/))
     wlsm(2, 9:12) = dble((/   9,  -56,    9,  -46/))
     wlsm(3, 9:12) = dble((/ -46,    9,  -56,    9/))
     wlsm(4, 9:12) = dble((/   9,  -46,    9,  -56/))
     !
     wlsm(1,13:16) = dble((/ -38,  -28,   17,    7/))
     wlsm(2,13:16) = dble((/   7,  -38,  -28,   17/))
     wlsm(3,13:16) = dble((/  17,    7,  -38,  -28/))
     wlsm(4,13:16) = dble((/ -28,   17,    7,  -38/))
     !
     wlsm(1,17:20) = dble((/ -18,  -18,   12,  -18/))
     wlsm(2,17:20) = dble((/ -18,  -18,  -18,   12/))
     wlsm(3,17:20) = dble((/  12,  -18,  -18,  -18/))
     wlsm(4,17:20) = dble((/ -18,   12,  -18,  -18/))
     !
     wlsm(1:4,1:20) = wlsm(1:4,1:20) / 1260_dp
     !
  else 
     !
     stop "tetra_type"
     !
  end if
  !
end subroutine tetra_type
!
! Integration for lambda
!
subroutine calc_usonic()
  !
  use mpi
  use ultrasonic_vals, only : my_rank, fstk, lstk, dgrid, beta, &
  &                  nkd, ngd, nbf, eig, vf, dltf, wlsm, ivvec, fstb, lstb, beta
  !use ultrasonic_routines, only : sort, tetrahedra2
  !
  integer :: it, ib, ii, jj, ikv1(3), ik, ik1, ierr
  real(8) :: ei1(nbf,4), vf1(3,nbf,4), dl1(nbf,4), vf2(3,3), dl2(3), &
  &          e(4), a(4,4), V, usonic0(2,3), usonic(2,3)
  real(8) :: tmp(5,4), tmp2(5,3), dmin, dmax
  !
if(my_rank == 0) write(*,*) minval(dltf(1:nbf,1:nkd)), maxval(dltf(1:nbf,1:nkd)), &
&                sum(dltf(1:nbf,1:nkd)) / dble(nbf * nkd)
  usonic(1:2,1:3) = 0d0
  !
  dmax = -1d10
  dmin = 1d10
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fstk,lstk,ngd,nbf,wlsm,eig,vf,dltf,usonic,fstb,lstb,dgrid,ivvec,dmin,dmax) &
  !$OMP & PRIVATE(it,ib,ii,jj,tmp,tmp2,a,e,ei1,vf1,dl1,vf2,dl2,V,usonic0,ikv1,ik1,ik)
  !
  !$OMP DO REDUCTION(+:usonic) REDUCTION(max:dmax) REDUCTION(min:dmin)
  do ik = fstk, lstk
     !
     do it = 1, 6
        !
        ei1(    1:nbf,1:4) = 0d0
        vf1(1:3,1:nbf,1:4) = 0d0
        dl1(    1:nbf,1:4) = 0d0
        !
        do ii = 1, 20
           !
           ikv1(1:3) = dgrid(1:3,ik) + ivvec(1:3,ii,it)
           ikv1(1:3) = modulo(ikv1(1:3), ngd(1:3))
           ik1 = 1 + ikv1(1) + ngd(1) * ikv1(2) + ngd(1) * ngd(2) * ikv1(3)
           ikv1(1:3) = ikv1(1:3) + 1
           !
           do jj = 1, 4
              ei1(    1:nbf,jj) = ei1(    1:nbf,jj) + wlsm(jj,ii) * eig(fstb:lstb,ikv1(1),ikv1(2),ikv1(3))
              vf1(1:3,1:nbf,jj) = vf1(1:3,1:nbf,jj) + wlsm(jj,ii) * vf(1:3,1:nbf, ik1)
              dl1(    1:nbf,jj) = dl1(    1:nbf,jj) + wlsm(jj,ii) * dltf(  1:nbf, ik1)
           end do
        end do
        !
        do ib = 1, nbf
           !
           tmp(1,   1:4) = ei1(     ib, 1:4)
           tmp(2:4, 1:4) = vf1(1:3, ib, 1:4)
           tmp(5,   1:4) = dl1(     ib, 1:4)
           !
           call sort(5, 4, tmp)
           !
           e(1:4) = tmp(1, 1:4)
           !
           do ii = 1, 4
              a(ii,1:4) = (0d0 - e(1:4)) / (e(ii) - e(1:4))
           end do
           !
           if(e(1) < 0d0 .and. 0d0 <= e(2)) then
              !
              ! A
              !
              !V = 3d0 * a(2,1) * a(3,1) * a(4,1) / (0d0 - e(1))
              V = 3d0 * a(2,1) * a(3,1)          / (e(4) - e(1))
              !
              tmp2(1:5,1) = tmp(1:5,1) * a(1,2) + tmp(1:5,2) * a(2,1)
              tmp2(1:5,2) = tmp(1:5,1) * a(1,3) + tmp(1:5,3) * a(3,1) 
              tmp2(1:5,3) = tmp(1:5,1) * a(1,4) + tmp(1:5,4) * a(4,1)
              ! 
              vf2(1:3,1:3) = tmp2(2:4,1:3)
              dl2(    1:3) = tmp2(  5,1:3)
              !
              dmin = min(dmin, minval(dl2(1:3)))
              dmax = max(dmax, maxval(dl2(1:3)))
              ! 
              call tetrahedra2(vf2,dl2,usonic0)
              usonic(1:2,1:3) = usonic(1:2,1:3) + V * usonic0(1:2,1:3)
              !
           else if(e(2) < 0d0 .and. 0d0 <= e(3)) then
              !
              ! B - 1
              !
              !V = 3d0 * a(3,1) * a(4,1) * a(2,4) / (0d0 - e(1))
              V = 3d0 *           a(4,1) * a(2,4) / (e(3) - e(1))
              !
              tmp2(1:5,1) = tmp(1:5,1) * a(1,3) + tmp(1:5,3) * a(3,1) 
              tmp2(1:5,2) = tmp(1:5,1) * a(1,4) + tmp(1:5,4) * a(4,1) 
              tmp2(1:5,3) = tmp(1:5,2) * a(2,4) + tmp(1:5,4) * a(4,2) 
              !
              vf2(1:3,1:3) = tmp2(2:4,1:3)
              dl2(    1:3) = tmp2(  5,1:3)
              ! 
              dmin = min(dmin, minval(dl2(1:3)))
              dmax = max(dmax, maxval(dl2(1:3)))
              !
              call tetrahedra2(vf2,dl2,usonic0)
              usonic(1:2,1:3) = usonic(1:2,1:3) + V * usonic0(1:2,1:3)
              !
              ! B - 2
              !
              !V = 3d0 * a(2,3) * a(3,1) * a(4,2) / (0d0 - e(1))
              V = 3d0 * a(2,3)           * a(4,2) / (e(3) - e(1))
              !
              tmp2(1:5,1) = tmp(1:5,1) * a(1,3) + tmp(1:5,3) * a(3,1) 
              tmp2(1:5,2) = tmp(1:5,2) * a(2,3) + tmp(1:5,3) * a(3,2) 
              tmp2(1:5,3) = tmp(1:5,2) * a(2,4) + tmp(1:5,4) * a(4,2) 
              !
              vf2(1:3,1:3) = tmp2(2:4,1:3)
              dl2(    1:3) = tmp2(  5,1:3)
              ! 
              dmin = min(dmin, minval(dl2(1:3)))
              dmax = max(dmax, maxval(dl2(1:3)))
              !
              call tetrahedra2(vf2,dl2,usonic0)
              usonic(1:2,1:3) = usonic(1:2,1:3) + V * usonic0(1:2,1:3)
              !
           else if(e(3) < 0d0 .and. 0d0 < e(4)) then
              !
              ! C
              !
              !V = 3d0 * a(1,4) * a(2,4) * a(3,4) / (e(4) - 0d0)
              V = 3d0 * a(1,4) * a(2,4)           / (e(4) - e(3))
              !
              tmp2(1:5,1) = tmp(1:5,1) * a(1,4) + tmp(1:5,4) * a(4,1) 
              tmp2(1:5,2) = tmp(1:5,2) * a(2,4) + tmp(1:5,4) * a(4,2) 
              tmp2(1:5,3) = tmp(1:5,3) * a(3,4) + tmp(1:5,4) * a(4,3) 
              !
              vf2(1:3,1:3) = tmp2(2:4,1:3)
              dl2(    1:3) = tmp2(  5,1:3)
              ! 
              dmin = min(dmin, minval(dl2(1:3)))
              dmax = max(dmax, maxval(dl2(1:3)))
              !
              call tetrahedra2(vf2,dl2,usonic0)
              usonic(1:2,1:3) = usonic(1:2,1:3) + V * usonic0(1:2,1:3)
              !
           end if
           !
        end do ! ib
        !
     end do ! it
     !
  end do ! ik
  !$OMP END DO
  !
  !$OMP END PARALLEL
  !
  call MPI_allREDUCE(MPI_IN_PLACE, usonic, 6, &
  &                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  !
  if(my_rank == 0) then
     write(*,*) "Min & Max : ", dmin / beta, dmax / beta
     write(*,*) ""
     write(*,*) "Ultrasonic attenuation coefficient"
     write(*,*) "q, a_N, a_S, a_S / a_N"
     write(*,*) "100 : ", usonic(1,1), usonic(2,1), usonic(2,1) / usonic(1,1)
     write(*,*) "110 : ", usonic(1,2), usonic(2,2), usonic(2,2) / usonic(1,2)
     write(*,*) "001 : ", usonic(1,3), usonic(2,3), usonic(2,3) / usonic(1,3)
     write(*,*) ""
  end if
  !
end subroutine calc_usonic
!
! 2nd step of tetrahedra method.
!
subroutine tetrahedra2(vf1,dl1,usonic)
  !
  !
  real(8),intent(in) :: vf1(3,3), dl1(3)
  real(8),intent(out) :: usonic(2,3)
  !
  integer :: id, ii
  real(8) :: tmp(2,3), e(3), a(3,3), dl2(3), dl3(2), V, usonic0
  !
  usonic(1:2,1:3) = 0d0
  !
  do id = 1, 3
     !
     if(maxval(abs(vf1(id,1:3))) < 1d-10) then
        !$OMP CRITICAL   
        write(*,*) vf1(id,1:3)
        stop "vf = 0 !!"
        !$OMP END CRITICAL
     end if
     !
     tmp(1, 1:3) = vf1(id, 1:3)
     tmp(2, 1:3) = dl1(    1:3)
     !
     call sort(2, 3, tmp)
     !
     e(  1:3) = tmp(1,1:3)
     dl2(1:3) = tmp(2,1:3)
     !
     do ii = 1, 3
        a(ii,1:3) = (0d0 - e(1:3)) / (e(ii) - e(1:3))
     end do
     !
     if((e(1) < 0d0 .and. 0d0 <= e(2)) .or. (e(1) <= 0d0 .and. 0d0 < e(2))) then
        !
        !V = 2d0 * a(2,1) * a(3,1) / (0d0 - e(1)) 
        V = 2d0 * a(2,1)           / (e(3) - e(1)) 
        !
        dl3(1) = dl2(1) * a(1,2) + dl2(2) * a(2,1)
        dl3(2) = dl2(1) * a(1,3) + dl2(3) * a(3,1)
        !
        if(abs(dl3(1) - dl3(2)) < 1d-10) then
           usonic0 = 0.5d0 * (1d0 - tanh(0.25d0 * abs(sum(dl3(1:2)))))
        else if(dl3(1) * dl3(2) >= 0d0) then
           usonic0 = (1d0 + tanh(0.5d0 * abs(dl3(1)))) / (1d0 + tanh(0.5d0 * abs(dl3(2))))
           usonic0 = log(usonic0) / (abs(dl3(1)) - abs(dl3(2)))
        else
           usonic0 = (1d0 + tanh(0.5d0 * abs(dl3(1)))) * (1d0 + tanh(0.5d0 * abs(dl3(2))))
           usonic0 = log(usonic0) / (abs(dl3(1)) + abs(dl3(2)))
        end if
        !
        usonic(1,id) = V
        usonic(2,id) = V * 2d0 * usonic0
        !
     else if((e(2) < 0d0 .and. 0d0 <= e(3)) .or. (e(2) <= 0d0 .and. 0d0 < e(3))) then
        !
        !V = 2d0 * a(1,3) * a(2,3) / (e(3) - 0d0)
        V = 2d0 * a(1,3)           / (e(3) - e(2))
        !
        dl3(1) = dl2(3) * a(3,1) + dl2(1) * a(1,3)
        dl3(2) = dl2(3) * a(3,2) + dl2(2) * a(2,3)
        !
        if(abs(dl3(1) - dl3(2)) < 1d-10) then
           usonic0 = 0.5d0 * (1d0 - tanh(0.25d0 * abs(sum(dl3(1:2)))))
        else if(dl3(1) * dl3(2) >= 0d0) then
           usonic0 = (1d0 + tanh(0.5d0 * abs(dl3(1)))) / (1d0 + tanh(0.5d0 * abs(dl3(2))))
           usonic0 = log(usonic0) / (abs(dl3(1)) - abs(dl3(2)))
        else
           usonic0 = (1d0 + tanh(0.5d0 * abs(dl3(1)))) * (1d0 + tanh(0.5d0 * abs(dl3(2))))
           usonic0 = log(usonic0) / (abs(dl3(1)) + abs(dl3(2)))
        end if
        !
        usonic(1,id) = V
        usonic(2,id) = V * 2d0 * usonic0
        !
     end if
     !
  end do ! id
  !
end subroutine tetrahedra2
!
! Simple sort
!
subroutine sort(n1,n2,a)
  !
  !
  integer,intent(in) :: n1, n2
  real(8),intent(inout) :: a(n1,n2) 
  !
  integer :: i, m
  real(8) :: am, atmp(n1)
  !
  do i = 1, n2 - 1
     am = minval(a(1,i+1:n2) )
     m  = minloc(a(1,i+1:n2),1) + i
     if(a(1,i) .gt. am) then
        atmp(1:n1) = a(1:n1, m)
        a(1:n1, m) = a(1:n1,i)
        a(1:n1, i) = atmp(1:n1)
     end if
  end do
  !
end subroutine sort
!
end module ultrasonic_routines
!
! Main routine
!
program ultrasonic
  !
  use mpi
  use omp_lib
  use ultrasonic_vals, only : petot, my_rank, Z, dltf
  use ultrasonic_routines, only : read_stdin, read_dat, read_elph, read_Coulomb, read_delta, make_Z, &
  &                               gapeq, average_matrix, tetra_type, calc_fvel, calc_usonic
  implicit none
  !
  integer :: date(8), ierr
  real(8) :: t0, t1, t2
  !
  ! MPI Initialize
  !
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, PETOT, ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, my_rank, ierr)
  !
  t0 = OMP_GET_WTIME()
  call date_and_time(values = date)
  if(my_rank == 0) then
     write(*,*) ""
     write(*,'(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') &
     &     "   Start time : ", &
     &     date(1), "/", date(2), "/", date(3), " ", &
     &     date(5), ":", date(6), ":", date(7)
     write(*,*) ""
  end if
  !
  if(my_rank == 0) write(*,*) "  # of PEs : ", petot
  !$OMP PARALLEL
  if(my_rank == 0 .and. OMP_GET_THREAD_NUM() == 0) &
  write(*,*) '  # of thread : ', OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  Read from STDIN  #####"
     write(*,*)
  end if
  t1 = OMP_GET_WTIME()
  call read_stdin()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,*) "  Time : ", t2 - t1, " sec"
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  Read from data-file.xml  #####"
     write(*,*)
  end if
  t1 = OMP_GET_WTIME()
  call read_dat()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,*) "  Time : ", t2 - t1, " sec"
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  Read from elph*.dat  #####"
     write(*,*)
  end if
  t1 = OMP_GET_WTIME()
  call read_elph()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,*) "  Time : ", t2 - t1, " sec"
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  Read from vc*.dat  #####"
     write(*,*)
  end if
  t1 = OMP_GET_WTIME()
  call read_Coulomb()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,*) "  Time : ", t2 - t1, " sec"
  !
  ! 
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  Average matrix in grid  #####"
     write(*,*)
  end if
  t1 = OMP_GET_WTIME()
  call average_matrix()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,*) "  Time : ", t2 - t1, " sec"
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  Set or read initial delta  #####"
     write(*,*)
  end if
  t1 = OMP_GET_WTIME()
  call read_delta()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,*) "  Time : ", t2 - t1, " sec"
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  Compute renormalization fuctor : Z  #####"
     write(*,*)
  end if
  t1 = OMP_GET_WTIME()
  call make_Z()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,*) "  Time : ", t2 - t1, " sec"
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  Integrate gapeq  #####"
     write(*,*)
  end if
  t1 = OMP_GET_WTIME()
  call gapeq()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,*) "  Time : ", t2 - t1, " sec"
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  Write delta.fs  #####"
     write(*,*)
  end if
  t1 = OMP_GET_WTIME()
  !
  call tetra_type()
  call calc_fvel()
  call calc_usonic()
  !
  t2 = OMP_GET_WTIME()
  !
  if(my_rank == 0) then
     write(*,*) ""
     write(*,*) "  Total time : ", t2 - t0, "sec"
     write(*,*) ""
     write(*,*) ''
     write(*,*) '###############  Done  ###############'  
     write(*,*) ''
  end if
  !
  call MPI_FINALIZE(ierr)
  !
end program ultrasonic


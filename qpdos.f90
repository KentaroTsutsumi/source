module qpdos_vals
  !
  implicit none
  !
  real(8),parameter :: &
  & pi = acos(-1d0) ! Pi = 3.14159265358979...
  !
  integer,save :: &
  & ncf,           & ! # of points for Chebyshev interpolation
  & nbee,          & ! # of Bands for Vc
  & fbee,          & ! First band for Vc
  & lbee,          & ! Last band for Vc
  & nbep,          & ! # of bands for elph
  & fbep,          & ! First band for elph
  & lbep,          & ! Last band for elph
  & fstk,          & ! First k of this PE
  & lstk,          & ! Last k of this PE
  & ivvec(3,20,6), & !
  & my_rank,       & ! Index of PE
  & petot,         & ! # of PEs
  & ne,            & !
  & nx,            & ! = 2 * ne - 1
  & nmf,           & ! # of Matsubara frequencoes
  & nt,            & ! Total # of Delta, Xi, ...
  & nb,            & ! # of bands
  & nbf,           & ! # of bands contais FS
  & fstb,          & ! first band contais FS
  & lstb,          & ! last band contais FS
  & ng(3),         & ! k grid for matrix elements
  & nk,            & ! ng(1) * ng(2) * ng(3). # of total k points for matrix elements
  & nk0,           & ! # of irreducible k points
  & ngd(3),        & ! k grid for DOS
  & nkd,           & ! ngd(1) * ngd(2) * ngd(3). # of total k points for DOS
  & nm,            & ! # of modes
  & nsym             ! # of symmetries
  !
  real(8),save :: &
  & emax,       & !
  & xic,        & !
  & wlsm(4,20), &
  & bvec(3,3),  & ! Reciplocal lattice vector
  & beta          ! inversed temperature [Ry]
  !
  integer,allocatable,save :: &
  & kindx(:),    & ! (nt) k point for gap equation
  & bindx(:),    & ! (nt) band index for gap equation
  & dgrid(:,:),  & ! (3,nkd) k points
  & sym(:,:,:)     ! (3,3,nsym) Symmetrical operators
  !
  real(8),allocatable,save :: &
  & rkv(:,:),       & ! (3,nk0) irreducible k points
  & grid(:,:),      & ! (3,nk) k points
  & dosk(:,:,:),    & ! (nx,nb,nkd)
  & xi0(:),         & ! (nx)
  & e0(:),          & ! (ne)
  & sdos(:),        & ! (ne)
  & mf(:),          & ! (nmf) Matsubara frequencies
  & wmf(:),         & ! (nmf) Weight for frequency integration
  & Z(:,:,:),       & ! (nx,nbf,nkd) Renormalization fuctor
  & omg0(:,:),      & ! (nm,nk0) Phonon frequencies [Ry]
  & omgf(:,:,:),    & ! (nm,nk,nkd) Phonon frequencies [Ry]
  & gg0(:,:,:,:,:), & ! (nm,nb,nb,nk,nk0) El-Ph matrix element [Ry]
  & Vc0(:,:,:,:,:), & ! (0:nmf,nb,nb,nk,nk0) Screened Coulomb matrix element [Ry]
  & ggf(:,:,:,:,:), & ! (nm,nb,nk,nbf,nkd) El-Ph matrix element [Ry]
  & Vcf(:,:,:,:,:), & ! (0:nmf,nb,nk,nbf,nkd) Screened Coulomb matrix element [Ry]
  & delta(:),       & ! (nt) Kohn-Sham gap functions [Ry]
  & dltf(:,:,:),    & ! (nx,nbf,nkd) Kohn-Sham gap functions [Ry]
  & xi(:),          & ! (nt) Kohn-Sham energy [Ry]
  & dk(:),          & ! (nt) Weight of k
  & eig(:,:)          ! (nb,nkd) Kohn-Sham energy [Ry]
  !
end module qpdos_vals
!
module qpdos_routines
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
  use qpdos_vals, only : my_rank, beta, fbee, lbee, nbee, xic, nmf, ne, emax
  !
  integer :: ierr
  !
  namelist /input/ beta, xic, fbee, lbee, xic, nmf, ne, emax
  !
  xic = -1d0
  !
  if(my_rank == 0) then
     !
     read(*,input,err=100)
     !
     write(*,*) '              Temparature[K] : ', beta
     beta = 157887d0 / beta
     write(*,*) '     Inverse temparature[/Ry] : ', beta
     write(*,*) '                    Xi cutoff : ', xic
     write(*,*) '                   First band : ', fbee
     write(*,*) '                    Last band : ', lbee
     write(*,*) "   # of Matsubara frequencies : ", nmf
     write(*,*) "       # of energy for QP-DOS : ", ne
     write(*,*) "  Max energy for QP-DOS [meV] : ", emax
     !
  end if
  !
  call MPI_BCAST(beta, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(xic,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(fbee, 1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(lbee, 1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  nbee = lbee - fbee + 1
  call MPI_BCAST(nmf,  1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ne,   1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(emax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  !
  return
  !
100 write(*,*) "stop. reading namelist file"
  !
  write(*,*) '               Temparature[K] : ', beta
  beta = 157887d0 / beta
  write(*,*) '     Inverse temparature[/Ry] : ', beta
  write(*,*) '                    Xi cutoff : ', xic
  write(*,*) '                   First band : ', fbee
  write(*,*) '                    Last band : ', lbee
  write(*,*) "   # of Matsubara frequencies : ", nmf     
  write(*,*) "       # of energy for QP-DOS : ", ne
  write(*,*) "  Max energy for QP-DOS [meV] : ", emax
  !
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
  use qpdos_vals, only : my_rank, nsym, sym, nb, ngd, nkd, eig, bvec, nbf, &
  &                  fstb, lstb, fstk, lstk
  !use qpdos_routines, only : cnt_and_dsp
  !
  integer :: isym, fi = 10, ib
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
     allocate(eig(nb, nkd))
     open(fi, file = 'eigval.dat')
     write(*,*) "  read from eigval.dat"
     read(fi,*) eig(1:nb,1:nkd)
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
  if(my_rank /= 0) allocate(eig(nb, nkd))
  call MPI_BCAST(eig, nb * nkd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, isym)
  call MPI_BCAST(bvec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, isym)
  !
  ! Which band contains Fermi level ?
  !
  do ib = 1, nb
     !
     if(minval(eig(ib,1:nkd)) < 0d0 .and. &
     &  maxval(eig(ib,1:nkd)) > 0d0) then
        fstb = ib
        exit
     end if
     !
  end do
  !
  do ib = fstb, nb
     !
     if(minval(eig(ib,1:nkd)) < 0d0 .and. &
     &  maxval(eig(ib,1:nkd)) > 0d0) then
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
  use qpdos_vals, only : my_rank, ng, nk, nk0, nm, gg0, omg0, rkv, nbep, fbep, lbep
  !use qpdos_routines, only : cnt_and_dsp
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
  use qpdos_vals, only : petot, my_rank
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
  use qpdos_vals, only : my_rank, nk, nk0, ng, rkv, grid, nsym, sym
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
  use qpdos_vals, only : my_rank, nk0, nk, nb, Vc0, rkv, ng, ncf, nmf
  !use qpdos_routines, only : cnt_and_dsp
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
  use qpdos_vals, only : petot, my_rank, nk0, nk, nkd, ngd, fstk, lstk, &
  &                  fstb, lstb, nbf, nbep, fbep, lbep, nm, &
  &                  ncf, Vc0, gg0, omg0, Vcf, ggf, omgf, nbee, fbee, lbee
  !use qpdos_routines, only : cnt_and_dsp, interpol_indx_weight, &
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
  use qpdos_vals, only : nk, nk0, nsym, sym, rkv, grid, ng, my_rank
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
  use qpdos_vals, only : nk, nk0, fstk, lstk, nsym, grid, ng, sym, nsym, my_rank, &
  &                  nkd, ngd
  !use qpdos_routines, only : interpol_indx_weight
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
  use qpdos_vals, only : my_rank, nt, xi, delta, dk, kindx, bindx
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
!
!
subroutine egrid()
  !
  use mpi
  use qpdos_vals, only : my_rank, ne, nx, e0, sdos, xi0, emax
  !
  integer :: ie, ierr
  real(8) :: de
  !
  emax = emax / 13.6057d3
  nx = ne * 2 - 1
  !
  allocate(e0(ne), sdos(ne), xi0(nx))
  de = emax / dble(ne)
  !
  do ie = 1, ne
     e0(ie) = de * dble(ie)
  end do
  !
  do ie = 1, nx
     xi0(ie) = de * dble(ie - ne)
  end do
  !
  sdos(1:ne) = 0d0
  !
end subroutine egrid
!
! Calc Z_{n k}
!
subroutine make_Z()
  !
  use mpi
  use qpdos_vals, only : nx, nkd, nbf, ngd, nm, nt, xi, dk, kindx, bindx, omgf, ggf, Z, xi0, &
  &                  fstk, lstk, fbep, lbep, beta
  !use qpdos_routines, only : Zweight
  !
  integer :: ik, ib, jt, jk, jb, im, ierr, ix
  real(8) :: x, xp, om, tx, txp, tom
  !
  allocate(Z(nx,nbf,nkd))
  Z(1:nx,1:nbf,1:nkd) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nx,fstk,lstk,nbf,nm,nt,xi,dk,kindx,bindx,Z,ggf,omgf,xi0,fbep,lbep,beta) &
  !$OMP & PRIVATE(ik, ib, jt, jk, jb, ix, im, x, xp, om, tx, txp, tom)
  !
  !$OMP DO
  do ik = fstk, lstk
     !
     do ib = 1, nbf
        !
        do ix = 1, nx
           !
           x = abs(0.5d0 * beta * xi0(ix))
           tx = tanh(x)
           !
           do jt = 1, nt
              !
              xp = abs(xi(jt) * 0.5d0 * beta)
              txp = tanh(xp)
              jk = kindx(jt)
              jb = bindx(jt)
              !
              if(jb < fbep .or. lbep < jb) cycle
              !
              do im = 1, nm
                 !
                 om = abs(omgf(im,jk,ik) * 0.5d0 * beta)
                 tom = tanh(om)
                 !    
                 Z(ix,ib,ik) = Z(ix,ib,ik) &
                 &        + dk(jt) * ggf(im,jb,jk,ib,ik) &
                 &        * beta**2 * Zweight(x, xp, om, tx, txp, tom)
                 !
              end do ! im
              !
           enddo ! jt
           !
           Z(ix,ib,ik) = - Z(ix,ib,ik)
           !
        end do
        !
     end do ! ib
     !
  end do ! ik
  !$OMP END DO
  !
  !$OMP END PARALLEL
  !
  call MPI_allREDUCE(MPI_IN_PLACE, Z, nx * nbf * nkd, &
  &                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  !
end subroutine make_Z
!
! Weight for Z
!
function Zweight(x,y,z,tx,ty,tz) result(Wz)
  !
  !
  real(8),intent(in) :: x, y, z, tx, ty, tz
  real(8) :: Wz
  !
  real(8) :: thr = 1d-8, txmz, txpz, mm, mp, pm, pp
  !
  if(abs(x) < thr) then
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
  else
     !
     mm = 1d0 / (x - y - z)
     mp = 1d0 / (x - y + z)
     pm = 1d0 / (x + y - z)
     pp = 1d0 / (x + y + z)
     !
     txpz = tanh(x + z)
     txmz = tanh(x - z)
     !
     if(abs(x - y - z) < thr) then
        !
        Wz = 0.0625d0 * ( &
        &   (1d0 + ty) * (tx - 1d0) * (1d0 + tz) * ty / (tz * tx) &
        & - (1d0 - 1d0 / (- tz * tx)) * (-1d0 + txpz**2 + (- ty + txpz) * mp) * mp &
        & + (1d0 - 1d0 / (  tz * tx)) * (-1d0 + txmz**2 + (  ty + txmz) * pm) * pm &
        & - (1d0 - 1d0 / (- tz * tx)) * (-1d0 + txpz**2 + (  ty + txpz) * pp) * pp &
        & )
        !
     else if(abs(x - y + z) < thr) then
        !
        Wz = 0.0625d0 * ( &
        &   (1d0 - 1d0 / (  tz * tx)) * (-1d0 + txmz**2 + (- ty + txmz) * mm) * mm &
        & - (1d0 + ty) * (tx - 1d0) * (1d0 - tz) * ty / (- tz * tx) &
        & + (1d0 - 1d0 / (  tz * tx)) * (-1d0 + txmz**2 + (  ty + txmz) * pm) * pm &
        & - (1d0 - 1d0 / (- tz * tx)) * (-1d0 + txpz**2 + (  ty + txpz) * pp) * pp &
        & )
        !
     else if(abs(x + y - z) < thr) then
        !
        Wz = 0.0625d0 * ( &
        &   (1d0 - 1d0 / (  tz * tx)) * (-1d0 + txmz**2 + (- ty + txmz) * mm) * mm &
        & - (1d0 - 1d0 / (- tz * tx)) * (-1d0 + txpz**2 + (- ty + txpz) * mp) * mp &
        & + (1d0 - ty) * (tx - 1d0) * (1d0 + tz) * (- ty) / (tz * tx) &
        & - (1d0 - 1d0 / (- tz * tx)) * (-1d0 + txpz**2 + (  ty + txpz) * pp) * pp &
        & )
        !
     else if(abs(x + y + z) < thr) then
        !
        Wz = 0.0625d0 * ( &
        &   (1d0 - 1d0 / (  tz * tx)) * (-1d0 + txmz**2 + (- ty + txmz) * mm) * mm &
        & - (1d0 - 1d0 / (- tz * tx)) * (-1d0 + txpz**2 + (- ty + txpz) * mp) * mp &
        & + (1d0 - 1d0 / (  tz * tx)) * (-1d0 + txmz**2 + (  ty + txmz) * pm) * pm &
        & - (1d0 - ty) * (tx - 1d0) * (1d0 - tz) * (- ty) / (- tz * tx) &
        & )
        !
     else 
        !
        Wz = 0.0625d0 * ( &
        &   (1d0 - 1d0 / (  tz * tx)) * (-1d0 + txmz**2 + (- ty + txmz) * mm) * mm &
        & - (1d0 - 1d0 / (- tz * tx)) * (-1d0 + txpz**2 + (- ty + txpz) * mp) * mp &
        & + (1d0 - 1d0 / (  tz * tx)) * (-1d0 + txmz**2 + (  ty + txmz) * pm) * pm &
        & - (1d0 - 1d0 / (- tz * tx)) * (-1d0 + txpz**2 + (  ty + txpz) * pp) * pp &
        & )
        !
     end if
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
  use qpdos_vals, only : pi, my_rank, nx, nkd, nt, ngd, nb, nbf, nm, ncf, omgf, ggf, Vcf, &
  &                  xi, dk, delta, kindx, bindx, dltf, Z, beta, xi0, xic, nmf, mf, wmf, &
  &                  fstk, lstk, fbep, lbep
  !use qpdos_routines, only : Kweight, calc_Kel, weightspoints_gl
  !
  integer :: ik, ib, jt, jk, jb, im, ierr, ix, nh
  real(8) :: x, xp, om, bx, bxp, tx, txp, tom, eqp, Kph, Kel, Kave, dlth, dosh, xmax
  !
  allocate(dltf(nx,nbf,nkd))
  dltf(1:nx,1:nbf,1:nkd) = 0d0
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
  if(my_rank == 0) write(*,*) nh, dlth, xmax, dosh
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(xi0, nx, fstk, lstk, nt, nb, nbf, nm, ncf, omgf, ggf, Vcf, &
  !$OMP &        xi, dk, delta, kindx, bindx, dltf, Z, beta, fbep, lbep, &
  !$OMP &        nh, dlth, dosh, xic, xmax) &
  !$OMP & PRIVATE(ik, ib, jt, jk, jb, im, ix, &
  !$OMP &         x, xp, om, bx, bxp, tx, txp, tom, eqp, Kph, Kel, Kave)
  !
  !$OMP DO
  do ik = fstk, lstk
     !
     do ib = 1, nbf
        !
        do ix = 1, nx 
           !
           Kave = 0d0
           x = abs(xi0(ix))
           bx = abs(0.5d0 * beta * x)
           tx = tanh(bx)
           !
           do jt = 1, nt
              !
              xp = abs(xi(jt))
              bxp = abs(0.5d0 * beta * xp)
              txp = tanh(bxp)
              jk = kindx(jt)
              jb = bindx(jt)
              !
              Kph = 0d0
              !
              if(fbep <= jb .and. jb <= lbep) then
                 !
                 do im = 1, nm
                    !
                    om = abs(omgf(im,jk,ik) * 0.5d0 * beta)
                    tom = tanh(om)
                    !
                    Kph = Kph + ggf(im,jb,jk,ib,ik) &
                    &   * beta * Kweight(bx, bxp, om, tx, txp, tom)
                    !
                 end do ! im
                 !               
                 Kph = Kph * 2d0
                 !
              end if ! (fbep <= jb .and. jb <= lbep)
              !
              Kel = calc_Kel(x,xp,Vcf(1:ncf,jb,jk,ib,ik))
              !
              eqp = sqrt(delta(jt)**2 + xp**2)
              !
              dltf(ix,ib,ik) = dltf(ix,ib,ik) &
              &                 - 0.5d0 * dk(jt) * (Kph + Kel) / (1d0 + Z(ix,ib,ik)) &
              &                   * delta(jt) / eqp * tanh(0.5d0 * beta * eqp)
              !
              if(xi(jt) > xic) Kave = Kave + (Kph + Kel)
              !
           enddo ! jt
           !
           if(xic > 0) then
              !
              Kave = Kave / dble(nh)
              dltf(ix,ib,ik) = dltf(ix,ib,ik) - 0.5d0 * Kave / (1d0 + Z(ix,ib,ik)) &
              &                          * dosh * dlth / xmax
              !
           end if ! (xic > 0)
           !
           dltf(ix,ib,ik) = sqrt(xi0(ix)**2 + dltf(ix,ib,ik)**2)
           !
        end do ! ix
        !
     end do ! ib
     !
  end do ! ik
  !$OMP END DO
  !
  !$OMP END PARALLEL
  ! 
  call MPI_allREDUCE(MPI_IN_PLACE, dltf, nx * nbf * nkd, &
  &                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  !
  deallocate(ggf, Vcf, omgf)
  !
end subroutine gapeq
!
! Weights & Points for Gauss-Legendre method
!
subroutine weightspoints_gl(n,x,w)
  !
  use mpi
  use qpdos_vals, only: pi
  !use qpdos_routines, only : legendre
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
function Kweight(x,y,z,tx,ty,tz) result(Wk)
  !
  !
  real(8),intent(in) :: x, y, z, tx, ty, tz
  real(8) :: Wk
  !
  real(8) :: thr = 1d-1, zp, zm
  !
  if(abs(x) < thr) then
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
  else if(abs(y) < thr) then
     !
     zp = 1d0 / (  x + z)
     zm = 1d0 / (- x + z)
     !
     if(abs(x + z) < thr) then
        !
        Wk = ( &
        &      (-1d0 + tx * tz + (- tx + tz) * zm) * (- zm) &
        &    ) / (4d0 * tx * tz)
        !
     else if(abs(- x + z) < thr) then
        !
        Wk = ( &
        &        (-1d0 - tx * tz + (  tx + tz) * zp) * (  zp) &
        &    ) / (4d0 * tx * tz)
        !
     else
        !
        Wk = ( &
        &        (-1d0 - tx * tz + (  tx + tz) * zp) * (  zp) &
        &      + (-1d0 + tx * tz + (- tx + tz) * zm) * (- zm) &
        &    ) / (4d0 * tx * tz)
        !
     end if
     !
  else
     !
     if(abs(x - y - z) < thr) then
        !
        Wk =  ( &
        &         (  (1d0 - tx) * (  tz + 1d0) * (1d0 + ty)) &
        &       - (- tx * ty * tz + tx - ty + tz) / (- (x - y + z)) &
        &       + (- tx * ty * tz + tx + ty - tz) / (- (x + y - z)) &
        &       - (  tx * ty * tz + tx + ty + tz) / (  (x + y + z)) &
        &     ) / (8d0 * tx * ty * tz)
        !
     else if(abs(x - y + z) < thr) then
        !
        Wk =  ( &
        &         (  tx * ty * tz + tx - ty - tz) / (  (x - y - z)) &
        &       - (- (1d0 - tx) * (- tz + 1d0) * (1d0 + ty)) &
        &       + (- tx * ty * tz + tx + ty - tz) / (- (x + y - z)) &
        &       - (  tx * ty * tz + tx + ty + tz) / (  (x + y + z)) &
        &     ) / (8d0 * tx * ty * tz)
        !
     else if(abs(x + y - z) < thr) then
        !
        Wk =  ( &
        &         (  tx * ty * tz + tx - ty - tz) / (  (x - y - z)) &
        &       - (- tx * ty * tz + tx - ty + tz) / (- (x - y + z)) &
        &       + (- (1d0 - tx) * (  tz + 1d0) * (1d0 - ty)) &
        &       - (  tx * ty * tz + tx + ty + tz) / (  (x + y + z)) &
        &     ) / (8d0 * tx * ty * tz)
        !
     else if(abs(x + y + z) < thr) then
        !
        Wk =  ( &
        &         (  tx * ty * tz + tx - ty - tz) / (  (x - y - z)) &
        &       - (- tx * ty * tz + tx - ty + tz) / (- (x - y + z)) &
        &       + (- tx * ty * tz + tx + ty - tz) / (- (x + y - z)) &
        &       - (  (1d0 - tx) * (- tz + 1d0) * (1d0 - ty)) &
        &     ) / (8d0 * tx * ty * tz)
        !
     else
        !
        Wk =  ( &
        &         (  tx * ty * tz + tx - ty - tz) / (  (x - y - z)) &
        &       - (- tx * ty * tz + tx - ty + tz) / (- (x - y + z)) &
        &       + (- tx * ty * tz + tx + ty - tz) / (- (x + y - z)) &
        &       - (  tx * ty * tz + tx + ty + tz) / (  (x + y + z)) &
        &     ) / (8d0 * tx * ty * tz)
        !
     end if
     !
  end if
  !
end function Kweight
!
! Compute Kel
!
function calc_Kel(x,xp,Vci) result(Kel)
  !
  use qpdos_vals, only : pi, ncf, nmf, mf, wmf
  !
  real(8),intent(in) :: x, xp, Vci(ncf)
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
  xpxp = abs(x) + abs(xp)
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
! define shortest diagonal line & define type of tetragonal
!
subroutine tetra_type()
  !
  use qpdos_vals, only : ngd, bvec, ivvec, wlsm
  !
  integer :: itype, i1, i2, i3, itet, divvec(4,4), ivvec0(4)
  real(8) :: l(4), bvec2(3,3), bvec3(3,4)
  !
  do i1 = 1, 3
     bvec2(1:3,i1) = bvec(1:3,i1) / dble(ngd(i1))
  end do
  !
  bvec3(1:3,1) = -bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
  bvec3(1:3,2) =  bvec2(1:3,1) - bvec2(1:3,2) + bvec2(1:3,3)
  bvec3(1:3,3) =  bvec2(1:3,1) + bvec2(1:3,2) - bvec2(1:3,3)
  bvec3(1:3,4) =  bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
  !
  ! length of delta bvec
  !
  do i1 = 1, 4
     l(i1) = dot_product(bvec3(1:3,i1),bvec3(1:3,i1))
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
  itet = 0
  do i1 = 1, 3
     do i2 = 1, 3
        if(i2 == i1) cycle
        do i3 = 1, 3
           if(i3 == i1 .or. i3 == i2) cycle
           !
           itet = itet + 1
           !
           ivvec(1:3,1,itet) = ivvec0(1:3)
           ivvec(1:3,2,itet) = ivvec(1:3,1,itet) + divvec(1:3,i1)
           ivvec(1:3,3,itet) = ivvec(1:3,2,itet) + divvec(1:3,i2)
           ivvec(1:3,4,itet) = ivvec(1:3,3,itet) + divvec(1:3,i3)
           !
        end do
     end do
  end do
  !
  ivvec(1:3, 5,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6)
  ivvec(1:3, 6,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,3,1:6)
  ivvec(1:3, 7,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,4,1:6)
  !
  ivvec(1:3, 8,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,1,1:6)
  ivvec(1:3, 9,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6)
  ivvec(1:3,10,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,4,1:6)
  !
  ivvec(1:3,11,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,1,1:6)
  ivvec(1:3,12,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6)
  ivvec(1:3,13,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,2,1:6)
  !
  ivvec(1:3,14,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6)
  ivvec(1:3,15,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,2,1:6)
  ivvec(1:3,16,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,3,1:6)
  !
  ivvec(1:3,17,1:6) = -ivvec(1:3,1,1:6) + ivvec(1:3,2,1:6) + ivvec(1:3,4,1:6)
  ivvec(1:3,18,1:6) =  ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6) + ivvec(1:3,3,1:6)
  ivvec(1:3,19,1:6) =  ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6) + ivvec(1:3,4,1:6)
  ivvec(1:3,20,1:6) =  ivvec(1:3,1,1:6) + ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6)
  !
  wlsm(1, 1:10) = dble((/1440,    0,   30,    0, -38, -56, -38, -28,  7,     9/))
  wlsm(2, 1:10) = dble((/   0, 1440,    0,   30, -28,   9,   7, -38, -38,  -56/))
  wlsm(3, 1:10) = dble((/  30,    0, 1440,    0,  17, -46,  17,   7, -28,    9/))
  wlsm(4, 1:10) = dble((/   0,   30,    0, 1440,   7,   9, -28,  17,  17,  -46/))
  !
  wlsm(1,11:20) = dble((/ -46,   17,   17,  -28,   9,   7, -18, -18,  12,  -18/))
  wlsm(2,11:20) = dble((/   9,    7,  -28,   17, -46,  17, -18, -18, -18,   12/))
  wlsm(3,11:20) = dble((/ -56,  -38,  -38,    7,   9, -28,  12, -18, -18,  -18/))
  wlsm(4,11:20) = dble((/   9,  -28,    7,  -38, -56, -38, -18,  12, -18,  -18/))
  !
  wlsm(1:4,1:20) = wlsm(1:4,1:20) / 1260
  !
  wlsm(1:4,1:20) = 0d0
  wlsm(1,1:4) = dble((/1, 0, 0, 0/))
  wlsm(2,1:4) = dble((/0, 1, 0, 0/))
  wlsm(3,1:4) = dble((/0, 0, 1, 0/))
  wlsm(4,1:4) = dble((/0, 0, 0, 1/))
  !
end subroutine tetra_type
!
! Compute Dos_k, xmin, xmax
! 
subroutine calc_dosk()
  !
  use mpi
  use qpdos_vals, only : ngd, nkd, nbf, fstb, ivvec, dgrid, nx, dosk, dgrid, &
  &                  xi0, eig, wlsm, fstk, lstk
  !use qpdos_routines, only : sort
  !
  integer :: ik, ik2, ib, it, ii, jj, ix, ikv1(3), i1, i2, i3, ierr
  real(8) :: ei(4,nbf), e(4), a(4,4), V, tmp(5,4), w0(4,4), w1(4)
  !
  allocate(dgrid(3,nkd), dosk(nx,nbf,nkd))
  !
  w0(1:4,1:4) = 0d0
  do ii = 1, 4
     w0(ii,ii) = 1d0
  end do
  !
  ! Calc. Dos_k
  !
  dosk(1:nx, 1:nbf, 1:nkd) = 0d0
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
  dosk(1:nx,1:nbf,1:nkd) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fstk,lstk,nbf,fstb,nx,ivvec,w0,ngd,dgrid,wlsm,eig,xi0,dosk) &
  !$OMP & PRIVATE(ik,ik2,it,ii,jj,ib,ikv1,ix, &
  !$OMP &         ei,e,V,a,w1,tmp)
  !
  !$OMP DO REDUCTION(+: dosk)
  do ik = fstk, lstk
     !
     do it = 1, 6
        !
        ei(1:4,1:nbf) = 0d0
        !
        do ii = 1, 20
           !
           ikv1(1:3) = dgrid(1:3,ik) + ivvec(1:3, ii, it)
           ikv1(1:3) = modulo(ikv1(1:3), ngd(1:3))
           ik2 = 1 + ikv1(1) + ngd(1) * ikv1(2) + ngd(1) * ngd(2) * ikv1(3)
           !
           do ib = 1, nbf
              ei(1:4,ib) = ei(1:4,ib) &
              &          + wlsm(1:4,ii) * eig(fstb + ib - 1,ik2)
           end do
           !
        end do
        !
        do ib = 1, nbf
           !
           tmp(  1,1:4) = ei(    1:4,ib)
           tmp(2:5,1:4) = w0(1:4,1:4)
           !
           call sort(5,4,tmp)
           !
           e(1:4) = tmp(1,1:4)
           !
           do ix = 1, nx
              !
              w1(1:4) = 0d0
              !
              do ii = 1, 4
                 do jj = 1, 4
                    a(jj,ii) = (xi0(ix) - e(ii)) / (e(jj) - e(ii))
                 end do
              end do
              !
              if(e(1) <= xi0(ix) .and. xi0(ix) < e(2)) then
                 !
                 V = a(2,1) * a(3,1) * a(4,1) / (xi0(ix) - e(1))
                 !
                 w1(1:4) = w1(1:4) + V * (tmp(2:5,1) * a(1,2) + tmp(2:5,2) * a(2,1))
                 w1(1:4) = w1(1:4) + V * (tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1))
                 w1(1:4) = w1(1:4) + V * (tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1))
                 !
              else if(e(2) <= xi0(ix) .and. xi0(ix) < e(3)) then
                 !
                 V = a(3,1) * a(4,1) * a(2,4) / (xi0(ix) - e(1))
                 !
                 w1(1:4) = w1(1:4) + V * (tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1))
                 w1(1:4) = w1(1:4) + V * (tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1))
                 w1(1:4) = w1(1:4) + V * (tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2))
                 !
                 V = a(2,3) * a(3,1) * a(4,2) / (xi0(ix) - e(1))
                 !
                 w1(1:4) = w1(1:4) + V * (tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1))
                 w1(1:4) = w1(1:4) + V * (tmp(2:5,2) * a(2,3) + tmp(2:5,3) * a(3,2))
                 w1(1:4) = w1(1:4) + V * (tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2))
                 !
              else if(e(3) <= xi0(ix) .and. xi0(ix) < e(4)) then
                 !
                 V = a(1,4) * a(2,4) * a(3,4) / (e(4) - xi0(ix))
                 !
                 w1(1:4) = w1(1:4) + V * (tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1))
                 w1(1:4) = w1(1:4) + V * (tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2))
                 w1(1:4) = w1(1:4) + V * (tmp(2:5,3) * a(3,4) + tmp(2:5,4) * a(4,3))
                 !
              end if
              !
              do ii = 1, 20
                 !
                 ikv1(1:3) = dgrid(1:3,ik) + ivvec(1:3, ii, it)
                 ikv1(1:3) = modulo(ikv1(1:3), ngd(1:3))
                 ik2 = 1 + ikv1(1) + ngd(1) * ikv1(2) + ngd(1) * ngd(2) * ikv1(3)
                 !
                 dosk(ix,ib,ik2) = dosk(ix,ib,ik2) &
                 &              + sum(wlsm(1:4,ii) * w1(1:4))
                 !
              end do ! ii
              !
           end do ! ix
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
  dosk(1:nx, 1:nbf, 1:nkd) = dosk(1:nx, 1:nbf, 1:nkd) / dble(6 * nkd)
  !
  call MPI_allREDUCE(MPI_IN_PLACE, dosk, nx * nbf * nkd, &
  &                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  !
end subroutine calc_dosk
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
        atmp(:) = a(:,m)
        a(:,m) = a(:,i)
        a(:,i) = atmp(:)
     end if
  end do
  !
end subroutine sort
!
!
!
subroutine calc_sdos()
  !
  use mpi
  use qpdos_vals, only : nkd, nbf, nx, ne, ngd, dgrid, dltf, dosk, e0, sdos, fstk, lstk
  !use qpdos_routines, only : penta_4d, sort
  !
  integer :: ik, ib, ix, ie, ip, ii, jj, ik2, ix2, ikv1(3), iv4d(4,5,24), ierr
  real(8) :: e(5), dos2(5), tmp(2,5), a(5,5)
  !
  call penta_4d(iv4d)
  !
  sdos = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fstk,lstk,nkd,nbf,nx,ne,ngd,dgrid,iv4d,dltf,dosk,e0,sdos) &
  !$OMP & PRIVATE(ik,ib,ix,ip,ii,jj,ie,ik2,ix2,e,a,dos2,tmp,ikv1)
  !
  !$OMP DO REDUCTION(+: sdos)
  do ik = fstk, lstk
     do ib = 1, nbf
        do ix = 1, nx - 1
           !
           do ip = 1, 24
              !
              do ii = 1, 5
                 !
                 ikv1(1:3) = dgrid(1:3,ik) + iv4d(1:3,ii,ip)
                 ikv1(1:3) = modulo(ikv1(1:3), ngd(1:3))
                 ik2 = 1 + ikv1(1) + ngd(1) * ikv1(2) + ngd(1) * ngd(2) * ikv1(3)
                 !
                 ix2 = ix + iv4d(4,ii,ip)
                 !
                 e(   ii) = dltf(ix2,ib,ik2)
                 dos2(ii) = dosk(ix2,ib,ik2)
                 !
              end do
              !
              tmp(1,1:5) = e(   1:5)
              tmp(2,1:5) = dos2(1:5)
              call sort(2, 5, tmp)
              e(   1:5) = tmp(1,1:5)
              dos2(1:5) = tmp(2,1:5)
              !
              do ie = 1, ne
                 !
                 do ii = 1, 5
                    do jj = 1, 5
                       a(jj,ii) = (e0(ie) - e(ii)) / (e(jj) - e(ii))
                    end do
                 end do
                 !
                 if(e(1) < e0(ie) .and. e0(ie) <= e(2)) then
                    !
                    sdos(ie) = sdos(ie) &
                    &        + a(2,1) * a(3,1) * a(4,1) * a(5,1) / (e0(ie) - e(1)) &
                    &        * ( dos2(1) * a(1,2) + dos2(2) * a(2,1) & 
                    &          + dos2(1) * a(1,3) + dos2(3) * a(3,1) &
                    &          + dos2(1) * a(1,4) + dos2(4) * a(4,1) &
                    &          + dos2(1) * a(1,5) + dos2(5) * a(5,1) )
                    !
                 else if(e(2) < e0(ie) .and. e0(ie) <= e(3)) then
                    !
                    sdos(ie) = sdos(ie) &
                    &        + a(3,1) * a(4,1) * a(5,1) * a(2,3) / (e0(ie) - e(1)) &
                    &        * ( dos2(1) * a(1,3) + dos2(3) * a(3,1) &
                    &          + dos2(1) * a(1,4) + dos2(4) * a(4,1) &
                    &          + dos2(1) * a(1,5) + dos2(5) * a(5,1) &
                    &          + dos2(2) * a(2,3) + dos2(3) * a(3,2) )
                    !
                    sdos(ie) = sdos(ie) &
                    &        + a(4,1) * a(5,1) * a(3,2) * a(2,4) / (e0(ie) - e(1)) &
                    &        * ( dos2(1) * a(1,4) + dos2(4) * a(4,1) &
                    &          + dos2(1) * a(1,5) + dos2(5) * a(5,1) &
                    &          + dos2(2) * a(2,3) + dos2(3) * a(3,2) &
                    &          + dos2(2) * a(2,4) + dos2(4) * a(4,2) )
                    !
                    sdos(ie) = sdos(ie) &
                    &        + a(5,1) * a(3,2) * a(4,2) * a(2,5) / (e0(ie) - e(1)) &
                    &        * ( dos2(1) * a(1,5) + dos2(5) * a(5,1) &
                    &          + dos2(2) * a(2,3) + dos2(3) * a(3,2) &
                    &          + dos2(2) * a(2,4) + dos2(4) * a(4,2) & 
                    &          + dos2(2) * a(2,5) + dos2(5) * a(5,2) )
                    !
                else if(e(3) < e0(ie) .and. e0(ie) <= e(4)) then
                    !
                    sdos(ie) = sdos(ie) &
                    &        + a(3,5) * a(2,5) * a(1,5) * a(4,3) / (e(5) - e0(ie)) &
                    &        * ( dos2(5) * a(5,3) + dos2(3) * a(3,5) &
                    &          + dos2(5) * a(5,2) + dos2(2) * a(2,5) &
                    &          + dos2(5) * a(5,1) + dos2(1) * a(1,5) &
                    &          + dos2(4) * a(4,3) + dos2(3) * a(3,4) )
                    !
                    sdos(ie) = sdos(ie) &
                    &        + a(2,5) * a(1,5) * a(3,4) * a(4,2) / (e(5) - e0(ie)) &
                    &        * ( dos2(5) * a(5,2) + dos2(2) * a(2,5) &
                    &          + dos2(5) * a(5,1) + dos2(1) * a(1,5) &
                    &          + dos2(4) * a(4,3) + dos2(3) * a(3,4) &
                    &          + dos2(4) * a(4,2) + dos2(2) * a(2,4) )
                    !
                    sdos(ie) = sdos(ie) &
                    &        + a(1,5) * a(3,4) * a(2,4) * a(4,1) / (e(5) - e0(ie)) &
                    &        * ( dos2(5) * a(5,1) + dos2(1) * a(1,5) &
                    &          + dos2(4) * a(4,3) + dos2(3) * a(3,4) &
                    &          + dos2(4) * a(4,2) + dos2(2) * a(2,4) &
                    &          + dos2(4) * a(4,1) + dos2(1) * a(1,4) )
                    !
                 else if(e(4) < e0(ie) .and. e0(ie) <= e(5)) then
                    !
                    sdos(ie) = sdos(ie) &
                    &        + a(4,5) * a(3,5) * a(2,5) * a(1,5) / (e(5) - e0(ie)) &
                    &        * ( dos2(5) * a(5,4) + dos2(4) * a(4,5) &
                    &          + dos2(5) * a(5,3) + dos2(3) * a(3,5) &
                    &          + dos2(5) * a(5,2) + dos2(2) * a(2,5) &
                    &          + dos2(5) * a(5,1) + dos2(1) * a(1,5) )
                    !
                 end if
                 !
              end do ! ie
              !
           end do ! ip
           !
        end do ! ix
     end do ! ib
  end do ! ik
  !$OMP END DO
  !$OMP END PARALLEL
  !
  sdos(1:ne) = sdos(1:ne) / dble(24) * (e0(2) - e0(1))
  !
  call MPI_allREDUCE(MPI_IN_PLACE, sdos, ne, &
  &                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  !
end subroutine calc_sdos
!
!
!
subroutine penta_4d(iv4d)
  !
  !
  integer,intent(out) :: iv4d(4,5,24)
  !
  integer :: i1, i2, i3, i4, ip
  !
  ip = 0
  do i1 = 1, 4
     do i2 = 1, 4
        if(i2 == i1) cycle
        do i3 = 1, 4
           if(i3 == i1 .or. i3 == i2) cycle
           do i4 = 1, 4
              if(i4 == i1 .or. i4 == i2 .or. i4 == i3) cycle
              !
              ip = ip + 1
              !
              iv4d(1:4,1,ip) = 0
              !
              iv4d(1:4,2,ip) = iv4d(1:4,1,ip)
              iv4d( i1,2,ip) = iv4d( i1,2,ip) + 1
              !
              iv4d(1:4,3,ip) = iv4d(1:4,2,ip)
              iv4d( i2,3,ip) = iv4d( i2,3,ip) + 1
              !
              iv4d(1:4,4,ip) = iv4d(1:4,3,ip)
              iv4d( i3,4,ip) = iv4d( i3,4,ip) + 1
              !
              iv4d(1:4,5,ip) = iv4d(1:4,4,ip)
              iv4d( i4,5,ip) = iv4d( i4,5,ip) + 1
              !
           end do
        end do
     end do
  end do
  !
end subroutine penta_4d
!
! write QPDOS
!
subroutine write_dos()
  !
  use qpdos_vals, only : ne, e0, sdos
  !
  integer :: ie, fo = 20
  !
  open(fo, file = 'qpdos.dat')
  !
  do ie = 1, ne
     write(fo,*) e0(ie) * 13.60569228d3, sdos(ie)
  end do
  !
  close(fo)
  !
end subroutine write_dos
!
end module qpdos_routines
!
! Main routine
!
program qpdos
  !
  use mpi
  use omp_lib
  use qpdos_vals, only : petot, my_rank, Z, dltf
  use qpdos_routines,only : read_stdin, read_dat, read_elph, read_Coulomb, &
  &                         read_delta, make_Z, gapeq, average_matrix, egrid, &
  &                         tetra_type, calc_dosk, calc_sdos, write_dos
  !
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
  call egrid()
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
     write(*,*) "#####  Compute qpdos  #####"
     write(*,*)
  end if
  t1 = OMP_GET_WTIME()
  call tetra_type()
  call calc_dosk()
  call calc_sdos()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,*) "  Time : ", t2 - t1, " sec"
  !
  if(my_rank == 0) call write_dos()
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
end program qpdos


module scdft_vals
  !
  implicit none
  !
  integer,parameter :: &
  & itrmax = 100, & ! Max # of iteration
  & DP = 8
  !
  real(8),parameter :: &
  & er0 = 1d-15, &  ! Convergence threshold
  & pi = acos(-1d0) ! Pi = 3.14159265358979...
  !
  integer,save :: &
  & ncf,          & ! # of points for Chebyshev interpolation
  & fbee,         & ! First band for gap eq
  & lbee,         & ! Last band for gap eq
  & nbee,         & ! = lbee - fbee + 1
  & nbep,         & ! # of bands for elph
  & fbep,         & ! First band for elph
  & lbep,         & ! Last band for elph
  & nbf,          & ! # of bands contais FS
  & fstb,         & ! first band contais FS
  & lstb,         & ! last band contais FS
  & nth,          & ! # of threads
  & ltetra,       & ! Type of tetra
  & petot,        & ! # of PEs
  & my_rank,      & ! Index of PE
  & nmf,          & ! # of Matsubara frequencies
  & nt1,          & ! Total # of Delta, Xi, ... for grid 1 (with Gamma)
  & nt2,          & ! Total # of Delta, Xi, ... for grid 2 (w/o Gamma)
  & nt,           & ! max(nt1,nt2)
  & nb,           & ! # of bands
  & ng(3),        & ! k grid for matrix elements
  & nk,           & ! ng(1) * ng(2) * ng(3). # of total k points for matrix elements
  & fstk,         & ! First k
  & lstk,         & ! Last k
  & nk0,          & ! # of irreducible k points
  & ngd(3),       & ! k grid for DOS
  & nkd,          & ! ngd(1) * ngd(2) * ngd(3). # of total k points for DOS
  & nm,           & ! # of modes
  & nsym,         & ! # of symmetries
  & nx,           & ! # of energy scale
  & ivvec(3,20,6)   ! points for tetrahedron method
  !
  real(8),save :: &
  & xic,       & ! Cut off Xi
  & bvec(3,3), & ! Reciplocal lattice vector
  & beta,      & ! inversed temperature [Ry]
  & emin,      & ! Minimum energy scale [Ry]
  & wlsm(4,20)   ! Weights for tetrahedron method
  !
  integer,allocatable,save :: &
  & kindx(:,:),  & ! (nt,2) k point for gap equation for grid 1
  & bindx(:,:),  & ! (nt,2) band index for gap equation for grid 1
  & sym(:,:,:)     ! (3,3,nsym) Symmetrical operators
  !
  real(8),allocatable,save :: &
  & rkv(:,:),       & ! (3,nk0) irreducible k points
  & grid(:,:),      & ! (3,nk) k points
  & mf(:),          & ! (nmf) Matsubara frequency
  & wmf(:),         & ! (nmf) Weight for frequency integration
  & effint(:,:),    & ! (nt1,nt2*) Effective interaction    
  & Z(:,:),         & ! (nt,2) Renormalization fuctor for grid 1      
  & omg0(:,:),      & ! (nm,nk0) Phonon frequencies [Ry]
  & gg0(:,:,:,:,:), & ! (nm,nb,nb,nk,nk0) El-Ph matrix element [Ry]
  & Vc0(:,:,:,:,:), & ! (0:nmf,nb,nb,nk,nk0) Screened Coulomb matrix element [Ry]
  & omg(:,:,:),     & ! (nm,nk,nk*) Phonon frequencies [Ry]
  & gg(:,:,:,:,:),  & ! (nm,nb,nk,nb,nk*) El-Ph matrix element [Ry]
  & Vc(:,:,:,:,:),  & ! (0:nmf,nb,nk,nb,nk*) Screened Coulomb matrix element [Ry]
  & delta(:,:),     & ! (nt) Kohn-Sham gap functions [Ry]
  & xi(:,:),        & ! (nt) Kohn-Sham energy [Ry]
  & dk(:,:),        & ! (nt) Weight of k
  & eig(:,:,:,:),   & ! (nb,ngd(1),ngd(2),ngd(3)) Kohn-Sham energy [Ry]
  & xi0(:),         & ! (nx) energy scale [Ry]
  & dx0(:)            ! (nx) weight for energy
  !
end module scdft_vals
!
! Routines for scdft
!
module scdft_routines
  !
  implicit none
  !
  interface
     !
     ! BLAS
     !
     subroutine daxpy(n,da,dx,incx,dy,incy)
       double precision dx(*),dy(*),da
       integer          incx,incy,n
     END subroutine daxpy
     !
     double precision function ddot(n,dx,incx,dy,incy)
       double precision dx(*),dy(*)
       integer          incx,incy,n
     END function ddot
     !
     subroutine dgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )
       double precision   alpha, beta
       integer            incx, incy, lda, m, n
       character*1        trans
       double precision   a( lda, * ), x( * ), y( * )
     END subroutine dgemv
     !
  end interface
  !
contains
  !
 ! Standard input
!
subroutine read_stdin()
  !
  use mpi
  use scdft_vals, only : my_rank, beta, nx, emin, ltetra, xic, fbee, lbee, nbee, nmf
  !
  integer :: ierr
  namelist /input/ beta, nx, emin, ltetra, xic, fbee, lbee, nmf
  !
  if(my_rank == 0) then
     !
     xic = -1d0
     ltetra = 1
     !
     write(*,*) ''
     write(*,*) '###############  Standard Input  ###############'
     write(*,*) ''
     !
     read(*,input,err=100)
     write(*,*) '              Temparature[K] : ', beta
     beta = 157887d0 / beta
     write(*,*) '    Inverse temparature[/Ry] : ', beta
     write(*,*) '                     # of xi : ', nx
     write(*,*) '        Minimum energy scale : ', emin
     write(*,*) "                   Xi cutoff : ", xic
     write(*,*) "                  First band : ", fbee
     write(*,*) "                   Last band : ", lbee
     write(*,*) "  # of Matsubara frequencies : ", nmf     
     !
  end if
  !
  call MPI_BCAST(nx,     1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ltetra, 1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(fbee,   1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(lbee,   1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(beta,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(emin,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(xic,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nmf,    1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  nbee = lbee - fbee + 1
  !
  return
  !
100 write(*,*) "Stop in read_stdin. reading namelist file"
  write(*,*) '           Temparature[K] : ', beta
  beta = 157887d0 / beta
  write(*,*) 'Inverse temparature[/Ry] : ', beta
  write(*,*) '                  # of xi : ', nx
  write(*,*) '     Minimum energy scale : ', emin
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
  use scdft_vals, only : my_rank, nsym, sym, nb, ngd, nkd, eig, bvec, nbf, fstb, lstb
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
     allocate(eig(nb, ngd(1), ngd(2), ngd(3)))
     open(fi, file = 'eigval.dat')
     write(*,*) "  read from eigval.dat"
     read(fi,*) eig(1:nb,1:ngd(1),1:ngd(2),1:ngd(3))
     close(fi)
     !
  end if
  !
  call MPI_BCAST(nsym, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, isym)
  if(my_rank /= 0) allocate(sym(3,3,nsym))
  call MPI_BCAST(sym,  9 * nsym, MPI_INTEGER, 0, MPI_COMM_WORLD, isym)
  call MPI_BCAST(nb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, isym)
  call MPI_BCAST(ngd, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, isym)
  call MPI_BCAST(nkd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, isym)
  if(my_rank /= 0) allocate(eig(nb, ngd(1), ngd(2), ngd(3)))
  call MPI_BCAST(eig, nb * nkd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, isym)
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
  use scdft_vals, only : my_rank, ng, nk, nk0, nm, gg0, omg0, rkv, &
  &                  nbep, fbep, lbep
  !use scdft_routines : cnt_and_dsp 
  use omp_lib
  !
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
  use scdft_vals, only : petot, my_rank
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
! Compute cnt and dsp
!
subroutine cnt_and_dsp_omp(n,cnt1,dsp1)
  !
  use scdft_vals, only : nth
  use omp_lib
  implicit none
  !
  integer,intent(in) :: n
  integer,intent(out) :: cnt1, dsp1
  !
  integer :: ii, ith
  integer :: cnt(0:nth-1), dsp(0:nth-1)
  !
  ith = OMP_GET_THREAD_NUM()
  !
  cnt(0:nth-1)        = n / nth
  cnt(0:mod(n,nth)-1) = n / nth + 1
  dsp(0) = 0
  do ii = 1, nth - 1
     dsp(ii) = dsp(ii - 1) + cnt(ii - 1)
  end do
  !
  cnt1 = cnt(ith)
  dsp1 = dsp(ith)
  !
end subroutine cnt_and_dsp_omp
!
! Irreducible Brillouin zone
!
subroutine irr_bz()
  !
  use scdft_vals, only : my_rank, nk, nk0, ng, rkv, grid, nsym, sym
  !
  integer :: i1, i2, i3, ik, isym, nkk
  real(8) :: kv0(3,nk), kv1(3), kv2(3)
  !
  allocate(grid(3,nk))
  !
  nk0 = 0
  nkk = 0
  !
  do i1 = 1, ng(1)
     do i2 = 1, ng(2)
        do i3 = 1, ng(3)
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
  use scdft_vals, only : my_rank, nk0, nk, nb, Vc0, rkv, ng, ncf, &
  &                  nmf
  !scdft_routines, only : cnt_and_dsp
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
! Define initial delta
!
subroutine ini_delta()
  !
  use mpi
  use scdft_vals, only : my_rank, nk0, nm, nk0, omg0, nt, nt1, nt2, fstk, lstk,&
  &                  xi, delta, dk, kindx, bindx, emin, nx, xi0, dx0
  !use scdft_routines, only : cnt_and_dsp, energy_grid
  !
  integer :: it, fi = 10, is, cnt, dsp, ierr
  real(8) :: thr, Z0, dosf
  character(1) :: tmp
  !
  if(my_rank == 0) open(fi, file = 'delta.dat',status="old", action = 'read',iostat = is)
  !
  call MPI_BCAST(is, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  !
  if(is == 0) then
     !
     if(my_rank == 0) then
        !
        read(fi,*) tmp, nt2, nt1
        !
        nt = max(nt1, nt2)
        write(*,*) "  # of total points for gap equation : ", nt1, nt2
        allocate(xi(nt,2), delta(nt,2), dk(nt,2), kindx(nt,2), bindx(nt,2))
        xi(   1:nt,1:2) = 0d0
        delta(1:nt,1:2) = 0d0
        dk(   1:nt,1:2) = 0d0
        kindx(1:nt,1:2) = 0
        bindx(1:nt,1:2) = 0
        !
        do it = 1, nt2
           read(fi,*) xi(it,2), delta(it,2), Z0, dk(it,2), kindx(it,2), bindx(it,2)
        enddo
        !
        do it = 1, nt1
           read(fi,*) xi(it,1), delta(it,1), Z0, dk(it,1), kindx(it,1), bindx(it,1)
        enddo
        close(fi)
        !
     end if ! (my_rank == 0)
     !
     call MPI_BCAST(nt,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(nt1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(nt2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     if(my_rank /= 0) allocate(xi(nt,2), delta(nt,2),    dk(nt,2), &
     &                         kindx(nt,2), bindx(nt,2))
     call MPI_BCAST(xi,    nt * 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(delta, nt * 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(dk,    nt * 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(kindx, nt * 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(bindx, nt * 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     !
     if(maxval(abs(delta(1:nt1,1))) > 1d-8) then
        !
        if(my_rank == 0) write(*,*) '  Initial delta is read from file (delta.dat) '
        !
        call energy_grid()
        !
        dosf = sum(dk(1:nt1,1), abs(xi(1:nt1,1)) < 1d-12)
        dosf = dosf / dx0(minloc(abs(xi0(1:nx)), 1))
        if(my_rank == 0) write(*,*) "  DOS(E_F)[Ry^-1/cell/spin] : ", dosf
        !
        dosf = sum(dk(1:nt2,2), abs(xi(1:nt2,2)) < 1d-12)
        dosf = dosf / dx0(minloc(abs(xi0(1:nx)), 1))
        if(my_rank == 0) write(*,*) "  DOS(E_F)[Ry^-1/cell/spin] : ", dosf
        !
        goto 5
        !
     end if ! (maxval(abs(delta(1:nt1,1))) > 1d-8)
     !
     ! If delta read from file is too small
     !
     if(my_rank == 0) write(*,*) 'Delta from file is too small !'
     !
     deallocate(xi, delta, dk, kindx, bindx)
     !
  end if ! (is == 0)
  !
  if(my_rank == 0) write(*,*) '  Initial delta is theta function '
  !
  call energy_grid()
  !
  call tetra_type()
  !
  call compute_d3k()
  !
  allocate(delta(nt,2))
  !
  call cnt_and_dsp(nk0,cnt,dsp)
  thr = maxval(omg0(1:nm,1:cnt))
  call MPI_allREDUCE(MPI_IN_PLACE, thr, 1, MPI_DOUBLE_PRECISION, &
  &                  MPI_MAX, MPI_COMM_WORLD,ierr)
  !
  if(my_rank == 0) then
     call random_seed()
     call random_number(delta(1:nt,1:2))
     delta(1:nt,1:2) = delta(1:nt,1:2) * thr
  end if
  call MPI_BCAST(delta, nt * 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  !
  !delta(1:nt,1:2) = thr**2 / (thr + xi(1:nt,1:2))
  !
5 continue
  !
  call cnt_and_dsp(nt1,cnt,dsp)
  fstk = minval(kindx(dsp + 1:dsp + cnt,1))
  lstk = maxval(kindx(dsp + 1:dsp + cnt,1))
  !
end subroutine ini_delta
!
! Compute energy groid
!
subroutine energy_grid()
  !
  use scdft_vals, only : eig, ngd, xi0, dx0, nx, emin, fbep, lbep
  !use scdft_routines, only : weightspoints_gl
  !
  !
  !
  integer :: ix
  real(8) :: xmax, xmin, xx, ww, rhs(nx)
  !
  ! Make logscale
  !
  xmax = maxval(eig(fbep:lbep,1:ngd(1),1:ngd(2),1:ngd(3)))
  xmin = minval(eig(fbep:lbep,1:ngd(1),1:ngd(2),1:ngd(3)))
  !
  allocate(xi0(nx), dx0(nx))
  call weightspoints_gl(nx,xi0,dx0)
  !
  do ix = 1, nx
     rhs(ix) = (1d0 + xi0(ix) - 2d0 / dble(nx)) &
     &         * log(  xmax / (0.5d0 * dble(nx) * emin * (1d0 - xi0(ix)))) &
     &       - (1d0 - xi0(ix) - 2d0 / dble(nx)) &
     &         * log(- xmin / (0.5d0 * dble(nx) * emin * (1d0 + xi0(ix))))
     rhs(ix) = abs(rhs(ix))
  end do
  !
  xx = xi0(minloc(rhs(1:nx),1))
  ww = max(log(- xmin / (0.5d0 * dble(nx) * emin * (1d0 + xx))) / (1d0 + xx - 2d0 / dble(nx)), &
  &        log(  xmax / (0.5d0 * dble(nx) * emin * (1d0 - xx))) / (1d0 - xx - 2d0 / dble(nx))  )
  !
  do ix = 1, nx
     dx0(ix) = dx0(ix) * (1d0 + ww * abs(xi0(ix) - xx)) &
     &                        * emin * 0.5d0 * dble(nx) * exp(ww * (abs(xi0(ix) - xx) - 2d0 / dble(nx))) 
     xi0(ix) = (xi0(ix) - xx) * emin * 0.5d0 * dble(nx) * exp(ww * (abs(xi0(ix) - xx) - 2d0 / dble(nx)))
  end do
  !
end subroutine energy_grid
!
! Weights & Points for Gauss-Legendre method
!
subroutine weightspoints_gl(n,x,w)
  !
  use mpi
  use scdft_vals, only: pi
  !use scdft_routines, only : legendre
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
! define shortest diagonal line & define type of tetragonal
!
subroutine tetra_type()
  !
  use mpi
  use scdft_vals, only : DP, ng, bvec, ivvec, wlsm, ltetra, my_rank
  !
  integer :: itype, i1, i2, i3, itet, divvec(4,4), ivvec0(4), ierr
  real(DP) :: l(4), bvec2(3,3), bvec3(3,4)
  !
  do i1 = 1, 3
     bvec2(1:3,i1) = bvec(1:3,i1) / real(ng(i1), DP)
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
  if(ltetra == 1) then
     !
     if(my_rank == 0) write(*,*) "Linear tetrahedron method is used"
     !
     wlsm(1:4,1:20) = 0.0_dp
     wlsm(1,1:4) = real((/1, 0, 0, 0/), DP)
     wlsm(2,1:4) = real((/0, 1, 0, 0/), DP)
     wlsm(3,1:4) = real((/0, 0, 1, 0/), DP)
     wlsm(4,1:4) = real((/0, 0, 0, 1/), DP)
     !
  else if(ltetra == 2) then
     !
     if(my_rank == 0) write(*,*) "Improved tetrahedron method is used"
     !
     wlsm(1, 1:10) = real((/1440,    0,   30,    0, -38, -56, -38, -28,  7,     9/), DP)
     wlsm(2, 1:10) = real((/   0, 1440,    0,   30, -28,   9,   7, -38, -38,  -56/), DP)
     wlsm(3, 1:10) = real((/  30,    0, 1440,    0,  17, -46,  17,   7, -28,    9/), DP)
     wlsm(4, 1:10) = real((/   0,   30,    0, 1440,   7,   9, -28,  17,  17,  -46/), DP)
     !
     wlsm(1,11:20) = real((/ -46,   17,   17,  -28,   9,   7, -18, -18,  12,  -18/), DP)
     wlsm(2,11:20) = real((/   9,    7,  -28,   17, -46,  17, -18, -18, -18,   12/), DP)
     wlsm(3,11:20) = real((/ -56,  -38,  -38,    7,   9, -28,  12, -18, -18,  -18/), DP)
     wlsm(4,11:20) = real((/   9,  -28,    7,  -38, -56, -38, -18,  12, -18,  -18/), DP)
     !
     wlsm(1:4,1:20) = wlsm(1:4,1:20) / 1260
     !
  else
     !
     write(*,*) "ltetra is wrong"
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     call MPI_FINALIZE(ierr)
     stop
     !
  end if
  !
end subroutine tetra_type
!
! Compute integration weights
!
subroutine compute_d3k()
  !
  use mpi
  use scdft_vals, only : nt, nt1, nt2, nx, nk, ng, ngd, &
  &                  kindx, bindx, grid, my_rank, fbep, lbep, &
  &                  xi, dk, xi0, eig, fbee, lbee
  !use scdft_routines, only : calc_dosk, symm_dosk, store_d3k
  !
  integer :: ik, ib, ix, shift(3)
  real(8) :: thr = 1d-12, dos(nx,fbep:lbep,nk,2)
  !
  call calc_dosk(dos)
  !
  call symm_dosk(dos)
  !
  !############  For debug  ###################
  !
  if(my_rank == 0) then
     do ix = 2, nx - 1
        write(90,*) xi0(ix), sum(dos(ix,fbep:lbep,1:nk,1)), sum(dos(ix,fbep:lbep,1:nk,2))
     end do ! ix
  end if
  !
  !############  End for debug  ###############
  !
  ix = minloc(abs(xi0(1:nx)),1)
  if(my_rank == 0) then
     write(*,*) "  DOS(E_F)[Ry^-1/cell/spin] : ", sum(dos(ix,fbep:lbep,1:nk,1)) 
     write(*,*) "  DOS(E_F)[Ry^-1/cell/spin] : ", sum(dos(ix,fbep:lbep,1:nk,2)) 
  end if
  !
  ! Query of nt1
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! modified by Kentaro Tsutsumi !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(223, file = "nt1.dat")
  nt1 = 0
  do ik = 1, nk
     do ib = fbee, lbee
        !
        if(fbep <= ib .and. ib <= lbep) then
           nt1 = nt1 + max(1, count(dos(1:nx,ib,ik,1) > thr))
           write(223, *) nt1, ib, ik
        else
           nt1 = nt1 + 1
           write(223, *) nt1, ib, ik
        end if
        !
     end do
  end do
  write(223, *) "nt1 =", nt1
  close(223)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! modified by Kentaro Tsutsumi !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Query of nt2
  !
  nt2 = 0
  do ik = 1, nk
     do ib = fbee, lbee
        !
        if(fbep <= ib .and. ib <= lbep) then
           nt2 = nt2 + max(1, count(dos(1:nx,ib,ik,2) > thr))
        else  
           nt2 = nt2 + 1
        end if
        !
     end do
  end do
  !
  if(my_rank == 0) write(*,*) "  # of total points for gap equation : ", nt1, nt2
  nt = max(nt1, nt2)
  allocate(xi(nt,2), dk(nt,2), kindx(nt,2), bindx(nt,2))
  !
  xi(   1:nt,1:2) = 0d0
  dk(   1:nt,1:2) = 0d0
  kindx(1:nt,1:2) = 0
  bindx(1:nt,1:2) = 0
  !
  ! Map xi, d3k, dosk, indx, bindx
  !
  shift(1:3) = 0
  call store_d3k(dos(1:nx,fbep:lbep,1:nk,1),shift,nt1,xi(1:nt1,1), &
  &              dk(1:nt1,1),kindx(1:nt1,1),bindx(1:nt1,1))
  shift(1:3) = 1
  call store_d3k(dos(1:nx,fbep:lbep,1:nk,2),shift,nt2,xi(1:nt2,2), &
  &              dk(1:nt2,2),kindx(1:nt2,2),bindx(1:nt2,2))
  !
end subroutine compute_d3k
!
! Compute Dos_k, xmin, xmax
! 
subroutine calc_dosk(dos)
  !
  use mpi
  use scdft_vals, only : DP, ngd, ng, nkd, nk, ivvec, nx, xi0, eig, wlsm, fbep, lbep, nbep
  !use scdft_routines, only : interpol_weight, cnt_and_dsp, cnt_and_dsp_omp, sort
  !
  real(8),intent(out) :: dos(nx,fbep:lbep,nk,2)
  !
  integer :: ik, ib, i1, i2, i3, it, ii, ix, ikv1(3), dgrid(3,nkd), &
  &          ik2, dsp, cnt, ierr
  real(8) :: ei(4,fbep:lbep), e(4), dosd(nx,fbep:lbep,nkd), a(4,4), w0(4,4), w1(4,nx), tmp(5,4), V, kv(3)
  !
  ! Dense grid
  !
  ik = 0
  do i3 = 1, ngd(3)
     do i2 = 1, ngd(2)
        do i1 = 1, ngd(1)
           ik = ik + 1
           dgrid(1:3,ik) = (/i1, i2, i3/) - 1
        end do ! i1
     end do ! i2
  end do ! i3
  !
  ! Calc. Dosd_k
  !
  dosd(1:nx, fbep:lbep, 1:nkd) = 0d0
  w0(1:4,1:4) = 0d0
  do ii = 1, 4
     w0(ii,ii) = 1d0
  end do
  !
  call cnt_and_dsp(nkd,cnt,dsp)
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(dsp,cnt,nx,ivvec,ngd,dgrid,wlsm,eig,xi0,dosd,w0,fbep,lbep) &
  !$OMP & PRIVATE(ik,it,ii,ib,ix,ikv1,ei,w1,a,e,tmp,V,ik2)
  !
  do ik = dsp + 1, dsp + cnt
     !
     do it = 1, 6
        !
        ei(1:4,fbep:lbep) = 0.0_dp
        !
        do ii = 1, 20
           !
           ikv1(1:3) = dgrid(1:3,ik) + ivvec(1:3, ii, it)
           ikv1(1:3) = modulo(ikv1(1:3), ngd(1:3)) + 1
           !
           do ib = fbep, lbep
              ei(1:4,ib) = ei(1:4,ib) &
              &          + wlsm(1:4,ii) * eig(ib,ikv1(1),ikv1(2),ikv1(3))
           end do
           !
        end do
        !
        do ib = fbep, lbep
           !
           w1(1:4,1:nx) = 0d0
           !
           tmp(  1,1:4) = ei(1:4,ib)
           tmp(2:5,1:4) = w0(1:4,1:4)
           call sort(5, 4, tmp)
           e(1:4) = tmp(1,1:4)
           !
           !$OMP DO
           do ix = 1, nx
              !
              do ii = 1, 4
                 a(ii,1:4) = (xi0(ix) - e(1:4)) / (e(ii) - e(1:4))
              end do
              !
              if(e(1) < xi0(ix) .and. xi0(ix) <= e(2)) then
                 !
                 ! A
                 !
                 V = a(2,1) * a(3,1) * a(4,1) / (xi0(ix) - e(1))
                 !
                 w1(1:4,ix) = V * ( tmp(2:5,1) * a(1,2) + tmp(2:5,2) * a(2,1) &
                 &                + tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1) &
                 &                + tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1) )
                 !
              else if( e(2) < xi0(ix) .and. xi0(ix) <= e(3)) then
                 !
                 ! B - 1
                 !
                 V = a(3,1) * a(4,1) * a(2,4) / (xi0(ix) - e(1))
                 !
                 w1(1:4,ix) = V * ( tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1) &
                 &                + tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1) &
                 &                + tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) )
                 !
                 ! B - 2
                 !
                 V = a(2,3) * a(3,1) * a(4,2) / (xi0(ix) - e(1))
                 !
                 w1(1:4,ix) = w1(1:4,ix)                                      &
                 &          + V * ( tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1) &
                 &                + tmp(2:5,2) * a(2,3) + tmp(2:5,3) * a(3,2) &
                 &                + tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) )
                 !
              else if(e(3) < xi0(ix) .and. xi0(ix) < e(4)) then
                 !
                 ! C
                 !
                 V = a(1,4) * a(2,4) * a(3,4) / (e(4) - xi0(ix))
                 !
                 w1(1:4,ix) = V * ( tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1) &
                 &                + tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) &
                 &                + tmp(2:5,3) * a(3,4) + tmp(2:5,4) * a(4,3) )
                 !
              end if
              !
           end do ! ix
           !$OMP END DO NOWAIT
           !
           do ii = 1, 20
              !
              ikv1(1:3) = dgrid(1:3,ik) + ivvec(1:3, ii, it)
              ikv1(1:3) = modulo(ikv1(1:3), ngd(1:3))
              ik2 = 1 + ikv1(1) + ngd(1) * ikv1(2) + ngd(1) * ngd(2) * ikv1(3)
              !
              !$OMP DO
              do ix = 1, nx
                 dosd(ix,ib,ik2) = dosd(ix,ib,ik2) + sum(wlsm(1:4,ii) * w1(1:4,ix))
              end do ! ix
              !$OMP END DO NOWAIT
              !
           end do ! ii
           !
        end do ! ib
        !
      end do ! it
     !
  end do ! ik
  !$OMP END PARALLEL
  !
  dosd(1:nx, fbep:lbep, 1:nkd) = dosd(1:nx, fbep:lbep, 1:nkd) / dble(6 * nkd)
  !
  ! Interpolation of weight
  !
  dos(1:nx,fbep:lbep,1:nk,1:2) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(dsp,cnt,nkd,dgrid,ngd,nk,nx,ng,dosd,dos,fbep,lbep,nbep) &
  !$OMP & PRIVATE(ik,kv)
  !
  !$OMP DO REDUCTION(+: dos)
  do ik = 1, nkd
     !
     kv(1:3) = dble(dgrid(1:3,ik)) / dble(ngd(1:3))
     call interpol_weight(nk,nx * nbep, ng, kv, dosd(1:nx,fbep:lbep,ik), dos(1:nx,fbep:lbep,1:nk,1))
     !
     kv(1:3) = dble(dgrid(1:3,ik)) / dble(ngd(1:3)) - 0.5d0 / dble(ng(1:3))
     call interpol_weight(nk,nx * nbep, ng, kv, dosd(1:nx,fbep:lbep,ik), dos(1:nx,fbep:lbep,1:nk,2))
     !
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  !
  call MPI_allREDUCE(MPI_IN_PLACE, dos, nx * nbep * nk * 2, MPI_DOUBLE_PRECISION, &
  &                  MPI_SUM,MPI_COMM_WORLD,ierr)
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
        atmp(1:n1) = a(1:n1, m)
        a(1:n1, m) = a(1:n1,i)
        a(1:n1, i) = atmp(1:n1)
     end if
  end do
  !
end subroutine sort
!
! first or third order interpolation of weights
!
subroutine interpol_weight(nk,nb,ng,ko,wi,wo)
  !
  use mpi
  use scdft_vals, only : ivvec, ltetra
  !
  integer,intent(in)  :: nk, nb, ng(3)
  real(8),intent(in)  :: ko(3)
  real(8),intent(in) ::wi(nb)
  real(8),intent(inout) :: wo(nb,nk)
  !
  integer :: ikv(3), ikv1(3), ik(20), ii, it, it0, ierr
  real(8) :: rot(3,3), res(3), prod(3), u, x, y, z, thr = 1d-10
  !
  rot(1:3,1) = (/  2d0, - 1d0,   0d0/)
  rot(1:3,2) = (/- 1d0,   2d0, - 1d0/)
  rot(1:3,3) = (/  0d0, - 1d0,   1d0/)
  !
  ! Search nearest neighbor grid points.
  !
  res(1:3) = ko(1:3) * dble(ng(1:3))
  ikv(1:3) = floor(res(1:3))
  res(1:3) = res(1:3) - dble(ikv(1:3))
  !
  do it = 1, 6
     !
     do ii = 1, 3
        prod(ii) = dot_product(dble(ivvec(1:3,1 + ii,it) - ivvec(1:3,1,it)), &
        &                                  res(1:3) - dble(ivvec(1:3,1,it))  )
     end do
     !
     prod(1:3) = matmul(rot(1:3,1:3), prod(1:3))
     !
     if(minval(prod(1:3)) > - thr .and. sum(prod(1:3)) < 1d0 + thr) then
        it0 = it
        goto 10
     end if
     !
  end do
  !
  write(*,*) "Stop in interpol_weight."
  call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  call MPI_FINALIZE(ierr)
  stop
  !
10 continue
  !
  x = prod(1)
  y = prod(2)
  z = prod(3)
  u = 1d0 - x - y - z
  !
  do ii = 1, 20
     !
     ikv1(1:3) = ikv(1:3) + ivvec(1:3,ii,it0)
     ikv1(1:3) = modulo(ikv1(1:3), ng(1:3))
     ik(ii) = 1 + ikv1(1) + ikv1(2) * ng(1) + ikv1(3) * ng(1) * ng(2)
     !
  end do
  !
  if(ltetra == 0 .or. ltetra == 1) then
     !
     wo(1:nb,ik(1)) = wo(1:nb,ik(1)) + wi(1:nb) * u
     wo(1:nb,ik(2)) = wo(1:nb,ik(2)) + wi(1:nb) * x
     wo(1:nb,ik(3)) = wo(1:nb,ik(3)) + wi(1:nb) * y
     wo(1:nb,ik(4)) = wo(1:nb,ik(4)) + wi(1:nb) * z
     !
  else if(ltetra == 2) then
     !
     wo(1:nb,ik( 1)) = wo(1:nb,ik( 1)) + wi(1:nb) * 0.5d0 * u * (2d0 + u * (1d0 - u) + 2d0 * y * (x + z))
     wo(1:nb,ik( 2)) = wo(1:nb,ik( 2)) + wi(1:nb) * 0.5d0 * x * (2d0 + x * (1d0 - x) + 2d0 * z * (u + y))
     wo(1:nb,ik( 3)) = wo(1:nb,ik( 3)) + wi(1:nb) * 0.5d0 * y * (2d0 + y * (1d0 - y) + 2d0 * u * (x + z))
     wo(1:nb,ik( 4)) = wo(1:nb,ik( 4)) + wi(1:nb) * 0.5d0 * z * (2d0 + z * (1d0 - z) + 2d0 * x * (u + y))
     wo(1:nb,ik( 5)) = wo(1:nb,ik( 5)) + wi(1:nb) * x * u * (2d0 * y - u - 1d0) / 6d0
     wo(1:nb,ik( 6)) = wo(1:nb,ik( 6)) + wi(1:nb) * x * y * (2d0 * z - x - 1d0) / 6d0
     wo(1:nb,ik( 7)) = wo(1:nb,ik( 7)) + wi(1:nb) * y * z * (2d0 * u - y - 1d0) / 6d0
     wo(1:nb,ik( 8)) = wo(1:nb,ik( 8)) + wi(1:nb) * z * u * (2d0 * x - z - 1d0) / 6d0
     wo(1:nb,ik( 9)) = wo(1:nb,ik( 9)) + wi(1:nb) * y * u * (2d0 * y + u - 3d0) / 6d0
     wo(1:nb,ik(10)) = wo(1:nb,ik(10)) + wi(1:nb) * x * z * (2d0 * z + x - 3d0) / 6d0
     wo(1:nb,ik(11)) = wo(1:nb,ik(11)) + wi(1:nb) * y * u * (2d0 * u + y - 3d0) / 6d0
     wo(1:nb,ik(12)) = wo(1:nb,ik(12)) + wi(1:nb) * x * z * (2d0 * x + z - 3d0) / 6d0
     wo(1:nb,ik(13)) = wo(1:nb,ik(13)) + wi(1:nb) * z * u * (2d0 * y - u - 1d0) / 6d0
     wo(1:nb,ik(14)) = wo(1:nb,ik(14)) + wi(1:nb) * x * u * (2d0 * z - x - 1d0) / 6d0
     wo(1:nb,ik(15)) = wo(1:nb,ik(15)) + wi(1:nb) * x * y * (2d0 * u - y - 1d0) / 6d0
     wo(1:nb,ik(16)) = wo(1:nb,ik(16)) + wi(1:nb) * y * z * (2d0 * x - z - 1d0) / 6d0
     wo(1:nb,ik(17)) = wo(1:nb,ik(17)) + wi(1:nb) * (- x * z * u)
     wo(1:nb,ik(18)) = wo(1:nb,ik(18)) + wi(1:nb) * (- x * y * u)
     wo(1:nb,ik(19)) = wo(1:nb,ik(19)) + wi(1:nb) * (- x * y * z)
     wo(1:nb,ik(20)) = wo(1:nb,ik(20)) + wi(1:nb) * (- y * z * u)
     !
  else
     !
     write(*,*) "Stop in interpol_weight. ltetra is wrong."
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     call MPI_FINALIZE(ierr)
     stop
     ! 
  end if
  !
end subroutine interpol_weight
!
! Map d3k, dosk, xi, kindx, bindx
!
subroutine store_d3k(dos,shift,nt,xi,dk,kindx,bindx)
  !
  use mpi
  use scdft_vals, only : nx, nb, nk, nkd, ng, ngd, eig, xi0, dx0, grid, &
  &                  fbee, lbee, fbep, lbep
  !use scdft_routines, only : interpol
  !
  integer,intent(in) :: shift(3), nt
  real(8),intent(in) :: dos(nx,fbep:lbep,nk)
  real(8),intent(out) :: xi(nt), dk(nt)
  integer,intent(out) :: kindx(nt), bindx(nt)
  !
  integer :: nt0, ik, ib, ix, ierr
  real(8) :: thr = 1d-12, eig1(ngd(1),ngd(2),ngd(3)), kv(3)
  !
  nt0 = 0
  do ik = 1, nk
     !
     do ib = fbee, lbee
        !
        eig1(1:ngd(1),1:ngd(2),1:ngd(3)) = eig(ib,1:ngd(1),1:ngd(2),1:ngd(3))
        !
        if(ib < fbep .or. lbep < ib) then
           !
           nt0 = nt0 + 1
           !
           kv(1:3) = grid(1:3,ik) + dble(shift(1:3)) / dble(2 * ng(1:3))
           call interpol(nkd,1,ngd,kv,eig1,xi(nt0))
           !
           dk(nt0) = 1d0 / dble(nk)
           kindx(nt0) = ik
           bindx(nt0) = ib
           !
        else if(count(dos(1:nx,ib,ik) > thr) < 2) then
           !
           nt0 = nt0 + 1
           !
           kv(1:3) = grid(1:3,ik) + dble(shift(1:3)) / dble(2 * ng(1:3))
           call interpol(nkd,1,ngd,kv,eig1,xi(nt0))
           !
           dk(nt0) = 1d0 / dble(nk)
           kindx(nt0) = ik
           bindx(nt0) = ib
           !
        else
           !
           do ix = 1, nx                        
              !
              if(dos(ix,ib,ik) > thr) then
                 !
                 nt0 = nt0 + 1
                 !
                 xi(nt0) = xi0(ix)
                 dk(nt0) = dos(ix,ib,ik) * dx0(ix)
                 kindx(nt0) = ik
                 bindx(nt0) = ib
                 !
              end if
              !
           end do ! ix
           !
        end if ! count(dos(1:nx,ib,jk) > thr) == 0
        !
     end do ! ib
     !
  end do ! ik
  !
  if(nt0 /= nt) then
     write(*,*) "Stop in store_d3k. nt0 /= nt"
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     call MPI_FINALIZE(ierr)
     stop
  end if
  !
end subroutine store_d3k
!
! Linear interpolation
!
subroutine interpol(nk,nb,ng,ko,ei,eo)
  !
  use mpi
  use scdft_vals, only : ivvec, ltetra
  !
  integer,intent(in)  :: nk, nb, ng(3)
  real(8),intent(in)  :: ko(3), ei(nb,nk)
  real(8),intent(out) :: eo(nb)
  !
  integer :: ikv(3), ikv1(3), ik(20), ii, it, it0, ierr
  real(8) :: rot(3,3), res(3), prod(3), u, x, y, z, thr = 1d-10
  !
  rot(1:3,1) = (/  2d0, - 1d0,   0d0/)
  rot(1:3,2) = (/- 1d0,   2d0, - 1d0/)
  rot(1:3,3) = (/  0d0, - 1d0,   1d0/)
  !
  ! Search nearest neighbor grid points.
  !
  res(1:3) = ko(1:3) * dble(ng(1:3))
  ikv(1:3) = floor(res(1:3))
  res(1:3) = res(1:3) - dble(ikv(1:3))
  !
  do it = 1, 6
     !
     do ii = 1, 3
        prod(ii) = dot_product(dble(ivvec(1:3,1 + ii,it) - ivvec(1:3,1,it)), &
        &                                  res(1:3) - dble(ivvec(1:3,1,it))  )
     end do
     !
     prod(1:3) = matmul(rot(1:3,1:3), prod(1:3))
     !
     if(minval(prod(1:3)) > - thr .and. sum(prod(1:3)) < 1d0 + thr) then
        it0 = it
        goto 10
     end if
     !
  end do
  !
  write(*,*) "Stop in interpol"
  call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  call MPI_FINALIZE(ierr)
  stop
  !
10 continue
  !
  x = prod(1)
  y = prod(2)
  z = prod(3)
  u = 1d0 - x - y - z
  !
  do ii = 1, 20
     !
     ikv1(1:3) = ikv(1:3) + ivvec(1:3,ii,it0)
     ikv1(1:3) = modulo(ikv1(1:3), ng(1:3))
     ik(ii) = 1 + ikv1(1) + ikv1(2) * ng(1) + ikv1(3) * ng(1) * ng(2)
     !
  end do
  !
  if(ltetra == 1) then
     !
     eo(1:nb) = ei(1:nb,ik(1)) * u &
     &        + ei(1:nb,ik(2)) * x &
     &        + ei(1:nb,ik(3)) * y &
     &        + ei(1:nb,ik(4)) * z
     !
  else if(ltetra == 2) then
     !
     eo(1:nb) = ei(1:nb,ik( 1)) * 0.5d0 * u * (2d0 + u * (1d0 - u) + 2d0 * y * (x + z)) &
     &        + ei(1:nb,ik( 2)) * 0.5d0 * x * (2d0 + x * (1d0 - x) + 2d0 * z * (u + y)) &
     &        + ei(1:nb,ik( 3)) * 0.5d0 * y * (2d0 + y * (1d0 - y) + 2d0 * u * (x + z)) &
     &        + ei(1:nb,ik( 4)) * 0.5d0 * z * (2d0 + z * (1d0 - z) + 2d0 * x * (u + y)) &
     &        + ei(1:nb,ik( 5)) * x * u * (2d0 * y - u - 1d0) / 6d0 &
     &        + ei(1:nb,ik( 6)) * x * y * (2d0 * z - x - 1d0) / 6d0 &
     &        + ei(1:nb,ik( 7)) * y * z * (2d0 * u - y - 1d0) / 6d0 &
     &        + ei(1:nb,ik( 8)) * z * u * (2d0 * x - z - 1d0) / 6d0 &
     &        + ei(1:nb,ik( 9)) * y * u * (2d0 * y + u - 3d0) / 6d0 &
     &        + ei(1:nb,ik(10)) * x * z * (2d0 * z + x - 3d0) / 6d0 &
     &        + ei(1:nb,ik(11)) * y * u * (2d0 * u + y - 3d0) / 6d0 &
     &        + ei(1:nb,ik(12)) * x * z * (2d0 * x + z - 3d0) / 6d0 &
     &        + ei(1:nb,ik(13)) * z * u * (2d0 * y - u - 1d0) / 6d0 &
     &        + ei(1:nb,ik(14)) * x * u * (2d0 * z - x - 1d0) / 6d0 &
     &        + ei(1:nb,ik(15)) * x * y * (2d0 * u - y - 1d0) / 6d0 &
     &        + ei(1:nb,ik(16)) * y * z * (2d0 * x - z - 1d0) / 6d0 &
     &        + ei(1:nb,ik(17)) * (- x * z * u) &
     &        + ei(1:nb,ik(18)) * (- x * y * u) &
     &        + ei(1:nb,ik(19)) * (- x * y * z) &
     &        + ei(1:nb,ik(20)) * (- y * z * u)
     !
  else
     ! 
     write(*,*) "Stop in interpol. ltetra is wrong."
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     call MPI_FINALIZE(ierr)
     stop
     ! 
  end if
  !
end subroutine interpol
!
! Symmetrize Dos
!
subroutine symm_dosk(dos)
  !
  use scdft_vals, only : nk, fbep, lbep, nx, nsym, sym, ng, grid
  !
  real(8),intent(inout) :: dos(nx,fbep:lbep,nk,2)
  !
  integer :: ik, ik2, isym, ikv2(3), nks(nk)
  real(8) :: dostmp(nx,fbep:lbep,nk), kv1(3), kv2(3)
  !
  ! 1
  !
  nks(1:nk) = 0
  dostmp(1:nx,fbep:lbep,1:nk) = 0d0
  !
  do ik = 1, nk
     !
     kv1(1:3) = grid(1:3,ik)
     !
     do isym = 1, nsym
        !
        kv2(1:3) = matmul(dble(sym(1:3,1:3,isym)),kv1(1:3)) * dble(ng(1:3))
        ikv2(1:3) = nint(kv2(1:3))
        !
        if(any(abs(kv2(1:3) - dble(ikv2(1:3))) > 1d-8)) cycle
        !
        ikv2(1:3) = modulo(ikv2(1:3), ng(1:3))
        ik2 = 1 + ikv2(1) + ng(1) * ikv2(2) + ng(1) * ng(2) * ikv2(3)
        !
        nks(ik2) = nks(ik2) + 1
        dostmp(1:nx,fbep:lbep,ik2) = dostmp(1:nx,fbep:lbep,ik2) + dos(1:nx,fbep:lbep,ik,1)
        !
     end do
     !
  end do
  !
  do ik = 1, nk
     dos(1:nx,fbep:lbep,ik,1) = dostmp(1:nx,fbep:lbep,ik) / dble(nks(ik))
  end do
  !
  ! 2
  !
  nks(1:nk) = 0
  dostmp(1:nx,fbep:lbep,1:nk) = 0d0
  !
  do ik = 1, nk
     !
     kv1(1:3) = grid(1:3,ik) + 0.5d0 / dble(ng(1:3))
     !
     do isym = 1, nsym
        !
        kv2(1:3) = matmul(dble(sym(1:3,1:3,isym)), kv1(1:3)) * dble(ng(1:3))
        kv2(1:3) = kv2(1:3) - 0.5d0
        ikv2(1:3) = nint(kv2(1:3))
        !
        if(any(abs(kv2(1:3) - dble(ikv2(1:3))) > 1d-8)) cycle
        !
        ikv2(1:3) = modulo(ikv2(1:3), ng(1:3))
        ik2 = 1 + ikv2(1) + ng(1) * ikv2(2) + ng(1) * ng(2) * ikv2(3)
        !
        nks(ik2) = nks(ik2) + 1
        dostmp(1:nx,fbep:lbep,ik2) = dostmp(1:nx,fbep:lbep,ik2) + dos(1:nx,fbep:lbep,ik,2)
        !
     end do
     !
  end do
  !
  do ik = 1, nk
     dos(1:nx,fbep:lbep,ik,2) = dostmp(1:nx,fbep:lbep,ik) / dble(nks(ik))
  end do
  !
end subroutine symm_dosk
!
! Average matrix in grid
!
subroutine average_matrix()
  !
  use omp_lib
  use mpi
  use scdft_vals, only : nk0, nk, nm, fstk, lstk, ncf, nsym, fbep, lbep, nbep, &
  &                  vc0, gg0, omg0, vc, gg, omg, petot, my_rank, fbee, lbee, nbee
  !use scdft_routines, only : cnt_and_dsp, average_matrix_symmetrize
  !
  integer :: iq, ik, jk, ik2, ii, cnt, dsp, cntmax, ierr, ipe, &
  &          nqs(nk0, fstk:lstk), iks(2, nsym, nk0, fstk:lstk)
  real(8) :: t1, t2
  integer,allocatable :: qindx(:)
  real(8),allocatable :: gg1(:,:,:,:,:), vc1(:,:,:,:,:), omg1(:,:)
  !
  t1 = MPI_WTIME()
  allocate(vc(ncf,fbee:lbee,nk,fbee:lbee,fstk:lstk), &
  &        gg(nm,fbep:lbep,nk,fbep:lbep,fstk:lstk), omg(nm,nk,fstk:lstk))
  !
  call cnt_and_dsp(nk0,cnt,dsp)
  !
  cntmax = cnt
  call MPI_allREDUCE(MPI_IN_PLACE, cntmax, 1, MPI_INTEGER, &
  &                  MPI_MAX, MPI_COMM_WORLD,ierr)
  !
  allocate(gg1(nm,fbep:lbep,fbep:lbep,nk,cntmax), Vc1(ncf,fbee:lbee,fbee:lbee,nk,cntmax), &
  &       omg1(nm,cntmax), qindx(cntmax))
  !
  ! #####  Symmetrize  #####
  !
  call average_matrix_symmetrize(nqs,iks)
  !
  t2 = MPI_WTIME()
  if(my_rank == 0) write(*,*) "  Time (symmetrize) : ", t2 - t1, " sec." 
  t1 = MPI_WTIME()
  !
  ! #####  Communicate  #####
  !
!  Vc(1:ncf,fbee:lbee,1:nk,fbee:lbee,fstk:lstk) = 0d0
  gg( 1:nm,fbep:lbep,1:nk,fbep:lbep,fstk:lstk) = 0d0
  omg(1:nm,          1:nk,          fstk:lstk) = 0d0
  do ipe = 0, petot - 1
     !
!     Vc1(1:ncf,fbee:lbee,fbee:lbee,1:nk,1:cntmax) = 0d0
     gg1( 1:nm,fbep:lbep,fbep:lbep,1:nk,1:cntmax) = 0d0
     omg1(1:nm,                         1:cntmax) = 0d0
     qindx(1:cntmax) = 1
     !
     if(ipe == my_rank) then
        !
!        Vc1(1:ncf,fbee:lbee,fbee:lbee,1:nk,1:cnt) = Vc0(1:ncf,fbee:lbee,fbee:lbee,1:nk,1:cnt)
        gg1( 1:nm,fbep:lbep,fbep:lbep,1:nk,1:cnt) = gg0( 1:nm,fbep:lbep,fbep:lbep,1:nk,1:cnt)
        omg1(1:nm,                         1:cnt) = omg0(1:nm,                         1:cnt)
        !        
        do iq = 1, cnt
           qindx(iq) = dsp + iq
        end do
        !
     end if
     !
!     call MPI_BCAST(Vc1, ncf * nbee * nbee * nk * cntmax, MPI_DOUBLE_PRECISION, &
!     &              ipe, MPI_COMM_WORLD, ierr)
     !
     call MPI_BCAST(gg1, nm * nbep * nbep * nk * cntmax, MPI_DOUBLE_PRECISION, &
     &              ipe, MPI_COMM_WORLD, ierr)
     !
     call MPI_BCAST(omg1, nm * cntmax, MPI_DOUBLE_PRECISION, &
     &              ipe, MPI_COMM_WORLD, ierr)
     !
     call MPI_BCAST(qindx, cntmax, MPI_INTEGER, &
     &              ipe, MPI_COMM_WORLD, ierr)
     !
     !$OMP PARALLEL DEFAULT(NONE) &
     !$OMP & SHARED(nqs,iks,fstk,lstk,cntmax,qindx,ncf,nm, &
     !$OMP &        fbep, lbep, Vc1,Vc,omg1,omg,gg1,gg,fbee,lbee) &
     !$OMP & PRIVATE(ik,jk,iq,ik2,ii)
     !
     !$OMP DO
     do ik = fstk, lstk
        do iq = 1, cntmax
           !
           do ii = 1, nqs(qindx(iq),ik)
              !
              jk  = iks(1, ii, qindx(iq), ik)
              ik2 = iks(2, ii, qindx(iq), ik)
              !
!              Vc(1:ncf,fbee:lbee,jk,fbee:lbee,ik) = Vc(1:ncf,fbee:lbee,jk,fbee:lbee,ik) &
!              &                   + Vc1(1:ncf,fbee:lbee,fbee:lbee,ik2,iq)
              !
              gg(1:nm,fbep:lbep,jk,fbep:lbep,ik) = gg(1:nm,fbep:lbep,jk,fbep:lbep,ik) &
              &        + gg1(1:nm,fbep:lbep,fbep:lbep,ik2,iq)
              !
              omg(1:nm,jk,ik) = omg(1:nm,jk,ik) &
              &        + omg1(1:nm,iq)
              !
           end do ! ii
           !
        end do ! ik2
     end do ! iq
     !$OMP END DO
     !$OMP END PARALLEL
     !
  end do ! ipe
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(fstk,lstk,ncf,nm,nk0,nk,Vc,gg,omg,iks,nsym,fbep,lbep,fbee,lbee) &
  !$OMP & PRIVATE(ik,jk)
  !
  !$OMP DO
  do ik = fstk, lstk
     do jk = 1, nk
!        Vc(1:ncf,fbee:lbee,jk,fbee:lbee,ik) = Vc(1:ncf,fbee:lbee,jk,fbee:lbee,ik) &
!        &                         / dble(count(iks(1,1:nsym,1:nk0,ik) == jk))
        gg( 1:nm,fbep:lbep,jk,fbep:lbep,ik) = gg( 1:nm,fbep:lbep,jk,fbep:lbep,ik) &
        &                                   / dble(count(iks(1,1:nsym,1:nk0,ik) == jk))
        omg(1:nm,jk,ik) = omg(1:nm,jk,ik) &
        &       / dble(count(iks(1,1:nsym,1:nk0,ik) == jk))
     end do ! jk
  end do ! ik
  !$OMP END DO
  !$OMP END PARALLEL
  !
!  deallocate(gg0, vc0, omg0, gg1, vc1, omg1)
  deallocate(gg0, omg0, gg1, omg1)
  !
end subroutine average_matrix
!
! Symmetrize
!
subroutine average_matrix_symmetrize(nqs,iks)
  !
  use mpi
  use omp_lib
  use scdft_vals, only : nk, nk0, fstk, lstk, nsym, grid, rkv, ng, sym, nsym
  !
  integer,intent(out) :: nqs(nk0,fstk:lstk)
  integer,intent(out) :: iks(2,nsym,nk0,fstk:lstk)
  !
  integer :: iq, ik, isym, ikv1(3), ik2, jk2
  real(8) :: kv1(3)
  !
  nqs(           1:nk0,fstk:lstk) = 0
  iks(1:2,1:nsym,1:nk0,fstk:lstk) = 0
  !
  do iq = 1, nk0
     !
     do ik = 1, nk
        !
        do isym = 1, nsym
           !
           kv1(1:3) = matmul(dble(sym(1:3,1:3,isym)), grid(1:3,ik)) * dble(ng(1:3))
           ikv1(1:3) = nint(kv1(1:3))
           !
           if(any(abs(kv1(1:3) - dble(ikv1(1:3))) > 1d-5)) cycle
           !
           ikv1(1:3) = modulo(ikv1(1:3), ng(1:3))
           ik2 = 1 + ikv1(1) + ng(1) * ikv1(2) + ng(1) * ng(2) * ikv1(3)
           !
           if(ik2 < fstk .or. lstk < ik2) cycle
           !
           kv1(1:3) = matmul(sym(1:3,1:3,isym), grid(1:3,ik) + rkv(1:3,iq))
           kv1(1:3) = kv1(1:3) * dble(ng(1:3)) - 0.5d0
           ikv1(1:3) = nint(kv1(1:3))
           !
           if(any(abs(kv1(1:3) - dble(ikv1(1:3))) > 1d-5)) cycle
           !
           ikv1(1:3) = modulo(ikv1(1:3), ng(1:3))
           jk2 = 1 + ikv1(1) + ng(1) * ikv1(2) + ng(1) * ng(2) * ikv1(3)
           !
           nqs(iq,ik2) = nqs(iq,ik2) + 1
           iks(1:2,nqs(iq,ik2),iq,ik2) = (/jk2, ik/)
           !
        end do ! jk
        !
     end do ! ik
     !
  end do ! isym
  !
end subroutine average_matrix_symmetrize
!
! Calc Z_{n k}
!
subroutine make_Z()
  !
  use mpi
  use scdft_vals, only : nt, nt1, nt2, nm, Z, xi, dk, kindx, bindx, gg, omg, beta, &
  &                      fbep, lbep, emin, my_rank
  !use scdft_routines, only : Zweight, cnt_and_dsp
  !
  integer :: it, ik, ib, jt, jk, jb, im, cnt, dsp, ierr
  real(8) :: x, xp, om, zave(2), ave(2), tx, txp, tom
  !
  allocate(Z(nt,2))
  Z(1:nt,1:2) = 0d0
  !
  call cnt_and_dsp(nt1,cnt,dsp)
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(cnt,dsp,nt2,nm,kindx,bindx,beta, &
  !$OMP &        fbep,lbep,xi,dk,Z,gg,omg) &
  !$OMP & PRIVATE(it,ik,ib,jt,jk,jb,im,x,xp,om,tx,txp,tom)
  !
  !$OMP DO REDUCTION(+: Z)
  do it = dsp + 1, dsp + cnt
     !
     x = abs(xi(it,1) * beta * 0.5d0)
     tx = tanh(x)
     ik = kindx(it,1)
     ib = bindx(it,1)
     !
     if(ib < fbep .or. lbep < ib) cycle
     !
     do jt = 1, nt2
        !
        xp = abs(xi(jt,2) * beta * 0.5d0)
        txp = tanh(xp)
        jk = kindx(jt,2)
        jb = bindx(jt,2)
        !
        if(jb < fbep .or. lbep < jb) cycle
        !
        do im = 1, nm
           !
           om = abs(omg(im,jk,ik) * beta * 0.5d0)
           tom = tanh(om)
           !
           Z(it,1) = Z(it,1) + dk(jt,2) * gg(im,jb,jk,ib,ik) * beta**2 &
           &                 * Zweight(x, xp, om, tx, txp, tom)
           Z(jt,2) = Z(jt,2) + dk(it,1) * gg(im,jb,jk,ib,ik) * beta**2 &
           &                 * Zweight(xp, x, om, txp, tx, tom)
           !
        end do ! im
        !
     enddo ! jt
     !
  enddo ! it
  !$OMP END DO
  !
  !$OMP END PARALLEL
  !
  Z(1:nt,1:2) = - Z(1:nt,1:2)
  !
  call MPI_allREDUCE(MPI_IN_PLACE, Z, 2 * nt, MPI_DOUBLE_PRECISION, &
  &                  MPI_SUM,MPI_COMM_WORLD,ierr)
  !
  zave(1) = sum(Z(1:nt,1) * dk(1:nt,1), abs(xi(1:nt,1)) < emin * 1d-2)
  zave(2) = sum(Z(1:nt,2) * dk(1:nt,2), abs(xi(1:nt,2)) < emin * 1d-2)
  ave(1) = sum(dk(1:nt,1), abs(xi(1:nt,1)) < emin * 1d-2)
  ave(2) = sum(dk(1:nt,2), abs(xi(1:nt,2)) < emin * 1d-2)
  !
  zave(1:2) = zave(1:2) / ave(1:2)
  !
  if(my_rank == 0)  write(*,*) "  Averaged Z_{FS} : ", zave(1:2)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! modified by Kentaro Tsutsumi !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(my_rank == 0) then
    open(223, file = "Znk_ttm.dat")
    do it = 1, nt
      write(223, *) it, Z(it, 1), Z(it, 2)
    end do
    close(223)
  end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! modified by Kentaro Tsutsumi !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
! Make Effective interaction
!
subroutine make_effint()
  !
  use mpi
  use scdft_vals, only : pi, my_rank, nt, nt1, nt2, nm, ncf, kindx, bindx, beta, &
  &                  xi, dk, Z, gg, omg, Vc, effint, fbep, nmf, lbep, mf, wmf
  !use scdft_routines, only : cnt_and_dsp, Kweight, calc_Kel, weightspoints_gl
  !
  integer :: it, ik, ib, jt, jk, jb, im, cnt, dsp
  real(8) :: x, xp, om, Kph, Kel, bx, bxp, tx, txp, tom
  !
  call cnt_and_dsp(nt1,cnt,dsp)
  allocate(effint(nt2,dsp + 1:dsp + cnt))
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
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(cnt, dsp, nt1, nt2, nm, ncf, kindx, bindx, beta, &
  !$OMP &        xi, dk, omg, gg, Vc, effint, fbep, lbep) &
  !$OMP & PRIVATE(it, jt, ik, jk, ib, jb, im, x, xp, om, Kph, Kel, &
  !$OMP &         bx, bxp, tx, txp, tom)
  !
  !$OMP DO
  do it = dsp + 1, dsp + cnt
     !
     x = xi(it,1)
     bx = 0.5d0 * beta * x
     tx = tanh(bx)
     ik = kindx(it,1)
     ib = bindx(it,1)
     !
     do jt = 1, nt2
        !
        xp = xi(jt,2)
        bxp = 0.5d0 * beta * xp
        txp = tanh(bxp)
        jk = kindx(jt,2)
        jb = bindx(jt,2)
        !
        Kph = 0d0
        !
        if(all((/fbep <= ib, ib <= lbep, fbep <= jb, jb <= lbep/))) then
           !
           do im = 1, nm
              !
              om = omg(im,jk,ik) * beta * 0.5d0
              tom = tanh(om)
              !
              Kph = Kph + gg(im,jb,jk,ib,ik) * beta &
              &         * Kweight(bx, bxp, om, tx, txp, tom)
              !       
           end do ! im
           Kph = Kph * 2d0
           !
        end if
        !
        Kel = calc_Kel(x,xp,Vc(1:ncf,jb,jk,ib,ik))
        !
        effint(jt,it) = Kph + Kel
        !
     enddo ! jt
     !
  enddo ! it
  !$OMP END DO
  !$OMP END PARALLEL
  !
  deallocate(gg, vc, omg)
  !
end subroutine make_effint
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
  use scdft_vals, only : pi, ncf, nmf, mf, wmf
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
  xcf = acos(- x0)
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
! Solve gap equation
!
subroutine gapeq()
  !
  use mpi
  use omp_lib
  use scdft_vals, only : my_rank, nt, nt1, nt2, itrmax, er0, xi, delta, emin, dk
  !use scdft_routines, only : make_rhs, ddot, daxpy, out_delta, cnt_and_dsp
  !
  integer :: itr, jtr, cnt, dsp, ierr
  real(8) :: res, alpha = 0.2d0, dd(nt,2), rhs(nt,2), rhs0(nt,2), drhs(nt,2), &
  &          jacob1(nt,2,itrmax), jacob2(nt,2,itrmax), ave(2), dave(2), t2, t1
  logical :: lfermi(nt)
  !
  call cnt_and_dsp(nt,cnt,dsp)
  !
  if(my_rank == 0) then
     write(*,*) '  Convergense threshold : ', er0
     write(*,*) '           Max itration : ',itrmax
     write(*,*) ""
  end if
  !
  t1 = OMP_GET_WTIME()
  !
  itr = 0
  if(my_rank == 0) write(*,*) "  Iteration ", itr
  !
  call make_rhs(rhs)
  res = ddot(nt * 2, rhs, 1, rhs, 1)
  res = sqrt(res) / dble(nt * 2)
  !
  ! $$$$$  Information of conversience  $$$$$
  !
  t2 = OMP_GET_WTIME() - t1
  t1 = OMP_GET_WTIME()
  !
  ave(1:2) = 0d0
  dave(1:2) = 0d0
  !
  lfermi(1:nt1) = abs(xi(1:nt1,1)) < emin * 1d-2
  dave(1) = sum(delta(1:nt1,1) * dk(1:nt1,1), lfermi(1:nt1))
  ave(1)  = sum(                 dk(1:nt1,1), lfermi(1:nt1))
  !
  lfermi(1:nt2) = abs(xi(1:nt2,2)) < emin * 1d-2
  dave(2) = sum(delta(1:nt2,2) * dk(1:nt2,2), lfermi(1:nt2))
  ave(2)  = sum(                 dk(1:nt2,2), lfermi(1:nt2))
  !
  dave(1:2) = dave(1:2) / ave(1:2) * 13.6057d3
  !
  if(my_rank == 0) then
     write(*,*) "      Residual[Ry] : ", res
     write(*,*) "       Delta [meV] : ", dave(1:2)
     write(*,*) "    Rap time [sec] : ", t2
     write(*,*) ""
  end if ! (my_rank == 0)
  !
  ! $$$$$  End information of conversience  $$$$$
  !
  if(res < er0) goto 5
  !
  dd(1:nt,1:2) = - alpha * rhs(1:nt,1:2)
  !
  do itr = 1, itrmax
     !
     if(my_rank == 0) write(*,*) "  Iteration ", itr
     !
     delta(1:nt,1:2) = delta(1:nt,1:2) + dd(1:nt,1:2)
     !
     rhs0(1:nt,1:2) = rhs(1:nt,1:2)
     call make_rhs(rhs)
     res = ddot(nt, rhs, 1, rhs, 1)
     res = sqrt(res) / dble(nt)
     !
     ! $$$$$  Information of conversience  $$$$$
     !
     t2 = OMP_GET_WTIME() - t1
     t1 = OMP_GET_WTIME()
     !
     ave(1:2) = 0d0
     dave(1:2) = 0d0
     !
     lfermi(1:nt1) = abs(xi(1:nt1,1)) < emin * 1d-2
     dave(1) = sum(delta(1:nt1,1) * dk(1:nt1,1), lfermi(1:nt1))
     ave(1)  = sum(                 dk(1:nt1,1), lfermi(1:nt1))
     !
     lfermi(1:nt2) = abs(xi(1:nt2,2)) < emin * 1d-2
     dave(2) = sum(delta(1:nt2,2) * dk(1:nt2,2), lfermi(1:nt2))
     ave(2)  = sum(                 dk(1:nt2,2), lfermi(1:nt2))
     !
     dave(1:2) = dave(1:2) / ave(1:2) * 13.6057d3
     !
     if(my_rank == 0) then
        write(*,*) "      Residual[Ry] : ", res
        write(*,*) "       Delta [meV] : ", dave(1:2)
        write(*,*) "    Rap time [sec] : ", t2
        write(*,*) ""
     end if ! (my_rank == 0)
     !
     ! $$$$$  End information of conversience  $$$$$
     !
     if(res < er0) then
        !       
        delta(1:nt,1) = delta(1:nt,1) * sign(1d0, dave(1))
        delta(1:nt,2) = delta(1:nt,2) * sign(1d0, dave(2))
        !
        goto 5
        !
     end if
     !
     if(modulo(itr, 2) == 0) then
        call out_delta("delta.evn")
     else
        call out_delta("delta.odd")
     end if
     !
     ! Update Jacobian with drhs
     !
     drhs(1:nt,1:2) = rhs(1:nt,1:2) - rhs0(1:nt,1:2)
     !
     jacob1(1:nt,1:2,itr) = - alpha * drhs(1:nt,1:2)
     do jtr = 1, itr - 1
        jacob1(1:nt,1:2,itr) = jacob1(1:nt,1:2,itr) - jacob1(1:nt,1:2,jtr) &
        &          * ddot(nt, jacob2(1:nt,1:2,jtr), 1, drhs(1:nt,1:2), 1)
     end do
     jacob1(1:nt,1:2,itr) = dd(1:nt,1:2) + jacob1(1:nt,1:2,itr)
     jacob2(1:nt,1:2,itr) = drhs(1:nt,1:2) / ddot(nt, drhs(1:nt,1:2), 1, drhs(1:nt,1:2), 1)
     !
     ! Compute dd with new Jacobian & rhs
     !
     dd(1:nt,1:2) = - alpha * rhs(1:nt,1:2)
     do jtr = 1, itr
        dd(1:nt,1:2) = dd(1:nt,1:2) - jacob1(1:nt,1:2,jtr) &
        &        * ddot(nt, jacob2(1:nt,1:2,jtr), 1, rhs(1:nt,1:2), 1)
     end do
     !
  end do ! itr
  !
  call out_delta("delta.dat")
  !
  if(my_rank == 0) write(*,*) ""
  if(my_rank == 0) write(*,*) '  Not converged! res = ',res
  return
  !
5 continue
  !
  call out_delta("delta.dat")
  !
  if(my_rank == 0) write(*,*) ""
  if(my_rank == 0) write(*,*) '  Converged! iter = ',itr
  !
end subroutine gapeq
!
! Make RHS vector for linear probrem.
!
subroutine make_rhs(veco)
  !
  use mpi
  use scdft_vals, only : nt, nt1, nt2, beta, nt, xi, delta, effint, Z, dk, &
  &                  nth, xic, my_rank
  !use scdft_routines, only : cnt_and_dsp_omp, dgemv
  !
  real(8),intent(out) :: veco(nt,2)
  !
  integer :: it, jt, cnt, dsp, cnt1, dsp1, nh
  real(8) :: chi(nt,2), dosh, Kh, dlth, xmax
  !
  call cnt_and_dsp(nt1, cnt, dsp)
  !
  veco(1:nt,1:2) = 0d0
  chi(1:nt,1:2) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nt1, nt2, cnt, dsp, beta, xi, delta, veco, effint, chi, dk, nth) &
  !$OMP & PRIVATE(it, jt,cnt1, dsp1)
  !
  !$OMP DO
  do it =1, nt1
     chi(it,1) = 0.5d0 * dk(it,1) * delta(it,1) &
     &         * tanh(0.5d0 * beta * sqrt(xi(it,1)**2 + delta(it,1)**2)) &
     &         / sqrt(xi(it,1)**2 + delta(it,1)**2)
  end do
  !$OMP END DO
  !
  !$OMP DO
  do it = 1, nt2
     chi(it,2) = 0.5d0 * dk(it,2) * delta(it,2) &
     &         * tanh(0.5d0 * beta * sqrt(xi(it,2)**2 + delta(it,2)**2)) &
     &         / sqrt(xi(it,2)**2 + delta(it,2)**2)
  end do
  !$OMP END DO
  !
  call cnt_and_dsp_omp(cnt, cnt1,dsp1)
  call dgemv("T", nt2, cnt1, 1d0, effint(1:nt2,dsp + dsp1 + 1:dsp + dsp1 + cnt1), nt2, &
  &    chi(1:nt2,2), 1, 1d0, veco(dsp + dsp1 + 1:dsp + dsp1 + cnt1,1), 1)
  !
  !$OMP DO REDUCTION(+: veco)
  do it = 1, nth
     call dgemv("N", nt2, cnt1, 1d0, effint(1:nt2,dsp+dsp1+1:dsp+dsp1+cnt1), nt2, &
     &    chi(dsp+dsp1+1:dsp+dsp1+cnt1,1), 1, 1d0, veco(1:nt2,2), 1)
  end do
  !$OMP END DO
  !
  !$OMP END PARALLEL
  !
  ! High energy region
  !
  if(xic > 0) then
     !
     ! 1, 
     !
     nh = count(xi(1:nt1,1) > xic)
     dlth = sum(delta(1:nt1,1) / xi(1:nt1,1), xi(1:nt1,1) > xic) &
     &    / sum(1d0 / xi(1:nt1,1)**2,         xi(1:nt1,1) > xic)
     !
     xmax = maxval(xi(1:nt1,1))
     dosh = sum(dk(1:nt1,1), xi(1:nt1,1) > xic) / (xmax - xic)
     !
     if(my_rank == 0) then
        write(*,*) "    Extrapol 1 : ", dlth / xic * 13.6d3
     end if ! (my_rank == 0)
     !
     !$OMP PARALLEL DEFAULT(NONE) &
     !$OMP & SHARED(veco,nt2,effint,dsp,cnt,xi,xic,xmax,nh,dosh,dlth) &
     !$OMP & PRIVATE(it,Kh)
     !$OMP DO
     do it = 1, nt2
        Kh = sum(effint(it,dsp + 1:dsp + cnt), xi(dsp + 1:dsp + cnt,1) > xic) / dble(nh)
        veco(it,2) = veco(it,2) + 0.5d0 * dosh * Kh * dlth / xmax
     end do
     !$OMP END DO
     !$OMP END PARALLEL
     !
     ! 2, 
     !
     nh = count(xi(1:nt2,2) > xic)
     dlth = sum(delta(1:nt2,2) / xi(1:nt2,2), xi(1:nt2,2) > xic) &
     &    / sum(1d0 / xi(1:nt2,2)**2,         xi(1:nt2,2) > xic)
     !
     xmax = maxval(xi(1:nt2,2))
     dosh = sum(dk(1:nt2,2), xi(1:nt2,2) > xic) / (xmax - xic)
     !
     if(my_rank == 0) then
        write(*,*) "    Extrapol 2 : ", dlth / xic * 13.6d3
     end if ! (my_rank == 0)
     !
     !$OMP PARALLEL DEFAULT(NONE) &
     !$OMP & SHARED(veco,nt2,effint,dsp,cnt,xi,xic,xmax,nh,dosh,dlth) &
     !$OMP & PRIVATE(it,Kh)
     !$OMP DO
     do it = dsp + 1, dsp + cnt
        Kh = sum(effint(1:nt2,it), xi(1:nt2,2) > xic) / dble(nh)
        veco(it,2) = veco(it,2) + 0.5d0 * dosh * Kh * dlth / xmax
     end do
     !$OMP END DO
     !$OMP END PARALLEL
     !
  end if ! (xic > 0)
  !
  call MPI_allREDUCE(MPI_IN_PLACE, veco, nt * 2, MPI_DOUBLE_PRECISION, &
  &                  MPI_SUM,MPI_COMM_WORLD,it)
  !
  veco(1:nt,1:2) = delta(1:nt,1:2) + veco(1:nt,1:2) / (1d0 + Z(1:nt,1:2)) 
  !
end subroutine make_rhs
!
! Output to file (delta)
!
subroutine out_delta(fname)
  !
  use scdft_vals, only : nt1, nt2, xi, delta, dk, kindx, bindx, Z
  !
  character(*),intent(in) :: fname
  !
  integer :: it, fo = 20
  !
  open(fo, file = fname)
  !
  write(fo,*) "#", nt2, nt1
  !
  write(fo,*) ""
  !
  do it = 1, nt2
     write(fo,'(4e25.15,2i8)') xi(it,2), delta(it,2), Z(it,2), dk(it,2), &
     &                         kindx(it,2), bindx(it,2)
  enddo
  !
  write(fo,*) ""
  !
  do it = 1, nt1
     write(fo,'(4e25.15,2i8)') xi(it,1), delta(it,1), Z(it,1), dk(it,1), &
     &                       kindx(it,1), bindx(it,1)
  enddo
  !
  close(fo)
  !
end subroutine out_delta
!
end module scdft_routines
!
! Main routine
!
program scdft
  !
  use scdft_vals, only : petot, my_rank, nth
  use scdft_routines, only : read_stdin, read_dat, read_elph, average_matrix, read_Coulomb, &
  &                          ini_delta, make_Z, make_effint, gapeq
  use mpi
  use omp_lib
  !
  implicit none
  !
  integer :: ierr, date(8)
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
  !$OMP MASTER
  nth = OMP_GET_NUM_THREADS()
  !$OMP END MASTER
  !$OMP END PARALLEL
  if(my_rank == 0) write(*,*) '  # of thread : ', nth
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
!  call read_Coulomb()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,*) "  Time : ", t2 - t1, " sec"
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  Set or read initial delta  #####"
     write(*,*)
  end if
  t1 = OMP_GET_WTIME()
  call ini_delta()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,*) "  Time : ", t2 - t1, " sec"
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
     write(*,*) "#####  Compute effective interaction  #####"
     write(*,*)
  end if
  t1 = OMP_GET_WTIME()
!  call make_effint()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,*) "  Time : ", t2 - t1, " sec"
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  Solve gap equation  #####"
     write(*,*)
  end if
  t1 = OMP_GET_WTIME()
!  call gapeq()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,*) "  Time : ", t2 - t1, " sec"
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
end program


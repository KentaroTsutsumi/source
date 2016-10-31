module rpa_el_vals
  !
  use,intrinsic :: iso_c_binding
  !
  real(8),parameter :: pi = acos(-1d0)
  !
  integer,save :: &
  & iqdo,          &
  & nkpe,          &
  & ngv,           & ! # of G-vector in wcut
  & gminxc(3),     & ! Min. G of dmxc
  & gmaxxc(3),     & ! Max. G of dmxc
  & laddxc,        & ! Switch for XC
  & ltetra,        & ! Switch for k integration scheme
  & ivvec(4,20,6), & ! points for tetrahedron method
  & nmf,           & ! # of Matsubara freq.
  & petot,         & ! Total # of PEs
  & my_rank,       & ! Index of PE
  & nsym,          & ! # of symmetries
  & ng(3),         & ! k point grid
  & nk,            & ! ng(1)*ng(2)*ng(3). Total # of k points
  & ngd(3),        & ! Dense k point grid
  & nkd,           & ! ngd(1)*ngd(2)*ngd(3). Total # of k points
  & nk0,           & ! # of irreducible k points
  & nb,            & ! # of bands
  & nf(3),         & ! # grid for FFT
  & nftot,         & ! Total # of G
  & igmin(3),      & ! Min. indices of G vector
  & npwmax           ! Max. # of PWs
  !
  real(8),save :: &
  & wcut,       & ! Cutoff energy for chi [Ry]
  & wlsm(4,20), & ! Weight for tetrahedron method
  & cellV,      & ! Unit cell volume [a.u.^3]
  & alat,       & ! Lattice constant [a.u.]
  & bvec(3,3)     ! Reciplocal latticee vectors [/a.u.]
  !
  integer,allocatable,save :: &
  & gindx(:),   & ! (nftot) G-vector in wcut
  & sym(:,:,:), & ! (3,3,nsym). Symmetry operators
  & iqv(:,:),   & ! (3,nk0). Irreducible k points
  & grid(:,:),  & ! (3,nk). All k points
  & npw(:,:),   & ! (nk,2). # of PWs
  & igv(:,:,:,:)  ! (3,npwmax,nk,2). G points
  !
  real(8),allocatable,save :: &
  & eig(:,:),       & ! (nb,nkd) Kohn-Sham eigenvaluues
  & gq2(:),         & ! (nftot) |G+q|^2
  & mf(:),          & ! (nmf) Matsubara frequencies
  & wmf(:),         & ! (nmf) Weights for Matsubara frequency
  & Kel(:,:,:,:)      ! (nmf,nb,nb,nk). Coulomb Kernel
  !
  complex(8),allocatable,save :: &
  & dmxc(:,:,:),    & ! (gminxc(1):gmaxxc(1), ...) derivative of XC potential
  & wscr(:,:,:),    & ! (nftot,nftot,0:nmf) Screened interaction
  & wfc1(:,:,:),    & ! (nftot,nb,nk) wfc(ib,ik)
  & wfc2(:,:,:),    & ! (nftot,nb,nk) wfc(jb,jk)
  & wfc(:,:,:,:)      ! (npwmax,nb,nk,2) wfc(G)
  !
  include 'fftw3.f'
  !
end module rpa_el_vals
!
module rpa_el_routines
  !
  implicit none
  !
  ! For blas
  !
  INTERFACE
     !
     double complex function zdotc(n,zx,incx,zy,incy)
       double complex zx(*),zy(*)
       integer        incx,incy,n
     END function zdotc
     !
     subroutine  dscal(n,da,dx,incx)
       double precision da,dx(*)
       integer          incx,n
     END subroutine dscal
     !
     subroutine  zscal(n,za,zx,incx)
       double complex za,zx(*)
       integer        incx,n
     END subroutine zscal
     !
     subroutine  dcopy(n,dx,incx,dy,incy)
       double precision dx(*),dy(*)
       integer          incx,incy,n
     END subroutine dcopy
     !
     subroutine  zcopy(n,zx,incx,zy,incy)
       double complex zx(*),zy(*)
       integer        incx,incy,n
     END subroutine zcopy
     !
     subroutine daxpy(n,da,dx,incx,dy,incy)
       double precision dx(*),dy(*),da
       integer          incx,incy,n
     END subroutine daxpy
     !
     subroutine zaxpy(n,za,zx,incx,zy,incy)
       double complex zx(*),zy(*),za
       integer        incx,incy,n
     END subroutine zaxpy
     !
     subroutine zgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
       character*1        transa, transb
       integer            m, n, k, lda, ldb, ldc
       complex*16         alpha, beta
       complex*16         a( lda, * ), b( ldb, * ), c( ldc, * )
     END subroutine zgemm
     !
  END INTERFACE
  !
contains
  !
 ! Read from STDIN
!
subroutine stdin()
  !
  use mpi, only : MPI_INTEGER, MPI_COMM_WORLD, MPI_DOUBLE_PRECISION
  use rpa_el_vals, only : my_rank, ng, iqdo, nmf, ngd, nkd, ltetra, laddxc, wcut
  !
  integer :: ierr
  namelist /input/ ng, iqdo, nmf, ngd, ltetra, laddxc, wcut
  !
  if(my_rank == 0) then
     !
     wcut = 0d0
     ltetra = 0
     laddxc = 0
     !
     read(*,input,err=100)
     write(*,*) "                 k grid : ", ng(1:3)
     write(*,*) "           Dense k grid : ", ngd(1:3)
     write(*,*) "   # of Matsubara freq. : ", nmf
     write(*,*) "                q index : ", iqdo
     write(*,*) "                 ltetra : ", ltetra
     write(*,*) "                 laddxc : ", laddxc
     write(*,*) "  Cutoff kinetic energy : ", wcut
     !
  end if
  !
  call MPI_BCAST(ng,     3, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ngd,    3, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(iqdo,   1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nmf,    1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ltetra, 1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(laddxc, 1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(wcut,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  !
  nkd = product(ngd(1:3))
  !
  return
  !
100 write(*,*) "Stop in stdin. reading namelist file"
  !
  call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  call MPI_FINALIZE(ierr)
  stop
  !
end subroutine stdin
!
! Read from data-file.xml
!
subroutine read_file()
  !
  use mpi, only : MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_INTEGER
  use iotk_module, only : iotk_open_read, iotk_scan_begin, iotk_scan_end, &
  &                       iotk_scan_dat, iotk_close_read
  use rpa_el_vals, only : my_rank, alat, bvec, nb, nk, ng, npwmax, &
  &                  cellV, nsym, sym, pi, eig, nkd
  !
  integer :: fi = 10, ierr
  real(8) :: avec(3,3)
  !
  if(my_rank == 0) then
     !
     call iotk_open_read(fi, file = 'data-file.xml')
     !
     ! Read lattice parameter
     !
     call iotk_scan_begin(fi,"CELL")
     !
     call iotk_scan_dat(fi,"LATTICE_PARAMETER",alat)
     write(*,*) "  lattice constant[a.u] : ", alat
     call iotk_scan_begin(fi,"DIRECT_LATTICE_VECTORS")
     call iotk_scan_dat(fi,"a1", avec(1:3,1))
     call iotk_scan_dat(fi,"a2", avec(1:3,2))
     call iotk_scan_dat(fi,"a3", avec(1:3,3))
     !
     ! Calcurate unit-cell volume
     !
     cellV = avec(1,1) * (avec(2,2) * avec(3,3) - avec(3,2) * avec(2,3)) &
     &     + avec(2,1) * (avec(3,2) * avec(1,3) - avec(1,2) * avec(3,3)) &
     &     + avec(3,1) * (avec(1,2) * avec(2,3) - avec(2,2) * avec(1,3)) 
     cellV = abs(cellV)
     !
     write(*,*) "  Cell volume[a.u.^3] : ", cellV
     !
     avec(1:3,1:3) = avec(1:3,1:3) / alat
     write(*,*) "  Direct lattice vector[a] : "
     write(*,*) avec(1:3,1)
     write(*,*) avec(1:3,2)
     write(*,*) avec(1:3,3)
     write(*,*) ""
     call iotk_scan_end(fi,"DIRECT_LATTICE_VECTORS")
     !
     call iotk_scan_begin(fi,"RECIPROCAL_LATTICE_VECTORS")
     call iotk_scan_dat(fi,"b1",bvec(1:3,1))
     call iotk_scan_dat(fi,"b2",bvec(1:3,2))
     call iotk_scan_dat(fi,"b3",bvec(1:3,3))
     call iotk_scan_end(fi,"RECIPROCAL_LATTICE_VECTORS")
     write(*,*) "  Reciprocal lattice vector[2p/a] : "
     write(*,*) bvec(1:3,1)
     write(*,*) bvec(1:3,2)
     write(*,*) bvec(1:3,3)
     write(*,*) ""
     !
     call iotk_scan_end(fi,"CELL")
     !
     call iotk_scan_begin(fi,"PLANE_WAVES")
     call iotk_scan_dat(fi,"MAX_NUMBER_OF_GK-VECTORS",npwmax)
     call iotk_scan_end(fi,"PLANE_WAVES")
     !
     write(*,*) "  Max G-vector # : ", npwmax
     !
     ! Read # of band & k-vector
     !
     call iotk_scan_begin(fi,"BAND_STRUCTURE_INFO")
     !
     call iotk_scan_dat(fi,"NUMBER_OF_K-POINTS",nk)
     if(nk /= 2 * product(ng(1:3))) then
        write(*,*) "Stop in read_file. nk /= product(ng(1:3))"
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        call MPI_FINALIZE(ierr)
        stop
     end if
     nk = nk / 2
     write(*,*) "  # of k-vector : ", nk, " * 2"
     !
     call iotk_scan_dat(fi,"NUMBER_OF_BANDS",nb)
     write(*,*) "  Total # of bands : ", nb
     !
     call iotk_scan_end(fi,"BAND_STRUCTURE_INFO")
     !
     call iotk_close_read(fi)
     !
     ! Read symmetry
     !
     write(*,*) "    open symm.dat"
     open(fi, file = 'symm.dat')
     !
     read(fi,*) nsym
     !
     write(*,*) "  # of symmetries : ", nsym
     !
     allocate(sym(3,3,nsym))
     read(fi,*) sym(1:3,1:3,1:nsym)
     !
     close(fi)
     !
     ! Read eigenvalues
     !
     allocate(eig(nb,nkd))
     !
     open(fi, file = "eigval.dat")
     read(fi,*) eig(1:nb,1:nkd)
     close(fi)
     !
  end if
  !
  call MPI_BCAST(alat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(bvec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nk, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(npwmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(cellV, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nsym, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if(my_rank /= 0) allocate(sym(3,3,nsym))
  call MPI_BCAST(sym,  9 * nsym, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if(my_rank /= 0) allocate(eig(nb,nkd))
  call MPI_BCAST(eig, nb * nkd, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  !
end subroutine read_file
!
! Irreducible Brillouin zone
!
subroutine irr_bz()
  !
  use rpa_el_vals, only : my_rank, nk, ng, nk0, iqv, grid, nsym, sym
  !
  integer :: i1, i2, i3, ik, nkk, isym
  real(8) :: qv0(3,nk), qv1(3), qv2(3)
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
           grid(1:3,nkk) = 2 * ((/i1, i2, i3/) - 1)
           grid(1:3,nkk) = modulo(grid(1:3,nkk), 2 * ng(1:3))
           !
           qv1(1:3) = (dble((/i1, i2, i3/)) - 0.5d0) / dble(ng(1:3))
           qv1(1:3) = qv1(1:3) - dble(floor(qv1(1:3) + 0.5d0 + 1d-4))
           !
           do isym = 1, nsym
              !
              qv2(1:3) = matmul(dble(sym(1:3,1:3,isym)), qv1(1:3))
              qv2(1:3) = qv2(1:3) - dble(floor(qv2(1:3) + 0.5d0 + 1d-4))
              !
              do ik = 1, nk0
                 if(all(abs(qv2(1:3) - qv0(1:3,ik)) < 1d-8)) goto 10
              end do
              !
           end do
           !
           nk0 = nk0 + 1
           qv0(1:3,nk0) = qv1(1:3)
           !
10         continue
           !
        end do
     end do
  end do
  !
  allocate(iqv(3,nk0))
  !
  do ik = 1, nk0
     !
     qv0(1:3,ik) = qv0(1:3,ik) * dble(2 * ng(1:3))
     iqv(1:3,ik) = nint(qv0(1:3,ik))
     !
  end do
  !
  if(my_rank == 0) write(*,*) "  # of irreducible k points : ", nk0
  !
end subroutine irr_bz
!
! read g-vector & wfc
!
subroutine get_wfcg()
  !
  use mpi, only : MPI_IN_PLACE, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD
  use iotk_module, only : iotk_namlenx, iotk_open_read, iotk_scan_dat, iotk_close_read
  use rpa_el_vals, only : my_rank, pi, alat, bvec, nk, npwmax, nb, npw, igv, wfc, &
  &                  grid, ng, nkpe
  !use rpa_routines, only : cnt_and_dsp
  !
  integer :: fi = 10, ik, ib, ikv(3), cnt, dsp, ierr, ii
  real(8) :: kvec(3), kvec2(3)
  character(iotk_namlenx) :: attr
  !
  call cnt_and_dsp(nk,cnt,dsp)
  nkpe = cnt
  call MPI_allREDUCE(MPI_IN_PLACE, nkpe, 1, &
  &                  MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
  if(my_rank == 0) write(*,*) "  # of k par PE : ", nkpe
  !
  allocate(npw(          1:nkpe,2), &
  &        igv(3,npwmax, 1:nkpe,2), &
  &        wfc(npwmax,nb,1:nkpe,2)  )
  !
  npw(              1:nkpe,1:2) = 0
  igv(1:3,1:npwmax, 1:nkpe,1:2) = 0
  wfc(1:npwmax,1:nb,1:nkpe,1:2) = cmplx(0d0, 0d0)
  !
  ! Read Wfc(k)
  !
  do ii = 1, 2
     !
     do ik = 1, cnt
        !
        ! Read G vector file
        !
        write(attr,'(a,i5.5,a)') "K", nk * (ii - 1) + dsp + ik, "/gkvectors.dat"
        call iotk_open_read(fi, trim(attr), binary=.true.)
        call iotk_scan_dat(fi,"NUMBER_OF_GK-VECTORS",npw(ik,ii))
        call iotk_scan_dat(fi,"K-POINT_COORDS",kvec(1:3))
        call iotk_scan_dat(fi,"GRID",igv(1:3,1:npw(ik,ii),ik,ii))
        !
        call iotk_close_read(fi)
        !
        ikv(1:3) = modulo(grid(1:3,dsp + ik) + ii - 1 + ng(1:3), 2 * ng(1:3)) - ng(1:3)
        kvec2(1:3) = dble(ikv(1:3)) / dble(2 * ng(1:3))
        kvec2(1:3) = matmul(bvec(1:3,1:3), kvec2(1:3))
        !
        if(all(abs(kvec2(1:3) - kvec(1:3)) > 1d-8)) then
           write(*,*) "Stop in get_wfcg. ik = ", dsp + ik, " ii = ", ii
           write(*,*) nk * (ii - 1) + dsp + ik
           write(*,*) kvec(1:3)
           write(*,*) kvec2(1:3)
           write(*,*) grid(1:3,dsp + ik) + ii - 1
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           call MPI_FINALIZE(ierr)
           stop
        end if
        !
        ! Read WFC file
        !
        write(attr,'(a,i5.5,a)') "K", nk * (ii - 1) + dsp + ik, "/evc.dat"
        call iotk_open_read(fi, trim(attr), binary=.true.)
        !
        do ib = 1, nb
           !
           write(attr,*) ib
           write(attr,'(a,a)') "evc.", trim(adjustl(attr))
           call iotk_scan_dat(fi, trim(attr), wfc(1:npw(ik,ii),ib,ik,ii))
           !
        end do
        !
        call iotk_close_read(fi)
        !
     end do ! ik
     !
  end do ! ii
  !
  bvec(1:3,1:3) = bvec(1:3,1:3) * 2d0 * pi / alat
  !
end subroutine get_wfcg
!
! Compute cnt and dsp
!
subroutine cnt_and_dsp(n,cnt1,dsp1)
  !
  use rpa_el_vals, only : petot, my_rank
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
subroutine cnt_and_dsp_full(n,cnt,dsp)
  !
  use rpa_el_vals, only : petot, my_rank
  !
  integer,intent(in) :: n
  integer,intent(out) :: cnt(0:petot - 1), dsp(0:petot - 1)
  !
  integer :: ii
  !
  cnt(0:petot-1)        = n / petot
  cnt(0:mod(n,petot)-1) = n / petot + 1
  dsp(0) = 0
  do ii = 1, petot - 1
     dsp(ii) = dsp(ii - 1) + cnt(ii - 1)
  end do
  !
end subroutine cnt_and_dsp_full
!
! Make fft mesh
!
subroutine fft_mesh()
  !
  use mpi, only : MPI_IN_PLACE, MPI_INTEGER, MPI_MAX, MPI_MIN, MPI_COMM_WORLD
  use rpa_el_vals, only : my_rank, nk, npwmax, igv, nf, nftot, igmin, nkpe
  !use rpa_routines, only : fft_base
  !
  integer :: ii, itr, igmax(3), ierr
  logical :: lbase
  !
  do ii = 1, 3
     !
     igmax(ii) = maxval(igv(ii,1:npwmax,1:nkpe,1:2) )
     igmin(ii) = minval(igv(ii,1:npwmax,1:nkpe,1:2) )
     !
     call MPI_allREDUCE(MPI_IN_PLACE, igmax(ii), 1, &
     &                  MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
     call MPI_allREDUCE(MPI_IN_PLACE, igmin(ii), 1, &
     &                  MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
     !
  end do
  !
  nf(1:3) = igmax(1:3) - igmin(1:3) + 1
  !
  do ii = 1, 3
     !   
     do itr = 1, 100
        call fft_base(nf(ii), lbase)
        if(lbase) exit
        nf(ii) = nf(ii) + 1
     enddo
     !
     if(itr >= 100) then
        write(*,*) 'Stop in fft_mesh.'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        call MPI_FINALIZE(ierr)
        stop 'FFT mesh'
     end if
     !
  enddo
  nftot = product(nf(1:3))
  !
  if(my_rank == 0) then
     write(*,*) "     igmin : ", igmin(1:3)
     write(*,*) "     igmax : ", igmax(1:3)
     write(*,*) "  FFT grid : ", nf(1:3)
  end if
  !
end subroutine fft_mesh
!
! Check fft base is 2, 3, 5, 7, or 11
!
subroutine fft_base(n,lbase)
  !
  use mpi, only : MPI_COMM_WORLD
  !
  integer,intent(in) :: n
  logical,intent(out) :: lbase
  !
  integer :: nmod, quot, itrmax = 100, itr, ierr
  !
  quot = n
  !   
  ! Devide by 2
  !
  do itr = 1, itrmax
     nmod = mod(quot, 2)
     if(nmod /= 0) exit
     quot = quot / 2
  end do
  if(itr >= itrmax) then
     write(*,*) 'Stop in fft_base. factor 2'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     call MPI_FINALIZE(ierr)
     stop 'factor 2'
  end if
  !   
  ! Devide by 3
  !
  do itr = 1, itrmax
     nmod = mod(quot,3)
     if(nmod /= 0) exit
     quot = quot / 3
  end do
  if(itr >= itrmax) then
     write(*,*) 'Stop in fft_base. factor 3'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     call MPI_FINALIZE(ierr)
     stop
  end if
  !
  ! Devide by 5
  !
  do itr = 1, itrmax
     nmod = mod(quot,5)
     if(nmod /= 0) exit
     quot = quot / 5
  end do
  if(itr >= itrmax) then
     write(*,*) 'Stop in fft_base. factor 5'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     call MPI_FINALIZE(ierr)
     stop
  end if
  !   
  ! Devide by 7
  !
  do itr = 1, itrmax
     nmod = mod(quot,7)
     if(nmod /= 0) exit
     quot = quot / 7
  end do
  if(itr >= itrmax) then
     write(*,*) 'Stop in fft_base. factor 7'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     call MPI_FINALIZE(ierr)
     stop
  end if
  !
  ! Devide by 11
  !
  do itr = 1, itrmax
     nmod = mod(quot,11)
     if(nmod /= 0) exit
     quot = quot / 11
  end do
  if(itr >= itrmax) then
     write(*,*) 'Stop in fft_base. factor 11'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     call MPI_FINALIZE(ierr)
     stop
  end if
  !
  if(quot == 1) then
     lbase = .true.
  else
     lbase = .false.
  end if
  !
end subroutine fft_base
!
! FFT wavefunctions
!
subroutine fft_wfc()
  !
  use rpa_el_vals, only : my_rank, nk, nb, npw, wfc, wfc1, wfc2, igv, igmin, nf, nkpe, nftot, &
  &                  FFTW_BACKWARD, FFTW_ESTIMATE
  !
  integer :: ig, ik, ib, igv2(3), ig2, cnt, dsp
  integer(8) :: plan
  complex(8) :: wfcin(nftot), wfcout(nftot)
  !
  call cnt_and_dsp(nk,cnt,dsp)
  !
  allocate(wfc1(nftot, nb, nkpe), wfc2(nftot, nb, nkpe))
  wfc1(1:nftot,1:nb,1:nkpe) = cmplx(0d0, 0d0)
  wfc2(1:nftot,1:nb,1:nkpe) = cmplx(0d0, 0d0)
  !
  call dfftw_plan_dft_3d(plan, nf(1), nf(2), nf(3), wfcin, wfcout, &
  &                      FFTW_BACKWARD, FFTW_ESTIMATE)
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(plan,cnt,dsp,nf,nb,npw,igmin,wfc,wfc1,wfc2,igv,nftot) &
  !$OMP & PRIVATE(wfcin,wfcout,ik,ig,ig2,igv2,ib)
  !
  ! Fourier trans. wfc(G) -> wfc(r)
  ! 
  !$OMP DO
  do ik = 1, cnt
     !
     do ig = 1, npw(ik,1)
        igv2(1:3) = igv(1:3,ig,ik,1) - igmin(1:3) + 1
        ig2 = 1 + (igv2(1) - 1) + nf(1) * (igv2(2) - 1) + nf(1) * nf(2) * (igv2(3) - 1)
        wfc1(ig2, 1:nb, ik) = wfc(ig, 1:nb, ik, 1)
     end do
     !
     do ig = 1, npw(ik,2)
        igv2(1:3) = igv(1:3,ig,ik,2) - igmin(1:3) + 1
        ig2 = 1 + (igv2(1) - 1) + nf(1) * (igv2(2) - 1) + nf(1) * nf(2) * (igv2(3) - 1)
        wfc2(ig2, 1:nb, ik) = wfc(ig, 1:nb, ik, 2)
     end do
     !
     do ib = 1, nb
        !
        wfcin(1:nftot) = wfc1(1:nftot, ib, ik)
        call dfftw_execute_dft(plan, wfcin, wfcout)
        wfc1(1:nftot, ib, ik) = wfcout(1:nftot)
        !
        wfcin(1:nftot) = wfc2(1:nftot, ib, ik)
        call dfftw_execute_dft(plan, wfcin, wfcout)
        wfc2(1:nftot, ib, ik) = wfcout(1:nftot)
        !
     end do
     !
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  !
  call dfftw_destroy_plan(plan)
  !
  deallocate(wfc, igv, npw)
  !
end subroutine fft_wfc
!
! Read dmxc(G)
!
subroutine read_dmxc()
  !
  use mpi, only : MPI_INTEGER, MPI_COMM_WORLD, MPI_DOUBLE_COMPLEX
  use rpa_el_vals, only : my_rank, dmxc, gminxc, gmaxxc
  !
  !
  integer :: fi = 10, ierr
  !
  if(my_rank == 0) then
     !
     open(fi, file = "dmuxc.dat")
     !
     read(fi,*) gminxc(1:3), gmaxxc(1:3)
     write(*,*) gminxc(1:3)
     write(*,*) gmaxxc(1:3)
     allocate(dmxc(gminxc(1):gmaxxc(1), gminxc(2):gmaxxc(2), gminxc(3):gmaxxc(3)))
     !
     read(fi,'(2e25.15)') dmxc(gminxc(1):gmaxxc(1), gminxc(2):gmaxxc(2), gminxc(3):gmaxxc(3))
     !
     close(fi)
     !
  end if
  !
  call MPI_BCAST(gminxc, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(gmaxxc, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if(my_rank /= 0) &
  &  allocate(dmxc(gminxc(1):gmaxxc(1), gminxc(2):gmaxxc(2), gminxc(3):gmaxxc(3)))
  call MPI_BCAST(dmxc, product(gmaxxc(1:3) - gminxc(1:3) + 1),&
  &              MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
  !
end subroutine read_dmxc
!
! Allocate Kel & ikq
!
subroutine alloc_Kel()
  !
  use rpa_el_vals, only : nk, nb, nf, nmf, mf, nftot, gq2, gindx, &
  &                  Kel, pi, my_rank
  !
  integer :: imf
  real(8) :: x0
  !
  allocate(mf(nmf), gq2(nftot), gindx(nftot))
  !
  x0 = cos(pi / dble(2 * (nmf + 2)))
  !
  if(my_rank == 0) write(*,*) "  x0 : ", x0
  !
  do imf = 1, nmf
     !
     mf(imf) = cos(dble(2 * imf + 1) * pi / dble(2 * (nmf + 2)))
     mf(imf) = (x0 + mf(imf)) / (x0 - mf(imf))
     !
  end do
  !
end subroutine alloc_Kel
!
! define shortest diagonal line & define type of tetragonal
!
subroutine tetra_type()
  !
  use mpi, only : MPI_COMM_WORLD
  USE rpa_el_vals, ONLY : my_rank, wlsm, ivvec, ng, ltetra, bvec
  !
  !
  integer :: itype, i1, i2, i3, it, &
  &          divvec(4,4), ivvec0(4), ierr
  real(8) :: l(4), bvec2(3,3), bvec3(3,4)
  !
  do i1 = 1, 3
     bvec2(1:3,i1) = bvec(1:3,i1) / dble(ng(i1))
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
     if(my_rank == 0) write(*,*) "  Linear tetrahedron method is used."
     !
     wlsm(1:4,1:20) = 0.0d0
     wlsm(1,1) = 1.0d0
     wlsm(2,2) = 1.0d0
     wlsm(3,3) = 1.0d0
     wlsm(4,4) = 1.0d0
     !
  else if(ltetra == 2) then
     !
     if(my_rank == 0) write(*,*) "  Optimized tetrahedron method is used."
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
     wlsm(1:4,1:20) = wlsm(1:4,1:20) / 1260d0
     !
  else
     !
     write(*,*) "Stop in tetra_type. ltetra = ", ltetra
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     call MPI_FINALIZE(ierr)
     stop
     !
  end if
  !
end subroutine tetra_type
!
! Screened coulomb interaction
!
subroutine prepare_q(iq)
  !
  use mpi, only : MPI_STATUS_SIZE, MPI_INTEGER, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD
  use rpa_el_vals, only : pi, nk, nb, grid, igmin, ng, bvec, wcut, ngv, gindx, nkpe, &
  &                  nf, nftot, iqv, wfc1, wfc2, gq2, my_rank, petot
  !use rpa_routines, only : cnt_and_dsp_full
  !
  integer,intent(in) :: iq
  !
  integer :: ik, jk, ib, i1, i2, i3, ikv(3), jkv(3), g0(3), ir(3), ifft, ig, dst, src, ierr, &
  &          cnt(0:petot - 1), dsp(0:petot - 1), jkindx(nk), status(MPI_STATUS_SIZE), org, ipe
  real(8) :: gv(3), qv(3), theta, gq20
  complex(8) :: phase(nftot), wfctmp(nftot,nb,nkpe)
  !
  ! |G+q|^2
  !
  qv(1:3) = dble(iqv(1:3,iq)) / dble(2 * ng(1:3))
  qv(1:3) = matmul(bvec(1:3,1:3), qv(1:3))
  !
  ifft = 0
  ngv = 0
  do i3 = 1, nf(3)
     do i2 = 1, nf(2)
        do i1 = 1, nf(1)
           !
           ifft = ifft + 1
           !
           gv(1:3) = dble((/i1, i2, i3/) - 1 + igmin(1:3))
           gv(1:3) = matmul(bvec(1:3,1:3), gv(1:3))
           gv(1:3) = gv(1:3) - qv(1:3)
           !
           gq20 = dot_product(gv(1:3), gv(1:3))
           !
           if(wcut < 1d-10 .or. gq20 < wcut) then
              !
              ngv = ngv + 1
              gq2(ngv) = gq20 / (8d0 * pi)
              gindx(ngv) = ifft
              !
           end if
           !
        end do ! i1 = 1, nf(1)
     end do ! i2 = 1, nf(2)
  end do ! i3 = 1, nf(3)
  !
  if(my_rank == 0) write(*,*) "    # of PWs for W : ", ngv
  !
  ! Prepare wave functions with pahse shift
  !
  call cnt_and_dsp_full(nk, cnt, dsp)
  !
  do ik = 1, nk
     !
     ikv(1:3) = modulo(grid(1:3,ik) + ng(1:3), 2 * ng(1:3)) - ng(1:3)
     !
     jkv(1:3) = grid(1:3,ik) + iqv(1:3,iq)
     jkv(1:3) = modulo(jkv(1:3) + ng(1:3), 2 * ng(1:3)) - ng(1:3)
     !
     g0(1:3) = (jkv(1:3) - ikv(1:3) - iqv(1:3,iq)) / (2 * ng(1:3))
     !
     jkv(1:3) = (jkv(1:3) - 1) / 2
     jkv(1:3) = modulo(jkv(1:3), ng(1:3))
     jk = 1 + jkv(1) + jkv(2) * ng(1) + jkv(3) * ng(1) * ng(2)
     jkindx(ik) = jk
     !
     if(ik <= dsp(my_rank) .or. dsp(my_rank) + cnt(my_rank) < ik) cycle
     !
     ig = 0
     do i3 = 1, nf(3)
        do i2 = 1, nf(2)
           do i1 = 1, nf(1)
              !
              ig = ig + 1
              !
              ir(1:3) = (/i1, i2, i3/) - 1
              ir(1:3) = ir(1:3) * (g0(1:3) + igmin(1:3))
              !             
              theta = sum(dble(ir(1:3)) / dble(nf(1:3)))
              theta = - 2d0 * pi * theta
              !
              phase(ig) = cmplx(cos(theta), sin(theta))
              !
           end do ! i1
        end do ! i2
     end do ! i3
     !
     do ib = 1, nb
        !
        wfc1(1:nftot,ib,ik - dsp(my_rank)) =  wfc1(1:nftot,ib,ik - dsp(my_rank)) &
        &                                  * phase(1:nftot)  / dble(nftot)
        wfc2(1:nftot,ib,ik - dsp(my_rank)) = conjg(wfc2(1:nftot,ib,ik - dsp(my_rank)))
        !
     end do ! ib
     !
  end do ! ik
  !
  !
  !
  wfctmp(1:nftot,1:nb,1:nkpe) = wfc2(1:nftot,1:nb,1:nkpe)
  wfc2(1:nftot,1:nb,1:nkpe) = cmplx(0d0, 0d0)
  !
  dst = modulo(my_rank + 1, petot)
  src  = modulo(my_rank - 1, petot)
  !
  do ipe = 1, petot
     !
     call MPI_SENDRECV_REPLACE(wfctmp, nftot * nb * nkpe, MPI_DOUBLE_COMPLEX, dst, 1, src, &
     &   1, MPI_COMM_WORLD, STATUS, ierr)
     !
     org = modulo(my_rank - ipe, petot)
     !
     do ik = 1, cnt(my_rank)
        !
        if(jkindx(dsp(my_rank) + ik) <= dsp(org) .or. &
        &  dsp(org) + cnt(org) < jkindx(dsp(my_rank) + ik)) cycle
        !
        wfc2(1:nftot,1:nb,ik) = wfctmp(1:nftot,1:nb,jkindx(dsp(my_rank) + ik) - dsp(org))
        !
     end do
     !
  end do
  !
end subroutine prepare_q
!
! Compute screened interaction
!
subroutine make_scrn(iq)
  !
  use mpi, only : MPI_STATUS_SIZE, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD
  use omp_lib, only : OMP_GET_WTIME, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
  use rpa_el_vals, only : nftot, nf, nk, nb, nmf, wfc1, wfc2, wscr, gq2, nkpe, &
  &                  FFTW_FORWARD, FFTW_ESTIMATE, my_rank, petot, ngv, gindx
  !use rpa_routines, only : zgemm, cnt_and_dsp, cnt_and_dsp_full
  !
  integer,intent(in) :: iq
  !
  integer :: ik, ib, jb, imf, cnt, dsp, fstg, lstg, org, dst, src, ierr, ipe, &
  &          kcnt(0:petot - 1), kdsp(0:petot - 1), status(MPI_STATUS_SIZE)
  integer(8) :: plan
  real(8) :: t1, t2
  !
  complex(8) :: wght(0:nmf,nb,nb,nk), rho1(ngv,nb), rhin(nftot), rhout(nftot)
  complex(8) :: one = cmplx(1d0, 0d0)
  complex(8),allocatable :: rho2(:,:)
  !
  call cnt_and_dsp(ngv, cnt, dsp)
  call cnt_and_dsp_full(nk, kcnt, kdsp)
  !
  allocate(wscr(ngv, dsp + 1:dsp + cnt, 0:nmf))
  wscr(1:ngv,dsp + 1:dsp + cnt, 0:nmf) = cmplx(0d0, 0d0)
  !
  ! Calc f * (1 - f') / (e - e' + iw)
  !
  t1 = OMP_GET_WTIME()
  call fermi_fuctor(iq,wght)
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,*) "    fermi_fuctor : ", t2 -t1
  !
  ! Calc. Chi
  !
  call dfftw_plan_dft_3d(plan, nf(1), nf(2), nf(3), rhin(1:nftot), &
  &                      rhout(1:nftot), FFTW_FORWARD, FFTW_ESTIMATE )
  !
  dst = modulo(my_rank + 1, petot)
  src  = modulo(my_rank - 1, petot)
  !
  do ipe = 1, petot
     !
     call MPI_SENDRECV_REPLACE(wfc1, nftot * nb * nkpe, MPI_DOUBLE_COMPLEX, &
     &                         dst, 1, src, 1, MPI_COMM_WORLD, STATUS, ierr  )
     !
     call MPI_SENDRECV_REPLACE(wfc2, nftot * nb * nkpe, MPI_DOUBLE_COMPLEX, &
     &                         dst, 1, src, 1, MPI_COMM_WORLD, STATUS, ierr  )
     !
     org = modulo(my_rank - ipe, petot)
     !
     !$OMP PARALLEL DEFAULT(NONE) &
     !$OMP & SHARED(cnt,dsp,nb,nk,nftot,nmf,one,wfc1,wfc2,plan,wscr,wght,ngv,gindx, &
     !$OMP &        kcnt,kdsp,org) &
     !$OMP & PRIVATE(ik,ib,jb,imf,rhin,rhout,rho1,rho2,fstg,lstg)
     !
     fstg = cnt / OMP_GET_NUM_THREADS() *  OMP_GET_THREAD_NUM() + 1
     lstg = cnt / OMP_GET_NUM_THREADS() * (OMP_GET_THREAD_NUM() + 1)
     !
     fstg = fstg + min(OMP_GET_THREAD_NUM(),     mod(cnt, OMP_GET_NUM_THREADS()))
     lstg = lstg + min(OMP_GET_THREAD_NUM() + 1, mod(cnt, OMP_GET_NUM_THREADS()))
     !
     fstg = fstg + dsp
     lstg = lstg + dsp
     !
     allocate(rho2(fstg:lstg,nb))
     !
     do ik = 1, kcnt(org)
        !
        do ib = 1, nb
           !
           do jb = 1, nb
              !
              rhin(1:nftot) = wfc1(1:nftot,ib,ik) * wfc2(1:nftot,jb,ik)
              call dfftw_execute_dft(plan, rhin(1:nftot), rhout(1:nftot))
              !
              rho1(1:ngv,jb) = rhout(gindx(1:ngv))
              !
           end do ! jb
           !
           do imf = 0, nmf
              !
              do jb = 1, nb
                 !
                 rho2(fstg:lstg,jb) = conjg(wght(imf,jb,ib,kdsp(org) + ik)) * rho1(fstg:lstg,jb)
                 !
              end do ! jb = 1, nb
              !
              call zgemm("N", "C", ngv, lstg - fstg + 1, nb, &
              &  one, rho1(1:ngv,           1:nb), ngv, &
              &       rho2(       fstg:lstg,1:nb), lstg - fstg + 1, &
              &  one, wscr(1:ngv, fstg:lstg, imf), ngv  )
              !
           end do ! imf = 0, nmf
           !
        end do ! ib = 1, nb
        !
     end do ! ik = 1, kcnt(org)
     !
     deallocate(rho2)
     !
     !$OMP END PARALLEL
     !
  end do ! ipe = 1, petot
  !
  call dfftw_destroy_plan(plan)
  !
  wscr(1:ngv,dsp + 1:dsp + cnt, 0:nmf) = &
  &  cmplx(dble(wscr(1:ngv,dsp + 1:dsp + cnt, 0:nmf)), 0d0)
  !
end subroutine make_scrn
!
! Calculation of weight function (f(1-f')) / (e - e')
!
subroutine fermi_fuctor(iq,wght)
  !
  use mpi, only : MPI_IN_PLACE, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD
  use rpa_el_vals, only : nb, nk, nkd, nmf, ng, ngd, iqv, ltetra, ivvec
  !use rpa_routines, only : tetraweight, interpol_weight
  !
  integer,intent(in) :: iq
  complex(8),intent(out) :: wght(0:nmf,nb,nb,nk)
  !
  integer :: cnt, dsp, nkd0, nt, it, ik, i1, i2, i3, ii, iqvd(3), ierr, ikv(3), &
  &          dgrid(3,nkd), indx1(2,20, 6 * nkd), indx2(20, 6 * nkd), indx3(20 * 6 * nkd)
  real(8) :: kv(3)
  complex(8),allocatable :: wghtd(:,:,:,:)
  !
  iqvd(1:3) = iqv(1:3,iq) * ngd(1:3) / (2 * ng(1:3))
  !
  nt = 0
  ik = 0
  do i3 = 1, ngd(3)
     do i2  = 1, ngd(2)
        do i1 = 1, ngd(1)
           !
           ik = ik + 1
           dgrid(1:3,ik) = (/i1, i2, i3/) - 1
           !
           do it = 1, 6
              !
              nt = nt + 1
              !
              do ii = 1, 20
                 !
                 ikv(1:3) = dgrid(1:3,ik) + ivvec(1:3,ii,it)
                 ikv(1:3) = modulo(ikv(1:3), ngd(1:3))
                 !
                 indx1(1,ii,nt) = 1 + ikv(1) + ngd(1) * ikv(2) + ngd(1) * ngd(2) * ikv(3)
                 !
                 ikv(1:3) = dgrid(1:3,ik) + ivvec(1:3,ii,it) + iqvd(1:3)
                 ikv(1:3) = modulo(ikv(1:3), ngd(1:3))
                 !
                 indx1(2,ii,nt) = 1 + ikv(1) + ngd(1) * ikv(2) + ngd(1) * ngd(2) * ikv(3)
                 !
              end do
              !
           end do
           !
        end do
     end do
  end do
  !
  indx2(1:20,1:6 * nkd) = 0
  indx3(1:20 * 6 * nkd) = 0
  !
  call cnt_and_dsp(nkd * 6,cnt,dsp)
  !
  nkd0 = 0
  do it = dsp + 1, dsp + cnt
     !
     do ii = 1, 20
        !
        do ik = 1, nkd0
           !
           if(indx1(1,ii,it) == indx3(ik)) then
              !
              indx2(ii,it) = ik
              goto 10
              !
           end if
           !
        end do
        !
        nkd0 = nkd0 + 1
        indx2(ii,it) = nkd0
        indx3(nkd0) = indx1(1,ii,it)
        !
10      continue
        !
     end do
     !
  end do
  !
  allocate(wghtd(0:nmf,nb,nb,nkd0))
  !
  if(ltetra == 1 .or. ltetra == 2) then
     call tetraweight(nkd0, dgrid,indx1,indx2,iqvd,wghtd)
  else
     write(*,*) "Stop in fermi fuctor. ltetra = ", ltetra
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     call MPI_FINALIZE(ierr)
     stop
  end if
  !
  ! Interpolation of weight
  !
  wght(0:nmf,1:nb,1:nb,1:nk) = 0d0
  do ik = 1, nkd0
     kv(1:3) = dble(dgrid(1:3,indx3(ik))) / dble(ngd(1:3))
     call interpol_weight(nk, (nmf + 1) * nb * nb, ng, kv(1:3), &
     &                     wghtd(0:nmf, 1:nb,1:nb,ik), wght)
  end do
  !
  call MPI_allREDUCE(MPI_IN_PLACE, wght, (nmf + 1) * nb * nb * nk, &
  &                  MPI_DOUBLE_COMPLEX, MPI_SUM,MPI_COMM_WORLD,ierr)
  !
  deallocate(wghtd)
  !
end subroutine fermi_fuctor
!
! Integration weight with tetrahedron method
!
subroutine tetraweight(nkd0,dgrid,indx1,indx2,iqvd,wghtd)
  !
  USE rpa_el_vals, ONLY : nkd, nb, nmf, ngd, ivvec, wlsm, eig, cellV
  !use rpa_routines, only : cnt_and_dsp, sort, tetra2
  !
  integer,intent(in) :: nkd0, dgrid(3,nkd), indx1(2,20,6 * nkd), indx2(20,6 * nkd), iqvd(3)
  complex(8),intent(out) :: wghtd(0:nmf,nb,nb,nkd0)
  !
  integer :: ik, it, ib, jb, imf, ii, cnt, dsp
  real(8) :: thr = 1d-8, V
  real(8) :: e(4), a(4,4), ei(nb,4), ej(nb,4), ei2(4), ej2(nb,4), &
  &        tmp(10,0:nmf,nb,4), tmp2(10,0:nmf,nb,4), &
  &         w0(4,2,0:nmf,nb,4), w1(4,2,0:nmf,nb), w2(4,2,0:nmf,nb,4)
  !
  call cnt_and_dsp(nkd * 6,cnt,dsp)
  !
  wghtd(0:nmf,1:nb,1:nb,1:nkd0) = 0d0
  !
  w0(1:4,1:2,0:nmf,1:nb,1:4) = 0d0
  do ii = 1, 4
     w0(ii,1:2,0:nmf,1:nb,ii) = 1d0
  end do
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(nkd,cnt,dsp,nb,dgrid,ivvec,iqvd,ngd,wlsm,nmf,eig,wghtd,thr, &
  !$OMP &        indx1,indx2,w0) &
  !$OMP & PRIVATE(ik,it,ii,ib,jb,imf, &
  !$OMP &         ei,ej,ei2,ej2,w1,w2,tmp,tmp2,a,e,V)
  !
  do it = dsp + 1, dsp + cnt
     !
     ei(1:nb, 1:4) = 0d0
     ej(1:nb, 1:4) = 0d0
     do ii = 1, 20
        !
        do ib = 1, nb
           !
           ei(ib, 1:4) = ei(ib, 1:4) + wlsm(1:4,ii) * eig(ib, indx1(1,ii,it))
           ej(ib, 1:4) = ej(ib, 1:4) + wlsm(1:4,ii) * eig(ib, indx1(2,ii,it))
           !               
        end do
        !
     end do
     !
     !$OMP DO
     do ib = 1, nb
        !
        w1(1:4,1:2,0:nmf,1:nb) = 0d0
        !
        tmp(1:10, 0:nmf, 1:nb, 1:4) = 0d0
        tmp(   1,     0,    1, 1:4) = ei(                 ib, 1:4)
        tmp(   2,     0, 1:nb, 1:4) = ej(               1:nb, 1:4)
        tmp(3: 6, 0:nmf, 1:nb, 1:4) = w0(1:4, 1, 0:nmf, 1:nb, 1:4)
        tmp(7:10, 0:nmf, 1:nb, 1:4) = w0(1:4, 2, 0:nmf, 1:nb, 1:4)
        !
        call sort(10 * (nmf + 1) * nb, 4, tmp)
        !
        e(1:4) = tmp(1, 0, 1, 1:4)
        !
        do ii = 1, 4
           a(ii,1:4) = ( 0d0 - e(1:4) ) / (e(ii) - e(1:4))
        end do
        !
        if(e(1) <= 0d0 .and. 0d0 < e(2) ) then
           !
           ! A - 1
           !
           V = a(2,1) * a(3,1) * a(4,1)
           !
           if(V > thr) then
              !
              tmp2(1:10,0:nmf,1:nb,1) = tmp(1:10,0:nmf,1:nb,1)
              tmp2(1:10,0:nmf,1:nb,2) = tmp(1:10,0:nmf,1:nb,1) * a(1,2) &
              &                       + tmp(1:10,0:nmf,1:nb,2) * a(2,1) 
              tmp2(1:10,0:nmf,1:nb,3) = tmp(1:10,0:nmf,1:nb,1) * a(1,3) &
              &                       + tmp(1:10,0:nmf,1:nb,3) * a(3,1)
              tmp2(1:10,0:nmf,1:nb,4) = tmp(1:10,0:nmf,1:nb,1) * a(1,4) &
              &                       + tmp(1:10,0:nmf,1:nb,4) * a(4,1)
              !
              ei2(                    1:4) = tmp2(   1,     0,    1, 1:4)
              ej2(              1:nb, 1:4) = tmp2(   2,     0, 1:nb, 1:4)
              w2(1:4, 1, 0:nmf, 1:nb, 1:4) = tmp2(3: 6, 0:nmf, 1:nb, 1:4)
              w2(1:4, 2, 0:nmf, 1:nb, 1:4) = tmp2(7:10, 0:nmf, 1:nb, 1:4)
              !
              call tetra2(ei2,ej2,w2)
              !
              do ii = 1, 4
                 w1(1:4,1:2,0:nmf,1:nb) = w1(1:4,1:2,0:nmf,1:nb) &
                 &                      + w2(1:4,1:2,0:nmf,1:nb,ii) * V
              end do
              !
           end if
           !
        else if(e(2) <= 0d0 .and. 0d0 < e(3)) then
           !
           ! B - 1
           !
           V = a(3,1) * a(4,1) * a(2,4)
           !
           if(V > thr) then
              !
              tmp2(1:10,0:nmf,1:nb,1) = tmp(1:10,0:nmf,1:nb,1)
              tmp2(1:10,0:nmf,1:nb,2) = tmp(1:10,0:nmf,1:nb,1) * a(1,3) &
              &                       + tmp(1:10,0:nmf,1:nb,3) * a(3,1) 
              tmp2(1:10,0:nmf,1:nb,3) = tmp(1:10,0:nmf,1:nb,1) * a(1,4) &
              &                       + tmp(1:10,0:nmf,1:nb,4) * a(4,1) 
              tmp2(1:10,0:nmf,1:nb,4) = tmp(1:10,0:nmf,1:nb,2) * a(2,4) &
              &                       + tmp(1:10,0:nmf,1:nb,4) * a(4,2) 
              !
              ei2(                    1:4) = tmp2(   1,     0,    1, 1:4)
              ej2(              1:nb, 1:4) = tmp2(   2,     0, 1:nb, 1:4)
              w2(1:4, 1, 0:nmf, 1:nb, 1:4) = tmp2(3: 6, 0:nmf, 1:nb, 1:4)
              w2(1:4, 2, 0:nmf, 1:nb, 1:4) = tmp2(7:10, 0:nmf, 1:nb, 1:4)
              ! 
              call tetra2(ei2,ej2,w2)
              !
              do ii = 1, 4
                 w1(1:4,1:2,0:nmf,1:nb) = w1(1:4,1:2,0:nmf,1:nb) &
                 &                      + w2(1:4,1:2,0:nmf,1:nb,ii) * V
              end do
              !
           end if
           !
           ! B - 2
           !
           V = a(3,2) * a(4,2)
           !
           if(V > thr) then
              !
              tmp2(1:10,0:nmf,1:nb,1:2) = tmp(1:10,0:nmf,1:nb,1:2)
              tmp2(1:10,0:nmf,1:nb,  3) = tmp(1:10,0:nmf,1:nb,2) * a(2,3) &
              &                         + tmp(1:10,0:nmf,1:nb,3) * a(3,2) 
              tmp2(1:10,0:nmf,1:nb,  4) = tmp(1:10,0:nmf,1:nb,2) * a(2,4) &
              &                         + tmp(1:10,0:nmf,1:nb,4) * a(4,2) 
              !
              ei2(                    1:4) = tmp2(   1,     0,    1, 1:4)
              ej2(              1:nb, 1:4) = tmp2(   2,     0, 1:nb, 1:4)
              w2(1:4, 1, 0:nmf, 1:nb, 1:4) = tmp2(3: 6, 0:nmf, 1:nb, 1:4)
              w2(1:4, 2, 0:nmf, 1:nb, 1:4) = tmp2(7:10, 0:nmf, 1:nb, 1:4)
              ! 
              call tetra2(ei2,ej2,w2)
              !
              do ii = 1, 4
                 w1(1:4,1:2,0:nmf,1:nb) = w1(1:4,1:2,0:nmf,1:nb) &
                 &                      + w2(1:4,1:2,0:nmf,1:nb,ii) * V
              end do
              !
           end if
           !
           ! B - 3
           !
           V = a(2,3) * a(3,1) * a(4,2)
           !
           if(V > thr) then
              !
              tmp2(1:10,0:nmf,1:nb,1) = tmp(1:10,0:nmf,1:nb,1)
              tmp2(1:10,0:nmf,1:nb,2) = tmp(1:10,0:nmf,1:nb,1) * a(1,3) &
              &                       + tmp(1:10,0:nmf,1:nb,3) * a(3,1) 
              tmp2(1:10,0:nmf,1:nb,3) = tmp(1:10,0:nmf,1:nb,2) * a(2,3) &
              &                       + tmp(1:10,0:nmf,1:nb,3) * a(3,2) 
              tmp2(1:10,0:nmf,1:nb,4) = tmp(1:10,0:nmf,1:nb,2) * a(2,4) &
              &                       + tmp(1:10,0:nmf,1:nb,4) * a(4,2) 
              !
              ei2(                    1:4) = tmp2(   1,     0,    1, 1:4)
              ej2(              1:nb, 1:4) = tmp2(   2,     0, 1:nb, 1:4)
              w2(1:4, 1, 0:nmf, 1:nb, 1:4) = tmp2(3: 6, 0:nmf, 1:nb, 1:4)
              w2(1:4, 2, 0:nmf, 1:nb, 1:4) = tmp2(7:10, 0:nmf, 1:nb, 1:4)
              ! 
              call tetra2(ei2,ej2,w2)
              !
              do ii = 1, 4
                 w1(1:4,1:2,0:nmf,1:nb) = w1(1:4,1:2,0:nmf,1:nb) &
                 &                      + w2(1:4,1:2,0:nmf,1:nb,ii) * V
              end do
              !
           end if
           !
        else if(e(3) <= 0d0 .and. 0d0 < e(4)) then
           !
           ! C - 1
           !
           V = a(4,3)
           !
           if(V > thr) then
              !
              tmp2(1:10,0:nmf,1:nb,1:3) = tmp(1:10,0:nmf,1:nb,1:3)
              tmp2(1:10,0:nmf,1:nb,  4) = tmp(1:10,0:nmf,1:nb,3) * a(3,4) &
              &                         + tmp(1:10,0:nmf,1:nb,4) * a(4,3) 
              !
              ei2(                    1:4) = tmp2(   1,     0,    1, 1:4)
              ej2(              1:nb, 1:4) = tmp2(   2,     0, 1:nb, 1:4)
              w2(1:4, 1, 0:nmf, 1:nb, 1:4) = tmp2(3: 6, 0:nmf, 1:nb, 1:4)
              w2(1:4, 2, 0:nmf, 1:nb, 1:4) = tmp2(7:10, 0:nmf, 1:nb, 1:4)
              ! 
              call tetra2(ei2,ej2,w2)
              !
              do ii = 1, 4
                 w1(1:4,1:2,0:nmf,1:nb) = w1(1:4,1:2,0:nmf,1:nb) &
                 &                      + w2(1:4,1:2,0:nmf,1:nb,ii) * V
              end do
              !
           end if
           !
           ! C - 2
           !
           V = a(3,4) * a(4,2)
           !
           if(V > thr) then
              !
              tmp2(1:10,0:nmf,1:nb,1:2) = tmp(1:10,0:nmf,1:nb,1:2)
              tmp2(1:10,0:nmf,1:nb,  3) = tmp(1:10,0:nmf,1:nb,2) * a(2,4) &
              &                         + tmp(1:10,0:nmf,1:nb,4) * a(4,2) 
              tmp2(1:10,0:nmf,1:nb,  4) = tmp(1:10,0:nmf,1:nb,3) * a(3,4) &
              &                         + tmp(1:10,0:nmf,1:nb,4) * a(4,3) 
              !
              ei2(                    1:4) = tmp2(   1,     0,    1, 1:4)
              ej2(              1:nb, 1:4) = tmp2(   2,     0, 1:nb, 1:4)
              w2(1:4, 1, 0:nmf, 1:nb, 1:4) = tmp2(3: 6, 0:nmf, 1:nb, 1:4)
              w2(1:4, 2, 0:nmf, 1:nb, 1:4) = tmp2(7:10, 0:nmf, 1:nb, 1:4)
              ! 
              call tetra2(ei2,ej2,w2)
              !
              do ii = 1, 4
                 w1(1:4,1:2,0:nmf,1:nb) = w1(1:4,1:2,0:nmf,1:nb) &
                 &                      + w2(1:4,1:2,0:nmf,1:nb,ii) * V
              end do
              !
           end if
           !
           ! C - 3
           !
           V = a(3,4) * a(2,4) * a(4,1)
           !
           if(V > thr) then
              !
              tmp2(1:10,0:nmf,1:nb,1) = tmp(1:10,0:nmf,1:nb,1)
              tmp2(1:10,0:nmf,1:nb,2) = tmp(1:10,0:nmf,1:nb,1) * a(1,4) &
              &                       + tmp(1:10,0:nmf,1:nb,4) * a(4,1) 
              tmp2(1:10,0:nmf,1:nb,3) = tmp(1:10,0:nmf,1:nb,2) * a(2,4) &
              &                       + tmp(1:10,0:nmf,1:nb,4) * a(4,2) 
              tmp2(1:10,0:nmf,1:nb,4) = tmp(1:10,0:nmf,1:nb,3) * a(3,4) &
              &                       + tmp(1:10,0:nmf,1:nb,4) * a(4,3) 
              !
              ei2(                    1:4) = tmp2(   1,     0,    1, 1:4)
              ej2(              1:nb, 1:4) = tmp2(   2,     0, 1:nb, 1:4)
              w2(1:4, 1, 0:nmf, 1:nb, 1:4) = tmp2(3: 6, 0:nmf, 1:nb, 1:4)
              w2(1:4, 2, 0:nmf, 1:nb, 1:4) = tmp2(7:10, 0:nmf, 1:nb, 1:4)
              ! 
              call tetra2(ei2,ej2,w2)
              !
              do ii = 1, 4
                 w1(1:4,1:2,0:nmf,1:nb) = w1(1:4,1:2,0:nmf,1:nb) &
                 &                      + w2(1:4,1:2,0:nmf,1:nb,ii) * V
              end do
              !
           end if
           !
        else if(e(4) <= 0d0 ) then
           !
           ! D - 1
           !
           V = 1d0
           !             
           tmp2(1:10,0:nmf,1:nb,1:4) = tmp(1:10,0:nmf,1:nb,1:4)
           !
           ei2(                    1:4) = tmp2(   1,     0,    1, 1:4)
           ej2(              1:nb, 1:4) = tmp2(   2,     0, 1:nb, 1:4)
           w2(1:4, 1, 0:nmf, 1:nb, 1:4) = tmp2(3: 6, 0:nmf, 1:nb, 1:4)
           w2(1:4, 2, 0:nmf, 1:nb, 1:4) = tmp2(7:10, 0:nmf, 1:nb, 1:4)
           ! 
           call tetra2(ei2,ej2,w2)
           !
           do ii = 1, 4
              w1(1:4,1:2,0:nmf,1:nb) = w1(1:4,1:2,0:nmf,1:nb) &
              &                      + w2(1:4,1:2,0:nmf,1:nb,ii) * V
           end do
           !
        end if
        !
        do ii = 1, 20
           !
           do jb = 1, nb
              !
              do imf = 0, nmf
                 !
                 wghtd(imf,jb,ib,indx2(ii,it)) = wghtd(imf,jb,ib,indx2(ii,it)) &
                 &          + cmplx(sum(wlsm(1:4,ii) * w1(1:4,1,imf,jb)), &
                 &                  sum(wlsm(1:4,ii) * w1(1:4,2,imf,jb))  )
                 !
              end do ! imf = 0, nmf
              !
           end do ! jb = 1, nb
           !               
        end do ! ii = 1, 20
        !
     end do ! ib
     !$OMP END DO NOWAIT
     !
  end do ! it
  !
  !$OMP END PARALLEL
  !
  wghtd(0:nmf,1:nb,1:nb,1:nkd0) = wghtd(0:nmf,1:nb,1:nb,1:nkd0) &
  &                            * 4d0 / (dble(6 * nkd) * cellV)
  !
end subroutine tetraweight
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
     if(a(1,i) > am) then
        atmp(1:n1) = a(1:n1,m)
        a(1:n1,m) = a(1:n1,i)
        a(1:n1,i) = atmp(1:n1)
     end if
  end do
  !
end subroutine sort
!
! Tetrahedra method for 1 - f(ep)
!
subroutine tetra2(ei,ej,w)
  !
  USE rpa_el_vals, ONLY : nb, nmf
  !use rpa_routines, only : sort, lindhard
  !
  real(8),intent(in) :: ei(4), ej(nb,4)
  real(8),intent(inout) :: w(4,2,0:nmf,nb,4)
  !
  integer :: ii, ib
  real(8) :: V, ei2(4), ej2(4), w2(4,2,0:nmf,4), thr = 1d-8
  real(8) :: tmp(10,0:nmf,4), tmp2(10,0:nmf,4), e(4), a(4,4)
  !
  do ib = 1, nb
     !
     tmp(1:10, 0:nmf, 1:4) = 0d0
     tmp(   1,     0, 1:4) = ej(           ib,1:4)
     tmp(   2,     0, 1:4) = ei(              1:4)
     tmp(3: 6, 0:nmf, 1:4) = w(1:4,1,0:nmf,ib,1:4)
     tmp(7:10, 0:nmf, 1:4) = w(1:4,2,0:nmf,ib,1:4)
     !
     call sort(10 * (nmf + 1), 4, tmp)
     !
     e(1:4) = tmp(1, 0, 1:4)
     !
     do ii = 1, 4
        a(ii,1:4) = ( 0d0 - e(1:4) ) / (e(ii) - e(1:4) )
     end do
     !
     w(1:4,1:2,0:nmf,ib,1:4) = 0d0
     !
     if(0d0 <= e(1) ) then
        !
        ! A - 1
        !
        V = 1d0
        !
        tmp2(1:10,0:nmf,1:4) = tmp(1:10,0:nmf,1:4)
        !
        ej2(             1:4) = tmp2(   1,     0, 1:4)
        ei2(             1:4) = tmp2(   2,     0, 1:4)
        w2(1:4, 1, 0:nmf,1:4) = tmp2(3: 6, 0:nmf, 1:4)
        w2(1:4, 2, 0:nmf,1:4) = tmp2(7:10, 0:nmf, 1:4)
        !
        call lindhard(ei2,ej2,w2)
        w(1:4,1:2,0:nmf,ib,1:4) = w(1:4,1:2,0:nmf,ib,1:4) + w2(1:4,1:2,0:nmf,1:4) * V
        !
     else if((e(1) < 0d0 .and. 0d0 <= e(2)) .or. (e(1) <= 0d0 .and. 0d0 < e(2))) then
        !
        ! B - 1
        !
        V = a(1,2)
        !
        if(V > thr) then
           !
           tmp2(1:10,0:nmf,1)   = tmp(1:10,0:nmf,1) * a(1,2) + tmp(1:10,0:nmf,2) * a(2,1)
           tmp2(1:10,0:nmf,2:4) = tmp(1:10,0:nmf,2:4)
           !
           ej2(             1:4) = tmp2(   1,     0, 1:4)
           ei2(             1:4) = tmp2(   2,     0, 1:4)
           w2(1:4, 1, 0:nmf,1:4) = tmp2(3: 6, 0:nmf, 1:4)
           w2(1:4, 2, 0:nmf,1:4) = tmp2(7:10, 0:nmf, 1:4)
           !
           call lindhard(ei2,ej2,w2)
           w(1:4,1:2,0:nmf,ib,1:4) = w(1:4,1:2,0:nmf,ib,1:4) + w2(1:4,1:2,0:nmf,1:4) * V
           !
        end if
        !
        ! B - 2
        !
        V = a(1,3) * a(2,1)
        !
        if(V > thr) then
           !
           tmp2(1:10,0:nmf,1) = tmp(1:10,0:nmf,1) * a(1,2) + tmp(1:10,0:nmf,2) * a(2,1)
           tmp2(1:10,0:nmf,2) = tmp(1:10,0:nmf,1) * a(1,3) + tmp(1:10,0:nmf,3) * a(3,1)
           tmp2(1:10,0:nmf,3:4) = tmp(1:10,0:nmf,3:4)
           !
           ej2(             1:4) = tmp2(   1,     0, 1:4)
           ei2(             1:4) = tmp2(   2,     0, 1:4)
           w2(1:4, 1, 0:nmf,1:4) = tmp2(3: 6, 0:nmf, 1:4)
           w2(1:4, 2, 0:nmf,1:4) = tmp2(7:10, 0:nmf, 1:4)
           !
           call lindhard(ei2,ej2,w2)
           w(1:4,1:2,0:nmf,ib,1:4) = w(1:4,1:2,0:nmf,ib,1:4) + w2(1:4,1:2,0:nmf,1:4) * V
           !
        end if
        !
        ! B - 3
        !
        V = a(1,4) * a(2,1) * a(3,1)
        !
        if(V > thr) then
           !
           tmp2(1:10,0:nmf,1) = tmp(1:10,0:nmf,1) * a(1,2) + tmp(1:10,0:nmf,2) * a(2,1)
           tmp2(1:10,0:nmf,2) = tmp(1:10,0:nmf,1) * a(1,3) + tmp(1:10,0:nmf,3) * a(3,1)
           tmp2(1:10,0:nmf,3) = tmp(1:10,0:nmf,1) * a(1,4) + tmp(1:10,0:nmf,4) * a(4,1)
           tmp2(1:10,0:nmf,4) = tmp(1:10,0:nmf,4)
           !
           ej2(             1:4) = tmp2(   1,     0, 1:4)
           ei2(             1:4) = tmp2(   2,     0, 1:4)
           w2(1:4, 1, 0:nmf,1:4) = tmp2(3: 6, 0:nmf, 1:4)
           w2(1:4, 2, 0:nmf,1:4) = tmp2(7:10, 0:nmf, 1:4)
           !
           call lindhard(ei2,ej2,w2)
           w(1:4,1:2,0:nmf,ib,1:4) = w(1:4,1:2,0:nmf,ib,1:4) + w2(1:4,1:2,0:nmf,1:4) * V
           !
        end if
        !          
     else if((e(2) < 0d0 .and. 0d0 <= e(3)) .or. (e(2) <= 0d0 .and. 0d0 < e(3))) then
        !          
        ! C - 1
        !
        V = a(2,4) * a(1,4) * a(3,1)
        !
        if(V > thr) then
           !
           tmp2(1:10,0:nmf,1) = tmp(1:10,0:nmf,1) * a(1,3) + tmp(1:10,0:nmf,3) * a(3,1)
           tmp2(1:10,0:nmf,2) = tmp(1:10,0:nmf,1) * a(1,4) + tmp(1:10,0:nmf,4) * a(4,1)
           tmp2(1:10,0:nmf,3) = tmp(1:10,0:nmf,2) * a(2,4) + tmp(1:10,0:nmf,4) * a(4,2)
           tmp2(1:10,0:nmf,4) = tmp(1:10,0:nmf,4)
           !
           ej2(             1:4) = tmp2(   1,     0, 1:4)
           ei2(             1:4) = tmp2(   2,     0, 1:4)
           w2(1:4, 1, 0:nmf,1:4) = tmp2(3: 6, 0:nmf, 1:4)
           w2(1:4, 2, 0:nmf,1:4) = tmp2(7:10, 0:nmf, 1:4)
           !
           call lindhard(ei2,ej2,w2)
           w(1:4,1:2,0:nmf,ib,1:4) = w(1:4,1:2,0:nmf,ib,1:4) + w2(1:4,1:2,0:nmf,1:4) * V
           !
        end if
        !
        ! C - 2
        !
        V = a(1,3) * a(2,3)
        !
        if(V > thr) then
           !
           tmp2(1:10,0:nmf,1) = tmp(1:10,0:nmf,1) * a(1,3) + tmp(1:10,0:nmf,3) * a(3,1)
           tmp2(1:10,0:nmf,2) = tmp(1:10,0:nmf,2) * a(2,3) + tmp(1:10,0:nmf,3) * a(3,2)
           tmp2(1:10,0:nmf,3:4) = tmp(1:10,0:nmf,3:4)
           !
           ej2(             1:4) = tmp2(   1,     0, 1:4)
           ei2(             1:4) = tmp2(   2,     0, 1:4)
           w2(1:4, 1, 0:nmf,1:4) = tmp2(3: 6, 0:nmf, 1:4)
           w2(1:4, 2, 0:nmf,1:4) = tmp2(7:10, 0:nmf, 1:4)
           !
           call lindhard(ei2,ej2,w2)
           w(1:4,1:2,0:nmf,ib,1:4) = w(1:4,1:2,0:nmf,ib,1:4) + w2(1:4,1:2,0:nmf,1:4) * V
           !
        end if
        !
        ! C - 3
        ! 
        V = a(1,3) * a(2,4) * a(3,2)
        !
        if(V > thr) then
           !
           tmp2(1:10,0:nmf,1) = tmp(1:10,0:nmf,1) * a(1,3) + tmp(1:10,0:nmf,3) * a(3,1)
           tmp2(1:10,0:nmf,2) = tmp(1:10,0:nmf,2) * a(2,3) + tmp(1:10,0:nmf,3) * a(3,2)
           tmp2(1:10,0:nmf,3) = tmp(1:10,0:nmf,2) * a(2,4) + tmp(1:10,0:nmf,4) * a(4,2)
           tmp2(1:10,0:nmf,4) = tmp(1:10,0:nmf,4)
           !
           ej2(             1:4) = tmp2(   1,     0, 1:4)
           ei2(             1:4) = tmp2(   2,     0, 1:4)
           w2(1:4, 1, 0:nmf,1:4) = tmp2(3: 6, 0:nmf, 1:4)
           w2(1:4, 2, 0:nmf,1:4) = tmp2(7:10, 0:nmf, 1:4)
           !
           call lindhard(ei2,ej2,w2)
           w(1:4,1:2,0:nmf,ib,1:4) = w(1:4,1:2,0:nmf,ib,1:4) + w2(1:4,1:2,0:nmf,1:4) * V
           !
        end if
        !          
     else if((e(3) < 0d0 .and. 0d0 <= e(4)) .or. (e(3) <= 0d0 .and. 0d0 < e(4))) then
        !
        ! D - 1
        !
        V = a(3,4) * a(2,4) * a(1,4) 
        !          
        if(V > thr) then
           !
           tmp2(1:10,0:nmf,1) = tmp(1:10,0:nmf,1) * a(1,4) + tmp(1:10,0:nmf,4) * a(4,1)
           tmp2(1:10,0:nmf,2) = tmp(1:10,0:nmf,2) * a(2,4) + tmp(1:10,0:nmf,4) * a(4,2)
           tmp2(1:10,0:nmf,3) = tmp(1:10,0:nmf,3) * a(3,4) + tmp(1:10,0:nmf,4) * a(4,3)
           tmp2(1:10,0:nmf,4) = tmp(1:10,0:nmf,4)
           !
           ej2(             1:4) = tmp2(   1,     0, 1:4)
           ei2(             1:4) = tmp2(   2,     0, 1:4)
           w2(1:4, 1, 0:nmf,1:4) = tmp2(3: 6, 0:nmf, 1:4)
           w2(1:4, 2, 0:nmf,1:4) = tmp2(7:10, 0:nmf, 1:4)
           !
           call lindhard(ei2,ej2,w2)
           w(1:4,1:2,0:nmf,ib,1:4) = w(1:4,1:2,0:nmf,ib,1:4) + w2(1:4,1:2,0:nmf,1:4) * V
           !
        end if
        !
     end if
     !
  end do
  !
end subroutine tetra2
!
! Tetarahedra method for delta(om - ep + e)
!
subroutine lindhard(ei,ej,w)
  !
  use mpi, only : MPI_COMM_WORLD
  USE rpa_el_vals, ONLY : nmf, mf
  !use rpa_routines, only : sort, &
  !&                        lindhard_1234,   lindhard_1231,   lindhard_1233, &
  !&                        lindhard_1221,   lindhard_1222,   lindhard_1211, &
  !&                        lindhard_1234_0, lindhard_1231_0, lindhard_1233_0, &
  !&                        lindhard_1221_0, lindhard_1222_0, lindhard_1211_0
  !
  !
  real(8),intent(in) :: ei(4), ej(4)
  real(8),intent(inout) :: w(4,2,0:nmf,4)
  !
  integer :: ii, imf, ierr
  real(8) :: tmp(9,0:nmf,4), w2(2,4), de(4), de0(4), lnd(4), thr, thr2
  !
  tmp(1:9, 0:nmf, 1:4) = 0d0
  tmp(  1,     0, 1:4) = ej(1:4) - ei(1:4)
  tmp(2:5, 0:nmf, 1:4) = w(1:4, 1, 0:nmf, 1:4)
  tmp(6:9, 0:nmf, 1:4) = w(1:4, 2, 0:nmf, 1:4)
  call sort(9 * (nmf + 1), 4, tmp)
  de0(          1:4) = tmp(  1,     0, 1:4)
  w(1:4,1,0:nmf,1:4) = tmp(2:5, 0:nmf, 1:4)
  w(1:4,2,0:nmf,1:4) = tmp(6:9, 0:nmf, 1:4)
  !
  thr2 = 1d-12
  !
  do ii = 1, 4
     if(de0(ii) < thr2) then
        if(ii == 3) then
           write(*,*) "Stop in lindhard. Nesting ! "
           write(*,'(9e15.5)') de0(1:4)
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           call MPI_FINALIZE(ierr)
           stop
        end if
        lnd(ii) = 0d0
        de0(ii) = 0d0
     else
        lnd(ii) = log(de0(ii))
     end if
  end do
  !
  de(1:4) = de0(1:4)
  thr = maxval(de(1:4)) * 1d-3
  !
  if(abs(de(4) - de(3)) < thr ) then
     if(abs(de(4) - de(2)) < thr ) then
        if(abs(de(4) - de(1)) < thr ) then
           !
           ! de(4) = de(3) = de(2) = de(1)
           !
           w2(1,4) = 0.25d0 / de(4)
           w2(1,3) = w2(1,4)
           w2(1,2) = w2(1,4)
           w2(1,1) = w2(1,4)
           !
        else
           !
           ! de(4) = de(3) = de(2)
           !
           w2(1,4) = lindhard_1211_0(de(4),de(1),lnd(4),lnd(1))
           w2(1,3) = w2(1,4)
           w2(1,2) = w2(1,4)
           w2(1,1) = lindhard_1222_0(de(1),de(4),lnd(1),lnd(4))
           !
           if(any(w2(1,1:4) < 0d0)) then
              write(*,*) "Stop in lindhard. weighting 4=3=2, imf = ", 0
              write(*,'(100e15.5)') de(1:4)
              write(*,'(100e15.5)') w2(1,1:4)
              call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
              call MPI_FINALIZE(ierr)
              stop
           end if
           !
        end if
     else if(abs(de(2) - de(1)) < thr ) then
        !
        ! de(4) = de(3), de(2) = de(1)
        !
        w2(1,4) = lindhard_1221_0(de(4),de(2), lnd(4),lnd(2))
        w2(1,3) = w2(1,4)
        w2(1,2) = lindhard_1221_0(de(2),de(4), lnd(2),lnd(4))
        w2(1,1) = w2(1,2)
        !
        if(any(w2(1,1:4) < 0d0)) then
           write(*,*) "Stop in lindhard. weighting 4=3 2=1, imf = ", 0
           write(*,'(100e15.5)') de(1:4)
           write(*,'(100e15.5)') w2(1,1:4)
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           call MPI_FINALIZE(ierr)
           stop
        end if
        !
     else
        !
        ! de(4) = de(3)
        !
        w2(1,4) = lindhard_1231_0(de(4),de(1),de(2),lnd(4),lnd(1),lnd(2))
        w2(1,3) = w2(1,4)
        w2(1,2) = lindhard_1233_0(de(2),de(1),de(4),lnd(2),lnd(1),lnd(4))
        w2(1,1) = lindhard_1233_0(de(1),de(2),de(4),lnd(1),lnd(2),lnd(4))
        !
        if(any(w2(1,1:4) < 0d0)) then
           write(*,*) "Stop in lindhard. weighting 4=3, imf = ", 0
           write(*,'(100e15.5)') de(1:4)
           write(*,'(100e15.5)') w2(1,1:4)
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           call MPI_FINALIZE(ierr)
           stop
        end if
        !
     end if
  else if(abs(de(3) - de(2)) < thr) then
     if(abs(de(3) - de(1)) < thr) then
        !
        ! de(3) = de(2) = de(1)
        !
        w2(1,4) = lindhard_1222_0(de(4),de(3), lnd(4),lnd(3))
        w2(1,3) = lindhard_1211_0(de(3),de(4), lnd(3),lnd(4))
        w2(1,2) = w2(1,3)
        w2(1,1) = w2(1,3)
        !
        if(any(w2(1,1:4) < 0d0)) then
           write(*,*) "Stop in lindhard. weighting 3=2=1, imf = ", 0
           write(*,'(100e15.5)') de(1:4)
           write(*,'(100e15.5)') w2(1,1:4)
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           call MPI_FINALIZE(ierr)
           stop
        end if
        !
     else
        !
        ! de(3) = de(2)
        !
        w2(1,4) = lindhard_1233_0(de(4),de(1),de(3),lnd(4),lnd(1),lnd(3))
        w2(1,3) = lindhard_1231_0(de(3),de(1),de(4),lnd(3),lnd(1),lnd(4))
        w2(1,2) = w2(1,3)
        w2(1,1) = lindhard_1233_0(de(1),de(4),de(3),lnd(1),lnd(4),lnd(3))
        !
        if(any(w2(1,1:4) < 0d0)) then
           write(*,*) "Stop in lindhard. weighting 3=2, imf = ", 0
           write(*,'(100e15.5)') de(1:4)
           write(*,'(100e15.5)') w2(1,1:4)
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           call MPI_FINALIZE(ierr)
           stop
        end if
        !
     end if
  else if(abs(de(2) - de(1)) < thr) then
     !
     ! de(2) = de(1)
     !
     w2(1,4) = lindhard_1233_0(de(4),de(3),de(2),lnd(4),lnd(3),lnd(2))
     w2(1,3) = lindhard_1233_0(de(3),de(4),de(2),lnd(3),lnd(4),lnd(2))
     w2(1,2) = lindhard_1231_0(de(2),de(3),de(4),lnd(2),lnd(3),lnd(4))
     w2(1,1) = w2(1,2)
     !
     if(any(w2(1,1:4) < 0d0)) then
        write(*,*) "Stop in lindhard. weighting 2=1, imf = ", 0
        write(*,'(100e15.5)') de(1:4)
        write(*,'(100e15.5)') w2(1,1:4)
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        call MPI_FINALIZE(ierr)
        stop
     end if
     !
  else
     !
     ! Different each other.
     !
     w2(1,4) = lindhard_1234_0(de(4),de(1),de(2),de(3),lnd(4),lnd(1),lnd(2),lnd(3))
     w2(1,3) = lindhard_1234_0(de(3),de(1),de(2),de(4),lnd(3),lnd(1),lnd(2),lnd(4))
     w2(1,2) = lindhard_1234_0(de(2),de(1),de(3),de(4),lnd(2),lnd(1),lnd(3),lnd(4))
     w2(1,1) = lindhard_1234_0(de(1),de(2),de(3),de(4),lnd(1),lnd(2),lnd(3),lnd(4))
     !      
     if(any(w2(1,1:4) < 0d0)) then
        write(*,*) "Stop in lindhard. weighting each other, imf = ", 0
        write(*,'(100e15.5)') de(1:4)
        write(*,'(100e15.5)') w2(1,1:4)
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        call MPI_FINALIZE(ierr)
        stop
     end if
     !
  end if
  !
  do ii = 1, 4
     w(1:4,1,0,ii) = w2(1,ii) * w(1:4,1,0,ii)
     w(1:4,2,0,ii) = 0d0
  end do ! ii
  !
  ! w /= 0 part
  !
  do imf = 1, nmf
     !
     de(1:4) = de0(1:4) / mf(imf)
     !thr = maxval(de(1:4)) * 1d-3
     thr = max(1d-3,  maxval(de(1:4)) * 1d-2)
     !
     if(abs(de(4) - de(3)) < thr ) then
        if(abs(de(4) - de(2)) < thr ) then
           if(abs(de(4) - de(1)) < thr ) then
              !
              ! de(4) = de(3) = de(2) = de(1)
              !
              w2(1,4) = 0.25d0 * de(4) / ((1d0 + de(4)**2))
              w2(2,4) = 0.25d0         / ((1d0 + de(4)**2))
              w2(1:2,3) = w2(1:2,4)
              w2(1:2,2) = w2(1:2,4)
              w2(1:2,1) = w2(1:2,4)
              !
           else
              !
              ! de(4) = de(3) = de(2)
              !
              w2(1:2,4) = lindhard_1211(de(4),de(1))
              w2(1:2,3) = w2(1:2,4)
              w2(1:2,2) = w2(1:2,4)
              w2(1:2,1) = lindhard_1222(de(1),de(4))
              !
              if(any(w2(1:2,1:4) < 0d0)) then
                 write(*,*) "Stop in lindhard. weighting 4=3=2. imf = ", imf
                 write(*,'(100e15.5)') de(1:4)
                 write(*,'(2e15.5)') w2(1:2,1:4)
                 call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
                 call MPI_FINALIZE(ierr)
                 stop
              end if
              !
           end if
        else if(abs(de(2) - de(1)) < thr ) then
           !
           ! de(4) = de(3), de(2) = de(1)
           !
           w2(1:2,4) = lindhard_1221(de(4),de(2))
           w2(1:2,3) = w2(1:2,4)
           w2(1:2,2) = lindhard_1221(de(2),de(4))
           w2(1:2,1) = w2(1:2,2)
           !
           if(any(w2(1:2,1:4) < 0d0)) then
              write(*,*) "Stop in lindhard. weighting 4=3, 2=1. imf = ", imf
              write(*,'(100e15.5)') de(1:4)
              write(*,'(2e15.5)') w2(1:2,1:4)
              call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
              call MPI_FINALIZE(ierr)
              stop
           end if
           !
        else
           !
           ! de(4) = de(3)
           !
           w2(1:2,4) = lindhard_1231(de(4),de(1),de(2))
           w2(1:2,3) = w2(1:2,4)
           w2(1:2,2) = lindhard_1233(de(2),de(1),de(4))
           w2(1:2,1) = lindhard_1233(de(1),de(2),de(4))
           !
           if(any(w2(1:2,1:4) < 0d0)) then
              write(*,*) "Stop in lindhard. weighting 4=3. imf = ", imf
              write(*,'(100e15.5)') de(1:4)
              write(*,'(2e15.5)') w2(1:2,1:4)
              call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
              call MPI_FINALIZE(ierr)
              stop
           end if
           !
        end if
     else if(abs(de(3) - de(2)) < thr) then
        if(abs(de(3) - de(1)) < thr) then
           !
           ! de(3) = de(2) = de(1)
           !
           w2(1:2,4) = lindhard_1222(de(4),de(3))
           w2(1:2,3) = lindhard_1211(de(3),de(4))
           w2(1:2,2) = w2(1:2,3)
           w2(1:2,1) = w2(1:2,3)
           !
           if(any(w2(1:2,1:4) < 0d0)) then
              write(*,*) "Stop in lindhard. weighting 3=2=1. imf = ", imf
              write(*,'(100e15.5)') de(1:4)
              write(*,'(2e15.5)') w2(1:2,1:4)
              call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
              call MPI_FINALIZE(ierr)
              stop
           end if
           !
        else
           !
           ! de(3) = de(2)
           !
           w2(1:2,4) = lindhard_1233(de(4),de(1),de(3))
           w2(1:2,3) = lindhard_1231(de(3),de(1),de(4))
           w2(1:2,2) = w2(1:2,3)
           w2(1:2,1) = lindhard_1233(de(1),de(4),de(3))
           !
           if(any(w2(1:2,1:4) < 0d0)) then
              write(*,*) "Stop in lindhard. weighting 3=2. imf = ", imf
              write(*,'(100e15.5)') de(1:4)
              write(*,'(2e15.5)') w2(1:2,1:4)
              call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
              call MPI_FINALIZE(ierr)
              stop
           end if
           !
        end if
     else if(abs(de(2) - de(1)) < thr) then
        !
        ! de(2) = de(1)
        !
        w2(1:2,4) = lindhard_1233(de(4),de(3),de(2))
        w2(1:2,3) = lindhard_1233(de(3),de(4),de(2))
        w2(1:2,2) = lindhard_1231(de(2),de(3),de(4))
        w2(1:2,1) = w2(1:2,2)
        !
        if(any(w2(1:2,1:4) < 0d0)) then
           write(*,*) "Stop in lindhard. weighting 2=1. imf = ", imf
           write(*,'(100e15.5)') de(1:4)
           write(*,'(2e15.5)') w2(1:2,1:4)
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           call MPI_FINALIZE(ierr)
           stop
        end if
        !
     else
        !
        ! Different each other.
        !
        w2(1:2,4) = lindhard_1234(de(4),de(1),de(2),de(3))
        w2(1:2,3) = lindhard_1234(de(3),de(1),de(2),de(4))
        w2(1:2,2) = lindhard_1234(de(2),de(1),de(3),de(4))
        w2(1:2,1) = lindhard_1234(de(1),de(2),de(3),de(4))
        !      
        if(any(w2(1:2,1:4) < 0d0)) then
           write(*,*) "Stop in lindhard. weighting each other. imf = ", imf
           write(*,'(100e15.5)') de(1:4)
           write(*,'(2e15.5)') w2(1:2,1:4)
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           call MPI_FINALIZE(ierr)
           stop
        end if
        !
     end if
     !
     do ii = 1, 4
        w(1:4,1,imf,ii) = w2(1,ii) * w(1:4,1,imf,ii) /    mf(imf)
        w(1:4,2,imf,ii) = w2(2,ii) * w(1:4,2,imf,ii) / (- mf(imf))
     end do ! ii
     !
  end do ! imf
  !
end subroutine lindhard
!
! Results of Integration (1-x-y-z)/(g0+(g1-g0)x+(g2-g0)y+(g3-g0))
!  for 0<x<1, 0<y<1-x, 0<z<1-x-y
!
! 1, Different each other
!
function lindhard_1234_0(g1,g2,g3,g4,lng1,lng2,lng3,lng4) result(w)
  !
  !
  real(8),intent(in) :: g1,g2,g3,g4,lng1,lng2,lng3,lng4
  real(8) :: w
  !
  real(8) :: w2, w3, w4
  !
  w2 = ((lng2 - lng1)/(g2 - g1)*g2 - 1d0)*g2/(g2 - g1)
  w3 = ((lng3 - lng1)/(g3 - g1)*g3 - 1d0)*g3/(g3 - g1)
  w4 = ((lng4 - lng1)/(g4 - g1)*g4 - 1d0)*g4/(g4 - g1)
  w2 = ((w2 - w3)*g2)/(g2 - g3)
  w4 = ((w4 - w3)*g4)/(g4 - g3)
  w = (w4 - w2)/(g4 - g2)
  !
end function lindhard_1234_0
!
! 1, Different each other
!
function lindhard_1234(g1,g2,g3,g4) result(w)
  !
  !
  real(8),intent(in) :: g1, g2, g3, g4
  real(8) :: w(2)
  !
  real(8) :: w2, w3, w4
  !
  ! Real
  !
  w2 = 2d0*(3d0*g2**2 - 1d0)*(atan(g2) - atan(g1)) + (g2**2 - &
  &      3d0)*g2*log((1d0 + g2**2)/( 1d0 + g1**2))
  w2 = -2d0*(g2**2 - 1d0) + w2/(g2 - g1 )
  w2 = w2/(g2 - g1 )
  w3 = 2d0*(3d0*g3**2 - 1d0)*(atan(g3) - atan(g1)) + (g3**2 -  &
  &      3d0)*g3*log((1d0 + g3**2)/( 1d0 + g1**2))
  w3 = -2d0*(g3**2 - 1d0) + w3/(g3 - g1 )
  w3 = w3/(g3 - g1 )
  w4 = 2d0*(3d0*g4**2 - 1d0)*(atan(g4) - atan(g1)) + (g4**2 -  &
  &      3d0)*g4*log((1d0 + g4**2)/( 1d0 + g1**2))
  w4 = -2d0*(g4**2 - 1d0) + w4/(g4 - g1 )
  w4 = w4/(g4 - g1 )
  w2 = (w2 - w3)/(g2 - g3)
  w4 = (w4 - w3)/(g4 - g3)
  w(1) = (w4 - w2)/(2d0*(g4 - g2))
  !
  ! Imaginal
  !
  w2 = 2d0*(3d0 - g2**2)* &
  &    g2*(atan(g2) - atan(g1)) + (3d0*g2**2 - 1d0)* &
  &    log((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 4d0*g2 - w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0*(3d0 - g3**2)* &
  &    g3*(atan(g3) - atan(g1)) + (3d0*g3**2 - 1d0)* &
  &    log((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 4d0*g3 - w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w4 = 2d0*(3d0 - g4**2)* &
  &    g4*(atan(g4) - atan(g1)) + (3d0*g4**2 - 1d0)* &
  &    log((1d0 + g4**2)/(1d0 + g1**2))
  w4 = 4d0*g4 - w4/(g4 - g1)
  w4 = w4/(g4 - g1)
  w2 = (w2 - w3)/(g2 - g3)
  w4 = (w4 - w3)/(g4 - g3)
  w(2) = (w4 - w2)/(2d0*(g4 - g2))
  !
end function lindhard_1234
!
! 2, g4 = g1
!
function lindhard_1231_0(g1,g2,g3,lng1,lng2,lng3) result(w)
  !
  !
  real(8),intent(in) :: g1,g2,g3,lng1,lng2,lng3
  real(8) :: w
  !
  real(8) :: w2, w3
  !
  w2 = ((lng2 - lng1)/(g2 - g1)*g2 - 1d0)*g2**2/(g2 - g1) - g1/( &
  &   2d0)
  w2 = w2/(g2 - g1)
  w3 = ((lng3 - lng1)/(g3 - g1)*g3 - 1d0)*g3**2/(g3 - g1) - g1/( &
  &   2d0)
  w3 = w3/(g3 - g1)
  w = (w3 - w2)/(g3 - g2)
  !
end function lindhard_1231_0
!
! 2, g4 = g1
!
function lindhard_1231(g1,g2,g3) result(w)
  !
  !
  real(8),intent(in) :: g1, g2, g3
  real(8) :: w(2)
  !
  real(8) :: w2, w3
  !
  ! Real
  !
  w2 = 2d0*(-1d0 + 3d0*g2**2)*(atan(g2) - atan(g1)) +  &
  &   g2*(-3d0 + g2**2)*log((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 2d0*(1d0 - g2**2) + w2/(g2 - g1)
  w2 = -g1 + w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0*(-1d0 + 3d0*g3**2)*(atan(g3) - atan(g1)) +  &
  &   g3*(-3d0 + g3**2)*log((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 2d0*(1 - g3**2) + w3/(g3 - g1)
  w3 = -g1 + w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w(1) = (w3 - w2)/(2d0*(g3 - g2))
  !
  ! Imaginal
  !
  w2 = 2d0* &
  &    g2*(3d0 - g2**2)*(atan(g2) - atan(g1)) + (-1d0 + 3d0*g2**2)* &
  &    log((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 4d0*g2 - w2/(g2 - g1)
  w2 = 1 + w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0* &
  &    g3*(3d0 - g3**2)*(atan(g3) - atan(g1)) + (-1d0 + 3d0*g3**2)* &
  &    log((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 4d0*g3 - w3/(g3 - g1)
  w3 = 1 + w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w(2) = (w3 - w2)/(2d0*(g3 - g2))
  !
end function lindhard_1231
!
! 3, g4 = g3
!
function lindhard_1233_0(g1,g2,g3,lng1,lng2,lng3) result(w)
  !
  !
  real(8),intent(in) :: g1,g2,g3,lng1,lng2,lng3
  real(8) :: w
  !
  real(8) :: w2, w3
  !
  w2 = (lng2 - lng1)/(g2 - g1)*g2 - 1d0
  w2 = (g2*w2)/(g2 - g1)
  w3 = (lng3 - lng1)/(g3 - g1)*g3 - 1d0
  w3 = (g3*w3)/(g3 - g1)
  w2 = (w3 - w2)/(g3 - g2)
  w3 = (lng3 - lng1)/(g3 - g1)*g3 - 1d0
  w3 = 1d0 - (2d0*w3*g1)/(g3 - g1)
  w3 = w3/(g3 - g1)
  w = (g3*w3 - g2*w2)/(g3 - g2)
  !
end function lindhard_1233_0
!
! 3, g4 = g3
!
function lindhard_1233(g1, g2, g3) result(w)
  !
  !
  real(8),intent(in) :: g1, g2, g3
  real(8) :: w(2)
  !
  real(8) :: w2, w3
  !
  ! Real
  !
  w2 = 2d0*(1d0 - 3d0*g2**2)*(atan(g2) - atan(g1)) +  &
  &   g2*(3d0 - g2**2)*log((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 2d0*(1 - g2**2) - w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0*(1d0 - 3d0*g3**2)*(atan(g3) - atan(g1)) +  &
  &   g3*(3d0 - g3**2)*log((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 2d0*(1 - g3**2) - w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w2 = (w3 - w2)/(g3 - g2)
  w3 = 4d0*(1d0 - 3d0*g1*g3)*(atan(g3) - atan(g1)) + (3d0*g1 +  &
  &      3d0*g3 - 3d0*g1*g3**2 + g3**3) * log((1d0 + g3**2)/( &
  &     1d0 + g1**2))
  w3 = -4d0*(1d0 - g1**2) + w3/(g3 - g1)
  w3 = 4d0*g1 + w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w(1) = (w3 - w2)/(2d0*(g3 - g2))
  !
  ! Imaginal
  !
  w2 = 2d0* &
  &    g2*(3d0 - g2**2)*(atan(g2) - atan(g1)) + (-1d0 + 3d0*g2**2)* &
  &    log((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 4d0*g2 - w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0* &
  &    g3*(3d0 - g3**2)*(atan(g3) - atan(g1)) + (-1d0 + 3d0*g3**2)* &
  &    log((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 4d0*g3 - w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w2 = (w3 - w2)/(g3 - g2)
  w3 = (3d0*g1 - 3d0*g1*g3**2 + 3d0*g3 + g3**3)*(atan(g3) -  &
  &      atan(g1)) + (3d0*g1*g3 - 1d0)* &
  &    log((1d0 + g3**2)/(1d0 + g1**2))
  w3 = w3/(g3 - g1) - 4d0*g1
  w3 = w3/(g3 - g1) - 2d0
  w3 = (2d0*w3)/(g3 - g1)
  w(2) = (w3 - w2)/(2d0*(g3 - g2))
  !
end function lindhard_1233
!
! 4, g4 = g1 and g3 = g2
!
function lindhard_1221_0(g1,g2,lng1,lng2) result(w)
  !
  !
  real(8),intent(in) :: g1, g2, lng1, lng2
  real(8) :: w
  !
  w = 1d0 - (lng2 - lng1)/(g2 - g1)*g1
  w = -1d0 + (2d0*g2*w)/(g2 - g1)
  w = -1d0 + (3d0*g2*w)/(g2 - g1)
  w = w/(2d0*(g2 - g1))
  !
end function lindhard_1221_0
!
! 4, g4 = g1 and g3 = g2
!
function lindhard_1221(g1,g2) result(w)
  !
  !
  real(8),intent(in) :: g1, g2
  real(8) :: w(2)
  !
  ! Real
  !
  w(1) = -2d0*(-1d0 + 2d0*g1*g2 + g2**2)*(atan(g2) -  &
  &      atan(g1)) + (g1 + 2d0*g2 - g1*g2**2)* &
  &    log((1d0 + g2**2)/(1d0 + g1**2))
  w(1) = 2d0*(-1d0 + g1**2) + w(1)/(g2 - g1)
  w(1) = 3d0*g1 + w(1)/(g2 - g1)
  w(1) = 2d0 + (3d0*w(1))/(g2 - g1)
  w(1) = w(1)/(2d0*(g2 - g1))
  !
  ! Imaginal
  !
  w(2) = 2d0*(g1 + 2d0*g2 - g1*g2**2)*(atan(g2) -  &
  &      atan(g1)) + (-1d0 + 2d0*g1*g2 + g2**2)* &
  &    log((1 + g2**2)/(1 + g1**2))
  w(2) = -4d0*g1 + w(2)/(g2 - g1)
  w(2) = -3d0 + w(2)/(g2 - g1)
  w(2) = (3d0*w(2))/(2d0*(g2 - g1)**2)
  !
end function lindhard_1221
!
! 5, g4 = g3 = g2
!
function lindhard_1222_0(g1,g2,lng1,lng2) result(w)
  !
  !
  real(8),intent(in) :: g1, g2, lng1, lng2
  real(8) :: w
  !
  w = (lng2 - lng1)/(g2 - g1)*g2 - 1d0
  w = (2d0*g1*w)/(g2 - g1) - 1d0
  w = (3d0*g1*w)/(g2 - g1) + 1d0
  w = w/(2d0*(g2 - g1))
  !
end function lindhard_1222_0
!
! 5, g4 = g3 = g2
!
function lindhard_1222(g1,g2) result(w)
  !
  !
  real(8),intent(in) :: g1, g2
  real(8) :: w(2)
  !
  ! Real
  !
  w(1) = 2d0*(-1d0 + g1**2 + 2d0*g1*g2)*(atan(g2) -  &
  &      atan(g1)) + (-2d0*g1 - g2 + g1**2*g2) * log((1d0 + g2**2)/( &
  &     1d0 + g1**2))
  w(1) = 2d0*(1d0 - g1**2) + w(1)/(g2 - g1)
  w(1) = g1 - w(1)/(g2 - g1)
  w(1) = 1d0 - (3d0*w(1))/(g2 - g1)
  w(1) = w(1)/(2d0*(g2 - g1))
  !
  ! Imaginal
  !
  w(2) = 2d0*(-2d0*g1 - g2 + g1**2*g2)*(atan(g2) - atan(g1)) + (1d0 - &
  &       g1**2 - 2d0*g1*g2) * log((1d0 + g2**2)/(1d0 + g1**2))
  w(2) = 4d0*g1 + w(2)/(g2 - g1)
  w(2) = 1d0 + w(2)/(g2 - g1)
  w(2) = (3d0*w(2))/(2d0*(g2 - g1)**2)
  !
end function lindhard_1222
!
! 6, g4 = g3 = g1
!
function lindhard_1211_0(g1,g2,lng1,lng2) result(w)
  !
  !
  real(8),intent(in) :: g1,g2,lng1,lng2
  real(8) :: w
  !
  w = -1d0 + (lng2 - lng1)/(g2 - g1)*g2
  w = -1d0 + (2d0*g2*w)/(g2 - g1)
  w = -1d0 + (3d0*g2*w)/(2d0*(g2 - g1))
  w = w/(3d0*(g2 - g1))
  !
end function lindhard_1211_0
!
! 6, g4 = g3 = g1
!
function lindhard_1211(g1,g2) result(w)
  !
  !
  real(8),intent(in) :: g1, g2
  real(8) :: w(2)
  !
  ! Real
  !
  w(1) = 2d0*(3d0*g2**2 - 1d0)*(atan(g2) - atan(g1)) +  &
  &   g2*(g2**2 - 3d0)*log((1d0 + g2**2)/(1d0 + g1**2))
  w(1) = 2d0*(1d0 - g1**2) + w(1)/(g2 - g1)
  w(1) = -5d0*g1 + w(1)/(g2 - g1)
  w(1) = -11d0 + (3d0*w(1))/(g2 - g1)
  w(1) = w(1)/(6d0*(g2 - g1))
  !
  ! Imaginal
  !
  w(2) = 2d0*g2*(-3d0 + g2**2)*(atan(g2) - atan(g1)) + (1d0 -  &
  &      3d0*g2**2)*log((1d0 + g2**2)/(1d0 + g1**2))
  w(2) = 4d0*g2 + w(2)/(g2 - g1)
  w(2) = 1d0 + w(2)/(g2 - g1)
  w(2) = w(2)/(2d0*(g2 - g1)**2)
  !
end function lindhard_1211
!
! first or third order interpolation of weights
!
subroutine interpol_weight(nk,nb,ng,ko,wi,wo)
  !
  use mpi, only : MPI_COMM_WORLD
  use rpa_el_vals, only : ivvec, ltetra
  !
  integer,intent(in)  :: nk, nb, ng(3)
  real(8),intent(in)  :: ko(3)
  complex(8),intent(in) ::wi(nb)
  complex(8),intent(inout) :: wo(nb,nk)
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
  write(*,*) "Stop in interpol weight."
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
     write(*,*) "Stop in interpol_weight. ltetra = ", ltetra
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     call MPI_FINALIZE(ierr)
     stop
     ! 
  end if
  !
end subroutine interpol_weight
!
! Apply XC term
!
subroutine apply_xc()
  !
  use rpa_el_vals, only : nmf, nftot, igmin, nf, wscr, gq2, &
  &                  dmxc, gminxc, gmaxxc, my_rank, laddxc, ngv, gindx
  !use rpa_routines, only : cnt_and_dsp
  !
  integer :: i1, i2, i3, ig, jg, imf, idg(3), igv0(3,nftot), cnt, dsp
  complex(8) :: vec(ngv), one = cmplx(1d0, 0d0), zero = cmplx(0d0, 0d0)
  complex(8),allocatable :: wscr0(:,:,:)
  !
  call cnt_and_dsp(ngv, cnt, dsp)
  allocate(wscr0(dsp + 1:dsp + cnt,0:nmf,ngv))
  !
  if(laddxc == 1) then
     !
     ig = 0
     do i3 = 1, nf(3)
        do i2 = 1, nf(2)
           do i1 = 1, nf(1)
              ig = ig + 1
              igv0(1:3,ig) = igmin(1:3) + (/i1, i2, i3/) - 1
           end do ! i1
        end do ! i2
     end do ! i3
     !
     !$OMP PARALLEL DEFAULT(NONE) &
     !$OMP & SHARED(ngv,igv0,gminxc,gmaxxc,gq2,dmxc,cnt, &
     !$OMP &        dsp,one,zero,wscr,wscr0,nmf,my_rank,gindx) &
     !$OMP & PRIVATE(ig,jg,idg,vec,imf)
     !
     !$OMP DO
     do ig = 1, ngv
        !
        do jg = 1, ngv
           !
           idg(1:3) = igv0(1:3,gindx(ig)) - igv0(1:3,gindx(jg))
           !
           if(all(gminxc(1:3) <= idg(1:3)) .and. all(idg(1:3) <= gmaxxc(1:3))) then
              vec(jg) = cmplx(gq2(ig), 0d0) * dmxc(idg(1), idg(2), idg(3))
           else
              vec(jg) = cmplx(0d0, 0d0)
           end if
           !
!if(my_rank == 0) write(90,'(2i7,2e15.5)') ig, jg, vec(jg)
           !
        end do
        vec(ig) = vec(ig) + cmplx(1d0, 0d0)
        !
        do imf = 0, nmf
           call zgemv("T", ngv, cnt, one, wscr(1:ngv, dsp + 1:dsp + cnt, imf), &
           &          ngv, vec, 1, zero, wscr0(dsp+1:dsp+cnt,imf,ig), 1)
        end do
        !
     end do ! ig
     !$OMP END DO
     !$OMP END PARALLEL
     !
     do ig = 1, ngv
        wscr(ig, dsp + 1:dsp + cnt, 0:nmf) = wscr0(dsp + 1:dsp + cnt, 0:nmf, ig)
     end do
     !
  end if
  !
  do ig = dsp + 1, dsp + cnt
     wscr(ig, ig, 0:nmf) = wscr(ig, ig, 0:nmf) + cmplx(gq2(ig), 0d0)
  end do
  !
  deallocate(wscr0)
  !
end subroutine apply_xc
!
! Invertion of matrices
!
subroutine invert()
  !
  use mpi, only : MPI_IN_PLACE, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, MPI_COMM_WORLD, &
  &               MPI_DOUBLE_COMPLEX, MPI_SUM
  use rpa_el_vals, only : wscr, ngv, nmf, pi, cellV
  !use rpa_routines, only : cnt_and_dsp, zscal, zcopy, zaxpy
  !
  integer :: cnt, dsp, imf, ig, jg, ierr
  real(8) :: maxpiv(2)
  complex(8) :: key1(ngv), key2(ngv), piv
  complex(8),allocatable :: wscr1(:,:), wscr2(:,:)
  !
  call cnt_and_dsp(ngv, cnt, dsp)
  allocate(wscr1(ngv,dsp + 1:dsp + cnt), wscr2(ngv,dsp + 1:dsp + cnt))
  !
  do imf = 0, nmf
     !
     wscr1(1:ngv, dsp + 1:dsp + cnt) = &
     &    wscr(1:ngv, dsp + 1:dsp + cnt,imf)
     !
     wscr2(1:ngv,dsp + 1:dsp + cnt) = cmplx(0d0, 0d0)
     do ig = dsp + 1, dsp + cnt
        wscr2(ig,ig) = cmplx(1d0, 0d0)
     end do
     !
     do ig = 1, ngv
        !
        ! Percial pivotting
        !
        jg = max(ig, dsp + 1)
        if(jg > dsp + cnt) then
           maxpiv(1) = - 1d10
           maxpiv(2) =   1d0
        else
           maxpiv(1) = maxval(dble(conjg(wscr1(ig, jg:dsp + cnt)) &
           &                           * wscr1(ig, jg:dsp + cnt) ))
           maxpiv(2) = maxloc(dble(conjg(wscr1(ig, jg:dsp + cnt)) &
           &                           * wscr1(ig, jg:dsp + cnt) ), 1) &
           &         + jg - 1
        end if
        !
        call MPI_allREDUCE(MPI_IN_PLACE, maxpiv, 1, &
        &                  MPI_2DOUBLE_PRECISION, MPI_MAXLOC,MPI_COMM_WORLD,ierr)
        !
        jg = nint(maxpiv(2))
        !
        ! Exchange wscr1
        !
        key1(1:ngv) = cmplx(0d0, 0d0)
        key2(1:ngv) = cmplx(0d0, 0d0)
        !
        if(dsp + 1 <= ig  .and. ig  <= dsp + cnt) key1(1:ngv) = wscr1(1:ngv,ig)
        if(dsp + 1 <= jg  .and. jg  <= dsp + cnt) key2(1:ngv) = wscr1(1:ngv,jg)
        !
        call MPI_allREDUCE(MPI_IN_PLACE, key1, ngv, &
        &                  MPI_DOUBLE_COMPLEX, MPI_SUM,MPI_COMM_WORLD,ierr)
        !
        call MPI_allREDUCE(MPI_IN_PLACE, key2, ngv, &
        &                  MPI_DOUBLE_COMPLEX, MPI_SUM,MPI_COMM_WORLD,ierr)
        !
        if(dsp + 1 <= ig  .and. ig  <= dsp + cnt) wscr1(1:ngv,ig) = key2(1:ngv)
        if(dsp + 1 <= jg  .and. jg  <= dsp + cnt) wscr1(1:ngv,jg) = key1(1:ngv)
        !
        ! Exchange wscr2
        !
        key1(1:ngv) = cmplx(0d0, 0d0)
        key2(1:ngv) = cmplx(0d0, 0d0)
        !
        if(dsp + 1 <= ig  .and. ig  <= dsp + cnt) key1(1:ngv) = wscr2(1:ngv,ig)
        if(dsp + 1 <= jg  .and. jg  <= dsp + cnt) key2(1:ngv) = wscr2(1:ngv,jg)
        !
        call MPI_allREDUCE(MPI_IN_PLACE, key1, ngv, &
        &                  MPI_DOUBLE_COMPLEX, MPI_SUM,MPI_COMM_WORLD,ierr)
        !
        call MPI_allREDUCE(MPI_IN_PLACE, key2, ngv, &
        &                  MPI_DOUBLE_COMPLEX, MPI_SUM,MPI_COMM_WORLD,ierr)
        !
        if(dsp + 1 <= ig  .and. ig  <= dsp + cnt) wscr2(1:ngv,ig) = key2(1:ngv)
        if(dsp + 1 <= jg  .and. jg  <= dsp + cnt) wscr2(1:ngv,jg) = key1(1:ngv)
        !
        ! Ordinally Gauss-Jordan
        !
        if(dsp + 1 <= ig .and. ig <= dsp + cnt) then
           !
           piv = cmplx(1d0, 0d0) / wscr1(ig, ig)
           !
           if(abs(piv) < 1d-12) then
              write(*,*) "Stop in invert. Singular imf = ", imf
              call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
              call MPI_FINALIZE(ierr)
              stop
           end if
           !
           call zscal(ngv, piv, wscr1(1:ngv, ig), 1)
           call zscal(ngv, piv, wscr2(1:ngv, ig), 1)
           call zcopy(ngv, wscr1(1:ngv, ig), 1, key1, 1)
           call zcopy(ngv, wscr2(1:ngv, ig), 1, key2, 1)
           !  
        else
           !
           key1(1:ngv) = cmplx(0d0, 0d0)
           key2(1:ngv) = cmplx(0d0, 0d0)
           !
        end if ! if(dsp + 1 <= ig .and. ig <= dsp + cnt) 
        !
        call MPI_allREDUCE(MPI_IN_PLACE, key1, ngv, &
        &                  MPI_DOUBLE_COMPLEX, MPI_SUM,MPI_COMM_WORLD,ierr)
        !
        call MPI_allREDUCE(MPI_IN_PLACE, key2, ngv, &
        &                  MPI_DOUBLE_COMPLEX, MPI_SUM,MPI_COMM_WORLD,ierr)
        !
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP & SHARED(dsp,cnt,ig,wscr1,wscr2,ngv,key1,key2) &
        !$OMP & PRIVATE(jg,piv)
        !
        !$OMP DO
        do jg = dsp + 1, dsp + cnt
           !
           if(jg == ig) cycle
           !
           piv = - wscr1(ig, jg)
           call zaxpy(ngv, piv, key1, 1, wscr1(1:ngv, jg), 1)
           call zaxpy(ngv, piv, key2, 1, wscr2(1:ngv, jg), 1)
           !
        end do ! jg
        !$OMP END DO
        !$OMP END PARALLEL
        !
     end do ! ig
     !
     wscr(1:ngv,dsp + 1:dsp + cnt, imf) = &
     &  wscr2(1:ngv,dsp + 1:dsp + cnt)
     !
  end do ! imf
  !
  deallocate(wscr1,wscr2)
  !
end subroutine invert
!
! Calc. Kel
!
subroutine make_Kel()
  !
  use mpi, only : MPI_STATUS_SIZE, MPI_DOUBLE_COMPLEX, MPI_IN_PLACE, MPI_COMM_WORLD, &
  &               MPI_DOUBLE_PRECISION, MPI_SUM
  use rpa_el_vals, only : nmf, nf, nk, nb, wfc1, wfc2, Kel, wscr, nftot, pi, gq2, cellV, &
  &                       my_rank, petot, ngv, gindx, FFTW_FORWARD, FFTW_ESTIMATE, nkpe
  !use rpa_routines, only : zgemm, zdotc, cnt_and_dsp, cnt_and_dsp_full
  !
  integer :: cnt, dsp, ik, ib, jb, imf, ierr, org, ipe, dst, src, &
  &          kcnt(0:petot - 1), kdsp(0:petot - 1), status(MPI_STATUS_SIZE)
  integer(8) :: plan
  complex(8) :: rhin(nftot), rhout(nftot), rho1(ngv,nb), rho3(ngv,nb), Kel0, &
  &             one = cmplx(1d0, 0d0), zero = cmplx(0d0, 0d0)
  complex(8),allocatable :: rho2(:,:)
  !
  call cnt_and_dsp(ngv, cnt, dsp)
  call cnt_and_dsp_full(nk, kcnt, kdsp)
  !
  allocate(Kel(0:nmf+1,nb,nb,nk))
  Kel(0:nmf + 1,1:nb,1:nb,1:nk) = 0d0
  !
  call dfftw_plan_dft_3d(plan, nf(1), nf(2), nf(3), rhin(1:nftot), &
  &                      rhout(1:nftot), FFTW_FORWARD, FFTW_ESTIMATE)
  !
  dst = modulo(my_rank + 1, petot)
  src  = modulo(my_rank - 1, petot)
  !
  do ipe = 1, petot
     !
     call MPI_SENDRECV_REPLACE(wfc1, nftot * nb * nkpe, MPI_DOUBLE_COMPLEX, &
     &                         dst, 1, src, 1, MPI_COMM_WORLD, STATUS, ierr  )
     !
     call MPI_SENDRECV_REPLACE(wfc2, nftot * nb * nkpe, MPI_DOUBLE_COMPLEX, &
     &                         dst, 1, src, 1, MPI_COMM_WORLD, STATUS, ierr  )
     !
     org = modulo(my_rank - ipe, petot)
     !
     !$OMP PARALLEL DEFAULT(NONE) &
     !$OMP & SHARED(my_rank,nk,nb,nftot,nmf,wfc1,wfc2,plan,cnt,dsp,zero,one,wscr,Kel,cellV, &
     !$OMP &        ngv,gindx,kcnt,kdsp,org,gq2) &
     !$OMP & PRIVATE(ik,ib,jb,imf,rhin,rhout,rho1,rho2,rho3,Kel0)
     !
     allocate(rho2(dsp + 1:dsp + cnt,nb))
     !
     !$OMP DO
     do ik = 1, kcnt(org)
        !
        do ib = 1, nb
           !
           do jb = 1, nb
              !
              rhin(1:nftot) = wfc1(1:nftot,ib,ik) * wfc2(1:nftot,jb,ik)
              call dfftw_execute_dft(plan, rhin(1:nftot), rhout(1:nftot))
              !
              rho1(1:ngv,            jb) = rhout(gindx(1:ngv))
              rho2(dsp + 1:dsp + cnt,jb) = rho1(dsp + 1:dsp + cnt,jb)
              !
           end do ! jb = 1, nb
           !
           do imf = 0, nmf
              !
              call zgemm("N", "N", ngv, nb, cnt, &
              &          one, wscr(1:ngv, dsp + 1:dsp + cnt,  imf), ngv, &
              &               rho2(       dsp + 1:dsp + cnt, 1:nb), cnt, &
              &         zero, rho3(1:ngv,                    1:nb), ngv  )
              !
              do jb = 1, nb
                 !
                 Kel0 = zdotc(ngv, rho1(1:ngv,jb), 1, &
                 &                 rho3(1:ngv,jb), 1)
                 !
                 Kel(imf,jb,ib,kdsp(org) + ik) = dble(Kel0) / cellV
                 !
              end do ! jb = 1, nb
              !
           end do ! imf = 0, nmf
           !
           ! Infinite frequency -> Bare Coulomb
           !
           do jb = 1, nb
              !
              Kel(nmf + 1,jb,ib,kdsp(org) + ik) = sum(dble(rho2(dsp + 1:dsp + cnt,jb) &
              &                                    * conjg(rho2(dsp + 1:dsp + cnt,jb))) &
              &                                    /        gq2(dsp + 1:dsp + cnt))  / cellV
              !
           end do ! jb = 1, nb
           !
        end do ! ib = 1, nb
        !
     end do ! ik = 1, kcnt(org)
     !$OMP END DO
     !
     deallocate(rho2)
     !
     !$OMP END PARALLEL
     !
  end do ! ipe = 1, petot
  !
  call dfftw_destroy_plan(plan)
  !
  call MPI_allREDUCE(MPI_IN_PLACE, Kel, (nmf + 2) * nb * nb * nk, &
  &                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  !
  deallocate(wscr)
  !
end subroutine make_Kel
!
! Compute parameter for Chebyshev interpolation
!
subroutine Chebyshev_interpol()
  !
  use mpi
  use rpa_el_vals, only : pi, nmf, nb, nk, Kel
  !use rpa_routines, only : cnt_and_dsp
  !
  integer :: ik, ib, jb, imf, jmf, cnt, dsp, ierr
  real(8) :: Chev(0:nmf + 1,0:nmf + 1), Kel0
  !
  call cnt_and_dsp(nk, cnt, dsp)
  !
  Kel(0:nmf + 1,1:nb,1:nb,            1:dsp) = 0d0
  Kel(0:nmf + 1,1:nb,1:nb,dsp + cnt + 1:nk ) = 0d0
  !
  do imf = 0, nmf + 1
     do jmf = 0, nmf + 1
        chev(jmf,imf) = 2d0 / dble(nmf + 2) &
        &  * cos(dble((2 * imf + 1) * jmf) * pi / dble(2 * (nmf + 2)))
     end do ! jmf
  end do ! imf
  chev(0,0:nmf + 1) = chev(0,0:nmf + 1) * 0.5d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(dsp, cnt, nb, nmf, Chev, Kel) &
  !$OMP & PRIVATE(ik,ib,jb,Kel0)
  !
  !$OMP DO
  do ik = dsp + 1, dsp + cnt
     do ib = 1, nb
        do jb = 1, nb
           !
           Kel0 = Kel(0,jb,ib,ik)
           Kel(0,jb,ib,ik) = Kel(nmf + 1,jb,ib,ik)
           Kel(nmf + 1,jb,ib,ik) = Kel0
           !
           Kel(0:nmf + 1,jb,ib,ik) = matmul(chev(0:nmf + 1,0:nmf + 1), Kel(0:nmf + 1,jb,ib,ik))
           !
        end do  ! jb
     end do ! ib
  end do ! ik
  !$OMP END DO
  !$OMP END PARALLEL
  !
  call MPI_allREDUCE(MPI_IN_PLACE, Kel, (nmf + 2) * nb * nb * nk, &
  &                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  !
end subroutine Chebyshev_interpol
!
! Output to file
!
subroutine write_data(iq)
  !
  use rpa_el_vals, only : my_rank, nmf, nk, nb, Kel, ng, iqv
  !
  integer,intent(in) :: iq
  !
  integer :: fo = 20
  character(100) :: fname
  !
  if(my_rank == 0) then
     !
     write(fname,*) iq
     write(fname,'(3a)') "vel", trim(adjustl(fname)), ".dat"
     open(fo, file = trim(fname), form = "unformatted")
     !
     write(fo) ng(1:3)
     write(fo) nb
     write(fo) dble(iqv(1:3,iq)) / dble(2 * ng(1:3))
     write(fo) nmf + 2
     !
     write(fo) Kel(0:nmf + 1,1:nb,1:nb,1:nk)
     !
     close(fo)
     !
  end if
  !
  deallocate(Kel)
  !
end subroutine write_data
!
end module rpa_el_routines
!
! Main routine
!
program rpa
  !
  use mpi, only : MPI_COMM_WORLD, MPI_INIT, MPI_COMM_SIZE, MPI_COMM_RANK
  use omp_lib, only : OMP_GET_WTIME, OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
  use rpa_el_vals, only : my_rank, petot, iqdo, ltetra, laddxc
  use rpa_el_routines, only : stdin, read_file, irr_bz, get_wfcg, fft_mesh, &
  &                           fft_wfc, read_dmxc, apply_xc, make_scrn, invert, make_kel, &
  &                           alloc_Kel, tetra_type, prepare_q, write_data, Chebyshev_interpol
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
  if(my_rank == 0 .and. OMP_GET_THREAD_NUM() ==0) &
  & write(*,*) '  # of thread : ', OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  Read from STDIN  #####"
     write(*,*)
  end if
  !
  t1 = OMP_GET_WTIME()
  call stdin()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,'(a,f15.3,a)') "  Time : ", t2 - t1, " sec."
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  Read from data-file.xml  #####"
     write(*,*)
  end if
  !
  t1 = OMP_GET_WTIME()
  call read_file()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,'(a,f15.3,a)') "  Time : ", t2 - t1, " sec."
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  Compute irreducible BZ  #####"
     write(*,*)
  end if
  !
  t1 = OMP_GET_WTIME()
  call irr_bz()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,'(a,f15.3,a)') "  Time : ", t2 - t1, " sec."
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  Read gvectors.dat & evc.dat  #####"
     write(*,*)
  end if
  !
  t1 = OMP_GET_WTIME()
  call get_wfcg()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,'(a,f15.3,a)') "  Time : ", t2 - t1, " sec."
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  Compute grid for FFT  #####"
     write(*,*)
  end if
  !
  t1 = OMP_GET_WTIME()
  call fft_mesh()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,'(a,f15.3,a)') "  Time : ", t2 - t1, " sec."
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  FFT wfc(G) -> wfc(r)  #####"
     write(*,*)
  end if
  !
  t1 = OMP_GET_WTIME()
  call fft_wfc()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,'(a,f15.3,a)') "  Time : ", t2 - t1, " sec."
  !
  if(laddxc == 1) then
     !
     if(my_rank == 0) then
        write(*,*)
        write(*,*) "#####  Read dmuxc.dat  #####"
        write(*,*)
     end if
     !
     t1 = OMP_GET_WTIME()
     call read_dmxc()
     t2 = OMP_GET_WTIME()
     if(my_rank == 0) write(*,'(a,f15.3,a)') "  Time : ", t2 - t1, " sec."
     !
  end if
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  Set frequency grid  #####"
     write(*,*)
  end if
  !
  t1 = OMP_GET_WTIME()
  call alloc_Kel()
  if(ltetra == 1 .or. ltetra == 2) call tetra_type()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,'(a,f15.3,a)') "  Time : ", t2 - t1, " sec."
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  Compute K_el  #####"
     write(*,*)
  end if
  !
  t1 = OMP_GET_WTIME()
  call prepare_q(iqdo)
  call make_scrn(iqdo)
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,'(a,f15.3,a)') "    Time (make_scrn) : ", t2 - t1
  !
  t1 = OMP_GET_WTIME()
  call apply_xc()
  call invert()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,'(a,f15.3,a)') "    Time (   invert) : ", t2 - t1
  !
  t1 = OMP_GET_WTIME()
  call make_kel()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,'(a,f15.3,a)') "    Time ( make_kel) : ", t2 - t1
  !
  call Chebyshev_interpol()
  call write_data(iqdo)
  !
  t2 = omp_get_wtime()
  !
  if(my_rank == 0) then
     write(*,*) ""
     write(*,'(a,f15.3,a)') "  Total time : ", t2 - t0
     write(*,*) ""
     call date_and_time(values = date)
     write(*,*) ""
     write(*,'(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') &
     &        "   End time : ", &
     &     date(1), "/", date(2), "/", date(3), " ", &
     &     date(5), ":", date(6), ":", date(7)
     write(*,*) ""
     write(*,*) "#####  Done  #####"
     write(*,*) ""
  end if
  !
  call MPI_FINALIZE(ierr)
  !
end program rpa

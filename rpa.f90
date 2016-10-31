module rpa_vals
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
  & nmf,           & ! # of Matsubara freq. (if nmf>0, set to # of lines of fort.135)
  & petot,         & ! Total # of PEs
  & my_rank,       & ! Index of PE
  & nsym,          & ! # of symmetries
  & ng(3),         & ! k point grid
  & nk,            & ! ng(1)*ng(2)*ng(3). Total # of k points
  & ngd(3),        & ! Dense k point grid
  & nkd,           & ! ngd(1)*ngd(2)*ngd(3). equal to k-grid in {pref}.nscf.in
  & nk0,           & ! # of irreducible k points
  & nb,            & ! # of bands
  & nf(3),         & ! # grid for FFT
  & nftot,         & ! Total # of G
  & igmin(3),      & ! Min. indices of G vector
  & npwmax,       & ! Max. # of PWs
  & gzero,        & ! index of G=0 vector (2016/01/18)
  & fbep,         & ! first band for el-ph
  & lbep,         & ! last band for el-ph
  & nbep,         & ! # of bands for el-ph
  & nx              ! # of energy scale, default value is 100.
  !
  real(8),save :: &
  & wcut,       & ! Cutoff energy for chi [Ry]
  & wlsm(4,20), & ! Weight for tetrahedron method
  & cellV,      & ! Unit cell volume [a.u.^3]
  & alat,       & ! Lattice constant [a.u.]
  & bvec(3,3),  & ! Reciplocal latticee vectors [/a.u.]
  & emin          ! Minimum energy scale [Ry]
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
  & eigval(:,:,:,:),& ! (nb,ngd(1),ngd(2),ngd(3)) Kohn-Sham energy [Ry]
  & gq2(:),         & ! (nftot) |G+q|^2
  & mf(:),          & ! (nmf) Matsubara frequencies
  & wmf(:),         & ! (nmf) Weights for Matsubara frequency
  & Kel(:,:,:,:),   & ! (nmf,nb,nb,nk). Coulomb Kernel
  & Ksf(:,:,:,:),   & ! (nmf, nb, nb, nk) Spin Fluctuation Kernel
  & xi0(:),         & ! (nx) energy scale [Ry]
  & dx0(:)            ! (nx) weight for energy
  !
  complex(8),allocatable,save :: &
  & dmxc(:,:,:),    & ! (gminxc(1):gmaxxc(1), ...) derivative of XC potential
  & wscr(:,:,:),    & ! (nftot,nftot,0:nmf) Screened interaction
  & wfc1(:,:,:),    & ! (nftot,nb,nk) wfc(ib,ik)
  & wfc2(:,:,:),    & ! (nftot,nb,nk) wfc(jb,jk)
  & wfc(:,:,:,:),   & ! (npwmax,nb,nk,2) wfc(G)
  !
  ! added by Kentaro Tsutsumi
  !
  & Ixc(:,:,:),     & ! (gminxc(1):gmaxxc(1), ...) derivative of XC potential 
                      ! with respect to local spin density
  & Chi_s(:,:,:),   & ! (nftot,nftot,0:nmf) spin susceptibility
  & Pi_s(:,:,:),      & ! (nftot,nftot,0:nmf) true polarization
  & Vscr(:,:,:),    & ! (nftot,nftot,0:nmf) screened Coulomb interaction
  & Wsf(:,:,:)        ! (nftot, nftot, 0:nmf) spin-mediated-interaction energy
  !
  include 'fftw3.f'
  !
end module rpa_vals
!
module rpa_routines
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
     subroutine zgetrf ( m, n, a, lda, ipiv, info)
       integer      m, n, lda, info, ipiv(*)
       complex*16   a(lda, *)
     end subroutine
  END INTERFACE
  !
contains
!
! Read from STDIN
!
subroutine stdin()
  !
  use mpi, only : MPI_INTEGER, MPI_COMM_WORLD, MPI_DOUBLE_PRECISION
  use rpa_vals, only : my_rank, ng, iqdo, nmf, ngd, nkd, ltetra, laddxc, wcut
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
     write(*,*) "                   k grid : ", ng(1:3)
     write(*,*) "             Dense k grid : ", ngd(1:3)
     write(*,*) "     # of Matsubara freq. : ", nmf
     write(*,*) "                  q index : ", iqdo
     write(*,*) "                   ltetra : ", ltetra
     write(*,*) "                   laddxc : ", laddxc
     write(*,*) "    Cutoff kinetic energy : ", wcut
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
  use rpa_vals, only : my_rank, alat, bvec, nb, nk, ng, npwmax, &
  &                    cellV, nsym, sym, pi, eig, nkd, eigval, ngd
  !
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
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !   modified by kentaro tsutsumi (2016/01/18)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     allocate(eigval(nb, ngd(1), ngd(2), ngd(3)))
     !
     open(fi, file = "eigval.dat")
     read(fi,*) eigval(1:nb, 1:ngd(1), 1:ngd(2), 1:ngd(3))
     close(fi)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !  end of modified region
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  use rpa_vals, only : my_rank, nk, ng, nk0, iqv, grid, nsym, sym
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
           if(ng(3) /= 1) then
             qv1(1:3) = (dble((/i1, i2, i3/)) - 0.5d0) / dble(ng(1:3))
           else
             qv1(1) = (dble(i1) - 0.5d0) / dble(ng(1))
             qv1(2) = (dble(i2) - 0.5d0) / dble(ng(2))
             qv1(3) = (dble(i3) - 1.0d0) / dble(ng(3))  ! modified to calc. in case of FeSe
           end if
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
!  allocate(iqv(3,nk0))
  allocate(iqv(3,nk0+1)) ! added by Kentaro Tsutsumi (2015/12/11)
  do ik = 1, nk0
     !
     qv0(1:3,ik) = qv0(1:3,ik) * dble(2 * ng(1:3))
     iqv(1:3,ik) = nint(qv0(1:3,ik))
     !
  end do
  iqv(1:3, nk0 + 1) = 0
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
  use rpa_vals, only : my_rank, pi, alat, bvec, nk, npwmax, nb, npw, igv, wfc, grid, ng, nkpe
  !use rpa_routines, only : cnt_and_dsp
  !
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
  use rpa_vals, only : petot, my_rank
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
  use rpa_vals, only : petot, my_rank
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
  use rpa_vals, only : my_rank, nk, npwmax, igv, nf, nftot, igmin, nkpe
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
  use rpa_vals, only : my_rank, nk, nb, npw, wfc, wfc1, wfc2, igv, igmin, nf, nkpe, nftot, &
  &                    FFTW_BACKWARD, FFTW_ESTIMATE
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
  use rpa_vals, only : my_rank, dmxc, gminxc, gmaxxc
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
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! modified in order to match the format of dmuxc.dat(2016/10/13)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !
     read(fi, *) 
     !
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! modified in order to match the format of dmuxc.dat(2016/10/13)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  use rpa_vals, only : nk, nb, nf, nmf, mf, nftot, gq2, gindx, Kel, wmf
  !use rpa_routines, only : weightspoints_gl
  !
  integer :: imf
  real(8) :: dummy, dummy1
  real(8),parameter :: H2Ry = 2.0d0
  !
  allocate(mf(nmf), wmf(nmf), gq2(nftot), gindx(nftot))
  !
!  call weightspoints_gl(nmf,mf,wmf)
  !
!  do imf = 1, nmf
!    !
!    wmf(imf) = wmf(imf) * 2d0 / (1d0 - mf(imf))**2
!    mf(imf) = (1d0 + mf(imf)) / (1d0 - mf(imf))
!    !
!  end do
  !
  open(10, file = "fort.135")
  read(10, *) dummy, dummy
  do imf = 1, nmf
    read(10, *) dummy, mf(imf)
  end do
  close(10)  
  !
  mf(1:nmf) = mf(1:nmf) * H2Ry ! set to Rydberg
  !
end subroutine alloc_Kel
!
! Weights & Points for Gauss-Legendre method
!
subroutine weightspoints_gl(n,x,w)
  !
  use mpi, only : MPI_COMM_WORLD
  use rpa_vals, only: pi
  !use rpa_routines, only : legendre
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
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        call MPI_FINALIZE(ierr)
        stop "error newton"
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
  use mpi, only : MPI_COMM_WORLD
  USE rpa_vals, ONLY : my_rank, wlsm, ivvec, ng, ltetra, bvec
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
  use rpa_vals, only : pi, nk, nb, grid, igmin, ng, bvec, wcut, ngv, gindx, nkpe, &
  &                  nf, nftot, iqv, wfc1, wfc2, gq2, my_rank, petot, gzero
  !use rpa_routines, only : cnt_and_dsp_full
  !
  integer,intent(in) :: iq
  !
  integer :: ik, jk, ib, i1, i2, i3, ikv(3), jkv(3), g0(3), ir(3), ifft, ig, dst, src, ierr, &
  &          cnt(0:petot - 1), dsp(0:petot - 1), jkindx(nk), status(MPI_STATUS_SIZE), org, ipe
  real(8) :: gv(3), qv(3), theta, gq20
  complex(8) :: phase(nftot), wfctmp(nftot,nb,nkpe)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! modified by Kentaro Tsutsumi. !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  character(100) :: fname
  real(8) :: tmp(3)
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
           tmp(1:3) = gv(1:3)
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
              if(dot_product(tmp(1:3), tmp(1:3)) .lt. 1.0d-10) then
                gzero = ngv
              end if
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
  use rpa_vals, only : my_rank, nftot, nf, nk, nb, nmf, wfc1, wfc2, wscr, gq2, nkpe, &
  &                    FFTW_FORWARD, FFTW_ESTIMATE, my_rank, petot, ngv, gindx, gzero
  !use rpa_routines, only : zgemm, cnt_and_dsp, cnt_and_dsp_full
  !
  integer,intent(in) :: iq
  !
  integer :: ik, ig, jg, ib, jb, imf, cnt, dsp, fstg, lstg, org, dst, src, ierr, ipe, &
  &          kcnt(0:petot - 1), kdsp(0:petot - 1), status(MPI_STATUS_SIZE)
  integer(8) :: plan
  real(8) :: t1, t2
  !
  complex(8) :: wght(0:nmf,nb,nb,nk), rho1(ngv,nb), rhin(nftot), rhout(nftot)
  complex(8) :: one = cmplx(1d0, 0d0)
  complex(8) :: zero = cmplx(0d0, 0d0)
  complex(8),allocatable :: rho2(:,:)
  character(40) :: fo
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
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! output Chi_0(G=G', q=0)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
!  write(fo, *) my_rank
!  write(fo, '(3a)')  "Chi_0_before", trim(adjustl(fo)), ".dat"
!  open(223, file = trim(fo))
!  do jg = dsp + 1, dsp + cnt
!    do ig = 1, ngv
!      if(ig == jg) then
!        write(223, *) gq2(ig), real(wscr(ig, jg, 0))
!      end if
!    end do
!  end do
!  close(223)
end subroutine make_scrn
!
! Calculation of weight function (f(1-f')) / (e - e')
!
subroutine fermi_fuctor(iq,wght)
  !
  use mpi, only : MPI_IN_PLACE, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD
  use rpa_vals, only : nb, nk, nkd, nmf, ng, ngd, iqv, ltetra, ivvec, xi0, nx, cellV, nk0, my_rank
  !use rpa_routines, only : tetraweight, interpol_weight
  !
  integer,intent(in) :: iq
  complex(8),intent(out) :: wght(0:nmf,nb,nb,nk)
  complex(8),allocatable :: dos(:,:,:,:) ! (nx, 1;nb, nk, 2)
  !
  integer :: cnt, dsp, nkd0, nt, it, ik, i1, i2, i3, ii, iqvd(3), ierr, ikv(3), &
  &          dgrid(3,nkd), indx1(2,20, 6 * nkd), indx2(20, 6 * nkd), indx3(20 * 6 * nkd)
  integer :: ix, ib, imf
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
  ! modify!!!
  ! calc. DOS [/Ryd/cell]
  !
!  wghtd(0:nmf, 1:nb, 1:nb, 1:nkd0) = cmplx(0d0, 0d0)
!  call calc_dosk(nkd0, indx1, indx2, wghtd)
!  open(10, file = "dos.dat")
!  write(10, *) "DOS =", sum(wghtd(0,1:nb, 1:nb, 1:nkd0)) * cellV, "[/Ry]"
!  close(10)
  !
  wghtd(0:nmf, 1:nb, 1:nb, 1:nkd0) = cmplx(0d0, 0d0)
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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! modified by kentaro tsutsumi
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  allocate(dosd(1:nb, 1:nkd))
!  !
!  call calc_dosk(dosd)
!  !
!  do ib = 1, nb
!    do imf = 0, nmf
!      wghtd(imf, ib, ib, 1:nkd0) = wghtd(imf, ib, ib, 1:nkd0) + dosd(ib, 1:nkd0)
!    end do
!  end do
  if(iq == nk0+1) then
    call calc_dosk(nkd0, indx1, indx2, wghtd)
  end if
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! end of modified region
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Interpolation of weight
  !
  wght(0:nmf,1:nb,1:nb,1:nk) = cmplx(0d0, 0d0)
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
!  if(my_rank == 0) then
!    write(*, *) "end of calling tetraweight"
!  end if
end subroutine fermi_fuctor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  codes given by kawamura-san
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! Compute Dos_k, xmin, xmax
      ! 
      subroutine calc_dosk(nkd0,indx1,indx2,wghtd)
        !
        use mpi, only : MPI_IN_PLACE, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD
        USE rpa_vals, ONLY : nkd, nmf, nb, wlsm, eig, cellV
        !use scdft_routines, only : cnt_and_dsp, sort
        !
        integer,intent(in) :: nkd0, indx1(2,20,6 * nkd), indx2(20,6 * nkd)
        complex(8),intent(inout) :: wghtd(0:nmf,nb,nb,nkd0)
        !
        integer :: ib, it, ii, dsp, cnt, ierr
        real(8) :: ei(4,nb), e(4), a(4,4), w0(4,4), w1(4), tmp(5,4), V, scal
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! modified by Kentaro Tsutsumi
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        complex(8) :: wghtdtmp(0:nmf,nb,nb,nkd0)
        !
        wghtdtmp(:,:,:,:) = cmplx(0d0, 0d0)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! end of modified region
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        scal = 2d0 / (dble(6 * nkd) * cellV)
        !
        call cnt_and_dsp(nkd * 6,cnt,dsp)
        !
        w0(1:4,1:4) = 0d0
        do ii = 1, 4
           w0(ii,ii) = 1d0
        end do
        !
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP & SHARED(scal,dsp,cnt,wlsm,eig,w0,nb,wghtdtmp,indx1, indx2) &
        !$OMP & PRIVATE(it,ii,ib,ei,w1,a,e,tmp,V)
        !
        do it = dsp + 1, dsp + cnt
           !
           ei(1:4,1:nb) = 0.0d0
           do ii = 1, 20
              !
              do ib = 1, nb
                 ei(1:4,ib) = ei(1:4,ib) + wlsm(1:4,ii) * eig(ib,indx1(1,ii,it))
              end do
              !
           end do
           !
           !$OMP DO
           do ib = 1, nb
              !
              w1(1:4) = 0d0
              !
              tmp(  1,1:4) = ei(1:4,ib)
              tmp(2:5,1:4) = w0(1:4,1:4)
              call sort(5, 4, tmp)
              e(1:4) = tmp(1,1:4)
              !
              do ii = 1, 4
                 a(ii,1:4) = (0d0 - e(1:4)) / (e(ii) - e(1:4))
              end do
              !
              if(e(1) < 0d0 .and. 0d0 <= e(2)) then
                 !
                 ! A
                 !
                 V = a(2,1) * a(3,1) * a(4,1) / (0d0 - e(1))
                 !
                 w1(1:4) = V * ( tmp(2:5,1) * a(1,2) + tmp(2:5,2) * a(2,1) &
                 &             + tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1) &
                 &             + tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1) )
                 !
              else if( e(2) < 0d0 .and. 0d0 <= e(3)) then
                 !
                 ! B - 1
                 !
                 V = a(3,1) * a(4,1) * a(2,4) / (0d0 - e(1))
                 !
                 w1(1:4) = V * ( tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1) &
                 &             + tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1) &
                 &             + tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) )
                 !
                 ! B - 2
                 !
                 V = a(2,3) * a(3,1) * a(4,2) / (0d0 - e(1))
                 !
                 w1(1:4) = w1(1:4)                                         &
                 &       + V * ( tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1) &
                 &             + tmp(2:5,2) * a(2,3) + tmp(2:5,3) * a(3,2) &
                 &             + tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) )
                 !
              else if(e(3) < 0d0 .and. 0d0 < e(4)) then
                 !
                 ! C
                 !
                 V = a(1,4) * a(2,4) * a(3,4) / (e(4) - 0d0)
                 !
                 w1(1:4) = V * ( tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1) &
                 &             + tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) &
                 &             + tmp(2:5,3) * a(3,4) + tmp(2:5,4) * a(4,3) )
                 !
              end if
              !
              do ii = 1, 20
                 !
                 wghtdtmp(0,ib,ib,indx2(ii,it)) = wghtdtmp(0,ib,ib,indx2(ii,it)) &
                 &          + cmplx(sum(wlsm(1:4,ii) * w1(1:4)) * scal, 0d0)
                 !
              end do ! ii = 1, 20
              !
           end do ! ib
           !$OMP END DO NOWAIT
           !
        end do ! it
        !$OMP END PARALLEL
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! modified by Kentaro Tsutsumi
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        call MPI_allREDUCE(MPI_IN_PLACE, wghtdtmp, (nmf + 1) * nb * nb * nkd0, &
!        &                  MPI_DOUBLE_COMPLEX, MPI_SUM,MPI_COMM_WORLD,ierr)
      !
        wghtd(0:nmf,1:nb,1:nb,1:nkd0) = wghtd(0:nmf,1:nb,1:nb,1:nkd0) + wghtdtmp(0:nmf,1:nb,1:nb,1:nkd0)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! end of modified region
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end subroutine calc_dosk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! end of codes given by kawamura-san
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Compute energy grid
!
subroutine energy_grid()
  !
  use rpa_vals, only : eigval, ngd, xi0, dx0, nx, emin, fbep, lbep, nb
  !use scdft_routines, only : weightspoints_gl
  !
  integer :: ix
  real(8) :: xmax, xmin, xx, ww, rhs(nx)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  modified by kentaro tsutsumi
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  emin = 1d-7
  fbep = 1
  lbep = 6
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  end of modified region
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Make logscale
  !
  xmax = maxval(eigval(1:nb,1:ngd(1),1:ngd(2),1:ngd(3)))
  xmin = minval(eigval(1:nb,1:ngd(1),1:ngd(2),1:ngd(3)))
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
! Integration weight with tetrahedron method
!
subroutine tetraweight(nkd0,dgrid,indx1,indx2,iqvd,wghtd)
  !
  USE rpa_vals, ONLY : nkd, nb, nmf, ngd, ivvec, wlsm, eig, cellV
  !use rpa_routines, only : cnt_and_dsp, sort, tetra2
  !
  !
  integer,intent(in) :: nkd0, dgrid(3,nkd), indx1(2,20,6 * nkd), indx2(20,6 * nkd), iqvd(3)
  complex(8),intent(out) :: wghtd(0:nmf,nb,nb,nkd0)
  !
  integer :: ik, it, ib, jb, imf, ii, cnt, dsp, ix
  real(8) :: thr = 1d-8, V
  real(8) :: e(4), a(4,4), ei(nb,4), ej(nb,4), ei2(4), ej2(nb,4), &
  &        tmp(10,0:nmf,nb,4), tmp2(10,0:nmf,nb,4), &
  &         w0(4,2,0:nmf,nb,4), w1(4,2,0:nmf,nb), w2(4,2,0:nmf,nb,4)
  !
  call cnt_and_dsp(nkd * 6,cnt,dsp)
  !
  wghtd(0:nmf,1:nb,1:nb,1:nkd0) = 0.0
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
  USE rpa_vals, ONLY : nb, nmf
  !use rpa_routines, only : sort, lindhard
  !
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
  USE rpa_vals, ONLY : nmf, mf
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
  use rpa_vals, only : ivvec, ltetra
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
  use rpa_vals, only : nmf, nftot, igmin, nf, wscr, gq2, &
  &                  dmxc, gminxc, gmaxxc, my_rank, laddxc, ngv, gindx
  !use rpa_routines, only : cnt_and_dsp
  !
  integer :: i1, i2, i3, ig, jg, imf, idg(3), igv0(3,nftot), cnt, dsp
  complex(8) :: vec(ngv), one = cmplx(1d0, 0d0), zero = cmplx(0d0, 0d0)
  complex(8),allocatable :: wscr0(:,:,:)
  !
  call cnt_and_dsp(ngv, cnt, dsp)
!  write(*,*) "ngv, cnt, dsp=", ngv, cnt, dsp
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              Added by Kentaro Tsutsumi.                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
!  Subroutines added by Kentaro Tsutsumi computes
!   - interacting spin susceptibility Chi_s
!   - spin-fluctuation-mediated interaction energy Wsf
!   - spin fluctuation kernel Ksf
!  using data written in Ixc.dat (generated by executing ph.x)
!
subroutine read_Ixc()
!
!  read Ixc_G from Ixc.dat
!
  use mpi, only : MPI_INTEGER, MPI_COMM_WORLD, MPI_DOUBLE_COMPLEX
  use rpa_vals, only : my_rank, Ixc, gminxc, gmaxxc
  !
  !
  integer :: fi = 10, ierr
  !
  if(my_rank == 0) then
     !
     open(fi, file = "Ixc.dat")
     !
     read(fi,*) gminxc(1:3), gmaxxc(1:3)
     write(*,*) gminxc(1:3)
     write(*,*) gmaxxc(1:3)
     allocate(Ixc(gminxc(1):gmaxxc(1), gminxc(2):gmaxxc(2), gminxc(3):gmaxxc(3)))
     !
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! modified in order to match the format of Ixc.dat (2016/10/13)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !
     read(fi, *) 
     !
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! modified in order to match the format of Ixc.dat (2016/10/13)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !
     read(fi,'(2e25.15)') Ixc(gminxc(1):gmaxxc(1), gminxc(2):gmaxxc(2), gminxc(3):gmaxxc(3))
     !
     write(*,*) "Ixc(G=G'=0) =", Ixc(0, 0, 0)
     !
     close(fi)
     !
  end if
  !
  call MPI_BCAST(gminxc, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(gmaxxc, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if(my_rank /= 0) &
  &  allocate(Ixc(gminxc(1):gmaxxc(1), gminxc(2):gmaxxc(2), gminxc(3):gmaxxc(3)))
  call MPI_BCAST(Ixc, product(gmaxxc(1:3) - gminxc(1:3) + 1),&
  &              MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! check whether Ixc was read correctly
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  if(my_rank == 0) then
     !
     open(223, file = "Ixc_read.dat")
     !
     write(223,*) gminxc(1:3)
     write(223,*) gmaxxc(1:3)
     !
     write(223,'(2e25.15)') Ixc(gminxc(1):gmaxxc(1), gminxc(2):gmaxxc(2), gminxc(3):gmaxxc(3))
     !
     close(223)
     !
  end if
end subroutine read_Ixc
!
!  prepare for calc. Chi_s using LU decomposition
!
subroutine LUdecomp(N, a, indx)
!
  integer, intent(in) :: N
  complex(8),intent(inout) :: a(N,N)
  integer, intent(out) :: indx(N)
  complex(8), parameter :: small =( 1.0d-20, 0d0)
  complex(8) :: s, s1, s2, d
  real(8) :: big, dum, vv(N)
  integer :: i, j, imax, k
  !
  !  prepare for pivoting
  !
  do i = 1, N
    big = 0.0d0
    do j = 1, N
      if ( abs(a(i, j)*conjg(a(i,j))) > big ) then
        big = sqrt(abs(a(i, j)*conjg(a(i,j))))
      end if
    end do
    if (big < abs(small)) then
      stop " Singular matrix "
    end if
    vv(i) = 1.0d0 / big ! scale factor
  end do
  !
  ! execute LU decomposition using Crout's algorithm
  !
  do j = 1, N
    do i = 1, j
      !
      ! computes Upper triangle matrix element
      !
        s1 = a(i, j)
        !if (i > 1) then
          do k= 1, i-1
            s1 = s1 - a(i, k) * a(k, j)
          end do
        !end if
        a(i, j) = s1
        !
    end do ! i
    big = 0.0d0 ! prepare for searching max pivot
    !
    ! computes Lower triangle matrix element
    !
    do i = j+1, N
      s2 = a(i, j) 
      !if (j > 1) then
        do k = 1, j-1
          s2 = s2 - a(i, k) * a(k, j)
        enddo
      !end if
      a(i, j) = s2
    ! 
    ! check the goodness of pivot
    !
      dum = vv(i) * sqrt(abs(a(i, j)*conjg(a(i,j))))
      if (dum >= big) then
        big = dum
        imax = i
      endif
      !
    end do ! i
    !
    ! exchange the column as necessary
    !
    if (j /= imax) then
      do k = 1, N
        dum = a(imax, k)
        a(imax, k) = a(j, k)
        a(j, k) = dum
      end do
      vv(imax) = vv(j)  ! exchange the scale factor
    endif
    !
    indx(j) = imax
    if (a(j, j) == 0.0d0) a(j, j) = small ! give small finite value
    !
    ! divide by the pivot 
    !
    if(j /= N) then
      d = 1.0d0 / a(j, j)
      do i = j+1, N
        a(i, j) = a(i, j) * d
      end do ! i
    end if
  !
  end do ! j
  !
end subroutine LUdecomp
!
! subroutine to calc Chi_s after LU decomposition
!
subroutine LUsolve(a, N, indx, b)
  implicit none
  !
  integer, intent(in) :: N
  complex(8), intent(in) :: a(N,N)
  integer, intent(in) :: indx(N)
  complex(8), intent(inout) :: b(N)
  !
  integer :: i, ip, j
  complex(8) :: s
  !
  do i = 1, N
    ip = indx(i)
    s = b(ip)
    b(ip) = b(i)
      if (i == 1) then 
        s = s
      else
        do j = 1, i-1
          s = s - a(i, j) * b(j) ! forward substitution
        end do
      end if
    b(i) = s
  end do ! i
  !
  do i = N, 1, -1
    s = b(i)
    do j = i+1, N
      s = s - a(i, j) * b(j) ! backward substitution
    end do
    b(i) = s/a(i, i)
  end do
  !
end subroutine LUsolve
!
!  computes Chi_s
!
subroutine make_Chi_s(iq)
!
! Computes spin susceptibility Chi_s by solving the matrix equation
!             Chi_s = Chi_0 + Chi_0 * Ixc * Chi_s
! using the subroutines "LUdecomp" and "LUsolve"
!
  use mpi, only : MPI_IN_PLACE, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD 
  use rpa_vals, only : nmf, nftot, igmin, nf, wscr, gq2, cellV, dmxc, gzero, &
  &                  Ixc, gminxc, gmaxxc, my_rank, laddxc, ngv, gindx, Chi_s, iqv, cellV, nk0
  !use rpa_routines, only : cnt_and_dsp
  !
  integer, intent(in) :: iq ! to compare the Chi_s to the exp. value
  !
  integer :: i1, i2, i3, ig,ig1,  jg, imf, idg(3), igv0(3,nftot), cnt, dsp, ierr, info
  complex(8),parameter :: one = cmplx(1d0, 0d0), zero = cmplx(0d0, 0d0)
  complex(8),allocatable :: wscr1(:,:,:), E(:,:), tmp(:), tmp0(:, :), vec(:), tmp4(:,:), wscrn2(:,:,:)
  integer,allocatable :: indx(:) !(ngv) used in the LU decomposition
  character(100) :: fo
  real(8) :: g2
  real(8) :: tmp1, tmp2, tmp3
  !
  call cnt_and_dsp(ngv, cnt, dsp)
  !
  allocate(Chi_s(1:ngv,dsp + 1:dsp + cnt,0:nmf))
  allocate(wscr1(1:ngv,1:ngv,0:nmf))
!  allocate(wscrn2(1:ngv,dsp+1:dsp+cnt,0:nmf))
  allocate(vec(dsp+1:dsp+cnt))
  allocate(tmp(ngv))
  allocate(indx(ngv))
  allocate(tmp0(ngv, ngv))
  !
!  wscrn2(:,:,:) = (-1d0) * wscr(:,:,:) ! modified (2016/09/07)
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
     wscr1(:,:,:) = zero
     do ig = 1, ngv ! index for G2
        !
        do jg = dsp+1, dsp+cnt ! index for G1
           !
           idg(1:3) = igv0(1:3,gindx(ig)) - igv0(1:3,gindx(jg)) ! G2 - G1
           !
           if(all(gminxc(1:3) <= idg(1:3)) .and. all(idg(1:3) <= gmaxxc(1:3))) then
             vec(jg) = Ixc(idg(1), idg(2), idg(3))
           else
             vec(jg) = zero
           end if
           !
        end do ! jg
!      end do ! ig
        !
        !  calc matrix Chi_0*Ixc
        !
        do imf = 0, nmf
!       call zgemm("N","N", ngv, ngv, cnt,                         &
!     &          one, wscr(1:ngv, dsp + 1:dsp + cnt, imf), ngv,   &
!     &                vec(       dsp + 1:dsp + cnt, 1:ngv),cnt, &
!     &        zero, wscr1(1:ngv, 1:ngv, imf), ngv)
          call zgemv("N", ngv, cnt, one, (-1d0)*wscr(1:ngv, dsp+1:dsp+cnt, imf), &
          &          ngv, vec(dsp+1:dsp+cnt), 1, zero, wscr1(1:ngv, ig, imf), 1)
        end do ! imf
      end do ! ig
      !
!      do imf = 0, nmf
!        tmp4(1:ngv, 1:ngv) = wscr1(1:ngv, 1:ngv, imf)
!      call MPI_allREDUCE(MPI_IN_PLACE, tmp4, ngv * ngv, &
!      &                  MPI_DOUBLE_COMPLEX, MPI_SUM,MPI_COMM_WORLD,ierr) ! nokosu
!        wscr1(1:ngv, 1:ngv, imf) = tmp4(1:ngv, 1:ngv)
!      end do
      call MPI_allREDUCE(MPI_IN_PLACE, wscr1, (nmf + 1) * ngv * ngv, &
      &                  MPI_DOUBLE_COMPLEX, MPI_SUM,MPI_COMM_WORLD,ierr) ! nokosu
           ! computes matrix multiple Pi_0*Ixc ( defined as wscr0 ).
           ! Here, "wscr" means Chi_0 because this subroutine is called
           ! soon after make_scrn(iq) is called and before apply_xc() is called.
     !
     allocate(E(ngv, ngv))
     E(:, :) = zero
     do ig = 1, ngv
       E(ig, ig) = one  ! make Unit matrix(ngv, ngv)
     end do
     do imf = 0, nmf
       wscr1(1:ngv, 1:ngv, imf) = E(1:ngv, 1:ngv) - wscr1(1:ngv, 1:ngv, imf)
     end do !imf
     !
     !
     tmp(1 : ngv) = zero
     do imf = 0, nmf
       !
       !  LU decomposition of wscr0(= E - Chi_0 * Ixc )
       !
       tmp0(1:ngv,1:ngv) = wscr1(1:ngv, 1:ngv, imf)
       !call LUdecomp(ngv, tmp0, indx)
       call zgetrf(ngv, ngv, tmp0, ngv, indx, info)
       !
       do i1 = 1, cnt
         !
         !  solve the matrix equation 
         !
         tmp(1:ngv) = (-1d0)*wscr(1:ngv, dsp + i1, imf)
         !tmp(1:ngv) = wscrn2(1:ngv, dsp + i1, imf)
         !
         call LUsolve(tmp0, ngv, indx, tmp)
         !
         Chi_s(1:ngv, dsp+i1, imf) = tmp(1:ngv)
       end do ! i1
     end do ! imf
     !
     !deallocate(vec,tmp, tmp0, indx, wscr1)
     deallocate(vec,tmp, tmp0, indx, wscr1)
     !
     !  check Chi_s is consistent to exp. value
     !
     if (iq == nk0+1) then ! gamma point
     do ig = 1, ngv
      if(ig == gzero) then
!       g2 = dble(dot_product(igv0(1:3,gindx(ig)), igv0(1:3, gindx(ig)))) 
!       if(g2 <= 1.0d-8 ) then  ! valid only when G = 0
        do i1 = dsp + 1, dsp + cnt
         if( i1 == ig ) then
          write(fo,*) gzero
          write(fo,'(3a)') "Chi_s", trim(adjustl(fo)), ".dat"
          open(10, file=trim(fo))
!!          write(10, *) "0 0 ", dble(iq)/dble(50)
!!          do imf = 0, nmf
           write(10, *) "Chi_s from matrix =", real(Chi_s(ig, i1, 0)) * cellV
           !write(10, *) "Chi_0 =", real(wscrn2(ig, i1, 0)) * cellV
           write(10, *) "Chi_0 =", real(wscr(ig, i1, 0)) * cellV * (-1d0)
           write(10, *) "Stoner factor S(matrix) =", (-1d0)*real(Chi_s(ig, i1, 0)) / real(wscrn2(ig, i1, 0))
!          tmp1 = real(wscr(ig, i1, 0))
!          tmp2 = real(Ixc(igv0(1,gindx(ig)), igv0(2,gindx(ig)), igv0(3,gindx(ig))))
!          tmp3 = tmp1 / (1d0 - tmp1 * tmp2)
!           write(10, *) "Chi_s from scalar =", tmp3
!           write(10, *) "Stoner factor S(scalar) =", 1d0/(1d0 - tmp1 * tmp2)
!           write(10, *) "Ixc [Ry] =", tmp2
!           write(10, *) "G =", igv0(1:3,gindx(ig))
!          write(10, *) "dmxc =", real(dmxc(igv0(1,gindx(ig)), igv0(2,gindx(ig)), igv0(3,gindx(ig))))
!          write(10, *) tmp3 ! Chi_s calculated from scalar
!          write(10, *) tmp1 ! Chi_0 
          !write(10, *) "                        Ixc  =", tmp2
          !write(10, *) igv0(1:3, gindx(ig))
!          end do ! imf
!         close(10)
         end if ! i1
       end do !i1
      end if ! g2
     end do ! ig
    end if ! iq
!    deallocate(wscrn2)
  end if ! laddxc
end subroutine make_Chi_s
!
!  Computes the spin-fluctuation-mediated interaction energy 
!       Wsf(q) = 3 * Ixc * Chi_s(q) * Ixc / cellV
!
subroutine make_Wsf()
  !
  use mpi, only : MPI_STATUS_SIZE, MPI_DOUBLE_COMPLEX, MPI_IN_PLACE, MPI_COMM_WORLD, &
  &               MPI_DOUBLE_PRECISION, MPI_SUM
  use rpa_vals, only : Chi_s, Ixc, gminxc, gmaxxc, Wsf, ngv, laddxc, nmf, nf, &
  &                     igmin, gindx, my_rank, nftot, cellV, mf
  !
  !integer, intent(in) :: iq
  integer :: i1, i2, i3, cnt, dsp, igp, ig2, ig1, ig, imf, ierr, counter
  complex(8) :: one = cmplx(1d0, 0d0), zero = cmplx(0d0, 0d0)
  complex(8), allocatable :: vec(:,:) ! (ngv, ngv)
  integer :: igv0(3,nftot)
  integer :: idg(3)
  complex(8),allocatable :: tmp1(:,:), tmp2(:,:), Chi_s_tmp(:,:,:)
  real(8) :: g2
  character(100) :: fname
  !
  call cnt_and_dsp(ngv, cnt, dsp)
  allocate(vec(1:ngv,1:ngv))
  allocate(tmp1(1:ngv, 1:ngv))
  allocate(tmp2(1:ngv, 1:ngv))
  allocate(Wsf(1:ngv, 1:ngv, 0:nmf))
  allocate(Chi_s_tmp(1:ngv, 1:ngv, 0:nmf))
  !
  Chi_s_tmp(:,:,:) = cmplx(0d0, 0d0)
  Chi_s_tmp(1:ngv, dsp + 1: dsp + cnt, 0:nmf) = Chi_s(1:ngv, dsp + 1 : dsp + cnt, 0:nmf)
  !
  call MPI_allREDUCE(MPI_IN_PLACE, Chi_s_tmp, ngv * ngv * (1+nmf), &
  &                   MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
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
     ! calc matrix multiple 
     ! tmp1(G1,G') = Chi_s_G1G2 * Ixc(G'-G2)
     !
     do ig2 = 1, ngv
        !
        do ig1 = 1, ngv
           !
           idg(1:3) = igv0(1:3,gindx(ig2)) - igv0(1:3,gindx(ig1)) 
           !
           ! make matrix Ixc
           !
           if(all(gminxc(1:3) <= idg(1:3)) .and. all(idg(1:3) <= gmaxxc(1:3))) then
              vec(ig1, ig2) = Ixc(idg(1), idg(2), idg(3))
           else
              vec(ig1, ig2) = zero
           end if
           !
        end do ! ig1
        !
     end do ! ig2
     !
     do imf = 0, nmf
        call zgemm("N","N", ngv, ngv, ngv,                              &
        &            one, Chi_s_tmp(1:ngv, 1:ngv,        imf),      &
        &            ngv,       vec(       1:ngv, 1:ngv     ), ngv, &
        &           zero,      tmp1(1:ngv,        1:ngv     ), ngv)
!         call zgemv("N", ngv, cnt, one, Chi_s(1:ngv, dsp+1:dsp+cnt, imf),  &
!         &          ngv, vec(dsp+1:dsp+cnt), 1, zero,                      &
!         &          tmp1(1:ngv, ig2, imf), 1)
!     end do
!     end do ! ig2
     !
!     call MPI_allREDUCE(MPI_IN_PLACE, tmp1, ngv * ngv , &
!     &                   MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
     !
     ! calc matrix multiple 
     ! Wsf(G,G') = Ixc(G1-G) * tmp1(G1,G')
     !
!     do ig1 = 1, ngv
!        !
!        do ig2 = dsp + 1, dsp + cnt
!           !
!           idg(1:3) = igv0(1:3,gindx(ig2)) - igv0(1:3,gindx(ig1)) 
!           !
!           ! make matrix Ixc
!           !
!           if(all(gminxc(1:3) <= idg(1:3)) .and. all(idg(1:3) <= gmaxxc(1:3))) then
!              vec(ig2) = Ixc(idg(1), idg(2), idg(3))
!           else
!              vec(ig2) = zero
!           end if
!           !
!        end do ! ig2
        !
!     do imf = 0, nmf
        call zgemm("N","N", ngv, ngv, ngv,                   &
        &            one, vec(1:ngv, 1:ngv            ), ngv,&
        &                tmp1(       1:ngv, 1:ngv     ), ngv,&
        &          zero,  Wsf(1:ngv,        1:ngv, imf), ngv)
!        call zgemm("N","N", ngv, ngv, cnt,                           &
!        &            one, vec(1:ngv, dsp+1:dsp+cnt       ), ngv,&
!        &                tmp1(       dsp+1:dsp+cnt, 1:ngv), cnt,&
!        &          zero, tmp2(1:ngv,                1:ngv), ngv)
!        call zgemv("T", ngv, cnt, one, tmp1(dsp+1:dsp+cnt, 1:ngv, imf),  &
!        &          ngv, vec(dsp+1:dsp+cnt), 1, zero,                      &
!        &          Wsf(ig1, 1:ngv, imf), 1)
!     end do
!     end do ! ig1
     !
!     call MPI_allREDUCE(MPI_IN_PLACE, tmp2, ngv * ngv , &
!     &                   MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
!     Wsf(1:ngv, 1:ngv, imf) = tmp2(1:ngv, 1:ngv)
     !
     end do ! imf
     !
     !Wsf(:,:,:) = 3.0d0 * Wsf(:,:,:)
     Wsf(:,:,:) = -3.0d0 * Wsf(:,:,:)
     !
     ! modified to correspond to Essenberger's formalism
     !
!     do imf = 0, nmf
!     do ig2 = 1, ngv
!        !
!        do ig1 = 1, ngv
!           !
!!           idg(1:3) = igv0(1:3,gindx(ig2)) - igv0(1:3,gindx(ig1)) 
!!           !
!!           ! make matrix Ixc
!!           !
!!           if(all(gminxc(1:3) <= idg(1:3)) .and. all(idg(1:3) <= gmaxxc(1:3))) then
!!              vec(ig1, ig2) = Ixc(idg(1), idg(2), idg(3))
!!           else
!!              vec(ig1, ig2) = zero
!!           end if
!           !
!           ! modify Wsf to correspond to Essenberger's formalism
!           !
!            Wsf(ig1, ig2, imf) = Wsf(ig1, ig2, imf) - vec(ig1, ig2)
!           !
!        end do ! ig1
!        !
!     end do ! ig2
!     end do ! imf
     deallocate(Chi_s)
     deallocate(Chi_s_tmp)
     deallocate(tmp1, tmp2, vec)
     !
!     write(fname,*) iq
!     write(fname,'(3a)') "Wsf", trim(adjustl(fname)), ".dat"
!     open(1000, file = trim(fname))
!     !
!     counter=0
!     do ig = 1, ngv
!       g2 = dble(dot_product(igv0(1:3,gindx(ig)), igv0(1:3, gindx(ig)))) 
!       if(g2 <= 1.0d-8 ) then  ! valid only when G = 0
!         counter = counter + 1
!         if( counter > 1) stop "Something is wrong ! "
!           write(1000, *) 0d0, real(Wsf(ig, ig, 0)), aimag(Wsf(ig, ig, 0))
!         do imf = 1, nmf
!           write(1000, *) mf(imf), real(Wsf(ig, ig, imf)), aimag(Wsf(ig, ig, imf))
!         end do
!       end if
!     end do
!     !
!     close(1000)
  end if ! laddxc
  !
!  call MPI_allREDUCE(MPI_IN_PLACE, Wsf, (nmf + 1) * ngv * ngv , &
!  &                  MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
!
end subroutine make_Wsf
!
!  Computes the kernel Ksf as following:
!        Ksf = rho * Wsf * rho
!
subroutine make_Ksf()
  !
  use mpi, only : MPI_STATUS_SIZE, MPI_DOUBLE_COMPLEX, MPI_IN_PLACE, MPI_COMM_WORLD, &
  &               MPI_DOUBLE_PRECISION, MPI_SUM
  use rpa_vals, only : nmf, nf, nk, nb, wfc1, wfc2, Wsf, Ksf, nftot, pi, gq2, cellV, my_rank, petot, ngv, gindx, &
  &                    FFTW_FORWARD, FFTW_ESTIMATE, nkpe
  !use rpa_routines, only : zgemm, zdotc, cnt_and_dsp, cnt_and_dsp_full
  !
  integer :: k, cnt, dsp, ig, ik, ib, jb, imf, ierr, org, ipe, dst, src, &
  &          kcnt(0:petot - 1), kdsp(0:petot - 1), status(MPI_STATUS_SIZE)
  integer(8) :: plan
  complex(8) :: rhin(nftot), rhout(nftot), rho1(ngv,nb), rho3(ngv,nb), Ksf0, &
  &             one = cmplx(1d0, 0d0), zero = cmplx(0d0, 0d0)
  complex(8),allocatable :: rho2(:,:)
  complex(8),allocatable :: Wsf_tmp(:,:,:)
  !
  call cnt_and_dsp(ngv, cnt, dsp)
  call cnt_and_dsp_full(nk, kcnt, kdsp)
  !
  allocate(Wsf_tmp(1:ngv,dsp+1:dsp+cnt,0:nmf))
  allocate(Ksf(0:nmf,nb,nb,nk))
  Wsf_tmp(1:ngv, dsp+1:dsp+cnt, 0:nmf) = Wsf(1:ngv, dsp+1:dsp+cnt, 0:nmf)
  Ksf(0:nmf,1:nb,1:nb,1:nk) = 0d0
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
     !$OMP & SHARED(my_rank,nk,nb,nftot,nmf,wfc1,wfc2,plan,cnt,dsp,zero,one,Wsf_tmp,Ksf,cellV, &
     !$OMP &        ngv,gindx,kcnt,kdsp,org) &
     !$OMP & PRIVATE(ik,ib,jb,imf,rhin,rhout,rho1,rho2,rho3,Ksf0)
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
              &          one, Wsf_tmp(1:ngv, dsp + 1:dsp + cnt,  imf), ngv, &
              &                  rho2(       dsp + 1:dsp + cnt, 1:nb), cnt, &
              &             zero,rho3(1:ngv,                    1:nb), ngv  )
  !            call zgemm("N", "N", ngv, nb, ngv, &
  !            &          one, Wsf(1:ngv, 1:ngv,  imf), ngv, &
  !            &               rho1(      1:ngv, 1:nb), ngv, &
  !            &         zero, rho3(1:ngv,       1:nb), ngv  )
              !
              do jb = 1, nb
                 !
                 Ksf0 = zdotc(ngv, rho1(1:ngv,jb), 1, &
                 &                 rho3(1:ngv,jb), 1)
                 !
                 Ksf(imf,jb,ib,kdsp(org) + ik) = dble(Ksf0) / cellV
                 !
              end do ! jb = 1, nb
              !
           end do ! imf = 0, nmf
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
  call MPI_allREDUCE(MPI_IN_PLACE, Ksf, (nmf + 1) * nb * nb * nk, &
  &                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  !
  deallocate(Wsf)
  deallocate(Wsf_tmp)
  !
end subroutine make_Ksf
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!             end of modified region. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Invertion of matrices
!
subroutine invert()
  !
  use mpi, only : MPI_IN_PLACE, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, MPI_COMM_WORLD, &
  &               MPI_DOUBLE_COMPLEX, MPI_SUM
  use rpa_vals, only : wscr, ngv, nmf, pi, cellV
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
  use rpa_vals, only : nmf, nf, nk, nb, wfc1, wfc2, Kel, wscr, nftot, pi, gq2, cellV, my_rank, petot, ngv, gindx, &
  &                    FFTW_FORWARD, FFTW_ESTIMATE, nkpe
  !use rpa_routines, only : zgemm, zdotc, cnt_and_dsp, cnt_and_dsp_full
  !
  integer :: cnt, dsp, ik, ig, ib, jb, imf, ierr, org, ipe, dst, k, src, &
  &          kcnt(0:petot - 1), kdsp(0:petot - 1), status(MPI_STATUS_SIZE)
  integer(8) :: plan
  complex(8) :: rhin(nftot), rhout(nftot), rho1(ngv,nb), rho3(ngv,nb), Kel0, &
  &             one = cmplx(1d0, 0d0), zero = cmplx(0d0, 0d0)
  complex(8),allocatable :: rho2(:,:)
  !
  call cnt_and_dsp(ngv, cnt, dsp)
  call cnt_and_dsp_full(nk, kcnt, kdsp)
  !
  allocate(Kel(0:nmf,nb,nb,nk))
  Kel(0:nmf,1:nb,1:nb,1:nk) = 0d0
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
     !$OMP &        ngv,gindx,kcnt,kdsp,org) &
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
  call MPI_allREDUCE(MPI_IN_PLACE, Kel, (nmf + 1) * nb * nb * nk, &
  &                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  !
  deallocate(wscr)
  !
end subroutine make_Kel
!
! Output to file
!
subroutine write_data(iq)
  !
  use rpa_vals, only : my_rank, nmf, nk, nb, Kel, Ksf, ng, nb, iqv, mf, wmf
  !
  integer,intent(in) :: iq
  !
  integer :: fo = 20, fo2 = 21,  imf
  character(100) :: fname
  !
  if(my_rank == 0) then
     !
     write(fname,*) iq
     write(fname,'(3a)') "vc", trim(adjustl(fname)), ".dat"
     !
     open(fo, file = trim(fname), form = "unformatted")
     !
     write(fo) ng(1:3)
     write(fo) nb
     write(fo) dble(iqv(1:3,iq)) / dble(2 * ng(1:3))
     write(fo) nmf
     !
     do imf = 1, nmf
        write(fo) mf(imf) !, wmf(imf)
     end do
     !
     write(fo) Kel(0:nmf,1:nb,1:nb,1:nk)
     !
     close(fo)
     !
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !      modified by Kentaro Tsutsumi.                   !
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     write(fname,*) iq
     write(fname,'(3a)') "Ksf", trim(adjustl(fname)), ".dat"
     !
     open(fo2, file = trim(fname), form = "unformatted")
     !
     write(fo2) ng(1:3)
     write(fo2) nb
     write(fo2) dble(iqv(1:3,iq)) / dble(2 * ng(1:3))
     write(fo2) nmf
     !
     do imf = 1, nmf
        write(fo2) mf(imf) !, wmf(imf)
     end do
     !
     write(fo2) Ksf(0:nmf,1:nb,1:nb,1:nk)
     !
     close(fo2)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !      end of modified region.                         !
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end if
  !
  deallocate(Kel)
  deallocate(Ksf)
  !
end subroutine write_data
!
end module rpa_routines
!
! Main routine
!
program rpa
  !
  use mpi, only : MPI_COMM_WORLD, MPI_INIT, MPI_COMM_SIZE, MPI_COMM_RANK
  use omp_lib, only : OMP_GET_WTIME, OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
  use rpa_vals, only : my_rank, petot, iqdo, ltetra, laddxc
  use rpa_routines,only : stdin, read_file, irr_bz, get_wfcg, fft_mesh, &
  &                       fft_wfc, read_dmxc, apply_xc, make_scrn, invert, make_kel, &
  &                       alloc_Kel, tetra_type, prepare_q, write_data, &
  &                       read_Ixc, make_Chi_s, make_Wsf, make_Ksf, &
  &                       energy_grid, calc_dosk
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
  if(laddxc == 1) then
     !
     if(my_rank == 0) then
        write(*,*)
        write(*,*) "#####  Read Ixc.dat  #####"
        write(*,*)
     end if
     !
     t1 = OMP_GET_WTIME()
     call read_Ixc()
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
  t1 = OMP_GET_WTIME()
  call prepare_q(iqdo)
  call make_scrn(iqdo)
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,'(a,f15.3,a)') "    Time (make_scrn) : ", t2 - t1
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      modified by Kentaro Tsutsumi.                   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  Compute K_sf  #####"
     write(*,*)
  end if
  !
  t1 = OMP_GET_WTIME()
  call make_Chi_s(iqdo)
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,'(a,f15.3,a)') "    Time (calc_Chi_s) : ", t2 - t1
  !
  t1 = OMP_GET_WTIME()
  call make_Wsf()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,'(a,f15.3,a)') "    Time (  make_Wsf) : ", t2 - t1
  !
  t1 = OMP_GET_WTIME()
  call make_Ksf()
  t2 = OMP_GET_WTIME()
  if(my_rank == 0) write(*,'(a,f15.3,a)') "    Time (  make_Ksf) : ", t2 - t1
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !       end of modified region.                        !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  if(my_rank == 0) then
     write(*,*)
     write(*,*) "#####  Compute K_el  #####"
     write(*,*)
  end if
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

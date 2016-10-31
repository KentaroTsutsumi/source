program bin2txt
  !
  implicit none
  !
  integer,parameter :: NTQ = 68
  integer :: fi = 10, ng(3), nk, nb, nmf, imf, iq
  real(8) :: qv(3)
  real(8),allocatable :: mf(:), wmf(:), Kel(:,:,:,:)
  character(100) :: fname, fname1, fname2
  !
  do iq = 1, NTQ
  !
  write(fname,*) iq
  !
  !  Kel.dat -> Kel.txt
  !
  write(fname1, '(3a)') "Kel", trim(adjustl(fname)), ".dat"
  write(fname2, '(3a)') "Kel", trim(adjustl(fname)), ".txt"
  !
  open(fi, file = trim(fname1), form = 'unformatted')
  !
  read(fi) ng(1:3)
  read(fi) nb
  read(fi) qv(1:3)
  read(fi) nmf
!write(*,*) ng(1:3)
!write(*,*) nb
!write(*,*) qv(1:3)
!write(*,*) nmf
  !
  nk = product(ng(1:3))
  allocate(mf(nmf), wmf(nmf), Kel(0:nmf,nb,nb,nk))
  !
  do imf = 1, nmf
     read(fi) mf(imf), wmf(imf)
  end do
  !
  read(fi) Kel(0:nmf,1:nb,1:nb,1:nk)
  !
  close(fi)
  !
  open(fi, file = trim(fname2))
!  open(fi, file = 'Kel_kwmr1.txt')
  !
  write(fi,*) ng(1:3)
  write(fi,*) nb
  write(fi,*) qv(1:3)
  write(fi,*) nmf
  !
  do imf = 1, nmf
     write(fi,*) mf(imf), wmf(imf)
  end do
  !
  write(fi,*) ""
  !
  write(fi,'(3e25.15)') Kel(0:nmf,1:nb,1:nb,1:nk)
  !
  deallocate(mf,wmf,Kel)
  close(fi)
  !
  !  Ksf.dat -> Ksf.txt
  !
  write(fname1, '(3a)') "Ksf", trim(adjustl(fname)), ".dat"
  write(fname2, '(3a)') "Ksf", trim(adjustl(fname)), ".txt"
  !
  open(fi, file = trim(fname1), form = 'unformatted')
  !
  read(fi) ng(1:3)
  read(fi) nb
  read(fi) qv(1:3)
  read(fi) nmf
!write(*,*) ng(1:3)
!write(*,*) nb
!write(*,*) qv(1:3)
!write(*,*) nmf
  !
  nk = product(ng(1:3))
  allocate(mf(nmf), wmf(nmf), Kel(0:nmf,nb,nb,nk))
  !
  do imf = 1, nmf
     read(fi) mf(imf), wmf(imf)
  end do
  !
  read(fi) Kel(0:nmf,1:nb,1:nb,1:nk)
  !
  close(fi)
  !
  open(fi, file = trim(fname2))
!  open(fi, file = 'Kel_kwmr1.txt')
  !
  write(fi,*) ng(1:3)
  write(fi,*) nb
  write(fi,*) qv(1:3)
  write(fi,*) nmf
  !
  do imf = 1, nmf
     write(fi,*) mf(imf), wmf(imf)
  end do
  !
  write(fi,*) ""
  !
  write(fi,'(3e25.15)') Kel(0:nmf,1:nb,1:nb,1:nk)
  !
  deallocate(mf,wmf,Kel)
  close(fi)
  end do ! iq
  !
!  write(*,*) " "
!  imf = 0
!  write(*,*) imf, 0d0, sum(Kel(imf,1:nb,1:nb,1:nk)) / dble(nb * nb * nk)
!  do imf = 1, nmf
!     write(*,*) imf, mf(imf), sum(Kel(imf,1:nb,1:nb,1:nk)) / dble(nb * nb * nk)
!  end do
  !
end program bin2txt

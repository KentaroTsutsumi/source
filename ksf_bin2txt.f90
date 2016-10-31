program ksf_bin2txt
  !
  implicit none
  !
  integer :: fi = 10, ng(3), nk, nb, nmf, imf, ik, ib, jb
  real(8) :: qv(3)
  real(8),allocatable :: mf(:), Ksf(:,:,:,:)
  integer :: iniq, lasq, iq
  character(100) :: fname, fname1, fname2, fname3, fname4
  !
  write(*,*) "input the start_q and last_q number"
  read(*,*) iniq
  read(*,*) lasq
  !
  do iq = iniq, lasq
    !
    write(fname,*) iq
    write(fname1,*) "vc", trim(adjustl(fname)), ".dat"
    write(fname2,*) "vc", trim(adjustl(fname)), ".txt"
    write(fname3,*) "Ksf", trim(adjustl(fname)), ".dat"
    write(fname4,*) "Ksf", trim(adjustl(fname)), ".txt"
  !  open(fi, file = 'Ksf69.dat', form = 'unformatted')
  !
  !  Kel${iq}.dat => Kel${iq}.txt
  !
    open(fi, file = trim(fname1), form = 'unformatted')
    !
    read(fi) ng(1:3)
    read(fi) nb
    read(fi) qv(1:3)
    read(fi) nmf
  write(*,*) ng(1:3)
  write(*,*) nb
  write(*,*) qv(1:3)
  write(*,*) nmf
    !
    nk = product(ng(1:3))
    allocate(mf(nmf), Ksf(0:nmf,nb,nb,nk))
    !
    do imf = 1, nmf
       read(fi) mf(imf) !, wmf(imf)
    end do
    !
    read(fi) Ksf(0:nmf,1:nb,1:nb,1:nk)
    !
    close(fi)
    !
  !  open(fi, file = 'Ksf69.txt')
    open(fi, file = trim(fname2))
    !
    write(fi,*) ng(1:3)
    write(fi,*) nb
    write(fi,*) qv(1:3)
    write(fi,*) nmf
    !
    do imf = 1, nmf
       write(fi,*) mf(imf) !, wmf(imf)
    end do
    !
    write(fi,*) ""
    !
    write(fi,'(3e25.15)') Ksf(0:nmf,1:nb,1:nb,1:nk)
    !
    close(fi)
    !
    write(*,*) " "
    imf = 0
    write(*,*) imf, 0d0, sum(Ksf(imf,1:nb,1:nb,1:nk)) / dble(nb * nb * nk)
    do imf = 1, nmf
       write(*,*) imf, mf(imf), sum(Ksf(imf,1:nb,1:nb,1:nk)) / dble(nb * nb * nk)
    end do
    deallocate(mf, Ksf)
    !
    !  Ksf${iq}.dat => Ksf${iq}.txt
    !
    open(fi, file = trim(fname3), form = 'unformatted')
    !
    read(fi) ng(1:3)
    read(fi) nb
    read(fi) qv(1:3)
    read(fi) nmf
  write(*,*) ng(1:3)
  write(*,*) nb
  write(*,*) qv(1:3)
  write(*,*) nmf
    !
    nk = product(ng(1:3))
    allocate(mf(nmf), Ksf(0:nmf,nb,nb,nk))
    !
    do imf = 1, nmf
       read(fi) mf(imf) !, wmf(imf)
    end do
    !
    read(fi) Ksf(0:nmf,1:nb,1:nb,1:nk)
    !
    close(fi)
    !
  !  open(fi, file = 'Ksf69.txt')
    open(fi, file = trim(fname4))
    !
    write(fi,*) ng(1:3)
    write(fi,*) nb
    write(fi,*) qv(1:3)
    write(fi,*) nmf
    !
    do imf = 1, nmf
       write(fi,*) mf(imf) !, wmf(imf)
    end do
    !
    write(fi,*) ""
    !
    write(fi,'(3e25.15)') Ksf(0:nmf,1:nb,1:nb,1:nk)
    !
    close(fi)
    !
    write(*,*) " "
    imf = 0
    write(*,*) imf, 0d0, sum(Ksf(imf,1:nb,1:nb,1:nk)) / dble(nb * nb * nk)
    do imf = 1, nmf
       write(*,*) imf, mf(imf), sum(Ksf(imf,1:nb,1:nb,1:nk)) / dble(nb * nb * nk)
    end do
    deallocate(mf, Ksf)
  end do
end program ksf_bin2txt

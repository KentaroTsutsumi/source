module verify_vals
  !
  implicit none
  !
  real(8),parameter :: pi = acos(-1d0)
  !
  integer,save :: &
  & ncf,           & ! # of points for spline
  & nmf,           & ! # of Matsubara freq.
  & ng(3),         & ! k point grid
  & nk,            & ! ng(1)*ng(2)*ng(3). Total # of k points
  & nb               ! # of bands
  !
  real(8),allocatable,save :: &
  & mf(:),           & ! (nmf) Matsubara frequencies
  & Kelc(:,:,:,:),   & ! (ncf,nb,nb,nk). Coulomb Kernel
  & Kel(:,:,:,:)       ! (0:nmf,nb,nb,nk). Coulomb Kernel
  !
end module verify_vals
!
module verify_routines
  !
  implicit none
  !
contains
  !
 ! Read vc.dat
!
subroutine read_vc()
  !
  use verify_vals, only : nk, nb, nmf, ncf, mf, Kel, Kelc, mf, ng
  !
  integer :: fi = 10, imf, ng0(3), nb0
  real(8) :: qv(3), qv0(3)
  !
  open(fi, file = "vc.dat", form= "unformatted")
  !
  read(fi) ng(1:3)
  nk = product(ng(1:3))
  !
  read(fi) nb
  read(fi) qv(1:3)
  read(fi) nmf
  allocate(Kel(0:nmf,nb,nb,nk), mf(nmf))
  !
  do imf = 1, nmf
     read(fi) mf(imf)
  end do
  !
  read(fi) Kel(0:nmf,1:nb,1:nb,1:nk)
  !
  close(fi)
  !
  open(fi, file = "vel.dat", form= "unformatted")
  !
  read(fi) ng0(1:3)
  if(.not. all(ng(1:3) == ng0(1:3))) stop "ng0 /= ng"
  !
  read(fi) nb0
  if(nb0 /= nb) stop "nb0 /= nb"
  !
  read(fi) qv0(1:3)
  qv(1:3) = qv(1:3) - qv0(1:3)
  if(dot_product(qv(1:3), qv(1:3)) > 1d-10) stop "qv0 /= qv"
  !
  read(fi) ncf
  allocate(Kelc(ncf,nb,nb,nk))
  !
  read(fi) Kelc(1:ncf,1:nb,1:nb,1:nk)
  !
  close(fi)
  !
  write(*,*) ""
  write(*,*) "     nk : ", nk
  write(*,*) "ng(1:3) : ", ng(1:3)
  write(*,*) "     nb : ", nb
  write(*,*) "    nmf : ", nmf
  write(*,*) "    ncf : ", ncf
  write(*,*) ""
  !
end subroutine read_vc
!
!
!
subroutine verify_vc()
  !
  use verify_vals, only : pi, nb, nk, nmf, ncf, mf, Kel, Kelc
  !
  integer :: ik, ib, jb, imf, icf
  real :: Kel0, aveer, maxer, mf0, err, x0
  !
  x0 = cos(pi / dble(2 * ncf))
  !
  aveer = 0d0
  maxer = 0d0
  !
  do ik = 1, nk
     do ib = 1, nb
        do jb = 1, nb
           !
           do imf = 1, nmf
              !
              mf0 = (mf(imf) - 1d0) / (mf(imf) + 1d0) * x0
              !
              Kel0 = 0
              do icf = 1, ncf
                 Kel0 = Kel0 + Kelc(icf,jb,ib,ik) * cos(dble(icf - 1) * acos(mf0))
              end do
              !
              err = abs((Kel0 - Kel(imf,jb,ib,ik)) / Kel(imf,jb,ib,ik))
              !
              aveer = aveer + err
              maxer = max(maxer, err)
              !
              Kel(imf,jb,ib,ik) = Kel0
              !
           end do ! imf
           !
        end do ! jb
     end do ! ib
  end do ! ik
  !
  aveer = aveer / dble(nmf * nb * nb * nk)
  !
  write(*,*) "maxer : ", 1d2 * maxer, "%"
  write(*,*) "aveer : ", 1d2 * aveer, "%"
  !
end subroutine verify_vc
!
!
!
subroutine write_vc()
  !
  use verify_vals, only : Kel, nb, nmf, nk
  !
  integer :: fo = 20
  !
  open(fo, file = "vc.new")
  !
  write(fo,'(3e25.15)') Kel(0:nmf,1:nb,1:nb,1:nk)
  !
  close(fo)
  !
end subroutine write_vc
!
end module verify_routines
!
program verify
  !
  use verify_routines, only : read_vc, verify_vc, write_vc
  !
  implicit none
  !
  call read_vc()
  call verify_vc()
  call write_vc()
  !
end program verify

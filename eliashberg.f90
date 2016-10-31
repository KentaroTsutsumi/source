module eliashberg_vals
  !
  implicit none
  !
  integer,parameter :: &
  & itrmax = 500 ! Max # of iteration
  !
  real(8),parameter :: &
  & alpha_mix = 0.2d0, &
  & pi = acos(-1d0), &
  & er0 = 1d-15   ! Convergence threshold
  !
  integer,save :: &
  & nmf  ! # of Matsubara frequency
  !
  real(8),save :: &
  & temp, & ! temperature [Ry]
  & omln, & ! Omega_ln [Ry]
  & lambda, & ! The Flochlich parameter
  & mu ! Coulomb potential
  !
  real(8),allocatable,save :: &
  & mf(:), & ! (nmf) Matsubara frequencies [Ry]
  & Z(:), &  ! (nmf) Renormalization
  & delta(:) ! (nmf) SC gap  
  !
end module eliashberg_vals
!
module eliashberg_routines
  !
  implicit none
  !
contains
  !
 ! Standard input
!
subroutine read_stdin()
  !
  use eliashberg_vals, only : nmf, temp, lambda, mu, omln
  !
  integer :: ierr
  namelist /input/ nmf, temp, lambda, mu, omln
  !
  write(*,*) ''
  write(*,*) '###############  Standard Input  ###############'
  write(*,*) ''
  !
  read(*,input,err=100)
  write(*,*) "              Temparature[K] : ", temp
  temp = temp / 1.578865419d5
  write(*,*) "            Temperature [Ry] : ", temp
  write(*,*) "  # of Matsubara frequencies : ", nmf     
  write(*,*) "                      lambda : ", lambda
  write(*,*) "               Omega_ln [Ry] : ", omln
  write(*,*) "           Coulomb potential : ", mu
  !
  return
  !
100 write(*,*) "Stop in read_stdin. reading namelist file"
  write(*,*) "              Temparature[K] : ", temp
  write(*,*) "  # of Matsubara frequencies : ", nmf     
  write(*,*) "                      Lambda : ", lambda
  write(*,*) "               Omega_ln [Ry] : ", omln
  write(*,*) "           Coulomb potential : ", mu
  stop
  !
end subroutine read_stdin
!
! Set Matsubara-frequency grid
!
subroutine set_matsubara_grid()
  !
  use eliashberg_vals, only : nmf, mf, Z, delta, temp, pi, lambda, omln
  !
  integer :: imf, is, fi = 11, nmf0
  real(8),allocatable :: mf0(:), Z0(:), delta0(:)
  !
  allocate(mf(nmf), Z(nmf), delta(nmf))
  !
  do imf = 1, nmf
     mf(imf) = (2d0 * dble(imf) - 1) * pi * temp
     Z(imf) = 1 + lambda / (1d0 + mf(imf) / omln)
     delta(imf) = lambda / (1d0 + mf(imf) / omln)
  end do
  !
  open(fi, file = 'delta.dat',status="old", action = 'read',iostat = is)
  !
  if(is == 0) then
     !
     write(*,*) "  Read delta.dat"
     !
     read(fi,*) nmf0
     !
     allocate(mf0(nmf0), Z0(nmf0), delta0(nmf0))
     !
     do imf = 1, nmf0
        read(fi,*) mf0(imf), Z0(imf), delta0(imf)
     end do
     !
     close(fi)
     !
     do imf = 1, nmf
        Z(imf) = Z0(minloc(abs(mf(imf) - mf0(1:nmf0)), 1))
        delta(imf) = delta0(minloc(abs(mf(imf) - mf0(1:nmf0)), 1))
     end do
     !
     deallocate(mf0,Z0,delta0)
     !
  end if
  !
end subroutine set_matsubara_grid
!
! Solve gap equation
!
subroutine gapeq()
  !
  use omp_lib
  use eliashberg_vals, only : alpha_mix, er0, itrmax, nmf, Z, delta
  !
  implicit none
  !
  integer :: itr, jtr
  real(8) :: res, dd(2 * nmf), rhs(2 * nmf), rhs0(2 * nmf), drhs(2 * nmf), &
  &          jacob1(2 * nmf,itrmax), jacob2(2 * nmf,itrmax), t2, t1
  !
  write(*,*) '  Convergense threshold : ', er0
  write(*,*) '           Max itration : ',itrmax
  write(*,*) '   Initial mixing alpha : ', alpha_mix
  write(*,*) ""
  !
  t1 = OMP_GET_WTIME()
  !
  itr = 0
  write(*,*) "  Iteration ", itr
  !
  call make_rhs(rhs)
  res = dot_product(rhs, rhs)
  res = sqrt(res) / dble(nmf * 2)
Z(1:nmf) = rhs(1:nmf)
delta(1:nmf) = rhs(nmf+1:nmf * 2)
return
  !
  ! $$$$$  Information of conversience  $$$$$
  !
  t2 = OMP_GET_WTIME() - t1
  t1 = OMP_GET_WTIME()
  !
  write(*,*) "      Residual[Ry] : ", res
  write(*,*) "     Delta_0 [meV] : ", delta(1) * 13.60569228d3
  write(*,*) "           Z_0 - 1 : ", Z(1) - 1
  write(*,*) "    Rap time [sec] : ", t2
  write(*,*) ""
  !
  ! $$$$$  End information of conversience  $$$$$
  !
  if(res < er0) goto 5
  !
  dd(1:2 * nmf) = - alpha_mix * rhs(1:2 * nmf)
  !
  do itr = 1, itrmax
     !
     write(*,*) "  Iteration ", itr
     !
     Z(    1:nmf) = Z(    1:nmf) + dd(      1 :     nmf)
     delta(1:nmf) = delta(1:nmf) + dd(nmf + 1 : 2 * nmf)
     !
     rhs0(1:2 * nmf) = rhs(1:2 * nmf)
     call make_rhs(rhs)
     res = dot_product(rhs, rhs)
     res = sqrt(res) / dble(2 * nmf)
     !
     ! $$$$$  Information of conversience  $$$$$
     !
     t2 = OMP_GET_WTIME() - t1
     t1 = OMP_GET_WTIME()
     !
     write(*,*) "      Residual[Ry] : ", res
     write(*,*) "     Delta_0 [meV] : ", delta(1) * 13.60569228d3
     write(*,*) "           Z_0 - 1 : ", Z(1) - 1
     write(*,*) "    Rap time [sec] : ", t2
     write(*,*) ""
     !
     ! $$$$$  End information of conversience  $$$$$
     !
     if(res < er0) then
        !       
        delta(1:nmf) = delta(1:nmf) * sign(1d0, delta(1))
        !
        goto 5
        !
     end if
     !
     ! Update Jacobian with drhs
     !
     drhs(1:2 * nmf) = rhs(1:2 * nmf) - rhs0(1:2 * nmf)
     !
     jacob1(1:2 * nmf, itr) = - alpha_mix * drhs(1:2 * nmf)
     do jtr = 1, itr - 1
        jacob1(1:2 * nmf,itr) = jacob1(1:2 * nmf,itr) - jacob1(1:2 * nmf,jtr) &
        &          * dot_product(jacob2(1:2 * nmf,jtr), drhs(1:2 * nmf))
     end do
     jacob1(1:2 * nmf,itr) = dd(1:2 * nmf) + jacob1(1:2 * nmf,itr)
     jacob2(1:2 * nmf,itr) = drhs(1:2 * nmf) / dot_product(drhs(1:2 * nmf), drhs(1:2 * nmf))
     !
     ! Compute dd with new Jacobian & rhs
     !
     dd(1:2 * nmf) = - alpha_mix * rhs(1:2 * nmf)
     do jtr = 1, itr
        dd(1:2 * nmf) = dd(1:2 * nmf) - jacob1(1:2 * nmf,jtr) &
        &        * dot_product(jacob2(1:2 * nmf,jtr), rhs(1:2 * nmf))
     end do
     !
  end do ! itr
  !
  write(*,*) ""
  write(*,*) '  Not converged! res = ',res
  return
  !
5 continue
  !
  write(*,*) ""
  write(*,*) '  Converged! iter = ',itr, ",     Delta_0 [meV] : ", delta(1) * 13.60569228d3
  !
end subroutine gapeq
!
! Make RHS vector for linear probrem.
!
subroutine make_rhs(veco)
  !
  use eliashberg_vals, only : pi, temp, nmf, mf, Z, delta, omln, lambda, mu
  implicit none
  !
  real(8),intent(out) :: veco(2 * nmf)
  !
  integer :: imf, jmf
  real(8) :: theta(nmf), Zrhs(nmf), Drhs(nmf), lambda2
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(temp,nmf,mf,Z,delta,omln,lambda,mu,theta,Zrhs,Drhs) &
  !$OMP & PRIVATE(imf,jmf,lambda2)
  !
  !$OMP DO
  do imf =1, nmf
     theta(imf) = 1d0 / sqrt(mf(imf)**2 + delta(imf)**2) 
  end do
  !$OMP END DO
  !
  !$OMP DO
  do imf = 1, nmf
     !
     Zrhs(imf) = 0d0
     Drhs(imf) = 0d0
     !
     do jmf = 1, nmf
        !
        lambda2 = lambda * omln**2 / ( ((mf(imf) - mf(jmf))**2 + omln**2) &
        &                            * ((mf(imf) + mf(jmf))**2 + omln**2) )
        !
!        Zrhs(imf) = Zrhs(imf) + mf(jmf)**2 * theta(jmf) * lambda2
!        Drhs(imf) = Drhs(imf) + delta(jmf) * theta(jmf) &
!        &         * ((mf(imf)**2 + mf(jmf)**2 + omln**2) * lambda2 - mu)
Zrhs(imf) = Zrhs(imf) + delta(jmf) * theta(jmf) &
&         * ((mf(imf)**2 + mf(jmf)**2 + omln**2) * lambda2)
Drhs(imf) = Drhs(imf) + delta(jmf) * theta(jmf) &
&         * (mu)
        ! 
     end do
     !
!     Zrhs(imf) = Z(imf) - 1d0 - 4d0 * pi * temp * Zrhs(imf)
!     Drhs(imf) = delta(imf) - 2d0 * pi * temp * Drhs(imf) / Z(imf)
Zrhs(imf) = 2d0 * pi * temp * Zrhs(imf) / Z(imf)
Drhs(imf) = 2d0 * pi * temp * Drhs(imf) / Z(imf)
     !
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  !
  veco(      1:    nmf) = Zrhs(1:nmf)
  veco(nmf + 1:2 * nmf) = Drhs(1:nmf)
  !
end subroutine make_rhs
!
! Output delta
!
subroutine out_delta()
  !
  use eliashberg_vals, only : nmf, mf, delta, Z
  !
  integer :: imf, fo = 20
  !
  open(fo, file = "delta.dat")
  !
  write(fo,*) nmf
  write(fo,*) ""
  !
  do imf = 1, nmf
     write(fo,*) mf(imf), Z(imf), delta(imf)
  end do
  !
  close(fo)
  !
end subroutine out_delta
!
end module eliashberg_routines
!
program eliashberg
  !
  use eliashberg_routines, only : read_stdin, set_matsubara_grid, gapeq, out_delta
  !
  implicit none
  !
  call read_stdin()
  !
  call set_matsubara_grid()
  !
  call gapeq()
  !
  call out_delta()
  !
end program eliashberg

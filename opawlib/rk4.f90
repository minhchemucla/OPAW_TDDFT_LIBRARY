subroutine rk4_prop_opaw(n,nocc,nstates,dt,p,ham)
  use opaw_mod, only : dv
  use opaw_ham_mod
  use mpi_lib_ours
  implicit none
  !integer :: is, jb, n, nb
  integer :: n,nocc,nstates
  integer :: is,ie,ns,i,st
  real*8 :: dt
  real*8 ::  nrm0(nstates), nrm1(nstates)
  !complex*16, dimension(n,ns) :: k1,k2,k4,y, p, pp, p2,pp2
  complex*16, dimension(n,nstates) :: p
  complex*16, allocatable, dimension(:,:) :: k1,k2,k4,y,pp
  complex*16, parameter :: ci = (0d0,1d0)
  type(opaw_ham_obj) :: ham
  !real*8 :: wtime, wtime2

  do i=1,nstates
    nrm0(i) = sqrt(sum(abs(p(:,i))**2d0)*dv)
  enddo

  call calc_is_ie(is,ie,nstates)
  ns = ie-is+1
  allocate(k1(n,ns),stat=st);if(st/=0) stop 'k1 rk4'
  allocate(k2,k4,y,pp,mold=k1,stat=st);if(st/=0) stop 'k2,k4,y rk4'

  call scatterv_c16(p,pp,size(p),size(pp),0)

  call rk4_step(pp,k1,y,dt/2d0,.false.)
  call rk4_step(y,k2,y,dt/2d0,.false.)
  call rk4_step(y,k4,y,dt,.false.) !k4=k3
  k2 = k2 + k4 !k2+k3 
  call rk4_step(y,k4,y,dt,.true.)
  pp = pp + dt/6d0*(k1+2d0*k2+k4)

  call gatherv_c16(pp,p,size(pp),size(p),0)
  !do i=1,nocc
  !  write(*,*) 'i, p rk4', i, sum(abs(p(:,i)))
  !enddo
  do i=1,nstates
    nrm1(i) = sqrt(sum(abs(p(:,i))**2d0)*dv)
    p(:,i) = p(:,i)*nrm0(i)/nrm1(i)
  enddo
  !do i=1,nocc
  !  write(*,*) 'i, norm p rk4', i, nrm0(i), nrm1(i), sum(abs(p(:,i)))
  !enddo

  contains 
    subroutine rk4_step(pin,k,y,dt,last)
      implicit none
      real*8     :: dt
      logical    :: last
      complex*16 :: k(n,ns), y(n,ns), pin(n,ns)

      do is=1,ns
        call shs(ham,pin(:,is),k(:,is))
        k(:,is) = -ci*k(:,is)
        if(.not. last) then
          y(:,is) = pp(:,is) + dt * k(:,is) !y1
        endif
      enddo

      if(.not.last) then
        call gatherv_c16(y,p,size(y),size(p),0)
        call opaw_make_hamiltonian(n,nocc,nstates,p,ham)
      endif
    end subroutine
end subroutine

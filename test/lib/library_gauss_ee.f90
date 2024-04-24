subroutine gauss_2b_damped(RA, na, za, qa, &
                           RB, nb, zb, qb, &
                           RC, nc, zc, qc, &
                           RD, nd, zd, qd , &
                           s, g )

  ! ( p_A(r1) p_B(r1) | erfc(s*r12)/r12 | p_A(r2) p_B(r2) ), p_A = sum_(ia=1,..na)qa(ia)* exp(-za(ia)*(r-RA)**2) 

  implicit none
  integer na, nb, nc, nd
  real*8 RA(3), za(na), qa(na)
  real*8 RB(3), zb(nb), qb(nb)
  real*8 RC(3), zc(nc), qc(nc)
  real*8 RD(3), zd(nd), qd(nd)
  real*8 s, g
  
  integer ia,ib,ic,id
  real*8 g1

  g = 0d0
  do ia=1,na
     do ib=1,nb
        do ic=1,nc
           do id=1,nd
              call gauss_1s(RA, RB, RC, RD, za(ia), zb(ib), zc(ic), zd(id), s, g1)
              g = g + qa(ia)*qb(ib)*qc(ic)*qd(id)* g1
           enddo
        enddo
     enddo
  enddo
end subroutine gauss_2b_damped


subroutine gauss_1s(RA, RB, RC, RD, a, b, c, d, s, govrlp) ! (AB| erfc(s |r12|)/|r12| |CD) where A,B,C,D un-normalized
  implicit none
  !integer, parameter :: QR = selected_real_kind(16)
  !real (kind=QR) :: gauss_1d

  real*8 :: ra(3), rb(3), rc(3), rd(3)
  real*8 :: rp(3), rq(3)
  real*8 :: a, b, c, d
  real*8 :: pi, p,  q, s, M, u, w, ws
  real*8 :: govrlp, g41
  real*8, parameter :: sml = 1d-10

  pi  = acos(-1d0)

  ! using Szabo and Ostlund, Eq. A.41, page 431, and, based on Eq. A.40, see also Word doc: Gaussian Integrlas.docx
  ! 
  
  p = a + b
  q = c + d

  rp = (a*ra+b*rb)/(a+b)
  rq = (c*rc+d*rd)/(c+d)

  u = max(sqrt(sum((rp-rq)**2)),sml)  ! x in the notes= |Rp-Rq|

  M = exp (-a*b/p * sum((ra-rb)**2) - c*d/q* sum((rc-rd)**2) )

  w  = (p+q)/(4d0*p*q)
  ws =  w + 1d0/(4d0*max(1d-20,s**2))

  govrlp = M * (pi**2d0/(p*q))**(3d0/2d0)* 1d0/u * (erf(u/(2d0*sqrt(w))) - erf(u/(2d0*sqrt(ws))) )

  !g41 = 2d0 * pi**2.5d0 /( (a+b)*(c+d)*sqrt(a+b+c+d) ) * M * &
  !     F0((a+b)*(c+d)/(a+b+c+d) * u**2 )
  !write(6,*)' govrlp ',govrlp
  !write(6,*)' g41 ',g41

  ! abandoned normalized gaussians.  conversion (not done) to normalized gaussians
  !govrlp = govrlp * (2d0/pi)**3d0 * (a*b*c*d)**0.75d0
contains
  real*8 function F0(t)
    implicit none
    real*8 t
    F0 = 0.5d0 * sqrt(pi/t)*erf(sqrt(t))
  end function F0
end subroutine gauss_1s


!call gauss_drvr
!end

subroutine gauss_drvr
  implicit none
  real*8 RA(3), RB(3), RC(3), RD(3)
  real*8, parameter :: a=0.73d0, b=0.5d0, c=0.35d0, d=1.09d0, fall=10,s=1d0/fall ! random values
  real*8 :: govrlp, gauss_1d
  RA = (/ -0.3d0, 1.5d0, 1.9d0 /)
  RB = (/  1.6d0,-0.7d0,-2.1d0 /)
  RC = (/ -1.3d0,-1.0d0, 1.3d0 /)
  RD = (/ -1.5d0, 0.6d0,-1.2d0 /)
  
  call gauss_1s(RA, RB, RC, RD, a, b, c, d, s, govrlp) 
  
  write(6,*)' overlap ',govrlp
end subroutine gauss_drvr
  

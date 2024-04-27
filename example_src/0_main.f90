module main_mod
  implicit none
  integer :: nx, ny, nz, nn
  integer :: funct
  integer :: nstates, nocc, nsp
  integer :: scale_vh
  integer :: funct_x, funct_c 
  integer :: nt
  integer :: ipol

  real*8  :: dx, dy, dz, dv
  real*8  :: ekcut
  real*8  :: rnel
  real*8  :: dt
  real*8 :: sm

  logical :: periodic
  logical :: flg_bin

  complex*16, allocatable :: wfs(:,:)
  real*8, allocatable :: eigs(:), vloc_tot(:)
  real*8, allocatable :: dens(:), vks(:), vxc(:)
  real*8, allocatable  :: vk(:)
end module

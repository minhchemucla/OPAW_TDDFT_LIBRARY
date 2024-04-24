module opaw_mod
!   use main_mod
!=======================================================================================
!     The following variables are things that I expect to be imported from
!     some main_mod in the code that this OPAW library is being integrated into.
!     I put in comments what the shape of the array should be
!=======================================================================================

! 100% need to import
  use main_mod, only : nx,ny,nz,nn ! nn=nx*ny*nz, number of grid points !integer
  use main_mod, only : dx,dy,dz,dv ! grid volume elements !real*8
  use main_mod, only : ekcut       ! kinetic energy cutoff. If ekcut <= 0d0 then code will use  !real*8
  use main_mod, only : scale_vh    ! use Martyna-Tucker (scale_vh=2) or not (scale_vh=1) !integer
  use main_mod, only : nstates     ! total number of states !integer
  use main_mod, only : periodic    !periodic or not !logical
  use main_mod, only : funct       !0=lda, 1=pbe xc functional to be used !integer
  use main_mod, only : funct_c     !correlation functional for LIBXC !integer
  use main_mod, only : funct_x     !exchange functional for LIBXC !integer
                                   !   abinit subroutines to calculate ekcut.
! can be moved to be internal variables if needed
  use main_mod, only : rnel        ! # of electrons !real*8
  use main_mod, only : nocc        ! number of occupied states !integer
  use main_mod, only : vloc_tot  !shape (nn) local part of the pseudopotential !real*8
  use main_mod, only : vk        !shape (nn) related to applying the hartree potential !real*8
  implicit none


!=======================================================================================
!       The following variables can be placed in this module or in main_mod
!=======================================================================================

!=======================================================================================
!                             particulars of OPAW code
!=======================================================================================
!  real*8  :: kkx,kky,kkz       !only gamma point use so these vars aren't used 
   integer :: xmax,ymax,zmax
   integer :: siz_dijall
   real*8, allocatable :: kpt(:,:)
   real*8, allocatable :: nhat(:,:,:)
   real*8, allocatable :: ncoret(:)
   real*8, allocatable :: dijall(:)
   logical :: ekread
   real*8, allocatable  :: ek3d(:,:,:)

!  parameters
   real*8  :: p_fg=0.15d0   !parameter for the fine-grid (Ono-Hirose) feel free to change
                            !nfovnr=n_fine/n_rough = max(dx/p_fg,dy/p_fg,dz/p_fg)
   integer :: nk_loc=1      !do not change this
   real*8  :: rpad=1.0d0
   real*8  :: rpad_r=2d0
   real*8  :: tollsij=0.01
   real*8  :: ek_factor=1d0
!
end module

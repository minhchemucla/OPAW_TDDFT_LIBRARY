!subroutine get_pot_opaw(nn,vks,vxc,vh)
subroutine get_pot_opaw(ham)
    use opaw_mod, only :  nn,nx,ny,nz,vloc_tot,&
        nhat,ncoret,scale_vh,funct,periodic, funct_x, funct_c
    use opaw_ham_mod
    use mpi_lib_ours
    implicit none
    type(opaw_ham_obj) :: ham


      if(rank==0) then
        call get_vxc
        call vh_sub_opaw(ham%dens+ham%nhat,ham%vh,scale_vh) !_opaw label just to avoid conflict
        if(periodic) then
         ham%vh=ham%vh-sum(ham%vh)/size(ham%vh)
        endif
        ham%vks=ham%vh+ham%vxc+vloc_tot
!        write(6,*) 'get_pot dens', sum(ham%dens), sum(ham%nhat)
!        write(6,*) 'get_pot vxc', sum(ham%vxc**3d0), sum(ham%vxc)
!        write(6,*) 'get_pot vh', sum(ham%vh**3d0), sum(ham%vh)
!        write(6,*) 'get_pot vloc_tot', sum(vloc_tot**3d0), sum(vloc_tot)
!        write(6,*) 'get_pot vks', sum(ham%vks**3d0), sum(ham%vks)
      endif
      call bcast_r8(ham%vks,size(ham%vks),0)
      call bcast_r8(ham%vxc,size(ham%vxc),0)
      call bcast_r8(ham%vh,size(ham%vh),0)
contains
    subroutine get_vxc
        implicit none

        if(funct==0) then
            funct_x = 1
            funct_c = 12
        else if(funct==1) then
            funct_x = 101
            funct_c = 130
        else
            write(*,*) 'funct should be 0 or 1'
            stop
        endif
        call vxc_libxc(ham%dens+ncoret,ham%vxc,nn,1)
    end subroutine get_vxc
end subroutine get_pot_opaw

subroutine get_pot_static_xc_opaw(ham)
    use opaw_mod, only :  nn,nx,ny,nz,vloc_tot,&
        nhat,ncoret,scale_vh,funct,periodic, funct_x, funct_c
    use opaw_ham_mod
    use mpi_lib_ours
    implicit none
    type(opaw_ham_obj) :: ham


    if(rank==0) then
      call vh_sub_opaw(ham%dens+ham%nhat,ham%vh,scale_vh) !_opaw label just to avoid conflict
      if(periodic) then
       ham%vh=ham%vh-sum(ham%vh)/size(ham%vh)
      endif
      ham%vks=ham%vh+ham%vxc+vloc_tot
    endif
    call bcast_r8(ham%vks,size(ham%vks),0)
    call bcast_r8(ham%vh,size(ham%vh),0)
end subroutine get_pot_static_xc_opaw

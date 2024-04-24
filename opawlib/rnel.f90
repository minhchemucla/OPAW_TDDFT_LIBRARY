subroutine get_rnel_opaw
      use opaw_mod, only : rnel, nocc!, is_opaw_nocc, ie_opaw_nocc
      use atom_mod
      use mpi_lib_ours
      implicit none
      integer :: ia,it,nds

      nds = max(nodes,1)

      if(rank==0) then
          rnel=0d0
          do ia=1,natom
            it=atom_map(ia)
            rnel=rnel+pawinfo(it)%val
          enddo
          write(*,*) 'total valence charge is',rnel
          nocc=int(rnel/2d0+1d-8)
          write(*,*) 'number of occupied states is',nocc
      endif

      call bcast_scalar_i(nocc)
      call bcast_scalar_r8(rnel)
      call sync_mpi
      
      !call calc_is_ie(is_opaw_nocc,ie_opaw_nocc,nocc)
end subroutine get_rnel_opaw


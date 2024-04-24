module rand_seed_mod
  implicit none; save
  integer, parameter :: nlseed=4,ndiffseeds=6
  integer            :: line_seed=1
  integer(kind=8)    :: seed_array(nlseed, ndiffseeds)
  integer(kind=8)    :: seed_read( nlseed, ndiffseeds)
  integer(kind=8)    :: seed_read_wcounter( nlseed, ndiffseeds)
  integer            :: rank_dependent(    ndiffseeds) 
  logical            :: first = .true. 
  integer            :: mdiffseeds,cnt
  integer(kind=8)    :: j1, j2, j3, j4, j5, cnt8, j8, big8
contains
  subroutine read_initiate_seeds
    use mpi_lib_ours, only : rank,nodes
    implicit none
    integer iseed,i

    if(.not.first) return
    first = .false.

    call read_counter(cnt)
    
    if(rank==0) then
       write(6,*) 
       write(6,*)' Give in random.inp between 1 and ',ndiffseeds,' lines, 5 integers each: '
       write(6,*)'  first 4 (integer*8): seeds; the last 1 or 0, depending on whether to assign '
       write(6,*)'  a rank-dependent seed or not '
       write(6,*)
    end if
    

    open(4,file='random.inp',status='old')
    rewind(4);      
    do iseed=1,ndiffseeds
       read(4,*,end=99)seed_read(:,iseed),rank_dependent(iseed)

       if(rank==0) &
            write(6,*)' iseed ',iseed,' seed_read ',     seed_read(:,iseed),&
                                      ' rank_dependent ',rank_dependent(iseed)
    end do
99  continue
    close(4)

    mdiffseeds = iseed-1
    call check_lele(1,mdiffseeds,ndiffseeds,' 1 mdiffseeds ndiffseeds ')

    do iseed=1,mdiffseeds
       if(rank_dependent(iseed)==1) then
          do i=1,nlseed 
             j1 = i
             j2 = 921778613
             j3 = rank
             j4 = iseed
             j5 = 6938938
             j8 = i+2*i**2
             cnt8 = cnt
             big8 = 5768301
         
             seed_read_wcounter(i,iseed) = &
                     seed_read( i,iseed)  + j4*j5 + j8*cnt8*big8
             seed_array(i,iseed) = seed_read_wcounter(i,iseed)+ j1*j2*j3

          enddo
       else
          do i=1,nlseed 
             j1 = i
             j2 = 693372617
             j4 = iseed
             cnt8 = cnt
             big8 = 41922341

             seed_read_wcounter(i,iseed) = &
             seed_read(i,iseed)+ j1*j2 + j1*cnt8*big8 !  int8(iseed)*int8(693)

             seed_array(i,iseed) = seed_read_wcounter(i,iseed)
          enddo
       end if
    end do

    line_seed = 1
    call ran_ps_putseed(seed_array(:,line_seed))

  end subroutine read_initiate_seeds

  subroutine read_counter(cnt)
    use mpi_lib_ours, only : bcast_scalar_i, rank
    implicit none
    integer cnt
    if(rank==0) then
       write(6,*)' reading counter from counter.inp ; ensure it exists! '
       open(2,file='counter.inp',status='old')
       rewind(2)
       read(2,*)cnt
       close(2)
       write(6,*)' counter= ',cnt
    end if
    call bcast_scalar_i(cnt)
  end subroutine read_counter
end module rand_seed_mod

function ran_ps()
  use mod_kiss,      only : kiss_uniform
  use mpi_lib_ours,  only : rank
  use rand_seed_mod, only : first,line_seed, mdiffseeds, seed_array, read_initiate_seeds
  implicit none
  integer       :: i
  real*8   r
  real*8   ran_ps
  
  if(first) then
     call read_initiate_seeds  
     first=.false. 
  end if

  call check_lele(1,line_seed, mdiffseeds,' 1 line_seed, mdiffseeds ')
  call ran_ps_putseed(seed_array(:,line_seed))

  call kiss_uniform(r)
  ran_ps = r
  call ran_ps_getseed(seed_array(:,line_seed))

end function ran_ps

subroutine ran_ps_getseed(seed)
  use mod_kiss, only : kiss_seed_extract_DN
  use rand_seed_mod, only : nlseed
  implicit none
  integer(kind=8) seed(nlseed)
  call check(nlseed, 4     ,' nseed, four   ')
  call kiss_seed_extract_DN(seed(1),seed(2),seed(3),seed(4))
end subroutine ran_ps_getseed

subroutine ran_ps_putseed(seed)
  use mod_kiss, only : kiss_seed
  use rand_seed_mod, only : nlseed
  implicit none
  integer(kind=8) seed(nlseed)
  call check(nlseed, 4     ,' nlseed, four   ')
  call kiss_seed(seed(1), seed(2), seed(3), seed(4))
end subroutine ran_ps_putseed

subroutine ran_ps_putseed_intoarray(seed, line)  ! sets line-seed to line, sets seed_array(:,line) too seed
  use mod_kiss, only : kiss_seed_extract_DN
  use rand_seed_mod, only : nlseed,line_seed,seed_array, mdiffseeds
  implicit none
  integer line
  integer(kind=8) seed(nlseed)
  call check(nlseed, 4     ,' nseed, four   ')
  call check_lele(1,line,mdiffseeds,' 1 line mdifseds ')
  line_seed=line
  seed_array(:,line)=seed
end subroutine ran_ps_putseed_intoarray

subroutine set_annweights(w,n_weight,seed)
    implicit none
    real(8), intent(inout):: w(n_weight)
    integer, intent(in)::  n_weight
    integer:: i,seed,clock
    real(8)::  tt,ampl_rand
    ampl_rand=2.d-1
    if (seed<1) then
        call random_number(tt)
        call system_clock(count=clock)    ! call itime(timeArray)     
        seed = clock*tt
    endif
    call srand(seed)
    write(*,*)"seed =" ,seed
    do i=1,n_weight
        call random_number(w(i))
        !w(i)=rand()
        tt=2.d0*ampl_rand
        w(i)=(w(i)-0.5d0)*tt
    enddo
    write(*,*) 'number of ANN wights:             ',n_weight
end subroutine set_annweights
!************************************************************************************
subroutine io_read_nconf(nat,nconf,file_input)
    implicit none
    integer:: stat, nat, nconf
    integer:: ios, iat, iconf
    character (15):: tmp
    character (16):: file_input

    open(unit=11,file=trim(file_input),status='old',iostat=ios)
    nconf=0
    do
        nconf=nconf+1
        read(11, *, IOSTAT=stat) nat 
        if(IS_IOSTAT_END(stat)) exit
        read(11, *, IOSTAT=stat)  tmp 
        do iat = 1,nat    
            read(11, *, IOSTAT=stat) 
        enddo
    end do
    nconf=nconf-1
    close (11)
end subroutine io_read_nconf
!*******************************************************************************
subroutine io_read_conf(nat,nconf,epot,typat,bonds,finput,bound,epot0,rat)
    implicit none
    integer:: ios, iconf,nconf,nat,iat
    integer:: ibond, nbond,jat
    real(8):: epot(nconf),rat(3,nat,nconf),epot0
    real(8):: bonds(nat*(nat-1)/2,nconf)    !scaled rat
    real(8):: bound(2) ,two_over_diff
    real(8):: dx, dy, dz 
    character(5)::tmp,typat(nat,nconf)
    character(5)::finput
    character(16)::file_name

    if (trim(finput)=="train" .or. trim(finput)=="valid") then
        write(file_name,"(a7,a5,a4)")'posinp_',trim(finput),'.xyz'
    else
        write(file_name,"(a10)")'posinp.xyz'
    endif
    open(unit=11,file=file_name,status='old',iostat=ios)
    do iconf=1, nconf
        read(11, *) nat,tmp,tmp,tmp,epot(iconf) 
        epot(iconf)=epot(iconf)-epot0
        read(11, *)  tmp 
        do iat = 1,nat    
            read(11, *)typat(iat,iconf),rat(1,iat,iconf),rat(2,iat,iconf),rat(3,iat,iconf) 
        enddo
    end do

end subroutine io_read_conf
!******************************************************************
subroutine set_input_layer(nat,nconf,rat,bonds,bound,epot0,finput)
    implicit none
    integer:: ios, iconf,nconf,nat,iat
    integer:: ibond, nbond,jat
    real(8):: epot(nconf),rat(3,nat,nconf),epot0
    real(8):: bonds(nat*(nat-1)/2,nconf)    !scaled rat
    real(8):: bound(2) ,two_over_diff
    real(8):: dx, dy, dz 
    character(5)::tmp,typat(nat,nconf)
    character(5)::finput

    nbond = nat*(nat-1)/2
    ibond=0
    do iconf=1, nconf
        ibond=0
        do iat = 1,nat
            do jat = iat+1, nat
                ibond=ibond+1
                dx= rat(1,iat,iconf)-rat(1,jat,iconf)
                dy= rat(2,iat,iconf)-rat(2,jat,iconf)
                dz= rat(3,iat,iconf)-rat(3,jat,iconf)
                bonds(ibond,iconf)=sqrt(dx**2+dy**2+dz**2)
            enddo
        enddo
    enddo
    if (trim(finput)=="train") then
        bound(2)=maxval(bonds)
        bound(1)=minval(bonds)
        two_over_diff=2.d0/(bound(1)-bound(2))
        do iconf=1, nconf
            do ibond=1,nbond
                bonds(ibond,iconf)=(bonds(ibond,iconf)-bound(2))*two_over_diff-1
            enddo
        end do
        open(unit=13,file="params.ann",status='unknown',iostat=ios)
            write(13,*) "epot:   ",epot0
            write(13,*) "bound:  ",bound(1),bound(2)
        close(13)
    else!if (finput=="valid") then
       two_over_diff=2.d0/(bound(1)-bound(2))
       do iconf=1, nconf
           do ibond=1,nbond
               bonds(ibond,iconf)=(bonds(ibond,iconf)-bound(2))*two_over_diff-1
           enddo
       end do
    endif
    close (11)
end subroutine set_input_layer

!*****************************************************************************************
subroutine read_params(nat,n_weight,nn,w,epot0,bound)
    integer:: ios, i,nat,n_weight,nn(0:3)
    real(8):: bound(2),epot0 
    real(8):: w(600)
    character(10) :: tmp
        open(unit=13,file="params.ann",status='unknown',iostat=ios)
            read(13,*) tmp ,epot0
            read(13,*) tmp ,bound(1),bound(2)
            read(13,*) tmp ,nn(1), nn(2)
            read(13,*) tmp  
            nn(0)=nat*(nat-1)/2
            nn(3)=1
            n_weight=nn(0)*nn(1) + nn(1)*nn(2) + nn(2) ! number of "a" parameters 
            n_weight=n_weight+nn(1) + nn(2) + 1               ! number of "b" parameters 
            write(*,*)nn(:)
            do i = 1,n_weight
                read(13,*) w(i)
            enddo
        close(13)
end subroutine read_params
!****************************************************************************
subroutine write_params(n_weight,nn,w)
    integer:: ios, i,nat,n_weight,nn(0:3)
    real(8):: bound(2),epot0 
    real(8):: w(n_weight)
    open(unit=13,file="params.ann",status='old',access='append',iostat=ios)
        write(13,*) "nodes:  ",nn(1), nn(2)
        write(13,*) "weights : "
        do i = 1,n_weight
            write(13,*) "  ",w(i)
        enddo
    close(13)
end subroutine write_params

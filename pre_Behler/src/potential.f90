subroutine potential(nn,task,nat,nconf)
    implicit none
    integer:: ios, m, iann,i
    integer, parameter :: nl=3 ! nl = Number of layers
    integer:: nat, nconf
    integer:: n_weight 
    integer:: iat, iconf
    integer:: nn(0:nl)  ! nn = Number of nodes 
    integer:: nbond
    real(8):: bound(2),epot0
    real(8):: wt(600)
    character (20):: task
    real(8) :: bonds(nat*(nat-1)/2,nconf), epot(nconf)
    real(8) :: rat(3,nat,nconf)
    real(8) :: fat(3,nat,nconf)
    real(8),allocatable :: epotd(:,:)
    real(8),allocatable :: w(:)
    character(5) :: typat(nat,nconf) 
    real(8):: a(140,140,nl), b(140,nl)

    nbond = nat*(nat-1)/2

    call read_params(nat,n_weight,nn,wt,epot0,bound)

    allocate(epotd(n_weight,nconf))
    allocate(w(n_weight))
    w(1:n_weight)=wt(1:n_weight)

    call io_read_conf(nat,nconf,epot,typat,bonds,"potential",bound,epot0,rat)
    call set_input_layer(nat,nconf,rat,bonds,bound,epot0,"potential")

    epot= 0.d0
    call convert_w_ann(n_weight,w,a,b,nl,nn)

    if (trim(task)=="single_point") then
        call force_energy(nl,nn,n_weight,nat,nconf,bonds,rat,bound,a,b,epot,fat)
    else if(trim(task)=="relax") then
        nconf=1
        call relaxation(nl,nn,n_weight,nat,1,bonds(:,1),rat(:,:,1),bound,a,b,epot(1),fat(:,:,1),epot0)
    endif


    open(unit=11,file='posout.xyz',status='unknown',iostat=ios)
    do iconf=1,nconf
        write(11,*)nat, "angstroem  struct00001  epot= " , epot(iconf)+epot0
        write(11,*)"free"
        do iat = 1,nat
            write(11,*)typat(iat,nconf), rat(:,iat,nconf)
        enddo
    enddo
    write(*,*)"Done. "
    write(*,*)"      results are printed in posout.xyz    "
    write(*,*)
end subroutine potential
!*****************************************************************************************   
subroutine relaxation(nl,nn,n_weight,nat,nconf,bonds,rat,bound,a,b,epot,fat,epot0)
    implicit none
    integer:: nl, nn(0:nl), n_weight,jat
    integer:: iat, nat, nconf,iter, ibond, nbond
    real(8)::  epot(1) ,epotd(n_weight)
    real(8):: a(140,140,nl), b(140,nl), x(140,nl), y(140,0:nl)
    real(8):: bonds(nat*(nat-1)/2,1)
    real(8):: fat(3,nat,1)
    real(8):: rat(3,nat,1)
    character(15):: event
    real(8):: rmse, errmax, tt,bound(2)
    integer:: iconf, ierrmax  
    real(8):: d(nat*(nat-1)/2)
    real(8):: dx, dy, dz ,r , de 
    real(8):: epot_old, fnrm, fmax, alpha,epot0 

    call force_energy(nl,nn,n_weight,nat,nconf,bonds,rat,bound,a,b,epot,fat)
    alpha=1.d-2
    iter=0
    fmax=1.d-5
    epot_old=epot(1)
    open(unit=21,file="monitor_sd",status='unknown')
    write(21,'(a80)')"#  iter         energy         energy difference   norm_forces    maximum force"
    do! while
         do iat=1,nat
             rat(:,iat,1)=rat(:,iat,1)+alpha*fat(:,iat,1)
        enddo
        call set_input_layer(nat,nconf,rat,bonds,bound,0.d0,"potential")
        call force_energy(nl,nn,n_weight,nat,nconf,bonds,rat,bound,a,b,epot,fat)
        iter=iter+1
        de=epot(1)-epot_old 
        epot_old=epot(1)

        fnrm=sqrt(sum(fat**2))

        if  (fnrm<fmax ) exit
        write(21,'(i5,es24.15,3es17.5)')iter,epot+epot0,de,fnrm,maxval(abs(fat))
    enddo
    write(*,*)"steepest decent"
    write(*,'(a7,a24,3a14)')"iter","epot      ","different ","fnrm   ","maxval(fat)"
    write(*,'(i7,es24.15,3es14.5)')iter,epot+epot0,de,fnrm,maxval(fat)
    close(21)

end subroutine relaxation
!***********************************************************************************
subroutine force_energy(nl,nn,n_weight,nat,nconf,bonds,rat,bound,a,b,epot,fat)
    implicit none
    integer:: nl, nn(0:nl), n_weight,jat
    integer:: iat, nat, nconf,iter, ibond, nbond
    real(8)::  epot(nconf) ,epotd(n_weight,nconf)
    real(8):: a(140,140,nl), b(140,nl), x(140,nl), y(140,0:nl)
    real(8):: bonds(nat*(nat-1)/2,nconf)
    real(8):: fat(3,nat,nconf)
    real(8):: rat(3,nat,nconf)
    character(15):: event
    real(8):: rmse, errmax, tt,bound(2)
    integer:: iconf, ierrmax  
    real(8):: d(nat*(nat-1)/2,nconf)
    real(8):: dx, dy, dz ,r 
    fat = 0.d0
    do iconf=1,nconf
        event='potential'
        call cal_ann_atombased(nl,nn,n_weight,nat,bonds(1,iconf),a,b,epot(iconf),epotd(1,iconf),d(1,iconf),event)
        tt=2.d0/(bound(1)-bound(2))
        d(:,iconf)=d(:,iconf)*tt
        ibond=0
        do iat = 1,nat
            do jat = iat+1, nat
                ibond=ibond+1
                dx= rat(1,iat,iconf)-rat(1,jat,iconf)
                dy= rat(2,iat,iconf)-rat(2,jat,iconf)
                dz= rat(3,iat,iconf)-rat(3,jat,iconf)
                r=sqrt(dx**2+dy**2+dz**2)
                fat(1,iat,iconf)=fat(1,iat,iconf)-d(ibond,iconf)*dx/r
                fat(2,iat,iconf)=fat(2,iat,iconf)-d(ibond,iconf)*dy/r
                fat(3,iat,iconf)=fat(3,iat,iconf)-d(ibond,iconf)*dz/r
                fat(1,jat,iconf)=fat(1,jat,iconf)+d(ibond,iconf)*dx/r
                fat(2,jat,iconf)=fat(2,jat,iconf)+d(ibond,iconf)*dy/r
                fat(3,jat,iconf)=fat(3,jat,iconf)+d(ibond,iconf)*dz/r
            enddo
        enddo
    enddo
end subroutine force_energy
!!************************************************************************************
subroutine cal_ann_atombased(nl,nn,n_weight,nat,bonds,a,b,epot,epotd,d,event)
    implicit none
    integer::  nat 
    integer:: nbond 
    integer:: nl, nn(0:nl), n_weight
    integer:: iat, jat,l,ib 
    real(8):: epot ,epotd(n_weight)
    real(8):: a(140,140,nl), b(140,nl), x(140,nl), y(140,0:nl)
    real(8):: ad(140*140,nl), bd(140,nl)
    real(8):: bonds(nat*(nat-1)/2),d(nat*(nat-1)/2)
    character(15):: event
    nbond = nat*(nat-1)/2
    do ib = 1,nbond
        y(ib,0) = bonds(ib)
    enddo
    if(trim(event)=='potential' .or. trim(event)=='evalu') then
        call cal_architecture_2hiddenlayer(nl,nn,a,b,x,y,epot,d)
    elseif(trim(event)=='train') then
        call cal_architecture_der_2hiddenlayer(nl,nn,a,b,x,y,ad,bd,epot)
        call convert_ann_epotd(nn,nl,ad,bd,n_weight,epotd(:))
    else
        stop 'ERROR: undefined content for event'
    endif
        
end subroutine cal_ann_atombased

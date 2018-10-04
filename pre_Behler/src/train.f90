!**********************************************************************************************************************
subroutine train(seed,epot0,nn,nat,nconf_train,nconf_valid) 

!> This subroutine works according to the training flowchart, that is: 
!>-----------------------------------------------------------------------------
!>       Input training data 
!>
!>            Read positions and DFT energies of atomic sructures (Here: H_2O):
!>                  By calling  "io_read_conf" and subroutine
!>-----------------------------------------------------------------------------
!>       Initialization
!>
!>            Construct input layer by determining bonds between atoms:
!>                  By calling  "set_input_layer" subroutine
!>            Set random weights:
!>                  By calling "set_annweights" subroutine
!>-----------------------------------------------------------------------------
!>       Computing 
!>            
!>            Decrease the root-mean-square error (RMSE) of energy by updating the weights using optimization methods:
!>                  By calling "prefit_optimizer_SD" and "optimizer_LM" subroutines 
!>                  both of them calls "cal_ann_atombased" subroutine in order to 
!>                  construct ANN architecture and compute E (potential energy)
!>-----------------------------------------------------------------------------                                    
!>        Get error value 
!>            
!>            Compute error value and print:
!>                  By calling "ann_evaluate" subroutine 
!>-----------------------------------------------------------------------------
!>        E < limit value 
!>            
!>            check if error is enough small:
!>                  By calling "ann_evaluate" subroutine
!>                  if yes => end and calling the "write_params".
!>                  if no  => change the weights using an algorithm and do previous steps iteratively.
!>-----------------------------------------------------------------------------
!*****************************************************************************************************
!   1- Define the variables that you need

    use mod_parlm, only: typ_parlm
    implicit none
    ! input variables
    integer, parameter :: nl=3                      ! Number of layers, 
    integer:: nat                                   ! Number of atom in each configuration
    integer:: nconf_train, nconf_valid              ! Number of configuration in training and validation data set
    integer:: seed, nn(0:nl)                        ! A random number, Number of nodes 
    real(8):: epot0                                 ! Energy   
    ! local variables 
    integer:: n_weight, m                           ! Number of weights, ...
    real(8):: rat_train(3,nat,nconf_train)          ! Positions of atoms for training data set
    real(8):: bonds_train(nat*(nat-1)/2,nconf_train)! Bond lengths of training data set
    real(8):: epot_train(nconf_train)               ! DFT energy of  training data set
    real(8):: rat_valid(3,nat,nconf_valid)          ! Positions of atoms for validating data set
    real(8):: bonds_valid(nat*(nat-1)/2,nconf_valid)! Bond lengths of validating data set
    real(8):: epot_valid(nconf_valid)               ! DFT energy of  validating data set
    integer:: nbond                                 ! Number of bonds ( H_H and H-O )
    character(5):: typat_train(nat,nconf_train)     ! The type of atoms (H or O) for training data set
    character(5):: typat_valid(nat,nconf_valid)     ! The type of atoms (H or O) for validating data set
    character (15):: event                          ! A string for choosing 
    real(8):: epot_ann(nconf_train)                 ! Energy computed from ANN process
    real(8):: a(140,140,nl), b(140,nl)              ! Weights "a" and "b"
    real(8),allocatable :: epotd(:,:),costd(:)      ! Derivatives of energy and cost function respect to weights
    real(8),allocatable :: w(:)                     ! An array for all weights, includes "a" and "b"
    type(typ_parlm):: parlm                         ! A derived type for Levenberg-Marquadt variables 
    integer:: iter                                  ! Number of iterations in optimization process
    real(8):: bound(2)                              ! An array for normalization of length of bonds.

    call random_seed()

    !********************************************************************************
!   2- Compute the number of bonds from number of atoms that it is defined 

    nbond = nat*(nat-1)/2
!   ********************************************************************************
!   3- Call io_read_conf for reading positions and DFT energies of atomic sructures, 
!      both training and validating data set 

    call io_read_conf(nat,nconf_train,epot_train,typat_train,bonds_train,"train",bound,epot0,rat_train)
    call io_read_conf(nat,nconf_valid,epot_valid,typat_valid,bonds_valid,"valid",bound,epot0,rat_valid)
    !********************************************************************************
!   4- call "set_input_layer" for setting bond lengths as input layer nodes instead of 
!      positions of atoms, both training and validating data set    

    call set_input_layer(nat,nconf_train,rat_train,bonds_train,bound,epot0,"train")
    call set_input_layer(nat,nconf_valid,rat_valid,bonds_valid,bound,epot0,"valid")
    !********************************************************************************
!   5- Compute the number of ANN weights using number of nodes and layers

    nn(0)=nbond
    nn(3)=1
    n_weight=nn(0)*nn(1) + nn(1)*nn(2) + nn(2)        ! number of "a" parameters 
    n_weight=n_weight+nn(1) + nn(2) + 1               ! number of "b" parameters 
    write(*,*) 'number of ANN wights:             ',n_weight
    !********************************************************************************
!   6- Do allocate the variables

    allocate(w(n_weight))
    allocate(epotd(n_weight,nconf_train),costd(n_weight))
    !********************************************************************************
!   7- Call "set_annweights"

    call set_annweights(w,n_weight,seed)

    !********************************************************************************
!   8- Call "prefit_optimizer_SD" to change the random weights properly 

    call prefit_optimizer_SD (n_weight,w,nl,nn,nconf_train,nat,bonds_train,epot_train)
    !********************************************************************************
!   9- Import a loop for doing optimization iteratively. 
!      In the loop you should call "optimizer_LM" and then call "ann_evaluate" to 
!      check if error is enough small 

    iter=0
    do !while
        call optimizer_LM (n_weight,w,nl,nn,nconf_train,nat,bonds_train,epot_train,iter,parlm,m,a,b)
        call ann_evaluate(nl,nn,n_weight,nconf_train,"train",nat,bonds_train,epot_train,a,b,iter)
        call ann_evaluate(nl,nn,n_weight,nconf_valid,"valid",nat,bonds_valid,epot_valid,a,b,iter)
        iter=iter+1
        if(parlm%finish) exit
        if (iter==201) exit
    enddo
!   ********************************************************************************
!   10- Call "write_params" to have final ANN weights and so the functionality of ANN. 

    call write_params(n_weight,nn,w)
    !*********************************************************************************
    write(*,*) 'iter= ',iter
    write(*,*) 'info = ',parlm%info
    call final_lmder_modified(parlm)

end subroutine train
!**********************************************************************************************************************
subroutine ann_evaluate(nl,nn,n_weight,nconf,data_set,nat,bonds,epot_ref,a,b,iter)
!> This subroutine compute the error value and check it if is enough small
    implicit none
    integer:: nl, nn(0:nl), n_weight
    integer:: iat, nat, nconf,iter , ifile
    real(8):: epot, epot_ref(nconf) ,epotd(n_weight,nconf)
    real(8):: a(140,140,nl), b(140,nl), x(140,nl), y(140,0:nl)
    real(8):: bonds(nat*(nat-1)/2,nconf)
    real(8):: fat(3,nat,nconf)
    character(15):: event
    real(8):: rmse, errmax, tt
    integer:: iconf, ierrmax  
    integer:: ilarge1, ilarge2, ilarge3, iunit, ios
    character(28):: frmt1='(a,i6,2f10.3,i7)'
    character(28):: frmt2='(a,i6,2e10.1,i7)'
    character(28):: frmt
    character(15):: filename, filename_err
    character(5):: data_set
    real(8):: d(nat*(nat-1)/2,nconf)
    rmse=0.d0
    errmax=0.d0
    ierrmax=0
    ilarge1=0
    ilarge2=0
    ilarge3=0
    event='evalu'
    ifile =11

!    write(filename,'(a12,i3.3)') 'detailed_err',iter
!    if(data_set=="train") then
!        open(unit=100,file=trim(filename),status='unknown')
!    elseif(data_set=="valid") then
!        open(unit=100,file=trim(filename),status='old',access='append')
!    endif
    write(filename_err,"(a4,a5)") "err_",data_set
    if (iter == 0) then
        open(unit=ifile,file=trim(filename_err),status='unknown')
    else
        open(unit=ifile,file=trim(filename_err),status='old',access='append')
    endif

    configuration: do iconf=1,nconf
        call cal_ann_atombased(nl,nn,n_weight,nat,bonds(:,iconf),a,b,epot,epotd(:,iconf),d(:,iconf),event)

        tt=abs(epot-epot_ref(iconf))/nat
!       write(100,'(i7,3es16.5,a10)') iconf,epot,epot_ref(iconf),tt,trim(data_set)
        if(tt>1.d-2) ilarge1=ilarge1+1
        if(tt>1.d-3) ilarge2=ilarge2+1
        if(tt>1.d-4) ilarge3=ilarge3+1
        if(tt>errmax) then
            errmax=tt
            ierrmax=iconf
        endif
        rmse=rmse+tt**2
    enddo configuration

    rmse=sqrt(rmse/real(nconf,8))
        rmse=rmse*1.d3
        errmax=errmax*1.d3
        if(rmse>99999.d0) then
            frmt=frmt2
        else
            frmt=frmt1
        endif
        write(ifile,frmt)"iter,RMSE(meV/atom),MAX_error,iconf", iter,rmse,errmax,ierrmax
!        close(100)
        close(ifile)
end subroutine ann_evaluate
!**********************************************************************************************************************
subroutine prefit_optimizer_SD (n_weight,w,nl,nn,nconf_train,nat,bonds_train,epot_train)
!> This subroutinre changes the random weights properly.
    implicit none
    integer::  i
    integer:: nl                  ! nl = Number of layers
    integer:: nat, nconf_train
    integer:: n_weight
    integer:: iconf
    integer:: nn(0:nl)            ! nn = Number of nodes 
    integer:: iter
    real(8):: de, fmax=0.1, fnrm, alpha= 1.d-2, cost0
    character (15):: event
    real(8):: bonds_train(nat*(nat-1)/2,nconf_train), epot_train(nconf_train)
    real(8):: epot_ann(nconf_train)
    real(8):: d_train(nat*(nat-1)/2,nconf_train)
    real(8):: epotd(n_weight,nconf_train),costd(n_weight),w(n_weight)
    real(8):: a(140,140,nl), b(140,nl)
    real(8):: cost 

    !> \f$ cost: \chi = \frac{1}{Nconf} \sum_i^{Nconf} (E_i^{ANN}-E_i^{DFT})^2  \f$ 


    !> \f$ costd: \frac{\partial \chi}{\partial W_i}= \frac{1}{Nconf} \sum_i^{Nconf} 2 (E_i^{ANN}-E_i^{DFT}) \frac{\partial
    !!E_i^{ANN}}{\partial W_i} \f$

    cost0= 0.d0
    iter=0
    do! while
        call convert_w_ann(n_weight,w,a,b,nl,nn)
        cost=0.d0
        costd=0.d0
        do iconf=1,nconf_train
            event='train'
            call cal_ann_atombased(nl,nn,n_weight,nat,bonds_train(:,iconf),a,b,epot_ann(iconf),&
                                   epotd(:,iconf),d_train(:,iconf),event)
            cost = cost+ (epot_train(iconf)-epot_ann(iconf))**2
            costd(1:n_weight)=costd(1:n_weight)-epotd(1:n_weight,iconf)*(epot_train(iconf)-epot_ann(iconf))*2.d0
        enddo
        cost=cost/nconf_train
        costd(1:n_weight)=costd(1:n_weight)/nconf_train
        de=cost-cost0 
        cost0=cost

        fnrm=sqrt(sum(costd**2))
        if  (fnrm<fmax .or. iter>30) exit

        do i=1,n_weight
             w(i)=w(i)-alpha*costd(i)
        enddo
        iter=iter+1
    enddo
end subroutine prefit_optimizer_SD 
!***********************************************************************************************************************
subroutine optimizer_LM (n_weight,w,nl,nn,nconf_train,nat,bonds_train,epot_train,iter,parlm,m,a,b)
!> This subroutine updates the weights to minimize the cost function
    use mod_parlm, only: typ_parlm
    implicit none
    type(typ_parlm):: parlm
    integer:: m,i
    integer:: nl          ! nl = Number of layers
    integer:: nat, nconf_train
    integer:: n_weight 
    integer:: iconf
    integer:: nn(0:nl)            ! nn = Number of nodes 
    integer:: iter
    character (15):: event
    real(8):: bonds_train(nat*(nat-1)/2,nconf_train), epot_train(nconf_train)
    real(8):: epot_ann(nconf_train)
    real(8):: d_train(nat*(nat-1)/2,nconf_train)
    real(8):: a(140,140,nl), b(140,nl)
    real(8):: epotd(n_weight,nconf_train),w(n_weight)

    if (iter ==0 .and. parlm%allocate_init) then
        m=nconf_train
        parlm%n=n_weight
        call init_lmder_modified(parlm,m,m)
        parlm%x(1:parlm%n)=w(1:parlm%n)

        parlm%allocate_init = .false.
    else 
        do
            call lmder_modified(parlm,m,m)
            if(parlm%finish) exit 
            if(parlm%icontinue==700) then
                call convert_w_ann(n_weight,parlm%wa2,a,b,nl,nn)
            else
                call convert_w_ann(n_weight,parlm%x,a,b,nl,nn)
            endif

            if(parlm%iflag==1) then
                event='evalu'
                do iconf=1,nconf_train
                    call cal_ann_atombased(nl,nn,n_weight,nat,bonds_train(:,iconf),a,b,epot_ann(iconf)&
                                          ,epotd(:,iconf),d_train(:,iconf),event)
                    if(parlm%icontinue==700) then
                        parlm%wa4(iconf)=(epot_train(iconf)-epot_ann(iconf))**2
                    else
                        parlm%fvec(iconf)=(epot_train(iconf)-epot_ann(iconf))**2
                    endif
                enddo
            elseif(parlm%iflag==2) then
                event='train'
                do iconf=1,nconf_train
                    call cal_ann_atombased(nl,nn,n_weight,nat,bonds_train(:,iconf),a,b,epot_ann(iconf)&
                                          ,epotd(:,iconf),d_train(:,iconf),event)
                    parlm%fjac(iconf,1:n_weight)=-epotd(1:n_weight,iconf)*(epot_train(iconf)-epot_ann(iconf))*2.d0
                enddo
            elseif(parlm%iflag==0) then
                exit
            endif
        enddo
    endif
    w(1:parlm%n)=parlm%x(1:parlm%n)
end subroutine optimizer_LM
!**********************************************************************************************************************

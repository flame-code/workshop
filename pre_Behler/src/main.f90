!****************************************************************************************
program ANN
    implicit none
    integer, parameter :: nl=3      ! nl = Number of layers
    integer:: seed, nn(0:nl)        ! nn = Number of nodes 
    character (20):: tmp, task, event
    real(8):: bound(2),epot0
    integer:: ios
    integer:: nat,nconf_train,nconf_valid,nconf 
    character (16):: file_input

    open(unit=1,file='input.ann',status='old',iostat=ios)
    read(1,*) tmp, nn(1), nn(2)
    read(1,*) tmp, task
    read(1,*) tmp, seed
    read(1,*) tmp, epot0
    close (1)

    if (trim(task)=="train") then
        call io_read_nconf(nat,nconf_valid,"posinp_valid.xyz")
        call io_read_nconf(nat,nconf_train,"posinp_train.xyz")
        write(*,*) "Number of training data points:     " ,nconf_train
        write(*,*) "Number of vlidating data points:    " ,nconf_valid
    else
        file_input= "posinp.xyz"
        call io_read_nconf(nat,nconf,file_input)
        write(*,*) "Number of data points:    " ,nconf
    end if

    if (trim(task)=="train") then
        call train(seed,epot0,nn,nat,nconf_train,nconf_valid) 
    else 
        call potential(nn,task,nat,nconf) 
    end if

end program ANN
!******************************************************************************************

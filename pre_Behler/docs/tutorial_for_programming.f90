!   PLEASE
!   Write a subroutine that its input arguments are : 
!   (a random number, energy reference, number of layers, number of nodes per layer,  
!   numbar of atoms, number of training set data, number of validating set data)
 
    subroutine train(seed,epot0,nn,nat,nconf_train,nconf_valid) 
!   
!  This subroutine works according to the training flowchart, that is: 
!
!       |----------------------------------------------------------------------------|
!       |       Initialization                                                       |
!       |                                                                            |
!       |            Set random weights:                                             |
!       |                  By calling "set_annweights" subroutine                    |
!       |            Get ann parameters include number of nodes and layers as input &|
!       |                  calling  "set_input_layer" subroutine                     |
!       |----------------------------------------------------------------------------|
!       |       Input training data                                                  |    
!       |                                                                            |
!       |            Read positions and DFT energies of atomic sructures(Here: H_2O) |
!       |                  By calling  "io_read_conf" subroutine                     |
!       |----------------------------------------------------------------------------|
!       |        Computing                                                           |        
!       |                                                                            |               
!       |           ANN architecture and compute E (potential energy):               |   
!       |                  By calling "prefit_optimizer_SD" and "optimizer_LM"       |   
!       |                  subroutines that call "cal_ann_atombased" subroutine      |
!       |                  to do this computation                                    |   
!       |----------------------------------------------------------------------------|                                                       
!       |        Get error value                                                     |
!       |                                                                            |
!       |            Compute error value:                                            |
!       |                  By calling "ann_evaluate" subroutine                      |
!       |----------------------------------------------------------------------------|
!       |        E < limit value                                                     |
!       |                                                                            |
!       |            check if error is enough small:                                 |
!       |                  By calling "ann_evaluate" subroutine                      |
!       |                  if yes => end and calling the "write_params".             | 
!       |                  if no  => change the weights using an algorithm and do    |
!       |                            previous steps iteratively.                     |
!       |----------------------------------------------------------------------------|
!       |        Change the weights                                                  |
!       |           Update the weights using one/two optimization method/methods:    | 
!       |                  By calling "prefit_optimizer_SD" and "optimizer_LM"       |
!       |                                                                            |
!       |----------------------------------------------------------------------------|
!   
!   You can do this using the following instructions step by step:
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
    
!   ********************************************************************************
!   2- Compute the number of bonds from number of atoms that it is defined 

    .......................

!   ********************************************************************************
!   3- Call io_read_conf for reading positions and DFT energies of atomic sructures, 
!      both training and validating data set 
    
    .........................................
    .........................................

!   ********************************************************************************
!   4- call "set_input_layer" for setting bond lengths as input layer nodes instead of 
!      positions of atoms, both training and validating data set    

    .........................................
    .........................................

!   ********************************************************************************
!   5- Compute the number of ANN weights using number of nodes and layers

    ...................
    ...................
    ..............................
    ....................................

!   ********************************************************************************
!   6- Do allocate the variables
    
    ..........................
    ...............................

!   ********************************************************************************
!   7- Call "set_annweights"

    ...........................

!   ********************************************************************************
!   8- Call "prefit_optimizer_SD" to change the random weights properly 
    
    ................................

!   ********************************************************************************
!   9- Import a loop for doing optimization iteratively. 
!      In the loop you should call "optimizer_LM" and then call "ann_evaluate" to 
!      check if error is enough small 
    
    ..........................
    ................................
    ......................................
    ......................................

!   ********************************************************************************
!   10- Call "write_params" to have final ANN weights and so the functionality of ANN. 

    ................................

!   *********************************************************************************
    
    call final_lmder_modified(parlm)
end subroutine train
!**********************************************************************************************************************

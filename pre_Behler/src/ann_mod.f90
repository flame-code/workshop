!*****************************************************************************************
module mod_parlm
    implicit none
    type typ_parlm
        real(8):: ftol=1.d-8
        real(8):: xtol=1.d-8
        real(8):: gtol=1.d-8
        real(8):: factor=100.d0
        integer:: maxfev=1000
        integer:: nprint=1
        integer:: n=0
        integer:: mode, info, nfev, njev
        integer:: iter
        integer:: icontinue
        integer:: iflag
        logical:: finish
        real(8):: epsmch
        real(8):: fnorm
        !real(8):: fnorm1
        real(8):: xnorm
        real(8):: gnorm
        real(8):: pnorm
        real(8):: par
        real(8):: delta
        real(8):: actred
        real(8):: prered
        real(8):: ratio
        real(8), allocatable:: wa1(:), wa2(:), wa3(:), wa4(:), qtf(:)
        real(8), allocatable:: x(:)
        real(8), allocatable:: fvec(:)
        real(8), allocatable:: fjac(:,:)
        real(8), allocatable:: diag(:)
        integer, allocatable:: ipvt(:)
        logical:: allocate_init = .true.
    end type typ_parlm
end module mod_parlm
!*****************************************************************************************

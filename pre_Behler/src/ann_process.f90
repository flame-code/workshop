subroutine cal_architecture_2hiddenlayer(nl,nn,a,b,x,y,epot,d)
    implicit none    
    integer, intent(in):: nl !number of hidden layer plus one
    integer, intent(in):: nn(0:nl) !number of nodes in each hiddenlayer
    real(8), intent(inout):: a(140,140,nl), b(140,nl), x(140,nl), y(140,0:nl)
    real(8), intent(inout):: epot
    real(8), intent(inout):: d(nn(0))
    !local variables
    real(8):: yd(140,0:nl)
    integer:: i, j, k 
    real(8):: tt
    real(8):: c(140)
    if(nl/=3) then
        write(*,*) 'ERROR: this routine works only for nl=3, while nl= ',nl
        stop
    endif
!> \f$ E = f_{1}^{3} \bigg( b_{1}^{3} + \sum_{k=1}^{P} a^{23}_{k1} .f_{k}^{2} \bigg( b_{k}^{2} + \sum_{j=1}^{L} a^{12}_{jk}. f_{j}^{1}
!>\bigg( b_{j}^{1} + \sum_{i=1}^{M} a^{01}_{ij}.G_i  \bigg)\bigg)\bigg) \f$ 

    !-------------------------------------------------------
    do j=1,nn(1) ! nn(1) = L
        tt=0.d0
        do i=1,nn(0)  ! nn(0) = M : number of atoms
            tt=tt+a(i,j,1)*y(i,0)
        enddo
        x(j,1)=b(j,1)+tt
        y(j,1)=tanh(x(j,1))
        yd(j,1)=1.d0/cosh(x(j,1))**2
    enddo
    !-------------------------------------------------------
    do k=1,nn(2) ! nn(2) = P
        tt=0.d0
        do j=1,nn(1) 
            tt=tt+a(j,k,2)*y(j,1)
        enddo
        x(k,2)=b(k,2)+tt
        y(k,2)=tanh(x(k,2))
        yd(k,2)=1.d0/cosh(x(k,2))**2
    enddo
    !-------------------------------------------------------
    tt=0.d0
    do k=1,nn(2) 
        tt=tt+a(k,1,3)*y(k,2)
    enddo
    x(1,3)=b(1,3)+tt
    y(1,3)=x(1,3)
    yd(1,3)=1.d0
    !-------------------------------------------------------
    epot=y(1,3)

    do j=1,nn(1)
        tt=0.d0
        do k=1,nn(2)
            tt=tt+a(k,1,3)*yd(k,2)*a(j,k,2)
        enddo
        c(j)=tt
    enddo
    do i=1,nn(0)
        tt=0.d0
        do j=1,nn(1)
            tt=tt+yd(j,1)*a(i,j,1)*c(j)
        enddo
        d(i)=tt
    enddo

end subroutine cal_architecture_2hiddenlayer

!*****************************************************************************************
subroutine cal_architecture_der_2hiddenlayer(nl,nn,a,b,x,y,ad,bd,epot)
    implicit none
    integer, intent(in):: nl !number of hidden layer plus one
    integer, intent(in):: nn(0:nl) !number of nodes in each hiddenlayer
    real(8), intent(inout):: a(140,140,nl), b(140,nl), x(140,nl), y(140,0:nl)
    real(8), intent(inout):: ad(140*140,nl), bd(140,nl)
    real(8), intent(inout):: epot
    !local variables
    integer:: i, j, k, l, lp
    real(8):: tt
    if(nl/=3) then
        write(*,'(a,i3)') 'ERROR: this routine works only for nl=3, while nl= ',nl
        stop
    endif
    !-------------------------------------------------------
    do j=1,nn(1)
        tt=0.d0
        do i=1,nn(0)
            tt=tt+a(i,j,1)*y(i,0)
        enddo
        x(j,1)=b(j,1)+tt
        y(j,1)=tanh(x(j,1))
    enddo
    !-------------------------------------------------------
    do k=1,nn(2)
        tt=0.d0
        do j=1,nn(1)
            tt=tt+a(j,k,2)*y(j,1)
        enddo
        x(k,2)=b(k,2)+tt
        y(k,2)=tanh(x(k,2))
    enddo
    !-------------------------------------------------------
    tt=0.d0
    do k=1,nn(2)
        tt=tt+a(k,1,3)*y(k,2)
    enddo
    x(1,3)=b(1,3)+tt
    epot=x(1,3)
    !-------------------------------------------------------
    do l=1,nn(1)
        tt=0.d0
        do k=1,nn(2)
            tt=tt+a(k,1,3)*a(l,k,2)/cosh(x(k,2))**2
        enddo
        tt=tt/(cosh(x(l,1)))**2
        bd(l,1)=tt
        do lp=1,nn(0)
            ad(lp+(l-1)*nn(0),1)=tt*y(lp,0)
        enddo
    enddo
    !-------------------------------------------------------
    do l=1,nn(2)
        tt=a(l,1,3)/(cosh(x(l,2)))**2
        bd(l,2)=tt
        do lp=1,nn(1)
            ad(lp+(l-1)*nn(1),2)=tt*y(lp,1)
        enddo
    enddo
    !-------------------------------------------------------
    tt=1.d0
    bd(1,3)=tt
    do l=1,nn(2)
        ad(l,3)=tt*y(l,2)
    enddo
    !-------------------------------------------------------
end subroutine cal_architecture_der_2hiddenlayer
!!!***************************************************************************************
subroutine convert_w_ann(n,w,a,b,nl,nn)
    implicit none
    integer, intent(in):: n
    real(8), intent(in):: w(n)
    integer, intent(in):: nl, nn(0:nl)
    real(8), intent(inout):: a(140,140,nl), b(140,nl)
    
    !local variables
    integer:: i, j, l, ialpha
    l=0
    do ialpha=1,nl
        do j=1,nn(ialpha)
            do i=1,nn(ialpha-1)
                l=l+1
                a(i,j,ialpha)=w(l)
            enddo
        enddo
        do i=1,nn(ialpha)
            l=l+1
            b(i,ialpha)=w(l)
        enddo
    enddo
    if(l/=n) stop 'ERROR: l/=n'
end subroutine convert_w_ann
!!***************************************************************************************
subroutine convert_ann_epotd(nn,nl,ad,bd,n,epotd)
    implicit none
    integer, intent(in):: nl
    integer, intent(in):: nn(0:nl)
    real(8), intent(in):: ad(140*140,nl), bd(140,nl)
    integer, intent(in):: n
    real(8), intent(inout):: epotd(n)
    !local variables
    integer:: i, ij, l, ialpha
    l=0
    do ialpha=1,nl
        do ij=1,nn(ialpha)*nn(ialpha-1)
            l=l+1
            epotd(l)=ad(ij,ialpha)
        enddo
        do i=1,nn(ialpha)
            l=l+1
            epotd(l)=bd(i,ialpha)
        enddo
    enddo
    if(l/=n) stop 'ERROR: l/=n'
end subroutine convert_ann_epotd


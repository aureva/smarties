      SUBROUTINE  wallstress2_ML_cal(ML_in,ML_out,B,
     +                       W1,W2,W3,W4,W5,W6,W7,
     +                       dpdx,dpdz)

        implicit none
        include 'dimen.h'
        real*8,dimension(18,1),intent(in)   :: ML_in
        real*8,dimension(2,1),intent(out)     :: ML_out
        real*8:: B(20,7),W1(20,18),W2(20,20),W3(20,20),W4(20,20),
     +           W5(20,20),W6(20,20),W7(2,20), temp(20,1),temp1(20,1),
     +     temp2(20,1),temp3(20,1),temp4(20,1),temp5(20,1),temp6(20,1)
          
        real*8,dimension(20,1) :: B1,B2,B3,B4,B5,B6
        real*8,dimension(2,1) :: B7
        integer*4 :: i
        real*8, dimension(nx,ny,nz2):: dpdx,dpdz

        do i = 1,20
          B1(i,1) = B(i,1)
          B2(i,1) = B(i,2)
          B3(i,1) = B(i,3)
          B4(i,1) = B(i,4)
          B5(i,1) = B(i,5)
          B6(i,1) = B(i,6)
        end do
        do i = 1,2
          B7(i,1) = B(i,7)
        end do

        temp(:,:) = 0d0

        ! Calculate the result
        temp1 = max(temp, matmul(W1,ML_in)+B1)
        temp2 = max(temp, matmul(W2,temp1)+B2)
        temp3 = max(temp, matmul(W3,temp2)+B3)
        temp4 = max(temp, matmul(W4,temp3)+B4)
        temp5 = max(temp, matmul(W5,temp4)+B5)
        temp6 = max(temp, matmul(W6,temp5)+B6)

        ML_out = matmul(W7,temp6)+B7

        RETURN
      
      END 


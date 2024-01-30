      subroutine divstress_uv (divt, tx, ty, tz,me,t)
      implicit none
      include 'dimen.h'
      integer*4 i,j,k,t
      real*8,dimension(nx,ny,nz2):: divt,tx,ty,tz,tdx,tdy,tdz
c...  Compute stress gradients
      if(S_flag.eq.1)then
         call ddx(tdx,tx,t,0)
         call ddy(tdy,ty,t,0)
      else
         call ddx(tdx,tx,t,1)
         call ddy(tdy,ty,t,1)
      endif
      call ddz_w(tdz, tz,me)

      do k=2,Nzb+1
         do j=1,Ny
            do i=1,Nx              
               divt(i,j,k)=tdx(i,j,k)+ tdy(i,j,k)+ tdz(i,j,k)
            end do
         end do
      end do

      return  
      end 

      subroutine divstress_w (divt, tx, ty, tdz, me,t)

      implicit none
      include 'dimen.h'
      integer*4 i,j,k,t
      real*8,dimension(nx,ny,nz2):: divt,tx,ty,tz,tdx,tdy,tdz
      
c...  Compute stress gradients      
      call ddx(tdx,tx,t,2)
      call ddy(tdy,ty,t,2)
c...  Note ddz_uv does not return values for k=1 level.
      do k=2,Nzb+1
         do j=1,Ny
            do i=1,Nx              
               divt(i,j,k)=tdx(i,j,k)+ tdy(i,j,k)+ tdz(i,j,k)
            end do
         end do

c...  At wall we have to assume that tdz(tzz)=0.0.  Any better ideas?
         if (me==0.and.k==2)then
            do j=1,Ny
               do i=1,Nx
                  divt(i,j,k)=tdx(i,j,k)+ tdy(i,j,k)
               end do
            end do
         end if
      end do
      return  
      end 

















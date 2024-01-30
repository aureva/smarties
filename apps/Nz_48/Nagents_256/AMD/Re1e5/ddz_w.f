      subroutine ddz_w (DFDZ, F,me)

c...F is on w nodes and dFdz is on uvp nodes
     
      implicit none
      include 'dimen.h'
      integer*4 i,j,k
      real*8 F(nx,ny,nz2), dfdz(nx,ny,nz2)
c...2*pi*aspect is the nondimensional vertical depth of the domain 
c...(Aspect is in dimen.h)
      
      do k=2,Nzb+1
         do j=1,Ny
            do i=1,Nx    
               dfdz(i,j,k)=(f(i,j,k+1)-f(i,j,k))*idz
            end do
         end do   
      end do  

c...  Stress free lid
	
      if (me==nprocs-1) then
         dfdz(:,:,Nzb+1)=0.0
      end if

      return
      end

      subroutine ddz_uv (DFDZ, F,me)
c
c  First deriv in z direction for boundary layer (2nd order numerics)
c...F is on UVP nodes and dFdz is on w nodes
c     
      implicit none
      include 'dimen.h'
      integer*4 i,j,k
      real*8 F(nx,ny,nz2), dfdz(nx,ny,nz2)
c...2*pi*aspect is the nondimensional vertical depth of the domain (Aspect is in dimen.h)     
c...Skip wall level, which will be done later with M.O.
      do k=2,Nzb+1
       do j=1,Ny
          do i=1,Nx    
          dfdz(i,j,k)=(f(i,j,k)-f(i,j,k-1))*idz
          end do
       end do   
      end do  

      if(me==0) then
         do j=1,ny
            do i=1,nx
               dfdz(i,j,2)=0.
            end do
         end do
      end if
      return
      end

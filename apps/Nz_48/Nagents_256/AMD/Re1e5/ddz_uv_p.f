      subroutine ddz_uv_p (DFDZ, F,me)

c  First deriv in z direction for boundary layer (2nd order numerics)
c...F is on UVP nodes and dFdz is on w nodes

      implicit none
      include 'dimen.h'
      integer*4 i,j,k
      real*8 F(nx,ny,nz2), dfdz(nx,ny,nz2)
      real*8 zz(nz2), d_avg
          
      do k=1,Nzb+2
         zz(k)=(k-1.+me*nzb-0.5)*DZ   
      end do
         
      do k=2,Nzb+1
         d_avg=0.
         do j=1,Ny
            do i=1,Nx    
               dfdz(i,j,k)=(f(i,j,k)-f(i,j,k-1))*idz
               d_avg=d_avg+dfdz(i,j,k)*inxny
            end do
         end do

c         if((k.eq.3).and.me==0)then
c            do j=1,ny
c               do i=1,nx
c                  dfdz(i,j,k)=d_avg*dz/((zz(k-1)+0.5*dz)*
c     +                 dlog(zz(k)/zz(k-1)))+(dfdz(i,j,k)-d_avg)
c               end do
c            end do
c         end if

         
      end do

      return
      end

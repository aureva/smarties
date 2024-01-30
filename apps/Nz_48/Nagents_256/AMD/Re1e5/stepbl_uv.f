      subroutine STEPBL_uv (ui, RHSi, RHSi_f, force,me)
c...  For staggered grid where U and V are not forced to 0 at k=1
      implicit none
      include 'dimen.h'
      integer*4 i,j,k
      real*8 ui(nx,ny,nz2),RHSi(nx,ny,nz2),RHSi_f(nx,ny,nz2),force(Nz2)
      real*8 u_meanP,u_meanL,uii(nx,ny,nz2)

      do k=2,Nzb+1
         do j=1,Ny
            do i=1,Nx
               ui(i,j,k)= ui(i,j,k)+DT*(1.5*RHSi(i,j,k)-
     +              0.5*RHSi_f(i,j,k)+force(k))
            end do
         end do
      end do     
      
c...No-stress top

      if(me==nprocs-1) then
         do j=1,Ny
            do i=1,Nx
               ui(i,j,Nzb+1)=ui(i,j,Nzb)   
            end do
         end do
      end if
      
      return            
      end 

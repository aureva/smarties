      subroutine STEP_S (S, RHS, RHS_f,flag,me)
      implicit none
      include 'dimen.h'
      integer*4 i,j,k,flag
      real*8,dimension(nx,ny,nz2):: s, RHS, RHS_f

      do k=2,Nzb+1
         do j=1,Ny
            do i=1,Nx
!               s(i,j,k)= s(i,j,k)+DT*(1.5*RHS(i,j,k)-0.5*RHS_f(i,j,k))
               s(i,j,k)= s(i,j,k)+DT*(1.5*RHS(i,j,k)-0.5*RHS_f(i,j,k)
     +                  -1*s_flux*z_i/l_z)				   
            end do
         end do
      end do     

c...  No-stress top
      if (me==nprocs-1) then
         if(flag.eq.0)then
            do j=1,Ny
               do i=1,Nx
                  s(i,j,Nzb+1)=s(i,j,Nzb)+inversion*dz
!                  s(i,j,Nzb+1)=268.0317d0				  
               end do
            end do
         else
            do j=1,Ny
               do i=1,Nx
                  s(i,j,Nzb+1)=s(i,j,Nzb)+qinversion*dz
               end do
            end do
         endif

      end if

      return            
      end 

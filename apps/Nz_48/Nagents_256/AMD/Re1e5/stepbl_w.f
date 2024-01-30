      subroutine STEPBL_w (ui,RHSi,RHSi_f,me)
      implicit none
  
      include 'dimen.h'
      integer*4 i,j,k
      real*8 ui(nx,ny,nz2),RHSi(nx,ny,nz2),RHSi_f(nx,ny,nz2),
     +       w_bar,uii(nx,ny,nz2)

      do k=2,Nzb+1

          do j=1,Ny
          do i=1,Nx
          ui(i,j,k)=ui(i,j,k)+DT*(1.5*RHSi(i,j,k)-0.5*RHSi_f(i,j,k))
	      end do
          end do
	  
c	  w_bar=0.0
c          do j=1,Ny
c              do i=1,Nx
c		w_bar= w_bar+ui(i,j,k)
c              end do
c          end do
c	  w_bar=w_bar*iNxNy
	  
c........And, remove the mean vertical velocity  
c          do j=1,Ny
c              do i=1,Nx
c		  ui(i,j,k)= ui(i,j,k)-w_bar
c	      end do
c          end do

      end do     

c...No slip wall and no flow through top

      if (me==0) then
         ui(:,:,2)=0.0
      end if

      if (me==nprocs-1) then
         ui(:,:,nzb+1)=0.0
      end if
	  
      return            
      end 

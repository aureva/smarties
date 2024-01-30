CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      Subroutine inflow_ts(theta,t,me,nall)
      Subroutine inflow_ts(theta,RHS_T_BF,t,me,nall)
      implicit none
      include 'dimen.h'
      
      integer*4 i,j,k,t,bx,fin,ii,jj,sn
      real*8, dimension(nx,ny,nz2)::theta
      real*8, dimension(ny,nz)::thetau

      real*8, dimension(Nxbe-Nxbm+1,ny,nz2)::tmptheta2
      real*8, dimension(Nxbe-Nxbm+1,ny,nz2)::theta_f,theta_b	  
	  
      real*8, dimension(nx,ny,nz2):: RHS_T_BF
      real*8  factor,cfrdmp

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      call read_instant2_parallel(1,tmpu2,t,me,nall)
!      call read_instant2_parallel(4,tmptheta2,t,me,nall)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc      ii=1
      do k=2,nzb+1
      do j=1,ny
      do i=1,Nxbe-Nxbm+1	  
!      theta_f(i,j,k)=tmptheta2(i,j,k)
      theta_f(i,j,k)=293.0d0
      end do
      end do
      enddo	 
	  
      sn=0
      do j=1,ny-sn
         theta_b(:,j+sn,:)=theta_b(:,j,:)
      end do
      do j=ny-sn+1,ny 
         theta_b(:,ny-j+1,:)=theta_b(:,ny-sn+ny-j+1,:)
      end do

!	   theta_f=theta_b

!      do k=2,nzb+1
!          do j=1,Ny
!             theta(Nxbm:Nxbe,j,k) = theta_f(:,j,k)
!          end do
!      end do	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RHS_T_BF=0.0d0	  	
      cfrdmp = 1./(rlx_time*u_star/z_i)	 

      do k=2,nzb+1
         do i=Nxbs,Nxbm
            do j=1,Ny
              factor = 0.5*(1.-cos(pi*1.*(i-Nxbs)/(Nxbm-Nxbs)))
              RHS_T_BF(i,j,k)=cfrdmp*factor*
     &                       (theta(i,j,k)-theta_f(1,j,k))			  
!            theta(i,j,k)=theta(Nxbs,j,k)+factor*
!     + 	                (theta(Nxbm,j,k)-theta(Nxbs,j,k))		
            end do
          end do
      end do
  
      return

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



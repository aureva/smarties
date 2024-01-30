CCCC NOTE: current version using linearized version of CCCC
CCCC virtual temp (see Garratt 92 page 22).            CCCC
      Subroutine CalcBeta (theta,q,b,me)
      implicit none
      include 'dimen.h'
      integer*4 i, j, k
      real*8,dimension(Nx,Ny,Nz2) :: b,theta,q
      real*8 theta_bar(Nz2),q_bar(Nz2),above, below
      
c...  Note Beta is stored on W nodes, but Theta is on UVP nodes
      
      theta_bar=0.d0
      do k=1,nzb+1
         do j=1,ny
            do i=1,nx
               theta_bar(k)=theta_bar(k)+theta(i,j,k)
            enddo
         enddo
         theta_bar(k)=theta_bar(k)*inxny
      enddo
      
      q_bar=0.d0
      if(q_flag.eq.1.and.qpass_flag.eq.0)then
         do k=1,nzb+1
            do j=1,ny
               do i=1,nx
                  q_bar(k)=q_bar(k)+q(i,j,k)
               enddo
            enddo
            q_bar(k)=q_bar(k)*inxny
         enddo
      endif
      
      do k=2,nzb+1
         if (Theta_bar(k).eq.(0.d0))then
            b(:,:,k) = 0.
         elseif(q_bar(k).eq.(0.d0).or.qpass_flag.eq.1)then
            do j=1,Ny
               do i=1,nx		    
                  above=(Theta(i,j,k)-Theta_bar(k))/(theta_0)
                  below=(Theta(i,j,k-1)-Theta_bar(k-1))/(theta_0)
                  b(i,j,k)=g_hat*(above + below)*0.5
               end do
            end do
         else
            do j=1,Ny
               do i=1,nx      		    
                  above=(Theta(i,j,k)-Theta_bar(k) + 
     +                 0.61d0*theta_0*(q(i,j,k)-q_bar(k)))/(theta_0)
                  below=(Theta(i,j,k-1)-Theta_bar(k-1) +
     +                 0.61d0*theta_0*(q(i,j,k-1)-q_bar(k-1))
     +                 )/(theta_0)
                  b(i,j,k)=g_hat*(above + below)*0.5
               end do
            end do
         endif
      end do
	
      if (me==0) then
         
         do j=1,Ny
            do i=1,nx
               b(i,j,2)=0
            enddo
         enddo
	   
      endif
      
      return
      
      end
				


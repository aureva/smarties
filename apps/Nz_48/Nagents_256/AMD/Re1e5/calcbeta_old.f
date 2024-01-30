	Subroutine CalcBeta (theta, b,me)
 	implicit none
 	include 'dimen.h'
 	integer*4 i, j, k
 	real*8 b(Nx,Ny,Nz2), theta(Nx,Ny,Nz2), theta_bar(Nz2),
     +	     above, below
	
c...    Note Beta is stored on W nodes, but Theta is on UVP nodes
	
	do k=1,nzb+1
	   theta_bar(k)=0.0
	   do j=1,ny
	      do i=1,nx
		 theta_bar(k)=theta_bar(k)+theta(i,j,k)
	      enddo
	   enddo
	   theta_bar(k)=theta_bar(k)*inxny
	enddo

	do k=2,nzb+1
	   do j=1,Ny
	      do i=1,nx
		 if (Theta_bar(k).eq.0) then
		    b(i,j,k) = 0.
		 else		    
		    above=(Theta(i,j,k)-Theta_bar(k))/(theta_0)
		    below=(Theta(i,j,k-1)-Theta_bar(k-1))/(theta_0)
		    b(i,j,k)=g_hat*(above + below)/2.
		 endif
	      end do
	   end do
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

				


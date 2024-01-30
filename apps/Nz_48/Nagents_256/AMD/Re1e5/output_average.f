ccc NOTE ave_x is expected to be of size anx,nz2 not nz ccc
ccc This routine does a reduce to me-hfact.eq.0 and sum ccc
ccc and then a reduce to me=0 for final stat output.    ccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine output_average(ave_x,nn,fnum,me,nall)

      implicit none
      include 'dimen.h'
      integer*4 fnum,ii,nn
      
      real*8,dimension(nn,nz2) :: ave_x,bar
      real*8,dimension(nn,nz) :: x_bar
         
ccc now send the sum to me=0 for output
         if(me.gt.0)then
            call MPI_SEND(ave_x(1,2),nn*nzb,MPI_DOUBLE_PRECISION,
     +           0,me,nall,ierr )
         else
            x_bar(:,1:nzb)=ave_x(:,2:nzb+1)
            do ii=1,nprocs-1
               call MPI_RECV(x_bar(1,ii*nzb+1),nzb*nn,
     +              MPI_DOUBLE_PRECISION,ii,ii,
     +              nall,status2,ierr)
            enddo
            
            write(fnum) x_bar
            call flush(fnum)

         endif

      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine BL_height(u,v,w,txz,tyz,t_total,t_t,me)
	  
      implicit none
      include 'dimen.h'
      integer*4 i,j,k,k_bl
      real*8,dimension(nx,ny,nz2):: u,v,w,txz,tyz,uw,vw
      real*8,dimension(nz2):: txz_total,tyz_total,t_total,txz_t,tyz_t
      real*8,dimension(nz2)::u_bar,v_bar,w_bar,t_t
      real*8 arg1,arg2,norm
	  
      norm =1.d0/(Ny*Nx)
      txz_total=0.0
      tyz_total=0.0
      t_total=0.0

ccc compute the plane averages of terms involved in products ccc
      txz_t=0.
      tyz_t=0.
      t_t=0.
      u_bar=0.d0
      v_bar=0.d0
      w_bar=0.d0
      do k=1,Nzb+1
         do j=1,Ny
            do i=1,Nx
               u_bar(k)=u_bar(k)+u(i,j,k)
               v_bar(k)=v_bar(k)+v(i,j,k)
            enddo
         enddo             
      enddo
      do k=2,Nzb+1
         do j=1,Ny
            do i=1,Nx
               w_bar(k)=w_bar(k)+w(i,j,k)
            enddo
         enddo             
      enddo
      u_bar=u_bar*inxny
      v_bar=v_bar*inxny
      w_bar=w_bar*inxny

      do k=2,Nzb+1
      do j=1,Ny
      do i=1,Nx	  
         if (k==2.and.me==0) then
           arg1=0.
           arg2=0.
         else
         arg1=0.5d0*(u(i,j,k)+u(i,j,k-1) - 
     +        u_bar(k) - u_bar(k-1))
         arg2=0.5d0*(v(i,j,k)+v(i,j,k-1) -
     +        v_bar(k) - v_bar(k-1))
         end if
               
         txz_t(k)=txz_t(k)+(w(i,j,k)-w_bar(k))*arg1+txz(i,j,k)
         tyz_t(k)=tyz_t(k)+(w(i,j,k)-w_bar(k))*arg2+tyz(i,j,k)
      enddo
      enddo
      enddo
	 
      txz_t=txz_t*inxny
      tyz_t=tyz_t*inxny

      t_t=(txz_t**2+tyz_t**2)**0.5	  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  
      do k=2,Nzb+1
      do j=1,Ny
      do i=1,Nx
      if (k==2.and.me==0) then
         arg1=0.
         arg2=0.
      else
         arg1=(u(i,j,k)+u(i,j,k-1))/2.   
         arg2=(v(i,j,k)+v(i,j,k-1))/2.
      end if
      uw(i,j,k)=w(i,j,k)*arg1  
      vw(i,j,k)=w(i,j,k)*arg2  
      end do
      end do
      end do	  
	  
      do k=2,Nzb+1
         do i=1,Nx
            do j=1,Ny	  
               txz_total(k)=txz_total(k)+(uw(i,j,k)+txz(i,j,k))
               tyz_total(k)=tyz_total(k)+(vw(i,j,k)+tyz(i,j,k))
            enddo
         enddo
      enddo
	  
      txz_total=txz_total*inxny
      tyz_total=tyz_total*inxny
	  
      t_total=(txz_total**2+tyz_total**2)**0.5
c      if (me==0) print*,'t_total(2)',t_total(3)**0.5,t_t(3)**0.5

      return
      end







	  
	  
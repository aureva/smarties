      subroutine rmsdiv(rms, du, dv, dw,me)
      implicit none
      include 'dimen.h'
      integer*4 i,j,k,kend
      real*8 du(nx,ny,nz2), dv(nx,ny,nz2), dw(nx,ny,nz2)
      real*8 rms
      
      if (me==nprocs-1) then
         kend=nzb
      else
         kend=nzb+1
      end if
      rms=0.0

      do k=2,kend
         do j=1,Ny
            do i=1,Nx
               rms=rms+abs(du(i,j,k)+dv(i,j,k)+dw(i,j,k))
            end do
         end do
      end do
      rms=rms/(nx*ny*nz)
      return
      end

	  
      subroutine Check_CFL(u,v,w,CFLx,CFLy,CFLz,me,nall)
      implicit none
      include 'dimen.h'
      integer*4 i,j,k,kend
      real*8 u(nx,ny,nz2),v(nx,ny,nz2),w(nx,ny,nz2)
      real*8 CFLx,CFLy,CFLz
               
      CFLx = maxval(abs(u))*dt/dx
      CFLy = maxval(abs(v))*dt/dy
      CFLz = maxval(abs(w))*dt/dz

      return
      end
	   
      subroutine derivwall2_ML (dtdz,dqdz,dudz,dvdz,u_hat,v_hat,fi,
     +     fi_h,t_flux,q_flux,ustar,M,t)
    
      implicit none
      include 'dimen.h'
      integer*4 i,j,t
      real*8,dimension(nx,ny,nz2):: dudz,dvdz,dtdz,dqdz
      real*8,dimension(nx,ny,2):: fi,fi_h
      real*8,dimension(nx,ny):: ustar,M,t_flux,q_flux,u_hat,v_hat

!c     Option 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !do j=1,ny
      !   do i=1,nx         
!            dvdz(i,j,2)=+0.
!     +           +0.5*(v(i,j,3)+v(i,j,2)-0.)/dz
      !      dvdz(i,j,2)=2.*dvdz(i,j,3)*fi(i,j,1)/fi(i,j,2)
      !   enddo
      !enddo
      do j=1,ny
         do i=1,nx
            !dudz(i,j,2)=2.*dudz(i,j,3)*fi(i,j,1)/fi(i,j,2)
            dtdz(i,j,2)=2.*dtdz(i,j,3)*fi_h(i,j,1)/fi_h(i,j,2)            
         end do
      end do
!c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
!c     Option 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       do j=1,ny
!          do i=1,nx
!             dudz(i,j,2)=fi(i,j,1)*ustar(i,j)*(u_hat(i,j)+Ugal)/
!     +            (M(i,j)*vonk*0.5*dz)
!             dvdz(i,j,2)=fi(i,j,1)*ustar(i,j)*(v_hat(i,j)+Vgal)/
!     +           (M(i,j)*vonk*0.5*dz)
!             if(s_flag.eq.1)then
!             dtdz(i,j,2)=fi_h(i,j,1)*(-t_flux(i,j)/ustar(i,j))/
!     +              (vonk*0.5*dz) 
!             endif
!             if(q_flag.eq.1)then
!             dqdz(i,j,2)=fi_h(i,j,1)*(-q_flux(i,j)/ustar(i,j))/
!     +              (vonk*0.5*dz)
!             endif
!           end do
!       end do
!c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      return
      end

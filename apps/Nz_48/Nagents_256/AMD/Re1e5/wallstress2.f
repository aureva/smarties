      subroutine wallstress2(txz,tyz,u_hat,v_hat,psi,psi0,zo,
     +   ustar,M,t,tauw,flag_RL)

      implicit none
      include 'dimen.h'
      integer*4 i,j,nx_decorr,t,flag_RL
      
      real*8,dimension(nx,ny,nz2) :: txz,tyz
      real*8,dimension(nx,ny) :: U_res,u_hat,v_hat,ustar,M,zo,tauw
      
      real*8 denom,ang_deg,factor,u_used,v_used,u_res_sum,psi(nx,ny,2),
     +   psi0(nx,ny,2),ustar_avg,U_mean,V_mean,tauw_avg
      
!c.......angle to account for inclined structure in velocity/shear correlation
!c     ang_deg=13.
!c     factor=(dz/2.)/(tan(ang_deg*Pi/180.)*(2.*Pi/Nx2))
      factor=0
      nx_decorr=0
      ustar_avg = 0d0
      tauw_avg = 0d0
      
      do j=1,Ny
         do i=1,Nx
            
            if(i.lt.Nx)then
               u_used=factor*(u_hat(i+1,j)+Ugal)+(1.-factor)*
     +              (u_hat(i,j)+Ugal)
               v_used=factor*(v_hat(i+1,j)+Vgal)+(1.-factor)*
     +              (v_hat(i,j)+Vgal)
            else
               u_used=factor*(u_hat(1,j)+Ugal)+(1.-factor)*
     +              (u_hat(Nx,j)+Ugal)
               v_used=factor*(v_hat(1,j)+Vgal)+(1.-factor)*
     +              (v_hat(Nx,j)+Vgal)
            end if
            
            denom=dlog((dz/2d0)/(zo(i,j)/z_i))+Psi(i,j,1)-Psi0(i,j,1)

            if (flag_RL == 1) then
              ustar(i,j)= dsqrt(dabs(tauw(i,j)))
              txz(i,j,2)=-tauw(i,j)*u_used/M(i,j)
              tyz(i,j,2)=-tauw(i,j)*v_used/M(i,j)
            else
              ustar(i,j)= M(i,j)*vonk/denom
              txz(i,j,2)=-ustar(i,j)**2*u_used/M(i,j)
              tyz(i,j,2)=-ustar(i,j)**2*v_used/M(i,j)
            end if
 
            tauw_avg = tauw_avg + tauw(i,j)
            ustar_avg=ustar_avg+ustar(i,j)   
                    
         enddo
      enddo
      tauw_avg=tauw_avg*inxny
      ustar_avg=ustar_avg*inxny
  
!      print *, 'tauw,ustar',tauw_avg,ustar_avg,txz(10,10,2),tyz(10,10,2)

!c     MKP model
!c      U_mean=(sum(u_hat))*inxny
!c      V_mean=(sum(v_hat))*inxny
!c      do j=1,Ny
!c         do i=1,Nx         
!c            txz(i,j,2)=-ustar_avg**2.-0.1*ustar_avg*(u_hat(i,j)-U_mean)
!c            tyz(i,j,2)=-0.1*ustar_avg*(v_hat(i,j)-V_mean)          
!c         enddo
!c      enddo			

!            if (mod(t,c_count)==0) then
!              write(1525,5111) ustar_avg
!              write(*,*) 'u*',ustar_avg  
!            endif

 5111 format (5(1x,f11.5))  
 
      return    
      
      end



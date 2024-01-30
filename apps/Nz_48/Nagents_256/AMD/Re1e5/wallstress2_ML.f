      subroutine wallstress2_ML(txz,tyz,u,v,w,psi,psi0,zo,
     +   ustar,M,t,B,W1,W2,W3,W4,W5,W6,W7,max_in,min_in,
     +   dpdx,dpdz)

      implicit none
      include 'dimen.h'
      integer*4 i,j,k,ind,nx_decorr,t,count
      
      real*8,dimension(nx,ny,nz2) :: txz,tyz,u,v,w
      real*8,dimension(nx,ny) :: U_res,ustar,M,zo
      
      real*8 :: denom,ang_deg,factor,u_used,v_used,u_res_sum,
     + psi(nx,ny,2),psi0(nx,ny,2),rho,dist,gradpx,up,uv,u_taup,y_star,
     + ustar_avg,hwm,h,u_used1,v_used1,ub

      real*8,dimension(18,1)   :: ML_in
      real*8,dimension(2,1)   :: ML_out
      real*8,dimension(20)   :: max_in, min_in
      real*8, dimension(20,7) :: B
      real*8, dimension(20,18) :: W1
      real*8, dimension(20,20) :: W2,W3,W4,
     +                            W5, W6
      real*8, dimension(2,20) :: W7
      real*8, dimension(nx,ny,nz2):: dpdx,dpdz
      
!c.......angle to account for inclined structure in velocity/shear correlation
!c     ang_deg=13.
!c     factor=(dz/2.)/(tan(ang_deg*Pi/180.)*(2.*Pi/Nx2))
      factor=0
      nx_decorr=0
      ustar_avg = 0.0

      !dist = 0.03d0/4.5d0*Pi/l_r*z_i
      gradpx = 1d0
      rho = 1d0
      up = dabs(nu/rho*gradpx)**(1d0/3d0)
!c    h = z_i/3.036
      h = L_z
                   
      do j=1,Ny
         do i=1,Nx
           
            do k = 1,3

              hwm = (dz/2d0 + (k-1)*dz)/h
 
              if(i.lt.Nx)then
                 u_used=factor*(u(i+1,j,k+1)+Ugal)+(1.-factor)*
     +                (u(i,j,k+1)+Ugal)
                 v_used=factor*(v(i+1,j,k+1)+Vgal)+(1.-factor)*
     +                (v(i,j,k+1)+Vgal)
              else
                 u_used=factor*(u(1,j,k+1)+Ugal)+(1.-factor)*
     +                (u(Nx,j,k+1)+Ugal)
                 v_used=factor*(v(1,j,k+1)+Vgal)+(1.-factor)*
     +                (v(Nx,j,k+1)+Vgal)
              end if
 
              if (k==1) then
                u_used1 = u_used
                v_used1 = v_used
              end if
              uv = dsqrt(dabs(nu*u_used/hwm))
              u_taup = dsqrt(uv**2 + up**2)
              y_star = nu/u_taup/h

              ub = 19.999d0 ! u_star/vonk*log(u_star*l_z/nu) + 5d0
 
 
              ind = 6*(k-1)

              ML_in(ind+1,1) = (dlog(hwm/y_star)-min_in(ind+1))
     +                            /(max_in(ind+1)-min_in(ind+1))  
              ML_in(ind+2,1) = (u_used/ub/hwm-min_in(ind+2))
     +                            /(max_in(ind+2)-min_in(ind+2))  
              ML_in(ind+3,1) = (w(i,j,k+1)/ub/hwm-min_in(ind+3))
     +                            /(max_in(ind+3)-min_in(ind+3))  
              ML_in(ind+4,1) = (v_used/ub/hwm-min_in(ind+4))
     +                            /(max_in(ind+4)-min_in(ind+4))  
              ML_in(ind+5,1) = (dpdx(i,j,k+1)*h/ub**2*hwm/h
     +                           -min_in(ind+5))
     +                            /(max_in(ind+5)-min_in(ind+5))  
              ML_in(ind+6,1) = (dpdz(i,j,k+1)*h/ub**2
     +                           *hwm/h-min_in(ind+6))
     +                            /(max_in(ind+6)-min_in(ind+6))  

            end do

            CALL wallstress2_ML_cal(ML_in,ML_out,B,
     +                       W1,W2,W3,W4,W5,W6,W7)


            tyz(i,j,2) = (ML_out(1,1)*(max_in(19)-min_in(19))
     +                   + min_in(19))*ub**2
            txz(i,j,2) = (ML_out(2,1)*(max_in(20)-min_in(20))
     +                   + min_in(20))*ub**2
            ustar(i,j) = dsqrt(dsqrt(txz(i,j,2)**2+tyz(i,j,2)**2)
     +                       *M(i,j)/(u_used1**2 + v_used1**2))

            if (1==0) then
            print *,'ub',ub,'i',i,',','j',j,',','in 1,1',ML_in(1,1),',',
     +             'in 2,1',ML_in(2,1),',','in 3,1',ML_in(3,1),',',
     +             'in 4,1',ML_in(4,1),',','in 5,1',ML_in(5,1),',',
     +             'in 6,1',ML_in(6,1),',', 'in 7,1',ML_in(7,1),',',
     +             'in 8,1',ML_in(8,1),',', 'in 9,1',ML_in(9,1),',',
     +             'in 10,1',ML_in(10,1),',',
     +             'in 11,1',ML_in(11,1),',','in 12,1',ML_in(12,1),',',
     +             'in 13,1',ML_in(13,1),',','in 14,1',ML_in(14,1),',',
     +             'in 15,1',ML_in(15,1),',', 'in 16,1',ML_in(16,1),',',
     +             'in 17,1',ML_in(17,1),',', 'in 18,1',ML_in(18,1),',',
     +             'txz XRWM',txz(i,j,2),',','txz ZWM',tyz(i,j,2),',',
     +             'tyz XRWM', tyz(i,j,2),',','tyz ZWM',txz(i,j,2),',',
     +                   'max out 1',max_in(19),',',
     +                   'min out 1',min_in(19),',',
     +                   'max out 2',max_in(20),',',
     +                   'min out 2',min_in(20),',',
     +                   'ML_out 1', ML_out(1,1),',',
     +                   'ML_out 2', ML_out(2,1),',',
     +                   'ML_out 1 first norm',
     +                   (ML_out(1,1)*(max_in(19)-min_in(19))
     +                   + min_in(19)),
     +                   'ML_out 2 first norm',
     +                   (ML_out(2,1)*(max_in(20)-min_in(20))
     +                   + min_in(20))

            end if


         enddo
      enddo

            if (mod(t,c_count)==0) then
              write(1525,5111) ustar_avg
              write(*,*) 'u*',ustar_avg  
            endif

 5111 format (5(1x,f11.5))  
 
      return    
      
      end



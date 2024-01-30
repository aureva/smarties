      subroutine surf_flux(theta,q,u,v,u_hat,v_hat,t_flux,q_flux,Psi,
     +     Psi0,fi,fiH,zo,t_s,q_s,coolrate,qcoolrate,ilow,U_res,t)

      implicit none
      include 'dimen.h'
      integer*4 i,j,k,kk,t
      integer*4 iii,ddd,swt,ssi,ilow(nx),il

      real*8, dimension(nx,ny,nz2):: u,v,theta,q
      real*8, dimension(nx,ny,2) :: Psi,PsiH,Psi0,PsiH0,fi,fiH 
      real*8, dimension(nx,ny) :: zo,t_s,ustar,q_s,coolrate,qcoolrate,
     +     t_flux,OB_L,u_res,t_res,u_hat,v_hat,t_sn,clrtn,qclrtn,zon,
     +     q_sn,q_res,q_flux,t_hat
      real*8 t_res_sum,u_res_sum,denom,t_flux_avg,z,x,y,denomH,
     +     OBL_avg,ustar_avg,t_bar,q_bar,q_flux_avg
      real*8 mhf,ssr,neg
      save swt

      if(t.eq.1)then
         mhf = ceiling(dx/(dt*Ugal)/2)
         swt = nint(10*dmod((mhf*dt*Ugal),dx)/dx)
         write(*,*) 'switching every -- ',swt
      endif

      call Filter_2dsl(u_hat,u(:,:,2),t,1)
      if(S_flag.eq.0)then
         call Filter_2dsl(v_hat,v(:,:,2),t,2)
      else
         call Filter_2dsl(v_hat,v(:,:,2),t,0)
      endif

!      u_hat=u(:,:,2)
!      v_hat=v(:,:,2)

!       Option: Heterogeneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      print*, Vgal,Ugal
      do j=1,Ny
         do i=1,Nx
         U_res(i,j)=((u_hat(i,j)+Ugal)**2.+(v_hat(i,j)+Vgal)**2.)**0.5
         enddo
      enddo
!       Option: Homogeneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!c$$$         U_res_sum = 0.
!c$$$         do j=1,Ny
!c$$$            do i=1,Nx
!c$$$               U_res_sum = U_res_sum + ((u_hat(i,j)+Ugal)**2.+
!c$$$     +              (v_hat(i,j)+Vgal)**2.)**0.5
!c$$$            end do
!c$$$         end do
!c$$$         U_res(:,:)=U_res_sum*iNxNy

      if(s_flag.eq.0)then
         
         do k=1,2
            Psi(:,:,k)=0.
            fi(:,:,k)=1.
            PsiH(:,:,k)=0.
            fiH(:,:,k)=1. 
         enddo
         
      else

         if(surf_flag.eq.0)then

!          t_flux=0.05*sin((t-0.5)*dt*(z_i/u_star)/3600.0d0
!     +          /24.0*2*Pi-Pi/6.0)+0.025d0 
            t_flux=s_flux 
            if(qsurf_flag.eq.0) q_flux=q_flux*10**(-4.)

            t_bar = 0.0
            do j=1,Ny
              do i=1,Nx
                t_res(i,j) = theta(i,j,2)
                t_bar = t_bar+t_res(i,j)
              end do
            end do
            t_bar = t_bar*inxny
 
!            t_s(:,:) = 3.d0*theta(:,:,2)-2.d0*theta(:,:,3)
            do j=1,Ny
            do i=1,Nx
              denom = dlog((dz/2.)/(zo(i,j)/z_i))+Psi(i,j,1)-Psi0(i,j,1)
              ustar(i,j) = (u_res(i,j)*vonk/denom)
              denomH = dlog((dz/2.)/(zo(i,j)/z_i))
     +                 +PsiH(i,j,1)-PsiH0(i,j,1)
              t_s(i,j) = t_flux(i,j)/(ustar(i,j)*vonk/denomH)+t_res(i,j)
            end do
            end do			
		
         else

            ssr = floor((dble(t+nrsub)*dt*Ugal)/dx)
            ssi = int(ssr)
            iii = mod(ssi,nx)
            ddd = nint(10*dmod((dble(t+nrsub)*dt*Ugal),dx)/dx)
            
            do i=1,nx
               il=i-iii
               if(il.lt.1) il=il+nx
               if(ddd.eq.swt)then
                  if(i.eq.nx)then
                     t_sn(i,:)=t_s(1,:)
                     clrtn(i,:)=coolrate(1,:)
                     zon(i,:)=zo(1,:)
                     if(qsurf_flag.eq.1.and.q_flag.eq.1) then
                          q_sn(i,:)=q_s(1,:)					 
                          qclrtn(i,:)=qcoolrate(1,:)
					 endif
                  else
                     t_sn(i,:)=t_s(i+1,:)
                     clrtn(i,:)=coolrate(i+1,:)
                     zon(i,:)=zo(i+1,:)
                     if(qsurf_flag.eq.1.and.q_flag.eq.1) then
                          q_sn(i,:)=q_s(i+1,:)					 
                          qclrtn(i,:)=qcoolrate(i+1,:)
					 endif
                  endif
               else
                  t_sn(i,:)=t_s(i,:)
                  clrtn(i,:)=coolrate(i,:)
                  zon(i,:)=zo(i,:)
                  if(qsurf_flag.eq.1.and.q_flag.eq.1) then
                       q_sn(i,:)=q_s(i,:)				  
                       qclrtn(i,:)=qcoolrate(i,:)
				  endif
               endif
               
               ilow(i)=il
               
            end do
            
            if(t+nrsub.gt.nint(0.*3600.*u_star/z_i/dt))then
               do j=1,ny
                  do i=1,nx
                     t_s(i,j)=t_s(i,j)-(c_coolrate/T_scale/3600)*
     +                    (z_i/u_star)*dt
                  end do
               end do
            else
              do j=1,ny
                  do i=1,nx
                     t_s(i,j)=t_sn(i,j)-(clrtn(i,j)/T_scale/3600)*
     +                    (z_i/u_star)*dt
                  end do
               end do
            endif
            
            if(q_flag.eq.1)then

               if(qsurf_flag.eq.1)then
                  if(t+nrsub.gt.nint(0.*3600.*u_star/z_i/dt))then
                     do j=1,ny
                        do i=1,nx
                           q_s(i,j)=q_s(i,j)-(cq_coolrate/Q_scale
     +                          /3600)*(z_i/u_star)*dt
                        end do
                     end do
                  else
                     do j=1,ny
                        do i=1,nx
                           q_s(i,j)=q_sn(i,j)-(qclrtn(i,j)/Q_scale
     +                          /3600)*(z_i/u_star)*dt
                        end do
                     end do
                  endif
                  qcoolrate = qclrtn
               else
                  q_flux=q_flux*10**(-4.)
               endif

            endif
            
            coolrate = clrtn
            zo       = zon
            
         endif
         
!c       Option: Homogeneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!c$$$         t_res_sum = 0.
!c$$$         do j=1,Ny
!c$$$            do i=1,Nx
!c$$$               t_res_sum = t_res_sum + theta(i,j,2)
!c$$$            end do
!c$$$         end do
!c$$$         t_res(:,:)=t_res_sum*iNxNy
!c$$$         if(q_flag.eq.1)then
!c$$$         q_res_sum = 0.
!c$$$         do j=1,Ny
!c$$$            do i=1,Nx
!c$$$               q_res_sum = q_res_sum + q(i,j,2)
!c$$$            end do
!c$$$         end do
!c$$$         q_res(:,:)=q_res_sum*iNxNy

!c       Option: Heterogeneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if(q_flag.eq.1)then
            call Filter_2dsl(t_res,theta(:,:,2),t,0)
         else
            call Filter_2dsl(t_res,theta(:,:,2),t,2)
         endif
!c$$$         t_res=theta(:,:,2)
         t_bar = 0.0
         do j=1,Ny
            do i=1,Nx
               t_bar      = t_bar+theta(i,j,2)
            end do
         end do
         t_bar = t_bar*inxny
         if(q_flag.eq.1)then
            call Filter_2dsl(q_res,q(:,:,2),t,2)
!c$$$            q_res=q(:,:,2)
            q_bar = 0.0
            do j=1,Ny
               do i=1,Nx
                  q_bar = q_bar+q(i,j,2)
               end do
            end do
            q_bar = q_bar*inxny
         endif
!c       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do kk=1,2
         
           t_flux_avg = 0.0
           q_flux_avg = 0.0
           ustar_avg = 0.0
           do j=1,Ny
              do i=1,Nx
                 denom = dlog((dz/2.*1.)/(zo(i,j)/z_i))+Psi(i,j,1)
     +                 -Psi0(i,j,1)
                 ustar(i,j) = (u_res(i,j)*vonk/denom)
                 if(surf_flag.eq.1)then
                    denomH = dlog((dz/2.*1.)/(zo(i,j)/z_i))+
     +                   PsiH(i,j,1)-PsiH0(i,j,1)
                    t_flux(i,j) = (t_s(i,j)-t_res(i,j))*ustar(i,j)*
     +                   vonk/denomH
!                    t_flux(i,j) = 0.0
                 endif
                 if(qsurf_flag.eq.1)then
                    q_flux(i,j) = (q_s(i,j)-q_res(i,j))*ustar(i,j)*
     +                   vonk/denomH
                 endif
                 q_flux_avg = q_flux_avg+q_flux(i,j)
                 t_flux_avg = t_flux_avg+t_flux(i,j)
                 ustar_avg=ustar_avg+ustar(i,j)
              end do
           end do

           t_flux_avg=t_flux_avg*inxny
           ustar_avg=ustar_avg*inxny
           q_flux_avg=q_flux_avg*inxny

           do j=1,ny
              do i=1,nx
!c$$$                 OB_L(i,j)=-ustar(i,j)**3.*t_res(i,j)*
!c$$$     +                (1.0/(vonk*g_hat*t_flux_avg))
                 OB_L(i,j)=-ustar(i,j)**3.*theta_0*
     +                (1.0/(vonk*g_hat*t_flux(i,j)))
!cc                 if (OB_L(i,j)*z_i<5.0) then
!cc                     OB_L(i,j)=5.0/z_i
!cc					 neg=neg+1.0
!cc                 endif
					 
                 do k=1,2
                    if(pass_flag.eq.0)then
                       z=+k*0.5*dz

                       if(z/OB_L(i,j).gt.5.)then
                          OB_L(i,j)=z/5.
                       elseif(z/OB_L(i,j).lt.-5.)then
                        OB_L(i,j)=-z/5.
                       endif
					 
                       if ((t_flux(i,j)).gt.0.) then
                          x=+(1.-(15.*z/OB_L(i,j)))**0.25
                          y=+(1.-(15.*(zo(i,j)/z_i)/OB_L(i,j)))**0.25
                          Psi(i,j,k)=-2.*dlog(0.5*(1+x))-
     +                         dlog(0.5*(1+x**2.))+2.*atan(x)-Pi/2.
                          Psi0(i,j,k)=-2.*dlog(0.5*(1+y))-
     +                       dlog(0.5*(1+y**2.))+2.*atan(y)-Pi/2.
                          fi(i,j,k)=1./x
!c     Equation 11.9 and 11.14 from Arya
                          psiH(i,j,k)=-2*dlog(0.5*(1+x**2.))
                          psiH0(i,j,k)=-2*dlog(0.5*(1+y**2.))
                          fiH(i,j,k)=fi(i,j,k)**2.0
                       else if ((t_flux(i,j)).lt.(-0.)) then
                          Psi(i,j,k)=+4.8*z/OB_L(i,j)
                          Psi0(i,j,k)=+4.8*(zo(i,j)/z_i)/OB_L(i,j)
                          fi(i,j,k)=1.+psi(i,j,k)
                          PsiH(i,j,k)=7.8*z/OB_L(i,j)
                          PsiH0(i,j,k)=7.8*(zo(i,j)/z_i)/OB_L(i,j)
                          fiH(i,j,k)=1.+psiH(i,j,k)
                       else
                          Psi(i,j,k)=0.
                          Psi0(i,j,k)=0.
                          fi(i,j,k)=1.
                          PsiH(i,j,k) = Psi(i,j,k)
                          PsiH0(i,j,k)= Psi0(i,j,k) 
                          fiH(i,j,k) = fi(i,j,k)
                       end if
                       
                    else
                       
                       Psi(i,j,k)=0.
                       Psi0(i,j,k)=0.
                       fi(i,j,k)=1.
                       PsiH(i,j,k)=0.
                       PsiH0(i,j,k)=0.
                       fiH(i,j,k)=1. 
!                       fiH(i,j,k)=0.74    
                       
                    endif
                 end do
               
              enddo
           enddo
           
        end do

        OBL_avg=0.0
        do j=1,Ny
           do i=1,Nx
              OBL_avg=OBL_avg+OB_L(i,j)
           enddo
        enddo

!c        write(1525,*) OBL_avg*z_i*inxny,ustar_avg,q_flux_avg
!c        write (1999,*) -t_flux_avg/ustar_avg,-q_flux_avg/ustar_avg
!c        write(*,*) 'tflux,u*,qflux',t_flux_avg,ustar_avg,q_flux_avg         

         if (mod(t,c_count)==0) then
!cc           write(1525,5111) ustar_avg,q_flux_avg,t_flux_avg,t_bar
           write (1999,5111) -t_flux_avg/ustar_avg,t_flux_avg,t_bar
     +           ,OBL_avg*z_i*inxny,ustar_avg !,neg
	  
          write(*,5555) 'tflux,u*,qflux',t_flux_avg,ustar_avg,q_flux_avg  
         endif
		 
      endif

 5111 format (6(1x,f15.5))  	  
 5555 format(1(1x,a14),3(1x,f15.8)) 
      
      return
      
      end

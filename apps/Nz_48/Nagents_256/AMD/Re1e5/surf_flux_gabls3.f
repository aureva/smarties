      subroutine surf_flux_gabls3(theta,q,u,v,u_hat,v_hat,t_flux,q_flux,
     +     zo,t_s2,q_s2,U_res,t,Psi,Psi0,fi,fiH)
                                  
      implicit none
      include 'dimen.h'
      integer*4 i,j,k,kk,t,cnt1,cnt2,cnt3
      integer*4 hourL

      real*8, dimension(nx,ny,nz2):: u,v,theta,q
      real*8, dimension(nx,ny,2) :: Psi,PsiH1,Psi0,PsiH2,fi,fiH 
      real*8, dimension(nx,ny) :: zo,ustar,
     +     t_flux,OB_L,u_res,t_res,u_hat,v_hat,t_sn,clrtn,qclrtn,zon,
     +     q_sn,q_res,q_flux
      real*8, dimension(nhrs) :: q_s2,t_s2

      real*8 t_res_sum,u_res_sum,denom,t_flux_avg,z,x,y,y2,denomH,
     +     OBL_avg,ustar_avg,q_flux_avg,OBL_max,OBL_min,q_res_sum
      real*8 Ptime,Pweig,q_s2used,t_s2used

      Ptime = dble(t+nrsub)*dt*z_i/u_star/3600.d0
      hourL = floor(Ptime)+1
      Pweig = dmod(Ptime,1.d0)
      if(Ptime.lt.9.0)then
         q_s2used = (1-Pweig)*q_s2(hourL)+Pweig*q_s2(hourL+1)
         t_s2used = (1-Pweig)*t_s2(hourL)+Pweig*t_s2(hourL+1)
      else
         q_s2used = q_s2(hourL)
         t_s2used = t_s2(hourL)
      endif
      
      call Filter_2dsl(u_hat,u(:,:,2),t,1)
      if(S_flag.eq.0)then
         call Filter_2dsl(v_hat,v(:,:,2),t,2)
      else
         call Filter_2dsl(v_hat,v(:,:,2),t,0)
      endif
      
c      u_hat=u(:,:,2)
c      v_hat=v(:,:,2)
c     Option: Heterogeneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      do j=1,Ny
         do i=1,Nx
            U_res(i,j)=((u_hat(i,j)+Ugal)**2.+
     +           (v_hat(i,j)+Vgal)**2.)**0.5
         enddo
      enddo
c     Option: Homogeneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      U_res_sum = 0.
c$$$      do j=1,Ny
c$$$         do i=1,Nx
c$$$            U_res_sum = U_res_sum + ((u_hat(i,j)+Ugal)**2.+
c$$$     +           (v_hat(i,j)+Vgal)**2.)**0.5
c$$$         end do
c$$$      end do
c$$$      U_res(:,:)=U_res_sum*iNxNy
c     Option: Heterogeneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(q_flag.eq.1)then
         call Filter_2dsl(t_res,theta(:,:,2),t,0)
      else
         call Filter_2dsl(t_res,theta(:,:,2),t,2)
      endif
c$$$  t_res=theta(:,:,2)
      if(q_flag.eq.1)then
         call Filter_2dsl(q_res,q(:,:,2),t,2)
c$$$  q_res=q(:,:,2)
      endif
c     Option: Homogeneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      t_res_sum = 0.
c$$$      do j=1,Ny
c$$$         do i=1,Nx
c$$$            t_res_sum = t_res_sum + theta(i,j,2)
c$$$         end do
c$$$      end do
c$$$      t_res_sum=t_res_sum*inxny
c$$$      t_res(:,:)=t_res_sum
c$$$      if(q_flag.eq.1)then
c$$$         q_res_sum = 0.
c$$$         do j=1,Ny
c$$$            do i=1,Nx
c$$$               q_res_sum = q_res_sum + q(i,j,2)
c$$$            end do
c$$$         end do
c$$$         q_res_sum=q_res_sum*inxny
c$$$         q_res(:,:)=q_res_sum
c$$$      endif
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      do kk=1,4
         
         t_flux_avg = 0.0
         q_flux_avg = 0.0
         ustar_avg = 0.0
         do j=1,Ny
            do i=1,Nx
               denom = dlog((dz/2.)/(zo(i,j)/z_i))+Psi(i,j,1)-
     +              Psi0(i,j,1)
               ustar(i,j) = (u_res(i,j)*vonk/denom)

c               denomH = dlog((dz/2.)/(zo(i,j)/z_i))+PsiH1(i,j,1)-
c     +              PsiH2(i,j,1)
               denomH = dlog((dz/2.)/(0.25d0/z_i))+PsiH1(i,j,1)-
     +              PsiH2(i,j,1)
               
               t_flux(i,j) = (t_s2used-t_res(i,j))*ustar(i,j)*
     +              vonk/denomH
               
               q_flux(i,j) = (q_s2used-q_res(i,j))*ustar(i,j)*
     +              vonk/denomH
               
               q_flux_avg = q_flux_avg+q_flux(i,j)
               t_flux_avg = t_flux_avg+t_flux(i,j)
               ustar_avg=ustar_avg+ustar(i,j)
            end do
         end do
         
         t_flux_avg=t_flux_avg*inxny
         ustar_avg=ustar_avg*inxny
         q_flux_avg=q_flux_avg*inxny

         cnt1=0
         cnt2=0
         cnt3=0
         
         do j=1,ny
            do i=1,nx

               if(qpass_flag.eq.0.and.q_flag.eq.1)then
                  OB_L(i,j)=-ustar(i,j)**3.*t_res(i,j)*
     +                 (1.0d0/(vonk*g_hat*(t_flux(i,j) +
     +                 0.61d0*theta_0*q_flux(i,j))))
               else
                  OB_L(i,j)=-ustar(i,j)**3.*t_res(i,j)*
     +                 (1.0d0/(vonk*g_hat*t_flux(i,j)))

c$$$                  OB_L(i,j)=-ustar_avg**3.*t_res_sum*
c$$$     +                 (1.0d0/(vonk*g_hat*t_flux_avg))
               endif
               
               do k=1,1
                  if(pass_flag.eq.0)then
                     z=+k*0.5d0*dz

		     if(z/OB_L(i,j).gt.5.d0)then 
                        OB_L(i,j)=z/5.d0
                        cnt1=cnt1+1
		     elseif(z/OB_L(i,j).lt.-1.d0)then 
                        OB_L(i,j)=-z/1.d0
                        cnt2=cnt2+1
                     endif

	             if(z/OB_L(i,j).lt.0.d0) cnt3=cnt3+1

                     if ((t_flux(i,j)).gt.0.) then
                        x=+(1.-(15.d0*z/OB_L(i,j)))**0.25
                        y=+(1.-(15.d0*(zo(i,j)/z_i)/OB_L(i,j)))**0.25
                        y2=+(1.-(15.d0*(0.25d0/z_i)/OB_L(i,j)))**0.25
                        Psi(i,j,k)=-2.d0*dlog(0.5*(1+x))-
     +                       dlog(0.5d0*(1+x**2.))+2.*atan(x)-Pi/2.
                        Psi0(i,j,k)=-2.d0*dlog(0.5*(1+y))-
     +                       dlog(0.5d0*(1+y**2.))+2.d0*atan(y)-Pi/2.d0
                        fi(i,j,k)=1./x
c     Equation 11.9 and 11.14 from Arya
                        psiH1(i,j,k)=-2.d0*dlog(0.5d0*(1+x**2.))
                        psiH2(i,j,k)=-2.d0*dlog(0.5d0*(1+y2**2.))
c                        psiH2(i,j,k)=-2.d0*dlog(0.5d0*(1+y**2.))
                        fiH(i,j,k)=fi(i,j,k)**2.0
                     else if ((t_flux(i,j)).lt.(-0.)) then
                        Psi(i,j,k)=+5.0d0*z/OB_L(i,j)
                        Psi0(i,j,k)=+5.0d0*(zo(i,j)/z_i)/OB_L(i,j)
                        fi(i,j,k)=1.d0+psi(i,j,k)
                        PsiH1(i,j,k)=5.0d0*z/OB_L(i,j)
                        PsiH2(i,j,k)=5.0d0*(0.25d0/z_i)/OB_L(i,j)
c                        PsiH2(i,j,k)=5.0d0*(zo(i,j)/z_i)/OB_L(i,j)
                        fiH(i,j,k)=1.d0+psiH1(i,j,k)
                     else
                        Psi(i,j,k)=0.d0
                        Psi0(i,j,k)=0.d0
                        fi(i,j,k)=1.d0
                        PsiH1(i,j,k) = Psi(i,j,k)
                        PsiH2(i,j,k)= Psi(i,j,k) 
                        fiH(i,j,k) = fi(i,j,k)
                     end if
                     
                  else
                     
                     Psi(i,j,k)=0.d0
                     Psi0(i,j,k)=0.d0
                     fi(i,j,k)=1.d0
                     PsiH1(i,j,k)=0.d0
                     PsiH2(i,j,k)=0.d0
                     fiH(i,j,k)=0.74d0 
                     
                  endif
                  
               end do
               
            enddo
         enddo
         
      end do
      
      OBL_avg=0.0
      OBL_min=1000
      OBL_max=-1000
c      do j=1,Ny
c         do i=1,Nx
c            OBL_avg=OBL_avg+0.5*dz/OB_L(i,j)
c	     if(0.5*dz/OB_L(i,j).gt.OBL_max) OBL_max=0.5*dz/OB_L(i,j)
c            if(0.5*dz/OB_L(i,j).lt.OBL_min) OBL_min=0.5*dz/OB_L(i,j)
c         enddo
c      enddo
      
c     write(1525,*) OBL_avg*z_i*inxny,ustar_avg,q_flux_avg
c     write (1999,*) -t_flux_avg/ustar_avg,-q_flux_avg/ustar_avg
c      write(*,*) 't,z/L,z/Lmax,z/Lmin',t,OBL_avg*inxny,
c     +     OBL_max,OBL_min,cnt1,cnt2,cnt3
c      write(*,*) 't,z/L,qflux',t,OBL_avg*inxny,q_flux_avg         

c      write(*,*) 'nrsub,t,t_flux,q_flux',nrsub,t,t_flux_avg,q_flux_avg
      
      return
      
      end

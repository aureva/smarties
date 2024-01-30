      subroutine surf_flux2(theta,q,u,v,u_hat,v_hat,t_flux,q_flux,Psi,
     +     Psi0,fi,fiH,zo,t_s,q_s,coolrate,ilow,U_res,t)
      implicit none
      include 'dimen.h'
      integer*4 i,j,k,kk,t,pos_cnt,neg_cnt,il,rlr,ip,ilow(nx),iii,
     +     ddd,swt,ssi,i_case,rule
      real*8, dimension(nx,ny,nz2):: u,v,theta,q
      real*8, dimension(nx,ny,2) :: Psi,PsiH,Psi0,PsiH0,fi,fiH
      real*8, dimension(nx,ny) :: zo,t_s,ustar,t_sn,coolrate,
     +     t_flux,OB_L,u_res,t_res,u_hat,v_hat,t_hat,clrtn,
     +     q_res,q_flux,q_s	 
      real*8 t_res_sum,u_res_sum,denom,t_flux_avg,z,x,y,denomH,
     +     OBL_avg,ustar_avg,t_bar,mhf,ssr,tmp1,tmp2,fac1,
     +     q_bar,q_flux_avg	 
      save swt     
 
      real*8,save:: tg(11),tg_old(11),t_bottom,h(11),cg,c_air,cgr,nug
      real*8,save:: tg3d(nx,ny,11),tg3d_old(nx,ny,11),hc,s_flux_dc
      real*8,save:: z0_ratio
     
      character*99 filename
      logical file_exists
 
      i_case=12
      if (t.eq.1) then
         mhf = ceiling(dx/(dt*Ugal)/2)
         swt = nint(10*dmod((mhf*dt*Ugal),dx)/dx)
         write(*,*) 'switching every -- ',swt
         do k=1,2
            Psi(:,:,k)=0.
            fi(:,:,k)=1.
            PsiH(:,:,k)=0.
            fiH(:,:,k)=1.0
         enddo
! surface roughness ratio between temperature and momentum
         z0_ratio = 0.1d0
         if (i_case .eq. 11) z0_ratio = 0.1d0
        
         hc = 0.05
! logarithmic spacing for land
         h(1) = 0.0
         h(2) = 0.01432*2.0
         h(3) = 0.03438*2.0
         h(4) = 0.06245*2.0
         h(5) = 0.10176*2.0
         h(6) = 0.15678*2.0
         h(7) = 0.23382*2.0
         h(8) = 0.34167*2.0
         h(9) = 0.49266*2.0
         h(10) = 0.70405*2.0
         h(11) = 1.0*2.0
         endif

c        call Filter_2dsl(t_hat,theta(:,:,2),t,1)
c        call Filter_2dsl(u_hat,u(:,:,2),t,0)
c        call Filter_2dsl(v_hat,v(:,:,2),t,0)
      
 
      if (i_case .eq. 0 .or. s_flag.eq.0)then
         u_hat=u(:,:,2)
         v_hat=v(:,:,2)
         u_res=dsqrt((u_hat)**2+(v_hat)**2)
         do k=1,2
            Psi(:,:,k)=0.
            fi(:,:,k)=1.
            PsiH(:,:,k)=0.
            fiH(:,:,k)=1.0
         enddo
 
         if (surf_flag .eq. 0) then
            t_flux=s_flux
         end if
      else
 
         if (surf_flag.eq.0) then
            t_flux=s_flux
 
! evaluate surface temperature profile using similarity theory, Hao Lu
            t_hat=theta(:,:,2)
            u_hat=u(:,:,2)
            v_hat=v(:,:,2)
 
            t_bar = 0.0
            do j=1,Ny
            do i=1,Nx
            u_res(i,j)=dsqrt((u_hat(i,j)+Ugal)**2+(v_hat(i,j)+Vgal)**2)
            t_res(i,j) = t_hat(i,j)
            t_bar = t_bar+t_res(i,j)
            end do
            end do
            t_bar = t_bar*inxny
 
!            t_s(:,:) = 3.d0*theta(:,:,2)-2.d0*theta(:,:,3)
            do j=1,Ny
            do i=1,Nx
              denom = dlog((dz/2.)/(zo(i,j)/z_i))+Psi(i,j,1)-Psi0(i,j,1)
              ustar(i,j) = (u_res(i,j)*vonk/denom)
              denomH = dlog((dz/2.)/(zo(i,j)*z0_ratio/z_i))
     +                 +PsiH(i,j,1)-PsiH0(i,j,1)
              t_s(i,j) = t_flux(i,j)/(ustar(i,j)*vonk/denomH)+t_res(i,j)
!              t_flux(i,j) = (t_s(i,j)-t_res(i,j))*ustar(i,j)*vonk/denomH
            end do
            end do
 
         elseif (surf_flag.eq.3) then
            ssr = floor((dble(t)*dt*Ugal)/dx)
            ssi = int(ssr)
            iii = mod(ssi,nx)
            ddd = nint(10*dmod((dble(t)*dt*Ugal),dx)/dx)
 
            rlr = t/2
            if(t.eq.1)then
               do i=1,nx   
                  ilow(i)=i
               enddo
            endif
           
            do i=1,nx
               il=i-iii
               if(il.lt.1) il=il+nx
               if(ddd.eq.swt)then
c                  if(i.eq.1) write(*,*) 'I moved...'
                  if(i.eq.nx)then
                     t_sn(i,:)=t_s(1,:)
                     clrtn(i,:)=coolrate(1,:)
                  else
                     t_sn(i,:)=t_s(i+1,:)
                     clrtn(i,:)=coolrate(i+1,:)
                  endif
               else
                  t_sn(i,:)=t_s(i,:)
                  clrtn(i,:)=coolrate(i,:)
               endif
              
c               if(t.eq.rlr*2) write(*,*) 'ts',i,il,
c     +              t_sn(i,ny/2)
 
               ilow(i)=il
 
            end do
 
!            if (t.gt.nint(9.*3600.*u_scale/z_i/dt))then
!               do j=1,ny
!                  do i=1,nx
!                     t_s(i,j)=t_sn(i,j)
!                  end do
!               end do
!            else
 
            tmp1 = (c_coolrate/T_scale/3600)*(z_i/u_scale)*dt
            do j=1,ny
            do i=1,nx
               t_s(i,j) = t_sn(i,j)-tmp1
            end do
            end do
 
!            endif
 
!            if(t.gt.nint(8.*3600.*u_scale/z_i/dt))then
!               do j=1,ny
!                  do i=1,nx
!                     t_s(i,j)=t_sn(i,j)-(c_coolrate/T_scale/3600)*
!     +                    (z_i/u_scale)*dt
!                  end do
!               end do
!            else
!               do j=1,ny
!                  do i=1,nx
!                     t_s(i,j)=t_sn(i,j)-(c_coolrate/T_scale/3600)*
!     +                    (z_i/u_scale)*dt
!c       Option: Heterogeneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!                     t_s(i,j)=t_sn(i,j)-(clrtn(i,j)/T_scale/3600)*
!!     +                    (z_i/u_scale)*dt
!                  end do
!               end do
!            endif
 
            coolrate=clrtn
 
         elseif (surf_flag .eq. 4) then
! fixed temperature underground 1m
            t_bottom = 265.0d0
            tg(11) = t_bottom
 
! heat capacity 1.3d6 J/m^3-K, diffusivity 5.0d-7 m^2/s
! air heat capacity is 1344 J/m^3-K
            cg = 1.3d6
            nug = 5.0d-7
            c_air = 1344.0d0
            cgr = cg/c_air
 
! evaluate surface temperature profile using similarity theory, Hao Lu
            t_hat=theta(:,:,2)
            u_hat=u(:,:,2)
            v_hat=v(:,:,2)
 
            t_bar = 0.0
            ustar_avg = 0.0d0
            do j=1,Ny
            do i=1,Nx
            u_res(i,j)=dsqrt((u_hat(i,j)+Ugal)**2+(v_hat(i,j)+Vgal)**2)
            denom = dlog((dz/2.)/(zo(i,j)/z_i))+Psi(i,j,1)-Psi0(i,j,1)
            ustar(i,j) = (u_res(i,j)*vonk/denom)
            ustar_avg = ustar_avg+ustar(i,j)
            t_res(i,j) = t_hat(i,j)
            t_bar = t_bar+t_res(i,j)
!              denomH = dlog((dz/2.)/(zo(i,j)/z_i))
!     +                 +PsiH(i,j,1)-PsiH0(i,j,1)
!              t_s(i,j) = t_flux(i,j)/(ustar(i,j)*vonk/denomH)+t_res(i,j)
!!              t_flux(i,j) = (t_s(i,j)-t_res(i,j))*ustar(i,j)*vonk/denomH
            end do
            end do
            ustar_avg = ustar_avg*inxny
            t_bar = t_bar*inxny
 
            if (t .eq. 1) then
               do i = 2, 10
                  tg(i) = t_bottom
               end do
               tg(1) = t_bar
            end if
 
            do i = 1, 11
               tg_old(i) = tg(i)
            end do
           
            fac1 = dt/(u_scale/z_i)*nug
            do i = 2, 10
               tg(i) = tg_old(i)+fac1*(tg_old(i+1)*(h(i)-h(i-1))+
     +          tg_old(i-1)*(h(i+1)-h(i))-tg_old(i)*(h(i+1)-h(i-1)))/
     +          (0.5d0*(h(i+1)-h(i-1))*(h(i+1)-h(i))*(h(i)-h(i-1)))
            end do
            tg(11) = t_bottom
           
            fac1 = 0.0d0
            do j=1,Ny
            do i=1,Nx
            denomH = dlog((dz/2.)/(zo(i,j)*z0_ratio/z_i))
     +                 +PsiH(i,j,1)-PsiH0(i,j,1)
            fac1 = fac1+denomH
            end do
            end do
            fac1 = fac1*inxny
            fac1 = vonk/fac1*ustar_avg
            tg(1) = (s_flux+fac1*t_bar+tg(2)*cgr*nug/(h(2)-h(1)))/
     +              (cgr*nug/(h(2)-h(1))+fac1)
            t_s = tg(1)
 
!            if (tg(1) .le. tg(2) .or. tg(1) .le. t_bar) then
!               write(*,*) 'ts=',tg(1),'tg=',tg(2),'tair=',t_bar
!               write(*,*) 'Qg=',(tg(1)-tg(2))*cgr*nug/h,
!     +                    'HF=',fac1*(tg(1)-t_bar)
!            end if
            if (mod(t,c_count) .eq. 0) then
               write(*,*) 'Qg=',(tg(1)-tg(2))*cgr*nug/(h(2)-h(1)),
     +                    'HF=',fac1*(tg(1)-t_bar)
               do i = 1, 11
                  write(*,*) h(i),tg(i)
               end do
            end if
           
            do j=1,Ny
            do i=1,Nx
              denomH = dlog((dz/2.)/(zo(i,j)*z0_ratio/z_i))
     +                 +PsiH(i,j,1)-PsiH0(i,j,1)
              t_s(i,j) = tg(1)
              t_flux(i,j) = (t_s(i,j)-t_res(i,j))*ustar(i,j)*vonk/denomH
            end do
            end do
 
         elseif (surf_flag .eq. 1) then
	 
c         s_flux_dc=0.1*(dexp(-(t*dt*z_i/3600-12.)**2/2/10)+
c     +                 dexp(-(t*dt*z_i/3600-36.)**2/2/10))-0.025		 
cc         s_flux_dc=0.1
         s_flux_dc=-0.1*dcos((t*dt*z_i/3600-12.0)*Pi/12.0)		 
         if (s_flux_dc .le. -0.025) 	s_flux_dc=-0.025
		 
c         s_flux_dc=0.04
! fixed temperature underground 1m
            t_bottom = 293.0d0
 
! heat capacity 1.3d6 J/m^3-K, diffusivity 5.0d-7 m^2/s
! air heat capacity is 1344 J/m^3-K
c            cg = 4.18d6
c            nug = 1.43d-7	

            cg = 1.3d6
            nug = 5.0d-7		
            c_air = 1344.0d0
            cgr = cg/c_air
 
            if (t .eq. 1) then
cc               filename = 'tg3d_ini.bin'
cc               inquire(FILE=filename, EXIST=file_exists)
cc               if (file_exists) then
cc                  call load_matrix_1P(tg3d,Nx,Ny,11,1.0d0,filename)
cc               else

cc               print*, 'yes,its done'
               open(unit=120,file='input/tg3d_ini.out',
     +                    status='unknown')
                do k=1,11
                  do j=1,Ny
                    do i=1,Nx                
                      read (120,*) tg3d(i,j,k)
                    end do
                  end do
                end do
                close(120)
				
cc                  t_s=t_bottom
cc                  tg3d(:,:,1) = t_s
cc                  fac1 = dsqrt(Pi/(24*3600.0d0*nug))
cc                  do k = 2, 10
cc                  do j = 1, ny
cc                  do i = 1, nx
cc                     tmp1 = (t_s(i,j)-t_bottom)
cc                     tg3d(i,j,k) = t_bottom+tmp1*dexp(-fac1*h(k))
cc                  end do
cc                  end do
cc                  end do
                  tg3d(:,:,11) = t_bottom
                  tg3d(:,:,:) = t_bottom	  
cc               end if
            end if
 
            tg3d_old = tg3d
 
            do k = 2, 10
               fac1 = dt/(u_scale/z_i)*nug
!               fac1 = dt/(u_scale/z_i)*nug/(hc**2)
               do j = 1, ny
               do i = 1, nx
                  tg3d(i,j,k) = tg3d_old(i,j,k)+fac1*(tg3d_old(i,j,k+1)*
     +             (h(k)-h(k-1))+tg3d_old(i,j,k-1)*(h(k+1)-h(k))-
     +              tg3d_old(i,j,k)*(h(k+1)-h(k-1)))/
     +              (0.5d0*(h(k+1)-h(k-1))*(h(k+1)-h(k))*(h(k)-h(k-1)))
!                  tg3d(i,j,k) = tg3d_old(i,j,k)+fac1*(tg3d_old(i,j,k+1)+
!     +             tg3d_old(i,j,k-1)-2*tg3d_old(i,j,k))
               end do
               end do
 
               fac1 = dt/(u_scale/z_i)*nug/((dx*z_i)**2)
               i = 1
               do j = 1, ny
                  tg3d(i,j,k) = tg3d(i,j,k)+fac1*(tg3d_old(i+1,j,k)-
     +                 tg3d_old(i,j,k)-tg3d_old(i,j,k)+tg3d_old(nx,j,k))
               end do
               do j = 1, ny
               do i = 2, nx-1
                  tg3d(i,j,k) = tg3d(i,j,k)+fac1*(tg3d_old(i+1,j,k)-
     +                tg3d_old(i,j,k)-tg3d_old(i,j,k)+tg3d_old(i-1,j,k))
               end do
               end do
               i = nx
               do j = 1, ny
                  tg3d(i,j,k) = tg3d(i,j,k)+fac1*(tg3d_old(1,j,k)-
     +                tg3d_old(i,j,k)-tg3d_old(i,j,k)+tg3d_old(i-1,j,k))
               end do
 
               fac1 = dt/(u_scale/z_i)*nug/((dy*z_i)**2)
               j = 1
               do i = 1, nx
                  tg3d(i,j,k) = tg3d(i,j,k)+fac1*(tg3d_old(i,j+1,k)-
     +                 tg3d_old(i,j,k)-tg3d_old(i,j,k)+tg3d_old(i,ny,k))
               end do
               do j = 2, ny-1
               do i = 1, nx
                  tg3d(i,j,k) = tg3d(i,j,k)+fac1*(tg3d_old(i,j+1,k)-
     +                tg3d_old(i,j,k)-tg3d_old(i,j,k)+tg3d_old(i,j-1,k))
               end do
               end do
               j = ny
               do i = 1, nx
                  tg3d(i,j,k) = tg3d(i,j,k)+fac1*(tg3d_old(i,1,k)-
     +                tg3d_old(i,j,k)-tg3d_old(i,j,k)+tg3d_old(i,j-1,k))
               end do
            end do
            tg3d(:,:,11) = t_bottom
 
! evaluate surface temperature profile using similarity theory, Hao Lu
            t_hat=theta(:,:,2)
            u_hat=u(:,:,2)
            v_hat=v(:,:,2)
 
            tmp1 = cgr*nug/((h(2)-h(1))*(h(3)-h(1))*(h(3)-h(2)))
            tmp2 = tmp1*((h(3)-h(1))**2-(h(2)-h(1))**2)
            do j=1,Ny
            do i=1,Nx
             u_res(i,j)=dsqrt((u_hat(i,j)+Ugal)**2+(v_hat(i,j)+Vgal)**2)
             denom = dlog((dz/2.)/(zo(i,j)/z_i))+Psi(i,j,1)-Psi0(i,j,1)
             ustar(i,j) = (u_res(i,j)*vonk/denom)
             t_res(i,j) = t_hat(i,j)
             denomH = dlog((dz/2.)/(zo(i,j)*z0_ratio/z_i))
     +                 +PsiH(i,j,1)-PsiH0(i,j,1)
             fac1 = vonk/denomH*ustar(i,j)
             
             if(q_flag.eq.1)then
                q_res=q(:,:,2) 	
c                q_s(i,j)=3.801664*exp(17.67*(t_s(i,j)-273.16)
c     +                    /(t_s(i,j)-29.66))				
                q_flux(i,j) = (q_s(i,j)-q_res(i,j))*ustar(i,j)*
     +                   vonk/denomH	
             endif	 
			 
! second order
!             tg3d(i,j,1) = (s_flux_dc+fac1*t_res(i,j)+tmp1*(-tg3d(i,j,3)*
!     +           (h(2)-h(1))**2+tg3d(i,j,2)*(h(3)-h(1))**2))/(tmp2+fac1)
! first order            
             tg3d(i,j,1)=(s_flux_dc+fac1*t_res(i,j)+tg3d(i,j,2)*cgr*nug/
     +                    (h(2)-h(1)))/(cgr*nug/(h(2)-h(1))+fac1)
             t_s(i,j) = tg3d(i,j,1)
             t_flux(i,j) = (t_s(i,j)-t_res(i,j))*ustar(i,j)*vonk/denomH
            end do
            end do
 
            if (mod(t,c_count) .eq. 0) then
               t_bar = 0.0d0
               do j = 1, Ny
               do i = 1, Nx
                  t_bar = t_bar+t_res(i,j)
               end do
               end do
               t_bar = t_bar*inxny
              
               do k = 1, 11
                  tg(k) = 0.0d0
                  do j = 1, Ny
                  do i = 1, Nx
                     tg(k) = tg(k)+tg3d(i,j,k)
                  end do
                  end do
                  tg(k) = tg(k)*inxny
               end do
c               do i = 1, 11
c                  write(*,*) h(i),tg(i)
c               end do
 
               filename = 'output/tg3d'
cc               call strdotint(filename,t,7)
cc               call save_matrix_1P(tg3d,Nx,Ny,11,1.0d0,filename)
            end if
 
         end if
 
!         if (i_case .eq. 1 .and. t .eq. 1) then
!            filename = 't_s_ini.bin'
!            call load_matrix_1P(t_s,nx,ny,1,1.0d0,filename)
!         end if
                 
c       Option: Homogeneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!         U_res_sum = 0.
!         t_res_sum = 0.
!         do j=1,Ny
!            do i=1,Nx
!               U_res_sum = U_res_sum + ((u(i,j,2)+Ugal)**2.+
!     +              v(i,j,2)**2.)**0.5
!               t_res_sum = t_res_sum + theta(i,j,2)
!            end do
!         end do
!         do j=1,Ny
!            do i=1,Nx
!               U_res(i,j)=U_res_sum*iNxNy
!               t_res(i,j)=t_res_sum*iNxNy
!            end do
!         end do
!        t_bar = t_res_sum*iNxNy
c       Option: Heterogeneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        call Filter_2dsl(t_hat,theta(:,:,2),t,1)
!        call Filter_2dsl(u_hat,u(:,:,2),t,0)
!        call Filter_2dsl(v_hat,v(:,:,2),t,0)
 
        t_hat=theta(:,:,2)
        u_hat=u(:,:,2)
        v_hat=v(:,:,2)
 
        t_bar = 0.0
        q_bar =0.0	
		
        do j=1,Ny
        do i=1,Nx
         U_res(i,j) = dsqrt((u_hat(i,j)+Ugal)**2+(v_hat(i,j)+Vgal)**2)
         t_res(i,j) = t_hat(i,j)
         t_bar = t_bar+theta(i,j,2)
		 
        q_bar = q_bar+q(i,j,2)		 
        end do
        end do
        t_bar = t_bar*inxny
        q_bar = q_bar*inxny		
c       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        do kk=1,2
        
c        t_flux_avg = 0.0
         do j=1,Ny
            do i=1,Nx
               denom = dlog((dz/2.)/(zo(i,j)/z_i))+Psi(i,j,1)
     +                 -Psi0(i,j,1)
               ustar(i,j) = (u_res(i,j)*vonk/denom)
               if(surf_flag.eq.1)then
                  denomH = dlog((dz/2.)/(zo(i,j)*z0_ratio/z_i))+
     +                 PsiH(i,j,1)-PsiH0(i,j,1)
                  t_flux(i,j) = (t_s(i,j)-t_res(i,j))*ustar(i,j)*
     +                 vonk/denomH
               endif
c               t_flux_avg = t_flux_avg+t_flux(i,j)
c                     ustar_avg=ustar_avg+ustar(i,j)
            end do
         end do
  
! the origional definition is more steady, by Hao Lu
c         ustar_avg=0.0d0
c         t_flux_avg=0.0d0
c         OBL_avg=0.0d0
c         do j=1, Ny
c         do i=1, Nx
c            ustar_avg = ustar_avg+ustar(i,j)
c            t_flux_avg = t_flux_avg+t_flux(i,j)
c            OBL_avg = OBL_avg+t_s(i,j)
c         end do
c         end do
c         ustar_avg = ustar_avg/(Nx*Ny)
c         t_flux_avg = t_flux_avg/(Nx*Ny)
c         OBL_avg = OBL_avg/(Nx*Ny)
c         OBL_avg = -(ustar_avg**3)*OBL_avg/(vonk*g_hat*t_flux_avg)
c         OB_L = OBL_avg
 
c         t_flux_avg=t_flux_avg*inxny
c                ustar_avg=ustar_avg*inxny
 
         pos_cnt=0
         neg_cnt=0
 
         do j=1,ny
            do i=1,nx
c               OB_L(i,j)=-ustar(i,j)**3.*t_res(i,j)/
c     +               (vonk*g_hat*t_flux(i,j))
               OB_L(i,j)=-ustar(i,j)**3.*theta_0/
     +               (vonk*g_hat*t_flux(i,j))
     
               do k=1,2
                  if(pass_flag.eq.0)then
                     z=+k*0.5*dz
                     if(z/OB_L(i,j).gt.5.)then
                        OB_L(i,j)=z/5.
                     elseif(z/OB_L(i,j).lt.-5.)then
                        OB_L(i,j)=-z/5.
                     endif
                     if ((t_flux(i,j)).gt.0.) then
                        if(k.eq.1) pos_cnt=pos_cnt+1
                        x=+(1.-(15.*z/OB_L(i,j)))**0.25
                        y=+(1.-(15.*(zo(i,j)/z_i)/OB_L(i,j)))**0.25
! see eqn 9.7.5i in Stull's "an intro. to blm" book
                        Psi(i,j,k)=-2.*dlog(0.5*(1+x))-
     +                       dlog(0.5*(1+x**2))+2.*atan(x)-Pi/2.
                        Psi0(i,j,k)=-2.*dlog(0.5*(1+y))-
     +                       dlog(0.5*(1+y**2))+2.*atan(y)-Pi/2.
! see eqn 9.7.5c in Stull's "an intro. to blm" book
                        fi(i,j,k)=1./x
! Eqn. 11.14 from Arya "Introduction2Micrometeorology" orginally from
! Paulson 1970 "the mathematical representation of wind speed ..."
                        psiH(i,j,k)=-2*dlog(0.5*(1+x**2))
                        psiH0(i,j,k)=-2*dlog(0.5*(1+y**2))
! Eqn. 11.9 from Arya "Introduction2Micrometeorology"
                        fiH(i,j,k)=fi(i,j,k)**2.0
! see eqn 9.6.1f, 9.7.5F (correct the power) Stull's "an intro. to blm"
c                        fiH(i,j,k)=0.74/((1.-(9.*z/OB_L(i,j)))**0.5)
                     else if ((t_flux(i,j)).lt.(-0.)) then
                        if(k.eq.1) neg_cnt=neg_cnt+1
                        Psi(i,j,k)=+4.7*z/OB_L(i,j)
                        Psi0(i,j,k)=+4.7*(zo(i,j)/z_i)/OB_L(i,j)
                        fi(i,j,k)=1.+psi(i,j,k)
! according to GABLS description
                        PsiH(i,j,k)=7.8*z/OB_L(i,j)
                        PsiH0(i,j,k)=7.8*(zo(i,j)/z_i)/OB_L(i,j)
                        fiH(i,j,k)=1.0+psiH(i,j,k)
c                        PsiH(i,j,k)=4.7*z/OB_L(i,j)
c                        PsiH0(i,j,k)=4.7*(zo(i,j)/z_i)/OB_L(i,j)
c                        fiH(i,j,k)=0.74+psiH(i,j,k)
                     else
                        Psi(i,j,k)=0.
                        Psi0(i,j,k)=0.
                        fi(i,j,k)=1.
                        PsiH(i,j,k) = Psi(i,j,k)
                        PsiH0(i,j,k)= Psi0(i,j,k)
                        fiH(i,j,k) = 1.
                     end if
                    
                  else
 
                     Psi(i,j,k)=0.
                     Psi0(i,j,k)=0.
                     fi(i,j,k)=1.
                     PsiH(i,j,k)=0.
                     PsiH0(i,j,k)=0.
                     fiH(i,j,k)=1.
 
                  endif
              end do
              
            enddo
         enddo
        
      end do
  
         rule=t/ruler    
         IF (rule*ruler.eq.t) then 
               open(unit=120,file='input/tg3d_ini.out',
     +              status='unknown')
	        do k=1,11
               do j=1,ny
                  do i=1,nx
                     write(120,*) tg3D(i,j,k)
                     write(897,*) tg3D(i,j,k) 					 
                  enddo
               enddo
            enddo   
               close(120)		 			   
         end if

! the origional definition is more steady, by Hao Lu
         ustar_avg=0.0d0
         t_flux_avg=0.0d0
         q_flux_avg=0.0d0		 
c         OBL_avg=0.0d0
         do j=1, Ny
         do i=1, Nx
            ustar_avg = ustar_avg+ustar(i,j)
            t_flux_avg = t_flux_avg+t_flux(i,j)
            q_flux_avg = q_flux_avg+q_flux(i,j)			
!            OBL_avg = OBL_avg+t_s(i,j)
         end do
         end do
         ustar_avg = ustar_avg/(Nx*Ny)
         t_flux_avg = t_flux_avg/(Nx*Ny)
         q_flux_avg = q_flux_avg/(Nx*Ny)		 
!         OBL_avg = OBL_avg/(Nx*Ny)
 
! set the reference temperature as the averaged surface temperature
c         OBL_avg = OBL_avg*z_i
!         OBL_avg = -(ustar_avg**3)*OBL_avg/(vonk*g_hat*t_flux_avg)*z_i
 
         OBL_avg=0.0
         do j=1,Ny
            do i=1,Nx
               OBL_avg=OBL_avg+OB_L(i,j)
            enddo
         enddo
         OBL_avg=OBL_avg*z_i/(Nx*Ny)

         if (mod(t,c_count)==0) then 
c         write(1525,*) OBL_avg*z_i*inxny,ustar_avg
           write (1999,5111) -t_flux_avg/ustar_avg,t_flux_avg,t_bar
     +                       ,OBL_avg,ustar_avg,tg(1),tg(10),s_flux_dc
     +                       ,q_flux_avg*2.5,q_s(1,1),q_bar 
cc            write(*,5555) 'tflux,u*',t_flux_avg,ustar_avg,t_bar,OBL_avg 
         endif			  
         
      endif
 
 5111 format (15(1x,f15.5))  	  
 5555 format(1(1x,a14),5(1x,f15.8))  
      return
 
      end
	  
	  
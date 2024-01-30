      Subroutine SGS_STAG (t_flux,qx,qy,qz,sgs_q1,sgs_q2,sgs_q3,
     +     u,v,w,dudx,dudy,dudz,dvdx,dvdy,
     +     dvdz,dwdx,dwdy,dwdz,dtdx,dtdy,dtdz,dqdx,dqdy,dqdz,t,
     +     Cs2,Pr2,Sc2,cst,txx,txy,txz,tyy,tyz,tzz,theta,q,
     +     a1_old,b1_old,c1_old,d1_old,e1_old,
     +     a2_old,b2_old,c2_old,d2_old,e2_old,
     +     a4_old,b4_old,c4_old,d4_old,e4_old,
     +     a8_old,b8_old,c8_old,d8_old,e8_old,
     +     a5_old,b5_old,c5_old,d5_old,e5_old,
     +     a9_old,b9_old,c9_old,d9_old,e9_old,
     +     beta1,beta2,beta3,ESGS3D,ET3D,EQ3D,ncs2g,ncs2l,
     +     nx,ny,nz,nx2,ny2,
     +  nz2,nzb,nsteps,inxny,averaging,Co,cs_count,delta,dx,dy,dz,
     + fgr,g_hat,me,model,mom_nodes,nall,nnn,nprocs,pass_flag,
     +q_flag,ri_flag,s_flag,sc,scl_nodes,t_scale,theta_0,
     +vonk,inx2ny2,nu,plan_f,plan_b,plan_ff,plan_bb)
      
      implicit none
!      include 'dimen.h'

      integer*4 :: nx,ny,nz,nx2,ny2,nz2,nzb,nsteps,averaging
      integer*4 :: cs_count,me,model,mom_nodes,nall,nprocs,pass_flag
      integer*4 :: q_flag,ri_flag,s_flag,scl_nodes
      real*8 :: inxny,Co,delta,dx,dy,dz,fgr,g_hat,nnn,sc,t_scale
      real*8 :: theta_0, vonk,inx2ny2,nu
      integer*8 :: plan_f,plan_b,plan_ff,plan_bb

      integer*4 i,j,k,t,pt,cst,ccount                    
      real*8,dimension(nx,ny,nz2):: dudx,dudy,dudz,dvdx,dvdy,dvdz,
     +     dwdx,dwdy,dwdz,ux,uy,uz,vx,vy,vz,wx,wy,wz,u,v,w,u_,v_,w_,S,
     +     txx,txy,txz,tyy,tyz,tzz,S11,S12,S13,S22,S23,S33,S_hat,S_hatd,
     +     dtdx,dtdy,dtdz,theta,tx,ty,tz,qx,qy,qz,t_,Rip,
     +     dqdx,dqdy,dqdz,dqx,dqy,dqz,q,q_,sgs_q1,sgs_q2,sgs_q3

      real*8,dimension(nx,ny,nz2)::
     +     a1_old,b1_old,c1_old,d1_old,e1_old,
     +     a2_old,b2_old,c2_old,d2_old,e2_old,Cs2,beta1,
     +     a4_old,b4_old,c4_old,d4_old,e4_old,
     +     a8_old,b8_old,c8_old,d8_old,e8_old,Pr2,beta2,
     +     a5_old,b5_old,c5_old,d5_old,e5_old,
     +     a9_old,b9_old,c9_old,d9_old,e9_old,Sc2,beta3,
     +     ESGS3D,ET3D,EQ3D
      
      real*8,dimension(nx2,ny2,nz2):: S11_m,S12_m,S13_m,S22_m,S23_m,
     +     S33_m,S_m,txx_m,txy_m,txz_m,tyy_m,tyz_m,tzz_m,qx_m,qy_m,qz_m,
     +     tx_m,ty_m,tz_m,S11_mu,S12_mu,S13_mu,S22_mu,S23_mu,
     +     S33_mu,S_mu,Cs2_m,Pr2_m,
     +     sgs_q1_m,sgs_q2_m,sgs_q3_m,Sc2_m,dqx_m,dqy_m,dqz_m

      real*8,dimension(nx2,ny2):: Cs2_L,Pr2_L,Sc2_L

      real*8 l(Nz2),
     +     ct,txzp(nx,ny),tyzp(nx,ny),zz(nz2),
     +     cs_ave,b1_ave,pr_ave,b2_ave,sc_ave

      real*8 t_flux(nx,ny),ncs2g,ncs2l  
!      save Cs2_m,Pr2_m,Sc2_m

      if (model .eq. 100 .and. S_FLAG .ne. 1) then
         call AMD_Const(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,
     +      txx,txy,txz,tyy,tyz,tzz,t,me,nall,
     +      txx_m,txy_m,txz_m,tyy_m,tyz_m,tzz_m,
     +      S11,S22,S33,S12,S13,S23,
     +      S11_m,S12_m,S13_m,S22_m,S23_m,S33_m,
     +      ux,uy,uz,vx,vy,vz,wx,wy,wz,Cs2,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + dx,dy,dz,nu,plan_f,plan_b,plan_ff,plan_bb)
         return
      end if	  
      if (model .eq. 100 .and. S_FLAG .eq. 1) then
!      if (model .eq. 10) then	  
         call AMD_Scalar(dudx,dudy,dudz,dvdx,dvdy,dvdz,
     +   dwdx,dwdy,dwdz,t_flux,
     +   u,v,w,txx,txy,txz,tyy,tyz,tzz,
     +   dtdx,dtdy,dtdz,qx,qy,qz,t,me,nall,
     +   txx_m,txy_m,txz_m,tyy_m,tyz_m,tzz_m,
     +   S11,S12,S13,S22,S23,S33,
     +   S11_m,S12_m,S13_m,S22_m,S23_m,S33_m,
     +   ux,uy,uz,vx,vy,vz,wx,wy,wz,
     +   tx,ty,tz,tx_m,ty_m,tz_m,qx_m,qy_m,qz_m,Cs2,Pr2)
        return
      end if


	  
      if (model .eq. 10 .and. S_FLAG .ne. 1) then
         call MGM_Const(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,
     +      txx,txy,txz,tyy,tyz,tzz,t,me,nall,
     +      txx_m,txy_m,txz_m,tyy_m,tyz_m,tzz_m,
     +      S11,S22,S33,S12,S13,S23,
     +      S11_m,S12_m,S13_m,S22_m,S23_m,S33_m,
     +      ux,uy,uz,vx,vy,vz,wx,wy,wz)
         return
      end if

      if (model .eq. 10 .and. S_FLAG .eq. 1) then
!      if (model .eq. 10) then	  
         call MGM_Scalar(dudx,dudy,dudz,dvdx,dvdy,dvdz,
     +   dwdx,dwdy,dwdz,t_flux,
     +   u,v,w,txx,txy,txz,tyy,tyz,tzz,
     +   dtdx,dtdy,dtdz,qx,qy,qz,t,me,nall,
     +   txx_m,txy_m,txz_m,tyy_m,tyz_m,tzz_m,
     +   S11,S12,S13,S22,S23,S33,
     +   S11_m,S12_m,S13_m,S22_m,S23_m,S33_m,
     +   ux,uy,uz,vx,vy,vz,wx,wy,wz,
     +   tx,ty,tz,tx_m,ty_m,tz_m,qx_m,qy_m,qz_m)
        return
      end if	  
	  
c...  Cs is Smagorinsky's constant. l is a filter size (non-dim.)  

      if (me==0) then
         do j=1,ny
            do i=1,nx
               txzp(i,j)=txz(i,j,2)
               tyzp(i,j)=tyz(i,j,2)
            end do
         end do
      end if
c.....For traditional Smagorinsky

      if(model.eq.1)then

         do k=1,Nzb+2
            do j=1,ny
               do i=1,nx
                  Cs2(i,j,k)=Co**2.
                  Pr2(i,j,k)=Co**2./Sc
                  Sc2(i,j,k)=Co**2./Sc
               end do
            end do
            do j=1,ny2
               do i=1,nx2
                  Cs2_m(i,j,k)=Co**2.
                  Pr2_m(i,j,k)=Co**2./Sc
                  Sc2_m(i,j,k)=Co**2./Sc
               end do
            end do

            if(me==0) then
               if(k==2) then
                  zz(k)=(k-1.-0.5)*dz
               else
                  zz(k)=(k-1.-1.)*DZ
               end if
            else
               zz(k)=(k-1.+me*nzb-1.)*DZ
            end if
            
            l(k)=(Co**(nnn)*(vonk*zz(k))**(-nnn)+
     +           (delta)**(-nnn))**(-1./nnn)
c            l(k)=delta
         end do

      else

         do k=1,Nzb+2
            l(k)=fgr*(DX*DY*DZ)**(1./3.)
         end do

      end if
      
      if(mom_nodes.eq.0)then
         do k=2,nzb+1 
            u_(:,:,k)=0.5*(u(:,:,k)+u(:,:,k-1))
            v_(:,:,k)=0.5*(v(:,:,k)+v(:,:,k-1))
            w_(:,:,k)=w(:,:,k)

            ux(:,:,k)=0.5*(dudx(:,:,k)+dudx(:,:,k-1)) 
            uy(:,:,k)=0.5*(dudy(:,:,k)+dudy(:,:,k-1)) 
            uz(:,:,k)=     dudz(:,:,k)
            vx(:,:,k)=0.5*(dvdx(:,:,k)+dvdx(:,:,k-1))
            vy(:,:,k)=0.5*(dvdy(:,:,k)+dvdy(:,:,k-1))
            vz(:,:,k)=     dvdz(:,:,k)
            wx(:,:,k)=     dwdx(:,:,k)
            wy(:,:,k)=     dwdy(:,:,k)
            wz(:,:,k)=0.5*(dwdz(:,:,k)+dwdz(:,:,k-1))
         end do
         
      else

         do k=2,nzb+1 
            u_(:,:,k)=u(:,:,k)
            v_(:,:,k)=v(:,:,k)
            w_(:,:,k)=0.5*(w(:,:,k)+w(:,:,k+1))
            
            ux(:,:,k)=     dudx(:,:,k) 
            uy(:,:,k)=     dudy(:,:,k) 
            uz(:,:,k)=0.5*(dudz(:,:,k)+dudz(:,:,k+1))
            vx(:,:,k)=     dvdx(:,:,k)
            vy(:,:,k)=     dvdy(:,:,k)
            vz(:,:,k)=0.5*(dvdz(:,:,k)+dvdz(:,:,k+1))
            wx(:,:,k)=0.5*(dwdx(:,:,k)+dwdx(:,:,k+1))
            wy(:,:,k)=0.5*(dwdy(:,:,k)+dwdy(:,:,k+1))
            wz(:,:,k)=     dwdz(:,:,k)
         end do
         
      endif

      if (me.eq.0) then
         u_(:,:,2)=u(:,:,2)
         v_(:,:,2)=v(:,:,2)
         w_(:,:,2)=w(:,:,2+1)/2.

         ux(:,:,2)=     dudx(:,:,2)
         uy(:,:,2)=     dudy(:,:,2) 
         uz(:,:,2)=     dudz(:,:,2)
         vx(:,:,2)=     dvdx(:,:,2)
         vy(:,:,2)=     dvdy(:,:,2)
         vz(:,:,2)=     dvdz(:,:,2)
         wx(:,:,2)=0.5*(dwdx(:,:,3)+dwdx(:,:,2))
         wy(:,:,2)=0.5*(dwdy(:,:,3)+dwdy(:,:,2))
         wz(:,:,2)=     dwdz(:,:,2)   
      endif

      if(t.eq.1)then
         cst=cs_count
      end if

      do k=2,Nzb+1
         S11(:,:,k)=0.5*(ux(:,:,k)+ux(:,:,k))
         S33(:,:,k)=0.5*(wz(:,:,k)+wz(:,:,k))
         S22(:,:,k)=0.5*(vy(:,:,k)+vy(:,:,k))
         S12(:,:,k)=0.5*(uy(:,:,k)+vx(:,:,k))
         S13(:,:,k)=0.5*(uz(:,:,k)+wx(:,:,k))
         S23(:,:,k)=0.5*(vz(:,:,k)+wy(:,:,k))
      end do

      if(mom_nodes.eq.0)then

         call dealias1(S11,S11_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(S12,S12_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(S13,S13_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(S22,S22_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(S23,S23_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(S33,S33_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)

         do k=2,nzb+1
            S_m(:,:,k)=sqrt(2*(S11_m(:,:,k)**2+S22_m(:,:,k)**2+
     +           S33_m(:,:,k)**2+2.*S12_m(:,:,k)**2+2.*S13_m(:,:,k)**2+
     +           2.*S23_m(:,:,k)**2))
         enddo

      else

         call dealias1(S11,S11_mu,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(S12,S12_mu,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(S13,S13_mu,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(S22,S22_mu,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(S23,S23_mu,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(S33,S33_mu,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         
         do k=2,nzb+1
            S_mu(:,:,k)=sqrt(2*(S11_mu(:,:,k)**2+S22_mu(:,:,k)**2+
     +           S33_mu(:,:,k)**2+2.*S12_mu(:,:,k)**2+
     +           2.*S13_mu(:,:,k)**2+2.*S23_mu(:,:,k)**2))
         enddo

      endif

c.....##### SMAGORINSKY: skip this (not call Optim_cs -only for dynamic!)

      if(model.gt.1)then
      
         if (cst.eq.cs_count) then
            
            if (model.eq.2) then
               
               if(averaging.eq.1)then

                  call update1(a1_old,me,nall)
                  call update1(b1_old,me,nall)

                  call optim_lag_dyn(Cs2,S11,S33,S22,S12,S13,S23,S,
     +                 S_hat,u_,v_,w_,l,t,me,a1_old,b1_old)

               endif

            elseif(model.eq.3)then

               if (averaging.eq.1) then

                  call update9(a1_old,b1_old,c1_old,d1_old,e1_old,
     +                 a2_old,b2_old,c2_old,d2_old,me,nall)
                  call update1(e2_old,me,nall)

               call optim_lag(Cs2,S11,S33,S22,S12,S13,S23,S,S_hat,
     +                 S_hatd,u_,v_,w_,l,t,me,beta1,a1_old,b1_old,
     +                 c1_old,d1_old,e1_old,a2_old,b2_old,c2_old,d2_old,
     +                 e2_old,ncs2g,ncs2l)
                 
               else if (averaging.eq.0) then
                  call optim_pln (Cs2,S11,S33,S22,S12,S13,S23,S,S_hat,
     +                 S_hatd,u_,v_,w_,L,t,me,beta1)
                  
               else if (averaging.eq.2) then
                  
c     call optim_loc (Cs2,S11,S33,S22,S12,S13,S23,S,S_hat,
c     +              S_hatd,u_,v_,w_,L,t,me,beta1)
                  
               else if (averaging.eq.3) then
                  
c     call optim_wong (Cs2,S11,S33,S22,S12,S13,S23,S,S_hat,
c     +              S_hatd,u_,v_,w_,L,t,me,beta1)
                  
               end if
               
            end if
            
            if(S_FLAG.eq.1)then
               
               if(scl_nodes.eq.0)then
                  do k=2,nzb+1
                     t_(:,:,k)=0.5*(theta(:,:,k)+theta(:,:,k-1))
                     tx(:,:,k)= 0.5*(dtdx(:,:,k)+dtdx(:,:,k-1))
                     ty(:,:,k)= 0.5*(dtdy(:,:,k)+dtdy(:,:,k-1))
                     tz(:,:,k)=      dtdz(:,:,k)
                  enddo
               else
                  do k=2,nzb+1
                     t_(:,:,k)=     theta(:,:,k)
                     tx(:,:,k)=      dtdx(:,:,k)
                     ty(:,:,k)=      dtdy(:,:,k)
                     tz(:,:,k)= 0.5*(dtdz(:,:,k)+dtdz(:,:,k+1))
                  enddo
               endif
               if(me.eq.0)then
                  t_(:,:,2)=theta(:,:,2)
                  tx(:,:,2)= dtdx(:,:,2)
                  ty(:,:,2)= dtdy(:,:,2)
                  tz(:,:,2)= dtdz(:,:,2)
               endif
               
               if(scl_nodes.ne.mom_nodes)then
                  
                  if(scl_nodes.eq.0)then
                     do k=2,nzb+1 
                        u_(:,:,k)=0.5*(u(:,:,k)+u(:,:,k-1))
                        v_(:,:,k)=0.5*(v(:,:,k)+v(:,:,k-1))
                        w_(:,:,k)=w(:,:,k)
                        
                        ux(:,:,k)=0.5*(dudx(:,:,k)+dudx(:,:,k-1)) 
                        uy(:,:,k)=0.5*(dudy(:,:,k)+dudy(:,:,k-1)) 
                        uz(:,:,k)=     dudz(:,:,k)
                        vx(:,:,k)=0.5*(dvdx(:,:,k)+dvdx(:,:,k-1))
                        vy(:,:,k)=0.5*(dvdy(:,:,k)+dvdy(:,:,k-1))
                        vz(:,:,k)=     dvdz(:,:,k)
                        wx(:,:,k)=     dwdx(:,:,k)
                        wy(:,:,k)=     dwdy(:,:,k)
                        wz(:,:,k)=0.5*(dwdz(:,:,k)+dwdz(:,:,k-1))
                     end do
                  else
                     do k=2,nzb+1 
                        u_(:,:,k)=u(:,:,k)
                        v_(:,:,k)=v(:,:,k)
                        w_(:,:,k)=0.5*(w(:,:,k)+w(:,:,k+1))
                        
                        ux(:,:,k)=     dudx(:,:,k) 
                        uy(:,:,k)=     dudy(:,:,k) 
                        uz(:,:,k)=0.5*(dudz(:,:,k)+dudz(:,:,k+1))
                        vx(:,:,k)=     dvdx(:,:,k)
                        vy(:,:,k)=     dvdy(:,:,k)
                        vz(:,:,k)=0.5*(dvdz(:,:,k)+dvdz(:,:,k+1))
                        wx(:,:,k)=0.5*(dwdx(:,:,k)+dwdx(:,:,k+1))
                        wy(:,:,k)=0.5*(dwdy(:,:,k)+dwdy(:,:,k+1))
                        wz(:,:,k)=     dwdz(:,:,k)
                     end do
                  endif
                  
                  if (me.eq.0) then           
                     u_(:,:,2)=    u(:,:,2)
                     v_(:,:,2)=    v(:,:,2)
                     w_(:,:,2)=0.5*w(:,:,3)
                     
                     ux(:,:,2)=     dudx(:,:,2)
                     uy(:,:,2)=     dudy(:,:,2) 
                     uz(:,:,2)=     dudz(:,:,2)
                     vx(:,:,2)=     dvdx(:,:,2)
                     vy(:,:,2)=     dvdy(:,:,2)
                     vz(:,:,2)=     dvdz(:,:,2)
                     wx(:,:,2)=0.5*(dwdx(:,:,3)+dwdx(:,:,2))
                     wy(:,:,2)=0.5*(dwdy(:,:,3)+dwdy(:,:,2))
                     wz(:,:,2)=     dwdz(:,:,2)   
                  endif
                  
                  do k=2,Nzb+1
                     S11(:,:,k)=0.5*(ux(:,:,k)+ux(:,:,k))
                     S33(:,:,k)=0.5*(wz(:,:,k)+wz(:,:,k))
                     S22(:,:,k)=0.5*(vy(:,:,k)+vy(:,:,k))
                     S12(:,:,k)=0.5*(uy(:,:,k)+vx(:,:,k))
                     S13(:,:,k)=0.5*(uz(:,:,k)+wx(:,:,k))
                     S23(:,:,k)=0.5*(vz(:,:,k)+wy(:,:,k))
                  end do
                     
               endif

               if(model.eq.2)then

                  if(averaging.eq.1)then
                     call update1(a4_old,me,nall)
                     call update1(b4_old,me,nall)

                     call optim_scl_lag_dyn(Pr2,S11,S33,S22,S12,S13,S23,
     +                    S_hat,u_,v_,w_,L,t,t_,tx,ty,tz,me,
     +                    a4_old,b4_old,0)
                  endif

               elseif(model.eq.3)then
                  
                  if (averaging.eq.1) then

                     call update9(a4_old,b4_old,c4_old,d4_old,e4_old,
     +                    a8_old,b8_old,c8_old,d8_old,me,nall)
                     call update1(e8_old,me,nall)
                     
                     call optim_scl_lag(Pr2,S11,S33,S22,S12,S13,S23,
     +                    S_hat,S_hatd,u_,v_,w_,L,t,t_,tx,ty,tz,me,
     +                    beta2,a4_old,b4_old,c4_old,d4_old,e4_old,
     +                    a8_old,b8_old,c8_old,d8_old,e8_old,0)
                  
                  else if (averaging.eq.0) then
                     
                     call optim_scl_pln(Pr2,S11,S33,S22,S12,S13,S23,
     +                    S_hat,S_hatd,u_,v_,w_,L,t,t_,tx,ty,tz,me,
     +                    beta2)
                  
                  else if (averaging.eq.2) then
                  
c     call optim_scl_loc(Pr2,S11,S33,S22,S12,S13,S23,S_hat,
c     +                 S_hatd,u_,v_,w_,L,t,t_,tx,ty,tz,me,beta2)
                  
                  else if (averaging.eq.3) then
                  
c     call optim_scl_wong(Pr2,S11,S33,S22,S12,S13,S23,S_hat,
c     +                 S_hatd,u_,v_,w_,L,t,t_,tx,ty,tz,me,beta2)
                  
                  end if

               end if

            endif
            
            if(Q_FLAG.eq.1)then

               if(scl_nodes.eq.0)then
                  do k=2,nzb+1
                     q_(:,:,k)=0.5*(q(:,:,k)+q(:,:,k-1))
                     dqx(:,:,k)= 0.5*(dqdx(:,:,k)+dqdx(:,:,k-1))
                     dqy(:,:,k)= 0.5*(dqdy(:,:,k)+dqdy(:,:,k-1))
                     dqz(:,:,k)=      dqdz(:,:,k)
                  enddo
               else
                  do k=2,nzb+1
                     q_(:,:,k) =      q(:,:,k)
                     dqx(:,:,k)=      dqdx(:,:,k)
                     dqy(:,:,k)=      dqdy(:,:,k)
                     dqz(:,:,k)= 0.5*(dqdz(:,:,k)+dqdz(:,:,k+1))
                  enddo
               endif
               if(me.eq.0)then
                  q_(:,:,2)= q(:,:,2)
                  dqx(:,:,2)= dqdx(:,:,2)
                  dqy(:,:,2)= dqdy(:,:,2)
                  dqz(:,:,2)= dqdz(:,:,2)
               endif
                  
               if(model.eq.2)then
                  
                  if(averaging.eq.1)then

                     call update1(a5_old,me,nall)
                     call update1(b5_old,me,nall)

                     call optim_scl_lag_dyn(Sc2,S11,S33,S22,S12,S13,S23,
     +                    S_hat,u_,v_,w_,L,t,q_,dqx,dqy,dqz,me,
     +                    a5_old,b5_old,1)

                  endif
                  
               elseif(model.eq.3)then
                  
                  if (averaging.eq.1) then
                     call update9(a5_old,b5_old,c5_old,d5_old,e5_old,
     +                    a9_old,b9_old,c9_old,d9_old,me,nall)
                     call update1(e9_old,me,nall)
                  end if
                  
                  if (averaging.eq.1) then
                     
                     call optim_scl_lag(Sc2,S11,S33,S22,S12,S13,S23,
     +                    S_hat,S_hatd,u_,v_,w_,L,t,q_,dqx,dqy,dqz,me,
     +                    beta3,a5_old,b5_old,c5_old,d5_old,e5_old,
     +                    a9_old,b9_old,c9_old,d9_old,e9_old,1)
                     
                  else if (averaging.eq.0) then
                     
                     call optim_scl_pln(Sc2,S11,S33,S22,S12,S13,S23,
     +                    S_hat,S_hatd,u_,v_,w_,L,t,q_,dqx,dqy,dqz,me,
     +                    beta3)
                     
                  else if (averaging.eq.2) then
                     
c     call optim_scl_loc(Sc2,S11,S33,S22,S12,S13,S23,S_hat,
c     +                 S_hatd,u_,v_,w_,L,t,q_,dqx,dqy,dqz,me,beta3)
                     
                  else if (averaging.eq.3) then
                     
c     call optim_scl_wong(Sc2,S11,S33,S22,S12,S13,S23,S_hat,
c     +                 S_hatd,u_,v_,w_,L,t,q_,dqx,dqy,dqz,me,beta3)
                     
                  end if
               
               endif
            
            endif
            
c     ... Richardson number criteria***************
            if (Ri_flag.eq.1) then
               if (pass_flag.eq.0.and.S_flag.eq.1) then
                  do k=2,Nzb+1
                     do j=1,Ny
                        do i=1,Nx
                           Rip(i,j,k)=(g_hat/(Theta_0/T_scale))*
     +                          tz(i,j,k)/(S(i,j,k)**2.)
                           if(Rip(i,j,k).gt.(Cs2(i,j,k)/Pr2(i,j,k)))then 
                              Cs2(i,j,k)=0.
                              Pr2(i,j,k)=0.
                              if(q_flag.eq.1) Sc2(i,j,k)=0.
                           end if               
                        end do
                     end do
                  end do
               end if
            end if
c     *******************************************
            call dealias1(Cs2,Cs2_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
            Cs2_m=abs(Cs2_m)
            call update1_m(Cs2_m,me,nall)
            if(S_flag.eq.1)then
               call dealias1(Pr2,Pr2_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
			   Pr2_m=abs(Pr2_m)
               call update1_m(Pr2_m,me,nall)
            endif
            if(Q_flag.eq.1)then
               call dealias1(Sc2,Sc2_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
			   Sc2_m=abs(Sc2_m)
               call update1_m(Sc2_m,me,nall)
            endif
            
            cst=1

cc            do k=2,nzb+1
c               cs_ave=0
c               pr_ave=0
c               sc_ave=0
cc               do j=1,ny2
cc                  do i=1,nx2
cc				      if (Cs2_m(i,j,k).lt.0) Cs2_m(i,j,k)=0.				  
cc					  if (S_flag.eq.1 .and. Pr2_m(i,j,k).lt.0) Pr2_m(i,j,k)=0.
cc					  if (Q_flag.eq.1 .and. Sc2_m(i,j,k).lt.0) Sc2_m(i,j,k)=0.
c                     cs_ave=Cs2(i,j,k)+cs_ave*inxny
c                     pr_ave=Pr2(i,j,k)+Pr_ave*inxny
c                     sc_ave=Sc2(i,j,k)+sc_ave*inxny
cc                  enddo
cc               enddo
c               write(*,*) me*nzb+k-1,sqrt(cs_ave),cs_ave/pr_ave,
c     +              cs_ave/sc_ave
cc            enddo

         else
            
            cst=cst+1
            
         end if
         
      else
         
c     ... Richardson number criteria***************
         if (Ri_flag.eq.1) then
            if (pass_flag.eq.0.and.S_flag.eq.1) then
               do k=2,Nzb+1
                  do j=1,Ny
                     do i=1,Nx
                        Rip(i,j,k)=(g_hat/(Theta_0/T_scale))*
     +                       tz(i,j,k)/(S(i,j,k)**2.)
                        if(Rip(i,j,k).gt.(Cs2(i,j,k)/Pr2(i,j,k))) then 
                           Cs2(i,j,k)=0.
                           Pr2(i,j,k)=0.
                           if(q_flag.eq.1) Sc2(i,j,k)=0.
                        end if               
                     end do
                  end do
               end do
            end if
         end if
c     *******************************************
         
      endif
      
      if(mom_nodes.eq.1)then
         do k=2,nzb+1 
            ux(:,:,k)=0.5*(dudx(:,:,k)+dudx(:,:,k-1)) 
            uy(:,:,k)=0.5*(dudy(:,:,k)+dudy(:,:,k-1)) 
            uz(:,:,k)=     dudz(:,:,k)
            vx(:,:,k)=0.5*(dvdx(:,:,k)+dvdx(:,:,k-1))
            vy(:,:,k)=0.5*(dvdy(:,:,k)+dvdy(:,:,k-1))
            vz(:,:,k)=     dvdz(:,:,k)
            wx(:,:,k)=     dwdx(:,:,k)
            wy(:,:,k)=     dwdy(:,:,k)
            wz(:,:,k)=0.5*(dwdz(:,:,k)+dwdz(:,:,k-1))
         end do
      else
         do k=2,nzb+1 
            ux(:,:,k)=     dudx(:,:,k) 
            uy(:,:,k)=     dudy(:,:,k) 
            uz(:,:,k)=0.5*(dudz(:,:,k)+dudz(:,:,k+1))
            vx(:,:,k)=     dvdx(:,:,k)
            vy(:,:,k)=     dvdy(:,:,k)
            vz(:,:,k)=0.5*(dvdz(:,:,k)+dvdz(:,:,k+1))
            wx(:,:,k)=0.5*(dwdx(:,:,k)+dwdx(:,:,k+1))
            wy(:,:,k)=0.5*(dwdy(:,:,k)+dwdy(:,:,k+1))
            wz(:,:,k)=     dwdz(:,:,k)
         end do
      endif
      
      if (me==0) then 
         ux(:,:,2)=     dudx(:,:,2)
         uy(:,:,2)=     dudy(:,:,2) 
         uz(:,:,2)=     dudz(:,:,2)
         vx(:,:,2)=     dvdx(:,:,2)
         vy(:,:,2)=     dvdy(:,:,2)
         vz(:,:,2)=     dvdz(:,:,2)
         wx(:,:,2)=0.5*(dwdx(:,:,3)+dwdx(:,:,2))
         wy(:,:,2)=0.5*(dwdy(:,:,3)+dwdy(:,:,2))
         wz(:,:,2)=     dwdz(:,:,2)   
      endif
      
      do k=2,Nzb+1
         S11(:,:,k)=0.5*(ux(:,:,k)+ux(:,:,k))
         S33(:,:,k)=0.5*(wz(:,:,k)+wz(:,:,k))
         S22(:,:,k)=0.5*(vy(:,:,k)+vy(:,:,k))
         S12(:,:,k)=0.5*(uy(:,:,k)+vx(:,:,k))
         S13(:,:,k)=0.5*(uz(:,:,k)+wx(:,:,k))
         S23(:,:,k)=0.5*(vz(:,:,k)+wy(:,:,k))
      end do
      
      if(mom_nodes.eq.1)then
         
         call dealias1(S11,S11_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(S12,S12_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(S13,S13_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(S22,S22_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(S23,S23_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(S33,S33_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         
         do k=2,nzb+1
            S_m(:,:,k)=sqrt(2*(S11_m(:,:,k)**2+S22_m(:,:,k)**2+
     +           S33_m(:,:,k)**2+2.*S12_m(:,:,k)**2+2.*S13_m(:,:,k)**2+
     +           2.*S23_m(:,:,k)**2))
         enddo
         
      else
         
         call dealias1(S11,S11_mu,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(S12,S12_mu,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(S13,S13_mu,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(S22,S22_mu,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(S23,S23_mu,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(S33,S33_mu,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         
         do k=2,nzb+1
            S_mu(:,:,k)=sqrt(2*(S11_mu(:,:,k)**2+S22_mu(:,:,k)**2+
     +           S33_mu(:,:,k)**2+2.*S12_mu(:,:,k)**2+
     +           2.*S13_mu(:,:,k)**2+2.*S23_mu(:,:,k)**2))
         enddo
         
      endif
      
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(S_flag.eq.1)then
         call dealias1(dtdx,tx_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(dtdy,ty_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(dtdz,tz_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      endif
      if(Q_flag.eq.1)then
         call dealias1(dqdx,dqx_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(dqdy,dqy_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(dqdz,dqz_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      do k=2,Nzb+1

         if(me.eq.0.and.k.eq.2)then
            
c     ... Wong-Lilly Model************************
            if(averaging.eq.3) then
               
               txx_m(:,:,k)=-2.*Cs2_m(:,:,k)*
     +              l(k)**(4./3.)*S11_mu(:,:,k)
               tyy_m(:,:,k)=-2.*Cs2_m(:,:,k)*
     +              l(k)**(4./3.)*S22_mu(:,:,k)
               tzz_m(:,:,k)=-2.*Cs2_m(:,:,k)*
     +              l(k)**(4./3.)*S33_mu(:,:,k)
               txy_m(:,:,k)=-2.*Cs2_m(:,:,k)*
     +              l(k)**(4./3.)*S12_mu(:,:,k)
c     ... Smagorinsky Model***********************
            else 
               txx_m(:,:,k)=-2.*Cs2_m(:,:,k)*
     +              l(k)**2.*S_mu(:,:,k)*S11_mu(:,:,k)
               tyy_m(:,:,k)=-2.*Cs2_m(:,:,k)*
     +              l(k)**2.*S_mu(:,:,k)*S22_mu(:,:,k)
               tzz_m(:,:,k)=-2.*Cs2_m(:,:,k)*
     +              l(k)**2.*S_mu(:,:,k)*S33_mu(:,:,k)
               txy_m(:,:,k)=-2.*Cs2_m(:,:,k)*
     +              l(k)**2.*S_mu(:,:,k)*S12_mu(:,:,k)
               
            end if
            
            if(S_flag.eq.1)then
               
c     ... Wong-Lilly Model************************
               if(averaging.eq.3) then
                  
                  qx_m(:,:,k)=-(Pr2_m(:,:,k))*
     +                 (l(k)**(4./3.))*tx_m(:,:,k)
                  qy_m(:,:,k)=-(Pr2_m(:,:,k))*
     +                 (l(k)**(4./3.))*ty_m(:,:,k)

c     ... Smagorinsky Model***********************
               else              

                  qx_m(:,:,k)=-(Pr2_m(:,:,k))*
     +                 (l(k)**2.)*(S_mu(:,:,k)*tx_m(:,:,k))
                  qy_m(:,:,k)=-(Pr2_m(:,:,k))*
     +                 (l(k)**2.)*(S_mu(:,:,k)*ty_m(:,:,k))
               end if
               
            endif
            if(Q_flag.eq.1)then
               
c     ... Wong-Lilly Model************************
               if(averaging.eq.3) then
                  
                  sgs_q1_m(:,:,k)=-(Sc2_m(:,:,k))*
     +                 (l(k)**(4./3.))*dqx_m(:,:,k)
                  sgs_q2_m(:,:,k)=-(Sc2_m(:,:,k))*
     +                 (l(k)**(4./3.))*dqy_m(:,:,k)

c     ... Smagorinsky Model***********************
               else              

                  sgs_q1_m(:,:,k)=-(Sc2_m(:,:,k))*
     +                 (l(k)**2.)*(S_mu(:,:,k)*dqx_m(:,:,k))
                  sgs_q2_m(:,:,k)=-(Sc2_m(:,:,k))*
     +                 (l(k)**2.)*(S_mu(:,:,k)*dqy_m(:,:,k))
               end if
               
            endif
            
         else

            if(mom_nodes.eq.0)then
               Cs2_L=0.5*(Cs2_m(:,:,k)+Cs2_m(:,:,k+1))
            else
               Cs2_L=Cs2_m(:,:,k)
            endif

c     ... Wong-Lilly Model************************
            if(averaging.eq.3) then

               txx_m(:,:,k)=-2.*Cs2_L*
     +              0.5*(l(k)**(4./3.)+l(k+1)**(4./3.))*
     +              S11_mu(:,:,k)
               tyy_m(:,:,k)=-2.*Cs2_L*
     +              0.5*(l(k)**(4./3.)+l(k+1)**(4./3.))*
     +              S22_mu(:,:,k)
               tzz_m(:,:,k)=-2.*Cs2_L*
     +              0.5*(l(k)**(4./3.)+l(k+1)**(4./3.))*
     +              S33_mu(:,:,k)
               txy_m(:,:,k)=-2.*Cs2_L*
     +              0.5*(l(k)**(4./3.)+l(k+1)**(4./3.))*
     +              S12_mu(:,:,k)
               
c     ... Smagorinsky Model***********************
            else
               txx_m(:,:,k)=-2.*Cs2_L*
     +              0.5*(l(k)**2.+l(k+1)**2.)*
     +              S_mu(:,:,k)*S11_mu(:,:,k)
               tyy_m(:,:,k)=-2.*Cs2_L*
     +              0.5*(l(k)**2.+l(k+1)**2.)*
     +              S_mu(:,:,k)*S22_mu(:,:,k)
               tzz_m(:,:,k)=-2.*Cs2_L*
     +              0.5*(l(k)**2.+l(k+1)**2.)*
     +              S_mu(:,:,k)*S33_mu(:,:,k)
               txy_m(:,:,k)=-2.*Cs2_L*
     +              0.5*(l(k)**2.+l(k+1)**2.)*
     +              S_mu(:,:,k)*S12_mu(:,:,k)
               
            end if
            
            if(S_flag.eq.1)then
               if(scl_nodes.eq.0)then
                  Pr2_L=0.5*(Pr2_m(:,:,k)+Pr2_m(:,:,k+1))
               else
                  Pr2_L=Pr2_m(:,:,k)
               endif

c     ... Wong-Lilly Model************************
               if(averaging.eq.3) then

                  qx_m(:,:,k)=-Pr2_L*
     +                 0.5*(l(k)**(4./3.)+l(k+1)**(4./3.))*
     +                 tx_m(:,:,k)
                  qy_m(:,:,k)=-Pr2_L*
     +                 0.5*(l(k)**(4./3.)+l(k+1)**(4./3.))*
     +                 ty_m(:,:,k)
                  
c     ... Smagorinsky Model***********************
               else
                  qx_m(:,:,k)=-Pr2_L*
     +                 0.5*(l(k)**2.+l(k+1)**2.)*
     +                 S_mu(:,:,k)*tx_m(:,:,k)
                  qy_m(:,:,k)=-Pr2_L*
     +                 0.5*(l(k)**2.+l(k+1)**2.)*
     +                 S_mu(:,:,k)*ty_m(:,:,k)
               endif
               
            end if
            
            if(Q_flag.eq.1)then
               if(scl_nodes.eq.0)then
                  Sc2_L=0.5*(Sc2_m(:,:,k)+Sc2_m(:,:,k+1))
               else
                  Sc2_L=Sc2_m(:,:,k)
               endif
               
c     ... Wong-Lilly Model************************
               if(averaging.eq.3) then
                  
                  sgs_q1_m(:,:,k)=-Sc2_L*
     +                 0.5*(l(k)**(4./3.)+l(k+1)**(4./3.))*
     +                 dqx_m(:,:,k)
                  sgs_q2_m(:,:,k)=-Sc2_L*
     +                 0.5*(l(k)**(4./3.)+l(k+1)**(4./3.))*
     +                 dqy_m(:,:,k)
                  
c     ... Smagorinsky Model***********************
               else
                  sgs_q1_m(:,:,k)=-Sc2_L*
     +                 0.5*(l(k)**2.+l(k+1)**2.)*
     +                 S_mu(:,:,k)*dqx_m(:,:,k)
                  sgs_q2_m(:,:,k)=-Sc2_L*
     +                 0.5*(l(k)**2.+l(k+1)**2.)*
     +                 S_mu(:,:,k)*dqy_m(:,:,k)
               endif
               
            end if
            
         endif
         
      end do
      
      do k=2,nzb+1
         
         if(mom_nodes.eq.1)then
            Cs2_L=0.5*(Cs2_m(:,:,k)+Cs2_m(:,:,k-1))
         else
            Cs2_L=Cs2_m(:,:,k)
         endif

c     ... Wong-Lilly Model************************
         if(averaging.eq.3) then

            txz_m(:,:,k)=-2.*Cs2_L*l(k)**(4./3.)*
     +           S13_m(:,:,k)
            tyz_m(:,:,k)=-2.*Cs2_L*l(k)**(4./3.)*
     +           S23_m(:,:,k)
            
c     ... Smagorinsky Model***********************
         else
            txz_m(:,:,k)=-2.*Cs2_L*l(k)**2.*
     +           S_m(:,:,k)*S13_m(:,:,k)
            tyz_m(:,:,k)=-2.*Cs2_L*l(k)**2.*
     +           S_m(:,:,k)*S23_m(:,:,k)

         end if
         
         if(S_flag.eq.1)then
            if(scl_nodes.eq.1)then
               Pr2_L=0.5*(Pr2_m(:,:,k)+Pr2_m(:,:,k-1))
            else
               Pr2_L=Pr2_m(:,:,k)
            endif
            if (me==0.and.k.eq.2)then

c     ... Wong-Lilly Model************************
               if(averaging.eq.3) then
                  qz_m(:,:,k)=-Pr2_m(:,:,k)*l(k)**(4./3.)*
     +                 tz_m(:,:,k)
c     ... Smagorinsky Model***********************
               else
                  qz_m(:,:,k)=-Pr2_m(:,:,k)*l(k)**2*
     +                 S_mu(:,:,k)*tz_m(:,:,k)
               end if
            else
c     ... Wong-Lilly Model************************
               if(averaging.eq.3) then
                  qz_m(:,:,k)=-Pr2_L*(l(k)**(4./3.))*
     +                 tz_m(:,:,k)
c     ... Smagorinsky Model***********************
               else
                  qz_m(:,:,k)=-Pr2_L*(l(k)**2.)*
     +                 (S_m(:,:,k)*tz_m(:,:,k))
               end if
            end if
         endif   
         if(Q_flag.eq.1)then
            if(scl_nodes.eq.1)then
               Sc2_L=0.5*(Sc2_m(:,:,k)+Sc2_m(:,:,k-1))
            else
               Sc2_L=Sc2_m(:,:,k)
            endif
            if (me==0.and.k.eq.2)then
               
c     ... Wong-Lilly Model************************
               if(averaging.eq.3) then
                  sgs_q3_m(:,:,k)=-Sc2_m(:,:,k)*l(k)**(4./3.)*
     +                 dqz_m(:,:,k)
c     ... Smagorinsky Model***********************
               else
                  sgs_q3_m(:,:,k)=-Sc2_m(:,:,k)*l(k)**2*
     +                 S_mu(:,:,k)*dqz_m(:,:,k)
               end if
            else
c     ... Wong-Lilly Model************************
               if(averaging.eq.3) then
                  sgs_q3_m(:,:,k)=-Sc2_L*(l(k)**(4./3.))*
     +                 dqz_m(:,:,k)
c     ... Smagorinsky Model***********************
               else
                  sgs_q3_m(:,:,k)=-Sc2_L*(l(k)**2.)*
     +                 (S_m(:,:,k)*dqz_m(:,:,k))
               end if
            end if
            
         endif
         
      end do
         
      if (me==nprocs-1) then
         do k=(nzb+1),(nzb+2)
            txz_m(:,:,k)=0.
            tyz_m(:,:,k)=0.
            tzz_m(:,:,k)=0.
            txx_m(:,:,k)=0.
            txy_m(:,:,k)=0.
            tyy_m(:,:,k)=0.
            if(S_flag.eq.1)then
               qz_m(:,:,k)=0.
            endif
            if(Q_flag.eq.1)then
               sgs_q3_m(:,:,k)=0.
            endif
         end do
      endif
      
      call dealias2(txx,txx_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)
      call dealias2(txy,txy_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)
      call dealias2(txz,txz_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)
      call dealias2(tyy,tyy_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)
      call dealias2(tyz,tyz_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)
      if(Q_flag.eq.1)then
         call dealias2(sgs_q1,sgs_q1_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)
         call dealias2(sgs_q2,sgs_q2_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)
         call dealias2(sgs_q3,sgs_q3_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)
      endif
      
      if(S_flag.eq.1)then
         call dealias2(qx,qx_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)
         call dealias2(qy,qy_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)
         call dealias2(qz,qz_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)
         call dealias2(tzz,tzz_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)
      else
         call dealias2(tzz,tzz_m,t,2,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)
      endif
      
      call update3(txz,tyz,tzz,me,nall)
      
      if (me==0) then
         txz(:,:,2)=txzp
         tyz(:,:,2)=tyzp
      end if

c     SGS Kinetic Energy Computation**************************
c     NOTE: THIS IS ONLY CORRECT IF scalar and momentum are on the same nodes
               do k=2,nzb+1
                  do j=1,ny
                     do i=1,nx
                        if (averaging.eq.3) then
                           ESGS3D(i,j,k)=l(k)**(4./3.)*S(i,j,k)*
     +                          Cs2(i,j,k)/0.3
                        else
                           ESGS3D(i,j,k)=((l(k)*S(i,j,k))**2.0)*
     +                          Cs2(i,j,k)/0.3d0
                        end if
                     end do
                  end do
               end do
c     ********************************************************

c     SGS Scalar Variance computation*****************      
      
      if (S_flag.eq.1) then
         call update3(qx,qy,qz,me,nall)
         do k=2,nzb+1
            do j=1,ny
               do i=1,nx
                  ET3D(i,j,k)=sqrt((0.5*(qx(i,j,k)+qx(i,j,k-1)))**2.+
     +                 (0.5*(qy(i,j,k)+qy(i,j,k-1)))**2.+qz(i,j,k)**2.)
               end do
            end do
         end do
      end if

      if (Q_flag.eq.1) then
         call update3(sgs_q1,sgs_q2,sgs_q3,me,nall)
         do k=2,nzb+1
            do j=1,ny
               do i=1,nx
                  EQ3D(i,j,k)=sqrt((0.5*(sgs_q1(i,j,k)+
     +                 sgs_q1(i,j,k-1)))**2.+(0.5*(sgs_q2(i,j,k)+
     +                 sgs_q2(i,j,k-1)))**2.+sgs_q3(i,j,k)**2.)
               end do
            end do
         end do
      end if
c     ************************************************

 5333 format (i5,501(1x,e11.4))
 5166 format (1x,i5,100(1x,f10.6))
 6623 format (100(f6.3,1x))
 6624 format (100(f8.5,1x))
      
 300  return					
      
      end









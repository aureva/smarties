      Subroutine optim_scl_lag_dyn(Pr2,S11,S33,S22,S12,S13,S23,S_hat,
     +     u_,v_,w_,L,t,t_,tx,ty,tz,me,KX_old,XX_old,flg)
      
      implicit none
      include 'dimen.h'
      integer*4 i,j,k,t,m,flg
      real*8,dimension(nx,ny,nz2):: tz,tx,ty,tx_hat,ty_hat,
     +     tz_hat,S,S11,S22,S33,S12,S13,S23,S_hat,S11_hat,S22_hat,
     +     S33_hat,S12_hat,S13_hat,S23_hat,u_hat,v_hat,w_hat,t_,t_hat,
     +     w_,u_,v_,s_tx,s_ty,s_tz,ut,vt,wt,ut_hat,vt_hat,wt_hat,
     +     s_tx_hat,s_ty_hat,s_tz_hat

      real*8,dimension(nx,ny,nz2):: Pr2,KX,XX,KX_old,XX_old

      real*8,dimension(nz2):: KX_avg,XX_avg,sig_t,avg_t

      real*8 a2,b2,c2,d2,e2,Tn,b_ave,l(nz2),
     +     pr_ave,uin,vin,win
     
      do k=2,nzb+1
         do j=1,ny
            do i=1,nx

               S(i,j,k)=sqrt(2.*(S11(i,j,k)**2.+S22(i,j,k)**2.+
     +              S33(i,j,k)**2.+2.*S12(i,j,k)**2.+
     +              2.*S13(i,j,k)**2.+2.*S23(i,j,k)**2.))

               S_tx(i,j,k)=S(i,j,k)*tx(i,j,k)
               S_ty(i,j,k)=S(i,j,k)*ty(i,j,k)
               S_tz(i,j,k)=S(i,j,k)*tz(i,j,k)
               ut(i,j,k)=u_(i,j,k)*t_(i,j,k)
               vt(i,j,k)=v_(i,j,k)*t_(i,j,k)
               wt(i,j,k)=w_(i,j,k)*t_(i,j,k)

            enddo
         enddo
      enddo

CC...Filtering to get the _hat variables (coarser resolution)

      Call Filter_La(u_hat,u_,t,0)
      Call Filter_La(v_hat,v_,t,0)
      Call Filter_La(w_hat,w_,t,0)
      Call Filter_La(t_hat,t_,t,0)
      Call Filter_La(ut_hat,ut,t,0)
      Call Filter_La(vt_hat,vt,t,0)
      Call Filter_La(wt_hat,wt,t,0)
      Call Filter_La(tx_hat,tx,t,0)
      Call Filter_La(ty_hat,ty,t,0)
      Call Filter_La(tz_hat,tz,t,0)
      
      if(scl_nodes.ne.mom_nodes)then

         Call Filter_La(S11_hat,S11,t,0)
         Call Filter_La(S22_hat,S22,t,0)
         Call Filter_La(S33_hat,S33,t,0)
         Call Filter_La(S12_hat,S12,t,0)
         Call Filter_La(S13_hat,S13,t,0)
         Call Filter_La(S23_hat,S23,t,0)

         do k=2,nzb+1
            do j=1,ny
               do i=1,nx
                  
                  S_hat(i,j,k)=sqrt(2.*(S11_hat(i,j,k)**2.+
     +                 S22_hat(i,j,k)**2.+
     +                 S33_hat(i,j,k)**2.+
     +                 2.*S12_hat(i,j,k)**2.+
     +                 2.*S13_hat(i,j,k)**2.+
     +                 2.*S23_hat(i,j,k)**2.))

               end do
            end do
         end do

      endif

      Call Filter_La(s_tx_hat,s_tx,t,0)
      Call Filter_La(s_ty_hat,s_ty,t,0)
      if(flg.eq.0.and.Q_flag.eq.1)then
         Call Filter_La(s_tz_hat,s_tz,t,0)
      else
         Call Filter_La(s_tz_hat,s_tz,t,2)
      endif

      do k=2,nzb+1
         
         do j=1,ny
            do i=1,nx

               a2=+L(k)**2.*
     +              ((ut_hat(i,j,k)-u_hat(i,j,k)*t_hat(i,j,k))*
     +              s_tx_hat(i,j,k)+(vt_hat(i,j,k)-v_hat(i,j,k)*
     +              t_hat(i,j,k))*s_ty_hat(i,j,k)+(wt_hat(i,j,k)-
     +              w_hat(i,j,k)*t_hat(i,j,k))*s_tz_hat(i,j,k))

               b2=-(tfr*L(k))**2.*((ut_hat(i,j,k)-u_hat(i,j,k)*
     +              t_hat(i,j,k))*s_hat(i,j,k)*tx_hat(i,j,k)
     +              +(vt_hat(i,j,k)-v_hat(i,j,k)*t_hat(i,j,k))*
     +              s_hat(i,j,k)*ty_hat(i,j,k)+(wt_hat(i,j,k)-
     +              w_hat(i,j,k)*t_hat(i,j,k))*s_hat(i,j,k)*
     +              tz_hat(i,j,k))

               c2=+L(k)**4.*((s_tx_hat(i,j,k))**2.
     +              +(s_ty_hat(i,j,k))**2.+(s_tz_hat(i,j,k))**2.)

               d2=-2.*L(k)**4.*(tfr**2.)*
     +              (s_tx_hat(i,j,k)*s_hat(i,j,k)*tx_hat(i,j,k)
     +              +s_ty_hat(i,j,k)*s_hat(i,j,k)*ty_hat(i,j,k)
     +              +s_tz_hat(i,j,k)*s_hat(i,j,k)*tz_hat(i,j,k))

               e2=+(tfr*L(k))**4.*
     +              ((s_hat(i,j,k)*tx_hat(i,j,k))**2.
     +              +(s_hat(i,j,k)*ty_hat(i,j,k))**2.
     +              +(s_hat(i,j,k)*tz_hat(i,j,k))**2.)

               KX(i,j,k)=a2+b2
               XX(i,j,k)=c2+d2+e2

            end do
         end do
      end do

c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c     **************initilization for lagrangian averaging***********   
      if(t.le.50.and.INITS.eq.0)then
         do k=2,nzb+1
            KX_avg(k)=0.
            XX_avg(k)=0.
            
            do j=1,Ny
               do i=1,Nx
                  
                  KX_avg(k)=KX_avg(k)+KX(i,j,k)/(Nx*Ny)
                  XX_avg(k)=XX_avg(k)+XX(i,j,k)/(Nx*Ny)
                  
               enddo
            enddo
         enddo
         
         do k=2,nzb+1
            do j=1,Ny
               do i=1,Nx
                  
                  KX(i,j,k)=KX_avg(k)
                  XX(i,j,k)=XX_avg(k)
                  
                  KX_old(i,j,k)=KX_avg(k)
                  XX_old(i,j,k)=XX_avg(k)
                  
               end do
            end do
         end do
      endif
c     **************end of the lagrangian initilization**************
  
      do k=2,nzb+1
         avg_t(k)=0.0
         do j=1,ny
            do i=1,nx
               avg_t(k)=avg_t(k)+t_(i,j,k)
            enddo
         enddo
         avg_t(k)=avg_t(k)*inxny
      enddo

      do k=2,nzb+1
         sig_t(k)=0.0
         do j=1,ny
            do i=1,nx
               sig_t(k)=sig_t(k)+(t_(i,j,k)-avg_t(k))**2.
            enddo
         enddo
         sig_t(k)=(sig_t(k)/dble(nx*ny-1))**0.5
      enddo

      do k=2,nzb+1

         do j=1,ny
            do i=1,nx
               
               if(me.eq.0.and.k.le.2)then
                  uin = u_(i,j,k)
                  vin = v_(i,j,k)
                  win = w_hat(i,j,k)
               else
                  uin = u_(i,j,k)
                  vin = v_(i,j,k)
                  win = w_hat(i,j,k)
               endif

               call lagrng_dyn_s(KX(i,j,k),XX(i,j,k),KX_old,XX_old,
     +              uin,vin,win,i,j,k,sig_t(k),me,t)
               
               if(KX(i,j,k).lt.0.or.XX(i,j,k).lt.(1e-10))then

                  Pr2(i,j,k)=0.0
                  
                  KX(i,j,k)=0.0
c     XX(i,j,k)=0.0
                  
               else

c     Note: Pr is actually Cs2/Prt
                  
                  Pr2(i,j,k)=KX(i,j,k)/XX(i,j,k)
                  
               endif
               
            enddo
         enddo
      enddo

      do k=2,nzb+1
         do j=1,ny
            do i=1,nx

               KX_old(i,j,k)=KX(i,j,k)
               XX_old(i,j,k)=XX(i,j,k)

            enddo
         enddo
      enddo
       
c     **************end of betaa calculation*************************

      return
      end












      Subroutine optim_scl_loc(Pr2,S11,S33,S22,S12,S13,S23,S_hat,S_hatd,
     +     u_,v_,w_,L,t,t_,tx,ty,tz,me,betaa1,bad_pr1)
      
      implicit none
      include 'dimen.h'
      integer*4 i,j,k,t,ii,jj,mX,mY
      real*8,dimension(nx,ny,nz2):: tz,tx,ty,tx_hat,ty_hat,
     +     tz_hat,S,S11,S22,S33,S12,S13,S23,S_hat,S11_hat,S22_hat,
     +     S33_hat,S12_hat,S13_hat,S23_hat,u_hat,v_hat,w_hat,t_,t_hat,
     +     S_hatd,S11_hatd,S22_hatd,S33_hatd,S12_hatd,S13_hatd,
     +     S23_hatd,u_hatd,v_hatd,w_hatd,w_,u_,v_,betaa1,s_tx,s_ty,s_tz,
     +     ut,vt,wt,ut_hat,vt_hat,wt_hat,ut_hatd,vt_hatd,wt_hatd,t_hatd,
     +     tx_hatd,ty_hatd,tz_hatd,s_tx_hat,s_ty_hat,s_tz_hat,
     +     s_tx_hatd,s_ty_hatd,s_tz_hatd,Pr2

      real*8 a2,b2,c2,d2,e2,a4,b4,c4,d4,e4,betaa1D(nz2)

      real*8 T_up,T_down,aa,bb,cc,dd,ee,ff,rtnewt,l(nz2),bad_pr1(Nz2)
     
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

      Call Filter_2Laa(u_hat,u_hatd,u_,t,0)
      Call Filter_2Laa(v_hat,v_hatd,v_,t,0)
      Call Filter_2Laa(w_hat,w_hatd,w_,t,0)
      Call Filter_2Laa(t_hat,t_hatd,t_,t,0)
      Call Filter_2Laa(ut_hat,ut_hatd,ut,t,0)
      Call Filter_2Laa(vt_hat,vt_hatd,vt,t,0)
      Call Filter_2Laa(wt_hat,wt_hatd,wt,t,0)
      Call Filter_2Laa(tx_hat,tx_hatd,tx,t,0)
      Call Filter_2Laa(ty_hat,ty_hatd,ty,t,0)
      Call Filter_2Laa(tz_hat,tz_hatd,tz,t,0)
      
      if(scl_nodes.ne.mom_nodes)then

         Call Filter_2Laa(S11_hat,S11_hatd,S11,t,0)
         Call Filter_2Laa(S22_hat,S22_hatd,S22,t,0)
         Call Filter_2Laa(S33_hat,S33_hatd,S33,t,0)
         Call Filter_2Laa(S12_hat,S12_hatd,S12,t,0)
         Call Filter_2Laa(S13_hat,S13_hatd,S13,t,0)
         Call Filter_2Laa(S23_hat,S23_hatd,S23,t,0)

         do k=2,nzb+1
            do j=1,ny
               do i=1,nx
                  
                  S_hat(i,j,k)=sqrt(2.*(S11_hat(i,j,k)**2.+
     +                 S22_hat(i,j,k)**2.+
     +                 S33_hat(i,j,k)**2.+
     +                 2.*S12_hat(i,j,k)**2.+
     +                 2.*S13_hat(i,j,k)**2.+
     +                 2.*S23_hat(i,j,k)**2.))

                  S_hatd(i,j,k)=sqrt(2.*(S11_hatd(i,j,k)**2.
     +                 +S22_hatd(i,j,k)**2.
     +                 +S33_hatd(i,j,k)**2.
     +                 +2.*S12_hatd(i,j,k)**2.
     +                 +2.*S13_hatd(i,j,k)**2.
     +                 +2.*S23_hatd(i,j,k)**2.))

               end do
            end do
         end do

      endif

      Call Filter_2Laa(s_tx_hat,s_tx_hatd,s_tx,t,0)
      Call Filter_2Laa(s_ty_hat,s_ty_hatd,s_ty,t,0)
      Call Filter_2Laa(s_tz_hat,s_tz_hatd,s_tz,t,2)

      do k=2,nzb+1                  
         a2=0.
         a4=0.
         b2=0.
         b4=0.
         c2=0.
         c4=0.
         d2=0.
         d4=0.
         e2=0.
         e4=0.

         do j=1,ny
            do i=1,nx

               a2=a2+L(k)**2.*
     +              ((ut_hat(i,j,k)-
     +              u_hat(i,j,k)*t_hat(i,j,k))*s_tx_hat(i,j,k)+
     +              (vt_hat(i,j,k)-
     +              v_hat(i,j,k)*t_hat(i,j,k))*s_ty_hat(i,j,k)+
     +              (wt_hat(i,j,k)-
     +              w_hat(i,j,k)*t_hat(i,j,k))*s_tz_hat(i,j,k))
               
               b2=b2-(L(k)**2.)*tfr**2.*
     +              ((ut_hat(i,j,k)-u_hat(i,j,k)*t_hat(i,j,k))*
     +              s_hat(i,j,k)*tx_hat(i,j,k)
     +              +(vt_hat(i,j,k)-v_hat(i,j,k)*t_hat(i,j,k))*
     +              s_hat(i,j,k)*ty_hat(i,j,k)
     +              +(wt_hat(i,j,k)-w_hat(i,j,k)*t_hat(i,j,k))*
     +              s_hat(i,j,k)*tz_hat(i,j,k))
               
               c2=c2+L(k)**4.*
     +              ((s_tx_hat(i,j,k))**2.
     +              +(s_ty_hat(i,j,k))**2.
     +              +(s_tz_hat(i,j,k))**2.)
               
               d2=d2-(2.*L(k)**4.)*(tfr**2.)*
     +              (s_tx_hat(i,j,k)*s_hat(i,j,k)*tx_hat(i,j,k)
     +              +s_ty_hat(i,j,k)*s_hat(i,j,k)*ty_hat(i,j,k)
     +              +s_tz_hat(i,j,k)*s_hat(i,j,k)*tz_hat(i,j,k))
               
               e2=e2+(L(k)**4.)*(tfr**4.)*
     +              ((s_hat(i,j,k)*tx_hat(i,j,k))**2.
     +              +(s_hat(i,j,k)*ty_hat(i,j,k))**2.
     +              +(s_hat(i,j,k)*tz_hat(i,j,k))**2.)
               
               a4=a4+L(k)**2.*
     +              ((ut_hatd(i,j,k)-
     +              u_hatd(i,j,k)*t_hatd(i,j,k))*s_tx_hatd(i,j,k)
     +              +(vt_hatd(i,j,k)-
     +              v_hatd(i,j,k)*t_hatd(i,j,k))*s_ty_hatd(i,j,k)
     +              +(wt_hatd(i,j,k)-
     +              w_hatd(i,j,k)*t_hatd(i,j,k))*s_tz_hatd(i,j,k))
               
               b4=b4-(L(k)**2.)*(tfr**4.)*
     +              ((ut_hatd(i,j,k)-u_hatd(i,j,k)*t_hatd(i,j,k))*
     +              s_hatd(i,j,k)*tx_hatd(i,j,k)
     +              +(vt_hatd(i,j,k)-v_hatd(i,j,k)*t_hatd(i,j,k))*
     +              s_hatd(i,j,k)*ty_hatd(i,j,k)
     +              +(wt_hatd(i,j,k)-w_hatd(i,j,k)*t_hatd(i,j,k))*
     +              s_hatd(i,j,k)*tz_hatd(i,j,k))
               
               c4=c4+L(k)**4.*
     +              ((s_tx_hatd(i,j,k))**2.
     +              +(s_ty_hatd(i,j,k))**2.
     +              +(s_tz_hatd(i,j,k))**2.)
               
               d4=d4-(2.*L(k)**4.)*(tfr**4.)*
     +              (s_tx_hatd(i,j,k)*s_hatd(i,j,k)*tx_hatd(i,j,k)
     +              +s_ty_hatd(i,j,k)*s_hatd(i,j,k)*ty_hatd(i,j,k)
     +              +s_tz_hatd(i,j,k)*s_hatd(i,j,k)*tz_hatd(i,j,k))
               
               e4=e4+(L(k)**4.)*(tfr**8.)*
     +              ((s_hatd(i,j,k)*tx_hatd(i,j,k))**2.
     +              +(s_hatd(i,j,k)*ty_hatd(i,j,k))**2.
     +              +(s_hatd(i,j,k)*tz_hatd(i,j,k))**2.)
               

            end do
         end do

         aa= a2*c4-a4*c2
         bb=-a4*d2+b2*c4
         cc=-c2*b4+a2*d4-a4*e2
         dd= b2*d4-b4*d2
         ee= a2*e4-b4*e2
         ff= b2*e4

         if (t.le.1000) then
            call root8(rtnewt,aa,bb,cc,dd,ee,ff)
         else
            call newroots(rtnewt,aa,bb,cc,dd,ee,ff)
         end if
         
         betaa1D(k)=rtnewt
 
         if(betaa1D(k).le.0.or.betaa1D(k).gt.1.2) then
            betaa1D(k)=1.0
         end if
         if(model.eq.2)then
            betaa1D(k)=1.
         end if
      end do

      do k=2,nzb+1
         bad_pr1(k)=0.0
         do j=1,ny
            do i=1,nx

            T_up = 0.
            T_down = 0.
               
            do mX=-1,1
               do mY=-1,1

                  if(i+mX.gt.Nx)then
                        ii=i+mX-Nx
                     elseif(i+mX.lt.1)then
                        ii=i+mX+Nx
                     else
                        ii=i+mX
                     endif

                     if(j+mY.gt.Ny)then
                        jj=j+mY-Ny
                     elseif(j+mY.lt.1)then
                        jj=j+mY+Ny
                     else
                        jj=j+mY
                     endif

               a2=L(k)**2.*
     +              ((ut_hat(ii,jj,k)-
     +              u_hat(ii,jj,k)*t_hat(ii,jj,k))*s_tx_hat(ii,jj,k)+
     +              (vt_hat(ii,jj,k)-
     +              v_hat(ii,jj,k)*t_hat(ii,jj,k))*s_ty_hat(ii,jj,k)+
     +              (wt_hat(ii,jj,k)-
     +              w_hat(ii,jj,k)*t_hat(ii,jj,k))*s_tz_hat(ii,jj,k))
               
               b2=-(L(k)**2.)*tfr**2.*
     +              ((ut_hat(ii,jj,k)-u_hat(ii,jj,k)*t_hat(ii,jj,k))*
     +              s_hat(ii,jj,k)*tx_hat(ii,jj,k)
     +                 +(vt_hat(ii,jj,k)-v_hat(ii,jj,k)*t_hat(ii,jj,k))*
     +              s_hat(ii,jj,k)*ty_hat(ii,jj,k)
     +              +(wt_hat(ii,jj,k)-w_hat(ii,jj,k)*t_hat(ii,jj,k))*
     +              s_hat(ii,jj,k)*tz_hat(ii,jj,k))
               
               c2=L(k)**4.*
     +              ((s_tx_hat(ii,jj,k))**2.
     +              +(s_ty_hat(ii,jj,k))**2.
     +              +(s_tz_hat(ii,jj,k))**2.)

               d2=-(2.*L(k)**4.)*(tfr**2.)*
     +              (s_tx_hat(ii,jj,k)*s_hat(ii,jj,k)*tx_hat(ii,jj,k)
     +              +s_ty_hat(ii,jj,k)*s_hat(ii,jj,k)*ty_hat(ii,jj,k)
     +              +s_tz_hat(ii,jj,k)*s_hat(ii,jj,k)*tz_hat(ii,jj,k))
               
               e2=(L(k)**4.)*(tfr**4.)*
     +              ((s_hat(ii,jj,k)*tx_hat(ii,jj,k))**2.
     +              +(s_hat(ii,jj,k)*ty_hat(ii,jj,k))**2.
     +              +(s_hat(ii,jj,k)*tz_hat(ii,jj,k))**2.)
               
               
               T_UP=T_up+a2+betaa1D(k)*b2
               T_DOWN=T_down+c2+betaa1D(k)*d2+e2*betaa1D(k)**2
               
            end do
            end do

            if ((T_UP/T_DOWN).lt.0) then
               bad_pr1(k) = bad_pr1(k)+1.0
            end if

           if(abs(T_DOWN).lt.(1e-10).or.(T_UP/T_DOWN).lt.0
     +         .or.(T_UP/T_DOWN).gt.1.0) then

             Pr2(i,j,k) = 0.0
           else
             Pr2(i,j,k) = T_UP/T_DOWN
           end if

        end do
      end do
      
      betaa1(:,:,k) = betaa1D(k)
      
      end do
      
      return
      end












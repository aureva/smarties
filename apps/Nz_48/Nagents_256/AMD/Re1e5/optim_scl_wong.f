      Subroutine optim_scl_wong(Pr2,S11,S33,S22,S12,S13,S23,S_hat,
     +     S_hatd,u_,v_,w_,L,t,t_,tx,ty,tz,me,betaa1)
      
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

      real*8 a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,betaa1D(nz2)

      real*8 aa,bb,cc,dd,ee,ff,L1,L2,L3,M1,M2,M3,Q1,Q2,Q3,
     +     LM,MM,rtnewt,l(nz2)
     
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

         a1=0.
         a2=0.
         a3=0.
         a4=0.
         a5=0.
         a6=0.
         a7=0.
         a8=0.
         a9=0.
         a10=0.

         do j=1,ny
            do i=1,nx

               L1=(ut_hat(i,j,k))-(u_hat(i,j,k))*(t_hat(i,j,k))
               L2=(vt_hat(i,j,k))-(v_hat(i,j,k))*(t_hat(i,j,k))
               L3=(wt_hat(i,j,k))-(w_hat(i,j,k))*(t_hat(i,j,k))
               
               Q1=(ut_hatd(i,j,k))-(u_hatd(i,j,k))*(t_hatd(i,j,k))
               Q2=(vt_hatd(i,j,k))-(v_hatd(i,j,k))*(t_hatd(i,j,k))
               Q3=(wt_hatd(i,j,k))-(w_hatd(i,j,k))*(t_hatd(i,j,k))
       
               a1=a1+(Q1*tx_hatd(i,j,k)+ 
     +              Q2*ty_hatd(i,j,k)+
     +              Q3*tz_hatd(i,j,k))
 
               a2=a1*(-tfr**(8./3.))

               a3=a3+(tx_hat(i,j,k)**2.+
     +              ty_hat(i,j,k)**2.+
     +              tz_hat(i,j,k)**2.)

               a4=(-2.*tfr**(4./3.))*a3
     
               a5=(tfr**(8./3.))*a3

               a6=a6+(L1*tx_hat(i,j,k)+ 
     +              L2*ty_hat(i,j,k)+
     +              L3*tz_hat(i,j,k))

               a7=a6*(-tfr**(4./3.))
               
               a8=a8+(tx_hatd(i,j,k)**2.+
     +              ty_hatd(i,j,k)**2.+
     +              tz_hatd(i,j,k)**2.)

               a9=(-2.*tfr**(8./3.))*a8

               a10=(tfr**(16./3.))*a8

            end do
         end do

         aa = a1*a3-a6*a8
         bb = a1*a4-a7*a8         
         cc = a2*a3+a1*a5-a6*a9
         dd = a2*a4-a7*a9
         ee = a2*a5-a6*a10
         ff = -a7*a10

         if (t.le.1000) then
            call root8(rtnewt,aa,bb,cc,dd,ee,ff)
         else 
            call newroots(rtnewt,aa,bb,cc,dd,ee,ff)
         end if
         
         betaa1D(k)=rtnewt
 
         if(betaa1D(k).le.0.or.betaa1D(k).gt.5.0) then
            betaa1D(k)=1.0
         end if

         if(model.eq.2)then
            betaa1D(k)=1.
         end if
      end do

      do k=2,nzb+1
         do j=1,ny
            do i=1,nx

            LM = 0.
            MM = 0.
               
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
                  
                  L1=(ut_hat(ii,jj,k))-(u_hat(ii,jj,k))*(t_hat(ii,jj,k))
                  L2=(vt_hat(ii,jj,k))-(v_hat(ii,jj,k))*(t_hat(ii,jj,k))
                  L3=(wt_hat(ii,jj,k))-(w_hat(ii,jj,k))*(t_hat(ii,jj,k))
                  
                  M1=L(k)**(4./3.)*(1.-tfr**(4./3.)*betaa1D(k))*
     +                 tx_hat(ii,jj,k)

                  M2=L(k)**(4./3.)*(1.-tfr**(4./3.)*betaa1D(k))*
     +                 ty_hat(ii,jj,k)

                  M3=L(k)**(4./3.)*(1.-tfr**(4./3.)*betaa1D(k))*
     +                 tz_hat(ii,jj,k)
               
                  LM=LM + L1*M1+L2*M2+L3*M3
                  MM=MM + M1*M1+M2*M2+M3*M3
               
               end do
            end do

      if(abs(MM).lt.(1e-10).or.(LM/MM).lt.0.or.(LM/MM).gt.10) then
        Pr2(i,j,k) = 0.0
      else
        Pr2(i,j,k) = LM/MM
      end if
            
         end do
      end do

      betaa1(:,:,k) = betaa1D(k)
      
      end do
      
      return
      end












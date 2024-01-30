      Subroutine OPTIM_LAG (Cs2,S11,S33,S22,S12,S13,S23,S,S_hat,S_hatd,
     +     u_,v_,w_,L,t,me,betaa,a1_old,b1_old,c1_old,d1_old,e1_old,
     +     a2_old,b2_old,c2_old,d2_old,e2_old,ncs2g,ncs2l)

      implicit none
      include 'dimen.h'
      integer*4 i,j,k,t  
      real*8,dimension(nx,ny,nz2)::
     +     S,S11,S22,S33,S12,S13,S23,S_hat,S11_hat,S22_hat,S33_hat,
     +     S12_hat,S13_hat,S23_hat,u_hat,v_hat,w_hat,uu_hat,uv_hat,
     +     uw_hat,vv_hat,vw_hat,ww_hat,S_hatd,S11_hatd,S22_hatd,
     +     S33_hatd,S12_hatd,S13_hatd,S23_hatd,u_hatd,v_hatd,w_hatd,
     +     uu_hatd,uv_hatd,uw_hatd,vv_hatd,vw_hatd,ww_hatd,
     +     SS11_hat,SS12_hat,SS13_hat,SS22_hat,SS23_hat,SS33_hat,
     +     SS11_hatd,SS12_hatd,SS13_hatd,SS22_hatd,SS23_hatd,SS33_hatd,
     +     uu,uv,uw,vv,vw,ww,SS11,SS22,SS33,SS12,SS13,SS23,w_,u_,v_

      real*8,dimension(nx,ny,nz2)::
     +     a1,b1,c1,d1,e1,a2,b2,c2,d2,e2,b,betaa,
     +     a1_old,b1_old,c1_old,d1_old,e1_old,a2_old,b2_old,c2_old,
     +     d2_old,e2_old,Cs2

      real*8 l(Nz2),M11,M12,M13,M22,M23,M33,
     +     LM,MM,L13,L23,L12,L22,L33,L11,Q13,Q23,Q12,Q22,Q33,Q11,
     +     ll,QN,NN,trace,Cs_ave,b_ave,uin,vin,win

      real*8,dimension(Nz2):: 
     +     a1_avg,b1_avg,c1_avg,d1_avg,e1_avg,a2_avg,b2_avg,c2_avg,
     +     d2_avg,e2_avg,Z

      real*8 X(nx+1),Y(ny+1),rtnewt,ncs2g,ncs2l 

      ncs2g=0.
      ncs2l=0.
	  
      do k=2,Nzb+1
         do j=1,Ny
            do i=1,Nx
 
               S(i,j,k)=sqrt(2.*(S11(i,j,k)**2.+S22(i,j,k)**2.+
     +              S33(i,j,k)**2.+2.*S12(i,j,k)**2.+
     +              2.*S13(i,j,k)**2.+
     +              2.*S23(i,j,k)**2.))

               SS11(i,j,k)=S(i,j,k)*S11(i,j,k)
               SS33(i,j,k)=S(i,j,k)*S33(i,j,k)
               SS22(i,j,k)=S(i,j,k)*S22(i,j,k)
               SS12(i,j,k)=S(i,j,k)*S12(i,j,k)
               SS13(i,j,k)=S(i,j,k)*S13(i,j,k)
               SS23(i,j,k)=S(i,j,k)*S23(i,j,k)
               
               uu(i,j,k)=u_(i,j,k)**2.
               vv(i,j,k)=v_(i,j,k)**2.
               ww(i,j,k)=w_(i,j,k)**2.
               uv(i,j,k)=u_(i,j,k)*v_(i,j,k)
               vw(i,j,k)=v_(i,j,k)*w_(i,j,k)
               uw(i,j,k)=u_(i,j,k)*w_(i,j,k)

            end do
         end do
      end do
      
CC...Filtering to get the _hat variables (coarser resolution)

      Call Filter_2Laa(u_hat,u_hatd,u_,t,1)
      Call Filter_2Laa(v_hat,v_hatd,v_,t,0)
      Call Filter_2Laa(w_hat,w_hatd,w_,t,0)
      Call Filter_2Laa(uu_hat,uu_hatd,uu,t,0)
      Call Filter_2Laa(vv_hat,vv_hatd,vv,t,0)
      Call Filter_2Laa(ww_hat,ww_hatd,ww,t,0)
      Call Filter_2Laa(uv_hat,uv_hatd,uv,t,0)
      Call Filter_2Laa(uw_hat,uw_hatd,uw,t,0)
      Call Filter_2Laa(vw_hat,vw_hatd,vw,t,0)
      Call Filter_2Laa(S11_hat,S11_hatd,S11,t,0)
      Call Filter_2Laa(S22_hat,S22_hatd,S22,t,0)
      Call Filter_2Laa(S33_hat,S33_hatd,S33,t,0)
      Call Filter_2Laa(S12_hat,S12_hatd,S12,t,0)
      Call Filter_2Laa(S13_hat,S13_hatd,S13,t,0)
      Call Filter_2Laa(S23_hat,S23_hatd,S23,t,0)
      Call Filter_2Laa(SS11_hat,SS11_hatd,SS11,t,0)
      Call Filter_2Laa(SS22_hat,SS22_hatd,SS22,t,0)
      Call Filter_2Laa(SS33_hat,SS33_hatd,SS33,t,0)
      Call Filter_2Laa(SS12_hat,SS12_hatd,SS12,t,0)
      Call Filter_2Laa(SS13_hat,SS13_hatd,SS13,t,0)
      if(S_flag.eq.1)then
         Call Filter_2Laa(SS23_hat,SS23_hatd,SS23,t,0)
      else
         Call Filter_2Laa(SS23_hat,SS23_hatd,SS23,t,2)
      endif
      
      do k=2,nzb+1
         do j=1,ny
            do i=1,nx
               S_hat(i,j,k)=sqrt(2.*(S11_hat(i,j,k)**2.+
     +              S22_hat(i,j,k)**2.+S33_hat(i,j,k)**2.+
     +              2.*S12_hat(i,j,k)**2.+2.*S13_hat(i,j,k)**2.+
     +              2.*S23_hat(i,j,k)**2.))
               S_hatd(i,j,k)=sqrt(2.*(S11_hatd(i,j,k)**2.+
     +              S22_hatd(i,j,k)**2.+S33_hatd(i,j,k)**2.+
     +              2.*S12_hatd(i,j,k)**2.+2.*S13_hatd(i,j,k)**2.+
     +              2.*S23_hatd(i,j,k)**2.))
            end do
         end do
      end do
      
      do k=2,Nzb+1
         do j=1,ny
            do i=1,nx
               L11=(uu_hat(i,j,k))-(u_hat(i,j,k))**2.
               L22=(vv_hat(i,j,k))-(v_hat(i,j,k))**2.
               L12=(uv_hat(i,j,k))-(u_hat(i,j,k)*v_hat(i,j,k))
               L13=(uw_hat(i,j,k))-(u_hat(i,j,k)*w_hat(i,j,k))
               L23=(vw_hat(i,j,k))-(v_hat(i,j,k)*w_hat(i,j,k))
               L33=(ww_hat(i,j,k))-(w_hat(i,j,k))**2.
               Q11=(uu_hatd(i,j,k))-(u_hatd(i,j,k))**2.
               Q22=(vv_hatd(i,j,k))-(v_hatd(i,j,k))**2.
               Q12=(uv_hatd(i,j,k))-(u_hatd(i,j,k)*v_hatd(i,j,k))
               Q13=(uw_hatd(i,j,k))-(u_hatd(i,j,k)*w_hatd(i,j,k))
               Q23=(vw_hatd(i,j,k))-(v_hatd(i,j,k)*w_hatd(i,j,k))
               Q33=(ww_hatd(i,j,k))-(w_hatd(i,j,k))**2.

               a1(i,j,k)=-2.*l(k)**2.*tfr**2.*S_hat(i,j,k)*
     +              (L11*S11_hat(i,j,k)+L22*S22_hat(i,j,k)+
     +              L33*S33_hat(i,j,k)+2.*(L12*S12_hat(i,j,k)+L13*
     +              S13_hat(i,j,k)+L23*S23_hat(i,j,k)))
               
               a2(i,j,k)=-2.*l(k)**2.*tfr**4.*S_hatd(i,j,k)*
     +              (Q11*S11_hatd(i,j,k)+Q22*S22_hatd(i,j,k)+Q33*
     +              S33_hatd(i,j,k)+2.*(Q12*S12_hatd(i,j,k)
     +              +Q13*S13_hatd(i,j,k)+Q23*S23_hatd(i,j,k)))
               
               b1(i,j,k)=-2.*l(k)**2.*(L11*SS11_hat(i,j,k)+L22*
     +              SS22_hat(i,j,k)+L33*SS33_hat(i,j,k)+
     +              2.*(L12*SS12_hat(i,j,k)+L13*SS13_hat(i,j,k)
     +              +L23*SS23_hat(i,j,k)))
               
               b2(i,j,k)=-2.*l(k)**2.*(Q11*SS11_hatd(i,j,k)+
     +              Q22*SS22_hatd(i,j,k)+Q33*SS33_hatd(i,j,k)+2.*
     +              (Q12*SS12_hatd(i,j,k)+Q13*SS13_hatd(i,j,k)+
     +              Q23*SS23_hatd(i,j,k)))
               
               c1(i,j,k)=(2.*l(k)**2.)**2.*
     +              (SS11_hat(i,j,k)**2.+SS22_hat(i,j,k)**2.
     +              +SS33_hat(i,j,k)**2.+2.*(SS12_hat(i,j,k)**2.+
     +              SS13_hat(i,j,k)**2.+SS23_hat(i,j,k)**2.))
               
               c2(i,j,k)=(2.*l(k)**2.)**2.*(SS11_hatd(i,j,k)**2.+
     +              SS22_hatd(i,j,k)**2.+SS33_hatd(i,j,k)**2.
     +              +2.*(SS12_hatd(i,j,k)**2.+SS13_hatd(i,j,k)**2.+
     +              SS23_hatd(i,j,k)**2.))
               
               d1(i,j,k)=(4.*l(k)**4.)*tfr**4.*S_hat(i,j,k)**2.*
     +              (S11_hat(i,j,k)**2.+S22_hat(i,j,k)**2.
     +              +S33_hat(i,j,k)**2.+2.*(S12_hat(i,j,k)**2.+
     +              S13_hat(i,j,k)**2.+S23_hat(i,j,k)**2.)) 
               
               d2(i,j,k)=(4.*l(k)**4.)*tfr**8.*S_hatd(i,j,k)**2.
     +              *(S11_hatd(i,j,k)**2.+S22_hatd(i,j,k)**2.+
     +              S33_hatd(i,j,k)**2.+2.*(S12_hatd(i,j,k)**2.+
     +              S13_hatd(i,j,k)**2.+S23_hatd(i,j,k)**2.)) 
               
               e1(i,j,k)=(8.*l(k)**4.)*tfr**2.*
     +              S_hat(i,j,k)*(S11_hat(i,j,k)*
     +              SS11_hat(i,j,k)+S22_hat(i,j,k)*
     +              SS22_hat(i,j,k)+S33_hat(i,j,k)*
     +              SS33_hat(i,j,k)+2.*(S12_hat(i,j,k)*
     +              SS12_hat(i,j,k)+S13_hat(i,j,k)*
     +              SS13_hat(i,j,k)+S23_hat(i,j,k)*
     +              SS23_hat(i,j,k)))             
               
               e2(i,j,k)=(8.*l(k)**4.)*tfr**4.*S_hatd(i,j,k)*
     +              (S11_hatd(i,j,k)*SS11_hatd(i,j,k)+
     +              S22_hatd(i,j,k)*SS22_hatd(i,j,k)+
     +              S33_hatd(i,j,k)*SS33_hatd(i,j,k)+2.
     +              *(S12_hatd(i,j,k)*SS12_hatd(i,j,k)+
     +              S13_hatd(i,j,k)*SS13_hatd(i,j,k)+
     +              S23_hatd(i,j,k)*SS23_hatd(i,j,k))) 
               
            enddo
         enddo
      enddo

c     **************initilization for lagrangian averaging***********
      if(t.le.50.and.INITU.ne.1)then
         do k=2,nzb+1
            a1_avg(k)=0.
            a2_avg(k)=0.
            b1_avg(k)=0.
            b2_avg(k)=0.
            c1_avg(k)=0.
            c2_avg(k)=0.
            d1_avg(k)=0.
            d2_avg(k)=0.
            e1_avg(k)=0.
            e2_avg(k)=0.	 
            do j=1,Ny
               do i=1,Nx
                  a1_avg(k)=a1_avg(k)+a1(i,j,k)/(Nx*Ny)
                  a2_avg(k)=a2_avg(k)+a2(i,j,k)/(Nx*Ny)
                  b1_avg(k)=b1_avg(k)+b1(i,j,k)/(Nx*Ny)
                  b2_avg(k)=b2_avg(k)+b2(i,j,k)/(Nx*Ny)
                  c1_avg(k)=c1_avg(k)+c1(i,j,k)/(Nx*Ny)
                  c2_avg(k)=c2_avg(k)+c2(i,j,k)/(Nx*Ny)
                  d1_avg(k)=d1_avg(k)+d1(i,j,k)/(Nx*Ny)
                  d2_avg(k)=d2_avg(k)+d2(i,j,k)/(Nx*Ny)
                  e1_avg(k)=e1_avg(k)+e1(i,j,k)/(Nx*Ny)
                  e2_avg(k)=e2_avg(k)+e2(i,j,k)/(Nx*Ny)
               enddo
            enddo
         enddo
         
         do k=2,nzb+1
            do j=1,Ny
               do i=1,Nx
                  
                  a1(i,j,k)=a1_avg(k)
                  b1(i,j,k)=b1_avg(k)
                  c1(i,j,k)=c1_avg(k)
                  d1(i,j,k)=d1_avg(k)
                  e1(i,j,k)=e1_avg(k)
                  a2(i,j,k)=a2_avg(k)
                  b2(i,j,k)=b2_avg(k)
                  c2(i,j,k)=c2_avg(k)
                  d2(i,j,k)=d2_avg(k)
                  e2(i,j,k)=e2_avg(k)

                  a1_old(i,j,k)=a1(i,j,k)
                  a2_old(i,j,k)=a2(i,j,k)
                  b1_old(i,j,k)=b1(i,j,k)
                  b2_old(i,j,k)=b2(i,j,k)
                  c1_old(i,j,k)=c1(i,j,k)
                  d1_old(i,j,k)=d1(i,j,k)
                  e1_old(i,j,k)=e1(i,j,k)
                  c2_old(i,j,k)=c2(i,j,k)
                  d2_old(i,j,k)=d2(i,j,k)
                  e2_old(i,j,k)=e2(i,j,k)
                  
                  if(t.le.5)then
                     betaa(i,j,k)=1.
                  endif

               enddo
            enddo
         enddo

      endif
c     **************end of the lagrangian initilization**************
      
ccccccccccccccccccccccccccccccccccccc
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

               call lagrng_sd(a1(i,j,k),b1(i,j,k),c1(i,j,k),
     +              d1(i,j,k),e1(i,j,k),a1_old,b1_old,c1_old,
     +              d1_old,e1_old,betaa(i,j,k),uin,vin,win,
     +              i,j,k,0,me,t)
               
c               call lagrng_sd(a2(i,j,k),b2(i,j,k),c2(i,j,k),
c     +              d2(i,j,k),e2(i,j,k),a2_old,b2_old,c2_old,
c     +              d2_old,e2_old,betaa(i,j,k),uin,vin,win,
c     +              i,j,k,1,me,t)
c
c               call root(rtnewt,a1(i,j,k),a2(i,j,k),b1(i,j,k),
c     +              b2(i,j,k),c1(i,j,k),c2(i,j,k),d1(i,j,k),
c     +              d2(i,j,k),e1(i,j,k),e2(i,j,k))
c               
c               if(rtnewt.le.0.or.rtnewt.gt.(100.0))then
c                                    
c                  betaa(i,j,k)=0.0
c                  cs2(i,j,k)=0.0
c                  betaa(i,j,k)=1.0
c                  cs2(i,j,k)=0.001				  
c                  
c                  a1(i,j,k)=0.
c                  b1(i,j,k)=0.
c                  c1(i,j,k)=0.
c                  d1(i,j,k)=0.
c                  e1(i,j,k)=0.
c                  
c                  a2(i,j,k)=0.
c                  b2(i,j,k)=0.
c                  c2(i,j,k)=0.
c                  d2(i,j,k)=0.
c                  e2(i,j,k)=0.
c                  
c               else
c                  
c                  if(t.le.1.and.INITU.ne.1) then
c                     betaa(i,j,k)=1.
c                  else	         
c                     betaa(i,j,k)=rtnewt
c                  end if
                  
                  betaa(i,j,k)=1-0.65*exp(-0.7*(me*nzb+k-1.5)*dz/delta)				  
                  LM=(a1(i,j,k)*betaa(i,j,k)-b1(i,j,k))
                  MM=(c1(i,j,k)+d1(i,j,k)*betaa(i,j,k)**2.-
     +                 e1(i,j,k)*betaa(i,j,k))
c                  QN=(a2(i,j,k)*betaa(i,j,k)**2.-b2(i,j,k))
c                  NN=(c2(i,j,k)+d2(i,j,k)*betaa(i,j,k)**4.-
c     +                 e2(i,j,k)*betaa(i,j,k)**2.)
c
                  if(LM.lt.0.or.MM.lt.(1.d-20))then 
                  
                     cs2(i,j,k)=0.0
                     
                     a1(i,j,k)=0.
                     b1(i,j,k)=0.
c                     c1(i,j,k)=0.
c                     d1(i,j,k)=0.
c                     e1(i,j,k)=0.
                     
c                     a2(i,j,k)=0.
c                     b2(i,j,k)=0.
c                     c2(i,j,k)=0.
c                     d2(i,j,k)=0.
c                     e2(i,j,k)=0.
                     
c                  elseif(QN.lt.(0.0).or.NN.lt.(1.d-20))then    
c                     
c                     cs2(i,j,k)=0.0
c                     
c                     a1(i,j,k)=0.
c                     b1(i,j,k)=0.
c                    c1(i,j,k)=0.
c                    d1(i,j,k)=0.
c                    e1(i,j,k)=0.
c                     
c                     a2(i,j,k)=0.
c                     b2(i,j,k)=0.
c                    c2(i,j,k)=0.
c                    d2(i,j,k)=0.
c                    e2(i,j,k)=0.
c                     
                  else
                     
                     cs2(i,j,k)=LM/MM

                  endif      
c               endif
            if (cs2(i,j,k) .gt. 0.09) then
                  cs2(i,j,k)=0.09
				  ncs2g=ncs2g+1
            endif
            if (cs2(i,j,k) .lt. 0.0001) then
                  cs2(i,j,k)=0.0001	
                  ncs2l=ncs2l+1				  
            endif				  
            enddo
         enddo
      enddo

      do k=2,nzb+1
         do j=1,ny
            do i=1,nx
               a1_old(i,j,k)=a1(i,j,k)
               b1_old(i,j,k)=b1(i,j,k)
               c1_old(i,j,k)=c1(i,j,k)
               d1_old(i,j,k)=d1(i,j,k)
               e1_old(i,j,k)=e1(i,j,k)
               a2_old(i,j,k)=a2(i,j,k)
               b2_old(i,j,k)=b2(i,j,k)
               c2_old(i,j,k)=c2(i,j,k)
               d2_old(i,j,k)=d2(i,j,k)
               e2_old(i,j,k)=e2(i,j,k)
            enddo
         enddo
      enddo
      
      return
      
      end
      















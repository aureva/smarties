      Subroutine OPTIM_LAG_DYN (Cs2,S11,S33,S22,S12,S13,S23,S,S_hat,
     +     u_,v_,w_,L,t,me,LM_old,MM_old)

      implicit none
      include 'dimen.h'
      integer*4 i,j,k,t   
      real*8,dimension(nx,ny,nz2)::
     +     S,S11,S22,S33,S12,S13,S23,S_hat,S11_hat,S22_hat,S33_hat,
     +     S12_hat,S13_hat,S23_hat,u_hat,v_hat,w_hat,uu_hat,uv_hat,
     +     uw_hat,vv_hat,vw_hat,ww_hat,
     +     SS11_hat,SS12_hat,SS13_hat,SS22_hat,SS23_hat,SS33_hat,
     +     uu,uv,uw,vv,vw,ww,SS11,SS22,SS33,SS12,SS13,SS23,w_,u_,v_

      real*8,dimension(nx,ny,nz2):: LM,MM,
     +     LM_old,MM_old,Cs2

      real*8 l(Nz2),M11,M12,M13,M22,M23,M33,
     +     L13,L23,L12,L22,L33,L11,
     +     ll,trace,Cs_ave,uin,vin,win

      real*8,dimension(Nz2):: 
     +     LM_avg,MM_avg,Z

      real*8 X(nx+1),Y(ny+1),a1,b1,c1,d1,e1

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

      Call Filter_La(u_hat,u_,t,1)
      Call Filter_La(v_hat,v_,t,0)
      Call Filter_La(w_hat,w_,t,0)
      Call Filter_La(uu_hat,uu,t,0)
      Call Filter_La(vv_hat,vv,t,0)
      Call Filter_La(ww_hat,ww,t,0)
      Call Filter_La(uv_hat,uv,t,0)
      Call Filter_La(uw_hat,uw,t,0)
      Call Filter_La(vw_hat,vw,t,0)
      Call Filter_La(S11_hat,S11,t,0)
      Call Filter_La(S22_hat,S22,t,0)
      Call Filter_La(S33_hat,S33,t,0)
      Call Filter_La(S12_hat,S12,t,0)
      Call Filter_La(S13_hat,S13,t,0)
      Call Filter_La(S23_hat,S23,t,0)
      Call Filter_La(SS11_hat,SS11,t,0)
      Call Filter_La(SS22_hat,SS22,t,0)
      Call Filter_La(SS33_hat,SS33,t,0)
      Call Filter_La(SS12_hat,SS12,t,0)
      Call Filter_La(SS13_hat,SS13,t,0)
      if(S_flag.eq.1)then
         Call Filter_La(SS23_hat,SS23,t,0)
      else
         Call Filter_La(SS23_hat,SS23,t,2)
      endif
      
      do k=2,nzb+1
         do j=1,ny
            do i=1,nx
               S_hat(i,j,k)=sqrt(2.*(S11_hat(i,j,k)**2.+
     +              S22_hat(i,j,k)**2.+S33_hat(i,j,k)**2.+
     +              2.*S12_hat(i,j,k)**2.+2.*S13_hat(i,j,k)**2.+
     +              2.*S23_hat(i,j,k)**2.))
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

               a1=-2.*l(k)**2.*tfr**2.*S_hat(i,j,k)*
     +              (L11*S11_hat(i,j,k)+L22*S22_hat(i,j,k)+
     +              L33*S33_hat(i,j,k)+2.*(L12*S12_hat(i,j,k)+L13*
     +              S13_hat(i,j,k)+L23*S23_hat(i,j,k)))
               
               b1=-2.*l(k)**2.*(L11*SS11_hat(i,j,k)+L22*
     +              SS22_hat(i,j,k)+L33*SS33_hat(i,j,k)+
     +              2.*(L12*SS12_hat(i,j,k)+L13*SS13_hat(i,j,k)
     +              +L23*SS23_hat(i,j,k)))
               
               c1=(2.*l(k)**2.)**2.*
     +              (SS11_hat(i,j,k)**2.+SS22_hat(i,j,k)**2.
     +              +SS33_hat(i,j,k)**2.+2.*(SS12_hat(i,j,k)**2.+
     +              SS13_hat(i,j,k)**2.+SS23_hat(i,j,k)**2.))
               
               d1=(4.*l(k)**4.)*tfr**4.*S_hat(i,j,k)**2.*
     +              (S11_hat(i,j,k)**2.+S22_hat(i,j,k)**2.
     +              +S33_hat(i,j,k)**2.+2.*(S12_hat(i,j,k)**2.+
     +              S13_hat(i,j,k)**2.+S23_hat(i,j,k)**2.)) 
               
               e1=(8.*l(k)**4.)*tfr**2.*
     +              S_hat(i,j,k)*(S11_hat(i,j,k)*
     +              SS11_hat(i,j,k)+S22_hat(i,j,k)*
     +              SS22_hat(i,j,k)+S33_hat(i,j,k)*
     +              SS33_hat(i,j,k)+2.*(S12_hat(i,j,k)*
     +              SS12_hat(i,j,k)+S13_hat(i,j,k)*
     +              SS13_hat(i,j,k)+S23_hat(i,j,k)*
     +              SS23_hat(i,j,k)))             

               LM(i,j,k) = a1 - b1
               MM(i,j,k) = c1 + d1 - e1

            enddo
         enddo
      enddo

c     **************initilization for lagrangian averaging***********
      if(t.le.50.and.INITU.ne.1)then
         do k=2,nzb+1
            LM_avg(k)=0.
            MM_avg(k)=0.
            	 
            do j=1,Ny
               do i=1,Nx

                  LM_avg(k)=LM_avg(k)+LM(i,j,k)/(Nx*Ny)
                  MM_avg(k)=MM_avg(k)+MM(i,j,k)/(Nx*Ny)
                  
               enddo
            enddo
         enddo
         
         do k=2,nzb+1
            do j=1,Ny
               do i=1,Nx
                  
                  LM(i,j,k)=LM_avg(k)
                  MM(i,j,k)=MM_avg(k)

                  LM_old(i,j,k)=LM(i,j,k)
                  MM_old(i,j,k)=MM(i,j,k)

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

               call lagrng_dyn(LM(i,j,k),MM(i,j,k),LM_old,MM_old,
     +              uin,vin,win,i,j,k,me,t)

               if(LM(i,j,k).lt.0.or.MM(i,j,k).lt.(1e-10))then 
                  
                  cs2(i,j,k)=0.0
                  
                  LM(i,j,k)=0.
c                 MM(i,j,k)=0.
                     
               else
                  
                  cs2(i,j,k)=LM(i,j,k)/MM(i,j,k)
                  
               endif      
                  
            enddo
         enddo
      enddo

      do k=2,nzb+1
         do j=1,ny
            do i=1,nx

               LM_old(i,j,k)=LM(i,j,k)
               MM_old(i,j,k)=MM(i,j,k)
               
            enddo
         enddo
      enddo
      
      return
      
      end
      















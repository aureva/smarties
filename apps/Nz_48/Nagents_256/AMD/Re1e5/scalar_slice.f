      subroutine scalar_slice(u,v,w,theta,sgs_t1,sgs_t2,sgs_t3,dtdx,
     +     dtdy,dtdz,Pr2,Cs2,beta2,at,t2,t3,asgs_t1,asgs_t2,asgs_t3,
     +     aut,avt,awt,adtdx,adtdy,adtdz,aPr,aCs2Pr,abeta2,aqz_s,
     +     ET3D,aET,ilow,wgx,t_s,ts_avg,me)

      implicit none
      include 'dimen.h'
      integer*4 i,j,k,l,kk,ilow(nx),h
      real*8,dimension(nx,ny,nz2):: u,v,w,theta,sgs_t1,sgs_t2,sgs_t3,
     +     dtdx,dtdy,dtdz,Pr2,Cs2,beta2,ET3D
      real*8,dimension(anx,nz2) :: at,t2,t3,asgs_t1,asgs_t2,asgs_t3,
     +     aut,avt,awt,adtdx,adtdy,adtdz,aPr,aCs2Pr,abeta2,aET
      real*8,dimension(nx,ny) :: aqz_s,t_s,ts_avg 
      
      real*8,dimension(anx,nz2) :: tu1,tt1,tt2,tt3,tsgst1,tsgst2,tsgst3,
     +     tut,tvt,twt,tdtdx,tdtdy,tCs2,tdtdz,tCs2Pr,tbeta2,tET,tw1
      real*8,dimension(nz2) :: u_bar,v_bar,w_bar,t_bar
      real*8 arg1,arg2,fr,norm,ftn,wgx

      fr=(1./p_count)*c_count
      norm=1.d0/(Ny*(1+Nx-aNx))
      ftn=fr*norm

ccc compute the plane averages of terms involved in products ccc
      u_bar=0.d0
      w_bar=0.d0
      t_bar=0.d0
      do k=1,Nzb+1
         do j=1,Ny
            do i=1,Nx
               t_bar(k)=t_bar(k)+theta(i,j,k)
            enddo
         enddo
      enddo
      do k=2,Nzb+1
         do j=1,Ny
            do i=1,Nx
               u_bar(k)=u_bar(k)+u(i,j,k)
               v_bar(k)=v_bar(k)+v(i,j,k)
               w_bar(k)=w_bar(k)+w(i,j,k)
            enddo
         enddo
      enddo
      u_bar=u_bar*inxny
      v_bar=v_bar*inxny
      t_bar=t_bar*inxny
      w_bar=w_bar*inxny
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      tu1=0.
      tt1=0.
      tt2=0.
      tt3=0.
      tsgst1=0.
      tsgst2=0.
      tsgst3=0.
      tut=0.
      tvt=0.
      twt=0.
      tdtdx=0.
      tdtdy=0.
      tdtdz=0.
      tCs2Pr=0.
      tCs2=0.
      tbeta2=0.
      tET=0.
      arg2=0.
      
      do k=2,Nzb+1
         do i=1,aNx

            do j=1,Ny  	
               do l=1,1+Nx-aNx

                  tt1(i,k)=tt1(i,k)+theta(i-1+l,j,k)                  
                  tu1(i,k)=tu1(i,k)+u(i-1+l,j,k)+Ugal
                  tw1(i,k)=tw1(i,k)+w(i-1+l,j,k)
                  tsgst1(i,k)=tsgst1(i,k)+sgs_t1(i-1+l,j,k)
                  tsgst2(i,k)=tsgst2(i,k)+sgs_t2(i-1+l,j,k)
                  tsgst3(i,k)=tsgst3(i,k)+sgs_t3(i-1+l,j,k)
                  tdtdx(i,k)=tdtdx(i,k)+dtdx(i-1+l,j,k)
                  tdtdy(i,k)=tdtdy(i,k)+dtdy(i-1+l,j,k)
                  tdtdz(i,k)=tdtdz(i,k)+dtdz(i-1+l,j,k)
                  tCs2Pr(i,k)=tCs2Pr(i,k)+Pr2(i-1+l,j,k)
                  tCs2(i,k)=tCs2(i,k)+Cs2(i-1+l,j,k)
                  tbeta2(i,k)=tbeta2(i,k)+beta2(i-1+l,j,k)
                  tET(i,k)=tET(i,k)+ET3D(i-1+l,j,k)
                  
                  tt2(i,k)=tt2(i,k)+(theta(i-1+l,j,k) - t_bar(k))**2
                  tt3(i,k)=tt3(i,k)+(theta(i-1+l,j,k) - t_bar(k))**3
                  tut(i,k)=tut(i,k)+(theta(i-1+l,j,k) - t_bar(k))*
     +                                  (u(i-1+l,j,k) - u_bar(k))
                  tvt(i,k)=tvt(i,k)+(theta(i-1+l,j,k) - t_bar(k))*
     +                                  (v(i-1+l,j,k) - v_bar(k))

                  if (k==2.and.me==0) then
                     arg1=0.
                  else
                     arg1=(theta(i-1+l,j,k)+theta(i-1+l,j,k-1) - 
     +                     t_bar(k) - t_bar(k-1))*0.5d0
                  end if

                  twt(i,k)=twt(i,k)+(w(i-1+l,j,k) - w_bar(k))*arg1

               enddo
            enddo
    
         enddo
      enddo



      do kk=2,Nzb+1
ccc NOTE rightnow just matching with the location of jump (no interpolation)
         do i=1,aNx
            if(aNx.ne.1)then
               l=ilow(i)
               h=ilow(i)
               if(l.eq.nx)then
                  h=nx
               endif
            else
               l=1
               h=1
            endif

            if(nprocs.eq.1)then
               k=kk-1
            else
               k=kk
            endif
			
            l=i
            h=i
	    
            at(i,k)=at(i,k)+ftn*((1-wgx)*tt1(l,k)+wgx*tt1(h,k))
            t2(i,k)=t2(i,k)+ftn*((1-wgx)*tt2(l,k)+wgx*tt2(h,k))
            t3(i,k)=t3(i,k)+ftn*((1-wgx)*tt3(l,k)+wgx*tt3(h,k))
            asgs_t1(i,k)=asgs_t1(i,k)+ftn*((1-wgx)*tsgst1(l,k)+
     +           wgx*tsgst1(h,k))
            asgs_t2(i,k)=asgs_t2(i,k)+ftn*((1-wgx)*tsgst2(l,k)+
     +           wgx*tsgst2(h,k))
            asgs_t3(i,k)=asgs_t3(i,k)+ftn*((1-wgx)*tsgst3(l,k)+
     +           wgx*tsgst3(h,k))
            aut(i,k)=aut(i,k)+ftn*((1-wgx)*tut(l,k)+
     +           wgx*tut(h,k))
            avt(i,k)=avt(i,k)+ftn*((1-wgx)*tvt(l,k)+
     +           wgx*tvt(h,k))
            awt(i,k)=awt(i,k)+ftn*((1-wgx)*twt(l,k)+
     +           wgx*twt(h,k))
            adtdx(i,k)=adtdx(i,k)+ftn*((1-wgx)*tdtdx(l,k)+
     +           wgx*tdtdx(h,k))
            adtdy(i,k)=adtdy(i,k)+ftn*((1-wgx)*tdtdy(l,k)+
     +           wgx*tdtdy(h,k))
            adtdz(i,k)=adtdz(i,k)+ftn*((1-wgx)*tdtdz(l,k)+
     +           wgx*tdtdz(h,k))
            arg2=((1-wgx)*tCs2Pr(l,k)+wgx*tCs2Pr(h,k))
        if(arg2.eq.0)then
               aPr(i,k)=aPr(i,k)
        else
               aPr(i,k)=aPr(i,k)+fr*((1-wgx)*tCs2(l,k)+
     +              wgx*tCs2(h,k))/arg2
            endif
            aCs2Pr(i,k)=aCs2Pr(i,k)+ftn*((1-wgx)*tCs2Pr(l,k)+
     +           wgx*tCs2Pr(h,k))
            abeta2(i,k)=abeta2(i,k)+ftn*((1-wgx)*tbeta2(l,k)+
     +           wgx*tbeta2(h,k))
            aET(i,k)=aET(i,k)+ftn*((1-wgx)*tET(l,k)+
     +           wgx*tET(h,k))

         enddo
      enddo

      if(me.eq.0)then
         do j=1,ny
            do i=1,nx
               aqz_s(i,j)=aqz_s(i,j)+fr*sgs_t3(i,j,2)
               ts_avg(i,j)=ts_avg(i,j)+fr*t_s(i,j)			   
            enddo
         enddo
      endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine scalar_slice1(u,v,w,theta,qx,qy,qz,dtdx,dtdy,dtdz,
     +            Pr2,Cs2,beta2,bat,bat2,baqx,baqy,baqz,baut,bavt,
     +            bawt,badtdx,badtdy,badtdz,baPr2,baBetaH,me)

      implicit none
      include 'dimen.h'
      integer*4 i,j,k
      real*8,dimension(nx,ny,nz2):: u,v,w,theta,qx,qy,qz,dtdx,
     +       dtdy,dtdz,Pr2,Cs2,beta2

      real*8,dimension(nx,ny,nz2):: bat,bat2,baqx,baqy,baqz,
     +       baut,bavt,bawt,badtdx,badtdy,badtdz,baPr2,baBetaH
      real*8 fr,arg1

      fr   =dble(c_count)/dble(p_count)
      
      bat=bat+fr*theta
      bat2=bat2+fr*theta*theta
      baqx=baqx+fr*qx
      baqy=baqy+fr*qy
      baqz=baqz+fr*qz

      baut=baut+fr*u*theta
      bavt=bavt+fr*v*theta

      badtdx=badtdx+fr*dtdx
      badtdy=badtdy+fr*dtdy
      badtdz=badtdz+fr*dtdz

      baPr2  =baPr2   + fr*Pr2
      baBetaH=baBetaH + fr*beta2

      do k=2,Nzb+1
      do j=1,Ny
      do i=1,Nx
      if (k==2.and.me==0) then
         arg1=0.
      else
         arg1=(theta(i,j,k)+theta(i,j,k-1))/2. 
      end if
      bawt(i,j,k)=bawt(i,j,k)+fr*w(i,j,k)*arg1
      end do
      end do
      end do

      return
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      fr   =(1./p_count)*c_count
      
      do k=2,Nzb+1
      do j=1,Ny
      do i=1,Nx
      bat (i,j,k)=bat (i,j,k)+fr*theta(i,j,k)
      bat2(i,j,k)=bat2(i,j,k)+fr*theta(i,j,k)*theta(i,j,k)
      baqx(i,j,k)=baqx(i,j,k)+fr*qx(i,j,k)
      baqy(i,j,k)=baqy(i,j,k)+fr*qy(i,j,k)
      baqz(i,j,k)=baqz(i,j,k)+fr*qz(i,j,k)

      baut(i,j,k)=baut(i,j,k)+fr*u(i,j,k)*theta(i,j,k)
      bavt(i,j,k)=bavt(i,j,k)+fr*v(i,j,k)*theta(i,j,k)

      if (k==2.and.me==0) then
         arg1=0.
      else
         arg1=(theta(i,j,k)+theta(i,j,k-1))/2. 
      end if
      bawt(i,j,k)=bawt(i,j,k)+fr*w(i,j,k)*arg1

      badtdx(i,j,k)=badtdx(i,j,k)+fr*dtdx(i,j,k)
      badtdy(i,j,k)=badtdy(i,j,k)+fr*dtdy(i,j,k)
      badtdz(i,j,k)=badtdz(i,j,k)+fr*dtdz(i,j,k)

      baPr2(i,j,k)=baPr2(i,j,k)+fr*Pr2(i,j,k)
      baBetaH(i,j,k)=baBetaH(i,j,k)+fr*beta2(i,j,k)
      end do
      end do
      end do

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
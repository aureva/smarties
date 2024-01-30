      subroutine avgslice(u,v,w,p,txx,txz,tyy,tyz,tzz,txy,dudz,dudx,
     +     dvdz,dwdz,dwdx,Cs2,beta1,au,av,aw,ap,u2,v2,w2,p2,w3,
     +     atxx,atxz,atyy,atyz,atzz,atxy,auw,avw,auv,adudz,adudx,
     +     advdz,adwdz,adwdx,e,aCs2,aCs,abeta1,atxz_s,ESGS3D,aESGS,
     +     aFSu1,aFSu2,aFSu3,aFSu4,aFSv1,aFSv2,aFSv3,aFSv4,	 
     +     tt,ilow,wgx,me)
      
      implicit none
      include 'dimen.h'
      integer*4 tt,i,j,k,l,kk,ilow(nx),h
      real*8,dimension(nx,ny,nz2):: u,v,w,p,txx,txz,tyy,tyz,tzz,txy,
     +     dudz,dudx,Cs2,beta1,dwdz,dwdx,dvdz,ESGS3D
      real*8,dimension(anx,nz2) :: au,av,aw,ap,u2,v2,w2,p2,w3,atxx,atxz,
     +     atyy,atyz,atzz,atxy,auw,avw,auv,adudz,adudx,e,aCs2,
     +     abeta1,aCs,adwdz,adwdx,advdz,aESGS,
     +     aFSu1,aFSu2,aFSu3,aFSu4,aFSv1,aFSv2,aFSv3,aFSv4 	 
      real*8,dimension(nx,ny) :: atxz_s
      real*8,dimension(anx,nz2) :: tu1,tv1,tw1,tp1,tu2,tv2,tw2,tp2,tw3,
     +     ttxx,ttxz,ttyy,ttyz,ttzz,ttxy,tuw,tvw,tuv,tdudz,tdudx,te,
     +     tCs2,tbeta1,tCs,tdwdz,tdwdx,tdvdz,tESGS,
     +     FSu1,FSu2,FSu3,FSu4,FSv1,FSv2,FSv3,FSv4	 
      real*8,dimension(nz2) :: u_bar,v_bar,w_bar
      real*8 fr,arg1,arg2,norm,DDD,xpnt,ftn,wgx     

      fr   =(1./p_count)*c_count
      norm =1.d0/(Ny*(1+Nx-aNx))
      ftn=fr*norm

      DDD = dmod(dble(tt)*dt*Ugal,dx)
      wgx=(DDD)/dx
c      if(me.eq.0) write(*,*) 'counting, wgx=',wgx

ccc compute the plane averages of terms involved in products ccc
      u_bar=0.d0
      v_bar=0.d0
      w_bar=0.d0
      do k=1,Nzb+1
         do j=1,Ny
            do i=1,Nx
               u_bar(k)=u_bar(k)+u(i,j,k)
               v_bar(k)=v_bar(k)+v(i,j,k)
            enddo
         enddo             
      enddo
      do k=2,Nzb+1
         do j=1,Ny
            do i=1,Nx
               w_bar(k)=w_bar(k)+w(i,j,k)
            enddo
         enddo             
      enddo
      u_bar=u_bar*inxny
      v_bar=v_bar*inxny
      w_bar=w_bar*inxny
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  
      tu1=0.
      tv1=0.
      tw1=0.
      tp1=0.
      tu2=0.
      tv2=0.
      tw2=0.
      tp2=0.
      tw3=0.
      ttxx=0.
      ttxz=0.
      ttyy=0.
      ttyz=0.
      ttzz=0.
      ttxy=0.
      tuw=0.
      tvw=0.
      tuv=0.
      tdudz=0.
      tdudx=0.
      tdvdz=0.
      tdwdz=0.
      tdwdx=0.
      te=0.
      tCs2=0.
      tCs=0.
      tbeta1=0.
      tESGS=0.

      FSu1=0.
      FSu2=0.
      FSu3=0.
      FSu4=0.	  
      FSv1=0.
      FSv2=0.
      FSv3=0.
      FSv4=0.	  
	  
      do k=2,Nzb+1

         do i=1,aNx

            do j=1,Ny  	
               do l=1,1+Nx-aNx
                  tu1(i,k)=tu1(i,k)+u(i-1+l,j,k)+Ugal
                  tv1(i,k)=tv1(i,k)+v(i-1+l,j,k)+Vgal
                  tw1(i,k)=tw1(i,k)+w(i-1+l,j,k)
                  tp1(i,k)=tp1(i,k)+p(i-1+l,j,k)
                  tp2(i,k)=tp2(i,k)+p(i-1+l,j,k)*p(i-1+l,j,k)
                  ttxx(i,k)=ttxx(i,k)+txx(i-1+l,j,k)
                  ttxz(i,k)=ttxz(i,k)+txz(i-1+l,j,k)
                  ttyy(i,k)=ttyy(i,k)+tyy(i-1+l,j,k)
                  ttyz(i,k)=ttyz(i,k)+tyz(i-1+l,j,k)
                  ttzz(i,k)=ttzz(i,k)+tzz(i-1+l,j,k)
                  ttxy(i,k)=ttxy(i,k)+txy(i-1+l,j,k)
                  tdudz(i,k)=tdudz(i,k)+dudz(i-1+l,j,k)
                  tdudx(i,k)=tdudx(i,k)+dudx(i-1+l,j,k)
                  tdvdz(i,k)=tdvdz(i,k)+dvdz(i-1+l,j,k)
                  tdwdz(i,k)=tdwdz(i,k)+dwdz(i-1+l,j,k)
                  tdwdx(i,k)=tdwdx(i,k)+dwdx(i-1+l,j,k)
                  tCs2(i,k)=tCs2(i,k)+Cs2(i-1+l,j,k)			  
                  tCs(i,k)=tCs(i,k)+abs(Cs2(i-1+l,j,k))**0.5
                  tbeta1(i,k)=tbeta1(i,k)+beta1(i-1+l,j,k)
                  tESGS(i,k)=tESGS(i,k)+ESGS3D(i-1+l,j,k)

                  if (k==2.and.me==0) then
                     arg1=0.
                     arg2=0.
                  else
                     arg1=0.5d0*(u(i-1+l,j,k)+u(i-1+l,j,k-1) - 
     +                           u_bar(k) - u_bar(k-1))
                     arg2=0.5d0*(v(i-1+l,j,k)+v(i-1+l,j,k-1) -
     +                           v_bar(k) - v_bar(k-1))
                  end if
                  
                  tuw(i,k)=tuw(i,k)+(w(i-1+l,j,k) - w_bar(k))*arg1
                  tvw(i,k)=tvw(i,k)+(w(i-1+l,j,k) - w_bar(k))*arg2
                  te(i,k)=te(i,k)+arg1*arg1 + arg2*arg2 +
     +                  (w(i-1+l,j,k) - w_bar(k))**2

                  tu2(i,k)=tu2(i,k)+(u(i-1+l,j,k) - u_bar(k))**2
                  tv2(i,k)=tv2(i,k)+(v(i-1+l,j,k) - v_bar(k))**2
                  tw2(i,k)=tw2(i,k)+(w(i-1+l,j,k) - w_bar(k))**2
                  tw3(i,k)=tw3(i,k)+(w(i-1+l,j,k) - w_bar(k))**3
                  tuv(i,k)=tuv(i,k)+(u(i-1+l,j,k) - u_bar(k))*
     +                              (v(i-1+l,j,k) - v_bar(k))
	 
                  if (arg1.gt.0. .and. (w(i-1+l,j,k)-w_bar(k)).gt.0.) 
     +            FSu1(i,k)=FSu1(i,k)+(w(i-1+l,j,k)-w_bar(k))*arg1	
                  if (arg1.lt.0. .and. (w(i-1+l,j,k)-w_bar(k)).gt.0.) 
     +            FSu2(i,k)=FSu2(i,k)+(w(i-1+l,j,k)-w_bar(k))*arg1	
                  if (arg1.lt.0. .and. (w(i-1+l,j,k)-w_bar(k)).lt.0.) 
     +            FSu3(i,k)=FSu3(i,k)+(w(i-1+l,j,k)-w_bar(k))*arg1	
                  if (arg1.gt.0. .and. (w(i-1+l,j,k)-w_bar(k)).lt.0.) 
     +            FSu4(i,k)=FSu4(i,k)+(w(i-1+l,j,k)-w_bar(k))*arg1		 
                  if (arg2.gt.0. .and. (w(i-1+l,j,k)-w_bar(k)).gt.0.) 
     +            FSv1(i,k)=FSv1(i,k)+(w(i-1+l,j,k)-w_bar(k))*arg2	
                  if (arg2.lt.0. .and. (w(i-1+l,j,k)-w_bar(k)).gt.0.) 
     +            FSv2(i,k)=FSv2(i,k)+(w(i-1+l,j,k)-w_bar(k))*arg2	
                  if (arg2.lt.0. .and. (w(i-1+l,j,k)-w_bar(k)).lt.0.) 
     +            FSv3(i,k)=FSv3(i,k)+(w(i-1+l,j,k)-w_bar(k))*arg2	
                  if (arg2.gt.0. .and. (w(i-1+l,j,k)-w_bar(k)).lt.0.) 
     +            FSv4(i,k)=FSv4(i,k)+(w(i-1+l,j,k)-w_bar(k))*arg2
	 
               enddo
            enddo

         enddo
      enddo

      do kk=2,Nzb+1
cccc NOTE WE ARE JUST USING STEPS HERE (NO INTERPOLATION)
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

            au(i,k)=au(i,k)+ftn*((1-wgx)*tu1(l,k)+wgx*tu1(h,k))
            av(i,k)=av(i,k)+ftn*((1-wgx)*tv1(l,k)+wgx*tv1(h,k))
            aw(i,k)=aw(i,k)+ftn*((1-wgx)*tw1(l,k)+wgx*tw1(h,k))
            ap(i,k)=ap(i,k)+ftn*((1-wgx)*tp1(l,k)+wgx*tp1(h,k))
            u2(i,k)=u2(i,k)+ftn*((1-wgx)*tu2(l,k)+wgx*tu2(h,k))
            v2(i,k)=v2(i,k)+ftn*((1-wgx)*tv2(l,k)+wgx*tv2(h,k))
            w2(i,k)=w2(i,k)+ftn*((1-wgx)*tw2(l,k)+wgx*tw2(h,k))
            p2(i,k)=p2(i,k)+ftn*((1-wgx)*tp2(l,k)+wgx*tp2(h,k))
            w3(i,k)=w3(i,k)+ftn*((1-wgx)*tw3(l,k)+wgx*tw3(h,k))

            atxx(i,k)=atxx(i,k)+ftn*((1-wgx)*ttxx(l,k)+
     +           wgx*ttxx(h,k))
            atxz(i,k)=atxz(i,k)+ftn*((1-wgx)*ttxz(l,k)+
     +           wgx*ttxz(h,k))
            atyy(i,k)=atyy(i,k)+ftn*((1-wgx)*ttyy(l,k)+
     +           wgx*ttyy(h,k))
            atyz(i,k)=atyz(i,k)+ftn*((1-wgx)*ttyz(l,k)+
     +           wgx*ttyz(h,k))
            atzz(i,k)=atzz(i,k)+ftn*((1-wgx)*ttzz(l,k)+
     +           wgx*ttzz(h,k))
            atxy(i,k)=atxy(i,k)+ftn*((1-wgx)*ttxy(l,k)+
     +           wgx*ttxy(h,k))
            auw(i,k)=auw(i,k)+ftn*((1-wgx)*tuw(l,k)+
     +           wgx*tuw(h,k))
            avw(i,k)=avw(i,k)+ftn*((1-wgx)*tvw(l,k)+
     +           wgx*tvw(h,k))
            auv(i,k)=auv(i,k)+ftn*((1-wgx)*tuv(l,k)+
     +           wgx*tuv(h,k))
            adudz(i,k)=adudz(i,k)+ftn*((1-wgx)*tdudz(l,k)+
     +           wgx*tdudz(h,k))
            adudx(i,k)=adudx(i,k)+ftn*((1-wgx)*tdudx(l,k)+
     +           wgx*tdudx(h,k))
            advdz(i,k)=advdz(i,k)+ftn*((1-wgx)*tdvdz(l,k)+
     +           wgx*tdvdz(h,k))
            adwdz(i,k)=adwdz(i,k)+ftn*((1-wgx)*tdwdz(l,k)+
     +           wgx*tdwdz(h,k))
            adwdx(i,k)=adwdx(i,k)+ftn*((1-wgx)*tdwdx(l,k)+
     +           wgx*tdwdx(h,k))
            aCs2(i,k)=aCs2(i,k)+ftn*((1-wgx)*tCs2(l,k)+
     +           wgx*tCs2(h,k))
            aCs(i,k)=aCs(i,k)+ftn*((1-wgx)*tCs(l,k)+
     +           wgx*tCs(h,k))
            abeta1(i,k)=abeta1(i,k)+ftn*((1-wgx)*tbeta1(l,k)+
     +           wgx*tbeta1(h,k))
            aESGS(i,k)=aESGS(i,k)+ftn*((1-wgx)*tESGS(l,k)+
     +           wgx*tESGS(h,k))
            e(i,k)=e(i,k)+ftn*((1-wgx)*te(l,k)+wgx*te(h,k))

			
            aFSu1(i,k)=aFSu1(i,k)+ftn*((1-wgx)*FSu1(l,k)+
     +           wgx*FSu1(h,k))
            aFSu2(i,k)=aFSu2(i,k)+ftn*((1-wgx)*FSu2(l,k)+
     +           wgx*FSu2(h,k))
            aFSu3(i,k)=aFSu3(i,k)+ftn*((1-wgx)*FSu3(l,k)+
     +           wgx*FSu3(h,k))
            aFSu4(i,k)=aFSu4(i,k)+ftn*((1-wgx)*FSu4(l,k)+
     +           wgx*FSu4(h,k))
            aFSv1(i,k)=aFSv1(i,k)+ftn*((1-wgx)*FSv1(l,k)+
     +           wgx*FSv1(h,k))
            aFSv2(i,k)=aFSv2(i,k)+ftn*((1-wgx)*FSv2(l,k)+
     +           wgx*FSv2(h,k))
            aFSv3(i,k)=aFSv3(i,k)+ftn*((1-wgx)*FSv3(l,k)+
     +           wgx*FSv3(h,k))
            aFSv4(i,k)=aFSv4(i,k)+ftn*((1-wgx)*FSv4(l,k)+
     +           wgx*FSv4(h,k))

			
         enddo
      enddo
      
      if(me.eq.0)then
         do j=1,ny
            do i=1,nx
               atxz_s(i,j)=atxz_s(i,j)+fr*
     +              (txz(i,j,2)**2+tyz(i,j,2)**2)**0.25
            enddo
         enddo
      endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine avgslice1(u,v,w,p,txx,txy,txz,tyy,tyz,tzz
     +      ,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
     +      ,dpdx,dpdy,dpdz,Cs2,beta1
     +      ,au1,av1,aw1,ap1,au2,av2,aw2,ap2,auv,auw,avw
     +      ,atxx,atxy,atxz,atyy,atyz,atzz,adpdx,adpdy,adpdz
     +      ,adudx,adudy,adudz,advdx,advdy,advdz,adwdx,adwdy,adwdz
     +      ,abeta1,aCs,aCs2,fx,fy,fz,afx,afy,afz
     +      ,fax,faz,AveWt,afax,afaz,aAveWt,aufx,avfy,awfz
     +      ,aTt,aTpz,aTsgsz,aTpx,aTsgsx,aTpy,aTsgsy 
     +      ,aS11,aS22,aS33,aS12,aS13,aS23
     +      ,aeps11,aeps22,aeps33,aeps12,aeps13,aeps23	 
     +      ,me)
      
      implicit none
      include 'dimen.h'
      integer*4 tt,i,j,k,l,kk,ilow(nx),h
      real*8,dimension(nx,ny,nz2):: u,v,w,p,txx,txy,txz,tyy,tyz,tzz
     +      ,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
     +      ,dpdx,dpdy,dpdz,Cs2,beta1,fx,fy,fz
      
      real*8,dimension(nx,ny,nz2) ::
     +       au1,av1,aw1,ap1,au2,av2,aw2,ap2,auv,auw,avw
     +      ,atxx,atxy,atxz,atyy,atyz,atzz,adpdx,adpdy,adpdz
     +      ,adudx,adudy,adudz,advdx,advdy,advdz,adwdx,adwdy,adwdz
     +      ,abeta1,aCs,aCs2,afx,afy,afz
     +      ,fax,faz,AveWt,afax,afaz,aAveWt,aufx,avfy,awfz
     +      ,aTt,aTpz,aTsgsz,aTpx,aTsgsx,aTpy,aTsgsy 
     +      ,S11,S22,S33,S12,S13,S23
     +      ,aS11,aS22,aS33,aS12,aS13,aS23	
     +      ,aeps11,aeps22,aeps33,aeps12,aeps13,aeps23
     +      ,u_new	 
      
      real*8 fr,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8
      real*8 u_line(1,ny,nz2), u_plane(nz2),u_span(nx,nz2)
      
      fr = dble(c_count)/dble(p_count)

      au1=au1+fr*u
      av1=av1+fr*v
      aw1=aw1+fr*w
      ap1=ap1+fr*p
      
      au2=au2+fr*u*u 
      av2=av2+fr*v*v
      aw2=aw2+fr*w*w 
      ap2=ap2+fr*p*p

      auv=auv+fr*u*v

      adpdx=adpdx+fr*dpdx
      adpdy=adpdy+fr*dpdy
      adpdz=adpdz+fr*dpdz
      
      atxx=atxx+fr*txx
      atxy=atxy+fr*txy
      atxz=atxz+fr*txz
      atyy=atyy+fr*tyy
      atyz=atyz+fr*tyz
      atzz=atzz+fr*tzz
      
      adudx=adudx+fr*dudx
      adudy=adudy+fr*dudy
      adudz=adudz+fr*dudz
      advdx=advdx+fr*dvdx
      advdy=advdy+fr*dvdy
      advdz=advdz+fr*dvdz
      adwdx=adwdx+fr*dwdx
      adwdy=adwdy+fr*dwdy
      adwdz=adwdz+fr*dwdz

      aCs=aCs+fr*dsqrt(abs(Cs2))
      aCs2=aCs2+fr*Cs2
      abeta1=abeta1+fr*beta1
      
      afx=afx+fr*fx
      afy=afy+fr*fy
      afz=afz+fr*fz
      afax=afax+fr*fax
      afaz=afaz+fr*faz
      
      aAveWt=aAveWt+fr*AveWt
      
      aufx=aufx+fr*u*fx
      avfy=avfy+fr*v*fy
      awfz=awfz+fr*w*fz	  

      S11=dudx
      S22=dvdy
      S33=dwdz
      S12=0.5*(dudy+dvdx)
      S13=0.5*(dudz+dwdx)
      S23=0.5*(dvdz+dwdy)
	  
      aS11=aS11+fr*S11
      aS22=aS22+fr*S22
      aS33=aS33+fr*S33
      aS12=aS12+fr*S12
      aS13=aS13+fr*S13
      aS23=aS23+fr*S23
	  
      aeps11=aeps11+fr*txx*S11
      aeps22=aeps22+fr*tyy*S22
      aeps33=aeps33+fr*tzz*S33
      aeps12=aeps12+fr*txy*S12
      aeps13=aeps13+fr*txz*S13
      aeps23=aeps23+fr*tyz*S23	  	  

      do k=2,Nzb+1
      do j=1,Ny
      do i=1,Nx
      if (k==2.and.me==0) then
         arg1=0.
         arg2=0.
         arg3=0.
         arg4=0.
         arg5=0.
         arg6=0.		 
         arg7=0.		 
      else
         arg1=(u(i,j,k)+u(i,j,k-1))/2.+Ugal   
         arg2=(v(i,j,k)+v(i,j,k-1))/2.
         arg3=(p(i,j,k)+p(i,j,k-1))/2.
         arg4=(tzz(i,j,k)+tzz(i,j,k-1))/2.
         arg5=(w(i,j,k)+w(i,j,k+1))/2.
         arg6=(txz(i,j,k)+txz(i,j,k+1))/2.	
         arg7=(tyz(i,j,k)+tyz(i,j,k+1))/2.		 
      end if
      auw(i,j,k)=auw(i,j,k)+fr*w(i,j,k)*arg1  
      avw(i,j,k)=avw(i,j,k)+fr*w(i,j,k)*arg2
      aTt(i,j,k)=aTt(i,j,k)+fr*(w(i,j,k)*arg1*arg1+
     +                          w(i,j,k)*arg2*arg2+	 
     +                          w(i,j,k)*w(i,j,k)*w(i,j,k))
      aTpz(i,j,k)=aTpz(i,j,k)+fr*w(i,j,k)*arg3
      aTpx(i,j,k)=aTpx(i,j,k)+fr*u(i,j,k)*p(i,j,k)
      aTpy(i,j,k)=aTpy(i,j,k)+fr*v(i,j,k)*p(i,j,k)		  
      aTsgsz(i,j,k)=aTsgsz(i,j,k)+fr*(txz(i,j,k)*arg1+
     +                              tyz(i,j,k)*arg2+
     +                              w(i,j,k)*arg4)
      aTsgsx(i,j,k)=aTsgsx(i,j,k)+fr*(txx(i,j,k)*u(i,j,k)+
     +                              txy(i,j,k)*v(i,j,k)+
     +                              arg5*arg6)
      aTsgsy(i,j,k)=aTsgsy(i,j,k)+fr*(txy(i,j,k)*u(i,j,k)+
     +                              tyy(i,j,k)*v(i,j,k)+
     +                              arg5*arg7) 	 
      end do
      end do
      end do
  
      return
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      fr = (1./p_count)*c_count
      
      do k=2,Nzb+1
      do j=1,Ny
      do i=1,Nx
      au1(i,j,k)=au1(i,j,k)+fr*u(i,j,k)
      av1(i,j,k)=av1(i,j,k)+fr*v(i,j,k)
      aw1(i,j,k)=aw1(i,j,k)+fr*w(i,j,k)
      ap1(i,j,k)=ap1(i,j,k)+fr*p(i,j,k)
      
      au2(i,j,k)=au2(i,j,k)+fr*(u(i,j,k)**2.) 
      av2(i,j,k)=av2(i,j,k)+fr*(v(i,j,k)**2.)
      aw2(i,j,k)=aw2(i,j,k)+fr*(w(i,j,k)**2.) 
      ap2(i,j,k)=ap2(i,j,k)+fr*(p(i,j,k)**2.)
      
      if (k==2.and.me==0) then
         arg1=0.
         arg2=0.
      else
         arg1=(u(i,j,k)+u(i,j,k-1))/2.+Ugal   
         arg2=(v(i,j,k)+v(i,j,k-1))/2.
      end if
      
      auv(i,j,k)=auv(i,j,k)+fr*u(i,j,k)*v(i,j,k)
      auw(i,j,k)=auw(i,j,k)+fr*w(i,j,k)*arg1  
      avw(i,j,k)=avw(i,j,k)+fr*w(i,j,k)*arg2  
      
      adpdx(i,j,k)=adpdx(i,j,k)+fr*dpdx(i,j,k)
      adpdy(i,j,k)=adpdy(i,j,k)+fr*dpdy(i,j,k)
      adpdz(i,j,k)=adpdz(i,j,k)+fr*dpdz(i,j,k)
      
      atxx(i,j,k)=atxx(i,j,k)+fr*txx(i,j,k)
      atxy(i,j,k)=atxy(i,j,k)+fr*txy(i,j,k)
      atxz(i,j,k)=atxz(i,j,k)+fr*txz(i,j,k)
      atyy(i,j,k)=atyy(i,j,k)+fr*tyy(i,j,k)
      atyz(i,j,k)=atyz(i,j,k)+fr*tyz(i,j,k)
      atzz(i,j,k)=atzz(i,j,k)+fr*tzz(i,j,k)
      
      adudx(i,j,k)=adudx(i,j,k)+fr*dudx(i,j,k)
      adudy(i,j,k)=adudy(i,j,k)+fr*dudy(i,j,k)
      adudz(i,j,k)=adudz(i,j,k)+fr*dudz(i,j,k)
      advdx(i,j,k)=advdx(i,j,k)+fr*dvdx(i,j,k)
      advdy(i,j,k)=advdy(i,j,k)+fr*dvdy(i,j,k)  
      advdz(i,j,k)=advdz(i,j,k)+fr*dvdz(i,j,k)
      adwdx(i,j,k)=adwdx(i,j,k)+fr*dwdx(i,j,k)
      adwdy(i,j,k)=adwdy(i,j,k)+fr*dwdy(i,j,k)
      adwdz(i,j,k)=adwdz(i,j,k)+fr*dwdz(i,j,k)
      
      aCs (i,j,k)=aCs (i,j,k)+fr*(Cs2(i,j,k)**0.5)
      aCs2(i,j,k)=aCs2(i,j,k)+fr*Cs2(i,j,k)
      abeta1(i,j,k)=abeta1(i,j,k)+fr*beta1(i,j,k)
      
      afx(i,j,k)=afx(i,j,k)+fr*fx(i,j,k)   
      afy(i,j,k)=afy(i,j,k)+fr*fy(i,j,k)   
      afz(i,j,k)=afz(i,j,k)+fr*fz(i,j,k)   
      afax(i,j,k)=afax(i,j,k)+fr*fax(i,j,k)
      afaz(i,j,k)=afaz(i,j,k)+fr*faz(i,j,k)
      
      aAveWt(i,j,k)=aAveWt(i,j,k)+fr*AveWt(i,j,k)   
      end do
      end do
      end do

      return
      
      end
	  
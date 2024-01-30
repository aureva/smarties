!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine AMD_Scalar(dudx,dudy,dudz,dvdx,dvdy,dvdz,
     &   dwdx,dwdy,dwdz,t_flux,
     &   u,v,w,txx,txy,txz,tyy,tyz,tzz,
     &   dtdx,dtdy,dtdz,qx,qy,qz,t,me,nall,
     &   txx_m,txy_m,txz_m,tyy_m,tyz_m,tzz_m,
     &   S11,S12,S13,S22,S23,S33,
     &   S11_m,S12_m,S13_m,S22_m,S23_m,S33_m,
     &   ux,uy,uz,vx,vy,vz,wx,wy,wz,
     &   tx,ty,tz,tx_m,ty_m,tz_m,qx_m,qy_m,qz_m,Cs2,Pr2)
      include 'dimen.h'

! these are inputs and outputs
      real*8,dimension(nx,ny,nz2)::
     &   dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,
     &   u,v,w,txx,txy,txz,tyy,tyz,tzz,dtdx,dtdy,dtdz,qx,qy,qz,
     &   Cs2,Pr2	 
      integer*4 t, countmgm1,countmgm2,countmgm3

      real*8 t_flux(nx,ny),tx_mean,ty_mean,tz_mean

! these are temps
      real*8,dimension(nx,ny,nz2)::
     &    S11,S12,S13,S22,S23,S33,tx,ty,tz
      real*8,dimension(nx2,ny2,nz2)::
     &    txx_m,txy_m,txz_m,tyy_m,tyz_m,tzz_m,
     &    S11_m,S12_m,S13_m,S22_m,S23_m,S33_m,
     &    tx_m,ty_m,tz_m,qx_m,qy_m,qz_m,
     &    Gkk_m,Gtt_m,Cs2_m,Pr2_m,bt_m

! these are memory temps
      real*8,dimension(nx,ny,nz2)::ux,uy,uz,vx,vy,vz,wx,wy,wz

      integer ix,iy,iz,izstart
      real*8 ux_,uy_,uz_,vx_,vy_,vz_,wx_,wy_,wz_,tx_,ty_,tz_
      
      real*8 factor, Ksgs,Cs
      real*8 tmp1, tmp2, tmp3, tmp4, tmp5, rzx2_, zz, usgs, ssgs, G_
      real*8 S, factor_s

      real*8, ALLOCATABLE:: txzp(:,:)
      real*8, ALLOCATABLE:: tyzp(:,:)

      real*8, save:: ryx2,rzx2,d_,Ce_const,Ces,eps,Cx,Cz,Cy
      integer, save:: i_ini
      
 
! 0. prepare
      if (i_ini .eq. 0) then
         d_ = (dx*dy*dz)**(1.0d0/3.0d0)

!         Cs = 0.1d0
         eps = 1.0d-6

         Ce_const = 0.132d0	
         Ces = 0.132d0		 

!         Cx= 0.10132d0
!         Cy= 0.10132d0		 
!         Cz= 0.33333d0
	
         Cx= 1./12.d0
         Cy= 1./12.d0	 
         Cz= 1./3.0d0
		 
         i_ini = 1
      end if

      if (me == 0) then
! Stress Boundary Condition at bottom has been given, 
! thus, it needs to be saved and restored. Here is for saving BC
         allocate(txzp(Nx,Ny))
         allocate(tyzp(Nx,Ny))
         do iy=1,Ny
         do ix=1,Nx
            txzp(ix,iy)=txz(ix,iy,2)
            tyzp(ix,iy)=tyz(ix,iy,2)
         end do
         end do

! the only reason doing that is: only at bottom(me=0,iz=2), 
! dudz and dvdz are computed at u layer
         iz = 2
         do iy = 1, Ny
         do ix = 1, Nx
            ux(ix,iy,iz) = dudx(ix,iy,iz)
            uy(ix,iy,iz) = dudy(ix,iy,iz)
            uz(ix,iy,iz) = dudz(ix,iy,iz)
            vx(ix,iy,iz) = dvdx(ix,iy,iz)
            vy(ix,iy,iz) = dvdy(ix,iy,iz)
            vz(ix,iy,iz) = dvdz(ix,iy,iz)
            wx(ix,iy,iz) = 0.5d0*(dwdx(ix,iy,iz)+dwdx(ix,iy,iz+1))
            wy(ix,iy,iz) = 0.5d0*(dwdy(ix,iy,iz)+dwdy(ix,iy,iz+1))
            wz(ix,iy,iz) = dwdz(ix,iy,iz)
!            wz(ix,iy,iz) = -(ux(ix,iy,iz)+vy(ix,iy,iz))
         end do
         end do
         
         izstart = 3
      else
         izstart = 2
      end if
      
      do iz = izstart, Nzb+1
         do iy = 1, Ny
         do ix = 1, Nx
            ux(ix,iy,iz) = dudx(ix,iy,iz)
            uy(ix,iy,iz) = dudy(ix,iy,iz)
            uz(ix,iy,iz) = 0.5d0*(dudz(ix,iy,iz)+dudz(ix,iy,iz+1))
            vx(ix,iy,iz) = dvdx(ix,iy,iz)
            vy(ix,iy,iz) = dvdy(ix,iy,iz)
            vz(ix,iy,iz) = 0.5d0*(dvdz(ix,iy,iz)+dvdz(ix,iy,iz+1))
            wx(ix,iy,iz) = 0.5d0*(dwdx(ix,iy,iz)+dwdx(ix,iy,iz+1))
            wy(ix,iy,iz) = 0.5d0*(dwdy(ix,iy,iz)+dwdy(ix,iy,iz+1))
            wz(ix,iy,iz) = dwdz(ix,iy,iz)
!            wz(ix,iy,iz) = -(ux(ix,iy,iz)+vy(ix,iy,iz))

         end do
         end do
      end do

! 2. get tau

! 2.1 u, v layer
      call dealias1(ux,S11_m,t,0)
      call dealias1(uy,S22_m,t,0)
      call dealias1(uz,S33_m,t,0)
      call dealias1(vx,S12_m,t,0)
      call dealias1(vy,S13_m,t,0)
      call dealias1(vz,S23_m,t,0)
      call dealias1(wx,txz_m,t,0)
      call dealias1(wy,tyz_m,t,0)
      call dealias1(wz,tzz_m,t,0)

      if (me .eq. 0) then
         iz = 2
         do iy = 1, Ny
         do ix = 1, Nx
            uz(ix,iy,iz) = dtdz(ix,iy,iz)
         end do
         end do

         do iz = 3, Nzb+1
         do iy = 1, Ny
         do ix = 1, Nx
            uz(ix,iy,iz) = 0.5d0*(dtdz(ix,iy,iz)+dtdz(ix,iy,iz+1))
         end do
         end do
         end do
      else
         do iz = 2, Nzb+1
         do iy = 1, Ny
         do ix = 1, Nx
            uz(ix,iy,iz) = 0.5d0*(dtdz(ix,iy,iz)+dtdz(ix,iy,iz+1))
         end do
         end do
         end do
      end if
      call dealias1(dtdx,tx_m,t,0)
      call dealias1(dtdy,ty_m,t,0)
      call dealias1(uz,tz_m,t,0)
  
      do iz = 2, Nzb+1
	  
      tx_mean=sum(tx_m(:,:,iz))/(Nx2*Ny2)
      ty_mean=sum(ty_m(:,:,iz))/(Nx2*Ny2)
      tz_mean=sum(tz_m(:,:,iz))/(Nx2*Ny2)

      do iy = 1, Ny2
      do ix = 1, Nx2
         ux_ = S11_m(ix,iy,iz)
         uy_ = S22_m(ix,iy,iz)
         uz_ = S33_m(ix,iy,iz)
         vx_ = S12_m(ix,iy,iz)
         vy_ = S13_m(ix,iy,iz)
         vz_ = S23_m(ix,iy,iz)
         wx_ = txz_m(ix,iy,iz)
         wy_ = tyz_m(ix,iy,iz)
         wz_ = tzz_m(ix,iy,iz)
!         wz_ = -(ux_+vy_)

         tx_ = tx_m(ix,iy,iz)
         ty_ = ty_m(ix,iy,iz)
         tz_ = tz_m(ix,iy,iz)

       txx_m(ix,iy,iz) = Cx*(dx*ux_)**2+Cy*(dy*uy_)**2+Cz*(dz*uz_)**2
       tyy_m(ix,iy,iz) = Cx*(dx*vx_)**2+Cy*(dy*vy_)**2+Cz*(dz*vz_)**2
       tzz_m(ix,iy,iz) = Cx*(dx*wx_)**2+Cy*(dy*wy_)**2+Cz*(dz*wz_)**2
       txy_m(ix,iy,iz) = Cx*(dx*ux_)*(dx*vx_)+Cy*(dy*uy_)*(dy*vy_)+
     +                   Cz*(dz*uz_)*(dz*vz_)	 
       txz_m(ix,iy,iz) = Cx*(dx*ux_)*(dx*wx_)+Cy*(dy*uy_)*(dy*wy_)+
     +                   Cz*(dz*uz_)*(dz*wz_)		   
       tyz_m(ix,iy,iz) = Cx*(dx*vx_)*(dx*wx_)+Cy*(dy*vy_)*(dy*wy_)+
     +                   Cz*(dz*vz_)*(dz*wz_)	
 
       Gkk_m(ix,iy,iz) = ux_**2+vx_**2+wx_**2+
     +                   uy_**2+vy_**2+wy_**2+
     +                   uz_**2+vz_**2+wz_**2 

         S11_m(ix,iy,iz) = ux_
         S22_m(ix,iy,iz) = vy_
         S33_m(ix,iy,iz) = wz_
         S12_m(ix,iy,iz) = 0.5d0*(uy_+vx_)
         S13_m(ix,iy,iz) = 0.5d0*(uz_+wx_)
         S23_m(ix,iy,iz) = 0.5d0*(vz_+wy_)

         qx_m(ix,iy,iz) =  Cx*(dx*ux_)*(dx*tx_)+Cy*(dy*uy_)*(dy*ty_)+
     +                     Cz*(dz*uz_)*(dz*tz_)
         qy_m(ix,iy,iz) =  Cx*(dx*vx_)*(dx*tx_)+Cy*(dy*vy_)*(dy*ty_)+
     +                     Cz*(dz*vz_)*(dz*tz_)
         qz_m(ix,iy,iz) =  Cx*(dx*wx_)*(dx*tx_)+Cy*(dy*wy_)*(dy*ty_)+
     +                     Cz*(dz*wz_)*(dz*tz_)
	 
         Gtt_m(ix,iy,iz) = tx_**2+ty_**2+tz_**2
		 
        bt_m(ix,iy,iz)=g_hat/theta_0*(
     +		Cx*(dx*wx_)*(dx*(tx_-tx_mean))+
     +      Cy*(dy*wy_)*(dy*(ty_-ty_mean))+
     +      Cz*(dz*wz_)*(dz*(tz_-tz_mean)))	
	 
!        bt_m(ix,iy,iz)=g_hat/theta_0*qz_m(ix,iy,iz)		
		
      end do
      end do
      end do

      tmp1 = (2.0d0*d_/Ce_const)**2
      do iz = 2, Nzb+1
      zz = (iz-1.5d0+me*Nzb)*dz
      do iy = 1, Ny2
      do ix = 1, Nx2
	  
         Gkk =  Gkk_m(ix,iy,iz)
         tmp3 = txx_m(ix,iy,iz)*S11_m(ix,iy,iz)+tyy_m(ix,iy,iz)*
     &          S22_m(ix,iy,iz)+tzz_m(ix,iy,iz)*S33_m(ix,iy,iz)+
     &          2.d0*(txy_m(ix,iy,iz)*S12_m(ix,iy,iz)+txz_m(ix,iy,iz)*
     &          S13_m(ix,iy,iz)+tyz_m(ix,iy,iz)*S23_m(ix,iy,iz))-
     &	        (1-pass_flag)*bt_m(ix,iy,iz) 

            Ksgs = max(-tmp3,0.0d0)

            if (Gkk .ne. 0.0d0) then			
            tmp3 = 1.0d0/Gkk	
            else			
            tmp3 = 0.0d0	
            end if
			
            factor = Ksgs*tmp3+nu
			
            S=S11_m(ix,iy,iz)**2+S22_m(ix,iy,iz)**2+S33_m(ix,iy,iz)**2+
     &      2*(S12_m(ix,iy,iz)**2+S13_m(ix,iy,iz)**2+S23_m(ix,iy,iz)**2)
            S = dsqrt(2*S)
            if (S.ne.0.d0) then 			
            Cs=(factor-nu)/((d_**2)*S)	
            else			 
            Cs=0.0d0
            endif
            if (Cs.gt.0.36d0) Cs=0.0d0	!Cs_max=0.6		
!            factor = Cs*(d_**2)*S+nu
				
            txx_m(ix,iy,iz) = -2.0d0*factor*S11_m(ix,iy,iz)
            tyy_m(ix,iy,iz) = -2.0d0*factor*S22_m(ix,iy,iz)
            tzz_m(ix,iy,iz) = -2.0d0*factor*S33_m(ix,iy,iz)
            txy_m(ix,iy,iz) = -2.0d0*factor*S12_m(ix,iy,iz)
			
!            Cs2_m(ix,iy,iz)=(factor-nu)/((d_**2)*S)	
            Cs2_m(ix,iy,iz)=(factor-nu)			  
			
            G_ =  Gtt_m(ix,iy,iz)	  
            tmp4 = qx_m(ix,iy,iz)*tx_m(ix,iy,iz)+
     &      qy_m(ix,iy,iz)*ty_m(ix,iy,iz)+qz_m(ix,iy,iz)*tz_m(ix,iy,iz)
	 
            ssgs = max(-tmp4,0.0d0)	

            if (G_ .ne. 0.0d0) then 			
            tmp3 = 1.0d0/G_	
            else			
            tmp3 = 0.0d0	
            end if
			
            factor_s = ssgs*tmp3+nu/Pr

            if (S.ne.0.d0) then 			
            Cs=(factor_s-nu/Pr)/((d_**2)*S)	
            else			 
            Cs=0.0d0
            endif
            if (Cs.gt.0.9d0) Cs=0.0d0	!Cs2Pr-1_max=0.9		
!            factor_s = Cs*(d_**2)*S+nu/Pr			
			
            qx_m(ix,iy,iz) = -1.0d0*factor_s*tx_m(ix,iy,iz)
            qy_m(ix,iy,iz) = -1.0d0*factor_s*ty_m(ix,iy,iz)	  

!            Pr2_m(ix,iy,iz)=(factor_s-nu/Pr)/((d_**2)*S)
            Pr2_m(ix,iy,iz)=(factor_s-nu/Pr)					
 
      end do
      end do
      end do

      call dealias2(txx,txx_m,t,0)
      call dealias2(tyy,tyy_m,t,0)
      call dealias2(tzz,tzz_m,t,0)
      call dealias2(txy,txy_m,t,0)
	  
      call dealias2(Cs2,Cs2_m,t,0)		  
      call dealias2(Pr2,Pr2_m,t,0)
	  
      call update3(txx,tyy,tzz,me,nall)
      call update1(txy,me,nall)

      call dealias2(qx,qx_m,t,0)
      call dealias2(qy,qy_m,t,0)

! 2.2 w layer
      do iz = 2, Nzb+1
      do iy = 1, Ny
      do ix = 1, Nx
         ux(ix,iy,iz) = 0.5d0*(dudx(ix,iy,iz-1)+dudx(ix,iy,iz))
         uy(ix,iy,iz) = 0.5d0*(dudy(ix,iy,iz-1)+dudy(ix,iy,iz))
         uz(ix,iy,iz) = dudz(ix,iy,iz)
         vx(ix,iy,iz) = 0.5d0*(dvdx(ix,iy,iz-1)+dvdx(ix,iy,iz))
         vy(ix,iy,iz) = 0.5d0*(dvdy(ix,iy,iz-1)+dvdy(ix,iy,iz))
         vz(ix,iy,iz) = dvdz(ix,iy,iz)
         wx(ix,iy,iz) = dwdx(ix,iy,iz)
         wy(ix,iy,iz) = dwdy(ix,iy,iz)
         wz(ix,iy,iz) = 0.5d0*(dwdz(ix,iy,iz-1)+dwdz(ix,iy,iz))
!         wz(ix,iy,iz) = -(ux_+vy_)
      end do
      end do
      end do

      call dealias1(ux,S11_m,t,0)
      call dealias1(uy,S22_m,t,0)
      call dealias1(uz,S33_m,t,0)
      call dealias1(vx,S12_m,t,0)
      call dealias1(vy,S13_m,t,0)
      call dealias1(vz,S23_m,t,0)
      call dealias1(wx,txz_m,t,0)
      call dealias1(wy,tyz_m,t,0)
      call dealias1(wz,tzz_m,t,0)

      do iz = 2, Nzb+1
      do iy = 1, Ny
      do ix = 1, Nx
         ux(ix,iy,iz) = 0.5d0*(dtdx(ix,iy,iz-1)+dtdx(ix,iy,iz))
         uy(ix,iy,iz) = 0.5d0*(dtdy(ix,iy,iz-1)+dtdy(ix,iy,iz))
      end do
      end do
      end do
      call dealias1(ux,tx_m,t,0)
      call dealias1(uy,ty_m,t,0)
      call dealias1(dtdz,tz_m,t,0)

      do iz = 2, Nzb+1
	  
      tx_mean=sum(tx_m(:,:,iz))/(Nx2*Ny2)
      ty_mean=sum(ty_m(:,:,iz))/(Nx2*Ny2)
      tz_mean=sum(tz_m(:,:,iz))/(Nx2*Ny2)	  
	  
      do iy = 1, Ny2
      do ix = 1, Nx2
         ux_ = S11_m(ix,iy,iz)
         uy_ = S22_m(ix,iy,iz)
         uz_ = S33_m(ix,iy,iz)
         vx_ = S12_m(ix,iy,iz)
         vy_ = S13_m(ix,iy,iz)
         vz_ = S23_m(ix,iy,iz)
         wx_ = txz_m(ix,iy,iz)
         wy_ = tyz_m(ix,iy,iz)
         wz_ = tzz_m(ix,iy,iz)
!         wz_ = -(ux_+vy_)

         tx_ = tx_m(ix,iy,iz)
         ty_ = ty_m(ix,iy,iz)
         tz_ = tz_m(ix,iy,iz)

       txx_m(ix,iy,iz) = Cx*(dx*ux_)**2+Cy*(dy*uy_)**2+Cz*(dz*uz_)**2
       tyy_m(ix,iy,iz) = Cx*(dx*vx_)**2+Cy*(dy*vy_)**2+Cz*(dz*vz_)**2
       tzz_m(ix,iy,iz) = Cx*(dx*wx_)**2+Cy*(dy*wy_)**2+Cz*(dz*wz_)**2
       txy_m(ix,iy,iz) = Cx*(dx*ux_)*(dx*vx_)+Cy*(dy*uy_)*(dy*vy_)+
     +                   Cz*(dz*uz_)*(dz*vz_)	 
       txz_m(ix,iy,iz) = Cx*(dx*ux_)*(dx*wx_)+Cy*(dy*uy_)*(dy*wy_)+
     +                   Cz*(dz*uz_)*(dz*wz_)		   
       tyz_m(ix,iy,iz) = Cx*(dx*vx_)*(dx*wx_)+Cy*(dy*vy_)*(dy*wy_)+
     +                   Cz*(dz*vz_)*(dz*wz_)	
 
       Gkk_m(ix,iy,iz) = ux_**2+vx_**2+wx_**2+
     +                   uy_**2+vy_**2+wy_**2+
     +                   uz_**2+vz_**2+wz_**2 

         S11_m(ix,iy,iz) = ux_
         S22_m(ix,iy,iz) = vy_
         S33_m(ix,iy,iz) = wz_
         S12_m(ix,iy,iz) = 0.5d0*(uy_+vx_)
         S13_m(ix,iy,iz) = 0.5d0*(uz_+wx_)
         S23_m(ix,iy,iz) = 0.5d0*(vz_+wy_)

         qx_m(ix,iy,iz) =  Cx*(dx*ux_)*(dx*tx_)+Cy*(dy*uy_)*(dy*ty_)+
     +                     Cz*(dz*uz_)*(dz*tz_)
         qy_m(ix,iy,iz) =  Cx*(dx*vx_)*(dx*tx_)+Cy*(dy*vy_)*(dy*ty_)+
     +                     Cz*(dz*vz_)*(dz*tz_)
         qz_m(ix,iy,iz) =  Cx*(dx*wx_)*(dx*tx_)+Cy*(dy*wy_)*(dy*ty_)+
     +                     Cz*(dz*wz_)*(dz*tz_)
	 
         Gtt_m(ix,iy,iz) = tx_**2+ty_**2+tz_**2
		 	
        bt_m(ix,iy,iz)=g_hat/theta_0*(
     +		Cx*(dx*wx_)*(dx*(tx_-tx_mean))+
     +      Cy*(dy*wy_)*(dy*(ty_-ty_mean))+
     +      Cz*(dz*wz_)*(dz*(tz_-tz_mean)))			 

!        bt_m(ix,iy,iz)=g_hat/theta_0*qz_m(ix,iy,iz)	 
	 
      end do
      end do
      end do

      tmp1 = (2.0d0*d_/Ce_const)**2
      do iz = 2, Nzb+1
      zz = (iz-1.5d0+me*Nzb)*dz
      do iy = 1, Ny2
      do ix = 1, Nx2
	  
         Gkk =  Gkk_m(ix,iy,iz)
         tmp3 = txx_m(ix,iy,iz)*S11_m(ix,iy,iz)+tyy_m(ix,iy,iz)*
     &          S22_m(ix,iy,iz)+tzz_m(ix,iy,iz)*S33_m(ix,iy,iz)+
     &          2.d0*(txy_m(ix,iy,iz)*S12_m(ix,iy,iz)+txz_m(ix,iy,iz)*
     &          S13_m(ix,iy,iz)+tyz_m(ix,iy,iz)*S23_m(ix,iy,iz))-
     &	        (1-pass_flag)*bt_m(ix,iy,iz)  

            Ksgs = max(-tmp3,0.0d0)
            if (Gkk .ne. 0.0d0) then			
            tmp3 = 1.0d0/Gkk	
            else			
            tmp3 = 0.0d0	
            end if
			
            factor = Ksgs*tmp3+nu
			
            S=S11_m(ix,iy,iz)**2+S22_m(ix,iy,iz)**2+S33_m(ix,iy,iz)**2+
     &      2*(S12_m(ix,iy,iz)**2+S13_m(ix,iy,iz)**2+S23_m(ix,iy,iz)**2)
            S = dsqrt(2*S)
            if (S.ne.0.d0) then 			
            Cs=(factor-nu)/((d_**2)*S)	
            else			 
            Cs=0.0d0
            endif
            if (Cs.gt.0.36d0) Cs=0.0d0	!Cs_max=0.6		
!            factor = Cs*(d_**2)*S+nu
			
            txz_m(ix,iy,iz) = -2.0d0*factor*S13_m(ix,iy,iz)
            tyz_m(ix,iy,iz) = -2.0d0*factor*S23_m(ix,iy,iz)
	  
            G_ =  Gtt_m(ix,iy,iz)	  
            tmp4 = qx_m(ix,iy,iz)*tx_m(ix,iy,iz)+
     &      qy_m(ix,iy,iz)*ty_m(ix,iy,iz)+qz_m(ix,iy,iz)*tz_m(ix,iy,iz)
	 
            ssgs = max(-tmp4,0.0d0)	
            if (G_ .ne. 0.0d0) then			
            tmp3 = 1.0d0/G_	
            else			
            tmp3 = 0.0d0	
            end if		
			
            factor_s = ssgs*tmp3+nu/Pr
			
            if (S.ne.0.d0) then 			
            Cs=(factor_s-nu/Pr)/((d_**2)*S)	
            else			 
            Cs=0.0d0
            endif
            if (Cs.gt.0.9d0) Cs=0.0d0	!Cs2Pr-1_max=0.9		
!            factor_s = Cs*(d_**2)*S+nu/Pr			
				
            qz_m(ix,iy,iz) = -1.0d0*factor_s*tz_m(ix,iy,iz)
 
      end do
      end do
      end do

      call dealias2(txz,txz_m,t,0)
      call dealias2(tyz,tyz_m,t,0)

      call dealias2(qz,qz_m,t,0)

      if (me==0) then
! Stress Boundary Condition at bottom has been given, 
! thus, it needs to be saved and restored. Here is for resotration
         do iy=1,Ny
         do ix=1,Nx
            txz(ix,iy,2)=txzp(ix,iy)
            tyz(ix,iy,2)=tyzp(ix,iy)
            qz(ix,iy,2) = t_flux(ix,iy)
         end do
         end do
         deallocate(txzp)
         deallocate(tyzp)
      end if

      call update1(txz,me,nall)
      call update1(tyz,me,nall)

      call update3(qx,qy,qz,me,nall)

!      write(*,*) 'countmgm1,2,3=',countmgm1,countmgm2,countmgm3        
      end
!-----------------------------------------------------------------------


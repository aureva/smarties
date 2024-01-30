!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine AMD_Const(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,
     &   txx,txy,txz,tyy,tyz,tzz,t,me,nall,
     &   txx_m,txy_m,txz_m,tyy_m,tyz_m,tzz_m,
     &   S11,S12,S13,S22,S23,S33,
     &   S11_m,S12_m,S13_m,S22_m,S23_m,S33_m,
     &   ux,uy,uz,vx,vy,vz,wx,wy,wz,Cs2,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + dx,dy,dz,nu,plan_f,plan_b,plan_ff,plan_bb)

      implicit none
!      include 'dimen.h'

! these are inputs and outputs
      integer*4 :: nx,ny,nzb,nx2,ny2,nz2,nsteps,me,nall
      real*8 :: inxny,inx2ny2,dx,dy,dz,nu
      integer*8 :: plan_f,plan_b,plan_ff,plan_bb

      real*8,dimension(nx,ny,nz2)::
     &   dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,
     &   txx,txy,txz,tyy,tyz,tzz,Cs2
      integer*4 t,countmgm1,countmgm2,countmgm3

! these are memory temps
      real*8,dimension(nx,ny,nz2)::
     &    S11,S12,S13,S22,S23,S33
      real*8,dimension(nx2,ny2,nz2)::
     &    txx_m,txy_m,txz_m,tyy_m,tyz_m,tzz_m,
     &    S11_m,S12_m,S13_m,S22_m,S23_m,S33_m,Gkk_m,Cs2_m
      real*8,dimension(nx,ny,nz2)::ux,uy,uz,vx,vy,vz,wx,wy,wz

      integer ix,iy,iz,izstart
      real*8 ux_,uy_,uz_,vx_,vy_,vz_,wx_,wy_,wz_
      
      real*8 factor, Ksgs, Cs
      real*8 tmp1, tmp2, tmp3, tmp4, tmp5, Gkk, S

      real*8, ALLOCATABLE:: txzp(:,:)
      real*8, ALLOCATABLE:: tyzp(:,:)

      real*8, save:: ryx2,rzx2,d_,Ce_const,eps,Cx,Cz,Cy
      integer, save:: i_ini

  

! 0. prepare
      if (i_ini .eq. 0) then
         ryx2 = (dy/dx)**2
         rzx2 = (dz/dx)**2

         d_ = (dx*dy*dz)**(1.0d0/3.0d0)

!         Cs = 0.1d0
         eps = 1.0d-6

         Ce_const = 0.132d0	

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

! the reason doing that is: only at bottom(me=0,iz=2), 
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
      call dealias1(ux,S11_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(uy,S22_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(uz,S33_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(vx,S12_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(vy,S13_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(vz,S23_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(wx,txz_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(wy,tyz_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(wz,tzz_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)

      do iz = 2, Nzb+1
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
      end do
      end do
      end do
  
      tmp1 = (2.0d0*d_/Ce_const)**2
      do iz = 2, Nzb+1
      do iy = 1, Ny2
      do ix = 1, Nx2
	  
         Gkk =  Gkk_m(ix,iy,iz)
         tmp3 = txx_m(ix,iy,iz)*S11_m(ix,iy,iz)+tyy_m(ix,iy,iz)*
     &          S22_m(ix,iy,iz)+tzz_m(ix,iy,iz)*S33_m(ix,iy,iz)+
     &          2.d0*(txy_m(ix,iy,iz)*S12_m(ix,iy,iz)+txz_m(ix,iy,iz)*
     &          S13_m(ix,iy,iz)+tyz_m(ix,iy,iz)*S23_m(ix,iy,iz))

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
            if (Cs.gt.0.36d0) Cs=0.0d0			
!            factor = Cs*(d_**2)*S+nu
		
            txx_m(ix,iy,iz) = -2.0d0*factor*S11_m(ix,iy,iz)
            tyy_m(ix,iy,iz) = -2.0d0*factor*S22_m(ix,iy,iz)
            tzz_m(ix,iy,iz) = -2.0d0*factor*S33_m(ix,iy,iz)
            txy_m(ix,iy,iz) = -2.0d0*factor*S12_m(ix,iy,iz) 			
		
            Cs2_m(ix,iy,iz)=(factor-nu)
      end do
      end do
      end do

      call dealias2(txx,txx_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)
      call dealias2(tyy,tyy_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)
      call dealias2(tzz,tzz_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)
      call dealias2(txy,txy_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)
      call dealias2(Cs2,Cs2_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)


      call update3(txx,tyy,tzz,me,nall)
      call update1(txy,me,nall)

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

      call dealias1(ux,S11_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(uy,S22_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(uz,S33_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(vx,S12_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(vy,S13_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(vz,S23_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(wx,txz_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(wy,tyz_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(wz,tzz_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)

      do iz = 2, Nzb+1
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
      end do
      end do
      end do

      tmp1 = (2.0d0*d_/Ce_const)**2
      do iz = 2, Nzb+1
      do iy = 1, Ny2
      do ix = 1, Nx2

         Gkk =  Gkk_m(ix,iy,iz)
         tmp3 = txx_m(ix,iy,iz)*S11_m(ix,iy,iz)+tyy_m(ix,iy,iz)*
     &          S22_m(ix,iy,iz)+tzz_m(ix,iy,iz)*S33_m(ix,iy,iz)+
     &          2.d0*(txy_m(ix,iy,iz)*S12_m(ix,iy,iz)+txz_m(ix,iy,iz)*
     &          S13_m(ix,iy,iz)+tyz_m(ix,iy,iz)*S23_m(ix,iy,iz)) 

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
            if (Cs.gt.0.36d0) Cs=0.0d0			
!            factor = Cs*(d_**2)*S+nu
			
            txz_m(ix,iy,iz) = -2.0d0*factor*S13_m(ix,iy,iz)
            tyz_m(ix,iy,iz) = -2.0d0*factor*S23_m(ix,iy,iz) 						
			
      end do
      end do
      end do

      call dealias2(txz,txz_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)
      call dealias2(tyz,tyz_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)

      if (me==0) then
! Stress Boundary Condition at bottom has been given, 
! thus, it needs to be saved and restored. Here is for resotration
         do iy=1,Ny
         do ix=1,Nx
            txz(ix,iy,2)=txzp(ix,iy)
            tyz(ix,iy,2)=tyzp(ix,iy)
         end do
         end do
         deallocate(txzp)
         deallocate(tyzp)
      end if

      call update1(txz,me,nall)
      call update1(tyz,me,nall)

!      write(*,*) 'countmgm1,2,3=',countmgm1,countmgm2,countmgm3  
      
      end
!-----------------------------------------------------------------------

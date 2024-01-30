      Subroutine calturbine(u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     turbine = 1 --> Actuator Disk model 
C     turbine = 2 --> An Improved Actuator Disk model 
C     turbine = 3 --> Actuator Line model 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'dimen.h'
      integer*4 t

      real*8, dimension(Nx,Ny,Nz2)::u,v,w,fx,fy,fz,fax,faz,AveWt

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Fx = 0.0d0
      Fy = 0.0d0
      Fz = 0.0d0

      if(nacelle_model.eq.1)then
      call CALNac   (u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
      end if

      if(tower_model.eq.1)then
c      call CALBar   (u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
      end if

      if(turbine_model.eq.1)then
      call CALADMNR (u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
      end if
      if(turbine_model.eq.2)then
      call CALADMR  (u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
      end if
      if(turbine_model.eq.3)then
      call CALALM   (u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
      end if

	  
      if(tower_model_vawt.eq.1)then
      call CALBarVAWT   (u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
      end if	  
      if(turbine_model_vawt.eq.1)then
      call CALVAWT   (u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
      end if	  
      if(turbine_model_vawt.eq.2)then
      call CALVAWT_ASSM   (u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
      end if	  
      if(turbine_model_vawt.eq.3)then
      call CALVAWT_FP   (u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
      end if
	  
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine CALBarVAWT(u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'dimen.h'
      integer*4 ii,i,j,k,t

      real*8, dimension(Nx,Ny,Nz2)::u,v,w,fx,fy,fz,fax,faz,AveWt
      real*8, dimension(Ny,Nz2):: CTangx
      real*8, dimension(num_turbine_vawt,Ny):: ker_y
      real*8, dimension(num_turbine_vawt,Nx,Ny,Nz2):: ker_u_bar_vawt  
      real*8, dimension(Ny,Nz2)::Nac_u

      real*8 a1,b1,at1,Tang(1:36),blade_jj,blade_kk,tmp,tmp1,tmp2
      real*8 mu_y,mu_z,sigma,hz1,hz2
      real*8 tmpy1,tmpz1,tmpy2,tmpz2
      real*8 ERRd1,ERRd2,tmpd1,tmpd2
      real*8 sum_ker_y,max_Nac_u,tmpb
      real*8 CTU,dA,CTFx,dr,CTV,CTM,CTFy
      real*8 sum_ker_u_bar(Nz2)	  

      real*8, dimension(num_turbine_vawt,Nx,Ny,Nz2):: CTy,CTx
      real*8, dimension(num_turbine_vawt,Nx,Ny,Nz2):: CTzu,CTzw,CTru,CTrw	
	  
      real*8, external :: interp2D,calang

      save ker_y,Nac_u,ker_u_bar_vawt
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (t.eq.1) then

        sum_ker_u_bar=0.0	
        ker_u_bar_vawt=0.0		
        do ii=1,num_turbine_vawt
        do k=2,nzb+1
        do j=1,ny
        do i=1,nx		
		
          CTy(ii,i,j,k)    = (wty_vawt(ii)-j)*dy+0.0d-50 
          CTx(ii,i,j,k)    = (wtx_vawt(ii)-i)*dx+0.0d-50		  
          CTru(ii,i,j,k)   = sqrt(CTy(ii,i,j,k)**2.+CTx(ii,i,j,k)**2.)

          if (CTru(ii,i,j,k).le.2.0*wtr_vawt(ii))then
          tmpb = CTru(ii,i,j,k)**2.0
          ker_u_bar_vawt(ii,i,j,k) = 1/(((1.*dx)**3)*(pi**1.5))
     +                           *dexp(-1.0*tmpb/((1.*dx)**2.0))
          sum_ker_u_bar(k) = sum_ker_u_bar(k)+ker_u_bar_vawt(ii,i,j,k)
          else
          ker_u_bar_vawt(ii,i,j,k) = 0.d0		  
          end if
		  
        end do
        end do
        end do
        end do		
	  
        do k=2,Nzb+1
        do j=1,Ny
        do i=1,Nx		  
          ker_u_bar_vawt(ii,i,j,k) = ker_u_bar_vawt(ii,i,j,k)
     +                               /sum_ker_u_bar(k)
        end do
        end do
        end do	

      end if

 
cc      if (me.le.15) then	  
	  
      do ii=1,num_turbine_vawt
      do k=2,nzb+1
 
          if ((me*nzb+k-1.5)*dz.gt.H1bar_vawt(ii) .and.
     +        (me*nzb+k-1.5)*dz.lt.H2bar_vawt(ii)) then	  
		  
          CTU = u(wtx_vawt(ii),wty_vawt(ii),k)
          CTV = v(wtx_vawt(ii),wty_vawt(ii),k)
          CTM = (CTU**2.0+CTV**2.0)**0.5		  
c         if(CTU.lt.0.d0)CTU=0.d0
          dA = Dbar_vawt(ii)*dz
          CTFx   = 0.5*(CTM*CTU)*CTB_vawt(ii)*dA/(dx*dy*dz)
          CTFy   = 0.5*(CTM*CTV)*CTB_vawt(ii)*dA/(dx*dy*dz)
		  
          do j=1,ny
          do i=1,nx
           if (ker_u_bar_vawt(ii,i,j,k).ge.1E-15)then
           Fx (i,j,k)=Fx (i,j,k) + CTFx  * ker_u_bar_vawt(ii,i,j,k)
           Fy (i,j,k)=Fy (i,j,k) + CTFy  * ker_u_bar_vawt(ii,i,j,k)
           end if
          end do
          end do

          end if		 

      end do
      end do
		
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine CALNac(u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'dimen.h'
      integer*4 ii,i,j,k,t,blade,N,me_k

      real*8, dimension(Nx,Ny,Nz2)::u,v,w,fx,fy,fz,fax,faz,AveWt
      real*8, dimension(num_turbine,Ny,Nz2):: CTy,CTzu,CTzw,CTru,CTrw
      real*8, dimension(num_turbine,Ny,Nz2):: CTangx,CTangz
      real*8, dimension(num_turbine,Nx):: ker_x
      real*8 tmp1,tmp2,tmp3

      real*8, dimension(num_turbine,N_ang,N_blade)::
     +uframe,vframe,wframe,meframe,u1frame,v1frame,w1frame,me1frame

      real*8, dimension(num_turbine,N_ang,N_blade,Ny,Nz2):: ker_u,ker_w

      real*8 a1,b1,at1,Tang(N_ang),blade_jj,blade_kk
      real*8 mu_y,mu_z,sigma,dr,omega,omega1,r,CTphi
      real*8 tmpy1,tmpz1,tmpy2,tmpz2,tmp
      real*8 ERRd1,ERRd2,tmpd1,tmpd2
      real*8 sum_ker_u,sum_ker_w,max_ker_u,max_ker_w
      real*8 CTU,CTU1,CTV,CTW,chordl,CTsoli
      real*8 CTV1,CTV2,CTUrel,CTFx,CTFt
      real*8 CTCL,CTCD,dang,Re,dA,RA,dL

      real*8, external :: chord,calang
      real*8, external :: itp2D,calU

      save CTy,CTzu,CTzw,CTangx,CTangz,ker_x,ker_u,ker_w
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      dang = 2.d0*pi/dble(N_ang)

      if (t.eq.1) then
        do ii=1,num_turbine
        do k=2,nzb+1
        do j=1,ny
          CTy(ii,j,k)    = (wty(ii)-j)*dy+0.0d-50 

          CTzu(ii,j,k)   = (me*nzb+k-1.5  )*dz-Zhub(ii)
          CTzw(ii,j,k)   = CTzu(ii,j,k)-0.5d0*dz

          CTru(ii,j,k)   = sqrt(CTy(ii,j,k)**2.+CTzu(ii,j,k)**2.)
          CTrw(ii,j,k)   = sqrt(CTy(ii,j,k)**2.+CTzw(ii,j,k)**2.)

          CTangx(ii,j,k) = calang(CTy(ii,j,k),CTzu(ii,j,k))
          CTangz(ii,j,k) = calang(CTy(ii,j,k),CTzw(ii,j,k))
        end do
        end do
        end do

        sigma = 1.d0
        tmp1  = 1/dsqrt(2.d0*pi*sigma*sigma)
        tmp2  = -1/(2.d0*sigma*sigma)
        do ii=1,num_turbine
        do i=1,Nx
          ker_x(ii,i) = tmp1*dexp(tmp2*(i-wtx(ii))**2.)
        end do
        tmp1 = sum(ker_x(ii,1:Nx))
        do i=1,Nx
          ker_x(ii,i) = ker_x(ii,i)/tmp1
        end do
        end do
      end if

      if (t.eq.1) then
      do i=1,N_ang	
      Tang(i) = (dble(i)-1)*dang
      end do

      do ii=1,num_turbine
        dr = wtr(ii)/dble(N_blade)
        do N=1,N_blade_s-1
        r = dble(N-0.5)*dr
        do blade=1,N_ang

          sum_ker_u  = 0.0
          sum_ker_w  = 0.0
          mu_y       = r*dcos(Tang(blade))  
          mu_z       = r*dsin(Tang(blade)) 
          sigma      = dsqrt((r*dang)**2.+dr**2.)
c          sigma      = dsqrt(dy*dy+dz*dz)*0.75d0

          tmp1       = 1/dsqrt(2.d0*pi*sigma*sigma)
          tmp2       = -1/(2.d0*sigma*sigma)

          do k=2,nzb+1
          do j=1,ny
          if (CTru(ii,j,k).le.2.0*wtr(ii))then
          tmp3 = (CTy(ii,j,k)-mu_y)**2.+(CTzu(ii,j,k)-mu_z)**2.
          ker_u(ii,blade,N,j,k) = tmp1*dexp(tmp2*tmp3)
          sum_ker_u = sum_ker_u+ker_u(ii,blade,N,j,k)
          else
          ker_u(ii,blade,N,j,k) = 0.d0
          end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
          if (CTrw(ii,j,k).le.2.0*wtr(ii))then
          tmp3 = (CTy(ii,j,k)-mu_y)**2.+(CTzw(ii,j,k)-mu_z)**2.
          ker_w(ii,blade,N,j,k) = tmp1*dexp(tmp2*tmp3)
          sum_ker_w = sum_ker_w+ker_w(ii,blade,N,j,k)
          else
          ker_w(ii,blade,N,j,k) = 0.d0
          end if
          end do
          end do

          call calsum(sum_ker_u,max_ker_u,me,nall)
          call calsum(sum_ker_w,max_ker_w,me,nall)

          do k=2,Nzb+1
          do j=1,Ny
          ker_u(ii,blade,N,j,k) = ker_u(ii,blade,N,j,k)/sum_ker_u
          ker_w(ii,blade,N,j,k) = ker_w(ii,blade,N,j,k)/sum_ker_w
          end do
          end do

        end do
        end do
      end do
      end if


      meframe=0.d0
      uframe=0.d0
      vframe=0.d0
      wframe=0.d0
      do ii=1,num_turbine
      dr = wtr(ii)/dble(N_blade)
      do N=1,N_blade_s-1
      r = dble(N-0.5)*dr
      do blade=1,N_ang
      blade_kk = dble(Zhub(ii)/dz+0.5)+r*dsin(Tang(blade))/dz-1.d0
      me_k=blade_kk/Nzb

      if (me.eq.me_k) then
      blade_jj = dble(wty(ii))-r*dcos(Tang(blade))/dy
      blade_kk = blade_kk-Nzb*me+2.d0
      
      meframe(ii,blade,N)=1.d0
      i = 2.d0*wtR(ii)/dx
      uframe(ii,blade,N)=itp2D(u,wtx(ii)-i,blade_jj,blade_kk,me,nall)
      i = 1
      vframe(ii,blade,N)=itp2D(v,wtx(ii)-i,blade_jj,blade_kk,me,nall)
      i = 1
      wframe(ii,blade,N)=itp2D(w,wtx(ii)-i,blade_jj,blade_kk,me,nall)
      end if
      end do
      end do
      end do

      IF (me>0) then
         call MPI_SEND(meframe(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(uframe(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(vframe(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(wframe(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
      ELSE
         do i=1,nprocs-1
         call MPI_RECV(me1frame(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         call MPI_RECV(u1frame(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         call MPI_RECV(v1frame(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         call MPI_RECV(w1frame(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )

         tmp1=sum(me1frame)
         if(tmp1.gt.0.5)then
           do ii=1,num_turbine
           do N=1,N_blade_s
           do blade=1,N_ang
           if(me1frame(ii,blade,N).gt.0.5d0)then
           uframe(ii,blade,N)=u1frame(ii,blade,N)
           vframe(ii,blade,N)=v1frame(ii,blade,N)
           wframe(ii,blade,N)=w1frame(ii,blade,N)
           end if
           end do
           end do
           end do
         end if
         end do

         do i=1,nprocs-1
         call MPI_SEND(uframe(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )
         call MPI_SEND(vframe(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )
         call MPI_SEND(wframe(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )
         end do
      END IF
      IF (me>0) then
          call MPI_RECV(uframe(1,1,1),num_turbine*N_ang*N_blade_s,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
          call MPI_RECV(vframe(1,1,1),num_turbine*N_ang*N_blade_s,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
          call MPI_RECV(wframe(1,1,1),num_turbine*N_ang*N_blade_s,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
      END IF


      do ii=1,num_turbine
      dr = wtr(ii)/dble(N_blade)
      do N=1,N_blade_s-1
      r = dble(N-0.5)*dr
      do blade=1,N_ang

      CTU = uframe(ii,blade,N)
cc      CTU = sum(uframe(ii,:,1))/dble(N_ang)
      tmp1=1.d0

      dA      = 0.5d0*(dble(N)**2-dble(N-1)**2)*(dr**2.)*dang

      CTFx    = 0.5*(CTU**2)*CTN(ii)*dA/(dx*dy*dz)*tmp1

      i=wtx(ii)
      do k=2,nzb+1
      do j=1,ny
         if (ker_u(ii,blade,N,j,k).ge.1E-15)then
         Fx (i,j,k)=Fx (i,j,k) + CTFx  * ker_u(ii,blade,N,j,k)
         end if
      end do
      end do

      end do
      end do
      end do

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine CALADMNR(u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'dimen.h'
      integer*4 ii,i,j,k,t,blade,N,me_k

      real*8, dimension(Nx,Ny,Nz2)::u,v,w,fx,fy,fz,fax,faz,AveWt
      real*8, dimension(num_turbine,Ny,Nz2):: CTy,CTzu,CTzw,CTru,CTrw
      real*8, dimension(num_turbine,Ny,Nz2):: CTangx,CTangz
      real*8, dimension(num_turbine,Nx):: ker_x
      real*8 tmp1,tmp2,tmp3

      real*8, dimension(num_turbine,N_ang,N_blade)::
     +uframe,vframe,wframe,meframe,u1frame,v1frame,w1frame,me1frame

      real*8, dimension(num_turbine,N_ang,N_blade,Ny,Nz2):: ker_u,ker_w

      real*8 a1,b1,at1,Tang(N_ang),blade_jj,blade_kk
      real*8 mu_y,mu_z,sigma,dr,omega,omega1,r,CTphi
      real*8 tmpy1,tmpz1,tmpy2,tmpz2,tmp
      real*8 ERRd1,ERRd2,tmpd1,tmpd2
      real*8 sum_ker_u,sum_ker_w,max_ker_u,max_ker_w
      real*8 CTU,CTU1,CTV,CTW,chordl,CTsoli
      real*8 CTV1,CTV2,CTUrel,CTFx,CTFt
      real*8 CTCL,CTCD,dang,Re,dA,RA,dL

      real*8, external :: chord,calang
      real*8, external :: itp2D,calU

      save CTy,CTzu,CTzw,CTangx,CTangz,ker_x,ker_u,ker_w
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      dang = 2.d0*pi/dble(N_ang)

      if (t.eq.1) then
        do ii=1,num_turbine
        do k=2,nzb+1
        do j=1,ny
          CTy(ii,j,k)    = (wty(ii)-j)*dy+0.0d-50 

          CTzu(ii,j,k)   = (me*nzb+k-1.5  )*dz-Zhub(ii)
          CTzw(ii,j,k)   = CTzu(ii,j,k)-0.5d0*dz

          CTru(ii,j,k)   = sqrt(CTy(ii,j,k)**2.+CTzu(ii,j,k)**2.)
          CTrw(ii,j,k)   = sqrt(CTy(ii,j,k)**2.+CTzw(ii,j,k)**2.)

          CTangx(ii,j,k) = calang(CTy(ii,j,k),CTzu(ii,j,k))
          CTangz(ii,j,k) = calang(CTy(ii,j,k),CTzw(ii,j,k))
        end do
        end do
        end do

        sigma = 1.d0
        tmp1  = 1/dsqrt(2.d0*pi*sigma*sigma)
        tmp2  = -1/(2.d0*sigma*sigma)
        do ii=1,num_turbine
        do i=1,Nx
          ker_x(ii,i) = tmp1*dexp(tmp2*(i-wtx(ii))**2.)
        end do
        tmp1 = sum(ker_x(ii,1:Nx))
        do i=1,Nx
          ker_x(ii,i) = ker_x(ii,i)/tmp1
        end do
        end do
      end if

      if (t.eq.1) then
      do i=1,N_ang	
      Tang(i) = (dble(i)-1)*dang
      end do

      do ii=1,num_turbine
        dr = wtr(ii)/dble(N_blade)
        do N=1,N_blade
        r = dble(N-0.5)*dr
        do blade=1,N_ang

          sum_ker_u  = 0.0
          sum_ker_w  = 0.0
          mu_y       = r*dcos(Tang(blade))  
          mu_z       = r*dsin(Tang(blade)) 
          sigma      = dsqrt((r*dang)**2.+dr**2.)

          tmp1       = 1/dsqrt(2.d0*pi*sigma*sigma)
          tmp2       = -1/(2.d0*sigma*sigma)

          do k=2,nzb+1
          do j=1,ny
          if (CTru(ii,j,k).le.2.0*wtr(ii))then
          tmp3 = (CTy(ii,j,k)-mu_y)**2.+(CTzu(ii,j,k)-mu_z)**2.
          ker_u(ii,blade,N,j,k) = tmp1*dexp(tmp2*tmp3)
          sum_ker_u = sum_ker_u+ker_u(ii,blade,N,j,k)
          else
          ker_u(ii,blade,N,j,k) = 0.d0
          end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
          if (CTrw(ii,j,k).le.2.0*wtr(ii))then
          tmp3 = (CTy(ii,j,k)-mu_y)**2.+(CTzw(ii,j,k)-mu_z)**2.
          ker_w(ii,blade,N,j,k) = tmp1*dexp(tmp2*tmp3)
          sum_ker_w = sum_ker_w+ker_w(ii,blade,N,j,k)
          else
          ker_w(ii,blade,N,j,k) = 0.d0
          end if
          end do
          end do

          call calsum(sum_ker_u,max_ker_u,me,nall)
          call calsum(sum_ker_w,max_ker_w,me,nall)

          do k=2,Nzb+1
          do j=1,Ny
          ker_u(ii,blade,N,j,k) = ker_u(ii,blade,N,j,k)/sum_ker_u
          ker_w(ii,blade,N,j,k) = ker_w(ii,blade,N,j,k)/sum_ker_w
          end do
          end do

        end do
        end do
      end do
      end if


      meframe=0.d0
      uframe=0.d0
      vframe=0.d0
      wframe=0.d0
      do ii=1,num_turbine
      dr = wtr(ii)/dble(N_blade)
      do N=1,N_blade
      r = dble(N-0.5)*dr
      do blade=1,N_ang
      blade_kk = dble(Zhub(ii)/dz+0.5)+r*dsin(Tang(blade))/dz-1.d0
      me_k=blade_kk/Nzb

      if (me.eq.me_k) then
      blade_jj = dble(wty(ii))-r*dcos(Tang(blade))/dy
      blade_kk = blade_kk-Nzb*me+2.d0
      
      meframe(ii,blade,N)=1.d0
      i = 2.d0*wtR(ii)/dx
      uframe(ii,blade,N)=itp2D(u,wtx(ii)-i,blade_jj,blade_kk,me,nall)
      i = 1
      vframe(ii,blade,N)=itp2D(v,wtx(ii)-i,blade_jj,blade_kk,me,nall)
      i = 1
      wframe(ii,blade,N)=itp2D(w,wtx(ii)-i,blade_jj,blade_kk,me,nall)
      end if
      end do
      end do
      end do

      IF (me>0) then
         call MPI_SEND(meframe(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(uframe(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(vframe(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(wframe(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
      ELSE
         do i=1,nprocs-1
         call MPI_RECV(me1frame(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         call MPI_RECV(u1frame(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         call MPI_RECV(v1frame(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         call MPI_RECV(w1frame(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )

         tmp1=sum(me1frame)
         if(tmp1.gt.0.5)then
           do ii=1,num_turbine
           do N=1,N_blade_s
           do blade=1,N_ang
           if(me1frame(ii,blade,N).gt.0.5d0)then
           uframe(ii,blade,N)=u1frame(ii,blade,N)
           vframe(ii,blade,N)=v1frame(ii,blade,N)
           wframe(ii,blade,N)=w1frame(ii,blade,N)
           end if
           end do
           end do
           end do
         end if
         end do

         do i=1,nprocs-1
         call MPI_SEND(uframe(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )
         call MPI_SEND(vframe(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )
         call MPI_SEND(wframe(1,1,1),num_turbine*N_ang*N_blade_s,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )
         end do
      END IF
      IF (me>0) then
          call MPI_RECV(uframe(1,1,1),num_turbine*N_ang*N_blade_s,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
          call MPI_RECV(vframe(1,1,1),num_turbine*N_ang*N_blade_s,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
          call MPI_RECV(wframe(1,1,1),num_turbine*N_ang*N_blade_s,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
      END IF


      do ii=1,num_turbine
      dr = wtr(ii)/dble(N_blade)
      do N=N_blade_s,N_blade
      r = dble(N-0.5)*dr
      do blade=1,N_ang

      CTU = sum(uframe(ii,:,1))/dble(N_ang)

      tmp1=1.d0

c      if(CTU*u_star.lt.8.)then
c      tmp1=(CTU*u_star-2.)/6.
c      end if
c      if(CTU*u_star.lt.2.)then
c      tmp1=0.0
c      end if

      dA      = 0.5d0*(dble(N)**2-dble(N-1)**2)*(dr**2.)*dang

      CTFx    = 0.5*(CTU**2)*CTT(ii)*dA/(dx*dy*dz)

      i=wtx(ii)
      do k=2,nzb+1
      do j=1,ny
         if (ker_u(ii,blade,N,j,k).ge.1E-15)then
         Fx (i,j,k)=Fx (i,j,k) + CTFx  * ker_u(ii,blade,N,j,k)
         end if
      end do
      end do

      end do
      end do
      end do

      return

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(me.eq.0)then
      do ii=1,num_turbine
      CTU = sum(uframe(ii,:,1))/dble(N_ang)*u_star
      write(570)CTU
      end do
      call flush(570)
      end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do ii=1,num_turbine
      do k=2,nzb+1
      do j=1,ny
      	tmp = Fx(wtx(ii),j,k)
        do i=1,nx
      	Fx(i,j,k)=tmp*ker_x(ii,i)
        end do

      	tmp = Fy(wtx(ii),j,k)
        do i=1,nx
      	Fy(i,j,k)=tmp*ker_x(ii,i)
        end do

      	tmp = Fz(wtx(ii),j,k)
        do i=1,nx
      	Fz(i,j,k)=tmp*ker_x(ii,i)
        end do

      end do
      end do
      end do

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine CALADMR(u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'dimen.h'
      integer*4 ii,i,j,k,t,blade,N,me_k,iter_omega,iter_pitch

      real*8, dimension(Nx,Ny,Nz2)::u,v,w,fx,fy,fz,fax,faz,AveWt
      real*8, dimension(num_turbine,Ny,Nz2):: CTy,CTzu,CTzw,CTru,CTrw
      real*8, dimension(num_turbine,Ny,Nz2):: CTangx,CTangz
      real*8, dimension(num_turbine,Nx):: ker_x
      real*8 tmp1,tmp2,tmp3

      real*8, dimension(num_turbine,N_ang,N_blade)::
     +uframe,vframe,wframe,meframe,u1frame,v1frame,w1frame,me1frame,
     +uframed,u1framed
      
      real*8, dimension(num_turbine,N_ang,N_blade,Ny,Nz2):: ker_u,ker_w

      real*8, dimension(num_turbine):: SUMCT,SUMa1,SUMdA,SUMCTU,SUMFx,
     +SUMb1,SUMCTUd,SUMQ,Power
      
      real*8 a1,b1,at1,Tang(N_ang),blade_jj,blade_kk
      real*8 mu_y,mu_z,sigma,dr,omega,omega1,r,CTphi
      real*8 tmpy1,tmpz1,tmpy2,tmpz2,tmp
      real*8 ERRd1,ERRd2,tmpd1,tmpd2
      real*8 sum_ker_u,sum_ker_w,max_ker_u,max_ker_w
      real*8 CTU,CTU1,CTV,CTW,chordl,CTsoli,CTUd
      real*8 CTV1,CTV2,CTUrel,CTFx,CTFt
      real*8 CTCL,CTCD,dang,Re,dA,RA,dL
      real*8 umean(num_turbine),wtomegaL(num_turbine),
     +       umeand(num_turbine),wtomegaLinf(num_turbine),
     +       Pitch_t(num_turbine),pitch_max(num_turbine),
     +       Power_max(num_turbine),wtomega_max(num_turbine)	 

      real*8, external :: chord,calang
      real*8, external :: itp2D,calU

      real*8 a1mean,b1mean,CTmean,CTUmean,CTUdmean,Fxmean,Powmean,
     +       omegamean,u_hub,v_hub

      save CTy,CTzu,CTzw,CTangx,CTangz,ker_x,ker_u,ker_w,dang
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      if (t.eq.1) then
	  
        do ii=1,num_turbine
        do k=2,nzb+1
        do j=1,ny
          CTy(ii,j,k)    = (wty(ii)-j)*dy+0.0d-50 

          CTzu(ii,j,k)   = (me*nzb+k-1.5  )*dz-Zhub(ii)
          CTzw(ii,j,k)   = CTzu(ii,j,k)-0.5d0*dz

          CTru(ii,j,k)   = sqrt(CTy(ii,j,k)**2.+CTzu(ii,j,k)**2.)
          CTrw(ii,j,k)   = sqrt(CTy(ii,j,k)**2.+CTzw(ii,j,k)**2.)

          CTangx(ii,j,k) = calang(CTy(ii,j,k),CTzu(ii,j,k))
          CTangz(ii,j,k) = calang(CTy(ii,j,k),CTzw(ii,j,k))
        end do
        end do
        end do

        sigma = 1.d0
        tmp1  = 1/dsqrt(2.d0*pi*sigma*sigma)
        tmp2  = -1/(2.d0*sigma*sigma)
        do ii=1,num_turbine
        do i=1,Nx
          ker_x(ii,i) = tmp1*dexp(tmp2*(i-wtx(ii))**2.)
        end do
        tmp1 = sum(ker_x(ii,1:Nx))
        do i=1,Nx
          ker_x(ii,i) = ker_x(ii,i)/tmp1
        end do
        end do
      end if

      if (t.eq.1) then
      dang = 2.d0*pi/dble(N_ang)
      do i=1,N_ang		
      Tang(i) = (dble(i)-1)*dang
      end do

      do ii=1,num_turbine
        dr = wtr(ii)/dble(N_blade)
        do N=1,N_blade
        r = dble(N-0.5)*dr
        do blade=1,N_ang

          sum_ker_u  = 0.0
          sum_ker_w  = 0.0
          mu_y       = r*dcos(Tang(blade))  
          mu_z       = r*dsin(Tang(blade)) 
          sigma      = dsqrt((r*dang)**2.+dr**2.)
c          sigma      = dsqrt(dy*dy+dz*dz)*0.75d0

          tmp1       = 1/dsqrt(2.d0*pi*sigma*sigma)
          tmp2       = -1/(2.d0*sigma*sigma)

          do k=2,nzb+1
          do j=1,ny
          if (CTru(ii,j,k).le.4.0*wtr(ii))then
          tmp3 = (CTy(ii,j,k)-mu_y)**2.+(CTzu(ii,j,k)-mu_z)**2.
          ker_u(ii,blade,N,j,k) = tmp1*dexp(tmp2*tmp3)
          sum_ker_u = sum_ker_u+ker_u(ii,blade,N,j,k)
          else
          ker_u(ii,blade,N,j,k) = 0.d0
          end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
          if (CTrw(ii,j,k).le.4.0*wtr(ii))then
          tmp3 = (CTy(ii,j,k)-mu_y)**2.+(CTzw(ii,j,k)-mu_z)**2.
          ker_w(ii,blade,N,j,k) = tmp1*dexp(tmp2*tmp3)
          sum_ker_w = sum_ker_w+ker_w(ii,blade,N,j,k)
          else
          ker_w(ii,blade,N,j,k) = 0.d0
          end if
          end do
          end do

          call calsum(sum_ker_u,max_ker_u,me,nall)
          call calsum(sum_ker_w,max_ker_w,me,nall)

          do k=2,Nzb+1
          do j=1,Ny
          ker_u(ii,blade,N,j,k) = ker_u(ii,blade,N,j,k)/sum_ker_u
          ker_w(ii,blade,N,j,k) = ker_w(ii,blade,N,j,k)/sum_ker_w
          end do
          end do

        end do
        end do
      end do
      end if


      meframe=0.d0
      uframe=0.d0
      uframed=0.d0
      vframe=0.d0
      wframe=0.d0

      do ii=1,num_turbine
      dr = wtr(ii)/dble(N_blade)
      do N=1,N_blade
      r = dble(N-0.5)*dr
      do blade=1,N_ang

      blade_kk = dble(Zhub(ii)/dz+0.5)+r*dsin(Tang(blade))/dz-1.d0
      me_k=blade_kk/Nzb

      if (me.eq.me_k) then
      blade_jj = dble(wty(ii))-r*dcos(Tang(blade))/dy
      blade_kk = blade_kk-Nzb*me+2.d0
      
      meframe(ii,blade,N)=1.d0
      i = 0*wtR(ii)/dx
      uframe(ii,blade,N)=itp2D(u,wtx(ii)-i,blade_jj,blade_kk,me,nall)
      i = 2.0*wtR(ii)/dx
      uframed(ii,blade,N)=itp2D(u,wtx(ii)-i,blade_jj,blade_kk,me,nall)	  
      i = 0*wtR(ii)/dx
      vframe(ii,blade,N)=itp2D(v,wtx(ii)-i,blade_jj,blade_kk,me,nall)
      i = 0*wtR(ii)/dx
      wframe(ii,blade,N)=itp2D(w,wtx(ii)-i,blade_jj,blade_kk,me,nall)
      end if
      end do
      end do
      end do

      IF (me>0) then
         call MPI_SEND(meframe(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(uframe(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(uframed(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	 
         call MPI_SEND(vframe(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(wframe(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
      ELSE
         do i=1,nprocs-1
         call MPI_RECV(me1frame(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         call MPI_RECV(u1frame(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         call MPI_RECV(u1framed(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	 
         call MPI_RECV(v1frame(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         call MPI_RECV(w1frame(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )

         tmp1=sum(me1frame)
         if(tmp1.gt.0.5)then
           do ii=1,num_turbine
           do N=1,N_blade
           do blade=1,N_ang
           if(me1frame(ii,blade,N).gt.0.5d0)then
           uframe(ii,blade,N)=u1frame(ii,blade,N)
           uframed(ii,blade,N)=u1framed(ii,blade,N)
           vframe(ii,blade,N)=v1frame(ii,blade,N)
           wframe(ii,blade,N)=w1frame(ii,blade,N)
           end if
           end do
           end do
           end do
         end if
         end do

         do i=1,nprocs-1
         call MPI_SEND(uframe(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )
         call MPI_SEND(uframed(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )	 
         call MPI_SEND(vframe(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )
         call MPI_SEND(wframe(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )
         end do
      END IF
      IF (me>0) then
          call MPI_RECV(uframe(1,1,1),num_turbine*N_ang*N_blade,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
          call MPI_RECV(uframed(1,1,1),num_turbine*N_ang*N_blade,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )	 
          call MPI_RECV(vframe(1,1,1),num_turbine*N_ang*N_blade,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
          call MPI_RECV(wframe(1,1,1),num_turbine*N_ang*N_blade,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
      END IF

      do ii=1,num_turbine

!      dir(ii)=1.0d0 
      dr = wtr(ii)/dble(N_blade)

      umean(ii)=sum(uframe(ii,:,:))/(dble(N_ang)*N_blade)
      umeand(ii)=sum(uframed(ii,:,:))/(dble(N_ang)*N_blade)

!      if (me.eq.0) write(*,*) 'ADMR yes'	  
	  
!      Power_max(ii)=0d0	 	  
!      wtomega(ii)=12.0d0
      do iter_omega=1,5
!      wtomega(ii)=wtomega(ii)+0.5d0	
      Pitch_t(ii)=0.0d0
!      do iter_pitch=1,21
!      Pitch_t(ii)=Pitch_t(ii)+0.5d0
      SUMQ(ii)=0.0d0
      do N=N_blade_s,N_blade
      r = dble(N-0.5)*dr 
      chordl  = chord(r*z_i)/z_i	 

!      if (me.eq.0 .and. N.eq.2) then
!           write(*,*) 'chord',chordl*z_i
!      endif		   
      CTsoli  = 3.0*chordl/(2.d0*pi*r)
      do blade=1,N_ang
      CTU = uframe(ii,blade,N)
      CTUd = uframed(ii,blade,N)
      CTV = vframe(ii,blade,N)
      CTW = wframe(ii,blade,N)
      omega1  = +0*CTW*dcos(Tang(blade))/r
     +          +0*CTV*dsin(Tang(blade))/r
     +          +wtomega(ii)*(2.d0*pi/60)*(z_i/u_star)
      call sol_CL_CD_1(ii,Pitch_t(ii),omega1,chordl,r,CTU,CTsoli,
     +     a1,b1,CTCL,CTCD,at1,CTphi,Re,me,nall)
!      CTCL=0.98*CTCL
!      CTCD=0.98*CTCD	 
      CTUrel  = CTU/dsin(CTphi)
!      call sol_CL_CD(ii,Pitch_t(ii),omega1,chordl,r,CTU,CTsoli,
!     +     a1,b1,CTCL,CTCD,at1,CTphi,Re,me,nall)
!      CTUrel  = CTU*(1-a1)/dsin(CTphi)

      dA      = 0.5d0*(dble(N)**2-dble(N-1)**2)*(dr**2.)*dang
      CTFx    = 0.5*(CTUrel**2)*CTsoli*dA/(dx*dy*dz)*
     +          (+CTCL*dcos(CTphi)+CTCD*dsin(CTphi))
      CTFt    = -dir(ii)*0.5*(CTUrel**2)*CTsoli*dA/(dx*dy*dz)*
     +          (+CTCL*dsin(CTphi)-CTCD*dcos(CTphi))
      SUMQ(ii)=SUMQ(ii)-dir(ii)*CTFt*(dx*dy*dz)*r*(z_i**3)
      end do
      end do  
      Power(ii)=SUMQ(ii)*wtomega(ii)*(2.d0*pi/60)	  
!      if (Power(ii).gt.Power_max(ii)) then
!      Power_max(ii)=Power(ii)
!          wtomega_max(ii)=wtomega(ii)
!          Pitch_max(ii)=Pitch_t(ii) 		  
!      endif	  
!      if (me.eq.0) print*, wtomega(ii),Pitch_t(ii),Power(ii)/10**6	  
!      wtomega(ii)=wtomega(ii)+0.5*(
!     +            +5.648D-11*(1.0*SUMQ(ii)/(10**3))**4
!     +            -1.365D-7*(1.0*SUMQ(ii)/(10**3))**3
!     +            +0.0001038*(1.0*SUMQ(ii)/(10**3))**2
!     +            -0.01895*(1.0*SUMQ(ii)/(10**3))
!     +            +13.23-wtomega(ii))	
      wtomega(ii)=wtomega(ii)+0.5*(
     +            +7.655D-11*(1.0*SUMQ(ii)/(10**3))**4
     +            -1.584D-7*(1.0*SUMQ(ii)/(10**3))**3
     +            +0.0001041*(1.0*SUMQ(ii)/(10**3))**2
     +            -0.01292*(1.0*SUMQ(ii)/(10**3))
     +            +12.37-wtomega(ii))	 
      if (wtomega(ii).gt.18.5) wtomega(ii)=18.5
      if (wtomega(ii).lt.11.5) wtomega(ii)=11.5
!      enddo
      enddo
!      wtomega(ii)=wtomega_max(ii) 	  
!      Pitch_t(ii)=Pitch_max(ii)
!      if (me.eq.0) print*, umean(ii),umeand(ii),wtomega(ii),Pitch_t(ii)
!      wtomega(ii)=11.9d0 	  
!      Pitch_t(ii)=0.0

!     wtomega(ii)=1.2d0*umean(ii)*TSR(ii)/(wtr(ii)*z_i)
!     +                    *60.d0/(2.0*pi)
  
      SUMdA(ii)=0.0
      SUMa1(ii)=0.0
      SUMb1(ii)=0.0	  
      SUMCT(ii)=0.0
      SUMCTU(ii)=0.0	  
      SUMCTUd(ii)=0.0
      SUMFx(ii)=0.0	
      SUMQ(ii)=0.0
	  
      do N=N_blade_s,N_blade
      r = dble(N-0.5)*dr 
      chordl  = chord(r*z_i)/z_i	  
      CTsoli  = 3.0*chordl/(2.d0*pi*r)
      do blade=1,N_ang
      CTU = uframe(ii,blade,N)
      CTUd = uframed(ii,blade,N)
      CTV = vframe(ii,blade,N)
      CTW = wframe(ii,blade,N)
      omega1  = +0*CTW*dcos(Tang(blade))/r
     +          +0*CTV*dsin(Tang(blade))/r
     +          +wtomega(ii)*(2.d0*pi/60)*(z_i/u_star)
cc      print*, 	  CTU,omega1,wtomegaL(ii)	 
cc      call sol_CL_CD(ii,Pitch_t(ii),omega1,chordl,r,CTU,CTsoli,
cc     +     a1,b1,CTCL,CTCD,at1,CTphi,Re,me,nall)
cc      CTUrel  = CTU*(1-a1)/dsin(CTphi)
cc      print*, 	  CTU	 

      call sol_CL_CD_1(ii,Pitch_t(ii),omega1,chordl,r,CTU,CTsoli,
     +     a1,b1,CTCL,CTCD,at1,CTphi,Re,me,nall)
!      CTCL=0.98*CTCL
!      CTCD=0.98*CTCD	  
      CTUrel  = 1.0*CTU/dsin(CTphi)	 
c      if (ii.eq.1 .and. N.eq.N_blade .and. blade.eq.1) then
c            print*, a1,b1
c      endif
cc      CTV1    =  CTU
cc      CTV2    = -CTW*dcos(Tang(blade))
cc     +            -CTV*dsin(Tang(blade))-omega1*r

	  
      dA      = 0.5d0*(dble(N)**2-dble(N-1)**2)*(dr**2.)*dang
      CTFx    = 0.5*(CTUrel**2)*CTsoli*dA/(dx*dy*dz)*
     +          (+CTCL*dcos(CTphi)+CTCD*dsin(CTphi))
      CTFt    = -dir(ii)*0.5*(CTUrel**2)*CTsoli*dA/(dx*dy*dz)*
     +          (+CTCL*dsin(CTphi)-CTCD*dcos(CTphi))

      SUMdA(ii) =SUMdA(ii)+dA*(z_i**2)
      SUMa1(ii) =SUMa1(ii)+a1*dA/(pi*wtR(ii)**2.d0)
      SUMb1(ii) =SUMb1(ii)+b1*dA/(pi*wtR(ii)**2.d0)
      SUMCT(ii) =SUMCT(ii)+CTFx/(0.5d0*((CTU/u_star)**2.)*	      
     +           dA/(dx*dy*dz))*dA/(pi*wtR(ii)**2.d0)
      SUMCTU(ii)=SUMCTU(ii)+CTU*dA/(pi*wtR(ii)**2.d0)
      SUMCTUd(ii)=SUMCTUd(ii)+CTUd*dA/(pi*wtR(ii)**2.d0)	  
      SUMFx(ii)=SUMFx(ii)+CTFx*(dx*dy*dz)
      SUMQ(ii)=SUMQ(ii)-dir(ii)*CTFt*(dx*dy*dz)*r*(z_i**3)

      i=wtx(ii)
      do k=2,nzb+1
      do j=1,ny
         if (ker_u(ii,blade,N,j,k).ge.1E-15)then
         Fx (i,j,k)=Fx (i,j,k) + CTFx  * ker_u(ii,blade,N,j,k)

         Fax(i,j,k)=Fax(i,j,k) + CTFt  * ker_u(ii,blade,N,j,k)
         Fy (i,j,k)=Fy (i,j,k) + CTFt  * ker_u(ii,blade,N,j,k)
     +             *dsin(CTangx(ii,j,k))
         end if

         if (ker_w(ii,blade,N,j,k).ge.1E-15)then
         Faz(i,j,k)=Faz(i,j,k) + CTFt  * ker_w(ii,blade,N,j,k)
         Fz (i,j,k)=Fz (i,j,k) + CTFt  * ker_w(ii,blade,N,j,k)
     +             *dcos(CTangz(ii,j,k))
         end if
      end do
      end do

      end do
      end do
      Power(ii)=SUMQ(ii)*wtomega(ii)*(2.d0*pi/60)	  
      end do
	  
      a1mean=sum(SUMa1)/(max(1,size(SUMa1)))
      b1mean=sum(SUMb1)/(max(1,size(SUMb1)))
      CTmean=sum(SUMCT)/(max(1,size(SUMCT)))
      CTUmean=sum(SUMCTU)/(max(1,size(SUMCTU)))
      CTUdmean=sum(SUMCTUd)/(max(1,size(SUMCTUd)))
      Fxmean=sum(SUMFx)/(max(1,size(SUMFx))) 
      Powmean=sum(Power)/(max(1,size(Power)))
      Omegamean=sum(wtomega)/(max(1,size(wtomega)))	  
	  
      if(me.eq.0)then
      if(t.eq.1)then
      Open (unit=1095,file='output/check_a_CT.out',position='append')
      end if
	  
      if (mod(t,c_count)==0) then

      write(1524,5544) Power/(10**6),SUMCTU
      write(1523,5544) wtomega,SUMCT

      write(134,5111) a1mean,b1mean,CTmean,
     +  Fxmean/(0.5d0*(pi*wtr(1)**2)*(CTUmean**2.)),
     +  Fxmean/(0.5d0*(pi*wtr(1)**2)*(CTUdmean**2.)),	 
     +  4*a1mean*(1-a1mean),
     +  Fxmean*(z_i**2)/(10**6),Powmean/(10**6),
     +  Fxmean*CTUmean*(z_i**2)/(10**6),	 
     +  v_hub,u_hub,
     +  CTUmean,CTUdmean,
     +  Omegamean	 
!     +  atan(v_hub/u_hub)*180/pi 	 
	  
      write(1095,5544) SUMa1,SUMdA

!c      write(*,5111) a1mean,b1mean,CTmean,
!c     +  Fxmean/(0.5d0*(pi*wtr(1)**2)*(CTUmean**2.)),
!c     +  Fxmean/(0.5d0*(pi*wtr(1)**2)*(CTUdmean**2.)),	 
!c     +  4*a1mean*(1-a1mean),
!c     +  Fxmean*(z_i**2)/(10**6),Powmean/(10**6),
!c     +  Fxmean*CTUmean*(z_i**2)/(10**6),	 
!c     +  v_hub,u_hub,
!c     +  CTUmean,CTUdmean,
!c     +  Omegamean,	 
!c     +  atan(v_hub/u_hub)*180/pi 	 	 
	  
      j=1
      write(*,5111) Power(j)/(10**6),SUMCTU(j),SUMCT(j),
!     +  4.0*SUMa1(j)*(1-SUMa1(j)),SUMCTU(j), 	  
!     +  SUMCTUd(j),1-(SUMCTU(j)/SUMCTUd(j)),
!     +  SUMFX(j)/(0.5*SUMCTU(j)**2*(pi*wtR(1)**2.d0)),
!     +  umean(j),wtomega(j),umeand(j),
     +  SUMdA(j),wtomega(j),SUMa1(j)

!      write(*,5111)SUMa1,SUMb1,SUMCT,1-(SUMCTU/SUMCTUd),
!     +                SUMCT*(1-SUMa1)**2,Power/(10**6)

!      print*, umean(j),umeand(j),wtomega(j),Power(j)

      endif 
	 
      end if

 5111 format (60(1x,f11.5)) 
 5544 format (60(f11.5,1x))

 
      return

      do ii=1,num_turbine
      do k=2,nzb+1
      do j=1,ny
      	tmp = Fx(wtx(ii),j,k)
        do i=1,nx
      	Fx(i,j,k)=tmp*ker_x(ii,i)
        end do

      	tmp = Fy(wtx(ii),j,k)
        do i=1,nx
      	Fy(i,j,k)=tmp*ker_x(ii,i)
        end do

      	tmp = Fz(wtx(ii),j,k)
        do i=1,nx
      	Fz(i,j,k)=tmp*ker_x(ii,i)
        end do

      end do
      end do
      end do



      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine CALALM(u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'dimen.h'
      integer*4 ii,i,j,k,t,M,blade,N,me_k

      real*8, dimension(Nx,Ny,Nz2)::u,v,w,fx,fy,fz,fax,faz,AveWt
      real*8, dimension(num_turbine,Ny,Nz2):: CTy,CTzu,CTzw,CTru,CTrw
      real*8, dimension(num_turbine,Ny,Nz2):: CTangx,CTangz
      real*8, dimension(num_turbine,Nx):: ker_x
      real*8 tmp1,tmp2,tmp3

      real*8, dimension(num_turbine,N_ang,N_blade)::
     +uframe,vframe,wframe,meframe,u1frame,v1frame,w1frame,me1frame,
     +uframed,u1framed	 
      
      real*8, dimension(num_turbine,N_ang,N_blade,Ny,Nz2):: ker_u,ker_w

      real*8 SUMCT,SUMa1,SUMdA,SUMFx,SUMQ,SUMb1,SUMCTd

      real*8 a1,b1,at1,Tang(N_ang),blade_jj,blade_kk
      real*8 mu_y,mu_z,sigma,dr,omega,omega1,r,CTphi
      real*8 tmpy1,tmpz1,tmpy2,tmpz2,tmp
      real*8 ERRd1,ERRd2,tmpd1,tmpd2
      real*8 sum_ker_u,sum_ker_w,max_ker_u,max_ker_w
      real*8 CTU,CTU1,CTV,CTW,chordl,CTsoli
      real*8 CTV1,CTV2,CTUrel,CTFx,CTFt
      real*8 CTCL,CTCD,dang,Re,dA,RA,dL
      real*8 umean(num_turbine),umeand(num_turbine)	  

      real*8, external :: chord,calang
      real*8, external :: itp2D,calU

      save CTy,CTzu,CTzw,CTangx,CTangz,ker_x,ker_u,ker_w
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (t.eq.1) then
        do ii=1,num_turbine
        do k=2,nzb+1
        do j=1,ny
          CTy(ii,j,k)    = (wty(ii)-j)*dy+0.0d-50 

          CTzu(ii,j,k)   = (me*nzb+k-1.5  )*dz-Zhub(ii)
          CTzw(ii,j,k)   = CTzu(ii,j,k)-0.5d0*dz

          CTru(ii,j,k)   = sqrt(CTy(ii,j,k)**2.+CTzu(ii,j,k)**2.)
          CTrw(ii,j,k)   = sqrt(CTy(ii,j,k)**2.+CTzw(ii,j,k)**2.)

          CTangx(ii,j,k) = calang(CTy(ii,j,k),CTzu(ii,j,k))
          CTangz(ii,j,k) = calang(CTy(ii,j,k),CTzw(ii,j,k))
        end do
        end do
        end do

        sigma = 1.d0
        tmp1  = 1/dsqrt(2.d0*pi*sigma*sigma)
        tmp2  = -1/(2.d0*sigma*sigma)
        do ii=1,num_turbine
        do i=1,Nx
          ker_x(ii,i) = tmp1*dexp(tmp2*(i-wtx(ii))**2.)
        end do
        tmp1 = sum(ker_x(ii,1:Nx))
        do i=1,Nx
          ker_x(ii,i) = ker_x(ii,i)/tmp1
        end do
        end do
      end if

      if (t.eq.1) then
      dang = 2.d0*pi/dble(N_ang)
      do i=1,N_ang		
      Tang(i) = (dble(i)-1)*dang
      end do

      do ii=1,num_turbine
        dr = wtr(ii)/dble(N_blade)
        do N=1,N_blade
        r = dble(N-0.5)*dr
        do blade=1,N_ang

          sum_ker_u  = 0.0
          sum_ker_w  = 0.0
          mu_y       = r*dcos(Tang(blade))  
          mu_z       = r*dsin(Tang(blade)) 
          sigma      = dsqrt((r*dang)**2.+dr**2.)*2.d0

          tmp1       = 1/dsqrt(2.d0*pi*sigma*sigma)
          tmp2       = -1/(2.d0*sigma*sigma)

          do k=2,nzb+1
          do j=1,ny
          if (CTru(ii,j,k).le.2.0*wtr(ii))then
          tmp3 = (CTy(ii,j,k)-mu_y)**2.+(CTzu(ii,j,k)-mu_z)**2.
          ker_u(ii,blade,N,j,k) = tmp1*dexp(tmp2*tmp3)
          sum_ker_u = sum_ker_u+ker_u(ii,blade,N,j,k)
          else
          ker_u(ii,blade,N,j,k) = 0.d0
          end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
          if (CTrw(ii,j,k).le.2.0*wtr(ii))then
          tmp3 = (CTy(ii,j,k)-mu_y)**2.+(CTzw(ii,j,k)-mu_z)**2.
          ker_w(ii,blade,N,j,k) = tmp1*dexp(tmp2*tmp3)
          sum_ker_w = sum_ker_w+ker_w(ii,blade,N,j,k)
          else
          ker_w(ii,blade,N,j,k) = 0.d0
          end if
          end do
          end do

          call calsum(sum_ker_u,max_ker_u,me,nall)
          call calsum(sum_ker_w,max_ker_w,me,nall)

          do k=2,Nzb+1
          do j=1,Ny
          ker_u(ii,blade,N,j,k) = ker_u(ii,blade,N,j,k)/sum_ker_u
          ker_w(ii,blade,N,j,k) = ker_w(ii,blade,N,j,k)/sum_ker_w
          end do
          end do

        end do
        end do
      end do
      end if

      meframe=0.d0
      uframe=0.d0
      vframe=0.d0
      wframe=0.d0
      uframed=0.d0	  

      do ii=1,num_turbine
      dr = wtr(ii)/dble(N_blade)
      do N=1,N_blade
      r = dble(N-0.5)*dr
      do blade=1,N_ang
      blade_kk = dble(Zhub(ii)/dz+0.5)+r*dsin(Tang(blade))/dz-1.d0
      me_k=blade_kk/Nzb

      if (me.eq.me_k) then
      blade_jj = dble(wty(ii))-r*dcos(Tang(blade))/dy
      blade_kk = blade_kk-Nzb*me+2.d0
      
      meframe(ii,blade,N)=1.d0
      i = 0.d0*wtR(ii)/dx
      uframe(ii,blade,N)=itp2D(u,wtx(ii)-i,blade_jj,blade_kk,me,nall)
      i = 6*wtR(ii)/dx
      uframed(ii,blade,N)=itp2D(u,wtx(ii)-i,blade_jj,blade_kk,me,nall)	  
      i = 0
      vframe(ii,blade,N)=itp2D(v,wtx(ii)-i,blade_jj,blade_kk,me,nall)
      i = 0
      wframe(ii,blade,N)=itp2D(w,wtx(ii)-i,blade_jj,blade_kk,me,nall)
      end if
      end do
      end do
      end do

      IF (me>0) then
         call MPI_SEND(meframe(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(uframe(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(uframed(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	 
         call MPI_SEND(vframe(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(wframe(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
      ELSE
         do i=1,nprocs-1
         call MPI_RECV(me1frame(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         call MPI_RECV(u1frame(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         call MPI_RECV(u1framed(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	 
         call MPI_RECV(v1frame(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         call MPI_RECV(w1frame(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )

         tmp1=sum(me1frame)
         if(tmp1.gt.0.5)then
           do ii=1,num_turbine
           do N=1,N_blade
           do blade=1,N_ang
           if(me1frame(ii,blade,N).gt.0.5d0)then
           uframe(ii,blade,N)=u1frame(ii,blade,N)
           uframed(ii,blade,N)=u1framed(ii,blade,N)
           vframe(ii,blade,N)=v1frame(ii,blade,N)
           wframe(ii,blade,N)=w1frame(ii,blade,N)
           end if
           end do
           end do
           end do
         end if
         end do

         do i=1,nprocs-1
         call MPI_SEND(uframe(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )
         call MPI_SEND(uframed(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )	 
         call MPI_SEND(vframe(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )
         call MPI_SEND(wframe(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )
         end do
      END IF
      IF (me>0) then
          call MPI_RECV(uframe(1,1,1),num_turbine*N_ang*N_blade,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
          call MPI_RECV(uframed(1,1,1),num_turbine*N_ang*N_blade,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )	 
          call MPI_RECV(vframe(1,1,1),num_turbine*N_ang*N_blade,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
          call MPI_RECV(wframe(1,1,1),num_turbine*N_ang*N_blade,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
      END IF

      dang = 2.d0*pi/dble(N_ang)

      SUMdA=0.0
      SUMa1=0.0
      SUMb1=0.0	  
      SUMCT=0.0
      SUMCTd=0.0	  
      SUMFx=0.0	
      SUMQ=0.0
	    
      do ii=1,num_turbine
	  
      umean(ii)=sum(uframe(ii,:,:))/(dble(N_ang)*N_blade)
      umeand(ii)=sum(uframed(ii,:,:))/(dble(N_ang)*N_blade)
	  
      dr = wtr(ii)/dble(N_blade)
      do N=N_blade_s,N_blade
      r = dble(N-0.5)*dr
      do M=1,3

      blade=t+(M-1)*N_ang/3
cc      blade=-t+(M-1)*N_ang/3+((t/N_ang)+1)*N_ang
      blade=mod(blade-1,N_ang)+1

      if(me.eq.0.and.N.eq.N_blade)then
      	write(*,*)ii,M,blade,N_ang
      end if

      CTU = uframe(ii,blade,N)
      CTV = vframe(ii,blade,N)
      CTW = wframe(ii,blade,N)

      chordl  = chord(r*z_i)/z_i
      CTsoli  = 3.0*chordl/(2.d0*pi*r)

      omega1  = +CTW*dcos(Tang(blade))/r
     +          +CTV*dsin(Tang(blade))/r
     +          +wtomega(ii)*(2.d0*pi/60)*(z_i/u_star)

      call sol_CL_CD_1(omega1,chordl,r,CTU,CTsoli,
     +     a1,b1,CTCL,CTCD,at1,CTphi,Re,me,nall)	  
      CTUrel  = CTU/dsin(CTphi)	  	
	  
c      call sol_CL_CD(omega1,chordl,r,CTU,CTsoli,
c     +     a1,b1,CTCL,CTCD,at1,CTphi,Re,me,nall)
c      CTV1    =  CTU
c      CTV2    = -CTW*dcos(Tang(blade))
c     +            -CTV*dsin(Tang(blade))-omega1*r
c      CTUrel  = CTU*(1-a1)/dsin(CTphi)

      dA      = dr*chordl

      CTFx    = 0.5*(CTUrel**2)*dA/(dx*dy*dz)*
     +          (+CTCL*dcos(CTphi)+CTCD*dsin(CTphi))
      CTFt    =-0.5*(CTUrel**2)*dA/(dx*dy*dz)*
     +      dabs(+CTCL*dsin(CTphi)-CTCD*dcos(CTphi))

      SUMdA =SUMdA+dA
      SUMa1 =SUMa1+a1*dA
      SUMCT =SUMCT+CTFx/(0.5d0*((umean(ii)/u_star)**2.)*dA/(dx*dy*dz))*
     +      dA/(pi*wtR(ii)**2.d0)
      SUMCTd =SUMCTd+CTFx/(0.5d0*((umeand(ii)/u_scale)**2.)*
     +dA/(dx*dy*dz))*dA/(pi*wtR(ii)**2.d0) 	 
      SUMFx=SUMFx+CTFx*(dx*dy*dz)
      SUMQ=SUMQ-CTFt*(dx*dy*dz)*r*
     +         wtomega(ii)*(2.d0*pi/60)*(z_i**3)

      i=wtx(ii)
      do k=2,nzb+1
      do j=1,ny
         if (ker_u(ii,blade,N,j,k).ge.1E-15)then
         Fx (i,j,k)=Fx (i,j,k) + CTFx  * ker_u(ii,blade,N,j,k)

         Fax(i,j,k)=Fax(i,j,k) + CTFt  * ker_u(ii,blade,N,j,k)
         Fy (i,j,k)=Fy (i,j,k) + CTFt  * ker_u(ii,blade,N,j,k)
     +             *dsin(CTangx(ii,j,k))
         end if

         if (ker_w(ii,blade,N,j,k).ge.1E-15)then
         Faz(i,j,k)=Faz(i,j,k) + CTFt  * ker_w(ii,blade,N,j,k)
         Fz (i,j,k)=Fz (i,j,k) + CTFt  * ker_w(ii,blade,N,j,k)
     +             *dcos(CTangz(ii,j,k))
         end if
      end do
      end do

      end do
      end do
      end do

      SUMa1=SUMa1/SUMdA
cc      SUMCT=SUMCT/SUMdA*(pi*wtR(1)**2.d0)
cc      SUMCT=SUMCT	  
      if(me.eq.0)then
      if(t.eq.1)then
      Open (unit=1095,file='output/check_a_CT.out',position='append')
      end if
      write(1095,5111)SUMa1,SUMCTd,SUMCT,umean,umeand,1-(umean/umeand),
     +       SUMQ/(10**6),SUMFx*(z_i**2)/(10**6)
c      write(*,*)SUMa1,SUMCT,SUMdA/(pi*wtR(1)**2.d0)/num_turbine

      end if

      return

      do ii=1,num_turbine
      do k=2,nzb+1
      do j=1,ny
      	tmp = Fx(wtx(ii),j,k)
        do i=1,nx
      	Fx(i,j,k)=tmp*ker_x(ii,i)
        end do

      	tmp = Fy(wtx(ii),j,k)
        do i=1,nx
      	Fy(i,j,k)=tmp*ker_x(ii,i)
        end do

      	tmp = Fz(wtx(ii),j,k)
        do i=1,nx
      	Fz(i,j,k)=tmp*ker_x(ii,i)
        end do

      end do
      end do
      end do

 5111 format (60(1x,f11.5)) 
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine CALVAWT(u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'dimen.h'
      integer*4 ii,i,j,k,t,M,blade,N,me_k

      real*8, dimension(Nx,Ny,Nz2)::u,v,w,fx,fy,fz,fax,faz,AveWt
      real*8, dimension(num_turbine_vawt,Nx,Ny,Nz2):: CTy,CTx
      real*8, dimension(num_turbine_vawt,Nx,Ny,Nz2):: CTzu,CTzw,CTru,CTrw	  
      real*8, dimension(num_turbine_vawt,Ny,Nz2):: CTangx,CTangz
      real*8, dimension(num_turbine_vawt,Nx):: ker_x
      real*8 tmp1,tmp2,tmp3

      real*8, dimension(num_turbine_vawt,N_ang_vawt,N_blade)::
     +uframe,vframe,wframe,meframe,u1frame,v1frame,w1frame,me1frame,
     +uframed,u1framed	 
      
      real*8, dimension(num_turbine_vawt,N_ang_vawt,Nx,Ny,Nz2):: 
     +        ker_u_vawt,ker_w_vawt
cc      real*8, dimension(num_turbine_vawt,N_ang,Nx,Ny):: ker_u,ker_w	  

      real*8, dimension(num_turbine_vawt)::SUMCT,SUMa1,SUMdA,SUMFx,
     +       SUMQ,SUMb1,SUMCTd,SUMQ1,SUMCTU,Power,Power1,
     +	     SUMCT1,SUMCTU1,SUMCTd1,SUMdA1,SUMFx1,
     +       SUMFr,SumFr1,SUMFt,SumFt1,Sumat,Sumat1,
     +       SumRe,SumRe1	 

      real*8 a1,b1,at1,Tang(N_ang_vawt),blade_jj,blade_kk
      real*8 mu_y,mu_x,sigma,dr,omega,omega1,r,CTphi
      real*8 tmpy1,tmpz1,tmpy2,tmpz2,tmp
      real*8 ERRd1,ERRd2,tmpd1,tmpd2
      real*8 sum_ker_u(Nz2),sum_ker_w,max_ker_u,max_ker_w
      real*8 CTU,CTU1,CTV,CTW,chordl,CTsoli
      real*8 CTV1,CTV2,CTUrel,CTFx,CTFt,CTFr,CTFy
      real*8 CTCL,CTCD,dang,Re,dA,RA,dL
      real*8 umean(num_turbine_vawt),umeand(num_turbine_vawt)	 
      real*8 AoA_old(num_turbine_vawt,3,Nz2),
     +       AoA_new(num_turbine_vawt,3,Nz2)
      real*8 AoA_rate(num_turbine_vawt,3,Nz2),
     +       AoA_ds(num_turbine_vawt,3,Nz2)	
      real*8 CL_new(num_turbine_vawt,3,Nz2),
     +       CD_new(num_turbine_vawt,3,Nz2)	  
      real*8 lambdar,at,CL,CD,phi,CL_max,beta,at_ss,CL_ss 	  
      real*8, external :: CALCL_vawt,CALCD_vawt	  

      real*8, external :: chord,calang
      real*8, external :: itp2D,calU

      real*8 rr(4),yy1,xx1,rr1(4),wp(4),atd,qr,q,Re_c,num	
      integer*4 n1y,n2y,n1x,n2x,ii1,i1,DS	  
	  
      save CTy,CTzu,CTzw,CTangx,CTangz,ker_x,ker_u_vawt,ker_w_vawt
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (t.eq.1) then
        do ii=1,num_turbine_vawt
        do k=2,nzb+1
        do j=1,ny
        do i=1,nx		
		
          CTy(ii,i,j,k)    = (wty_vawt(ii)-j)*dy+0.0d-50 
          CTx(ii,i,j,k)    = (wtx_vawt(ii)-i)*dx+0.0d-50		  
          CTru(ii,i,j,k)   = sqrt(CTy(ii,i,j,k)**2.+CTx(ii,i,j,k)**2.)

        end do
        end do
        end do
        end do		
      end if

      if (t.eq.1) then
	 	  
      dang = 2.d0*pi/dble(N_ang_vawt)
      do i=1,N_ang_vawt		
      Tang(i) = (dble(i)-1)*dang
      end do

      ker_u_vawt=0.0d0	  
      do ii=1,num_turbine_vawt  	  
        do blade=1,N_ang_vawt

          sum_ker_u  = 0.0
          mu_x = wtr_vawt(ii)*dcos(Tang(blade))  
          mu_y = wtr_vawt(ii)*dsin(Tang(blade)) 

          do k=2,nzb+1
          do j=1,ny
          do i=1,nx		  
          if (CTru(ii,i,j,k).le.2.0*wtr_vawt(ii))then
          tmp3 = (CTy(ii,i,j,k)-mu_y)**2.+(CTx(ii,i,j,k)-mu_x)**2.
          ker_u_vawt(ii,blade,i,j,k) = 1/(((2.5*dx)**3)*(pi**1.5))
     +                           *dexp(-1.0*tmp3/((2.5*dx)**2.0))
          sum_ker_u(k) = sum_ker_u(k)+ker_u_vawt(ii,blade,i,j,k)
          else
          ker_u_vawt(ii,blade,i,j,k) = 0.d0
          end if
          end do
          end do
          end do
	  	  
          do k=2,Nzb+1
          do j=1,Ny
          do i=1,Nx		  
          ker_u_vawt(ii,blade,i,j,k) = ker_u_vawt(ii,blade,i,j,k)
     +                                 /sum_ker_u(k)
          end do
          end do
          end do

        end do
      end do
      end if

	  
!      do k=2,nzb+1
!      do j=1,Ny	 
!      do i=5,5	  
!      dA=dz*dy 
!c      n1x=i
!      CTU=u(i,j,k)			
!      CTV=v(i,j,k)
	         
!cc      CALL RANDOM_NUMBER(num)
!c      print*, 'num', num	  
!cc      CTFx=0.5*num*4.0d0*(CTU**2)*dA/(dx*dy*dz) 
!cc      Fx(i,j,k)=Fx(i,j,k)+CTFx	

!      end do
!      end do	
!      end do	  
  
      dang = 2.d0*pi/dble(N_ang_vawt)

c      SUMdA=0.0
c      SUMa1=0.0
c      SUMb1=0.0	  
c      SUMCT=0.0
c      SUMCTU=0.0	  
c      SUMCTd=0.0	  
c      SUMFx=0.0	
c      SUMQ=0.0
c      Power=0.0	  
CCCCCCCCCCCCCCCCCCCCCCCCCCC
c      do ii=1,num_turbine_vawt
c      do k=2,nzb+1
c	  
c       if ((me*nzb+k-1.5)*dz.gt.(Zhub_vawt(ii)-0.5*wtH_vawt(ii)) .and.
c     +    (me*nzb+k-1.5)*dz.lt.(Zhub_vawt(ii)+0.5*wtH_vawt(ii))) then
c	                                                      	  
c       do blade=1,N_ang_vawt
c
c        xx1= wtx_vawt(ii)-(wtr_vawt(ii)*dcos(Tang(blade))/dx)  
c        yy1= wty_vawt(ii)-(wtr_vawt(ii)*dsin(Tang(blade))/dy) 
c        n1x=floor(xx1)
c        n2x=n1x+1
c        n1y=floor(yy1)
c        n2y=n1y+1
c        rr1(1)=(((yy1-n1y)*dy)**2+((xx1-n1x)*dx)**2)**0.5
c        rr1(2)=(((yy1-n1y)*dy)**2+((xx1-n2x)*dx)**2)**0.5
c        rr1(3)=(((yy1-n2y)*dy)**2+((xx1-n1x)*dx)**2)**0.5
c        rr1(4)=(((yy1-n2y)*dy)**2+((xx1-n2x)*dx)**2)**0.5
c        ii1=0
c        do i1=1,4
c        if (rr1(i1).le.1.0d-6) then
c        wp=0.0
c        wp(i1)=1.0
c        ii1=1 
c        endif
c        enddo
c        if (ii1==0) then
c         do i1=1,4	  
c         wp(i1)=(1/rr1(i1))/(1/rr1(1)+1/rr1(2)+1/rr1(3)+1/rr1(4))		 
c         enddo
c        endif	        
c        CTU=wp(1)*u(n1x,n1y,k)+wp(2)*u(n2x,n1y,k)+
c     +      wp(3)*u(n1x,n2y,k)+wp(4)*u(n2x,n2y,k)			
c        CTV=wp(1)*v(n1x,n1y,k)+wp(2)*v(n2x,n1y,k)+
c     +      wp(3)*v(n1x,n2y,k)+wp(4)*v(n2x,n2y,k)		
 
cc       	CTU=0.93d0*CTU 
cc       	CTV=0.93d0*CTV
		
c      dA      = wtr_vawt(ii)*dang*dz
c      SUMdA(ii) =SUMdA(ii)+dA*(z_i**2.)	
c      SUMCTU(ii)=SUMCTU(ii)+CTU*dA/(2.0*pi*wtR_vawt(ii)*wtH_vawt(ii))	  

c       end do
	  
c       endif
c
c      end do
c      IF (me>0) then
c         call MPI_SEND(SUMCTU(ii),1,
c     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	 
c      ELSE
c         do i=1,nprocs-1	  
c         call MPI_RECV(SUMCTU1(ii),1,
c     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	 
c         SUMCTU(ii)=SUMCTU(ii)+SUMCTU1(ii)	
c         enddo		 
c      ENDIF

c      IF (me.eq.0) then
cc         wtomega_vawt(ii)=33.5d0
c         wtomega_vawt(ii)=1.2d0*SUMCTU(ii)*3.8/(wtR_vawt(ii)*z_i)*
c     +                    60.d0/(2.0*pi)      
c         do i=1,nprocs-1
c          call MPI_SEND(wtomega_vawt(ii),1, MPI_DOUBLE_PRECISION,
c     +           i,i,nall,ierr ) 
c         enddo
c      ENDIF
c      if (me.ne.0) then	  
c      call MPI_RECV(wtomega_vawt(ii),1, MPI_DOUBLE_PRECISION,
c     +                  0,me,nall,status2,ierr )
c      endif	 
 		
c      end do 	  
CCCCCCCCCCCCCCCCCCCCCCCCCCC	  
      SUMdA=0.0
      SUMa1=0.0
      SUMb1=0.0	  
      SUMCT=0.0
      SUMCTU=0.0	  
      SUMCTd=0.0	  
      SUMFx=0.0	
      SUMQ=0.0
      SUMFr=0.0
      SUMFt=0.0	  
      Power=0.0
      SumRe=0.0
      Sumat=0.0 	  
CCCCCCCCCCCCCCCCCCCCCCCCCCC
	  
      if (t.eq.1) AoA_old=0.0 	  
  
cc      if (me.ge.7 .and. me.le.8) then
      do ii=1,num_turbine_vawt
      do k=2,nzb+1
	  
       if ((me*nzb+k-1.5)*dz.gt.(Zhub_vawt(ii)-0.5*wtH_vawt(ii)) .and.
     +    (me*nzb+k-1.5)*dz.lt.(Zhub_vawt(ii)+0.5*wtH_vawt(ii))) then
	                                                      	  
c       do blade=1,N_ang_vawt

      do M=1,N_blade_vawt
        if (dir_vawt(ii).gt.0.d0) then 
        blade=t+(M-1)*N_ang_vawt/N_blade_vawt
        blade=mod(blade-1,N_ang_vawt)+1
        else		
        blade=-t+(M-1)*N_ang_vawt/N_blade_vawt+
     +         ((t/N_ang_vawt)+1)*N_ang_vawt
        blade=mod(blade-1,N_ang_vawt)+1		
        endif 
		
        xx1= wtx_vawt(ii)-(wtr_vawt(ii)*dcos(Tang(blade))/dx)  
        yy1= wty_vawt(ii)-(wtr_vawt(ii)*dsin(Tang(blade))/dy) 
        n1x=floor(xx1)
        n2x=n1x+1
        n1y=floor(yy1)
        n2y=n1y+1
        rr1(1)=(((yy1-n1y)*dy)**2+((xx1-n1x)*dx)**2)**0.5
        rr1(2)=(((yy1-n1y)*dy)**2+((xx1-n2x)*dx)**2)**0.5
        rr1(3)=(((yy1-n2y)*dy)**2+((xx1-n1x)*dx)**2)**0.5
        rr1(4)=(((yy1-n2y)*dy)**2+((xx1-n2x)*dx)**2)**0.5
        ii1=0
        do i1=1,4
        if (rr1(i1).le.1.0d-6) then
        wp=0.0
        wp(i1)=1.0
        ii1=1 
        endif
        enddo
        if (ii1==0) then
         do i1=1,4	  
         wp(i1)=(1/rr1(i1))/(1/rr1(1)+1/rr1(2)+1/rr1(3)+1/rr1(4))		 
         enddo
        endif	        
        CTU=wp(1)*u(n1x,n1y,k)+wp(2)*u(n2x,n1y,k)+
     +      wp(3)*u(n1x,n2y,k)+wp(4)*u(n2x,n2y,k)			
        CTV=wp(1)*v(n1x,n1y,k)+wp(2)*v(n2x,n1y,k)+
     +      wp(3)*v(n1x,n2y,k)+wp(4)*v(n2x,n2y,k)
	 
c       	CTU=0.93d0*CTU 
c       	CTV=0.93d0*CTV
		
c      chordl  = chord(r*z_i)/z_i
       chordl  = chord_vawt(ii)	  
       CTsoli  = N_blade_vawt*chordl/(2.d0*pi*wtr_vawt(ii))

       omega1  = wtomega_vawt(ii)*(2.d0*pi/60)*(z_i/u_star)
c       lambdar=(wtr_vawt(ii)*omega1-dir_vawt(ii)*CTU*dsin(Tang(blade))+
c     +         dir_vawt(ii)*CTV*dcos(Tang(blade)))/
c     +        (CTU*dcos(Tang(blade))+CTV*dsin(Tang(blade)))
       CTphi = datan((CTU*dcos(Tang(blade))+1.0*CTV*dsin(Tang(blade)))/
     +	  (wtr_vawt(ii)*omega1-dir_vawt(ii)*CTU*dsin(Tang(blade))+
     +         dir_vawt(ii)*1.0*CTV*dcos(Tang(blade))))	

c       CTphi=0.9*CTphi	 
       at  = CTphi

       CTCL   = 1.0*CALCL_vawt(at*180.d0/pi,wtr_vawt(ii)*z_i)
       CTCD   = 1.0*CALCD_vawt(at*180.d0/pi,wtr_vawt(ii)*z_i)

       CTUrel=sqrt((CTU*dcos(Tang(blade))+1.0*CTV*dsin(Tang(blade)))**2.
     +       +(wtr_vawt(ii)*omega1-dir_vawt(ii)*CTU*dsin(Tang(blade))+
     +         dir_vawt(ii)*1.0*CTV*dcos(Tang(blade)))**2.)
 
       Re_c=chordl*CTUrel/nu  
  
       qr=(abs((me*nzb+k-1.5)*dz-Zhub_vawt(ii)))/
     +    (0.5*wtH_vawt(ii))+1.0d-50
cc       q=2./pi*dacos(
cc     +   dexp(-(N_blade_vawt/2.*(1.-qr))/(qr*abs(dsin(CTphi)))))
       q=2./pi*dacos(
     +   dexp(-(N_blade_vawt/2.*(1.-qr))/(abs(dsin(CTphi)))))	 
!       q=1.0d0
c       if (M.eq.1) print*, q,me
CCCCCCCCCCCCCCCCCCCCCCCCCCC DYNAMIC STALL
      beta=0.077d0	   
      at_ss=14.0d0
      CL_ss=1.08d0
CCCCCC Dynamic Stall CCCCCC  	  
       AoA_new(ii,M,k)=at*180.d0/pi 
       AoA_rate(ii,M,k)=(abs(AoA_new(ii,M,k))-abs(AoA_old(ii,M,k)))
     + 	                *pi/180.d0/dt*chordl/(2.0*CTUrel)	 
       if (t.eq.1) AoA_rate(ii,M,k)=0.0 

       AoA_ds(ii,M,k)=0.0d0
       if (AoA_rate(ii,M,k).gt.0.0d0) then
          if (AoA_new(ii,M,k).gt.0.0d0) then
       AoA_ds(ii,M,k)=at_ss+(1.0*(AoA_rate(ii,M,k))**0.5)*180.d0/pi
          else if (AoA_new(ii,M,k).lt.0.0d0) then
       AoA_ds(ii,M,k)=-at_ss-(1.0*(AoA_rate(ii,M,k))**0.5)*180.d0/pi		  
          endif
       endif

       DS=1
       if (DS.eq.1) then 
       if (AoA_rate(ii,M,k).gt.0.0d0) then
        if (AoA_new(ii,M,k).gt.at_ss .and. 
     +       AoA_new(ii,M,k).lt.AoA_ds(ii,M,k)) then	   
c           	 CTCL=CL_ss+beta*(AoA_new(ii,M,k)-at_ss)
c             CL_max=CL_ss+2.0d0
           	 CTCL=CL_ss/(at_ss*pi/180.d0)*dsin(at)
             CL_max=CL_ss+2.0d0
         if (40.*2.*AoA_rate(ii,M,k).lt.2.0) then
c       if ((CL_ss+40.*2.*AoA_rate(ii,M,k)).gt.
c     +     beta*(AoA_ds(ii,M,k)-at_ss)) then   
                CL_max=CL_ss+40.0*2.*AoA_rate(ii,M,k) 
c             else				
c                CL_max=CL_ss+beta*(AoA_ds(ii,M,k)-at_ss)				
c       endif
         endif			 
             if (CTCL .ge. CL_max) CTCL=CL_max         			 
        endif
       endif		 
       if (AoA_rate(ii,M,k).gt.0.0d0) then
         if (AoA_new(ii,M,k).gt.at_ss .and.
     +       AoA_new(ii,M,k).gt.AoA_ds(ii,M,k)) then
c           	 CTCL=CL_ss+beta*(AoA_new(ii,M,k)-at_ss)
c             CL_max=CL_ss+2.0d0
           	 CTCL=CL_ss/(at_ss*pi/180.d0)*dsin(at)
             CL_max=CL_ss+2.0d0
        if (40.*2.*AoA_rate(ii,M,k).lt.2.0) then
c        if (40.*2.*AoA_rate(ii,M,k).gt.beta*(AoA_ds(ii,M,k)-at_ss)) then   
                CL_max=CL_ss+40.0*2.*AoA_rate(ii,M,k) 
c             else				
c                CL_max=CL_ss+beta*(AoA_ds(ii,M,k)-at_ss)				
c        endif
        endif			
             if (CTCL .ge. CL_max) CTCL=CL_max
c             CTCD=abs(CTCL*dsin(AoA_new(ii,M,k)*pi/180.0d0))      
             CTCD=abs(CTCL*dtan(AoA_new(ii,M,k)*pi/180.0d0))  			 
         endif
       endif
	   
       if (AoA_rate(ii,M,k).gt.0.0d0) then
         if (AoA_new(ii,M,k).lt.-at_ss .and. 
     +       AoA_new(ii,M,k).gt.AoA_ds(ii,M,k)) then	   
c          	 CTCL=-CL_ss-beta*(-AoA_new(ii,M,k)-at_ss)
c             CL_max=-CL_ss-2.0d0
           	 CTCL=CL_ss/(at_ss*pi/180.d0)*dsin(at)
             CL_max=-CL_ss-2.0d0			 
			 
       if (40.*2.*AoA_rate(ii,M,k).lt.2.0) then
c       if (40.*2.*AoA_rate(ii,M,k).gt.beta*(-AoA_ds(ii,M,k)-at_ss)) then    
                CL_max=-CL_ss-40.0*2.*AoA_rate(ii,M,k) 
c             else				
c                CL_max=-CL_ss-beta*(-AoA_ds(ii,M,k)-at_ss)				
c       endif
       endif			
             if (CTCL .le. CL_max) CTCL=CL_max         			 
         endif
       endif		 
       if (AoA_rate(ii,M,k).gt.0.0d0) then
         if (AoA_new(ii,M,k).lt.-at_ss .and.
     +       AoA_new(ii,M,k).lt.AoA_ds(ii,M,k)) then
c          	 CTCL=-CL_ss-beta*(-AoA_new(ii,M,k)-at_ss)
c             CL_max=-CL_ss-2.0d0
           	 CTCL=CL_ss/(at_ss*pi/180.d0)*dsin(at)
             CL_max=-CL_ss-2.0d0	
			 
       if (40.*2.*AoA_rate(ii,M,k).lt.2.0) then
c       if (40.*2.*AoA_rate(ii,M,k).gt.beta*(-AoA_ds(ii,M,k)-at_ss)) then   
                CL_max=-CL_ss-40.0*2.*AoA_rate(ii,M,k) 
c             else				
c                CL_max=-CL_ss-beta*(-AoA_ds(ii,M,k)-at_ss)				
c       endif
       endif					
             if (CTCL .le. CL_max) CTCL=CL_max
c             CTCD=abs(CTCL*dsin(AoA_new(ii,M,k)*pi/180.0d0))   
             CTCD=abs(CTCL*dtan(AoA_new(ii,M,k)*pi/180.0d0))           			 
         endif
       endif

       CL_new(ii,M,k)=CTCL
       CD_new(ii,M,k)=CTCD	   
	   
       endif	
       AoA_old(ii,M,k)=AoA_new(ii,M,k)	   
CCCCCCCCCCCCCCCCCCCCCCCCCCC END OF DYNAMIC STALL
CCCCCC NEW Dynamic Stall CCCCCC  	  
       DS=0
       if (DS.eq.10) then 
       if (AoA_rate(ii,M,k).gt.0.0d0) then
         if (AoA_new(ii,M,k).gt.0. .and. 
     +       AoA_new(ii,M,k).lt.30.) then	   
           	 CTCL=AoA_new(ii,M,k)*0.073d0		 
         endif
       endif
       if (AoA_rate(ii,M,k).lt.0.0d0) then
         if (AoA_new(ii,M,k).gt.0. .and. 
     +       AoA_new(ii,M,k).lt.30.) then	   
c           	 CTCL=AoA_new(ii,M,k)*0.02d0		 
         endif
       endif
       if (AoA_rate(ii,M,k).gt.0.0d0) then
         if (AoA_new(ii,M,k).lt.(-0.) .and. 
     +       AoA_new(ii,M,k).gt.(-20.)) then	   
           	 CTCL=AoA_new(ii,M,k)*0.073d0		 
         endif
       endif
       if (AoA_rate(ii,M,k).lt.0.0d0) then
         if (AoA_new(ii,M,k).lt.(-0.) .and. 
     +       AoA_new(ii,M,k).gt.(-30.)) then	   
c           	 CTCL=AoA_new(ii,M,k)*0.02d0		 
         endif
       endif

       CL_new(ii,M,k)=CTCL
       CD_new(ii,M,k)=CTCD	   
	   
       endif	      
CCCCCCCCCCCCCCCCCCCCCCCCCCC END OF DYNAMIC STALL
      dA      = dz*chordl
      CTFr    = q*0.5*(CTUrel**2)*dA/(dx*dy*dz)*
     +          (+CTCL*dcos(CTphi)+CTCD*dsin(CTphi))
      CTFt    = q*dir_vawt(ii)*0.5*(CTUrel**2)*dA/(dx*dy*dz)*
     +          (+CTCL*dsin(CTphi)-CTCD*dcos(CTphi))	
  
cc     +      dabs(+CTCL*dsin(CTphi)-CTCD*dcos(CTphi))

c       dA      = wtr_vawt(ii)*dang*dz
c       CTFr    = 0.5*(CTUrel**2)*CTsoli*dA/(dx*dy*dz)*
c     +          (+CTCL*dcos(CTphi)+CTCD*dsin(CTphi))
c       CTFt    = dir_vawt(ii)*0.5*(CTUrel**2)*CTsoli*dA/(dx*dy*dz)*
c     +          (+CTCL*dsin(CTphi)-CTCD*dcos(CTphi))	 

       CTFx=CTFr*dcos(Tang(blade))+CTFt*dsin(Tang(blade))     	 
       CTFy=CTFr*dsin(Tang(blade))-CTFt*dcos(Tang(blade)) 

      SUMdA(ii) =SUMdA(ii)+dA*(z_i**2.)
      if (CTU.ne.0.0d0) then	  
      SUMCT(ii) =SUMCT(ii)+CTFx/(0.5d0*((CTU/u_star)**2.)*	      
     +           dA/(dx*dy*dz))*dA/(2.0*pi*wtR_vawt(ii)*wtH_vawt(ii))
      else
      SUMCT(ii) =SUMCT(ii)+0.0d0
      endif	  
      SUMCTd(ii) =SUMCTd(ii)+CTFx/(0.5d0*((7.04d0/u_star)**2.)*	      
     +           dA/(dx*dy*dz))*dA/(2.0*pi*wtR_vawt(ii)*wtH_vawt(ii))	 
      SUMCTU(ii)=SUMCTU(ii)+CTU*dA/(2.0*pi*wtR_vawt(ii)*wtH_vawt(ii))	  

      SUMQ(ii)=SUMQ(ii)+dir_vawt(ii)*CTFt*(dx*dy*dz)*wtr_vawt(ii)*
     +         (z_i**3)		

      SUMFx(ii)=SUMFx(ii)+CTFx*(dx*dy*dz)*(z_i**2.)	 	 
      if (M.eq.1) then 	   
      SUMFr(ii)=SUMFr(ii)+CTFr*(dx*dy*dz)*(z_i**2)	
      SUMFt(ii)=SUMFt(ii)+CTFt*(dx*dy*dz)*(z_i**2)	
      Sumat(ii)=Sumat(ii)+
     +          at*180.0/pi*dA/(chordl*wtH_vawt(ii))
      SumRe(ii)=SumRe(ii)+Re_c*dA/(chordl*wtH_vawt(ii))  
      endif 
	  
      Fx (n1x,n1y,k)=Fx (n1x,n1y,k) + CTFx*wp(1)	
      Fx (n2x,n1y,k)=Fx (n2x,n1y,k) + CTFx*wp(2)	
      Fx (n1x,n2y,k)=Fx (n1x,n2y,k) + CTFx*wp(3)	
      Fx (n2x,n2y,k)=Fx (n2x,n2y,k) + CTFx*wp(4)			 
      Fy (n1x,n1y,k)=Fy (n1x,n1y,k) + CTFy*wp(1)	
      Fy (n2x,n1y,k)=Fy (n2x,n1y,k) + CTFy*wp(2)	
      Fy (n1x,n2y,k)=Fy (n1x,n2y,k) + CTFy*wp(3)	
      Fy (n2x,n2y,k)=Fy (n2x,n2y,k) + CTFy*wp(4)	

      do j=1,ny
      do i=1,nx
         if (ker_u_vawt(ii,blade,i,j,k).ge.1E-15)then
c         Fx (i,j,k)=Fx (i,j,k) + 1.*CTFx  * ker_u_vawt(ii,blade,i,j,k)
c         Fy (i,j,k)=Fy (i,j,k) + 1.*CTFy  * ker_u_vawt(ii,blade,i,j,k)
        end if
      end do
      end do
   
      end do
	  
      endif

      end do

      Power(ii)=SUMQ(ii)*wtomega_vawt(ii)*(2.d0*pi/60)	
	  
      end do 	  

c      end if

      IF (me>0) then
         call MPI_SEND(SUMQ,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(SUMCT,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(SUMCTU,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	 
         call MPI_SEND(SUMCTd,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	
         call MPI_SEND(SUMdA,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	
         call MPI_SEND(SUMFx,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	
         call MPI_SEND(Power,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	 
         call MPI_SEND(SUMFr,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	
         call MPI_SEND(SUMFt,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(SUMRe,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	
         call MPI_SEND(SUMat,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	 
      ELSE
         do i=1,nprocs-1	  
         call MPI_RECV(SUMQ1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         call MPI_RECV(SUMCT1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         call MPI_RECV(SUMCTU1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	 
         call MPI_RECV(SUMCTd1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	
         call MPI_RECV(SUMdA1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	 
         call MPI_RECV(SUMFx1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	 
         call MPI_RECV(Power1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	
         call MPI_RECV(SUMFr1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	
         call MPI_RECV(SUMFt1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	
         call MPI_RECV(SUMRe1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	
         call MPI_RECV(SUMat1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )		 
         SUMQ=SUMQ+SUMQ1
         SUMCT=SUMCT+SUMCT1
         SUMCTU=SUMCTU+SUMCTU1	
         SUMCTd=SUMCTd+SUMCTd1		
         SUMdA=SUMdA+SUMdA1	
         SUMFx=SUMFx+SUMFx1			 
         Power=Power+Power1	
         SUMFr=SUMFr+SUMFr1
         SUMFt=SUMFt+SUMFt1		
         SUMRe=SUMRe+SUMRe1	
         SUMat=SUMat+SUMat1			 
         enddo		 
      ENDIF

	  
      if(me.eq.0)then
      if(t.eq.1)then
      Open (unit=1097,file='output/check_a_CT_vawt.out',
     +      position='append')
      end if
c      write(1097,5111) Power/1000.d0,SUMQ,SUMCTU,SUMCT,SUMCTd,SUMdA,
c     +                 SUMFx,wtomega_vawt  
      write(1097,5111) Power/1000.d0,SUMCTU*1.0d0,SUMCT,SUMCTd,
     +                 wtomega_vawt,SUMdA
c     +   ,SUMFx/(0.5*SUMdA*SUMCTU**2.0)
     +    ,SUMFx/(0.5*SUMdA*7.04**2.0)/1.0d0
c     +   ,Power/(0.5*SUMdA*SUMCTU**3.0)
     +    ,Power/(0.5*SUMdA*7.04**3.0)/1.0d0
     +    ,SUMFr/(0.5*SUMdA*7.04**2.0)/1.0d0
     +    ,SUMFt/(0.5*SUMdA*7.04**2.0)/1.0d0
     +    ,SUMRe/10**6,SUMat  	 
cc      write(*,*) SUMdA/(2*pi*wtR(1)*dz)/num_turbine
      end if	  

      if(me.eq.5)then
      if(t.eq.1)then
      Open (unit=1096,file='output/AoA_vawt.out',position='append')
      end if
      write(1096,5111) AoA_new(1,1,2),AoA_rate(1,1,2),AoA_ds(1,1,2),
     +                 CL_new(1,1,2),CD_new(1,1,2) 	  
cc      write(*,*) SUMdA/(2*pi*wtR(1)*dz)/num_turbine
      end if


	  
 5111 format (60(1x,f11.5)) 
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine CALVAWT_ASSM(u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'dimen.h'
      integer*4 ii,i,j,k,t,M,blade,N,me_k

      real*8, dimension(Nx,Ny,Nz2)::u,v,w,fx,fy,fz,fax,faz,AveWt
!      real*8, dimension(num_turbine_vawt,Nx,Ny,Nz2):: CTy,CTx
!      real*8, dimension(num_turbine_vawt,Nx,Ny,Nz2):: CTzu,CTzw,CTru,CTrw	  
!      real*8, dimension(num_turbine_vawt,Ny,Nz2):: CTangx,CTangz
!      real*8, dimension(num_turbine_vawt,Nx):: ker_x
      real*8 tmp1,tmp2,tmp3

!      real*8, dimension(num_turbine_vawt,N_ang_vawt,N_blade)::
!     +uframe,vframe,wframe,meframe,u1frame,v1frame,w1frame,me1frame,
!     +uframed,u1framed	 
      
!      real*8, dimension(num_turbine_vawt,N_ang_vawt,Nx,Ny,Nz2):: 
!     +        ker_u_vawt,ker_w_vawt

!!      real*8, dimension(num_turbine_vawt,N_ang,Nx,Ny):: ker_u,ker_w	  

      real*8, dimension(num_turbine_vawt)::SUMCT,SUMa1,SUMdA,SUMFx,
     +       SUMQ,SUMb1,SUMCTd,SUMQ1,SUMCTU,Power,Power1,
     +	     SUMCT1,SUMCTU1,SUMCTd1,SUMdA1,SUMFx1,
     +       SUMFr,SumFr1,SUMFt,SumFt1,Sumat,Sumat1,
     +       SumRe,SumRe1	 

      real*8 a1,b1,at1,Tang(N_ang_vawt),blade_jj,blade_kk
      real*8 mu_y,mu_x,sigma,dr,omega,omega1,r,CTphi
      real*8 tmpy1,tmpz1,tmpy2,tmpz2,tmp
      real*8 ERRd1,ERRd2,tmpd1,tmpd2
      real*8 sum_ker_u(Nz2),sum_ker_w,max_ker_u,max_ker_w
      real*8 CTU,CTU1,CTV,CTW,chordl,CTsoli
      real*8 CTV1,CTV2,CTUrel,CTFx,CTFt,CTFr,CTFy
      real*8 CTCL,CTCD,dang,Re,dA,RA,dL
      real*8 umean(num_turbine_vawt),umeand(num_turbine_vawt)	 
      real*8 AoA_old(num_turbine_vawt,3,Nz2),
     +       AoA_new(num_turbine_vawt,3,Nz2)
      real*8 AoA_rate(num_turbine_vawt,3,Nz2),
     +       AoA_ds(num_turbine_vawt,3,Nz2)	
      real*8 CL_new(num_turbine_vawt,3,Nz2),
     +       CD_new(num_turbine_vawt,3,Nz2)	  
      real*8 lambdar,at,CL,CD,phi,CL_max,beta,at_ss,CL_ss 	  
      real*8, external :: CALCL_vawt,CALCD_vawt	  

      real*8, external :: chord,calang
      real*8, external :: itp2D,calU

      real*8 rr(4),yy1,xx1,rr1(4),wp(4),atd,qr,q,Re_c,num	
      integer*4 n1y,n2y,n1x,n2x,ii1,i1,DS	  
	  
!      save CTy,CTzu,CTzw,CTangx,CTangz,ker_x,ker_u_vawt,ker_w_vawt
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      if (t.eq.1) then
!        do ii=1,num_turbine_vawt
!        do k=2,nzb+1
!        do j=1,ny
!        do i=1,nx		
		
!          CTy(ii,i,j,k)    = (wty_vawt(ii)-j)*dy+0.0d-50 
!          CTx(ii,i,j,k)    = (wtx_vawt(ii)-i)*dx+0.0d-50		  
!          CTru(ii,i,j,k)   = sqrt(CTy(ii,i,j,k)**2.+CTx(ii,i,j,k)**2.)

!        end do
!        end do
!        end do
!        end do		
!      end if

      if (t.eq.1) then
	 	  
      dang = 2.d0*pi/dble(N_ang_vawt)
      do i=1,N_ang_vawt		
      Tang(i) = (dble(i)-1)*dang
      end do

!      ker_u_vawt=0.0d0	  
!      do ii=1,num_turbine_vawt  	  
!        do blade=1,N_ang_vawt

!          sum_ker_u  = 0.0
!          mu_x = wtr_vawt(ii)*dcos(Tang(blade))  
!          mu_y = wtr_vawt(ii)*dsin(Tang(blade)) 

!          do k=2,nzb+1
!          do j=1,ny
!          do i=1,nx		  
!          if (CTru(ii,i,j,k).le.2.0*wtr_vawt(ii))then
!          tmp3 = (CTy(ii,i,j,k)-mu_y)**2.+(CTx(ii,i,j,k)-mu_x)**2.
!          ker_u_vawt(ii,blade,i,j,k) = 1/(((2.5*dx)**3)*(pi**1.5))
!     +                           *dexp(-1.0*tmp3/((2.5*dx)**2.0))
!          sum_ker_u(k) = sum_ker_u(k)+ker_u_vawt(ii,blade,i,j,k)
!          else
!          ker_u_vawt(ii,blade,i,j,k) = 0.d0
!          end if
!          end do
!          end do
!          end do
	  	  
!          do k=2,Nzb+1
!          do j=1,Ny
!          do i=1,Nx		  
!          ker_u_vawt(ii,blade,i,j,k) = ker_u_vawt(ii,blade,i,j,k)
!     +                                 /sum_ker_u(k)
!          end do
!          end do
!          end do

!        end do
!      end do
      end if

	  
!      do k=2,nzb+1
!      do j=1,Ny	 
!      do i=5,5	  
!      dA=dz*dy 
!c      n1x=i
!      CTU=u(i,j,k)			
!      CTV=v(i,j,k)
	         
!cc      CALL RANDOM_NUMBER(num)
!c      print*, 'num', num	  
!cc      CTFx=0.5*num*4.0d0*(CTU**2)*dA/(dx*dy*dz) 
!cc      Fx(i,j,k)=Fx(i,j,k)+CTFx	

!      end do
!      end do	
!      end do	  
  
      dang = 2.d0*pi/dble(N_ang_vawt)

      SUMdA=0.0
      SUMa1=0.0
      SUMb1=0.0	  
      SUMCT=0.0
      SUMCTU=0.0	  
      SUMCTd=0.0	  
      SUMFx=0.0	
      SUMQ=0.0
      Power=0.0	  
CCCCCCCCCCCCCCCCCCCCCCCCCCC
      do ii=1,num_turbine_vawt
       do k=2,nzb+1
	  
        if ((me*nzb+k-1.5)*dz.gt.(Zhub_vawt(ii)-0.5*wtH_vawt(ii)) .and.
     +     (me*nzb+k-1.5)*dz.lt.(Zhub_vawt(ii)+0.5*wtH_vawt(ii))) then
	                                                      	  
        do blade=1,N_ang_vawt

         xx1= wtx_vawt(ii)-(wtr_vawt(ii)*dcos(Tang(blade))/dx)  
         yy1= wty_vawt(ii)-(wtr_vawt(ii)*dsin(Tang(blade))/dy) 
         n1x=floor(xx1)
         n2x=n1x+1
         n1y=floor(yy1)
         n2y=n1y+1
         rr1(1)=(((yy1-n1y)*dy)**2+((xx1-n1x)*dx)**2)**0.5
         rr1(2)=(((yy1-n1y)*dy)**2+((xx1-n2x)*dx)**2)**0.5
         rr1(3)=(((yy1-n2y)*dy)**2+((xx1-n1x)*dx)**2)**0.5
         rr1(4)=(((yy1-n2y)*dy)**2+((xx1-n2x)*dx)**2)**0.5
         ii1=0
         do i1=1,4
          if (rr1(i1).le.1.0d-6) then
           wp=0.0
           wp(i1)=1.0
           ii1=1 
          endif
         enddo
         if (ii1==0) then
          do i1=1,4	  
          wp(i1)=(1/rr1(i1))/(1/rr1(1)+1/rr1(2)+1/rr1(3)+1/rr1(4))		 
          enddo
         endif	        
         CTU=wp(1)*u(n1x,n1y,k)+wp(2)*u(n2x,n1y,k)+
     +       wp(3)*u(n1x,n2y,k)+wp(4)*u(n2x,n2y,k)			
         CTV=wp(1)*v(n1x,n1y,k)+wp(2)*v(n2x,n1y,k)+
     +       wp(3)*v(n1x,n2y,k)+wp(4)*v(n2x,n2y,k)		
         
!            CTU=0.93d0*CTU 
!            CTV=0.93d0*CTV
                       
         dA      = wtr_vawt(ii)*dang*dz
         SUMdA(ii) =SUMdA(ii)+dA*(z_i**2.)	
         SUMCTU(ii)=SUMCTU(ii)+CTU*dA/(2.0*pi*wtR_vawt(ii)*wtH_vawt(ii))	  
         
        end do
	  
        endif

       end do
       IF (me>0) then
         call MPI_SEND(SUMCTU(ii),1,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	 
       ELSE
         do i=1,nprocs-1	  
         call MPI_RECV(SUMCTU1(ii),1,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	 
         SUMCTU(ii)=SUMCTU(ii)+SUMCTU1(ii)	
         enddo		 
       ENDIF

       IF (me.eq.0) then
!c         wtomega_vawt(ii)=33.5d0
       wtomega_vawt(ii)=1.2d0*SUMCTU(ii)*TSR_vawt(ii)/
     +                  (wtR_vawt(ii)*z_i)*60.d0/(2.0*pi)
!       print*, wtomega_vawt(ii) 	 
         do i=1,nprocs-1
          call MPI_SEND(wtomega_vawt(ii),1, MPI_DOUBLE_PRECISION,
     +           i,i,nall,ierr ) 
         enddo
       ENDIF
       if (me.ne.0) then	  
       call MPI_RECV(wtomega_vawt(ii),1, MPI_DOUBLE_PRECISION,
     +                  0,me,nall,status2,ierr )
       endif	 
 		
      end do 	  
CCCCCCCCCCCCCCCCCCCCCCCCCCC	  
      SUMdA=0.0
      SUMa1=0.0
      SUMb1=0.0	  
      SUMCT=0.0
      SUMCTU=0.0	  
      SUMCTd=0.0	  
      SUMFx=0.0	
      SUMQ=0.0
      SUMFr=0.0
      SUMFt=0.0	  
      Power=0.0
      SumRe=0.0
      Sumat=0.0 	  
CCCCCCCCCCCCCCCCCCCCCCCCCCC
	  
      if (t.eq.1) AoA_old=0.0 	  
  
!cc      if (me.ge.7 .and. me.le.8) then
      do ii=1,num_turbine_vawt
      do k=2,nzb+1
	  
       if ((me*nzb+k-1.5)*dz.gt.(Zhub_vawt(ii)-0.5*wtH_vawt(ii)) .and.
     +    (me*nzb+k-1.5)*dz.lt.(Zhub_vawt(ii)+0.5*wtH_vawt(ii))) then
	                                                      	  
       do blade=1,N_ang_vawt
!      do M=1,N_blade_vawt

!        if (dir_vawt(ii).gt.0.d0) then 
!        blade=t+(M-1)*N_ang_vawt/N_blade_vawt
!        blade=mod(blade-1,N_ang_vawt)+1
!        else		
!        blade=-t+(M-1)*N_ang_vawt/N_blade_vawt+
!     +         ((t/N_ang_vawt)+1)*N_ang_vawt
!        blade=mod(blade-1,N_ang_vawt)+1		
!        endif 
		
        xx1= wtx_vawt(ii)-(wtr_vawt(ii)*dcos(Tang(blade))/dx)  
        yy1= wty_vawt(ii)-(wtr_vawt(ii)*dsin(Tang(blade))/dy) 
        n1x=floor(xx1)
        n2x=n1x+1
        n1y=floor(yy1)
        n2y=n1y+1
        rr1(1)=(((yy1-n1y)*dy)**2+((xx1-n1x)*dx)**2)**0.5
        rr1(2)=(((yy1-n1y)*dy)**2+((xx1-n2x)*dx)**2)**0.5
        rr1(3)=(((yy1-n2y)*dy)**2+((xx1-n1x)*dx)**2)**0.5
        rr1(4)=(((yy1-n2y)*dy)**2+((xx1-n2x)*dx)**2)**0.5
        ii1=0
        do i1=1,4
        if (rr1(i1).le.1.0d-6) then
        wp=0.0
        wp(i1)=1.0
        ii1=1 
        endif
        enddo
        if (ii1==0) then
         do i1=1,4	  
         wp(i1)=(1/rr1(i1))/(1/rr1(1)+1/rr1(2)+1/rr1(3)+1/rr1(4))		 
         enddo
        endif	        
        CTU=wp(1)*u(n1x,n1y,k)+wp(2)*u(n2x,n1y,k)+
     +      wp(3)*u(n1x,n2y,k)+wp(4)*u(n2x,n2y,k)			
        CTV=wp(1)*v(n1x,n1y,k)+wp(2)*v(n2x,n1y,k)+
     +      wp(3)*v(n1x,n2y,k)+wp(4)*v(n2x,n2y,k)
	 
!       	CTU=0.93d0*CTU 
!       	CTV=0.93d0*CTV
		
!      chordl  = chord(r*z_i)/z_i
       chordl  = chord_vawt(ii)	  
       CTsoli  = N_blade_vawt*chordl/(2.d0*pi*wtr_vawt(ii))

       omega1  = wtomega_vawt(ii)*(2.d0*pi/60)*(z_i/u_star)
!       lambdar=(wtr_vawt(ii)*omega1-dir_vawt(ii)*CTU*dsin(Tang(blade))+
!     +         dir_vawt(ii)*CTV*dcos(Tang(blade)))/
!     +        (CTU*dcos(Tang(blade))+CTV*dsin(Tang(blade)))
       CTphi = datan((CTU*dcos(Tang(blade))+1.0*CTV*dsin(Tang(blade)))/
     +	  (wtr_vawt(ii)*omega1-dir_vawt(ii)*CTU*dsin(Tang(blade))+
     +         dir_vawt(ii)*1.0*CTV*dcos(Tang(blade))))	

!       CTphi=0.9*CTphi	 
       at  = CTphi

       CTCL   = 1.0*CALCL_vawt(at*180.d0/pi,wtr_vawt(ii)*z_i)
       CTCD   = 1.0*CALCD_vawt(at*180.d0/pi,wtr_vawt(ii)*z_i)

       CTUrel=sqrt((CTU*dcos(Tang(blade))+1.0*CTV*dsin(Tang(blade)))**2.
     +       +(wtr_vawt(ii)*omega1-dir_vawt(ii)*CTU*dsin(Tang(blade))+
     +         dir_vawt(ii)*1.0*CTV*dcos(Tang(blade)))**2.)
 
       Re_c=chordl*CTUrel/nu  
  
       qr=(abs((me*nzb+k-1.5)*dz-Zhub_vawt(ii)))/
     +    (0.5*wtH_vawt(ii))+1.0d-50
!!       q=2./pi*dacos(
!!     +   dexp(-(N_blade_vawt/2.*(1.-qr))/(qr*abs(dsin(CTphi)))))
       q=2./pi*dacos(
     +   dexp(-(N_blade_vawt/2.*(1.-qr))/(abs(dsin(CTphi)))))	 
!       q=1.0d0
!       if (M.eq.1) print*, q,me
CCCCCCCCCCCCCCCCCCCCCCCCCCC DYNAMIC STALL
      beta=0.077d0	   
      at_ss=14.0d0
      CL_ss=1.08d0
CCCCCC Dynamic Stall CCCCCC  	  
!       AoA_new(ii,M,k)=at*180.d0/pi 
!       AoA_rate(ii,M,k)=(abs(AoA_new(ii,M,k))-abs(AoA_old(ii,M,k)))
!     + 	                *pi/180.d0/dt*chordl/(2.0*CTUrel)	 
!       if (t.eq.1) AoA_rate(ii,M,k)=0.0 
!
!       AoA_ds(ii,M,k)=0.0d0
!       if (AoA_rate(ii,M,k).gt.0.0d0) then
!          if (AoA_new(ii,M,k).gt.0.0d0) then
!       AoA_ds(ii,M,k)=at_ss+(1.0*(AoA_rate(ii,M,k))**0.5)*180.d0/pi
!          else if (AoA_new(ii,M,k).lt.0.0d0) then
!       AoA_ds(ii,M,k)=-at_ss-(1.0*(AoA_rate(ii,M,k))**0.5)*180.d0/pi		  
!          endif
!       endif

!       DS=1
!       if (DS.eq.1) then 
!       if (AoA_rate(ii,M,k).gt.0.0d0) then
!        if (AoA_new(ii,M,k).gt.at_ss .and. 
!     +       AoA_new(ii,M,k).lt.AoA_ds(ii,M,k)) then	   
!c           	 CTCL=CL_ss+beta*(AoA_new(ii,M,k)-at_ss)
!c             CL_max=CL_ss+2.0d0
!           	 CTCL=CL_ss/(at_ss*pi/180.d0)*dsin(at)
!             CL_max=CL_ss+2.0d0
!         if (40.*2.*AoA_rate(ii,M,k).lt.2.0) then
!c       if ((CL_ss+40.*2.*AoA_rate(ii,M,k)).gt.
!c     +     beta*(AoA_ds(ii,M,k)-at_ss)) then   
!                CL_max=CL_ss+40.0*2.*AoA_rate(ii,M,k) 
!c             else				
!c                CL_max=CL_ss+beta*(AoA_ds(ii,M,k)-at_ss)				
!c       endif
!         endif			 
!             if (CTCL .ge. CL_max) CTCL=CL_max         			 
!        endif
!       endif		 
!       if (AoA_rate(ii,M,k).gt.0.0d0) then
!         if (AoA_new(ii,M,k).gt.at_ss .and.
!     +       AoA_new(ii,M,k).gt.AoA_ds(ii,M,k)) then
!c           	 CTCL=CL_ss+beta*(AoA_new(ii,M,k)-at_ss)
!c             CL_max=CL_ss+2.0d0
!           	 CTCL=CL_ss/(at_ss*pi/180.d0)*dsin(at)
!             CL_max=CL_ss+2.0d0
!        if (40.*2.*AoA_rate(ii,M,k).lt.2.0) then
!c        if (40.*2.*AoA_rate(ii,M,k).gt.beta*(AoA_ds(ii,M,k)-at_ss)) then   
!                CL_max=CL_ss+40.0*2.*AoA_rate(ii,M,k) 
!c             else				
!c                CL_max=CL_ss+beta*(AoA_ds(ii,M,k)-at_ss)				
!c        endif
!        endif			
!             if (CTCL .ge. CL_max) CTCL=CL_max
!c             CTCD=abs(CTCL*dsin(AoA_new(ii,M,k)*pi/180.0d0))      
!             CTCD=abs(CTCL*dtan(AoA_new(ii,M,k)*pi/180.0d0))  			 
!         endif
!       endif
!	   
!       if (AoA_rate(ii,M,k).gt.0.0d0) then
!         if (AoA_new(ii,M,k).lt.-at_ss .and. 
!     +       AoA_new(ii,M,k).gt.AoA_ds(ii,M,k)) then	   
!c          	 CTCL=-CL_ss-beta*(-AoA_new(ii,M,k)-at_ss)
!c             CL_max=-CL_ss-2.0d0
!           	 CTCL=CL_ss/(at_ss*pi/180.d0)*dsin(at)
!             CL_max=-CL_ss-2.0d0			 
!			 
!       if (40.*2.*AoA_rate(ii,M,k).lt.2.0) then
!c       if (40.*2.*AoA_rate(ii,M,k).gt.beta*(-AoA_ds(ii,M,k)-at_ss)) then    
!                CL_max=-CL_ss-40.0*2.*AoA_rate(ii,M,k) 
!c             else				
!c                CL_max=-CL_ss-beta*(-AoA_ds(ii,M,k)-at_ss)				
!c       endif
!       endif			
!             if (CTCL .le. CL_max) CTCL=CL_max         			 
!         endif
!       endif		 
!       if (AoA_rate(ii,M,k).gt.0.0d0) then
!         if (AoA_new(ii,M,k).lt.-at_ss .and.
!     +       AoA_new(ii,M,k).lt.AoA_ds(ii,M,k)) then
!c          	 CTCL=-CL_ss-beta*(-AoA_new(ii,M,k)-at_ss)
!c             CL_max=-CL_ss-2.0d0
!           	 CTCL=CL_ss/(at_ss*pi/180.d0)*dsin(at)
!             CL_max=-CL_ss-2.0d0	
			 
!       if (40.*2.*AoA_rate(ii,M,k).lt.2.0) then
!c       if (40.*2.*AoA_rate(ii,M,k).gt.beta*(-AoA_ds(ii,M,k)-at_ss)) then   
!                CL_max=-CL_ss-40.0*2.*AoA_rate(ii,M,k) 
!c             else				
!c                CL_max=-CL_ss-beta*(-AoA_ds(ii,M,k)-at_ss)				
!c       endif
!       endif					
!             if (CTCL .le. CL_max) CTCL=CL_max
!c             CTCD=abs(CTCL*dsin(AoA_new(ii,M,k)*pi/180.0d0))   
!             CTCD=abs(CTCL*dtan(AoA_new(ii,M,k)*pi/180.0d0))           			 
!         endif
!       endif
!
!       CL_new(ii,M,k)=CTCL
!       CD_new(ii,M,k)=CTCD	   
!	   
!       endif	
!       AoA_old(ii,M,k)=AoA_new(ii,M,k)	   
CCCCCCCCCCCCCCCCCCCCCCCCCCC END OF DYNAMIC STALL
CCCCCC NEW Dynamic Stall CCCCCC  	  
!       DS=0
!       if (DS.eq.10) then 
!       if (AoA_rate(ii,M,k).gt.0.0d0) then
!         if (AoA_new(ii,M,k).gt.0. .and. 
!     +       AoA_new(ii,M,k).lt.30.) then	   
!           	 CTCL=AoA_new(ii,M,k)*0.073d0		 
!         endif
!       endif
!       if (AoA_rate(ii,M,k).lt.0.0d0) then
!         if (AoA_new(ii,M,k).gt.0. .and. 
!     +       AoA_new(ii,M,k).lt.30.) then	   
!c           	 CTCL=AoA_new(ii,M,k)*0.02d0		 
!         endif
!       endif
!       if (AoA_rate(ii,M,k).gt.0.0d0) then
!         if (AoA_new(ii,M,k).lt.(-0.) .and. 
!     +       AoA_new(ii,M,k).gt.(-20.)) then	   
!           	 CTCL=AoA_new(ii,M,k)*0.073d0		 
!         endif
!       endif
!       if (AoA_rate(ii,M,k).lt.0.0d0) then
!         if (AoA_new(ii,M,k).lt.(-0.) .and. 
!     +       AoA_new(ii,M,k).gt.(-30.)) then	   
!c           	 CTCL=AoA_new(ii,M,k)*0.02d0		 
!         endif
!       endif

!       CL_new(ii,M,k)=CTCL
!       CD_new(ii,M,k)=CTCD	   
	   
!       endif	      
CCCCCCCCCCCCCCCCCCCCCCCCCCC END OF DYNAMIC STALL
!      dA      = dz*chordl
!      CTFr    = q*0.5*(CTUrel**2)*dA/(dx*dy*dz)*
!     +          (+CTCL*dcos(CTphi)+CTCD*dsin(CTphi))
!      CTFt    = q*dir_vawt(ii)*0.5*(CTUrel**2)*dA/(dx*dy*dz)*
!     +          (+CTCL*dsin(CTphi)-CTCD*dcos(CTphi))	
  
!cc     +      dabs(+CTCL*dsin(CTphi)-CTCD*dcos(CTphi))

       dA      = wtr_vawt(ii)*dang*dz
       CTFr    = 0.5*(CTUrel**2)*CTsoli*dA/(dx*dy*dz)*
     +          (+CTCL*dcos(CTphi)+CTCD*dsin(CTphi))
       CTFt    = dir_vawt(ii)*0.5*(CTUrel**2)*CTsoli*dA/(dx*dy*dz)*
     +          (+CTCL*dsin(CTphi)-CTCD*dcos(CTphi))	 

       CTFx=CTFr*dcos(Tang(blade))+CTFt*dsin(Tang(blade))     	 
       CTFy=CTFr*dsin(Tang(blade))-CTFt*dcos(Tang(blade)) 

      SUMdA(ii) =SUMdA(ii)+dA*(z_i**2.)
      if (CTU.ne.0.0d0) then	  
      SUMCT(ii) =SUMCT(ii)+CTFx/(0.5d0*((CTU/u_star)**2.)*	      
     +           dA/(dx*dy*dz))*dA/(2.0*pi*wtR_vawt(ii)*wtH_vawt(ii))
      else
      SUMCT(ii) =SUMCT(ii)+0.0d0
      endif	  
      SUMCTd(ii) =SUMCTd(ii)+CTFx/(0.5d0*((8.0d0/u_star)**2.)*	      
     +           dA/(dx*dy*dz))*dA/(2.0*pi*wtR_vawt(ii)*wtH_vawt(ii))	 
      SUMCTU(ii)=SUMCTU(ii)+CTU*dA/(2.0*pi*wtR_vawt(ii)*wtH_vawt(ii))	  

      SUMQ(ii)=SUMQ(ii)+dir_vawt(ii)*CTFt*(dx*dy*dz)*wtr_vawt(ii)*
     +         (z_i**3)		

      SUMFx(ii)=SUMFx(ii)+CTFx*(dx*dy*dz)*(z_i**2.)	 	 

!      if (M.eq.1) then 	   
!      SUMFr(ii)=SUMFr(ii)+CTFr*(dx*dy*dz)*(z_i**2)	
!      SUMFt(ii)=SUMFt(ii)+CTFt*(dx*dy*dz)*(z_i**2)	
!      Sumat(ii)=Sumat(ii)+
!     +          at*180.0/pi*dA/(chordl*wtH_vawt(ii))
!      SumRe(ii)=SumRe(ii)+Re_c*dA/(chordl*wtH_vawt(ii))  
!      endif 
	  
      Fx (n1x,n1y,k)=Fx (n1x,n1y,k) + CTFx*wp(1)	
      Fx (n2x,n1y,k)=Fx (n2x,n1y,k) + CTFx*wp(2)	
      Fx (n1x,n2y,k)=Fx (n1x,n2y,k) + CTFx*wp(3)	
      Fx (n2x,n2y,k)=Fx (n2x,n2y,k) + CTFx*wp(4)			 
      Fy (n1x,n1y,k)=Fy (n1x,n1y,k) + CTFy*wp(1)	
      Fy (n2x,n1y,k)=Fy (n2x,n1y,k) + CTFy*wp(2)	
      Fy (n1x,n2y,k)=Fy (n1x,n2y,k) + CTFy*wp(3)	
      Fy (n2x,n2y,k)=Fy (n2x,n2y,k) + CTFy*wp(4)	

!      do j=1,ny
!      do i=1,nx
!         if (ker_u_vawt(ii,blade,i,j,k).ge.1E-15)then
!         Fx (i,j,k)=Fx (i,j,k) + 1.*CTFx  * ker_u_vawt(ii,blade,i,j,k)
!         Fy (i,j,k)=Fy (i,j,k) + 1.*CTFy  * ker_u_vawt(ii,blade,i,j,k)
!        end if
!      end do
!      end do
   
      end do
	  
      endif

      end do

      Power(ii)=SUMQ(ii)*wtomega_vawt(ii)*(2.d0*pi/60)	
	  
      end do 	  

!      end if

      IF (me>0) then
         call MPI_SEND(SUMQ,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(SUMCT,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(SUMCTU,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	 
         call MPI_SEND(SUMCTd,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	
         call MPI_SEND(SUMdA,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	
         call MPI_SEND(SUMFx,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	
         call MPI_SEND(Power,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	 
         call MPI_SEND(SUMFr,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	
         call MPI_SEND(SUMFt,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(SUMRe,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	
         call MPI_SEND(SUMat,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	 
      ELSE
         do i=1,nprocs-1	  
         call MPI_RECV(SUMQ1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         call MPI_RECV(SUMCT1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         call MPI_RECV(SUMCTU1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	 
         call MPI_RECV(SUMCTd1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	
         call MPI_RECV(SUMdA1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	 
         call MPI_RECV(SUMFx1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	 
         call MPI_RECV(Power1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	
         call MPI_RECV(SUMFr1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	
         call MPI_RECV(SUMFt1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	
         call MPI_RECV(SUMRe1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	
         call MPI_RECV(SUMat1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )		 
         SUMQ=SUMQ+SUMQ1
         SUMCT=SUMCT+SUMCT1
         SUMCTU=SUMCTU+SUMCTU1	
         SUMCTd=SUMCTd+SUMCTd1		
         SUMdA=SUMdA+SUMdA1	
         SUMFx=SUMFx+SUMFx1			 
         Power=Power+Power1	
         SUMFr=SUMFr+SUMFr1
         SUMFt=SUMFt+SUMFt1		
         SUMRe=SUMRe+SUMRe1	
         SUMat=SUMat+SUMat1			 
         enddo		 
      ENDIF

	  
      if(me.eq.0)then
      if(t.eq.1)then
      Open (unit=1097,file='output/check_a_CT_vawt.out',
     +      position='append')
      end if
	  
      if (mod(t,c_count)==0) then	

      write(1528,5111) Power/(10**3),SUMCTU
      write(1529,5111) wtomega_vawt,SUMCT
	  
      write(*,5111) Power(1)/1000.d0,SUMCTU(1),SUMCT(1),
     +              SUMdA(1),wtomega_vawt(1)
  
      write(1097,5111) SUMdA
	 
!     +   ,SUMFx/(0.5*SUMdA*SUMCTU**2.0)
!     +    ,SUMFx/(0.5*SUMdA*8**2.0)/1.0d0
!     +   ,Power/(0.5*SUMdA*SUMCTU**3.0)
!     +    ,Power/(0.5*SUMdA*8**3.0)/1.0d0
!     +    ,SUMFr/(0.5*SUMdA*8**2.0)/1.0d0
!     +    ,SUMFt/(0.5*SUMdA*8**2.0)/1.0d0
!     +    ,SUMRe/10**6,SUMat  	
	 
!      write(*,*) SUMdA/(2*pi*wtR(1)*dz)/num_turbine
      end if	  

!      if(me.eq.5)then
!      if(t.eq.1)then
!      Open (unit=1096,file='output/AoA_vawt.out',position='append')
!      end if
!      write(1096,5111) AoA_new(1,1,2),AoA_rate(1,1,2),AoA_ds(1,1,2),
!     +                 CL_new(1,1,2),CD_new(1,1,2) 	  
!cc      write(*,*) SUMdA/(2*pi*wtR(1)*dz)/num_turbine
!      end if
       endif

	  
 5111 format (60(1x,f11.5)) 
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine CALVAWT_FP(u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'dimen.h'
      integer*4 ii,i,j,k,t,M,blade,N,me_k

      real*8, dimension(Nx,Ny,Nz2)::u,v,w,fx,fy,fz,fax,faz,AveWt
      real*8, dimension(num_turbine_vawt,Nx,Ny,Nz2):: CTy,CTx
      real*8, dimension(num_turbine_vawt,Nx,Ny,Nz2):: CTzu,CTzw,CTru,CTrw	  
      real*8, dimension(num_turbine_vawt,Ny,Nz2):: CTangx,CTangz
      real*8, dimension(num_turbine_vawt,Nx):: ker_x
      real*8 tmp1,tmp2,tmp3

      real*8, dimension(num_turbine_vawt,N_ang_vawt,N_blade)::
     +uframe,vframe,wframe,meframe,u1frame,v1frame,w1frame,me1frame,
     +uframed,u1framed	 
      
      real*8, dimension(num_turbine_vawt,N_ang_vawt,Nx,Ny,Nz2):: 
     +        ker_u_vawt,ker_w_vawt
cc      real*8, dimension(num_turbine_vawt,N_ang,Nx,Ny):: ker_u,ker_w	  

      real*8, dimension(num_turbine_vawt)::SUMCT,SUMa1,SUMdA,SUMFx,
     +       SUMQ,SUMb1,SUMCTd,SUMQ1,SUMCTU,Power,Power1,
     +	     SUMCT1,SUMCTU1,SUMCTd1,SUMdA1,SUMFx1,
     +       SUMFr,SumFr1,SUMFt,SumFt1,Sumat,Sumat1,
     +       SumRe,SumRe1,SumCTUd,SumCTUd1	 

      real*8 a1,b1,at1,Tang(N_ang_vawt),blade_jj,blade_kk
      real*8 mu_y,mu_x,sigma,dr,omega,omega1,r,CTphi
      real*8 tmpy1,tmpz1,tmpy2,tmpz2,tmp
      real*8 ERRd1,ERRd2,tmpd1,tmpd2
      real*8 sum_ker_u(Nz2),sum_ker_w,max_ker_u,max_ker_w
      real*8 CTU,CTU1,CTV,CTW,chordl,CTsoli,CTUd
      real*8 CTV1,CTV2,CTUrel,CTFx,CTFt,CTFr,CTFy
      real*8 CTCL,CTCD,dang,Re,dA,RA,dL
      real*8 umean(num_turbine_vawt),umeand(num_turbine_vawt)	 
      real*8 AoA_old(num_turbine_vawt,3,Nz2),
     +       AoA_new(num_turbine_vawt,3,Nz2)
      real*8 AoA_rate(num_turbine_vawt,3,Nz2),
     +       AoA_ds(num_turbine_vawt,3,Nz2)	
      real*8 CL_new(num_turbine_vawt,3,Nz2),
     +       CD_new(num_turbine_vawt,3,Nz2)	  
      real*8 lambdar,at,CL,CD,phi,CL_max,beta,at_ss,CL_ss 	  
      real*8, external :: CALCL_vawt,CALCD_vawt	  

      real*8, external :: chord,calang
      real*8, external :: itp2D,calU

      real*8 rr(4),yy1,xx1,rr1(4),wp(4),atd,qr,q,Re_c,num	
      integer*4 n1y,n2y,n1x,n2x,ii1,i1,DS	  

      real*8 CTT_x,CTT_y 	  
      save CTy,CTzu,CTzw,CTangx,CTangz,ker_x,ker_u_vawt,ker_w_vawt
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCC	  
      SUMdA=0.0
      SUMa1=0.0
      SUMb1=0.0	  
      SUMCT=0.0
      SUMCTU=0.0	  
      SUMCTd=0.0	  
      SUMFx=0.0	
      SUMQ=0.0
      SUMFr=0.0
      SUMFt=0.0	  
      Power=0.0
      SumRe=0.0
      Sumat=0.0 
      SumCTUd=0.0	  
CCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      do k=2,nzb+1
!      do j=1,Ny	 
!      do i=5,5	  
!      dA=dz*dy 
!c      n1x=i
!      CTU=u(i,j,k)			
!      CTV=v(i,j,k)
!	         
!c      CALL RANDOM_NUMBER(num)
!cc      print*, 'num', num	  
!c      CTFx=0.5*num*4.0d0*(CTU**2)*dA/(dx*dy*dz) 
!c      Fx(i,j,k)=Fx(i,j,k)+CTFx	
!      end do
!      end do	
!      end do	  
CCCCCCCCCCCCCCCCCCCCCCCCCCC
      do ii=1,num_turbine_vawt
	  
!      n1x=wtx_vawt(ii)	 
      n1x=wtx_vawt(ii)-3.*wtr_vawt(ii)/dx	  
	  
      do k=2,nzb+1
      if ((me*nzb+k-1.5)*dz.gt.(Zhub_vawt(ii)-0.5*wtH_vawt(ii)) .and.
     +    (me*nzb+k-1.5)*dz.lt.(Zhub_vawt(ii)+0.5*wtH_vawt(ii))) then

!      do i=1,5	 
      do j=1,Ny	 
      if (abs(j-wty_vawt(ii))*dy .lt. wtr_vawt(ii)) then
      if (j .lt. wty_vawt(ii)) then
          CTT_x=0.64d0
          CTT_y=0.0d0
      elseif (j .gt. wty_vawt(ii)) then 		  
          CTT_x=0.64d0
          CTT_y=0.0d0
      else
          CTT_x=0.64d0
          CTT_y=0.0d0
      endif

!      q=10.0*dcos(pi/2.*((j-wty_vawt(ii))*dy)/(wtr_vawt(ii)))/6.3138d0
      q=1.0d0
	  
      dA=dz*dy 
               
      n1x=4			   
      CTU=u(n1x,j,k)			
      CTV=v(n1x,j,k)
      CTUd=u(wtx_vawt(ii),j,k)
	        
!      if (mod(j+0*mod(t/100,2),2).eq.0) then 
!      CALL RANDOM_NUMBER(num)
!      print*, 'num', num	  
!!      CTFx=0.5*num*CTT_vawt(ii)*(CTU**2)*dA/(dx*dy*dz) 
!      else
      CTFx=0.5*CTT_x*q*(CTU**2)*dA/(dx*dy*dz)
      CTFy=0.0*CTT_y*q*(CTU**2)*dA/(dx*dy*dz)  	  
!      endif	  

      SUMdA(ii) =SUMdA(ii)+dA*(z_i**2.)	
      SUMCTU(ii) =SUMCTU(ii)+CTU*dA/(2.*wtr_vawt(ii)*wtH_vawt(ii))
      SUMCTUd(ii) =SUMCTUd(ii)+CTUd*dA/(2.*wtr_vawt(ii)*wtH_vawt(ii))	  
      Power(ii)=Power(ii)+CTFx*(dx*dy*dz)*CTUd*(z_i**2.)		
      SUMFx(ii)=SUMFx(ii)+CTFx*(dx*dy*dz)*(z_i**2.)	 	 
!      SUMFy(ii)=SUMFy(ii)+CTFx*(dx*dy*dz)*(z_i**2.)	
	  
      Fx(wtx_vawt(ii),j,k)=Fx(wtx_vawt(ii),j,k)+CTFx
      Fy(wtx_vawt(ii),j,k)=Fy(wtx_vawt(ii),j,k)+CTFy	  
!      Fx(i,j,k)=Fx(i,j,k)+CTFx	  
      endif	  
!      enddo
      enddo	  
	  
      endif
      end do
	  
      end do 	  

      IF (me>0) then
         call MPI_SEND(SUMCTU,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	 
         call MPI_SEND(SUMCTUd,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	
         call MPI_SEND(SUMdA,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	
         call MPI_SEND(SUMFx,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	
         call MPI_SEND(Power,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )	 
      ELSE
         do i=1,nprocs-1	  
         call MPI_RECV(SUMCTU1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	 
         call MPI_RECV(SUMCTUd1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	
         call MPI_RECV(SUMdA1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	 
         call MPI_RECV(SUMFx1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	 
         call MPI_RECV(Power1,num_turbine_vawt,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	
	 
         SUMCTU=SUMCTU+SUMCTU1	
         SUMCTUd=SUMCTUd+SUMCTUd1		
         SUMdA=SUMdA+SUMdA1	
         SUMFx=SUMFx+SUMFx1			 
         Power=Power+Power1	
         enddo		 
      ENDIF

	  
      if(me.eq.0)then
      if(t.eq.1)then
      Open (unit=1097,file='output/check_a_CT_vawt.out',
     +      position='append')
      end if
c      write(1097,5111) Power/1000.d0,SUMQ,SUMCTU,SUMCT,SUMCTd,SUMdA,
c     +                 SUMFx,wtomega_vawt  
      write(1097,5111) Power/1000.0d0,SUMCTU,SUMCTUd,
     +                 SUMdA
c     +   ,SUMFx/(0.5*SUMdA*SUMCTU**2.0)
     +   ,SUMFx/(0.5*SUMdA*7.04**2.0)
c     +   ,Power/(0.5*SUMdA*SUMCTU**3.0)
     +    ,Power/(0.5*SUMdA*7.04**3.0)  	 
cc      write(*,*) SUMdA/(2*pi*wtR(1)*dz)/num_turbine
      end if	  

      if(me.eq.40)then
      if(t.eq.1)then
      Open (unit=1096,file='output/AoA_vawt.out',position='append')
      end if
c      write(1096,5111) AoA_new(1,1,2),AoA_rate(1,1,2),AoA_ds(1,1,2),
c     +                 CL_new(1,1,2),CD_new(1,1,2) 	  
cc      write(*,*) SUMdA/(2*pi*wtR(1)*dz)/num_turbine
      end if


	  
 5111 format (60(1x,f11.5)) 
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function CALCL_vawt(x,r)
      include 'dimen.h'
!------------------
      real*8 x,r,CALCL_vawt
      integer*4  I,J
  
      I=floor(x)+181
      J=CEILING(x)+181
      CALCL_vawt=(x-floor(x))/1.D0*ALFA_CL_CD_vawt(2,J)+(CEILING(x)-x)/
     +      1.D0*ALFA_CL_CD_vawt(2,I)
!---------------------------------------------------------------------
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function CALCD_vawt(x,r)
      include 'dimen.h'
!------------------
      real*8 x,r,CALCD_vawt
      integer*4 I,J
   
      I=floor(x)+181
      J=CEILING(x)+181

      CALCD_vawt=(x-floor(x))/1.D0*ALFA_CL_CD_vawt(3,J)+(CEILING(x)-x)
     +      /1.D0*ALFA_CL_CD_vawt(3,I)
!------------------
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function CALCL(x,r)
      include 'dimen.h'
!------------------
      real*8 x,r,CALCL
      integer*4  I,J
  
      I=floor(x)+91
      J=CEILING(x)+91
      CALCL=(x-floor(x))/1.D0*ALFA_CL_CD(2,J)+(CEILING(x)-x)/
     +      1.D0*ALFA_CL_CD(2,I)
!---------------------------------------------------------------------
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function CALCD(x,r)
      include 'dimen.h'
!------------------
      real*8 x,r,CALCD
      integer*4 I,J
   
      I=floor(x)+91
      J=CEILING(x)+91

      CALCD=(x-floor(x))/1.D0*ALFA_CL_CD(3,J)+(CEILING(x)-x)
     +      /1.D0*ALFA_CL_CD(3,I)
!------------------
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine CALALM1(u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'dimen.h'
      integer*4 ii,jj,kk,i1,j1,k1,i,j,k,t,N
      integer*4 Hubk(num_turbine),Tkk1,Tkk2,mek,mek1,mek2,kk1,kk2

      real*8, dimension(Nx,Ny,Nz2)::u,v,w,fx,fy,fz,fax,faz,AveWt
      real*8, dimension(num_turbine,Ny,Nz2):: CTy,CTzu,CTzw,CTru,CTrw
      real*8, dimension(num_turbine,Ny,Nz2):: CTangx,CTangz,dru,drw
      real*8, dimension(num_turbine,Nx):: ker_x
      real*8 tmp1,tmp2,tmp3,Bang(num_turbine,1:3),omega(num_turbine)

      real*8, dimension(:,:,:), allocatable::
     +uframe,vframe,wframe
      real*8, dimension(:,:,:), allocatable::
     +tempu,tempv,tempw
      
      real*8, dimension(num_turbine,N_ang,N_blade,Ny,Nz2):: ker_u,ker_w

      real*8 a1,b1,at1,Tang(N_ang),blade_jj,blade_kk
      real*8 mu_y,mu_z,sigma,dr,omega1,r,CTphi
      real*8 tmpy1,tmpz1,tmpy2,tmpz2,tmp
      real*8 ERRd1,ERRd2,tmpd1,tmpd2
      real*8 sum_ker_u,sum_ker_w,max_ker_u,max_ker_w
      real*8 CTU,CTU1,CTV,CTW,chordl,CTsoli
      real*8 CTV1,CTV2,CTUrel,CTFx,CTFt
      real*8 CTCL,CTCD,dang,Re,dA,RA,dL

      real*8, external :: chord,calang
      real*8, external :: itp,calU

      save CTy,CTzu,CTzw,CTangx,CTangz,ker_x,dru,drw
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (t.eq.1) then
        do ii=1,num_turbine
        do k=2,nzb+1
        do j=1,ny
          CTy(ii,j,k)    = (wty(ii)-j)*dy+0.0d-50 

          CTzu(ii,j,k)   = (me*nzb+k-1.5  )*dz-Zhub(ii)
          CTzw(ii,j,k)   = CTzu(ii,j,k)-0.5d0*dz

          CTru(ii,j,k)   = sqrt(CTy(ii,j,k)**2.+CTzu(ii,j,k)**2.)
          CTrw(ii,j,k)   = sqrt(CTy(ii,j,k)**2.+CTzw(ii,j,k)**2.)

          CTangx(ii,j,k) = calang(CTy(ii,j,k),CTzu(ii,j,k))
          CTangz(ii,j,k) = calang(CTy(ii,j,k),CTzw(ii,j,k))

          dru(ii,j,k)    = dsqrt((dy*dcos(CTangx(ii,j,k)))**2.
     +                   +       (dz*dsin(CTangx(ii,j,k)))**2.)
          drw(ii,j,k)    = dsqrt((dy*dcos(CTangz(ii,j,k)))**2.
     +                   +       (dz*dsin(CTangz(ii,j,k)))**2.)
        end do
        end do
        end do

        sigma = 1.d0
        tmp1  = 1/dsqrt(2.d0*pi*sigma*sigma)
        tmp2  = -1/(2.d0*sigma*sigma)
        do ii=1,num_turbine
        do i=1,Nx
          ker_x(ii,i) = tmp1*dexp(tmp2*(i-wtx(ii))**2.)
        end do
        tmp1 = sum(ker_x(ii,1:Nx))
        do i=1,Nx
          ker_x(ii,i) = ker_x(ii,i)/tmp1
        end do
        end do

      end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      i1=num_turbine        
      tmp1=wtr(1)/dy
      j1=nint(tmp1)+2
      tmp1=wtr(1)/dz
      k1=nint(tmp1)+2

      do ii=1,num_turbine
      Hubk(ii)=dble(Zhub(ii)/dz+0.5d0)
      end do

      kk=(int((Hubk(1)+k1)/Nzb)+1)*Nzb
      mek=kk/Nzb-1

      if(t.eq.1)then
      allocate(uframe(2*j1+1,kk,i1))
      allocate(vframe(2*j1+1,kk,i1))
      allocate(wframe(2*j1+1,kk,i1))

      allocate(tempu(2*j1+1,Nzb,i1))
      allocate(tempv(2*j1+1,Nzb,i1))
      allocate(tempw(2*j1+1,Nzb,i1))
      end if

      uframe=0.d0
      vframe=0.d0
      wframe=0.d0

      if(me.le.mek.and.me.ne.0)then
      do ii=1,num_turbine
      do j=1,2*j1+1
      do k=2,Nzb+1
      i = 2.d0*wtR(ii)/dx
      tempu(j,k-1,ii)=u(wtx(ii)-i,wty(ii)-j1-1+j,k)
      i = 1
      tempv(j,k-1,ii)=v(wtx(ii)-i,wty(ii)-j1-1+j,k)
      i = 1
      tempw(j,k-1,ii)=w(wtx(ii)-i,wty(ii)-j1-1+j,k)
      end do
      end do
      end do

      call MPI_SEND(tempu(1,1,1),(2*j1+1)*Nzb*num_turbine,
     +     MPI_DOUBLE_PRECISION,0,me,nall, ierr )
      call MPI_SEND(tempv(1,1,1),(2*j1+1)*Nzb*num_turbine,
     +     MPI_DOUBLE_PRECISION,0,me,nall, ierr )
      call MPI_SEND(tempw(1,1,1),(2*j1+1)*Nzb*num_turbine,
     +     MPI_DOUBLE_PRECISION,0,me,nall, ierr )
      end if
    
      if(me.eq.0)then
      do i=1,mek
      call MPI_RECV(tempu(1,1,1),(2*j1+1)*Nzb*num_turbine,
     +     MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
      call MPI_RECV(tempv(1,1,1),(2*j1+1)*Nzb*num_turbine,
     +     MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
      call MPI_RECV(tempw(1,1,1),(2*j1+1)*Nzb*num_turbine,
     +     MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )

        do ii=1,num_turbine
        do j=1,2*j1+1
        do k=1,Nzb
          uframe(j,nzb*i+k,ii)=tempu(j,k,ii)
          vframe(j,nzb*i+k,ii)=tempv(j,k,ii)
          wframe(j,nzb*i+k,ii)=tempw(j,k,ii)
        end do
        end do
        end do
      end do

      k=Nzb*(mek+1)
      do i=1,mek
        call MPI_SEND(uframe(1,1,1),(2*j1+1)*k*num_turbine,
     +       MPI_DOUBLE_PRECISION,i,i,nall, ierr )
        call MPI_SEND(vframe(1,1,1),(2*j1+1)*k*num_turbine,
     +       MPI_DOUBLE_PRECISION,i,i,nall, ierr )
        call MPI_SEND(wframe(1,1,1),(2*j1+1)*k*num_turbine,
     +       MPI_DOUBLE_PRECISION,i,i,nall, ierr )
      end do
      end if

      if(me.le.mek.and.me.ne.0)then
      k=Nzb*(mek+1)
          call MPI_RECV(uframe(1,1,1),(2*j1+1)*k*num_turbine,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
          call MPI_RECV(vframe(1,1,1),(2*j1+1)*k*num_turbine,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
          call MPI_RECV(wframe(1,1,1),(2*j1+1)*k*num_turbine,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
      end if

      do ii=1,num_turbine
      omega(ii)=-wtomega(ii)*(2.d0*pi/60.d0)*(z_i/u_star)
      do N=1,3
      Bang(ii,N)=omega(ii)*dt*dble(t-1)+dble(N-1)*2./3.*pi
      Bang(ii,N)=dmod(Bang(ii,N),pi+pi)
      end do
      end do

      do ii=1,num_turbine
      do N=1,3
      do k=2,nzb+1
      do j=1,ny
      if(CTru(ii,j,k).le.wtr(ii))then
      r  = CTru(ii,j,k)
      dr = dru(ii,j,k)

      blade_jj=dble(j1+1)-r*dcos(Bang(ii,N))/dy
      blade_kk=(zhub(ii)+r*dsin(Bang(ii,N)))/dz+0.5

      CTU=itp(uframe(:,:,ii),blade_jj,blade_kk,2*j1+1,kk,me,nall)
      CTV=itp(vframe(:,:,ii),blade_jj,blade_kk,2*j1+1,kk,me,nall)
      blade_kk=(zhub(ii)+r*dsin(Bang(ii,N)))/dz+1.0
      CTW=itp(wframe(:,:,ii),blade_jj,blade_kk,2*j1+1,kk,me,nall)

c      write(*,*)r,dr,CTU,CTV,CTW

      chordl  = chord(r*z_i)/z_i

      omega1  = +CTW*dcos(Bang(ii,N))/r
     +          +CTV*dsin(Bang(ii,N))/r
     +          +omega(ii)
      CTsoli  = 3.0*chordl/(2.d0*pi*r)

      call sol_CL_CD(omega1,chordl,r,CTU,CTsoli,
     +     a1,b1,CTCL,CTCD,at1,CTphi,Re,me,nall)

      CTV1    =  CTU
c      CTV2    = -CTW*dcos(CTangx(ii,j,k))
c     +          -CTV*dsin(CTangx(ii,j,k))-omega1*r
      CTUrel  = CTU*(1-a1)/dsin(CTphi)

      dA      = dr*chordl

      CTFx    = 0.5*(CTUrel**2)*dA/(dx*dy*dz)*
     +          (+CTCL*dcos(CTphi)+CTCD*dsin(CTphi))
      CTFt    = 0.5*(CTUrel**2)*dA/(dx*dy*dz)*
     +          (+CTCL*dsin(CTphi)-CTCD*dcos(CTphi))

      tmp1=(dabs(Bang(ii,N)-CTangx(ii,j,k)))
      tmp1=dmod(tmp1,pi+pi)
      if(tmp1.gt.pi)tmp1=pi+pi-tmp1
      
      tmpd1=(tmp1-0.5*dt*dabs(omega(ii)))*r
      tmpd2=(tmp1+0.5*dt*dabs(omega(ii)))*r

      tmpd1=tmpd1/dsqrt(2.d0*chordl*chordl*1.)
      tmpd2=tmpd2/dsqrt(2.d0*chordl*chordl*1.)
      call ERROR(tmpd1,ERRd1)
      call ERROR(tmpd2,ERRd2)
      tmp1=ERRd2-ERRd1

      write(*,*)CTFx,tmp1,tmpd1,tmpd2

c        sigma = dx*4.0
c        tmp1  = 1/dsqrt(2.d0*pi*sigma*sigma)
c        tmp2  = -1/(2.d0*sigma*sigma)
c        tmpd1 = tmp1*dexp(tmp2*(tmpd1)**2.)
c        tmpd2 = tmp1*dexp(tmp2*(tmpd2)**2.)
c        tmp1=0.5*(tmpd1+tmpd2)

      Fx (i,j,k)=Fx (i,j,k)+CTFx*tmp1

c      write(*,*)CTFx,tmp1,tmpd1,tmpd2



      Fax(i,j,k)=Fax(i,j,k)+CTFt*tmp1
      Fy (i,j,k)=Fy (i,j,k)+CTFt*tmp1*dsin(CTangx(ii,j,k))
      end if

      if(CTrw(ii,j,k).le.wtr(ii))then
      r  = CTrw(ii,j,k)
      dr = drw(ii,j,k)

      blade_jj=dble(j1+1)-r*dcos(Bang(ii,N))/dy
      blade_kk=(zhub(ii)+r*dsin(Bang(ii,N)))/dz+0.5

      CTU=itp(uframe(:,:,ii),blade_jj,blade_kk,2*j1+1,kk,me,nall)
      CTV=itp(vframe(:,:,ii),blade_jj,blade_kk,2*j1+1,kk,me,nall)
      blade_kk=(zhub(ii)+r*dsin(Bang(ii,N)))/dz+1.0
      CTW=itp(wframe(:,:,ii),blade_jj,blade_kk,2*j1+1,kk,me,nall)

      chordl  = chord(r*z_i)/z_i

      omega1  = +CTW*dcos(Bang(ii,N))/r
     +          +CTV*dsin(Bang(ii,N))/r
     +          +omega(ii)
      CTsoli  = 3.0*chordl/(2.d0*pi*r)

c      call sol_CL_CD_1(omega1,chordl,r,CTU,CTsoli,
c     +     a1,b1,CTCL,CTCD,at1,CTphi,Re,me,nall)
c      CTUrel  = CTU/dsin(CTphi)	  
	  
      call sol_CL_CD(omega1,chordl,r,CTU,CTsoli,
     +     a1,b1,CTCL,CTCD,at1,CTphi,Re,me,nall)

      CTV1    =  CTU
c      CTV2    = -CTW*dcos(CTangx(ii,j,k))
c     +          -CTV*dsin(CTangx(ii,j,k))-omega1*r
      CTUrel  = CTU*(1-a1)/dsin(CTphi)

      dA      = dr*chordl

      CTFt    = 0.5*(CTUrel**2)*dA/(dx*dy*dz)*
     +          (+CTCL*dsin(CTphi)-CTCD*dcos(CTphi))

      tmp1=(dabs(Bang(ii,N)-CTangz(ii,j,k)))
      tmp1=dmod(tmp1,pi+pi)
      if(tmp1.gt.pi)tmp1=pi+pi-tmp1
      
      tmpd1=tmp1-0.5*dt*dabs(omega(ii))
      tmpd2=tmp1+0.5*dt*dabs(omega(ii))

      tmpd1=r*tmpd1/dsqrt(2.d0*chordl*chordl)
      tmpd2=r*tmpd2/dsqrt(2.d0*chordl*chordl)
      call ERROR(tmpd1,ERRd1)
      call ERROR(tmpd2,ERRd2)
      tmp1=ERRd2-ERRd1

      Faz(i,j,k)=Faz(i,j,k)+CTFt*tmp1
      Fz (i,j,k)=Fz (i,j,k)+CTFt*tmp1*dcos(CTangz(ii,j,k))
      end if

      end do
      end do
      end do
      end do

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function chord(r_)
!------------------
!------------------
      real*8 r_, r
      real*8 chord

      r=r_ 
      chord = - 2.8110671517D-8*r**6 + 3.8126820348D-6*r**5 
     +        - 2.2363395740D-4*r**4 + 7.2893308664D-3*r**3 
     +        - 1.3510840879D-1*r**2 + 1.1745120720D+0*r**1  
     +        - 9.4466401330D-2
      
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine sol_CL_CD_1(ii,Pitch_a,omega,chordl,r,bu,soli
     +           ,a1,b1,CL11,CD,at1,phi1,Re1,me,nall)
      implicit none
      include 'dimen.h'
      integer*4 it,i,ii
      real*8 lamdar,omega,r,bu,bu1,soli
      real*8 a1,a2,a3,b1,b2,b3,a1n
      real*8 phi1,phi2,phi3
      real*8 CL11,CL12,CL21,CL22,CL31,CL32,CD
      real*8 at1,at2,at3,res1,res2,Fcr1,Fcr2,Fcr3
      real*8 Urel1,Urel2,Urel3,Re1,Re2,Re3,chordl
      real*8, external :: pitch,CALCL,CALCD,CALCL1,CALCD1
      real*8 Pitch_a	  

      a1=0.0
      b1=0.0
	  
      a1n=a1
      lamdar=omega*r/bu
      phi1 = datan((bu)/(omega*r*(1+b1)))
                    
      Fcr1 = 1.0
 
        at1=phi1-(pitch(r*z_i,wtomega(ii))+Pitch_a)*pi/180.d0		  
cc      at1  = phi1-pitch(r*z_i,wtomega(ii))*pi/180.d0		
cc      if (ii.eq.2) at1=phi1-(pitch(r*z_i,wtomega(ii))+5.0)*pi/180.d0		  
cc      Re1  = Urel1*chordl/nu*1.d0
      CL11 = CALCL(at1*180.d0/pi,r*z_i)
      CD   = CALCD(at1*180.d0/pi,r*z_i)
      a1   = (soli*(CL11*dcos(phi1)+CD*dsin(phi1))/ 
     +       (4*(dsin(phi1))**2))/(1+
     +       (soli*(CL11*dcos(phi1)+CD*dsin(phi1))/ 
     +       (4*(dsin(phi1))**2)))
      b1   = soli*(CL11*dsin(phi1)-CD*dcos(phi1))/ 
     +       (4*lamdar*(dsin(phi1))**2)
      return

      end                           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine sol_CL_CD(omega,chordl,r,bu,soli,a1,b1,CL11,CD
     +           ,at1,phi1,Re1,me,nall)
      implicit none
      include 'dimen.h'
      integer*4 it
      real*8 lamdar,omega,r,bu,soli
      real*8 a1,a2,a3,b1,b2,b3
      real*8 phi1,phi2,phi3
      real*8 CL11,CL12,CL21,CL22,CL31,CL32,CD
      real*8 at1,at2,at3,res1,res2,Fcr1,Fcr2,Fcr3
      real*8 Urel1,Urel2,Urel3,Re1,Re2,Re3,chordl
      real*8, external :: pitch,CALCL,CALCD,CALCL1,CALCD1

      lamdar=omega*r/bu
                    
      a1=0.05
      a2=0.6
      a3=0.5*(a1+a2)

      b1=0
      b2=0
                    
      do it=1,50
         phi1 = datan((bu*(1-a1))/(-omega*r*(1+b1)+0))
         phi2 = datan((bu*(1-a2))/(-omega*r*(1+b1)+0))
         phi3 = datan((bu*(1-a3))/(-omega*r*(1+b1)+0))

c         Fcr1 = 2/pi*dacos(dexp(-3.0*(wtR-r)/(2*r*dsin(phi1))))
c         Fcr2 = 2/pi*dacos(dexp(-3.0*(wtR-r)/(2*r*dsin(phi2))))
c         Fcr3 = 2/pi*dacos(dexp(-3.0*(wtR-r)/(2*r*dsin(phi3))))

         Fcr1 = 1.0
         Fcr2 = 1.0
         Fcr3 = 1.0

         CL11 = 4.d0*Fcr1*dsin(phi1)*(dcos(phi1)-lamdar*dsin(phi1))
     +        / (soli*(dsin(phi1)+lamdar*dcos(phi1)))
         at1  = phi1-pitch(r*z_i,wtomega(1))*pi/180.d0
c         Urel1= bu*(1.d0-a1)/dsin(phi1)
         Urel1= bu/((soli*CL11/(4.*Fcr1))/dtan(phi1)+dsin(phi1))
         Re1  = Urel1*chordl/nu*1.d0
         CL12 = CALCL(at1*180.d0/pi,r*z_i)
c         CL12 = CALCL1(at1*180.d0/pi,Re1)
                        
         CL21 = 4.d0*Fcr2*dsin(phi2)*(dcos(phi2)-lamdar*dsin(phi2))
     +        / (soli*(sin(phi2)+lamdar*cos(phi2)))
         at2  = phi2-pitch(r*z_i,wtomega(1))*pi/180.d0
c         Urel2= bu*(1.d0-a2)/dsin(phi2)
         Urel2= bu/((soli*CL21/(4.*Fcr2))/dtan(phi2)+dsin(phi2))
         Re2  = Urel2*chordl/nu*1.d0
         CL22 = CALCL(at2*180.d0/pi,r*z_i)
c         CL22 = CALCL1(at2*180.d0/pi,Re2)
                        
         CL31 = 4.d0*Fcr3*dsin(phi3)*(dcos(phi3)-lamdar*dsin(phi3))
     +        / (soli*(dsin(phi3)+lamdar*dcos(phi3)))
         at3  = phi3-pitch(r*z_i,wtomega(1))*pi/180.d0
c         Urel3= bu*(1.d0-a3)/dsin(phi3)
         Urel3= bu/((soli*CL31/(4.*Fcr3))/dtan(phi3)+dsin(phi3))
         Re3  = Urel3*chordl/nu*1.d0
         CL32 = CALCL(at3*180.d0/pi,r*z_i)
c         CL32 = CALCL1(at3*180.d0/pi,Re3)
                        

      if(me.eq.0)then
c      write(*,500)it,phi1,phi2,phi3,soli,lamdar,bu,r,omega,Re1,Re2
c      write(*,500)it,a1,a2,a3,CL11,CL12,CL21,CL22,CL31,CL32,at1*180/pi
c      if(it.eq.50)write(*,500)it,a1,a2,a3,CL11,CL12,CL21,CL22,at1*180/pi
      end if

500   format(I2,10(3x,f10.4))

         if (CL11.gt.CL12.and.CL31.gt.CL32)then
            a1=a3
            a3=(a1+a2)/2.d0
         end if
         if (CL11.lt.CL12.and.CL31.lt.CL32)then
            a1=a3
            a3=(a1+a2)/2.d0
         end if                     
         if (CL21.gt.CL22.and.CL31.gt.CL32)then
            a2=a3 
            a3=(a1+a2)/2.d0
         end if
         if (CL21.lt.CL22.and.CL31.lt.CL32)then
            a2=a3
            a3=(a1+a2)/2.d0
         end if
                        
         res1 = abs(CL11-CL21)
         res2 = abs(CL31-CL21)

         if (res1.le.0.001.and.res2.le.0.001) then
            goto 1001
         end if
      enddo

1001     a1 = 1.d0/(1.d0+4.d0*Fcr1*dsin(phi1)**2/(soli*CL11*dcos(phi1)))
         b1 = 1.d0/(4.d0*Fcr1*dcos(phi1)/(soli*CL11)-1.d0)
         CD  = CALCD(at1*180.d0/pi,r*z_i)
c         if (Re1.gt.4000) CD = CD*(4000/Re1)
c         CD  = CALCD1(at1*180.d0/pi,Re1)
         phi1 = datan((bu*(1-a1))/(-omega*r*(1+b1)+0))
         CL11 = 4.d0*Fcr1*dsin(phi1)*(dcos(phi1)-lamdar*dsin(phi1))
     +        / (soli*(dsin(phi1)+lamdar*dcos(phi1)))
         if(me.eq.0)then
c         write(*,500)it,r*z_i,a1,b1,CL11,CD,at1*180/pi,phi1*180/pi
c         write(6855,500)it,r*z_i,a1,b1,CL11,at1*180/pi,phi1*180/pi
c         write(*,500)it,phi1,phi2,phi3,soli,lamdar,bu,r,omega,Re1,Re2
         end if

      return

      end                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function pitch(r_,RPM)
!------------------
!------------------
      real*8 r_,r
      real*8 Measured_Pitch, RPM, pitch

      r=r_
      pitch=-2.0921301914D-4*r**3+3.2105894326D-2*r**2  
     +      -1.5300370009D+0*r**1+2.3276875553D+1
      
      Measured_Pitch=-5.8320471905D-4*RPM**6+4.6419810970D-2*RPM**5  
     +  -1.5008145878D+0*RPM**4 + 2.4977428807D+1*RPM**3
     +  -2.2171588200D+2*RPM**2 + 9.5971735595D+2*RPM**1 
     +  - 1.4373108828D+3

cc      print*, r_,pitch,Measured_Pitch	 
      pitch = pitch + Measured_Pitch + 0.0D0
!---------------------------------------------------------------------
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function CALCD1(x,Re)
      include 'dimen.h'
      real*8 Re,x,CALCD1

      if(Re.le.3500)then      
CCCCCC NACA 2309 from Xfoil RE=2,000
      if(x.le.20.)then
      CALCD1 = 4.104787E-09*x**6 - 2.044210E-07*x**5 + 2.706574E-06*x**4
     +       - 2.812527E-06*x**3 + 4.233886E-04*x**2 + 6.675463E-04*x**1
     +       + 7.860550E-02
      else
      CALCD1 =-2.360187E-10*x**6 + 7.343527E-08*x**5 - 9.124715E-06*x**4
     +       + 5.744508E-04*x**3 - 1.926885E-02*x**2 + 3.431459E-01*x**1
     +       - 2.238189E+00
      end if
      elseif(Re.gt.3500.and.Re.le.6500)then
CCCCCC NACA 2309 from Xfoil RE=5,000
      if(x.le.20.)then
      CALCD1 = 1.484731E-08*x**6 - 7.426866E-07*x**5 + 1.133017E-05*x**4
     +       - 4.525399E-05*x**3 + 4.338993E-04*x**2 + 9.560049E-04*x**1
     +       + 4.975932E-02
      else
      CALCD1 =-2.300255E-10*x**6 + 7.156753E-08*x**5 - 8.894411E-06*x**4
     +       + 5.602613E-04*x**3 - 1.880877E-02*x**2 + 3.351148E-01*x**1
     +       - 2.192792E+00
      end if
      elseif(Re.gt.6500.and.Re.le.9500)then
CCCCCC NACA 2309 from Xfoil RE=8,000
      if(x.le.20.)then
      CALCD1 = 2.442313E-08*x**6 - 1.214148E-06*x**5 + 1.813410E-05*x**4
     +       - 5.854174E-05*x**3 + 2.294585E-04*x**2 + 1.044494E-03*x**1
     +       + 4.010342E-02
      else
      CALCD1 =-2.268082E-10*x**6 + 7.055546E-08*x**5 - 8.767933E-06*x**4
     +       + 5.522824E-04*x**3 - 1.854033E-02*x**2 + 3.303485E-01*x**1
     +       - 2.161266E+00
      end if
      elseif(Re.gt.9500.and.Re.le.12500)then
CCCCCC NACA 2309 from Xfoil RE=11,000
      if(x.le.20.)then
      CALCD1 = 3.956439E-08*x**6 - 1.995309E-06*x**5 + 3.055711E-05*x**4
     +       - 1.006242E-04*x**3 - 9.429684E-05*x**2 + 1.721290E-03*x**1
     +       + 3.552569E-02
      else
      CALCD1 =-2.314883E-10*x**6 + 7.208434E-08*x**5 - 8.967408E-06*x**4
     +       + 5.655857E-04*x**3 - 1.901677E-02*x**2 + 3.388775E-01*x**1
     +       - 2.220659E+00
      end if
      else
CCCCCC NACA 2309 from Xfoil RE=14,000
      if(x.le.20.)then
      CALCD1 = 4.255185E-08*x**6 - 2.199702E-06*x**5 + 3.484976E-05*x**4
     +       - 1.245836E-04*x**3 - 1.899481E-04*x**2 + 2.034694E-03*x**1
     +       + 3.181027E-02
      else
      CALCD1 =-2.236269E-10*x**6 + 6.956973E-08*x**5 - 8.646708E-06*x**4
     +       + 5.447765E-04*x**3 - 1.829609E-02*x**2 + 3.263042E-01*x**1
     +       - 2.137689E+00
      end if
      end if

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        SUBROUTINE ERROR(X,ERR)
C
C       =========================================
C       Purpose: Compute error function erf(x)
C       Input:   x   --- Argument of erf(x)
C       Output:  ERR --- erf(x)
C       =========================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EPS=1.0D-15
        PI=3.141592653589793D0
        X2=X*X
        IF (DABS(X).LT.3.5D0) THEN
           ER=1.0D0
           R=1.0D0
           DO 10 K=1,50
              R=R*X2/(K+0.5D0)
              ER=ER+R
              IF (DABS(R).LE.DABS(ER)*EPS) GO TO 15
10         CONTINUE
15         C0=2.0D0/DSQRT(PI)*X*DEXP(-X2)
           ERR=C0*ER
        ELSE
           ER=1.0D0
           R=1.0D0
           DO 20 K=1,12
              R=-R*(K-0.5D0)/X2
20            ER=ER+R
           C0=DEXP(-X2)/(DABS(X)*DSQRT(PI))
           ERR=1.0D0-C0*ER
           IF (X.LT.0.0) ERR=-ERR
        ENDIF
        RETURN
        END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine pass_slice(ac,at,me,nall)
      implicit none
      include 'dimen.h'

      integer*4 fid,i,k
      real*8, dimension(ny,nz2):: ac
      real*8, dimension(ny,nz)::  at
         
      IF (me>0) then
         call MPI_SEND(ac(1,2),nzb*ny,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
      ELSE
         do k=1,nzb
            at(:,k)=ac(:,k+1)
         enddo  
         do i=1,nprocs-1
               call MPI_RECV(at(1,i*nzb+1),nzb*ny,
     +              MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         end do
         do i=1,nprocs-1
         call MPI_SEND(at(1,1),ny*nz,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )
         end do
      END IF
      IF (me>0) then
          call MPI_RECV(at(1,1),ny*nz,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
      END IF

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine calsum(sum_kernal,max_kernal,me,nall)
      implicit none
      include 'dimen.h'

      real*8 sum_kernal,max_kernal,tmp
      integer*4 fid,i,k
         
      IF (me>0) then
         call MPI_SEND(sum_kernal,1,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(max_kernal,1,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
      ELSE
         tmp = sum_kernal
         do i=1,nprocs-1
               call MPI_RECV(sum_kernal,1,
     +              MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         tmp = tmp + sum_kernal
         end do
         sum_kernal=tmp

         tmp = max_kernal
         do i=1,nprocs-1
               call MPI_RECV(max_kernal,1,
     +              MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         tmp = dmax1(tmp, max_kernal)
         end do
         max_kernal=tmp


         do i=1,nprocs-1
         call MPI_SEND(sum_kernal,1,
     +        MPI_DOUBLE_PRECISION,i,me,nall, ierr )

         call MPI_SEND(max_kernal,1,
     +        MPI_DOUBLE_PRECISION,i,me,nall, ierr )
         end do
      END IF
      IF (me>0) then
          call MPI_RECV(sum_kernal,1,
     +         MPI_DOUBLE_PRECISION,0,0,nall,status2,ierr )

          call MPI_RECV(max_kernal,1,
     +         MPI_DOUBLE_PRECISION,0,0,nall,status2,ierr )
      END IF

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function calang(CTy,CTz)
      include 'dimen.h'
      real*8 CTy,CTz,gau,calang,RR
         RR = dsqrt(CTy**2+CTz**2)
         if (RR.le.1E-8) then
            calang = 0.d0
         else
            gau = abs(CTz/dsqrt(CTy**2+CTz**2))
            if (CTz.ge.0.0)then
               if(CTy.ge.0.0) calang = dasin(gau)
               if(CTy.lt.0.0) calang = PI-abs(dasin(gau))
            else
               if(CTy.lt.0.0) calang = PI+abs(dasin(gau))
               if(CTy.ge.0.0) calang = 2.0*PI-dasin(gau)
            end if   
         end if
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function itp2D(ac,ii,blade_jj,blade_kk,me,nall)
      include 'dimen.h'
      integer ii,blade_j,blade_k
      real*8, dimension(nx,ny,nz2)::ac
      real*8 blade_jj,blade_kk,ly1,ly2,lz1,lz2
      real*8 itp2D


      blade_j = int(blade_jj) 
      blade_k = int(blade_kk) 

      ly1 = dble(blade_jj-blade_j)
      ly2 = (1.0-ly1)
      lz1 = dble(blade_kk-blade_k)
      lz2 = (1.0-lz1)

      itp2D=ac(ii,blade_j  ,blade_k  )*ly2*lz2
     +     +ac(ii,blade_j+1,blade_k  )*ly1*lz2
     +     +ac(ii,blade_j  ,blade_k+1)*ly2*lz1
     +     +ac(ii,blade_j+1,blade_k+1)*ly1*lz1

c      if(me.eq.0)write(*,*)ac(ii,blade_j,blade_k),interp2D

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function itp(ac,blade_jj,blade_kk,jj,kk,me,nall)
      include 'dimen.h'
      integer ii,blade_j,blade_k
      real*8, dimension(jj,kk)::ac
      real*8 blade_jj,blade_kk,ly1,ly2,lz1,lz2
      real*8 itp

      blade_j = int(blade_jj) 
      blade_k = int(blade_kk) 

      ly1 = dble(blade_jj-blade_j)
      ly2 = (1.0-ly1)
      lz1 = dble(blade_kk-blade_k)
      lz2 = (1.0-lz1)

      itp  =ac(blade_j  ,blade_k  )*ly2*lz2
     +     +ac(blade_j+1,blade_k  )*ly1*lz2
     +     +ac(blade_j  ,blade_k+1)*ly2*lz1
     +     +ac(blade_j+1,blade_k+1)*ly1*lz1

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function interp2D_slice(ac,blade_jj,blade_kk)
      include 'dimen.h'
      integer blade_j,blade_k
      real*8, dimension(Ny,Nz)::ac
      real*8 blade_jj,blade_kk,interp2D_slice,ly1,ly2,lz1,lz2

      blade_j = int(blade_jj) 
      blade_k = int(blade_kk) 

      ly1 = dble(blade_jj-blade_j)
      ly2 = (1.0-ly1)
      lz1 = dble(blade_kk-blade_k)
      lz2 = (1.0-lz1)

      interp2D_slice=ac(blade_j  ,blade_k  )*ly2*lz2
     +              +ac(blade_j+1,blade_k  )*ly1*lz2
     +              +ac(blade_j  ,blade_k+1)*ly2*lz1
     +              +ac(blade_j+1,blade_k+1)*ly1*lz1
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function calU(x)
      real*8 x,calU
     
      calU = -13.077*x**2 + 6.4431*x + 1.6279
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

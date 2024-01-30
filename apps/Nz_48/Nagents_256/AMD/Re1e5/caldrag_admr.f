      Subroutine calturbine(alpha,u_hub,v_hub,u,v,w,
     +           Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     turbine = 1 --> Actuator Disk model 
C     turbine = 2 --> An Improved Actuator Disk model 
C     turbine = 3 --> Actuator Line model 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'dimen.h'
      integer*4 t,k

      real*8, dimension(Nx,Ny,Nz2)::u,v,w,fx,fy,fz,fax,faz,AveWt
      real*8 alpha,u_hub,v_hub

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Fx = 0.0d0
      Fy = 0.0d0
      Fz = 0.0d0
  
      if(nacelle_model.eq.1)then
      call CALNac   (u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
      end if

      if(tower_model.eq.1)then
      call CALBar   (u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
      end if

      if(turbine_model.eq.1)then
      call CALADMNR (alpha,u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
      end if

      if(turbine_model.eq.2)then
      call CALADMR  (u_hub,v_hub,u,v,w,Fx,Fy,Fz,
     +               Fax,Faz,AveWt,t,me,nall)
      end if

c      if(turbine_model.eq.3)then
c      call CALALM   (u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
c      end if

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine CALBar(u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'dimen.h'
      integer*4 ii,i,j,k,t

      real*8, dimension(Nx,Ny,Nz2)::u,v,w,fx,fy,fz,fax,faz,AveWt
      real*8, dimension(Ny,Nz2):: CTy,CTzu,CTru,CTangx
      real*8, dimension(num_turbine,Ny):: ker_y
      real*8, dimension(Ny,Nz2)::Nac_u

      real*8 a1,b1,at1,Tang(1:36),blade_jj,blade_kk,tmp,tmp1,tmp2
      real*8 mu_y,mu_z,sigma,hz1,hz2
      real*8 tmpy1,tmpz1,tmpy2,tmpz2
      real*8 ERRd1,ERRd2,tmpd1,tmpd2
      real*8 sum_ker_y,max_Nac_u
      real*8 CTU,dA,CTFx,dr

      real*8, external :: interp2D,calang

      save ker_y,Nac_u
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if (t.eq.1) then
        sigma = 0.5d0
        tmp1  = 1/dsqrt(2.d0*pi*sigma*sigma)
        tmp2  = -1/(2.d0*sigma*sigma)
        do ii=1,num_turbine
        do i=1,Ny
          ker_y(ii,i) = tmp1*dexp(tmp2*(i-wty(ii))**2.)
        end do
        tmp1 = sum(ker_y(ii,1:Ny))
        do i=1,Ny
          ker_y(ii,i) = ker_y(ii,i)/tmp1
        end do
        end do
      end if

      do ii=1,num_turbine
      do k=2,nzb+1
        hz1= dble(me*nzb+k-1.5-0.5)*dz
        hz2= dble(me*nzb+k-1.5+0.5)*dz

        if(Zhub(ii)*0.94d0.gt.hz1)then
          i = 2.d0*wtR(ii)/dx
          CTU = u(wtx(ii)-i,nint(wty(ii)),k)
          if(CTU.lt.0.d0)CTU=0.d0
          dA = (3.0d0/z_i)*dz
          if(Zhub(ii)*0.94d0.lt.hz2)then
      	    hz2=Zhub(ii)*0.94d0
            dA =(3.0d0/z_i)*(hz2-hz1)
          endif
        CTFx   = 0.5*(CTU**2)*CTB(ii)*dA/(dx*dy*dz)
        do j=1,ny
        if (ker_y(ii,j).ge.1E-15)then
          Fx (wtx(ii),j,k)=Fx (wtx(ii),j,k) + CTFx * ker_y(ii,j)
        end if
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
      i = 4.d0*wtR(ii)/dx
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
      Subroutine CALADMNR(angle,u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'dimen.h'
      integer*4 ii,i,j,k,t,blade,N,me_k

      real*8, dimension(Nx,Ny,Nz2)::u,v,w,fx,fy,fz,fax,faz,AveWt
      real*8, dimension(num_turbine,Ny,Nz2):: CTy,CTzu,CTzw,CTru,CTrw
      real*8, dimension(num_turbine,Ny,Nz2):: CTangx,CTangz
      real*8, dimension(num_turbine,Nx):: ker_x
      real*8 tmp1,tmp2,tmp3

      real*8, dimension(num_turbine)::
     +uframe,vframe,wframe,meframe,u1frame,v1frame,w1frame,me1frame

      real*8, dimension(num_turbine,N_ang,N_blade,Ny,Nz2):: ker_u,ker_w

      real*8 a1,b1,at1,Tang(N_ang),blade_jj,blade_kk
      real*8 mu_y,mu_z,sigma,dr,omega,omega1,r,CTphi
      real*8 tmpy1,tmpz1,tmpy2,tmpz2,tmp
      real*8 ERRd1,ERRd2,tmpd1,tmpd2
      real*8 sum_ker_u,sum_ker_w,max_ker_u,max_ker_w
      real*8 CTU1,CTV,CTW,chordl,CTsoli
      real*8 CTV1,CTV2,CTUrel,CTFx,CTFt
      real*8 CTCL,CTCD,dang,Re,dA,RA,dL
      real*8,dimension(num_turbine) :: Ua,Va,Uanew,Vanew
      
      real*8 Pow(num_turbine),CTU(num_turbine),CTUnew(num_turbine),
     +       anglein(num_turbine)

      real*8, external :: chord,calang
      real*8, external :: itp2D,calU

c     Mahdi 

      real*8,dimension(nx,ny,nz2):: u_hat,u_hatd,v_hat,v_hatd,
     +                              w_hat,w_hatd    
      real*8, dimension(num_turbine,Ny,Nz2):: gama
      real*8 rr(4),y1,z1,yy1,xx1,rr1(4),wp(4),f,angle,Powmean,CTUmean,
     +       CTUnewmean 
      integer*4 nn,jm,km,k1,j1,n1y,n2y,n1x,n2x,i1,ii1
      integer*4, dimension(num_turbine) :: wnpx,wnpy 
 	  
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
		  
          rr(1)=((CTy(ii,j,k)+dy/2)**2+(CTzu(ii,j,k)+dz/2)**2)**0.5
          rr(2)=((CTy(ii,j,k)+dy/2)**2+(CTzu(ii,j,k)-dz/2)**2)**0.5
          rr(3)=((CTy(ii,j,k)-dy/2)**2+(CTzu(ii,j,k)+dz/2)**2)**0.5
          rr(4)=((CTy(ii,j,k)-dy/2)**2+(CTzu(ii,j,k)-dz/2)**2)**0.5
		  
          if (maxval(rr)<wtr(ii)) then
           gama(ii,j,k)=1.
          else if (minval(rr)>wtr(ii)) then
           gama(ii,j,k)=0.
          else
          jm=50
          km=50
          nn=0
          do k1=1,km
          do j1=1,jm
                y1=(CTy(ii,j,k)-dy/2)+dy/(jm-1)*(j1-1)
                z1=(CTzu(ii,j,k)-dz/2)+dz/(km-1)*(k1-1)
                if (((y1**2+z1**2)**0.5) <wtr(ii)) nn=nn+1
          end do
          end do
          gama(ii,j,k)=1.0*nn/(1.*jm*km) 
          endif
        end do
        end do
        end do

      end if

      if (t.eq.1) then
      do i=1,N_ang	
      Tang(i) = (dble(i)-1)*dang
      end do
      wnpx=wtr/dx+1
      wnpy=wtr/dy+1
      end if

c      if (me==5) write(*,*) wnpx,wnpy 

      meframe=0.d0
      uframe=0.d0
      vframe=0.d0
      wframe=0.d0

      do ii=1,num_turbine

      blade_kk = dble(Zhub(ii)/dz-0.5)
      me_k=blade_kk/Nzb

      if (me.eq.me_k) then
      blade_kk = blade_kk-Nzb*me+2.d0
      blade_jj = dble(wty(ii))	  
      meframe(ii)=1.d0

      uframe(ii)=itp2D(u,wtx(ii),blade_jj,blade_kk,me,nall)
      vframe(ii)=itp2D(v,wtx(ii),blade_jj,blade_kk,me,nall)
      wframe(ii)=itp2D(w,wtx(ii),blade_jj,blade_kk,me,nall)
c       uframe(ii)=u(wtx(ii),wty(ii),blade_kk)
c       vframe(ii)=v(wtx(ii),wty(ii),blade_kk)
c       wframe(ii)=w(wtx(ii),wty(ii),blade_kk)
	  
c      if (mod(t,c_count)==0) then
c        print*,u(wtx(ii),wty(ii),blade_kk),
c     +         itp2D(u,wtx(ii),blade_jj,blade_kk,me,nall),
c     +         itp2D(u,wtx(ii),wty(ii),blade_kk,me,nall),
c     +         blade_kk
c      endif

      end if

      end do

      IF (me>0) then
         call MPI_SEND(meframe(1),num_turbine,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(uframe(1),num_turbine,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(vframe(1),num_turbine,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(wframe(1),num_turbine,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
      ELSE
         do i=1,nprocs-1
         call MPI_RECV(me1frame(1),num_turbine,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         call MPI_RECV(u1frame(1),num_turbine,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         call MPI_RECV(v1frame(1),num_turbine,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         call MPI_RECV(w1frame(1),num_turbine,
     +        MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )

         tmp1=sum(me1frame)
         if(tmp1.gt.0.5)then
           do ii=1,num_turbine
           if(me1frame(ii).gt.0.5d0)then
           uframe(ii)=u1frame(ii)
           vframe(ii)=v1frame(ii)
           wframe(ii)=w1frame(ii)
           end if
           end do
         end if
         end do

         do i=1,nprocs-1
         call MPI_SEND(uframe(1),num_turbine,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )
         call MPI_SEND(vframe(1),num_turbine,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )
         call MPI_SEND(wframe(1),num_turbine,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )
         end do
      END IF
      IF (me>0) then
          call MPI_RECV(uframe(1),num_turbine,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
          call MPI_RECV(vframe(1),num_turbine,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
          call MPI_RECV(wframe(1),num_turbine,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
      END IF

c      write(*,*) uframe,vframe,wframe

      do ii=1,num_turbine

      Uanew(ii)=uframe(ii)
      Vanew(ii)=vframe(ii)
      CTUnew(ii) =(Uanew(ii)**2.+Vanew(ii)**2.)**0.5

       if (t.eq.1) Ua(ii)=Uanew(ii)
       Ua(ii) = 1.0*Uanew(ii)+0.0*Ua(ii)
       if (t.eq.1) Va(ii)=Vanew(ii)
       Va(ii) = 1.0*Vanew(ii)+0.0*Va(ii)
       CTU(ii) =(Ua(ii)**2.+Va(ii)**2.)**0.5
       anglein(ii)=angle
c      anglein(ii)=atan(Va(ii)/Ua(ii))		

c      Va(ii)=0.0d0
c      anglein(ii)=0.0d0
c      CTU(ii) =Ua(ii) 
		
c      i=wtx(ii)
	  
      do k=2,nzb+1
      do j=wty(ii)-wnpy(ii),wty(ii)+wnpy(ii)
 
      F=0.5*(CTU(ii)**2)*CTT(ii)*gama(ii,j,k)/dx

      yy1=(j-wty(ii))*dy*cos(anglein(ii))
      xx1=-1.0*(j-wty(ii))*dy*sin(anglein(ii))
      n1y=floor(yy1/dy)
      n2y=n1y+1
      n1x=floor(xx1/dx)
      n2x=n1x+1
      rr1(1)=((yy1-n1y*dy)**2+(xx1-n1x*dx)**2)**0.5
      rr1(2)=((yy1-n1y*dy)**2+(xx1-n2x*dx)**2)**0.5
      rr1(3)=((yy1-n2y*dy)**2+(xx1-n1x*dx)**2)**0.5
      rr1(4)=((yy1-n2y*dy)**2+(xx1-n2x*dx)**2)**0.5
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
      
c      if (me==0) print*,wp	  

      Fx(n1x+wtx(ii),n1y+wty(ii),k)=Fx(n1x+wtx(ii),n1y+wty(ii),k) +
     +              F*wp(1)*Ua(ii)/CTU(ii)
      Fy(n1x+wtx(ii),n1y+wty(ii),k)=Fy(n1x+wtx(ii),n1y+wty(ii),k) + 
     +              F*wp(1)*Va(ii)/CTU(ii)	 
      Fx(n2x+wtx(ii),n1y+wty(ii),k)=Fx(n2x+wtx(ii),n1y+wty(ii),k) + 
     +              F*wp(2)*Ua(ii)/CTU(ii)	 
      Fy(n2x+wtx(ii),n1y+wty(ii),k)=Fy(n2x+wtx(ii),n1y+wty(ii),k) + 
     +              F*wp(2)*Va(ii)/CTU(ii)
      Fx(n1x+wtx(ii),n2y+wty(ii),k)=Fx(n1x+wtx(ii),n2y+wty(ii),k) + 
     +              F*wp(3)*Ua(ii)/CTU(ii)
      Fy(n1x+wtx(ii),n2y+wty(ii),k)=Fy(n1x+wtx(ii),n2y+wty(ii),k) + 
     +              F*wp(3)*Va(ii)/CTU(ii)
      Fx(n2x+wtx(ii),n2y+wty(ii),k)=Fx(n2x+wtx(ii),n2y+wty(ii),k) + 
     +              F*wp(4)*Ua(ii)/CTU(ii)
      Fy(n2x+wtx(ii),n2y+wty(ii),k)=Fy(n2x+wtx(ii),n2y+wty(ii),k) + 
     +              F*wp(4)*Va(ii)/CTU(ii)	 
	 
      end do
      end do
  
      Pow(ii)=0.5*CTT(ii)*CTU(ii)**3*pi*((wtR(ii)*z_i)**2)
      end do
	  
        CTUmean=sum(CTU)/(max(1,size(CTU)))
        CTUnewmean=sum(CTUnew)/(max(1,size(CTUnew)))
        Powmean=sum(Pow)/(max(1,size(Pow)))
		
        if (mod(t,c_count)==0 .and. me==0) then
        write(1524,5544) Pow/(10**6.)
        write(134,5544) CTUmean,Powmean/(10**6.)
		
c        write(*,5543) 'power',  Pow/(10**6.),Powmean/(10**6.)
c        write(*,5543) 'CTU', CTU,CTUmean
        write(*,5543) 'mean P&M',Powmean/(10**6.),CTUmean
        write(*,5543) 'angle', angle*180/pi		
        endif

5543   format (a10,1x,33(f11.5,1x))
5544   format (33(f11.5,1x))
		
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine CALADMR(u_hub,v_hub,u,v,w,Fx,Fy,Fz,
     +                   Fax,Faz,AveWt,t,me,nall)
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
     +       omegamean,u_hub,v_hub,dir(num_turbine)

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
      i = 4.0*wtR(ii)/dx
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

      dir(ii)=1.0d0 
      dr = wtr(ii)/dble(N_blade)

      umean(ii)=sum(uframe(ii,:,:))/(dble(N_ang)*N_blade)
      umeand(ii)=sum(uframed(ii,:,:))/(dble(N_ang)*N_blade)

	  
c      Power_max(ii)=0d0	 	  
c      wtomega(ii)=12.0d0
      do iter_omega=1,10
c      wtomega(ii)=wtomega(ii)+0.5d0	
      Pitch_t(ii)=0.0d0
c      do iter_pitch=1,21
c      Pitch_t(ii)=Pitch_t(ii)+0.5d0
      SUMQ(ii)=0.0d0
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
      call sol_CL_CD_1(ii,Pitch_t(ii),omega1,chordl,r,CTU,CTsoli,
     +     a1,b1,CTCL,CTCD,at1,CTphi,Re,me,nall)
      CTCL=0.9*CTCL
      CTCD=0.9*CTCD	 
      CTUrel  = CTU/dsin(CTphi)
c      call sol_CL_CD(ii,Pitch_t(ii),omega1,chordl,r,CTU,CTsoli,
c     +     a1,b1,CTCL,CTCD,at1,CTphi,Re,me,nall)
c      CTUrel  = CTU*(1-a1)/dsin(CTphi)

      dA      = 0.5d0*(dble(N)**2-dble(N-1)**2)*(dr**2.)*dang
      CTFx    = 0.5*(CTUrel**2)*CTsoli*dA/(dx*dy*dz)*
     +          (+CTCL*dcos(CTphi)+CTCD*dsin(CTphi))
      CTFt    = -dir(ii)*0.5*(CTUrel**2)*CTsoli*dA/(dx*dy*dz)*
     +          (+CTCL*dsin(CTphi)-CTCD*dcos(CTphi))
      SUMQ(ii)=SUMQ(ii)-dir(ii)*CTFt*(dx*dy*dz)*r*(z_i**3)
      end do
      end do  
      Power(ii)=SUMQ(ii)*wtomega(ii)*(2.d0*pi/60)	  
c      if (Power(ii).gt.Power_max(ii)) then
c      Power_max(ii)=Power(ii)
c          wtomega_max(ii)=wtomega(ii)
c          Pitch_max(ii)=Pitch_t(ii) 		  
c      endif	  
c      if (me.eq.0) print*, wtomega(ii),Pitch_t(ii),Power(ii)/10**6	  
c      wtomega(ii)=wtomega(ii)+0.5*(
c     +            +5.648D-11*(1.0*SUMQ(ii)/(10**3))**4
c     +            -1.365D-7*(1.0*SUMQ(ii)/(10**3))**3
c     +            +0.0001038*(1.0*SUMQ(ii)/(10**3))**2
c     +            -0.01895*(1.0*SUMQ(ii)/(10**3))
c     +            +13.23-wtomega(ii))	
      wtomega(ii)=wtomega(ii)+0.5*(
     +            +7.655D-11*(1.0*SUMQ(ii)/(10**3))**4
     +            -1.584D-7*(1.0*SUMQ(ii)/(10**3))**3
     +            +0.0001041*(1.0*SUMQ(ii)/(10**3))**2
     +            -0.01292*(1.0*SUMQ(ii)/(10**3))
     +            +12.37-wtomega(ii))	 
      if (wtomega(ii).gt.18.5) wtomega(ii)=18.5
      if (wtomega(ii).lt.11.5) wtomega(ii)=11.5
c      enddo
      enddo
c      wtomega(ii)=wtomega_max(ii) 	  
c      Pitch_t(ii)=Pitch_max(ii)
c      if (me.eq.0) print*, umean(ii),wtomega(ii),Pitch_t(ii)
c      wtomega(ii)=11.9d0 	  
c      Pitch_t(ii)=0.0
       	  
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
      CTCL=0.9*CTCL
      CTCD=0.9*CTCD	  
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

      SUMdA(ii) =SUMdA(ii)+dA
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

      write(1524,5544) Power/(10**6)
      write(1523,5544) wtomega

      write(134,5111) a1mean,b1mean,CTmean,
     +  Fxmean/(0.5d0*(pi*wtr(1)**2)*(CTUmean**2.)),
     +  Fxmean/(0.5d0*(pi*wtr(1)**2)*(CTUdmean**2.)),	 
     +  4*a1mean*(1-a1mean),
     +  Fxmean*(z_i**2)/(10**6),Powmean/(10**6),
     +  Fxmean*CTUmean*(z_i**2)/(10**6),	 
     +  v_hub,u_hub,
     +  CTUmean,CTUdmean,
     +  Omegamean,	 
     +  atan(v_hub/u_hub)*180/pi 	 
	  
      write(1095,5111)SUMa1,SUMb1,SUMCT,1-(SUMCTU/SUMCTUd),
     +                SUMCT*(1-SUMa1)**2

c      write(*,5111) a1mean,b1mean,CTmean,
c     +  Fxmean/(0.5d0*(pi*wtr(1)**2)*(CTUmean**2.)),
c     +  Fxmean/(0.5d0*(pi*wtr(1)**2)*(CTUdmean**2.)),	 
c     +  4*a1mean*(1-a1mean),
c     +  Fxmean*(z_i**2)/(10**6),Powmean/(10**6),
c     +  Fxmean*CTUmean*(z_i**2)/(10**6),	 
c     +  v_hub,u_hub,
c     +  CTUmean,CTUdmean,
c     +  Omegamean,	 
c     +  atan(v_hub/u_hub)*180/pi 	 	 
	  
c      j=1
c      write(*,5111)SUMa1(j),SUMb1(j),SUMCT(j),
c     +  4.0*SUMa1(j)/(1-SUMa1(j)),4.0*SUMa1(j)*(1-SUMa1(j)),SUMCTU(j), 	  
c     +  SUMFX(j),SUMCTUd(j),1-(SUMCTU(j)/SUMCTUd(j)),
c     +  SUMCT(j)*(1-SUMa1(j))**2,
c     +  SUMFX(j)/(0.5*SUMCTU(j)**2*(pi*wtR(1)**2.d0)),
c     +  umean(j),wtomegaL(j),umeand(j),
c     +  SUMQ(j)/(10**6),SUMFX(j)*SUMCTU(j)*(z_i**2)/(10**6)

      endif 
	 
      end if

 5111 format (18(1x,f11.5)) 
 5544 format (56(f11.5,1x))

 
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
      subroutine sol_CL_CD(ii,Pitch_a,omega,chordl,r,bu,soli,a1,b1,CL11
     +           ,CD,at1,phi1,Re1,me,nall)
      implicit none
      include 'dimen.h'
      integer*4 it,ii
      real*8 lamdar,omega,r,bu,soli
      real*8 a1,a2,a3,b1,b2,b3
      real*8 phi1,phi2,phi3
      real*8 CL11,CL12,CL21,CL22,CL31,CL32,CD
      real*8 at1,at2,at3,res1,res2,Fcr1,Fcr2,Fcr3
      real*8 Urel1,Urel2,Urel3,Re1,Re2,Re3,chordl
      real*8, external :: pitch,CALCL,CALCD,CALCL1,CALCD1
      real*8 Pitch_a	  

      lamdar=omega*r/bu
                    
      a1=0.05
      a2=0.6
      a3=0.5*(a1+a2)

      b1=0
      b2=0
                    
      do it=1,20
         phi1 = datan((bu*(1-a1))/(omega*r*(1+b1)+0))
         phi2 = datan((bu*(1-a2))/(omega*r*(1+b1)+0))
         phi3 = datan((bu*(1-a3))/(omega*r*(1+b1)+0))

c         Fcr1 = 2/pi*dacos(dexp(-3.0*(wtR-r)/(2*r*dsin(phi1))))
c         Fcr2 = 2/pi*dacos(dexp(-3.0*(wtR-r)/(2*r*dsin(phi2))))
c         Fcr3 = 2/pi*dacos(dexp(-3.0*(wtR-r)/(2*r*dsin(phi3))))

         Fcr1 = 1.0
         Fcr2 = 1.0
         Fcr3 = 1.0

         CL11 = 4.d0*Fcr1*dsin(phi1)*(dcos(phi1)-lamdar*dsin(phi1))
     +        / (soli*(dsin(phi1)+lamdar*dcos(phi1)))
         at1  = phi1-(pitch(r*z_i,wtomega(ii))+Pitch_a)*pi/180.d0
c         Urel1= bu*(1.d0-a1)/dsin(phi1)
         Urel1= bu/((soli*CL11/(4.*Fcr1))/dtan(phi1)+dsin(phi1))
         Re1  = Urel1*chordl/nu*1.d0
         CL12 = CALCL(at1*180.d0/pi,r*z_i)
c         CL12 = CALCL1(at1*180.d0/pi,Re1)
                        
         CL21 = 4.d0*Fcr2*dsin(phi2)*(dcos(phi2)-lamdar*dsin(phi2))
     +        / (soli*(sin(phi2)+lamdar*cos(phi2)))
         at2  = phi2-(pitch(r*z_i,wtomega(ii))+Pitch_a)*pi/180.d0
c         Urel2= bu*(1.d0-a2)/dsin(phi2)
         Urel2= bu/((soli*CL21/(4.*Fcr2))/dtan(phi2)+dsin(phi2))
         Re2  = Urel2*chordl/nu*1.d0
         CL22 = CALCL(at2*180.d0/pi,r*z_i)
c         CL22 = CALCL1(at2*180.d0/pi,Re2)
                        
         CL31 = 4.d0*Fcr3*dsin(phi3)*(dcos(phi3)-lamdar*dsin(phi3))
     +        / (soli*(dsin(phi3)+lamdar*dcos(phi3)))
         at3  = phi3-(pitch(r*z_i,wtomega(ii))+Pitch_a)*pi/180.d0
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
         phi1 = datan((bu*(1-a1))/(omega*r*(1+b1)+0))
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
      subroutine sol_CL_CD_2(ii,Pitch_a,omega,chordl,r,bu,soli,a1,b1,CL
     +           ,CD,at,phi,Re1,me,nall)
      implicit none
      include 'dimen.h'
      integer*4 it,ii
      real*8 lamdar,omega,r,bu,soli
      real*8 a1,a2,a3,b1,b2,b3
      real*8 phi
      real*8 CL,CD
      real*8 at
      real*8 Urel1,Urel2,Urel3,Re1,Re2,Re3,chordl
      real*8, external :: pitch,CALCL,CALCD,CALCL1,CALCD1
      real*8 Pitch_a
	  
      lamdar=omega*r/bu

      a1=0
      a2=0
      
      phi=2./3.*datan(1.0/lamdar);		
cc      print*, phi*180/pi	   
      do it=1,20      
      a2=a1
cc      at  = phi-pitch(r*z_i,wtomega(1))*pi/180.d0	 
      at=phi-(pitch(r*z_i,wtomega(ii))+Pitch_a)*pi/180.d0		  
cc      print*, pitch(r*z_i,wtomega(1))
      CL = CALCL(at*180.d0/pi,r*z_i)
      CD  = CALCD(at*180.d0/pi,r*z_i)                        
      a1   = (1+(soli*(CL*dcos(phi)+CD*dsin(phi))/ 
     +       (4*(dsin(phi))**2))**(-1.0))**(-1.0)
      b1   = (1-a1)*soli*(CL*dsin(phi)-CD*dcos(phi))/ 
     +       (4*lamdar*(dsin(phi))**2)
                        
      phi=datan((1-a1)/(lamdar*(1+b1))) 
cc      if (abs(a1-a2).le.10**(-5)) exit	  
      enddo
	  
cc      print*, r*z_i,a1,at*180.d0/pi	
cc      print*, at*180.d0/pi	  
	  
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
      function calU(x)
      real*8 x,calU
     
      calU = -13.077*x**2 + 6.4431*x + 1.6279
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

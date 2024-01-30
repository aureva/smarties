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
      call CALBar   (u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
      end if

      if(turbine_model.eq.1)then
      call CALADMNR (u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
      end if
      if(turbine_model.eq.2)then
      call CALADMR  (u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
      end if
      if(turbine_model.eq.3)then
c      call CALALM   (u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
      end if

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
          CTU = u(wtx(ii)-i,wty(ii),k)
          if(CTU.lt.0.d0)CTU=0.d0
          dA = (0.005d0/z_i)*dz
          if(Zhub(ii)*0.94d0.lt.hz2)then
      	    hz2=Zhub(ii)*0.94d0
            dA =(0.005d0/z_i)*(hz2-hz1)
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
      real*8 CTU1,CTV,CTW,chordl,CTsoli
      real*8 CTV1,CTV2,CTUrel,CTFx,CTFt
      real*8 CTCL,CTCD,dang,Re,dA,RA,dL
      real*8 Ua,Va
      
      real*8 Pow(num_turbine),CTU(num_turbine),CTUnew(num_turbine),
     +       anglein(num_turbine)

      real*8, external :: chord,calang
      real*8, external :: itp2D,calU

c     Mahdi      
      real*8, dimension(num_turbine,Ny,Nz2):: gama
      real*8 rr(4),y1,z1
      integer*4 nn,jm,km,k1,j1
 	  
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
c      i = 2.d0*wtR(ii)/dx
	  i = 0
      uframe(ii,blade,N)=itp2D(u,wtx(ii)-i,blade_jj,blade_kk,me,nall)
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
         call MPI_SEND(vframe(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )
         call MPI_SEND(wframe(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )
         end do
      END IF
      IF (me>0) then
          call MPI_RECV(uframe(1,1,1),num_turbine*N_ang*N_blade,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
          call MPI_RECV(vframe(1,1,1),num_turbine*N_ang*N_blade,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
          call MPI_RECV(wframe(1,1,1),num_turbine*N_ang*N_blade,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
      END IF


      do ii=1,num_turbine
c      dr = wtr(ii)/dble(N_blade)
c      do N=N_blade_s,N_blade
c      r = dble(N-0.5)*dr
c      do blade=1,N_ang

      Ua=sum(uframe(ii,:,1))/dble(N_ang)
      Va=sum(vframe(ii,:,1))/dble(N_ang)	
      CTU(ii) =(Ua**2.+Va**2.)**0.5
      anglein(ii)=atan(Va/Ua)

cc      CTU = sum(uframe(ii,:,1))/dble(N_ang)
c        CTUnew(ii)=sum(uframe(ii,:,:))/(dble(N_ang)*N_blade)
c        if (t.eq.1) CTU(ii)=CTUnew(ii)
c        CTU(ii) = 1.*CTUnew(ii)+0.0*CTU(ii)
c         CTU(ii)=sum(uframe(ii,:,:))/(dble(N_ang)*N_blade)	   
c         CTU(ii)=sum(uframe(ii,:,1))/(dble(N_ang))	

      tmp1=1.d0

c      if(CTU*u_star.lt.8.)then
c      tmp1=(CTU*u_star-2.)/6.
c      end if
c      if(CTU*u_star.lt.2.)then
c      tmp1=0.0
c      end if

c      dA      = 0.5d0*(dble(N)**2-dble(N-1)**2)*(dr**2.)*dang
c      CTFx    = 0.5*(CTU**2)*CTT(ii)*dA/(dx*dy*dz)

      i=wtx(ii)
      do k=2,nzb+1
      do j=1,ny
c         if (ker_u(ii,blade,N,j,k).ge.1E-15)then
c        Fx (i,j,k)=Fx (i,j,k) + CTFx  * ker_u(ii,blade,N,j,k)
c        Fx (i,j,k)=Fx (i,j,k) + CTFx  * ker_u(ii,blade,N,j,k)*Ua/CTU
c		 Fy (i,j,k)=Fy (i,j,k) + CTFx  * ker_u(ii,blade,N,j,k)*Va/CTU
        Fx (i,j,k)=Fx (i,j,k) + 0.5*(CTU(ii)**2)*CTT(ii)*gama(ii,j,k)/dx
c         end if
      end do
      end do

c     end do
c     end do

      Pow(ii)=0.5*(CTU(ii)**3)*CTT(ii)
      end do
	  
        if (mod(t,c_count)==0 .and. me==0) then
        write(1524,5544) Pow(1)
c		,Pow(2),Pow(3),Pow(4),Pow(5),Pow(6),
c     +               Pow(7),Pow(8)
c     +               Pow(9),Pow(10),Pow(11),Pow(12),Pow(13),Pow(14),
c     +               Pow(15),Pow(16),Pow(17),Pow(18)
        write(*,*) 'power',  Pow(1)
c        write(*,*) CTU(1),sum(uframe(1,:,1))/dble(N_ang)
        write(*,*) CTU(1),CTUnew(1),anglein(1)*180/pi		
        endif

5544   format (24(f11.5,1x))
		
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
      integer*4 ii,i,j,k,t,blade,N,me_k

      real*8, dimension(Nx,Ny,Nz2)::u,v,w,fx,fy,fz,fax,faz,AveWt
      real*8, dimension(num_turbine,Ny,Nz2):: CTy,CTzu,CTzw,CTru,CTrw
      real*8, dimension(num_turbine,Ny,Nz2):: CTangx,CTangz
      real*8, dimension(num_turbine,Nx):: ker_x
      real*8 tmp1,tmp2,tmp3

      real*8, dimension(num_turbine,N_ang,N_blade)::
     +uframe,vframe,wframe,meframe,u1frame,v1frame,w1frame,me1frame
      
      real*8, dimension(num_turbine,N_ang,N_blade,Ny,Nz2):: ker_u,ker_w

      real*8 SUMCT,SUMa1,SUMdA
      
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
         call MPI_SEND(meframe(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         call MPI_SEND(uframe(1,1,1),num_turbine*N_ang*N_blade,
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
         call MPI_SEND(vframe(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )
         call MPI_SEND(wframe(1,1,1),num_turbine*N_ang*N_blade,
     +        MPI_DOUBLE_PRECISION,i,i,nall, ierr )
         end do
      END IF
      IF (me>0) then
          call MPI_RECV(uframe(1,1,1),num_turbine*N_ang*N_blade,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
          call MPI_RECV(vframe(1,1,1),num_turbine*N_ang*N_blade,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
          call MPI_RECV(wframe(1,1,1),num_turbine*N_ang*N_blade,
     +         MPI_DOUBLE_PRECISION,0,me,nall,status2,ierr )
      END IF

      do ii=1,num_turbine
      dr = wtr(ii)/dble(N_blade)
      do N=N_blade_s,N_blade
      r = dble(N-0.5)*dr
      do blade=1,N_ang

      CTU = uframe(ii,blade,N)
      CTV = vframe(ii,blade,N)
      CTW = wframe(ii,blade,N)

      chordl  = chord(r*z_i)/z_i
      CTsoli  = 3.0*chordl/(2.d0*pi*r)

      omega1  = +CTW*dcos(Tang(blade))/r
     +          +CTV*dsin(Tang(blade))/r
     +          +wtomega(ii)*(2.d0*pi/60)*(z_i/u_star)

      call sol_CL_CD(omega1,chordl,r,CTU,CTsoli,
     +     a1,b1,CTCL,CTCD,at1,CTphi,Re,me,nall)

      CTV1    =  CTU
      CTV2    = -CTW*dcos(Tang(blade))
     +            -CTV*dsin(Tang(blade))-omega1*r
      CTUrel  = CTU*(1-a1)/dsin(CTphi)

      dA      = 0.5d0*(dble(N)**2-dble(N-1)**2)*(dr**2.)*dang

      CTFx    = 0.5*(CTUrel**2)*CTsoli*dA/(dx*dy*dz)*
     +          (+CTCL*dcos(CTphi)+CTCD*dsin(CTphi))
      CTFt    = 0.5*(CTUrel**2)*CTsoli*dA/(dx*dy*dz)*
     +          (+CTCL*dsin(CTphi)-CTCD*dcos(CTphi))

      SUMdA =SUMdA+dA
      SUMa1 =SUMa1+a1*dA/(pi*wtR(ii)**2.d0)
      SUMCT =SUMCT+CTFx/(0.5d0*((9.d0/u_star)**2.)*dA/(dx*dy*dz))*
     +dA/(pi*wtR(ii)**2.d0)

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

      SUMa1=SUMa1/SUMdA*(pi*wtR(1)**2.d0)
      SUMCT=SUMCT/SUMdA*(pi*wtR(1)**2.d0)
      if(me.eq.0)then
      if(t.eq.1)then
      Open (unit=1095,file='output/check_a_CT.out',position='append')
      end if
      write(1095,*)SUMa1,SUMCT
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



      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function chord(r)
      real*8 r,chord

      real*8 c3_data(3)
! SWT-2.3-93 blade: chord = 0.058065*(1.0-r)/(1.0-0.032258)+0.017204
      data c3_data/0.032258d0,0.06d0,0.017204d0/

      r=r/46.5d0

      if (r < c3_data(1) .or. r > 1.0d0) then
         chord = 0.0d0
      else
         chord = c3_data(2)*(1.0d0-r)+c3_data(3)
      end if

      chord=chord*46.5d0

      r=r*46.5d0

c      if((r.lt.0.01).or.(r.gt.0.07))then
c      chord=-1.648759985E+07*r**6+3.707847138E+06*r**5
c     +      -3.252428628E+05*r**4+1.409761773E+04*r**3
c     +      -3.167275003E+02*r**2+3.476124426E+00*r**1 
c     +      -1.881364199E-05      
c      elseif((r.ge.0.01).and.(r.le.0.07))then
c      chord=-2.958333336E+06*r**6+7.279166675E+05*r**5
c     +      -7.131250008E+04*r**4+3.530208335E+03*r**3
c     +      -9.378916663E+01*r**2+1.249200004E+00*r**1
c     +      +8.159999868E-03
c      end if
	
c      if(r.le.0.07)then
c        chord=23.889*r**3 - 5.106*r**2 + 0.212*r**1 + 0.0123
c      else
c        chord=-44.667*r**2 + 5.6347*r**1 - 0.1654
c      end if
c      if(r.gt.0.075)then
c      chord=0.0d0
c      end if

c      if (r.ge.6.) then
c         chord = 3.3-2.4*(r-6.d0)/24.56d0
c      else 
c         chord = 3.3d0
c      end if
      
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine sol_CL_CD_1(omega,chordl,r,bu,soli
     +           ,a1,b1,CL11,CD,at1,phi1,Re1,me,nall)
      implicit none
      include 'dimen.h'
      integer*4 it
      real*8 lamdar,omega,r,bu,bu1,soli
      real*8 a1,a2,a3,b1,b2,b3
      real*8 phi1,phi2,phi3
      real*8 CL11,CL12,CL21,CL22,CL31,CL32,CD
      real*8 at1,at2,at3,res1,res2,Fcr1,Fcr2,Fcr3
      real*8 Urel1,Urel2,Urel3,Re1,Re2,Re3,chordl
      real*8, external :: pitch,CALCL,CALCD,CALCL1,CALCD1

      lamdar=-omega*r/bu
      phi1 = datan((bu)/(-omega*r))
                    
c     Fcr1 = 2/pi*dacos(dexp(-3.0*(wtR-r)/(2*r*dsin(phi1))))
      Fcr1 = 1.0

      at1  = phi1-pitch(r*z_i)*pi/180.d0
      Re1  = Urel1*chordl/nu*1.d0
      CL11 = CALCL(at1*180.d0/pi,r*z_i)
      CD   = CALCD(at1*180.d0/pi,r*z_i)


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

      lamdar=-omega*r/bu
                    
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
         at1  = phi1-pitch(r*z_i)*pi/180.d0
c         Urel1= bu*(1.d0-a1)/dsin(phi1)
         Urel1= bu/((soli*CL11/(4.*Fcr1))/dtan(phi1)+dsin(phi1))
         Re1  = Urel1*chordl/nu*1.d0
         CL12 = CALCL(at1*180.d0/pi,r*z_i)
c         CL12 = CALCL1(at1*180.d0/pi,Re1)
                        
         CL21 = 4.d0*Fcr2*dsin(phi2)*(dcos(phi2)-lamdar*dsin(phi2))
     +        / (soli*(sin(phi2)+lamdar*cos(phi2)))
         at2  = phi2-pitch(r*z_i)*pi/180.d0
c         Urel2= bu*(1.d0-a2)/dsin(phi2)
         Urel2= bu/((soli*CL21/(4.*Fcr2))/dtan(phi2)+dsin(phi2))
         Re2  = Urel2*chordl/nu*1.d0
         CL22 = CALCL(at2*180.d0/pi,r*z_i)
c         CL22 = CALCL1(at2*180.d0/pi,Re2)
                        
         CL31 = 4.d0*Fcr3*dsin(phi3)*(dcos(phi3)-lamdar*dsin(phi3))
     +        / (soli*(dsin(phi3)+lamdar*dcos(phi3)))
         at3  = phi3-pitch(r*z_i)*pi/180.d0
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
      function pitch(r)
      real*8 r,pitch,gamma

      r=r/46.5d0
! SWT-2.3-93 blade
!      if (r < 0.032258d0 .or. r > 1.0d0) then
!         gamma = 0.0d0
!      else
! linear distribution from 10 to 0 (tip pitch)
!         gamma = 10.4d0*(1.0d0-r)
!         gamma = 10.0d0*(1.0d0-r)/(1.0d0-0.032258d0)
!      end if



! SWT-2.3-93 blade from "Optimization of Wind Turbines with respect to Noise"
      if (r < 0.032258d0 .or. r > 1.0d0) then
      gamma = 0.0d0
      else
c      tmp1 = 13.0d0
c      gamma = (3.65228d0+r*(-10.04934d0+r*(9.36866d0+r*(-2.98057d0))))
c      gamma = dabs(gamma)
c      gamma = gamma*tmp1
c      if(r<0.38)gamma = 1.d0*tmp1

      tmp1 = 10.3d0
      gamma=-2.85600E+01*r**6 + 1.15008E+02*r**5 - 1.84639E+02*r**4 
     +     + 1.46526E+02*r**3 - 5.40413E+01*r**2 + 3.04717E+00*r**1
     +     + 2.65404E+00

      end if
      gamma = gamma*tmp1

      r=r*46.5d0

      pitch = gamma

CCCCCC old
c      pitch = -3.66065E+09*r**6 + 7.26503E+08*r**5 - 5.17506E+07*r**4 
c     +      +  1.67158E+06*r**3 - 3.26723E+04*r**2 + 4.24808E+02*r**1 
c     +      +  1.82914E+01

CCCCCC New
c      pitch = -14469267270.1562*r**6 + 3117024854.7461*r**5 
c     +        - 256580382.0576*r**4 + 10211248.1846*r**3 
c     +        - 208430.3097*r**2 + 1967.0664*r**1 + 16.4945

CCCCCC From Jim
CCCCCC BACK BACK BACK BACK BACK BACK BACK BACK BACK BACK BACK BACK BACK BACK BACK
c      pitch =-7.107662E+09*r**6 + 2.370132E+09*r**5 - 3.316227E+08*r**4 
c     +      + 2.479045E+07*r**3 - 1.033707E+06*r**2 + 2.221309E+04*r**1
c     +      - 1.695696E+02
c      if(r.le.0.03)pitch = 19.63

CCCCCC FACE FACE FACE FACE FACE FACE FACE FACE FACE FACE FACE FACE FACE FACE FACE
c      pitch =-1.891237E+07*r**4 + 4.145681E+06*r**3 - 3.320177E+05*r**2
c     +      + 1.113371E+04*r**1 - 1.096808E+02
c      if(r.le.0.03)pitch = 22.5

CCCCCC FACE BACK FACE BACK FACE BACK FACE BACK FACE BACK FACE BACK FACE BACK FACE BACK
c      pitch =-1.271396E+07*r**4 + 2.735155E+06*r**3 - 2.131390E+05*r**2
c     +      + 6.788930E+03*r**1 - 5.425028E+01
c      if(r.le.0.03)pitch = 21.0625

c      if (r.ge.6.) then
c         pitch = 12.-(r-6.d0)/3.d0
c      else 
c         pitch = 12.d0
c      end if

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function CALCL(x,r)
      include 'dimen.h'
      real*8 r,x,CALCL,CL1,CL2

      real*8 attack_angle,cl
! Airfoil NACA 63-618 for SWT-2.3-xx blade
      real*8 cl_data(63)
      data cl_data/0.056d0,0.114d0,0.172d0,0.230d0,0.288d0,0.346d0,
     & 0.404d0,0.461d0,0.519d0,0.577d0,0.635d0,0.692d0,0.750d0,0.807d0,
     & 0.865d0,0.922d0,0.980d0,1.037d0,1.094d0,1.151d0,1.208d0,1.265d0,
     & 1.322d0,1.376d0,1.430d0,1.483d0,1.536d0,1.589d0,1.642d0,1.695d0,
     & 1.605d0,1.614d0,1.637d0,1.659d0,1.679d0,1.693d0,1.702d0,1.706d0,
     & 1.705d0,1.703d0,1.702d0,1.698d0,1.695d0,1.687d0,1.678d0,1.665d0,
     & 1.650d0,1.632d0,1.611d0,1.586d0,1.558d0,1.525d0,1.488d0,1.446d0,
     & 1.400d0,1.349d0,1.293d0,1.231d0,1.165d0,1.093d0,1.017d0,0.934d0,
     & 0.934d0/

      integer ix
      real*8 tmp
      attack_angle=x

! Riso-A1-18 (Riso-P)
      if (attack_angle <= 20.0d0 .and. attack_angle >= -5.0d0) then
         tmp = attack_angle+5.0d0
         ix = tmp+1
         cl = cl_data(ix)+(cl_data(ix+1)-cl_data(ix))*(tmp+1.0d0-ix)
!!         cl = cl_data(ix)+(cl_data(ix+1)-cl_data(ix))*
!!     &        (attack_angle-((ix-1)-5.0d0))
      elseif (attack_angle > 20.0d0) then
         cl = cl_data(27)
      elseif (attack_angle < -5.0d0) then
         cl = cl_data(1)
      end if

      CALCL=cl

      if(x.le.13.5)then
      CALCL= 4.4650E-08*x**6 - 7.9983E-07*x**5 - 4.4934E-06*x**4 
     +     - 1.8503E-05*x**3 - 7.8678E-04*x**2 + 1.2223E-01*x**1
     +     + 4.1655E-01
      else
      CALCL= 1.9724E-09*x**6 + 5.7082E-07*x**5 - 4.0601E-05*x**4
     +     - 1.7913E-04*x**3 + 6.0053E-02*x**2 - 1.3857E+00*x**1
     +     + 1.0923E+01
      end if

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      if (x.lt.6) then
c      CALCL = 0.1117*x - 1E-16
c      elseif (x.ge.6.and.x.lt.20) then
c      CALCL = 0.0114*x + 0.6014
c      elseif (x.ge.20.and.x.lt.42.5) then
c      CALCL = -7.3319E-05*x**3 + 5.2197E-03*x**2 - 8.1195E-02*x 
c     +      + 9.5255E-01
c      else
c      CALCL = -1.7651E-04*x**2 + 7.9061E-03*x + 7.1406E-01
c      end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      if (x.le.12) then
c      CALCL = 2*pi*dsin(x*pi/180.d0)
c      else
c      CALCL = 2*pi*x*pi/180.d0
c      end if 
c      if (x.gt.16) CALCL = 2*pi*16*pi/180.d0 


c      if (x.lt.6) then
c      CALCL = 2*pi*dsin(x*pi/180.d0)
c      else
c      CALCL = -3.4226E-05*x**3 + 1.3393E-03*x**2 - 4.7202E-03*x + 6.5750E-01
c      CALCL = -3.4226E-05*x**3 + 1.3393E-03*x**2 - 4.7202E-03*x 
c     +      + 6.5750E-0

c      CALCL = 5.3571E-05*x**3 - 1.9643E-03*x**2 + 2.8429E-02*x 
c     +      + 5.3857E-01
c      end if



c       CALCL=-1.6351E-04*x**3.d0+6.5557E-05*x**2.d0
c     +       +1.2171E-01*x+3.4492E-01

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function CALCD(x,r)
      include 'dimen.h'
      real*8 r,x,CALCD,CD1,CD2

      real*8 attack_angle
      real*8 cd

! Airfoil NACA 63-618 for SWT-2.3-xx blade
      real*8 cd_data(63)
      data cd_data/0.0067d0,0.0066d0,0.0065d0,0.0065d0,0.0064d0,
     & 0.0064d0,0.0063d0,0.0062d0,0.0064d0,0.0065d0,0.0066d0,0.0067d0,
     & 0.0068d0,0.0070d0,0.0072d0,0.0075d0,0.0078d0,0.0081d0,0.0086d0,
     & 0.0090d0,0.0099d0,0.0102d0,0.0106d0,0.0111d0,0.0115d0,0.0119d0,
     & 0.0124d0,0.0129d0,0.0134d0,0.0139d0,0.0148d0,0.0154d0,0.0160d0,
     & 0.0167d0,0.0175d0,0.0182d0,0.0190d0,0.0199d0,0.0208d0,0.0218d0,
     & 0.0228d0,0.0239d0,0.0249d0,0.0261d0,0.0273d0,0.0286d0,0.0300d0,
     & 0.0314d0,0.0330d0,0.0346d0,0.0363d0,0.0381d0,0.0400d0,0.0420d0,
     & 0.0441d0,0.0464d0,0.0488d0,0.0513d0,0.0540d0,0.0568d0,0.0598d0,
     & 0.0630d0,0.0630d0/

      integer ix
      real*8 tmp
      real*8 cl

      attack_angle=x

! Riso-A1-18 (Riso-P)
      if (attack_angle <= 10.0d0 .and. attack_angle >= -5.0d0) then
         tmp = attack_angle+5.0d0
         ix = tmp+1
         cd = cd_data(ix)+(cd_data(ix+1)-cd_data(ix))*(tmp+1.0d0-ix)
!!         cd = cd_data(ix)+(cd_data(ix+1)-cd_data(ix))*
!!     &        (attack_angle-((ix-1)-5.0d0))
      elseif (attack_angle > 10.0d0) then
         cd = cd_data(27)

c         cd = 3.09612E-06*x**3 + 8.86736E-05*x**2
c     +      + 1.11250E-02*x**1 - 9.13290E-02

          cd = 2.7385E-07*x**4 - 2.3052E-05*x**3 + 7.3966E-04*x**2
     +       + 8.3038E-03*x**1 - 7.1226E-02


      elseif (attack_angle < -5.0d0) then
         cd = cd_data(1)
      end if

      CALCD = cd



c      CALCD  = 2.0*dsin(x/180*pi)*dsin(x/180*pi)

c      CALCD = 6.77083E-06*x**4. - 3.10417E-04*x**3. + 4.47292E-03*x**2.
c     +      - 5.45833E-03*x + 5.00000E-02
	  
c      CALCD = 1.3713E-06*x**3 - 4.6243E-05*x**2 + 1.4764E-02*x
c     +      + 4.6131E-02

c      CALCD = 4.3750E-04*x**2 + 5.7500E-03*x + 1.0000E-02

c      CALCD = 0.0001*x**2 + 0.0123*x + 0.01

c      if (x.lt.20) then
c      CALCD = 0.015*x + 0.03
c      elseif (x.ge.20.and.x.lt.42.5) then
c      CALCD = -3.5097E-05*x**2 + 4.2109E-02*x - 5.1393E-01
c      else 
c      CALCD = -0.0002*x**2 + 0.0386*x - 0.6145
c      end if


c      CALCD = 0.01

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

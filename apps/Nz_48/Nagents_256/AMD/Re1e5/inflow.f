CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine post_instant2_parallel(fid,ac,t,me,nall)
      implicit none
      include 'dimen.h'

      integer*4 fid,t,flag,k,Nrec,j,i
      real*8 aaa,aaaa,rand1
      real*8, dimension(nx,ny,nz2):: ac
      real*4, dimension(nxbe-nxbm+1,ny,nzb):: ac2
cc      real*4, dimension(1,ny,nzb):: ac2
      real*4 u_line(1,ny,nz2), u_plane(nz2)
ccm      real*8, dimension(nx,ny) :: t_flux

      character(len=80) path
      character(len=20) fm1,fm2,fm3,fme,fnum,ff,pp,ss,aa

      integer*4 irecl
      real*4 recx
      inquire(iolength=irecl) recx

      flag = 4100+1000*fid+me
cc      resubnum=0      
	  
      if(t.eq.1)then

      Nrec = (Nxbe-Nxbm+1)*Ny*Nzb*irecl
cc      Nrec = (1)*Ny*Nzb*irecl
c      write(*,*)fid,irecl,Nrec,flag

      if(me.lt.10)then
      write (fme,'(I1)') me
      elseif(me.lt.100)then 
      write (fme,'(I2)') me
      elseif(me.lt.1000)then 
      write (fme,'(I3)') me
      elseif(me.lt.10000)then 
      write (fme,'(I4)') me
      end if

      if(resubnum.lt.10)then
      write (fnum,'(I1)') resubnum
      elseif(resubnum.lt.100)then 
      write (fnum,'(I2)') resubnum
      elseif(resubnum.lt.1000)then 
      write (fnum,'(I3)') resubnum
      elseif(resubnum.lt.10000)then 
      write (fnum,'(I4)') resubnum
      end if
      fnum='0'
	  
      ff  ='unformatted'
      pp  ='append'
      ss  ='unknown'
c     aa  ='sequential'
      aa  ='direct'
c      fm1 ='instant_field/'
      fm2 ='_'
      fm3 ='.bin'

      if(fid.eq.1) fm1 ='instant_field/Su'
      if(fid.eq.2) fm1 ='instant_field/Sv'
      if(fid.eq.3) fm1 ='instant_field/Sw'
      if(fid.eq.4) fm1 ='instant_field/ST'

      path=trim(fm1)//trim(fnum)//trim(fm2)//trim(fme)//trim(fm3)

c     Open(flag,file=path,access=aa,status=ss,form=ff,position=pp)
      Open(flag,file=path,access=aa,status=ss,form=ff,recl=Nrec)

!      write(*,*)fid,irecl,Nrec,flag,path
      end if

c      u_line=0.0
c      u_plane=0.0
c      do k=2,Nzb+1
c         do j=1,Ny
c            do i=1,Nx
c               u_line(1,j,k)=u_line(1,j,k)+ac(i,j,k)			
c               u_plane(k)=u_plane(k)+ac(i,j,k)
c            end do
c         end do
c      end do		 
c      u_line=u_line/Nx
c      u_plane=u_plane/(Nx*Ny)	
	  
      do k=2,Nzb+1
      ac2(1:Nxbe-Nxbm+1,:,k-1)=ac(Nxbm:Nxbe,:,k)
cc       ac2(1,:,k-1)=u_plane(k)+(ac(Nxbe,:,k)-u_line(1,:,k))
cc       ac2(1,:,k-1)=ac(Nxbe,:,k)   
	   
	   
ccm      do j=1,ny 
ccm      CALL RANDOM_NUMBER (rand1)
ccm      rand1=0.02*(rand1-0.5)
ccm      write(*,*), rand1
ccm      ac2(1,j,k-1)=Ts1-
ccm     +      0.02*(1./vonk)*log(((k-1-0.5)+me*Nzb)*dz/(zo1/10./z_i))+
ccm     +      rand1	 
ccm      end do
      end do
	  
c      end if

c     write(flag) ac2
      write(flag,rec=t) ac2
      call flush(flag)

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine read_instant2_parallel(fid,ac,t,me,nall)
      implicit none
      include 'dimen.h'

      integer*4 fid,t,t1,flag,k,Nrec,j,i
      real*8, dimension(nxbe-nxbm+1,ny,nz2):: ac
      real*4, dimension(nxbe-nxbm+1,ny,nzb):: ac1
	  
cc      real*8, dimension(1,ny,nz2):: ac
cc      real*4, dimension(1,ny,nzb):: ac1

      real*8 aaa,aaaa
	  
      character(len=80) path
      character(len=20) fm1,fm2,fm3,fme,fnum,ff,pp,ss,aa

      integer*4 irecl
      real*4 recx
      inquire(iolength=irecl) recx

      flag = 4100+1000*fid+me

c... control to pick up the initial data 
c     t1=mod(t-1,50000)+1
c     if(t1.gt.25001) t1=50001-t1

      t1=t

      if(t.eq.1)then

      Nrec = (Nxbe-Nxbm+1)*Ny*Nzb*irecl
!      Nrec = (1)*Ny*Nzb*irecl
!      write(*,*)fid,irecl,Nrec

      if(me.lt.10)then
      write (fme,'(I1)') me
      elseif(me.lt.100)then 
      write (fme,'(I2)') me
      elseif(me.lt.1000)then 
      write (fme,'(I3)') me
      elseif(me.lt.10000)then 
      write (fme,'(I4)') me
      end if

      if(resubnum.lt.10)then
      write (fnum,'(I1)') resubnum
      elseif(resubnum.lt.100)then 
      write (fnum,'(I2)') resubnum
      elseif(resubnum.lt.1000)then 
      write (fnum,'(I3)') resubnum
      elseif(resubnum.lt.10000)then 
      write (fnum,'(I4)') resubnum
      end if
      fnum='0'

      ff  ='unformatted'
c     pp  ='append'
      pp  ='rewind'
      ss  ='unknown'
c     aa  ='sequential'
      aa  ='direct'
c     fm1 ='instant_field/'
      fm2 ='_'
      fm3 ='.bin'

      if(fid.eq.1) fm1 ='instant_field/Su'
      if(fid.eq.2) fm1 ='instant_field/Sv'
      if(fid.eq.3) fm1 ='instant_field/Sw'
      if(fid.eq.4) fm1 ='instant_field/ST'

      path=trim(fm1)//trim(fnum)//trim(fm2)//trim(fme)//trim(fm3)

c     Open(flag,file=path,access=aa,status=ss,form=ff,position=pp)
      Open(flag,file=path,access=aa,status=ss,form=ff,recl=Nrec)
      end if

c     read(flag) ac1
      read(flag,rec=t1) ac1

      do k=2,Nzb+1
      ac(:,:,k)=ac1(:,:,k-1)
      end do
	  
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      Subroutine inflow_t(u,v,w,uinf,vinf,winf,
!     +                    u2inf,v2inf,w2inf,t,me,nall)
!      Subroutine inflow_t(u,v,w,RHSx_BF,RHSy_BF,RHSz_BF,t,me,nall)
      Subroutine inflow_t(u,v,w,
     +		 uinf,vinf,winf,u2inf,v2inf,w2inf,
     +		 RHSx_BF,RHSy_BF,RHSz_BF,t,me,nall)	  
	 
      implicit none
      include 'dimen.h'
      
      integer*4 i,j,k,t,bx,fin,ii,jj,ttt,sn
      real*8, dimension(nx,ny,nz2)::u,v,w  

      real*8, dimension(Nxbe-Nxbm+1,ny,nz2)::tmpu2,tmpv2,tmpw2
      real*8, dimension(Nxbe-Nxbm+1,ny,nz2)::u_f,v_f,w_f	
      real*8, dimension(Nxbe-Nxbm+1,ny,nz2)::u_b,v_b,w_b	  
	  
	  
      real*8, dimension(nxx,ny,nz2)::uinf,vinf,winf
      real*8, dimension(nxx,ny,nz2)::u2inf,v2inf,w2inf	
      real*8 avg_u_inf(nz2),avg_u2_inf(nz2)
      real*8 avg_v_inf(nz2),avg_v2_inf(nz2)
      real*8 avg_w_inf(nz2),avg_w2_inf(nz2)	 	  
    
      real*8, dimension(nx,ny,nz2):: RHSx_BF,RHSy_BF,RHSz_BF
      real*8  factor,cfrdmp	

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call read_instant2_parallel(1,tmpu2,t,me,nall)
      call read_instant2_parallel(2,tmpv2,t,me,nall)
      call read_instant2_parallel(3,tmpw2,t,me,nall)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do k=2,nzb+1
      do j=1,ny
      do i=1,Nxbe-Nxbm+1
	  
      u_f(i,j,k)=tmpu2(i,j,k)
      v_f(i,j,k)=tmpv2(i,j,k)
      w_f(i,j,k)=tmpw2(i,j,k)
!      u_f(i,j,k)=8.0d0
!      v_f(i,j,k)=0.0d0
!      w_f(i,j,k)=0.0d0  
	  
      end do
      end do
      end do	

      do k=2,nzb+1
        avg_u_inf(k) = sum(uinf(:,:,k))/(Nxx*Ny)		 
        avg_v_inf(k) = sum(vinf(:,:,k))/(Nxx*Ny)	
        avg_w_inf(k) = sum(winf(:,:,k))/(Nxx*Ny)	
        avg_u2_inf(k) = sum(u2inf(:,:,k))/(Nxx*Ny)	
        avg_v2_inf(k) = sum(v2inf(:,:,k))/(Nxx*Ny)
        avg_w2_inf(k) = sum(w2inf(:,:,k))/(Nxx*Ny) 
      end do	  

      do k=2,nzb+1
        do j=1,Ny
        if (avg_u2_inf(k) .gt. 0.005) then
        u_b(:,j,k) = avg_u_inf(k)+(u_f(:,j,k)-uinf(Nxx:Nxx,j,k))
!     +                     *(avg_u2_inf(k)/u2inf(Nxx:Nxx,j,k))**0.5
        else	 
        u_b(:,j,k) = avg_u_inf(k)+(u_f(:,j,k)-uinf(Nxx:Nxx,j,k))		
        endif		
        if (avg_v2_inf(k) .gt. 0.005) then
        v_b(:,j,k) = avg_v_inf(k)+(v_f(:,j,k)-vinf(Nxx:Nxx,j,k))
!     +                     *(avg_v2_inf(k)/v2inf(Nxx:Nxx,j,k))**0.5	
        else	 
        v_b(:,j,k) = avg_v_inf(k)+(v_f(:,j,k)-vinf(Nxx:Nxx,j,k))	
        endif		
        if (avg_w2_inf(k) .gt. 0.005) then
        w_b(:,j,k) = avg_w_inf(k)+(w_f(:,j,k)-winf(Nxx:Nxx,j,k))
!     +                     *(avg_w2_inf(k)/w2inf(Nxx:Nxx,j,k))**0.5
        else	 
        w_b(:,j,k) = avg_w_inf(k)+(w_f(:,j,k)-winf(Nxx:Nxx,j,k))		
        endif			
        end do
      end do	
  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      sn=45
!      do j=1,ny-sn
!         u_b(:,j+sn,:)=u_f(:,j,:)
!         v_b(:,j+sn,:)=v_f(:,j,:)
!         w_b(:,j+sn,:)=w_f(:,j,:)
!      end do
!      do j=ny-sn+1,ny 
!         u_b(:,ny-j+1,:)=u_f(:,ny-sn+ny-j+1,:)
!         v_b(:,ny-j+1,:)=v_f(:,ny-sn+ny-j+1,:)
!         w_b(:,ny-j+1,:)=w_f(:,ny-sn+ny-j+1,:)
!      end do

	   u_f=u_b
	   v_f=v_b
	   w_f=w_b

      do k=2,nzb+1
        do j=1,Ny
        u(Nxbm:Nxbe,j,k) = u_f(:,j,k)	
        v(Nxbm:Nxbe,j,k) = v_f(:,j,k)
        w(Nxbm:Nxbe,j,k) = w_f(:,j,k)			
        end do
      end do	
	  
      RHSx_BF=0.0d0	  
      RHSy_BF=0.0d0	 
      RHSz_BF=0.0d0	
      cfrdmp = 1./(rlx_time*u_star/z_i)	 
	  
      do k=2,nzb+1
       do i=Nxbs,Nxbm	
         factor = 0.5*(1.-cos(pi*1.*(i-Nxbs)/(Nxbm-Nxbs)))
         do j=1,Ny		 
!		   RHSx_BF(i,j,k)=cfrdmp*factor*(u(i,j,k)-u_f(1,j,k))
!		   RHSy_BF(i,j,k)=cfrdmp*factor*(v(i,j,k)-v_f(1,j,k))
!		   RHSz_BF(i,j,k)=cfrdmp*factor*(w(i,j,k)-w_f(1,j,k))	

          u(i,j,k)=u(Nxbs,j,k)+factor*(u_f(1,j,k)-u(Nxbs,j,k))
!          u(i,j,k)=u(i,j,k)+(factor*(1-factor))*(u(i,j,k)-uinf(i,j,k))	
          v(i,j,k)=v(Nxbs,j,k)+factor*(v_f(1,j,k)-v(Nxbs,j,k))	
!          v(i,j,k)=v(i,j,k)+(factor*(1-factor))*(v(i,j,k)-vinf(i,j,k))	
          w(i,j,k)=w(Nxbs,j,k)+factor*(w_f(1,j,k)-w(Nxbs,j,k))	
!          w(i,j,k)=w(i,j,k)+(factor*(1-factor))*(w(i,j,k)-winf(i,j,k))
		  		   
         end do
       end do
      end do
	  
!      if (me.eq.0) print*,avg_u2_inf(2) 	  
      return

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



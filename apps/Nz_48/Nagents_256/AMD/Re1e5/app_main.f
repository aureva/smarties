!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     c    This program is developed for parallel            c
!     c    computation of scalar transport in ABL            c
!     c                                                      c
!     c					    July, 2002       c
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      
      module app_main_module
 
      implicit none

      public app_main
 
      contains 

      function app_main(smarties_comm, f_mpicomm) result(result_value)
     +      bind(c, name= 'app_main')

      use, intrinsic :: iso_c_binding
      use smarties
      use smarties_stat

      implicit none
      include 'dimen.h'
      integer*4 t,i,j,k,ct,cst,pt,kk,t01,t02,ii,i1,iii,jj,tct         
      integer*4 rule,fst,ilow(nx),kstart,ios,
     +          ttt,oldsteps,RL_init
      integer*4 u_L10,u_L25,u_L50,u_L100,u_L200

      real*8, dimension(nx,ny,nz2):: u,v,w,P,cx,cy,cz,RHSx,RHSy,RHSz,
     +     RHSx_f,RHSy_f,RHSz_f,dudx,dudy,dudz,dvdx,dvdy,dvdz,
     +     dwdx,dwdy,dwdz,txx,txy,txz,tyy,tyz,tzz,divtx,divty,divtz,
     +     dpdx,dpdy,dpdz,dtdx,dtdy,Theta,Beta,qx,qy,qz,
     +     RHS_T,RHS_Tf,S,dtdz,ddtzz,ESGS3D,ET3D,EQ3D,
     +     q,RHS_Q,RHS_Qf,dqdx,dqdy,dqdz,sgs_q1,sgs_q2,sgs_q3,
     +     rdmp,u_n,v_n
	 
      real*8, dimension(nxx,ny,nz2)::uinf,vinf,winf,u2inf,v2inf,w2inf
	  
      real*8, dimension(nx2,ny2,nz2):: u_m,v_m,w_m

      real*8,dimension (nx,ny):: zo,t_s,ustar,M,t_flux,q_flux,atxz_s,
     +     aqz_s,aqsurf,coolrate,q_s,qcoolrate,uh,vh,ts_avg

      real*8,dimension(anx,nz2) :: au,av,aw,u2,v2,w2,w3,atxz,
     +     atyy,atxx,atyz,atzz,atxy,p2,auw,avw,auv,ap,adudz,adudx,
     +     abeta1,at,t2,t3,asgs_t1,asgs_t2,asgs_t3,aut,avt,awt,adtdx,
     +     adtdy,adtdz,aPr,abeta2,e,aCs2,aCs,aCs2Pr,adwdz,adwdx,
     +     advdz,aESGS,aET,aq,q2,q3,asgs_q1,asgs_q2,asgs_q3,auq,avq,awq,
     +     adqdx,adqdy,adqdz,aSc,aCs2Sc,abeta3,aEQ,
     +     aFSu1,aFSu2,aFSu3,aFSu4,aFSv1,aFSv2,aFSv3,aFSv4
      
      real*8,dimension(nx/2+1,nz2) :: specu,specv,specw,spect,specq

      real*8 sum_s,rmsdivvel,rmsdivve,kee,ke,ave_t,ave_q,
     +     fi(nx,ny,2),fi_h(nx,ny,2),Psi(nx,ny,2),Psi0(nx,ny,2),wgx,
     +     oldtime,force(nz2),cfrdmp,ztemp,Tbase(nz2),Qbase(nz2)

      real*8 q_s2(nhrs),t_s2(nhrs),Ugeo_g,Vgeo_g,Uadv_d,Vadv_d,
     +     u_lv,ave_txz,ave_tyz,Tadv_d,Qadv_d,
     +     ave_sgst3,ave_sgsq3,Ptime,ave_u,ave_v

      real*8,dimension(nz2) :: Tadv,Qadv,Uadv,Vadv,Ugeo,Vgeo,uu2,
     +                         Ugeo_n,Vgeo_n	  

      real*8,dimension(nx,ny,nz2-2):: u_frame,v_frame,w_frame,t_frame,
     +                             q_frame	 

      real*8  time1,time2,umean,vmean,wmean	 
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      real*8,dimension(nx,ny,nz2)::
     +     a1_old,b1_old,c1_old,d1_old,e1_old,
     +     a2_old,b2_old,c2_old,d2_old,e2_old,
     +     a4_old,b4_old,c4_old,d4_old,e4_old,
     +     a8_old,b8_old,c8_old,d8_old,e8_old,
     +     a5_old,b5_old,c5_old,d5_old,e5_old,
     +     a9_old,b9_old,c9_old,d9_old,e9_old,
     +     Pr2,Cs2,Sc2,beta1,beta2,beta3
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      real*8 uave,tave,ubar,wvar	  
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      real*8 CFLx,CFLy,CFLz,CFLx1,CFLy1,CFLz1,tau_total(nz),t_total(nz2)
      real*8 alpha,h_bl,us_top,t_t(nz2),tau_t(nz),us_top2,xx,
     +       u_hub,v_hub,alpha_n
      real*8 umean_top,tmean_top,u_p,u_l,ncs2g,ncs2l,ncs2g1,ncs2l1 
      integer*4 me_k,l_kk, kk_m,k_bl, k_top, iT, jT, kT

      real*8, dimension(nx,ny,nz2)::
     +       fx,fy,fz,bafx,bafy,bafz
     +      ,fax,faz,AveWt,Bfx,Bfy,Bfz

      real*8, dimension(nx,ny,nz2)::
     +       bau1,bav1,baw1,bap1,bau2,bav2,baw2,bap2,bauv,bauw,bavw
     +      ,batxx,batxy,batxz,batyy,batyz,batzz,badpdx,badpdy,badpdz
     +      ,badudx,badudy,badudz,badvdx,badvdy,badvdz,badwdx,badwdy
     +      ,badwdz,babeta1,baCs,baCs2
     +      ,bafax,bafaz,bAveWt,baufx,bavfy,bawfz
     +      ,baTt,baTpx,baTpy,baTpz,baTsgsx,baTsgsy,baTsgsz	 
     +      ,baS11,baS22,baS33,baS12,baS13,baS23
     +      ,baeps11,baeps22,baeps33,baeps12,baeps13,baeps23

      real*8, dimension(nx,ny,nz2)::
     +        bat,bat2,baqx,baqy,baqz,baut,bavt,
     +        bawt,badtdx,badtdy,badtdz,baPr2,baBetaH
	 
      real*8, dimension(nx,ny,nz2)::
     +        baqm,baqm2,baqmx,baqmy,baqmz,bauqm,bavqm,
     +        bawqm,badqmdx,badqmdy,badqmdz,baSc2,baBetaqm
	 
      real*8, dimension(nx,ny,nz2):: RHSx_BF,RHSy_BF,RHSz_BF,RHS_T_BF	 
 
      type(c_ptr), intent(in),value :: smarties_comm
      integer(c_int), intent(in), value :: f_mpicomm
      integer*4 :: num_actions
      integer*4:: state_size,nx_agents,ny_agents
      integer(c_int)  , dimension(6),  target :: b_observable
      real(c_double), dimension(1), target ::
     +                     upper_action_bound,lower_action_bound
      logical(c_bool)   :: bounded
      real*8,dimension(nx,ny) :: tauw
      CHARACTER(len=255) :: path
      integer*4 :: me_test,nprocs_test,train
      integer*8 plan_f,plan_b
      integer*8 plan_ff,plan_bb
      integer*4 :: flag_RL,id,RL_stop
      integer(c_int) :: result_value
      real*8, dimension(:),allocatable :: rewards_old
      integer*4, dimension(nx,ny,4) :: neighbor
      real*8, dimension(:,:),allocatable :: actions_2D
      real*8, dimension(:,:,:),allocatable :: states_2D
      


      flag_RL = 1
      num_actions = 1
      state_size = 6
      nx_agents = 16
      ny_agents = 16
      RL_stop=9
      train = 0

      allocate(rewards_old(ny_agents))
      allocate(actions_2D(nx_agents,ny_agents))
      allocate(states_2D(nx_agents,ny_agents,2))
!cccccccccccccccccccccccccccccccccccccccccccccccccc

      call MPI_INIT( ierr )
 
      call MPI_COMM_RANK(f_mpicomm, me, ierr)
      call MPI_COMM_SIZE(f_mpicomm, nprocs, ierr)
      nall=f_mpicomm

      write(*,*) 'rank #', me, ' of ', nprocs,' is alive in Fortran'
      if (me .LT. nprocs) then

      CALL getcwd(path)
      WRITE(*,*) TRIM(path)
      CALL chdir
     +("/home/aurelien.vadrot/TOOLS/smarties/apps/model3/Nz_"
     + //trim(str(nz))//
     + "/Nagents_"//trim(str(nx_agents*ny_agents)) 
     + //"/AMD/N256/Re1e5")
      CALL getcwd(path)
      WRITE(*,*) TRIM(path)

      if (train == 1) then
      if (me == 0) then 
        call execute_command_line ("rm -rf output/*.bin")
        call execute_command_line ("rm -rf output/*.dat")
        call execute_command_line (trim("cp /home/aurelien.vadrot/")//
     + trim("CFD_Codes/LES_CTRSP2022/Nz_48/Viscosity_on/XRWM/AMD/")//
     + trim("Re1e5/input/*.out input/"))
       print *, trim("cp /home/aurelien.vadrot/")//
     + trim("CFD_Codes/LES_CTRSP2022/Nz_48/Viscosity_on/XRWM/AMD/")//
     + trim("Re1e5/input/*.out input/")
      end if
      end if

      if(me.eq.0)then
         open(unit=61,file='RL_init',status='old',iostat=ios)
         if(ios.eq.0)then
           read(61,*) RL_init
         elseif(resub.ne.0)then
           RL_init = 0
         end if
         print *,"RL_init", RL_init
      end if

      call MPI_Bcast(RL_init,1,MPI_INTEGER,0,nall,ierr)

!      call MPI_COMM_RANK( MPI_COMM_WORLD, me, ierr )
!      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )
!      nall=MPI_COMM_WORLD

!cccc  calculated parameters ccccc
      nzb = nz/nprocs
      dtl = dt*cs_count
      dx = 2.0d0*Pi/Nx 
      dy=2.0d0*Pi/Ny/l_r
      dz=(L_z/z_i)/(Nz-1)
      delta=fgr*(dx*dy*dz)**(1./3.)  
      idz=1.d0/dz
      idx=1.d0/dx
      idy=1.d0/dy
      inxny=1.d0/(Nx*Ny)
      inx2ny2=1.d0/(Nx2*Ny2)
      fc=f_c*10.**(-4)*z_i/u_star
      g_hat=9.81d0*(z_i/(u_star**2))

      !Read bias and weights files for data-driven Wall model
 
!cccc  Check for Resubmit for initu inits
      if(me.eq.0)then
         open(unit=61,file='resubmit',status='old',iostat=ios)
         if(ios.eq.0)then
            read(61,*) resubnum
         elseif(resub.ne.0)then
            print *, 'resubmit is missing and resub/=0 
     +           in dimen, exiting'
            goto 999
         else
            resubnum=-1
         endif   
      endif
      call MPI_Bcast(resubnum,1,MPI_INTEGER,0,nall,ierr)

      if(resubnum.lt.0)then
         initu=inituA
         inits=initsA
      elseif(resub.eq.resubnum)then
         initu=0
         inits=0
      else
         initu=1
         inits=1
      endif
      if(me.eq.0)then
      write(*,*) 'initu,inits=',initu,',',inits
      endif
!cccccccccccccccccccccccccccccccc
      
      write(*,*) 'NZB =',nzb,'NZ2=',nz2

!     .. Output the error if the array size is improperly posed

      IF (nz2<nzb+2.or.nx2.ne.3*nx/2.or.ny2.ne.3*ny/2) then
         print *,'error, array size too small, terminate'
         goto 999
      ENDIF

      write(*,*) ' My rank =', me,' Tot procs=', nprocs
      write(*,*) ' NZB =', nzb,'INITU =', initu
      
      cst=0

      u(:,:,:) = 0d0
      v(:,:,:) = 0d0
      w(:,:,:) = 0d0
      
      u_m(:,:,:) =0d0
      v_m(:,:,:) =0d0
      w_m(:,:,:) =0d0


!     ... Open out files      

      IF (me==nprocs-1) then

!C     ... me=0 only opens vel.out, q.out and t.out      

         Open (unit=84,file='input/vel.out',status='unknown')
!c     ... file to initilize lagrangian averaging ...
         Open (unit=20,file='input/a1_e1.bin',status='unknown',
     +        form='unformatted')
         Open (unit=21,file='input/a2_e2.bin',status='unknown',
     +        form='unformatted')
!c     ... end file for lag ...

         Open (unit=25,file='input/RHS_mom.bin',status='unknown',
     +        form='unformatted')

         IF (S_Flag.eq.1) then
            Open (unit=108,file='input/t.out',status='unknown')
!c     ... file to initilize lagrangian averaging ...
            Open (unit=22,file='input/a4_e4.bin',status='unknown',
     +           form='unformatted')
            Open (unit=23,file='input/a8_e8.bin',status='unknown',
     +           form='unformatted')
!c     ... end file for lag ...
            
            Open (unit=26,file='input/RHS_T.bin',status='unknown',
     +           form='unformatted')
            
         ENDIF

         IF (Q_Flag.eq.1) then
            Open (unit=109,file='input/q.out',status='unknown')
!c     ... file to initilize lagrangian averaging ...
            Open (unit=24,file='input/a5_e5.bin',status='unknown',
     +           form='unformatted')
            Open (unit=27,file='input/a9_e9.bin',status='unknown',
     +           form='unformatted')
!c     ... end file for lag ...
            
            Open (unit=28,file='input/RHS_Q.bin',status='unknown',
     +           form='unformatted')
            
         ENDIF         
         
      ENDIF

      IF(ME.EQ.nprocs-1) THEN

!C     ... that means the processor nprocs-1 is used
!C     ... as the I/O processor      
         
!cc         read (84,*) oldsteps,oldtime
         read (84,*) 
         read (84,*) 		 
         read (84,*)	

         print*, 's',oldsteps,oldtime,zref
		 
         if(resubnum.ge.0.and.initu.eq.1)then
            nrsub=oldsteps
         else
            nrsub=0
         endif

!c     Runesha      
!c     Read only the portion of each processor and send it out
!c     to the processor
         IF(nprocs.gt.1)THEN
            
            DO II=0,nprocs-2
               
               do k=1,NZB
                  do j=1,Ny
                     do i=1,Nx
                        read (84,*) u(i,j,k),v(i,j,k),w(i,j,k)
                        
                        u(i,j,k)=u(i,j,k)-Ugal
                        v(i,j,k)=v(i,j,k)-Vgal
                     end do
                  end do
               end do
               
               call MPI_SEND(u(1,1,1),(nzb)*nx*ny, 
     +              MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
               call MPI_SEND(v(1,1,1),(nzb)*nx*ny, 
     +              MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
               call MPI_SEND(w(1,1,1),(nzb)*nx*ny, 
     +              MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
               
            END DO

 
         ENDIF
         
!c     read for himself

         do k=2,NZB+1
            do j=1,Ny
               do i=1,Nx
                  read (84,*) u(i,j,k),v(i,j,k),w(i,j,k)
                       u(i,j,k)=u(i,j,k)-Ugal
                       v(i,j,k)=v(i,j,k)-Vgal
               end do
            end do
         end do

!c     ... use the frames and temp0 to hold arrays to pass out ...
         IF(initu.eq.1)THEN
!            read (20) u_frame,v_frame,w_frame,t_frame,q_frame
			
            IF(nprocs.gt.1)THEN
               DO II=0,nprocs-2

               do k=1,NZB
                  do j=1,Ny
                     do i=1,Nx
			      read (20)    a1_old(i,j,k),b1_old(i,j,k),
     + 			               c1_old(i,j,k),d1_old(i,j,k),
     +                         e1_old(i,j,k)  
                     end do
                  end do
               end do
			   
                  call MPI_SEND(a1_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(b1_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(c1_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(d1_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(e1_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  
               ENDDO
            ENDIF
!c     read for himself
            do k=2,NZB+1
               do j=1,Ny
                  do i=1,Nx
                     
                  read (20) a1_old(i,j,k),b1_old(i,j,k),
     + 			            c1_old(i,j,k),d1_old(i,j,k),
     +                      e1_old(i,j,k) 
                     
                  enddo
               enddo
            enddo

!            read (21) u_frame,v_frame,w_frame,t_frame,q_frame
            IF(nprocs.gt.1)THEN
               DO II=0,nprocs-2
                  
               do k=1,NZB
                  do j=1,Ny
                     do i=1,Nx			  
                  read (21) a2_old(i,j,k),b2_old(i,j,k),
     + 			            c2_old(i,j,k),d2_old(i,j,k),
     +                      e2_old(i,j,k)  
                     end do
                  end do
               end do
			   
                  call MPI_SEND(a2_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(b2_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(c2_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(d2_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(e2_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  
               ENDDO
            ENDIF
!c     read for himself
            do k=2,NZB+1
               do j=1,Ny
                  do i=1,Nx
                    
                  read (21) a2_old(i,j,k),b2_old(i,j,k),
     + 			            c2_old(i,j,k),d2_old(i,j,k),
     +                      e2_old(i,j,k) 
                     
                  enddo
               enddo
            enddo
            
!            read (25) u_frame,v_frame,w_frame
            IF(nprocs.gt.1)THEN
               DO II=0,nprocs-2

               do k=1,NZB
                  do j=1,Ny
                     do i=1,Nx	   					 
                          read (25) RHSx(i,j,k),RHSy(i,j,k),
     + 			                    RHSz(i,j,k) 
                     end do
                  end do
               end do	 

			   
                  call MPI_SEND(RHSx(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(RHSy(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(RHSz(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )


               ENDDO
            ENDIF
!c     read for himself
            do k=2,NZB+1	
               do j=1,Ny
                  do i=1,Nx

                  read (25) RHSx(i,j,k),RHSy(i,j,k),
     + 			            RHSz(i,j,k) 

                  enddo
               enddo
            enddo

         endif
         
!C     ... IF scalar is computed

         IF(S_Flag.eq.1) then
               
            print *,'Reading Initial temp field' 

            IF(nprocs.gt.1)THEN
               DO II=0,nprocs-2
                  
                  do k=1,NZB
                     do j=1,Ny
                        do i=1,Nx
                           read (108,*) theta(i,j,k)
!                           theta(i,j,k)=293.0
                        end do
                     end do
                  end do
                  
                  call MPI_SEND(theta(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  
               END DO
            ENDIF
            
!c     read for himself
            
            do k=2,NZB+1
               do j=1,Ny
                  do i=1,Nx
                     
                     read (108,*) theta(i,j,k)
!                           theta(i,j,k)=293.0					 
                     
                  end do
               end do
            end do

!c     ... use the frames and temp0 to hold arrays to pass out ...
            IF(inits.eq.1)then
!               read (22) u_frame,v_frame,w_frame,t_frame,q_frame
               IF(nprocs.gt.1)THEN
                  DO II=0,nprocs-2

               do k=1,NZB
                  do j=1,Ny
                     do i=1,Nx				  
                         read (22) a4_old(i,j,k),b4_old(i,j,k),
     + 			                   c4_old(i,j,k),d4_old(i,j,k),
     +                             e4_old(i,j,k)  
                     end do
                  end do
               end do
			   
                  call MPI_SEND(a4_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(b4_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(c4_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(d4_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(e4_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                     
                  ENDDO
               ENDIF
!c     read for himself
               do k=2,NZB+1
                  do j=1,Ny
                     do i=1,Nx
                        
                  read (22) a4_old(i,j,k),b4_old(i,j,k),
     + 			            c4_old(i,j,k),d4_old(i,j,k),
     +                      e4_old(i,j,k) 
                        
                     enddo
                  enddo
               enddo
               
!               read (23) u_frame,v_frame,w_frame,t_frame,q_frame
               IF(nprocs.gt.1)THEN
                  DO II=0,nprocs-2

               do k=1,NZB
                  do j=1,Ny
                     do i=1,Nx				  
                          read (23) a8_old(i,j,k),b8_old(i,j,k),
     + 			                    c8_old(i,j,k),d8_old(i,j,k),
     +                              e8_old(i,j,k)  
                     end do
                  end do
               end do
			   
                  call MPI_SEND(a8_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(b8_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(c8_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(d8_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(e8_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                     
                  ENDDO
               ENDIF
!     read for himself
               do k=2,NZB+1
                  do j=1,Ny
                     do i=1,Nx
                  read (23) a8_old(i,j,k),b8_old(i,j,k),
     + 			            c8_old(i,j,k),d8_old(i,j,k),
     +                      e8_old(i,j,k)
                     enddo
                  enddo
               enddo
               
!               read (26) t_frame
               IF(nprocs.gt.1)THEN
                  DO II=0,nprocs-2

               do k=1,NZB
                  do j=1,Ny
                     do i=1,Nx				  
                          read (26) RHS_T(i,j,k)
                     end do
                  end do
               end do
			   
                  call MPI_SEND(RHS_T(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )


                  ENDDO
               ENDIF
!c     read for himself
               do k=2,NZB+1
                  do j=1,Ny
                     do i=1,Nx
                         read (26) RHS_T(i,j,k)
                     enddo
                  enddo
               enddo
         
           endif
            
         ENDIF

!c....If a second scalar is computed

         IF(Q_Flag.eq.1) then
               
            print *,'Reading Initial q field' 

            IF(nprocs.gt.1)THEN
               DO II=0,nprocs-2
                  
                  do k=1,NZB
                     do j=1,Ny
                        do i=1,Nx
                           read (109,*) q(i,j,k)
                        end do
                     end do
                  end do
                  
                  call MPI_SEND(q(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  
               END DO
            ENDIF
            
!c     read for himself
            
            do k=2,NZB+1
               do j=1,Ny
                  do i=1,Nx
                     
                     read (109,*) q(i,j,k)
                     
                  end do
               end do
            end do

!c     ... use the frames and temp0 to hold arrays to pass out ...
            IF(inits.eq.1)then
!               read (24) u_frame,v_frame,w_frame,t_frame,q_frame
               IF(nprocs.gt.1)THEN
                  DO II=0,nprocs-2

               do k=1,NZB
                  do j=1,Ny
                     do i=1,Nx				  
                         read (24) a5_old(i,j,k),b5_old(i,j,k),
     + 			                   c5_old(i,j,k),d5_old(i,j,k),
     +                             e5_old(i,j,k)  
                     end do
                  end do
               end do
			   
                  call MPI_SEND(a5_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(b5_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(c5_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(d5_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(e5_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                     
                  ENDDO
               ENDIF
!     read for himself
               do k=2,NZB+1
                  do j=1,Ny
                     do i=1,Nx
                        
                  read (24) a5_old(i,j,k),b5_old(i,j,k),
     + 			            c5_old(i,j,k),d5_old(i,j,k),
     +                      e5_old(i,j,k) 
                        
                     enddo
                  enddo
               enddo
               
!               read (27) u_frame,v_frame,w_frame,t_frame,q_frame
               IF(nprocs.gt.1)THEN
                  DO II=0,nprocs-2

               do k=1,NZB
                  do j=1,Ny
                     do i=1,Nx				  
                         read (27) a9_old(i,j,k),b9_old(i,j,k),
     + 			                   c9_old(i,j,k),d9_old(i,j,k),
     +                             e9_old(i,j,k)  
                     end do
                  end do
               end do
			   
                  call MPI_SEND(a9_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(b9_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(c9_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(d9_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                  call MPI_SEND(e9_old(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
                     
                  ENDDO
               ENDIF
!c     read for himself
               do k=2,NZB+1
                  do j=1,Ny
                     do i=1,Nx
                          read (27) a9_old(i,j,k),b9_old(i,j,k),
     + 			                    c9_old(i,j,k),d9_old(i,j,k),
     +                              e9_old(i,j,k)  
                     enddo
                  enddo
               enddo
               
!               read (28) t_frame
               IF(nprocs.gt.1)THEN
                  DO II=0,nprocs-2

               do k=1,NZB
                  do j=1,Ny
                     do i=1,Nx					  
                          read (28) RHS_Q(i,j,k)
                     end do
                  end do
               end do
			   
                  call MPI_SEND(RHS_Q(1,1,1),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )

                  ENDDO
               ENDIF
!c     read for himself
               do k=2,NZB+1
                  do j=1,Ny
                     do i=1,Nx
                          read (28) RHS_Q(i,j,k)
                     enddo
                  enddo
               enddo
         
           endif
            
        ENDIF

!c     close files so that me=0 can write to them in the end

         close(unit= 84)
         close(unit= 20)
         close(unit= 21)
         close(unit= 25)
         if(S_flag.eq.1)then
            close(unit=108)
            close(unit= 22)
            close(unit= 23)
            close(unit= 26)
         endif
         if(Q_flag.eq.1)then
            close(unit=109)
            close(unit= 24)
            close(unit= 27)
            close(unit= 28)
         endif
      
      END IF
 
!c     ... open the surface input files
      if (me.eq.0)then
         open(unit=110,file='input/zo.out',status='unknown')
         do j=1,ny
            do i=1,nx
               read(110,*) zo(i,j)
!            zo(i,j)=0.000001d0
            enddo
         enddo
      close(unit=110)
         if(S_flag.eq.1)then
            if(surf_flag.eq.1)then
               open(unit=781,file='input/coolrate.out',status='unknown')
               open(unit=111,file='input/t_surf.out',status='unknown')
               do j=1,ny
                  do i=1,nx
                     read(111,*) t_s(i,j)
                     read(781,*) coolrate(i,j)
!                     t_s(i,j)=293.0
!                     coolrate(i,j)=c_coolrate
                     t_s(i,j)=t_s(i,j)/T_scale
                  enddo
               enddo
               close(111)
               close(781)
            endif
         endif

         if(Q_flag.eq.1)then
            if(qsurf_flag.eq.1)then
               open(unit=783,file='input/qcoolrate.out',
     +              status='unknown')
               open(unit=782,file='input/q_surf.out',status='unknown')
               do j=1,ny
                  do i=1,nx
                     read(782,*) q_s(i,j)
                     read(783,*) qcoolrate(i,j)
                     q_s(i,j)=q_s(i,j)/Q_scale
                  enddo
               enddo
               close(782)
               close(783)
            endif
         endif

         call openfiles()
         call openfiles1()

         Print *,'Number of timesteps', Nsteps
         Print *,'DT = ',dt
         print *,'Nx, Ny, Nz= ',nx,ny,nz
         print *,'sampling stats every ',c_count,' timesteps'
         print *,'writing stats every ',p_count,' timesteps'
         print *,'writing u,v,w fields at end of run'

        ! write(fname,'(A)')'output/debug.dat'
        ! open(375,file=fname)
        ! do i=1,ny
        ! do j=1,nx
        !  write(375,*) ML_in(:,1), txz(i,j,2),tyz(i,j,2)
        ! end do
        ! end do
        ! close(375)

      endif


!     now all procs have to run through this initialization for statistics
      call zeroslice(au,av,aw,ap,u2,v2,w2,w3,atxx,atxz,atyy,
     +     atyz,atzz,atxy,p2,auw,avw,auv,adudz,adudx,advdz,adwdz,adwdx,
     +     e,aCs2,aCs,abeta1,atxz_s,aESGS,
     +     aFSu1,aFSu2,aFSu3,aFSu4,aFSv1,aFSv2,aFSv3,aFSv4)
      

      if(S_Flag.eq.1)then
         call s_zeroslice(at,t2,t3,asgs_t1,asgs_t2,asgs_t3,aut,avt,awt,
     +        adtdx,adtdy,adtdz,aPr,aCs2Pr,abeta2,aqz_s,aET,ts_avg)
      endif
      if(Q_flag.eq.1)then
         call s_zeroslice(aq,q2,q3,asgs_q1,asgs_q2,asgs_q3,auq,avq,awq,
     +        adqdx,adqdy,adqdz,aSc,aCs2Sc,abeta3,aqsurf,aEQ,ts_avg)
      endif

      IF(nprocs.gt.1)THEN

        IF (me<=nprocs-2) then
            
            call MPI_RECV(u(1,1,2),(nzb)*nx*ny, MPI_DOUBLE_PRECISION,
     +           nprocs-1,me,nall,status2, ierr )
            call MPI_RECV(v(1,1,2),(nzb)*nx*ny, MPI_DOUBLE_PRECISION,
     +           nprocs-1,me,nall,status2, ierr )
            call MPI_RECV(w(1,1,2),(nzb)*nx*ny, MPI_DOUBLE_PRECISION,
     +           nprocs-1,me,nall,status2, ierr )

            if(initu.eq.1)then
               call MPI_RECV(a1_old(1,1,2),(nzb)*nx*ny, 
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(b1_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(c1_old(1,1,2),(nzb)*nx*ny, 
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(d1_old(1,1,2),(nzb)*nx*ny, 
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(e1_old(1,1,2),(nzb)*nx*ny, 
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               
               call MPI_RECV(a2_old(1,1,2),(nzb)*nx*ny, 
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(b2_old(1,1,2),(nzb)*nx*ny, 
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(c2_old(1,1,2),(nzb)*nx*ny, 
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(d2_old(1,1,2),(nzb)*nx*ny, 
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(e2_old(1,1,2),(nzb)*nx*ny, 
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)

               call MPI_RECV(RHSx(1,1,2),(nzb)*nx*ny, 
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(RHSy(1,1,2),(nzb)*nx*ny, 
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(RHSz(1,1,2),(nzb)*nx*ny, 
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)

            endif

            IF (S_Flag.eq.1) then
            
               call MPI_RECV(theta(1,1,2),(nzb)*nx*ny, 
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               if(inits.eq.1)then
               call MPI_RECV(a4_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(b4_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(c4_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(d4_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(e4_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               
               call MPI_RECV(a8_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(b8_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(c8_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(d8_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(e8_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)

               call MPI_RECV(RHS_T(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)

               endif
               
            ENDIF

            IF (Q_Flag.eq.1) then
            
               call MPI_RECV(q(1,1,2),(nzb)*nx*ny, 
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               if(inits.eq.1)then
               call MPI_RECV(a5_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(b5_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(c5_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(d5_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(e5_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               
               call MPI_RECV(a9_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(b9_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(c9_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(d9_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)
               call MPI_RECV(e9_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)

               call MPI_RECV(RHS_Q(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,nprocs-1,me,nall,status2,ierr)

               endif
               
            ENDIF
            
         ENDIF
         
      ENDIF


      call MPI_Bcast(nrsub,1,MPI_INTEGER,nprocs-1,nall,ierr)         
      if(me.eq.nprocs-1) write(*,*) 'Previous timesteps = ',nrsub

      ct  = mod(nrsub,c_count)+1
      pt  = mod(nrsub,p_count)+1
      fst = mod(nrsub,framestep)+1
 
!     When Buffer is 1
!      if (buffer.eq.1) then 
!       do k=2,NZB+1
!        umean=sum(u(:,:,k))/(Nx*Ny)
!        u(:,:,k)=umean
!        vmean=sum(v(:,:,k))/(Nx*Ny)
!        v(:,:,k)=vmean
!        wmean=sum(w(:,:,k))/(Nx*Ny)
!        w(:,:,k)=wmean			
!       end do	
!      endif	  
!      if (me==0) print*, 'yesct',ct 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(turbine_model.gt.0)then
	  
      do i=1,num_turbine 
      wtz(i)=6
      Zhub(i)=70.0d0/z_i
      wtR(i)=40.0d0/z_i
      wtomega(i)=15.0d0
      CTN(i)=0.85d0
      CTB(i)=1.2d0
      CTT(i)=1.3413d0
      dir(i)=1.0 	  
      end do

      i=1
       wtx(i)=25
       wty(i)=81

!      iT=6
!      jT=3
!      do i=1,iT
!         do j=1,jT
!             kT=(i-1)*jT+j
!             wtx(kT)=41+(i-1)*56
!             wty(kT)=(j-0.5)*(1.0*ny/jT)+1
!         end do
!      end do

      if (me==0) print*, 'yesY-HAWT'
      if (me==0) print*, wty	  
      if (me==0) print*, 'yesX-HAWT'
      if (me==0) print*, wtx
!      if (me==0) print*, 'yesZ', Zhub	  
      
      endif
	  
      if(turbine_model_vawt.gt.0)then
	  
      do i=1,num_turbine_vawt 
      H1bar_vawt(i)=44.40d0/z_i
      H2bar_vawt(i)=68.40d0/z_i	  
      Dbar_vawt(i)=0.01d0/z_i	  
      Zhub_vawt(i)=40.8d0/z_i
      wtH_vawt(i)=24.0d0/z_i	  
      wtR_vawt(i)=13.0d0/z_i
!      wtH_vawt(i)=8.0d0/z_i	  
!      wtR_vawt(i)=5.0d0/z_i
      wtomega_vawt(i)=19.65d0
      CTB_vawt(i)=1.0d0
!      CTT_vawt(i)=0.8315d0
      CTT_vawt(i)=0.64d0	  
      dir_vawt(i)=1.0 
      chord_vawt(i)=0.75/z_i
      TSR_vawt(i)=3.8d0
      end do

!      iT=5
!      jT=6
!      do i=1,iT
!         do j=1,jT
!             kT=(i-1)*jT+j
!             wtx_vawt(kT)=69+(i-1)*56
!             wty_vawt(kT)=(j-0.5)*(1.0*ny/jT)+1
!             dir_vawt(kT)=1.0			 
!         end do
!      end do
	  
      i=1
      wtx_vawt(i)=16
      wty_vawt(i)=37
      dir_vawt(i)=1.0 	 

!      i=2
!      wtx_vawt(i)=16
!      wty_vawt(i)=24
!      dir_vawt(i)=-1.0 	  
	  
      if (me==0) print*, 'yesY-VAWT'
      if (me==0) print*,  wty_vawt	  
      if (me==0) print*, 'yesX-VAWT'
      if (me==0) print*,  wtx_vawt
      if (me==0) print*,  wtR_vawt*z_i,wtH_vawt*z_i	      
      endif	  
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF(me.eq.0) t01=MPI_WTIME()
      IF(me.eq.0) time1=MPI_WTIME()	  
!     ==============================S
!     ... Start Time Loop      ==========S
!     ... Start Time Loop      ===============S
!     ... Start Time Loop      ==========S
!     ==============================S
      ct=1 

      if (flag_RL == 1) then

        Do i = 1,nx
          innerx:do i1 = 1,nx_agents
            if (i1 .lt. nx_agents) then
              ii = 1 + nx/nx_agents*(i1-1)
              iii = 1 + nx/nx_agents*i1
            else if (i1.eq.nx_agents) then
              ii = 1 + nx/nx_agents*(i1-1)
              iii = nx
            end if
            if (i.eq.nx) then
              ii = nx
              iii = nx+1
            end if
            if (ii.le.i.and.i.lt.iii) then
              neighbor(i,:,1) = ii
              neighbor(i,:,2) = iii
              exit innerx
            endif
          end do innerx
        enddo

        Do i = 1,ny
          innery:do i1 = 1,ny_agents
            if (i1 .lt. ny_agents) then
              ii = 1 + ny/ny_agents*(i1-1)
              iii = 1 + ny/ny_agents*i1
            else if (i1.eq.ny_agents) then
              ii = 1 + ny/ny_agents*(i1-1)
              iii = ny
            end if
            if (i.eq.ny) then
              ii = ny
              iii = ny+1
            end if
            if (ii.le.i.and.i.lt.iii) then
              neighbor(:,i,3) = ii
              neighbor(:,i,4) = iii
              exit innery
            endif
          end do innery
        enddo

        if (RL_init == 0) then
          call smarties_setNumAgents(smarties_comm, nx_agents*ny_agents)
          print *, 'set nagents',smarties_comm,',',nx_agents*ny_agents
          id = 0
          call smarties_setStateActionDims(smarties_comm, state_size,
     +              num_actions, id)
          print *, 'set state actions dim',state_size
          bounded = .true.
          upper_action_bound = (/1.1/)
          lower_action_bound = (/0.9/)
          call smarties_setActionScales(smarties_comm,
     +        c_loc(upper_action_bound), c_loc(lower_action_bound),
     +        bounded, num_actions, id)
          print *,'set actions scales',',',num_actions
          do i = 1,2
            b_observable(i) = 1
          end do
          do i=3,6
            b_observable(i) = 0
          end do
          call smarties_setStateObservable(smarties_comm,
     +                    c_loc(b_observable),state_size, id)
        END IF
      END IF


      DO 728 T=1,Nsteps
      call post_dimen(me)
    

         ttt=t+nrsub

!...  Truncate wavenumbers and calculate ddxy          
         call FILT_DA(u,dudx,dudy,t,1,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,fgr,l_r)
         call FILT_DA(v,dvdx,dvdy,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,fgr,l_r)
         
         if(Q_flag.eq.1.and.t.eq.1)then
            call FILT_DA(q,dqdx,dqdy,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,fgr,l_r)
         endif

         if(S_FLAG.eq.1)then
            call FILT_DA(w,dwdx,dwdy,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,fgr,l_r)
            if(t.eq.1)then
               call FILT_DA(theta,dtdx,dtdy,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,fgr,l_r)
            endif
         else
            call FILT_DA(w,dwdx,dwdy,t,2,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,fgr,l_r)
         endif

!     mpi update ghostlayers
         call update3(u,v,w,me,nall)


         if(t.eq.1)then
            if(S_flag.eq.1)then
               call update1(theta,me,nall)
            endif
            if(Q_flag.eq.1)then
               call update1(q,me,nall)
            endif
         endif
         
!     ... Save previous time's right-hand-sides for Adams-Bashforth 
!     ... Integration (In subroutine "STEP" use first order time
!     ... advancement on first time step)
         
!         IF (t.gt.1) then    
            RHSx_f=RHSx
            RHSy_f=RHSy
            RHSz_f=RHSz
!         ENDIF

!     ... Compute spatial derivatives from u,v,w (phys. space) and 
!     ... place in phys space.  
         call DDZ_UV_p(dudz,u,me)
         call DDZ_UV_p(dvdz,v,me)
         call DDZ_W(dwdz,w,me)   
         
         IF (S_FLAG.eq.1) then
            call ddz_uv_p(dtdz,theta,me)
         ENDIF
         IF (Q_FLAG.eq.1) then
            call ddz_uv_p(dqdz,q,me)
         ENDIF
         
!     ... MPI UPDATE THE GHOSTLAYERS for the derivatives 
!     ... (this is important for the SGS model)
         
         call update9(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,
     +        me,nall)

         if(S_flag.eq.1)then 
            call update3(dtdx,dtdy,dtdz,me,nall)
         endif
         if(Q_flag.eq.1)then
            call update3(dqdx,dqdy,dqdz,me,nall)
         endif


!     Compute Surface Boundary Conditions*************************
!     NOTE: currently q is treated as a passive scalar in the BCs
         IF (me==0) then

!          call surf_flux2(theta,q,u,v,uh,vh,t_flux,q_flux,Psi,Psi0,			   
!     +         fi,fi_H,zo,t_s,q_s,coolrate,ilow,M,t)		 
!            call surf_flux2(theta,u,v,uh,vh,t_flux,Psi,Psi0,			   
!     +           fi,fi_H,zo,t_s,coolrate,ilow,M,t)

            call surf_flux(theta,q,u,v,uh,vh,t_flux,q_flux,Psi,Psi0,
     +              fi,fi_H,zo,t_s,q_s,coolrate,qcoolrate,ilow,M,t)
         end if

         if (t==1.and.flag_RL==1) then
           rewards_old = 0d0
         end if
         IF ((ct.eq.c_count).or.(T==1).or.(T==Nsteps)) then
         if (flag_RL == 1) then
            call send_recv_state_action(u,v,smarties_comm,T,
     +       num_actions,state_size,nx,ny,nz2,dx,dy,dz,vonk,nu,
     +       nsteps,nx_agents,ny_agents,tauw,rewards_old,zo,M,z_i,me,
     +       neighbor,ct,c_count,RL_init,nall,ierr,MPI_DOUBLE_PRECISION,
     +RL_stop,nprocs,status2,mpi_status_size,train,actions_2D,states_2D)
            call MPI_Barrier(  nall, ierr)
         endif
         endif

         IF (me ==0) then
            call wallstress2(txz,tyz,uh,vh,Psi,Psi0,zo,ustar,M,t,
     +          tauw, flag_RL)
            call derivwall2 (dtdz,dqdz,dudz,dvdz,uh,vh,fi,fi_H,
     +                          t_flux,q_flux,ustar,M,t)
         ENDIF
         IF (ct.eq.c_count) then 
		 
!              Zhub(1)=zref
              l_kk = int(zref/dz-0.5d0)
!              l_kk = nint(dz/2./dz-0.5d0)
              if (l_kk .le. 0.) l_kk=0.0d0  
              me_k= l_kk/Nzb
              if (me.eq.me_k) then
			  
!              print*,me,l_kk
              kk_m=l_kk- me*Nzb +2
              xx=(zref/dz-0.5d0)-l_kk			  
              u_hub=((1-xx)*sum(u(:,:,kk_m))+xx*sum(u(:,:,kk_m+1)))/
     +              (Nx*Ny)
              v_hub=((1-xx)*sum(v(:,:,kk_m))+xx*sum(v(:,:,kk_m+1)))/
     +              (Nx*Ny)
              alpha=atan(v_hub/u_hub) 
!              print*,kk_m,u_hub,v_hub,alpha*1e5/pi,xx			  
              do i=0,nprocs-1
              if (i.ne.me)	then		  
              call MPI_SEND(alpha,1, MPI_DOUBLE_PRECISION,
     +              i,i,nall,ierr ) 
              endif 	 
              enddo
!              print*,'send_alpha'
			  
              do i=0,nprocs-1
              if (i.ne.me)	then			  
              call MPI_SEND(u_hub,1, MPI_DOUBLE_PRECISION,
     +              i,i,nall,ierr ) 
              endif		 
              enddo	
!              print*,'send_u_huv'	
			  
              do i=0,nprocs-1
              if (i.ne.me)	then			  
              call MPI_SEND(v_hub,1, MPI_DOUBLE_PRECISION,
     +              i,i,nall,ierr )
              endif		 
              enddo	
!              print*,'send_v_hub'
			   
              endif

              if (me.ne.me_k) then			  
              call MPI_RECV(alpha,1, MPI_DOUBLE_PRECISION,
     +                     me_k,me,nall,status2,ierr )
              call MPI_RECV(u_hub,1, MPI_DOUBLE_PRECISION,
     +                     me_k,me,nall,status2,ierr )	 
              call MPI_RECV(v_hub,1, MPI_DOUBLE_PRECISION,
     +                     me_k,me,nall,status2,ierr ) 
              endif	 
!              print*,'recieve'
            call BL_height(u,v,w,txz,tyz,t_total,t_t,me)
            IF (me>0) then
               call MPI_SEND(t_total(2),nzb,
     +               MPI_DOUBLE_PRECISION,0,me,nall, ierr )
               call MPI_SEND(t_t(2),nzb,
     +               MPI_DOUBLE_PRECISION,0,me,nall, ierr )
            ELSE
               do k=1,nzb
                  tau_total(k)=t_total(k+1)
                  tau_t(k)=t_t(k+1)				  
               enddo  
               do i=1,nprocs-1
                 call MPI_RECV(tau_total(i*nzb+1),nzb,
     +              MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
                 call MPI_RECV(tau_t(i*nzb+1),nzb,
     +              MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )	 
              end do
!              print*, tau_total(1),tau_t(1),tau_total(nz/2),
!     +                tau_t(nz/2)

              do k_bl=2,nz
                  if (tau_t(k_bl)<0.05*tau_t(1)) then
				  h_bl=(k_bl-1)*dz*z_i
!                  do i=1,nprocs-1
!                  call MPI_SEND(h_bl,1, MPI_DOUBLE_PRECISION,
!     +              i,i,nall,ierr ) 
!                  enddo							  
!                  write(*,*),k_bl,tau_total(1),tau_total(k_bl),
!     + 				  tau_total(k_bl+1),h_bl				  
!                  exit
                  go to 101 
				  endif
              enddo
			  
101         END IF 


!            if (me.ne.0) then
!            call MPI_RECV(h_bl,1, MPI_DOUBLE_PRECISION,
!     +                    0,me,nall,status2,ierr )
!            endif  	 

            rmsdivve=0.
            kee=0.
            ke=0.0
            
            if(me==0) then
               kstart=3
            else
               kstart=2
            end if
            
            do k=kstart,Nzb+1
               do j=1,Ny
                  do i=1,Nx
                     ke=ke+u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+
     +                    w(i,j,k)*w(i,j,k)
                  end do
               end do
            end do
          
            ke=ke/(1.d0*(nx*ny*(nz-1)))

            call rmsdiv (rmsdivvel,dudx,dvdy,dwdz,me) 
            call Check_CFL(u,v,w,CFLx,CFLy,CFLz,me,nall)
            
            IF (me>0) THEN
               
               call MPI_SEND(rmsdivvel,1, MPI_DOUBLE_PRECISION,
     +              0,me,nall,ierr ) 
               call MPI_SEND(ke,1, MPI_DOUBLE_PRECISION,
     +              0,me,nall,ierr ) 
               call MPI_SEND(CFLx,1, MPI_DOUBLE_PRECISION,
     +              0,me,nall,ierr )               
               call MPI_SEND(CFLy,1, MPI_DOUBLE_PRECISION,
     +              0,me,nall,ierr )               
               call MPI_SEND(CFLz,1, MPI_DOUBLE_PRECISION,
     +              0,me,nall,ierr )   
!               call MPI_SEND(ncs2g,1, MPI_DOUBLE_PRECISION,
!     +              0,me,nall,ierr )
!               call MPI_SEND(ncs2l,1, MPI_DOUBLE_PRECISION,
!     +              0,me,nall,ierr )	 
            ELSE
               
!     ... Receive and write other data
               if(nprocs.gt.1)then
                  do i=1,nprocs-1
                     
                     call MPI_RECV(rmsdivvel,1, MPI_DOUBLE_PRECISION,
     +                    i,i,nall,status2,ierr ) 
                     call MPI_RECV(ke,1, MPI_DOUBLE_PRECISION,
     +                    i,i,nall,status2,ierr ) 
                     call MPI_RECV(CFLx1,1, MPI_DOUBLE_PRECISION,
     +                    i,i,nall,status2,ierr ) 
                     call MPI_RECV(CFLy1,1, MPI_DOUBLE_PRECISION,
     +                    i,i,nall,status2,ierr ) 
                     call MPI_RECV(CFLz1,1, MPI_DOUBLE_PRECISION,
     +                    i,i,nall,status2,ierr ) 
!                     call MPI_RECV(ncs2g1,1, MPI_DOUBLE_PRECISION,
!     +                    i,i,nall,status2,ierr ) 	
!                     call MPI_RECV(ncs2l1,1, MPI_DOUBLE_PRECISION,
!     +                    i,i,nall,status2,ierr ) 	 

                     rmsdivve=rmsdivve+rmsdivvel
                     kee=kee+ke
!					 ncs2g=ncs2g+ncs2g1
!                     ncs2l=ncs2l+ncs2l1					 

                     CFLx=max(CFLx,CFLx1)
                     CFLy=max(CFLy,CFLy1)
                     CFLz=max(CFLz,CFLz1)                     
                  end do
               endif

! end measure time per step
            time2 = MPI_WTIME()
            write(*,*) 'Elapsed time (s):', time2-time1	   
			
            write (11,5111) ttt,dt,kee,rmsdivve,CFLx,CFLy,CFLz,
     +             ttt*dt*z_i/3600,u_hub,v_hub,alpha*1e5/pi,h_bl,
     +      Ugeo(nz2),Vgeo(nz2),time2-time1	 
	 
            WRITE(*,5555) 't= kee CFL',ttt,dt,kee,rmsdivve, 
     +		CFLx,CFLy,CFLz,ttt*dt*z_i/3600,u_hub,v_hub,alpha*1e5/pi,h_bl,
     +      Ugeo(nz2),Vgeo(nz2) 
	 
! start measure time per step
            time1 = MPI_WTIME()	 
            ENDIF

            call MPI_Bcast(ilow,nx,MPI_INTEGER,0,nall,ierr)
            call avgslice(u,v,w,p,txx,txz,tyy,tyz,tzz,txy,dudz,dudx,
     +           dvdz,dwdz,dwdx,Cs2,beta1,au,av,aw,ap,u2,v2,w2,p2,w3,
     +           atxx,atxz,atyy,atyz,atzz,atxy,auw,avw,auv,adudz,adudx,
     +           advdz,adwdz,adwdx,e,aCs2,aCs,abeta1,atxz_s,ESGS3D,
     +           aESGS,aFSu1,aFSu2,aFSu3,aFSu4,aFSv1,aFSv2,aFSv3,aFSv4,
     +           t,ilow,wgx,me)

 
!          if (me.eq.0) write(*,*) 'avgslice1: Start' 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            call avgslice1(u,v,w,p,txx,txy,txz,tyy,tyz,tzz
     +      ,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
     +      ,dpdx,dpdy,dpdz,Cs2,beta1
     +      ,bau1,bav1,baw1,bap1,bau2,bav2,baw2,bap2,bauv,bauw,bavw
     +      ,batxx,batxy,batxz,batyy,batyz,batzz,badpdx,badpdy,badpdz
     +      ,badudx,badudy,badudz,badvdx,badvdy,badvdz,badwdx,badwdy
     +      ,badwdz,babeta1,baCs,baCs2
     +      ,fx,fy,fz,bafx,bafy,bafz
     +      ,fax,faz,AveWt,bafax,bafaz,bAveWt,baufx,bavfy,bawfz
     +      ,baTt,baTpz,baTsgsz,baTpx,baTsgsx,baTpy,baTsgsy	 
     +      ,baS11,baS22,baS33,baS12,baS13,baS23
     +      ,baeps11,baeps22,baeps33,baeps12,baeps13,baeps23
     +      ,me)
!          if (me.eq.0) write(*,*) 'avgslice1: End' 	 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            If (S_Flag.eq.1) then
	        call update1(Pr2,me,nall)
               call scalar_slice(u,v,w,theta,qx,qy,qz,dtdx,dtdy,dtdz,
     +              Pr2,Cs2,beta2,at,t2,t3,asgs_t1,asgs_t2,asgs_t3,aut,
     +              avt,awt,adtdx,adtdy,adtdz,aPr,aCs2Pr,abeta2,aqz_s,
     +              ET3D,aET,ilow,wgx,t_s,ts_avg,me)
	 
               call scalar_slice1(u,v,w,theta,qx,qy,qz,dtdx,dtdy,dtdz,
     +            Pr2,Cs2,beta2,bat,bat2,baqx,baqy,baqz,baut,bavt,
     +            bawt,badtdx,badtdy,badtdz,baPr2,baBetaH,me)

            ENDIF

            If (Q_Flag.eq.1) then
	        call update1(Sc2,me,nall)
               call scalar_slice(u,v,w,q,sgs_q1,sgs_q2,sgs_q3,dqdx,dqdy,
     +              dqdz,Sc2,Cs2,beta3,aq,q2,q3,asgs_q1,asgs_q2,asgs_q3,
     +              auq,avq,awq,adqdx,adqdy,adqdz,aSc,aCs2Sc,abeta3,
     +              aqsurf,EQ3D,aEQ,ilow,wgx,t_s,ts_avg,me)
	 
               call scalar_slice1(u,v,w,q,sgs_q1,sgs_q2,sgs_q3,dqdx,
     +            dqdy,dqdz,Sc2,Cs2,beta3,baqm,baqm2,baqmx,baqmy,baqmz,
     +            bauqm,bavqm,bawqm,badqmdx,badqmdy,badqmdz,baSc2,
     +            baBetaqm,me)   
            ENDIF
            
            ct=1
            
         ELSE
            
            ct=ct+1
            
         ENDIF	
         
!     ... Compute Convective term in physical space.                  

!          if (me.eq.0) write(*,*) 'dealias1: Start'	 



         call dealias1(u,u_m,t,1,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(v,v_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
         call dealias1(w,w_m,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)

!         DO k = 1,nz2
!           DO i=1,nx2
!             DO j=1,ny2 
!               if (isnan(u_m(i,j,k)).or.isnan(v_m(i,j,k)).or.
!     + isnan(w_m(i,j,k))) print *, 'NAN velocity M before',me,i,j,k,
!     + u_m(i,j,k),
!     + v_m(i,j,k),w_m(i,j,k),u(i,j,k),
!     + v(i,j,k),w(i,j,k)
!             end do
!           enddo
!         enddo


         call update3_m(u_m,v_m,w_m,me,nall)


!          if (me.eq.0) write(*,*) 'dealias1: end'	
!          if (me.eq.0) write(*,*) 'CONVEC: Start'			  

         Call CONVEC(Cx,Cy,Cz,u_m,v_m,w_m,dudy,dudz,dvdx,dvdz,dwdx,dwdy,
     +        t,me,nall,
     +  nx,ny,nx2,ny2,nz2,nzb,nprocs,inxny,nsteps,inx2ny2,
     + plan_f,plan_b,
     + plan_ff,plan_bb)          

!          if (me.eq.0) write(*,*) 'CONVEC: end'	
!          if (me.eq.0) write(*,*) 'SGS: Start'	 
         Call SGS_STAG (t_flux,qx,qy,qz,sgs_q1,sgs_q2,sgs_q3,
     +        u,v,w,dudx,dudy,dudz,dvdx,dvdy,
     +        dvdz,dwdx,dwdy,dwdz,dtdx,dtdy,dtdz,dqdx,dqdy,dqdz,
     +        t,Cs2,Pr2,Sc2,cst,txx,txy,
     +        txz,tyy,tyz,tzz,theta,q,
     +        a1_old,b1_old,c1_old,d1_old,e1_old,
     +        a2_old,b2_old,c2_old,d2_old,e2_old,
     +        a4_old,b4_old,c4_old,d4_old,e4_old,
     +        a8_old,b8_old,c8_old,d8_old,e8_old,
     +        a5_old,b5_old,c5_old,d5_old,e5_old,
     +        a9_old,b9_old,c9_old,d9_old,e9_old,
     +        beta1,beta2,beta3,ESGS3D,ET3D,EQ3D,ncs2g,ncs2l,
     +        nx,ny,nz,nx2,ny2,
     +  nz2,nzb,nsteps,inxny,averaging,Co,cs_count,delta,dx,dy,dz,
     + fgr,g_hat,me,model,mom_nodes,nall,nnn,nprocs,pass_flag,
     +q_flag,ri_flag,s_flag,sc,scl_nodes,t_scale,theta_0,
     +vonk,inx2ny2,nu,plan_f,plan_b,plan_ff,plan_bb)
 
!          if (me.eq.0) write(*,*) 'SGS: End'	 
!CCCCCCCCCCCC OUTPUT TIME STATS WHILE EVERYTHING IS SYNCED UP CCCCCCCCC
!     ... Compute Scalar Stats and advance scalar fields in time

         IF (S_Flag.eq.1) then

            call update1(qz,me,nall)
          
            RHS_Tf=RHS_T
            
            call Scalar_RHS(theta,u_m,v_m,w_m,qx,qy,qz,txz,tyz,dtdx,
     +           dtdy,dtdz,RHS_T,me,t_flux,t,nall)
	 
!            IF (t.eq.1.and.inits.eq.0) RHS_Tf=RHS_T
            
            if(Q_flag.eq.1)then

               call update1(sgs_q3,me,nall)
               
               RHS_Qf=RHS_Q
            
              call Scalar_RHS(q,u_m,v_m,w_m,sgs_q1,sgs_q2,sgs_q3,txz,
     +              tyz,dqdx,dqdy,dqdz,RHS_Q,me,q_flux,t,nall)
               
!               IF (t.eq.1.and.inits.eq.0) RHS_Qf=RHS_Q

            endif

!      ... Sponge Layer **********************************
            if (sponge.eq.1) then
 
               cfrdmp = 1./(rlx_time*u_star/z_i)
               do k=2,Nzb+1
                  ztemp=(me*Nzb+(k-2.))*dz*z_i
                  if (ztemp.ge.z_d.and.ztemp.le.l_z) then
                     do j=1,Ny
                        do i=1,Nx 
                           rdmp(i,j,k)=cfrdmp*0.5*
     +                          (1.0-COS(pi*(ztemp-z_d)/(l_z-z_d)))
                        end do
                     end do
                  else
                     do j=1,Ny
                        do i=1,Nx
                           rdmp(i,j,k)=0.0
                        end do
                     end do
                  end if
               end do            
               
               if (me.eq.(nprocs-1)) then
                  do j=1,Ny
                     do i=1,Nx
                        rdmp(i,j,Nzb+2)=rdmp(i,j,Nzb+1)
                     end do
                  end do
               end if  
               call update1(rdmp,me,nall)
!     *******************************************************
               do k=2,Nzb+1
                  Tbase(k) = 0.0
                  do j=1,Ny
                     do i=1,Nx
                        Tbase(k) = Tbase(k)+theta(i,j,k)
                     end do
                  end do
                  Tbase(k)=Tbase(k)*inxny
                  do j=1,Ny
                     do i=1,Nx
                        RHS_T(i,j,k)=RHS_T(i,j,k)-0.5*(rdmp(i,j,k)+
     +                       rdmp(i,j,k+1))*(theta(i,j,k)-Tbase(k))
                     end do
                  end do
               end do

               if(Q_flag.eq.1)then

                  do k=2,Nzb+1
                     Qbase(k) = 0.0
                     do j=1,Ny
                        do i=1,Nx
                           Qbase(k) = Qbase(k)+q(i,j,k)
                        end do
                     end do
                     Qbase(k)=Qbase(k)*inxny
                     do j=1,Ny
                        do i=1,Nx
                           RHS_Q(i,j,k)=RHS_Q(i,j,k)-0.5*(rdmp(i,j,k)+
     +                          rdmp(i,j,k+1))*(q(i,j,k)-Qbase(k))
                        end do
                     end do
                  end do
                  
               endif
            end if
!     ******End Sponge Layer for Scalars***************************
			
            if(Q_flag.eq.1)then
			
               IF (t.eq.1.and.inits.eq.0) RHS_Qf=RHS_Q			
			
               call Step_S (q, RHS_Q, RHS_Qf,1,me)
               call FILT_DA(q,dqdx,dqdy,t,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,fgr,l_r)
               call update1(q,me,nall)
            endif

!     ******Add in large scale tendency****************************

            if(bufferT.eq.1)then
			
!           call inflow_ts(theta,t,me,nall)
            call inflow_ts(theta,RHS_T_BF,t,me,nall)		 
		    RHS_T=RHS_T-RHS_T_BF
			
            end if

            IF (t.eq.1.and.inits.eq.0) RHS_Tf=RHS_T
			
            call Step_S (Theta, RHS_T, RHS_Tf,0,me)
            call FILT_DA(theta,dtdx,dtdy,t,2,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,fgr,l_r)
            call update1(theta,me,nall)

!            call calcbeta (theta,beta,me)
            call calcbeta (theta,q,beta,me)
!ccccccccccccccccccccccccccccccccccccccccccccccc
         ENDIF
         
!     ... Compute Divergence of SGS shear stresses              
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         call ddz_uv(ddtzz, tzz,me)
         call update1(ddtzz,me,nall)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         call divstress_UV (divtx, txx, txy, txz,me,t)
         call divstress_UV (divty, txy, tyy, tyz,me,t)
         call divstress_W (divtz, txz, tyz, ddtzz,me,t)

!     ... Calculate the proper Ugeo for forcing
         if(press_cor.eq.1)then
         if (t.eq.1) then
            Ugeo=Ugeo_inf
            Vgeo=Vgeo_inf
            alpha_n=0.0d0			
         endif	 
         elseif(press_cor.eq.2)then
            Ptime = dble(ttt)*dt*z_i/u_star/3600.d0
            if(Ptime.lt.3.d0)then
               Ugeo_g=-6.5d0/u_star+(Ptime+1.d0)/4.d0*
     +              (-5.d0+6.5d0)/u_star
               Vgeo_g=4.5d0/u_star+(Ptime+1.d0)/4.d0*
     +              (4.5d0-4.5d0)/u_star
            elseif(Ptime.lt.6.d0)then
               Ugeo_g=-5.d0/u_star+dmod(Ptime,3.d0)/3.d0*
     +              (-5.d0+5.d0)/u_star
               Vgeo_g=4.5d0/u_star+dmod(Ptime,3.d0)/3.d0*
     +              (4.5d0-4.5d0)/u_star
            else
               Ugeo_g=-5.d0/u_star+dmod(Ptime,6.d0)/6.d0*
     +              (-6.5d0+5.d0)/u_star
               Vgeo_g=4.5d0/u_star+dmod(Ptime,6.d0)/6.d0*
     +              (2.5d0-4.5d0)/u_star
            endif
            
            do k=2,nzb+1
               ztemp=(me*nzb+k-1)*dz*z_i-dz*z_i/2.d0
               Ugeo(k)=Ugeo_g*((Ugeo_inf/Ugeo_g-1.d0)*
     +              ztemp/2000.d0+1.d0)
               Vgeo(k)=Vgeo_g*((Vgeo_inf/Vgeo_g-1.d0)*
     +              ztemp/2000.d0+1.d0)
            enddo
         endif

		 
!     ... Calculate dynamic tendencies for momentum

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!....... Compute the drag force Fx
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!          if (me.eq.0) write(*,*) 'Turbine: Start'	
      if(turbine_model.gt.0 .or. turbine_model_vawt.gt.0)then  
        if (turbine_model.gt.0 .and. t.eq.1) then
          Open(unit=4444,file='input/CL_CD.DAT')
            DO k=1,181
            READ(4444,*)ALFA_CL_CD(1:3,k)
          END DO
!          if (me.eq.0) write(*,*) 'yesBladeHAWT'
          CLOSE(4444)	  
        endif	 

        if (turbine_model_vawt.gt.0 .and. t.eq.1) then
          Open(unit=4444,file='input/NACA0018_1e6_N.DAT')
            DO k=1,361
              READ(4444,*)ALFA_CL_CD_vawt(1:3,k)
            END DO
!          if (me.eq.0) write(*,*) 'yesBladeVAWT'			
          CLOSE(4444)	  
        endif	
		
!           call calturbine(alpha,u_hub,v_hub,u,v,w,Fx,Fy,Fz,
!     +                      Fax,Faz,AveWt,t,me,nall)
	 
         call calturbine(u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)	 

!       if(t.eq.1 .or. mod(t,c_count)==0)then
!          if (me.eq.0) write(*,*) 'WT Start'	   
!         call calturbine(u,v,w,Fx,Fy,Fz,Fax,Faz,AveWt,t,me,nall)
!          if (me.eq.0) write(*,*) 'WT End'		 
!       end if		 
!      print*,   ALFA_CL_CD(1:3,181)
!        if(t.eq.1 .or. mod(t,c_count)==0)then
!           call calturbine(alpha,u_hub,v_hub,u,v,w,Fx,Fy,Fz,
!     +                      Fax,Faz,AveWt,t,me,nall)
!        end if
!        Fx=Fx*(1-exp(-t*dt*z_i/3600.0))
!        Fy=Fy*(1-exp(-t*dt*z_i/3600.0))
!        Fz=Fz*(1-exp(-t*dt*z_i/3600.0))
      end if
		 
!     ... Compute preliminary RHS matrices for pressure calculation         

         do k=2,Nzb+1
            do j=1,Ny
               do i=1,Nx
                  if (press_cor.eq.0) then
                     RHSx(i,j,k)=-Cx(i,j,k)-divtx(i,j,k)+Uadv(k)
     +				 -fx(i,j,k)
                     RHSy(i,j,k)=-Cy(i,j,k)-divty(i,j,k)+Vadv(k)
     +				 -fy(i,j,k)
                  else
                     RHSx(i,j,k)=-Cx(i,j,k)-divtx(i,j,k)
     +                    -fc*(Vgeo(k)-Vgal-v(i,j,k))+Uadv(k)-fx(i,j,k)
!     +                    -1.0*alpha_n/dt*(v(i,j,k)-Vgeo(k))	 
                     RHSy(i,j,k)=-Cy(i,j,k)-divty(i,j,k)
     +                    +fc*(Ugeo(k)-Ugal-u(i,j,k))+Vadv(k)-fy(i,j,k)
!     +                    +1.0*alpha_n/dt*(u(i,j,k)-Ugeo(k)) 	 
                  end if
                  if (pass_flag.eq.0) then
                    RHSz(i,j,k)=-Cz(i,j,k)-divtz(i,j,k)+beta(i,j,k)
     +              -fz(i,j,k)
                  else
                    RHSz(i,j,k)=-Cz(i,j,k)-divtz(i,j,k)-fz(i,j,k)
                  end if
               end do
            end do
         end do
         
!     ********Sponge Layer for Momentum****************************
         if (sponge.eq.1) then
            do k=2,Nzb+1
               do j=1,Ny
                  do i=1,Nx
                     RHSx(i,j,k)=RHSx(i,j,k)-0.5*
     +                    (rdmp(i,j,k)+rdmp(i,j,k+1))*
     +                    (u(i,j,k)-(Ugeo(k)-Ugal))
                     RHSy(i,j,k)=RHSy(i,j,k)-0.5*
     +                    (rdmp(i,j,k)+rdmp(i,j,k+1))*
     +                    (v(i,j,k)-(Vgeo(k)-Vgal))
                  end do
               end do
            end do
            do k=2,Nzb+1
               do j=1,Ny
                  do i=1,Nx
                     RHSz(i,j,k)=RHSz(i,j,k)-rdmp(i,j,k)*w(i,j,k)
                  end do
               end do
            end do
         end if
!     *********End of Sponge Layer for Momentum********************
!Mahdi add to calculate buffer zone
         if (buffer.eq.1 .and. t.eq.1) then

!         if (me.eq.0) print*, 'buffer:start'		 
         Open (unit=200,file='input/vel_inflow.out',status='unknown')
         Open (unit=201,file='input/vel2_inflow.out',status='unknown')		
!         Open (unit=200,file='input/vel.out',status='unknown')
!         Open (unit=201,file='input/vel.out',status='unknown')		 
         IF(ME.EQ.nprocs-1) THEN
         read (200,*) 
         read (200,*) 		 
         read (200,*)	
         IF(nprocs.gt.1)THEN            
            DO II=0,nprocs-2               
               do k=1,NZB
                  do j=1,Ny
                     do i=1,Nxx
                        read (200,*) uinf(i,j,k),vinf(i,j,k),winf(i,j,k)
                     end do
                  end do
               end do               
               call MPI_SEND(uinf(1,1,1),(nzb)*nxx*ny, 
     +              MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
               call MPI_SEND(vinf(1,1,1),(nzb)*nxx*ny, 
     +              MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
               call MPI_SEND(winf(1,1,1),(nzb)*nxx*ny, 
     +              MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )               
            END DO            
         ENDIF       
!     read for himself
         do k=2,NZB+1
            do j=1,Ny
               do i=1,Nxx
                  read (200,*) uinf(i,j,k),vinf(i,j,k),winf(i,j,k)
               end do
            end do
         end do	
		 
         read (201,*) 
         read (201,*) 		 
         read (201,*)	
         IF(nprocs.gt.1)THEN            
            DO II=0,nprocs-2               
               do k=1,NZB
                  do j=1,Ny
                     do i=1,Nxx
          read (201,*) u2inf(i,j,k),v2inf(i,j,k),w2inf(i,j,k)
                     end do
                  end do
               end do               
               call MPI_SEND(u2inf(1,1,1),(nzb)*nxx*ny, 
     +              MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
               call MPI_SEND(v2inf(1,1,1),(nzb)*nxx*ny, 
     +              MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )
               call MPI_SEND(w2inf(1,1,1),(nzb)*nxx*ny, 
     +              MPI_DOUBLE_PRECISION,ii,ii,nall, ierr )               
            END DO            
         ENDIF       
!     read for himself
         do k=2,NZB+1
            do j=1,Ny
               do i=1,Nxx
          read (201,*) u2inf(i,j,k),v2inf(i,j,k),w2inf(i,j,k)
               end do
           end do
         end do			 
         ENDIF
		 
         IF(nprocs.gt.1)THEN
         IF (me<=nprocs-2) then           
          call MPI_RECV(uinf(1,1,2),(nzb)*nxx*ny, MPI_DOUBLE_PRECISION,
     +           nprocs-1,me,nall,status2, ierr )
          call MPI_RECV(vinf(1,1,2),(nzb)*nxx*ny, MPI_DOUBLE_PRECISION,
     +           nprocs-1,me,nall,status2, ierr )
          call MPI_RECV(winf(1,1,2),(nzb)*nxx*ny, MPI_DOUBLE_PRECISION,
     +           nprocs-1,me,nall,status2, ierr )
          call MPI_RECV(u2inf(1,1,2),(nzb)*nxx*ny, MPI_DOUBLE_PRECISION,
     +           nprocs-1,me,nall,status2, ierr )
          call MPI_RECV(v2inf(1,1,2),(nzb)*nxx*ny, MPI_DOUBLE_PRECISION,
     +           nprocs-1,me,nall,status2, ierr )
          call MPI_RECV(w2inf(1,1,2),(nzb)*nxx*ny, MPI_DOUBLE_PRECISION,
     +           nprocs-1,me,nall,status2, ierr )	 
         ENDIF
         ENDIF
         endif 

!         if (me.eq.0) print*, 'buffer:end'		 
!         if(buffer.eq.1)then		 
!          call inflow_t(u,v,w,uinf,vinf,winf,u2inf,v2inf,
!     +                  w2inf,t,me,nall)

!         if (me.eq.0) print*, 'buffer:no'
         if(buffer.eq.1)then	
!         if (me.eq.0) print*, 'buffer:yes'  		 
!         call inflow_t(u,v,w,uinf,vinf,winf,
!     +                 u2inf,v2inf,w2inf,t,me,nall)

!         call inflow_t(u,v,w,RHSx_BF,RHSy_BF,RHSz_BF,t,me,nall)
         call inflow_t(u,v,w,
     +		 uinf,vinf,winf,u2inf,v2inf,w2inf,
     +		 RHSx_BF,RHSy_BF,RHSz_BF,t,me,nall)
		 
         RHSx=RHSx-RHSx_BF
         RHSy=RHSy-RHSy_BF
         RHSz=RHSz-RHSz_BF		 	 
         end if
		 
!         if(bufferT.eq.1)then
!         call inflow_ts(theta,t,me,nall)
!         end if		         

!Mahdi end of buffer zone		 
!     ... Solve Poisson Equation for pressure.  

         call Press_STAG(P,RHSx,RHSy,RHSz,RHSx_f,RHSy_f,RHSz_f,
     +        u,v,w,dpdx,dpdy,divtz,me,nall,t) 

         call update1(P,me,nall)
         
         
         Call DDZ_UV (dpdz,P,me)

!...  Finalize right hand side for explicit time advancement
         
         do k=2,Nzb+1
            do j=1,Ny
               do i=1,Nx
                  RHSx(i,j,k)=  RHSx(i,j,k) - dpdx(i,j,k) 
                  RHSy(i,j,k)=  RHSy(i,j,k) - dpdy(i,j,k) 
                  RHSz(i,j,k)=  RHSz(i,j,k) - dpdz(i,j,k)
               end do
            end do
         end do
                  
!     ... Start first step using new C as former C   
         
         IF (t.eq.1.and.initu.eq.0) then 
            RHSx_f=RHSx
            RHSy_f=RHSy
            RHSz_f=RHSz
         ENDIF
         
!     ... Compute required mean forcing to keep total x momentum constant.
 1       do k=2,Nzb+1
!            if (press_cor.eq.0 .and. me.le.6) then
            if (press_cor.eq.0) then			
               force(k)=u2_surface*z_i/l_z
            else
               force(k)=0.0
            end if
         end do
         
         call stepbl_UV(u,RHSx,RHSx_f,force,me) 

         do k=2,Nzb+1
            force(k)=0.0
         end do

         Call STEPBL_UV (v,RHSy,RHSy_f,force,me)              
         Call STEPBL_W (w, RHSz,RHSz_f,me)              

         alpha_n=0.d0*alpha/1000.d0
!         alpha_n=(0.9*(500.0*alpha_n)+0.1*alpha)/500.d0		 
!         Ugeo_n=Ugeo*dcos(-alpha_n)-Vgeo*dsin(-alpha_n)
!         Vgeo_n=Ugeo*dsin(-alpha_n)+Vgeo*dcos(-alpha_n)	
!         Ugeo=Ugeo_n
!         Vgeo=Vgeo_n			 
!         u_n=u*dcos(-alpha_n)-v*dsin(-alpha_n)
!         v_n=u*dsin(-alpha_n)+v*dcos(-alpha_n)	 
!         u=u_n
!         v=v_n 		 

!        umean_top=0.0	
!         tmean_top=0.0		 
!         if (me==26) then
!         do j=1,Ny
!            do i=1,Nx
!              u(i,j,Nzb+1)=u(i,j,Nzb)
!               v(i,j,Nzb+1)=v(i,j,Nzb)
!               w(i,j,Nzb+1)=0.0	
!               theta(i,j,Nzb+1)=theta(i,j,Nzb)			   
!               umean_top=umean_top+u(i,j,Nzb+1)
!               tmean_top=tmean_top+theta(i,j,Nzb+1)			   
!            end do
!         end do		    
!         umean_top=umean_top*inxny
!         tmean_top=tmean_top*inxny		 
!         do i=0,nprocs-1
!           call MPI_SEND(umean_top,1, MPI_DOUBLE_PRECISION,
!     +              i,i,nall,ierr ) 
!         enddo	
!         do i=0,nprocs-1
!           call MPI_SEND(tmean_top,1, MPI_DOUBLE_PRECISION,
!     +              i,i,nall,ierr ) 
!         enddo		 
!        endif
!         call MPI_RECV(umean_top,1, MPI_DOUBLE_PRECISION,
!     +                     26,me,nall,status2,ierr )	
!         call MPI_RECV(tmean_top,1, MPI_DOUBLE_PRECISION,
!     +                     24,me,nall,status2,ierr )		 
!
!         if (me .ge. 27) then
!         do k=2,Nzb+1
!         do j=1,Ny
!            do i=1,Nx
!               u(i,j,k)=umean_top
!               v(i,j,k)=0.0
!               w(i,j,k)=0.0		
!               theta(i,j,k)=tmean_top+0.0*inversion*((me-25)*2+(k-1))*dz				   
!           end do
!         end do
!         end do		 
!         endif 		
	 
!     ... Print stats to file
         
         IF (pt.eq.p_count) then
            
            do k=2,nzb+1
               
               call spectrum (u,k,specu,t,1)
               call spectrum (v,k,specv,t,0)
               IF (S_Flag.eq.1) then
                  call spectrum (theta,k,spect,t,0)
               ENDIF
               IF (Q_Flag.eq.1) then
                  call spectrum (q,k,specq,t,0)
               ENDIF
               call spectrum (w,k,specw,t,2)
            end do
          if ((me==0).and.(t.ne.Nsteps)) call output_tauw(tauw,t,nx,ny) 
          if ((me==0).and.(t.ne.Nsteps)) call
     +  output_actions(states_2D,actions_2D,t,nx_agents,ny_agents)
            call output_average(specu,nx/2+1,299,me,nall)
            call output_average(specv,nx/2+1,298,me,nall)
            call output_average(specw,nx/2+1,297,me,nall)
            
            IF(S_flag.eq.1)then
              call output_average(spect,nx/2+1,300,me,nall)
            endif
            IF(Q_flag.eq.1)then
               call output_average(specq,nx/2+1,915,me,nall)
            endif
!ccc FORCING PARAMETERS ccc
            call output_average(Ugeo,1,61,me,nall)
            call output_average(Vgeo,1,62,me,nall)
            call output_average(Uadv,1,63,me,nall)
            call output_average(Vadv,1,64,me,nall)
            call output_average(Tadv,1,65,me,nall)
            call output_average(Qadv,1,66,me,nall)
            call output_average(au,anx,71,me,nall)
            call output_average(av,anx,72,me,nall)
            call output_average(aw,anx,73,me,nall)
            call output_average(ap,anx,85,me,nall)
            call output_average(u2,anx,74,me,nall)
            call output_average(v2,anx,75,me,nall)
            call output_average(w2,anx,76,me,nall)
            call output_average(w3,anx,40,me,nall)
            call output_average(p2,anx,81,me,nall)
            call output_average(atxx,anx,57,me,nall)
            call output_average(atxz,anx,78,me,nall)
            call output_average(atyy,anx,79,me,nall)
            call output_average(atyz,anx,80,me,nall)
            call output_average(atzz,anx,87,me,nall)
            call output_average(atxy,anx,90,me,nall)
            call output_average(auw,anx,82,me,nall)
            call output_average(avw,anx,83,me,nall)
            call output_average(auv,anx,77,me,nall)
            call output_average(adudx,anx,91,me,nall)
            call output_average(adudz,anx,88,me,nall)
            call output_average(advdz,anx,95,me,nall)
            call output_average(adwdx,anx,92,me,nall)
            call output_average(adwdz,anx,93,me,nall)
            call output_average(e,anx,39,me,nall)
            call output_average(aCs2,anx,398,me,nall)
            call output_average(aCs,anx,399,me,nall)
            call output_average(abeta1,anx,978,me,nall)
            call output_average(aESGS,anx,888,me,nall)
			
            call output_average(aFSu1,anx,889,me,nall)
            call output_average(aFSu2,anx,890,me,nall)
            call output_average(aFSu3,anx,891,me,nall)
            call output_average(aFSu4,anx,892,me,nall)	

            call output_average(aFSv1,anx,893,me,nall)
            call output_average(aFSv2,anx,894,me,nall)
            call output_average(aFSv3,anx,895,me,nall)
            call output_average(aFSv4,anx,896,me,nall)			
			

            write (132) atxz_s

            if(S_flag.eq.1)then
               call output_average(at,anx,100,me,nall)
               call output_average(t2,anx,102,me,nall)
               call output_average(t3,anx,41,me,nall)
               call output_average(asgs_t1,anx,112,me,nall)
               call output_average(asgs_t2,anx,116,me,nall)
               call output_average(asgs_t3,anx,104,me,nall)
               call output_average(aut,anx,113,me,nall)
               call output_average(avt,anx,117,me,nall)
               call output_average(awt,anx,106,me,nall)
               call output_average(adtdx,anx,114,me,nall)
               call output_average(adtdy,anx,115,me,nall)
               call output_average(adtdz,anx,89,me,nall)
               call output_average(aPr,anx,1526,me,nall)
               call output_average(aCs2Pr,anx,1527,me,nall)
               call output_average(abeta2,anx,979,me,nall)
               call output_average(aET,anx,999,me,nall)

               write (133) aqz_s
               write (135) ts_avg			   
            endif

            if(Q_flag.eq.1)then
               call output_average(aq,anx,906,me,nall)
               call output_average(q2,anx,907,me,nall)
               call output_average(q3,anx,908,me,nall)
               call output_average(asgs_q1,anx,909,me,nall)
               call output_average(asgs_q2,anx,910,me,nall)
               call output_average(asgs_q3,anx,911,me,nall)
               call output_average(auq,anx,913,me,nall)
               call output_average(avq,anx,914,me,nall)
               call output_average(awq,anx,912,me,nall)
               call output_average(adqdx,anx,904,me,nall)
               call output_average(adqdy,anx,905,me,nall)
               call output_average(adqdz,anx,903,me,nall)
               call output_average(aSc,anx,916,me,nall)
               call output_average(aCs2Sc,anx,917,me,nall)
               call output_average(abeta3,anx,918,me,nall)
               call output_average(aEQ,anx,919,me,nall)

               write (920) aqsurf
            endif
         
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

            call zeroslice(au,av,aw,ap,u2,v2,w2,w3,atxx,atxz,atyy,
     +           atyz,atzz,atxy,p2,auw,avw,auv,adudz,adudx,advdz,adwdz,
     +           adwdx,e,aCs2,aCs,abeta1,atxz_s,aESGS,
     +           aFSu1,aFSu2,aFSu3,aFSu4,aFSv1,aFSv2,aFSv3,aFSv4)
            
            if(S_Flag.eq.1)then
               call s_zeroslice(at,t2,t3,asgs_t1,asgs_t2,asgs_t3,aut,
     +              avt,awt,adtdx,adtdy,adtdz,aPr,aCs2Pr,abeta2,
     +              aqz_s,aET,ts_avg)
            endif
            
            if(Q_Flag.eq.1)then

               call s_zeroslice(aq,q2,q3,asgs_q1,asgs_q2,asgs_q3,auq,
     +              avq,awq,adqdx,adqdy,adqdz,aSc,aCs2Sc,abeta3,
     +              aqsurf,aEQ,ts_avg)

            endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            call post(501,bau1,me,nall) 
            call post(502,bav1,me,nall)
            call post(503,baw1,me,nall)
            call post(504,bap1,me,nall)
            
            call post(505,bau2,me,nall)
            call post(506,bav2,me,nall)
            call post(507,baw2,me,nall)   
            call post(508,bap2,me,nall)

            call post(509,bauv,me,nall) 
            call post(510,bauw,me,nall) 
            call post(511,bavw,me,nall)   
            
            call post(512,batxx,me,nall) 
            call post(513,batxy,me,nall) 
            call post(514,batxz,me,nall)
            call post(515,batyy,me,nall)
            call post(516,batyz,me,nall)
            call post(517,batzz,me,nall)
            
            call post(518,badudx,me,nall)
            call post(519,badudy,me,nall)
            call post(520,badudz,me,nall)
            
            call post(521,badvdx,me,nall)
            call post(522,badvdy,me,nall)
            call post(523,badvdz,me,nall)
            
            call post(524,badwdx,me,nall)
            call post(525,badwdy,me,nall)
            call post(526,badwdz,me,nall)

            call post(527,badpdx,me,nall)
            call post(528,badpdy,me,nall)
            call post(529,badpdz,me,nall)

            call post(530,baCs,me,nall)
            call post(531,babeta1,me,nall)

            call post(605,baS11,me,nall)
            call post(606,baS22,me,nall)
            call post(607,baS33,me,nall)
            call post(608,baS12,me,nall)
            call post(609,baS13,me,nall)
            call post(610,baS23,me,nall)
            call post(611,baeps11,me,nall)
            call post(612,baeps22,me,nall)
            call post(613,baeps33,me,nall)
            call post(614,baeps12,me,nall)
            call post(615,baeps13,me,nall)
            call post(616,baeps23,me,nall)			
            call post(617,baTt,me,nall)			
            call post(618,baTpz,me,nall)			
            call post(619,baTsgsz,me,nall)
            call post(621,baTpx,me,nall)			
            call post(623,baTsgsx,me,nall)
            call post(622,baTpy,me,nall)			
            call post(624,baTsgsy,me,nall)			

            call post_instant3(551,u,me,nall)
            call post_instant3(552,v,me,nall)
            call post_instant3(553,w,me,nall)
            call post_instant3(554,dudx,me,nall)
            call post_instant3(555,dudy,me,nall)
            call post_instant3(556,dudz,me,nall)
            call post_instant3(557,dvdx,me,nall)
            call post_instant3(558,dvdy,me,nall)
            call post_instant3(559,dvdz,me,nall)
            call post_instant3(560,dwdx,me,nall)
            call post_instant3(561,dwdy,me,nall)
            call post_instant3(562,dwdz,me,nall)
            call post_instant3(563,dpdx,me,nall)
            call post_instant3(564,dpdy,me,nall)
            call post_instant3(565,dpdz,me,nall)
!            call post_instant3(625,P,me,nall)	
!            call post_instant3(625,u,me,nall)			
			
            if(S_Flag.eq.1)then
!            call post_instant3(620,theta,me,nall)			
            call post(532,bat,me,nall)
            call post(533,bat2,me,nall)
            call post(534,baqx,me,nall)
            call post(535,baqy,me,nall)
            call post(536,baqz,me,nall)
            call post(537,baut,me,nall)
            call post(538,bavt,me,nall)
            call post(539,bawt,me,nall)
            call post(540,badtdx,me,nall)
            call post(541,badtdy,me,nall)
            call post(542,badtdz,me,nall)
            call post(543,baPr2,me,nall)
            call post(544,baBetaH,me,nall)
            end if

            if(Q_Flag.eq.1)then
            call post(589,baqm,me,nall)
            call post(590,baqm2,me,nall)
            call post(591,baqmx,me,nall)
            call post(592,baqmy,me,nall)
            call post(593,baqmz,me,nall)
            call post(594,bauqm,me,nall)
            call post(595,bavqm,me,nall)
            call post(596,bawqm,me,nall)
            call post(597,badqmdx,me,nall)
            call post(598,badqmdy,me,nall)
            call post(599,badqmdz,me,nall)
            call post(600,baSc2,me,nall)
            call post(601,baBetaqm,me,nall)
            end if

			
            if(turbine_model.gt.0 .or. turbine_model_vawt.gt.0)then
            call post(545,bafx,me,nall) 
            call post(546,bafy,me,nall)  
            call post(547,bafz,me,nall)   
            call post(602,baufx,me,nall) 
            call post(603,bavfy,me,nall) 
            call post(604,bawfz,me,nall) 			
            call post(548,bafax,me,nall)
            call post(549,bafaz,me,nall)
            call post(550,bAveWT,me,nall)
            end if
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            
            pt=1
            
         ELSE
            
            pt=pt+1
            
         ENDIF
         
!C     =========== C
          if(vel_instant.eq.1)then
          call post_instant2_parallel(1,u,t,me,nall)
          call post_instant2_parallel(2,v,t,me,nall)
          call post_instant2_parallel(3,w,t,me,nall)
          endif
          if(T_instant.eq.1)then
          call post_instant2_parallel(4,theta,t,me,nall)
          endif
!C     =========== C
         if(ttt.ge.sframe .and. nframe.gt.0 .and.
     +     ttt.le.sframe+nframe*framestep)then

            if((ttt.eq.sframe).or.(fst.eq.framestep))then	

!            call post_instant3(551,u,me,nall)
!            call post_instant3(552,v,me,nall)
!            call post_instant3(553,w,me,nall)
            call post_instant3(625,u,me,nall)			
			
            if(S_Flag.eq.1)then
            call post_instant3(620,theta,me,nall)			
            end if


			
!             IF (me>0) then
!                   call MPI_SEND(u(1,1,2),nzb*nx*ny,
!     +                    MPI_DOUBLE_PRECISION,0,me,nall,ierr)
!                   call MPI_SEND(v(1,1,2),nzb*nx*ny,
!     +                    MPI_DOUBLE_PRECISION,0,me,nall,ierr)
!                   call MPI_SEND(w(1,1,2),nzb*nx*ny,
!     +                    MPI_DOUBLE_PRECISION,0,me,nall,ierr)
!
!                   if(S_flag.eq.1)then
!                      call MPI_SEND(theta(1,1,2),nzb*nx*ny,
!     +                    MPI_DOUBLE_PRECISION,0,me,nall,ierr)
!                   endif
!                   if(Q_flag.eq.1)then
!                      call MPI_SEND(q(1,1,2),nzb*nx*ny,
!     +                    MPI_DOUBLE_PRECISION,0,me,nall,ierr)
!                   endif
!             ELSE
!             do k=1,nzb
!                    u_frame(:,:,k)=u(:,:,k+1)
!                    v_frame(:,:,k)=v(:,:,k+1)
!                    w_frame(:,:,k)=w(:,:,k+1)	
!                   if(S_flag.eq.1)then
!                    t_frame(:,:,k)=theta(:,:,k+1)  
!                   endif
!                   if(Q_flag.eq.1)then
!                    q_frame(:,:,k)=q(:,:,k+1)
!                   endif				   
!             enddo
!             write(994) u_frame(:,:,:)+Ugal
!             write(995) v_frame(:,:,:)+Vgal
!             write(998) w_frame(:,:,:)
!                   
!             if(s_flag.eq.1)then
!                 write(997) t_frame(:,:,:)
!             endif
!             if(s_flag.eq.1)then
!                 write(901) q_frame(:,:,:)
!             endif
! 				  
!             do ii=1,nprocs-1
!                 call MPI_RECV(u_frame(1,1,0*nzb+1),nzb*nx*ny, 
!     +               MPI_DOUBLE_PRECISION,ii,ii,nall,
!     +               status2,ierr)
!                 call MPI_RECV(v_frame(1,1,0*nzb+1),nzb*nx*ny, 
!     +               MPI_DOUBLE_PRECISION,ii,ii,nall,
!     +               status2,ierr)
!                 call MPI_RECV(w_frame(1,1,0*nzb+1),nzb*nx*ny, 
!     +               MPI_DOUBLE_PRECISION,ii,ii,nall,
!     +               status2,ierr)
!                         
!                 if(s_flag.eq.1)then
!                 call MPI_RECV(t_frame(1,1,0*nzb+1),
!     +               nzb*nx*ny,MPI_DOUBLE_PRECISION,
!     +               ii,ii,nall,status2,ierr)
!                 endif
!                 if(q_flag.eq.1)then
!                 call MPI_RECV(q_frame(1,1,0*nzb+1),
!     +               nzb*nx*ny,MPI_DOUBLE_PRECISION,
!     +               ii,ii,nall,status2,ierr)
!                 endif			
!             write(994) u_frame(:,:,:)+Ugal
!             write(995) v_frame(:,:,:)+Vgal
!             write(998) w_frame(:,:,:)
!                   
!             if(s_flag.eq.1)then
!                 write(997) t_frame(:,:,:)
!             endif
!             if(s_flag.eq.1)then
!                 write(901) q_frame(:,:,:)
!             endif
!				 
!             end do
!
!             END IF			
 
!cccc output the instantaneous stress and flux at the wall cccc
!                  write(132) ustar
!                  if(s_flag.eq.1)then
!                     write(133) t_flux
!                  endif                  
!                  if(q_flag.eq.1)then
!                     write(920) q_flux
!                  endif

               fst=1		   
               
            else
               
               fst=fst+1
               
            endif
            
         endif
		 
 
         rule=t/ruler
         
         IF (rule*ruler.eq.t) then

            if(me.eq.0)then
               open(unit=110,file='input/zo.out',
     +              status='unknown')
               do j=1,ny
                  do i=1,nx
                     write(110,*) zo(i,j)
                  enddo
               enddo
               close(110)
            endif
            if(S_flag.eq.1)then
               if(me.eq.0)then
                  if(surf_flag.eq.1)then
                     open(unit=111,file='input/t_surf.out',
     +                    status='unknown')
                     open(unit=781,file='input/coolrate.out',
     +                    status='unknown')
                     do j=1,ny
                        do i=1,nx
                           write(781,*) coolrate(i,j)
                           write(111,*) t_s(i,j)*T_scale
                        enddo
                     enddo
                     close(111)
                     close(781)
                  elseif(surf_flag.eq.2)then
                     open(unit=111,file='input/t_025meters.out',
     +                    status='unknown')
                     do i=1,nhrs
                        write(111,*) t_s2(i)*T_scale
                     enddo
                     close(111)
                  endif
               endif
            endif
            if(Q_flag.eq.1)then
               if(me.eq.0)then
                  if(qsurf_flag.eq.1)then
                     open(unit=782,file='input/q_surf.out',
     +                    status='unknown')
                     open(unit=783,file='input/qcoolrate.out',
     +                    status='unknown')
                     do j=1,ny
                        do i=1,nx
                           write(783,*) qcoolrate(i,j)
                           write(782,*) q_s(i,j)*Q_scale
                        enddo
                     enddo
                     close(782)
                     close(783)
                  elseif(qsurf_flag.eq.2)then
                     open(unit=782,file='input/q_025meters.out',
     +                    status='unknown')
                     do i=1,nhrs
                        write(782,*) q_s2(i)*Q_scale
                     enddo
                     close(782)
                  endif
               endif
            endif
            
            IF (me>0) then
               
               call MPI_SEND(u(1,1,2),(nzb)*nx*ny,MPI_DOUBLE_PRECISION,
     +              0,me,nall, ierr )
               call MPI_SEND(v(1,1,2),(nzb)*nx*ny,MPI_DOUBLE_PRECISION,
     +              0,me,nall, ierr )
               call MPI_SEND(w(1,1,2),(nzb)*nx*ny,MPI_DOUBLE_PRECISION,
     +              0,me,nall, ierr )
               
               call MPI_SEND(a1_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,0,me,nall, ierr )
               call MPI_SEND(b1_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,0,me,nall, ierr )
               call MPI_SEND(c1_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,0,me,nall, ierr )
               call MPI_SEND(d1_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,0,me,nall, ierr )
               call MPI_SEND(e1_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,0,me,nall, ierr )

               if(S_Flag.eq.1)then
                  call MPI_SEND(theta(1,1,2),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
               endif
               if(Q_Flag.eq.1)then
                  call MPI_SEND(q(1,1,2),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
               endif

            else
               
               Open (unit=84,file='input/vel.out',status='unknown')			   
               Open (unit=20,file='input/a1_e1.bin',status='unknown',
     +              form='unformatted')
               IF (S_Flag.eq.1) then
                  Open (unit=108,file='input/t.out',status='unknown')
               ENDIF
               IF (Q_Flag.eq.1) then
                  Open (unit=109,file='input/q.out',status='unknown')
               ENDIF
               
               write (84,*) ttt, dt
               write (84,*) Nx,Ny,Nz
               write (84,*)

               do k=2,Nzb+1
                  do j=1,Ny
                     do i=1,Nx
                        
                        write (84,6623) u(i,j,k)+Ugal,v(i,j,k)+Vgal,
     +                       w(i,j,k)
	 
                        write(20) a1_old(i,j,k),b1_old(i,j,k),
     +                            c1_old(i,j,k),d1_old(i,j,k),
     +                            e1_old(i,j,k)
			   
                        
                        IF (S_Flag.eq.1) then
                           write (108,*) theta(i,j,k)
                        ENDIF
                        IF (Q_Flag.eq.1) then
                           write (109,*) q(i,j,k)
                        ENDIF
                     end do
                  end do
               end do

               dudz=u
               dvdz=v
               dwdz=w
			   
               dudx=a1_old
               dvdx=b1_old
               dwdx=c1_old			   
               dudy=d1_old
               dvdy=e1_old			   
               
               IF (S_Flag.eq.1) then
                  dtdz=theta
               ENDIF
               IF (Q_Flag.eq.1) then
                  dqdz=q
               ENDIF
                            
               IF(nprocs.gt.1)THEN
                  
                  DO ii=1,nprocs-1
                     
                     call MPI_RECV(u(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(v(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(w(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     
                     call MPI_RECV(a1_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(b1_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(c1_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(d1_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(e1_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     IF (S_Flag.eq.1) then
                        call MPI_RECV(theta(1,1,2),(nzb)*nx*ny, 
     +                       MPI_DOUBLE_PRECISION,ii,ii,nall,
     +                       status2,ierr)
                     ENDIF
                     IF (Q_Flag.eq.1) then
                        call MPI_RECV(q(1,1,2),(nzb)*nx*ny, 
     +                       MPI_DOUBLE_PRECISION,ii,ii,nall,
     +                       status2,ierr)
                     ENDIF
                     
                     do k=2,Nzb+1
                        do j=1,Ny
                           do i=1,Nx
                              
                              write (84,6623) u(i,j,k)+Ugal,v(i,j,k)+
     +                             Vgal,w(i,j,k)
	 
                              write(20) a1_old(i,j,k),b1_old(i,j,k),
     +                                  c1_old(i,j,k),d1_old(i,j,k),
     +                                  e1_old(i,j,k)	 
                              
                              IF (S_Flag.eq.1) then
                                 write (108,*) theta(i,j,k)
                              ENDIF
                              IF (Q_Flag.eq.1) then
                                 write (109,*) q(i,j,k)
                              ENDIF
                           end do
                        end do
                     end do
                     
                  END DO
                  
               endif
               
!               write(20) u_frame,v_frame,w_frame,t_frame,q_frame
               close(20)
               close(84)

	       IF(S_flag.eq.1) close(108)
	       IF(Q_flag.eq.1) close(109)
               
c     restore u field
               u=dudz
               v=dvdz
               w=dwdz

               a1_old=dudx
               b1_old=dvdx
               c1_old=dwdx			   
               d1_old=dudy
               e1_old=dvdy				   
			   
			   
               IF (S_Flag.eq.1) then
                  theta=dtdz
               endif
               IF (Q_Flag.eq.1) then
                  q=dqdz
               endif
               
            endif
            
            IF (me>0) then

               call MPI_SEND(a2_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,0,me,nall, ierr )
               call MPI_SEND(b2_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,0,me,nall, ierr )
               call MPI_SEND(c2_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,0,me,nall, ierr )
               call MPI_SEND(d2_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,0,me,nall, ierr )
               call MPI_SEND(e2_old(1,1,2),(nzb)*nx*ny,
     +              MPI_DOUBLE_PRECISION,0,me,nall, ierr )

            else
               
               Open (unit=21,file='input/a2_e2.bin',status='unknown',
     +              form='unformatted')
               do k=2,Nzb+1
                  do j=1,Ny
                     do i=1,Nx
                              write(21) a2_old(i,j,k),b2_old(i,j,k),
     +                                  c2_old(i,j,k),d2_old(i,j,k),
     +                                  e2_old(i,j,k)	 
                     enddo
                  enddo
               enddo

               dudx=a2_old
               dvdx=b2_old
               dwdx=c2_old			   
               dudy=d2_old
               dvdy=e2_old				   
			   
			   
               IF(nprocs.gt.1)THEN
                  DO ii=1,nprocs-1
                     call MPI_RECV(a2_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(b2_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(c2_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(d2_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(e2_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)

                     do k=2,Nzb+1
                        do j=1,Ny
                           do i=1,Nx
						   
                              write(21) a2_old(i,j,k),b2_old(i,j,k),
     +                                  c2_old(i,j,k),d2_old(i,j,k),
     +                                  e2_old(i,j,k)                              

                           end do
                        end do
                     end do
			   
                  END DO
               ENDIF			   
			   
!               write(21) u_frame,v_frame,w_frame,t_frame,q_frame
               close(21)
			   
               a2_old=dudx
               b2_old=dvdx
               c2_old=dwdx			   
               d2_old=dudy
               e2_old=dvdy				   
			   
            endif
            
            IF (me>0) then

                  call MPI_SEND(RHSx(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  call MPI_SEND(RHSy(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  call MPI_SEND(RHSz(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )

               else

                  Open (unit=25,file='input/RHS_mom.bin',
     +                 status='unknown',form='unformatted')

                  do k=2,Nzb+1
                     do j=1,Ny
                        do i=1,Nx
                              write(25) RHSx(i,j,k),RHSy(i,j,k),
     +                                  RHSz(i,j,k)
                        enddo
                     enddo
                  enddo
				  
               dudz=RHSx
               dvdz=RHSy
               dwdz=RHSz
			   
                  IF(nprocs.gt.1)THEN
                     DO ii=1,nprocs-1
                        call MPI_RECV(RHSx(1,1,2),(nzb)*nx*ny,
     +                       MPI_DOUBLE_PRECISION,ii,ii,nall,
     +                       status2,ierr)
                        call MPI_RECV(RHSy(1,1,2),(nzb)*nx*ny,
     +                       MPI_DOUBLE_PRECISION,ii,ii,nall,
     +                       status2,ierr)
                        call MPI_RECV(RHSz(1,1,2),(nzb)*nx*ny,
     +                       MPI_DOUBLE_PRECISION,ii,ii,nall,
     +                       status2,ierr)
	 

                  do k=2,Nzb+1
                     do j=1,Ny
                        do i=1,Nx
                              write(25) RHSx(i,j,k),RHSy(i,j,k),
     +                                  RHSz(i,j,k)
                        enddo
                     enddo
                  enddo 
	 
	 
                     END DO
                  ENDIF

!                  write(25) u_frame,v_frame,w_frame
                  close(25)
				  
               RHSx=dudz
               RHSy=dvdz
               RHSz=dwdz				  
               endif

            IF (S_Flag.eq.1) then
               
               IF (me>0) then
                  
                  call MPI_SEND(a4_old(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  call MPI_SEND(b4_old(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  call MPI_SEND(c4_old(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  call MPI_SEND(d4_old(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  call MPI_SEND(e4_old(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  
               else
                  
                  Open (unit=22,file='input/a4_e4.bin',status='unknown',
     +                 form='unformatted')
                  
               do k=2,Nzb+1
                  do j=1,Ny
                     do i=1,Nx
                              write(22) a4_old(i,j,k),b4_old(i,j,k),
     +                                  c4_old(i,j,k),d4_old(i,j,k),
     +                                  e4_old(i,j,k)	 
                     enddo
                  enddo
               enddo

               dudx=a4_old
               dvdx=b4_old
               dwdx=c4_old			   
               dudy=d4_old
               dvdy=e4_old				   
			   
			   
               IF(nprocs.gt.1)THEN
                  DO ii=1,nprocs-1
                     call MPI_RECV(a4_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(b4_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(c4_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(d4_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(e4_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)

                     do k=2,Nzb+1
                        do j=1,Ny
                           do i=1,Nx
						   
                              write(22) a4_old(i,j,k),b4_old(i,j,k),
     +                                  c4_old(i,j,k),d4_old(i,j,k),
     +                                  e4_old(i,j,k)                              

                           end do
                        end do
                     end do
			   
                  END DO
               ENDIF			   
			   
!               write(22) u_frame,v_frame,w_frame,t_frame,q_frame
               close(22)
			   
               a4_old=dudx
               b4_old=dvdx
               c4_old=dwdx			   
               d4_old=dudy
               e4_old=dvdy	

               endif
               
               IF (me>0) then
                  
                  call MPI_SEND(a8_old(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  call MPI_SEND(b8_old(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  call MPI_SEND(c8_old(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  call MPI_SEND(d8_old(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  call MPI_SEND(e8_old(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  
               else
                  
                  Open (unit=23,file='input/a8_e8.bin',status='unknown',
     +                 form='unformatted')
                  
               do k=2,Nzb+1
                  do j=1,Ny
                     do i=1,Nx
                              write(23) a8_old(i,j,k),b8_old(i,j,k),
     +                                  c8_old(i,j,k),d8_old(i,j,k),
     +                                  e8_old(i,j,k)	 
                     enddo
                  enddo
               enddo

               dudx=a8_old
               dvdx=b8_old
               dwdx=c8_old			   
               dudy=d8_old
               dvdy=e8_old				   
			   
			   
               IF(nprocs.gt.1)THEN
                  DO ii=1,nprocs-1
                     call MPI_RECV(a8_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(b8_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(c8_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(d8_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(e8_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)

                     do k=2,Nzb+1
                        do j=1,Ny
                           do i=1,Nx
						   
                              write(23) a8_old(i,j,k),b8_old(i,j,k),
     +                                  c8_old(i,j,k),d8_old(i,j,k),
     +                                  e8_old(i,j,k)                              

                           end do
                        end do
                     end do
			   
                  END DO
               ENDIF			   
			   
!               write(23) u_frame,v_frame,w_frame,t_frame,q_frame
               close(23)
			   
               a8_old=dudx
               b8_old=dvdx
               c8_old=dwdx			   
               d8_old=dudy
               e8_old=dvdy	
               endif

               IF (me>0) then

                  call MPI_SEND(RHS_T(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )

               else

                  Open (unit=26,file='input/RHS_T.bin',status='unknown',
     +                 form='unformatted')

                  do k=2,Nzb+1
                     do j=1,Ny
                        do i=1,Nx
                              write(26) RHS_T(i,j,k)						
                        enddo
                     enddo
                  enddo
				  
               dudz=RHS_T				  
				  
                  IF(nprocs.gt.1)THEN
                     DO ii=1,nprocs-1
                        call MPI_RECV(RHS_T(1,1,2),(nzb)*nx*ny,
     +                       MPI_DOUBLE_PRECISION,ii,ii,nall,
     +                       status2,ierr)
	 
                  do k=2,Nzb+1
                     do j=1,Ny
                        do i=1,Nx
                              write(26) RHS_T(i,j,k)						
                        enddo
                     enddo
                  enddo	 
	 
                     END DO
                  ENDIF

!                  write(26) t_frame
                  close(26)
               RHS_T=dudz				  
				  
               endif

            ENDIF

            IF (Q_Flag.eq.1) then
               
               IF (me>0) then
                  
                  call MPI_SEND(a5_old(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  call MPI_SEND(b5_old(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  call MPI_SEND(c5_old(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  call MPI_SEND(d5_old(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  call MPI_SEND(e5_old(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  
               else
                  
                  Open (unit=24,file='input/a5_e5.bin',status='unknown',
     +                 form='unformatted')
                  
              do k=2,Nzb+1
                  do j=1,Ny
                     do i=1,Nx
                              write(24) a5_old(i,j,k),b5_old(i,j,k),
     +                                  c5_old(i,j,k),d5_old(i,j,k),
     +                                  e5_old(i,j,k)	 
                     enddo
                  enddo
               enddo

               dudx=a5_old
               dvdx=b5_old
               dwdx=c5_old			   
               dudy=d5_old
               dvdy=e5_old				   
			   
			   
               IF(nprocs.gt.1)THEN
                  DO ii=1,nprocs-1
                     call MPI_RECV(a5_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(b5_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(c5_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(d5_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(e5_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)

                     do k=2,Nzb+1
                        do j=1,Ny
                           do i=1,Nx
						   
                              write(24) a5_old(i,j,k),b5_old(i,j,k),
     +                                  c5_old(i,j,k),d5_old(i,j,k),
     +                                  e5_old(i,j,k)                              

                           end do
                        end do
                     end do
			   
                  END DO
               ENDIF			   
			   
!               write(24) u_frame,v_frame,w_frame,t_frame,q_frame
               close(24)
			   
               a5_old=dudx
               b5_old=dvdx
               c5_old=dwdx			   
               d5_old=dudy
               e5_old=dvdy	
               endif
               
               IF (me>0) then
                  
                  call MPI_SEND(a9_old(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  call MPI_SEND(b9_old(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  call MPI_SEND(c9_old(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  call MPI_SEND(d9_old(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  call MPI_SEND(e9_old(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )
                  
               else
                  
                  Open (unit=27,file='input/a9_e9.bin',status='unknown',
     +                 form='unformatted')
                  
              do k=2,Nzb+1
                  do j=1,Ny
                     do i=1,Nx
                              write(27) a9_old(i,j,k),b9_old(i,j,k),
     +                                  c9_old(i,j,k),d9_old(i,j,k),
     +                                  e9_old(i,j,k)	 
                     enddo
                  enddo
               enddo

               dudx=a9_old
               dvdx=b9_old
               dwdx=c9_old			   
               dudy=d9_old
               dvdy=e9_old				   
			   
			   
               IF(nprocs.gt.1)THEN
                  DO ii=1,nprocs-1
                     call MPI_RECV(a9_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(b9_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(c9_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(d9_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
                     call MPI_RECV(e9_old(1,1,2),(nzb)*nx*ny, 
     +                    MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)

                     do k=2,Nzb+1
                        do j=1,Ny
                           do i=1,Nx
						   
                              write(27) a9_old(i,j,k),b9_old(i,j,k),
     +                                  c9_old(i,j,k),d9_old(i,j,k),
     +                                  e9_old(i,j,k)                              

                           end do
                        end do
                     end do
			   
                  END DO
               ENDIF			   
			   
!               write(27) u_frame,v_frame,w_frame,t_frame,q_frame
               close(27)
			   
               a9_old=dudx
               b9_old=dvdx
               c9_old=dwdx			   
               d9_old=dudy
               e9_old=dvdy	
               endif

               IF (me>0) then

                  call MPI_SEND(RHS_Q(1,1,2),(nzb)*nx*ny,
     +                 MPI_DOUBLE_PRECISION,0,me,nall, ierr )

               else

                  Open (unit=28,file='input/RHS_Q.bin',status='unknown',
     +                 form='unformatted')

                  do k=2,Nzb+1
                     do j=1,Ny
                        do i=1,Nx
                              write(28) RHS_Q(i,j,k)						
                        enddo
                     enddo
                  enddo
				  
               dudz=RHS_Q				  
				  
                  IF(nprocs.gt.1)THEN
                     DO ii=1,nprocs-1
                        call MPI_RECV(RHS_Q(1,1,2),(nzb)*nx*ny,
     +                       MPI_DOUBLE_PRECISION,ii,ii,nall,
     +                       status2,ierr)
	 
                  do k=2,Nzb+1
                     do j=1,Ny
                        do i=1,Nx
                              write(28) RHS_Q(i,j,k)						
                        enddo
                     enddo
                  enddo	 
	 
                     END DO
                  ENDIF

!                  write(28) t_frame
                  close(28)
               RHS_Q=dudz	
               endif

            ENDIF
            
         ENDIF
         
!     ===========================================E
!     ... End Time Loop  ===================E
!     ... End Time Loop  ==============E
!     ... End Time Loop  ===================E          
!     ===========================================E
         
 728  ENDDO


 5166 format (5(1x,e15.6))    
 5167 format (7(1x,e15.6))      
 5111 format (1x,i10,15(1x,f11.5))  
 5555 format(1x,a10,1x,i10,15(1x,f11.5)) 
      
      if(me==0)then 
         t02=MPI_WTIME()
         call closefiles()
      endif

      if(me.eq.0)then
         open(unit=110,file='input/zo.out',
     +        status='unknown')
         do j=1,ny
            do i=1,nx
               write(110,*) zo(i,j)
            enddo
         enddo
         close(110)
      endif
      if(S_flag.eq.1)then
         if(me.eq.0)then
            if(surf_flag.eq.1)then
               open(unit=781,file='input/coolrate.out',status='unknown')
               open(unit=111,file='input/t_surf.out',status='unknown')
               do j=1,ny
                  do i=1,nx
                     write(781,*) coolrate(i,j)
                     write(111,*) t_s(i,j)*T_scale
                  enddo
               enddo
               close(111)
            elseif(surf_flag.eq.2)then
               open(unit=111,file='input/t_025meters.out',
     +              status='unknown')
               do i=1,nhrs
                  write(111,*) t_s2(i)*T_scale
               enddo
               close(111)
            endif
         endif
      endif
      if(Q_flag.eq.1)then
         if(me.eq.0)then
            if(qsurf_flag.eq.1)then
               open(unit=782,file='input/q_surf.out',
     +              status='unknown')
               open(unit=783,file='input/qcoolrate.out',
     +              status='unknown')
               do j=1,ny
                  do i=1,nx
                     write(783,*) qcoolrate(i,j)
                     write(782,*) q_s(i,j)*Q_scale
                  enddo
               enddo
               close(782)
               close(783)
            elseif(qsurf_flag.eq.2)then
               open(unit=782,file='input/q_025meters.out',
     +              status='unknown')
               do i=1,nhrs
                  write(782,*) q_s2(i)*Q_scale
               enddo
               close(782)
            endif
         endif
      endif

      IF (me>0) then
         
         call MPI_SEND(u(1,1,2),(nzb)*nx*ny,MPI_DOUBLE_PRECISION,
     +        0,me,nall, ierr )
         call MPI_SEND(v(1,1,2),(nzb)*nx*ny,MPI_DOUBLE_PRECISION,
     +        0,me,nall, ierr )
         call MPI_SEND(w(1,1,2),(nzb)*nx*ny,MPI_DOUBLE_PRECISION,
     +        0,me,nall, ierr )
         

         
         if(S_Flag.eq.1)then
            call MPI_SEND(theta(1,1,2),(nzb)*nx*ny, 
     +           MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         endif
         if(Q_Flag.eq.1)then
            call MPI_SEND(q(1,1,2),(nzb)*nx*ny, 
     +           MPI_DOUBLE_PRECISION,0,me,nall, ierr )
         endif

      else
         
         Open (unit=84,file='input/vel.out',status='unknown')
         IF (S_Flag.eq.1) then
            Open (unit=108,file='input/t.out',status='unknown')
         ENDIF
         IF (Q_Flag.eq.1) then
            Open (unit=109,file='input/q.out',status='unknown')
         ENDIF
         write (84,*) ttt, dt
         write (84,*) Nx,Ny,Nz
         write (84,*)
         
         do k=2,Nzb+1
            do j=1,Ny
               do i=1,Nx
                  
                  write (84,6623) u(i,j,k)+Ugal,v(i,j,k)+Vgal,w(i,j,k)
                  
                  IF (S_Flag.eq.1) then
                     write (108,*) theta(i,j,k)
                  ENDIF
                  IF (Q_Flag.eq.1) then
                     write (109,*) q(i,j,k)
                  ENDIF
               end do
            end do
         end do
         
         dudz=u
         dvdz=v
         dwdz=w
         
         IF (S_Flag.eq.1) then
            dtdz=theta
         ENDIF
         IF (Q_Flag.eq.1) then
            dqdz=q
         ENDIF

         
         IF(nprocs.gt.1)THEN
            
            DO ii=1,nprocs-1
               
               call MPI_RECV(u(1,1,2),(nzb)*nx*ny, 
     +              MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
               call MPI_RECV(v(1,1,2),(nzb)*nx*ny, 
     +              MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
               call MPI_RECV(w(1,1,2),(nzb)*nx*ny, 
     +              MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
               
               IF (S_Flag.eq.1) then
                  call MPI_RECV(theta(1,1,2),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
               ENDIF
               IF (Q_Flag.eq.1) then
                  call MPI_RECV(q(1,1,2),(nzb)*nx*ny, 
     +                 MPI_DOUBLE_PRECISION,ii,ii,nall,status2,ierr)
               ENDIF
            
               do k=2,Nzb+1
                  do j=1,Ny
                     do i=1,Nx
                        
                        write (84,6623) u(i,j,k)+Ugal,v(i,j,k)+Vgal,
     +                       w(i,j,k)
                        
                        IF (S_Flag.eq.1) then
                           write (108,*) theta(i,j,k)
                        ENDIF
                        IF (Q_Flag.eq.1) then
                           write (109,*) q(i,j,k)
                        ENDIF
                     end do
                  end do
               end do
            
            END DO

         endif
         
         close(108)
         close(84)
         
!     restore u field
         u=dudz
         v=dvdz
         w=dwdz

         IF (S_Flag.eq.1) then
            theta=dtdz
         endif
         IF (Q_Flag.eq.1) then
            q=dqdz
         endif
         
      endif

      deallocate(rewards_old)    
      deallocate(actions_2D)    
      deallocate(states_2D)    
            
 6623 format (24(f11.5,1x))
      
      if (me==0) then
         write(*,*) 'Total time = ', t02-t01
      end if

      endif 
      
 999  result_value = 0
      
      if (train==0) call MPI_FINALIZE(ierr)

      print *,"Fortran side ends",me

      if (me == 0) then
      if (RL_init == 0) then
        open(unit=55,file = 'RL_init',status="new")
        write(55,*) 1
        close(55)
      else if (RL_init.gt.0) then
        open(unit=55,file = 'RL_init',status="old")
        write(55,*) RL_init+1
        close(55)
      end if
      print *,"end-RL_init",RL_init
      end if  

      if (train == 0 ) then
        call exit(1)
      end if
       
      END function app_main
      
      character(len=20) function str(k)
!     "Convert an integer to string."
      integer, intent(in) :: k
      write (str, *) k
      str = adjustl(str)
      end function str

      END module app_main_module   

!=================================================================

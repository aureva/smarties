ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ------------------------
c Notes on the Code (general description)::
c ------------------------
c This code solves the filtered conservation of momentum and the filtered
c conservation of scalar concentration. Its basic form and numerics are
c described in Albertson 1996 (Phd thesis UC Davis), Albertson and Parlange 1999
c (Water Resources Res.) and Albertson and Parlange 1999 (Adv. in Water Res.).
c The following is a brief discription of some of the key points in the code.
c
c SGS model:     Scale-dependent dynamic model (scalars and momentum)
c                described in Porte-Agel et al 2000 (JFM), Porte-Agel 2004 (BLM). The
c                averaging scheme for SGS optimization is done with either 
c                Stoll and Porte-Agel 2005 (WRR) or Basu and Porte-Agel 2005 (JAS) 
c                Wong-Lilly base model is also incorporated as an option 
c                (Wong-Lilly, Phys Fluids, 1994)
c
c Momentum BC:   local log-law (see Stoll and Porte-Agel 2005 (BLM)
c 	         for description of Momentum BCs) at surface, stress free 
c                condition at domain top and periodic in horizontal directions.
c
c Scalar BC:     local log-law with constant cooling/heating rate or
c                prescribed constant flux at the surface, constant inversion at 
c                domain top and periodic in horizontal directions.
c
c Numerics:      Spectral methods are used in the horizontal directions and
c                for the SGS filtering operations.  2nd order central differences
c                are used in the vertical direction for derivatives.  The
c                Fourier transforms are computed by the fftw libraries version 3.0.1
c                (see www.fftw.org for algorithm details and syntax).
c
c -------------------------
c Folder and file Structure
c -------------------------
c File extensions are as follows: *.bin for all files in Fortran binary format
c and *.out for all files in ASCII format (this includes input and output files).
c Some input files may be labeled *.ini, these files are initialization files
c and should be copied to *.out before execution if desired (of course you don't 
c need to do this for continuation of simulations etc.).
c
c The code files are split into three sections. :
c
c 1) all of the *.f files, dimen.h, makefile, jobsub, resubmit, les.out and les.err
c files are contained in the main directory.  les.out is the standard output for IBM systems
c and les.err is the standard error file.  For the Altix the standard output is jobsub.ojobid
c and the standard error is jobsub.ejobid where jobid is the system assigned jobid number.
c
c 2) the input files are contained in the directory input\ 
c
c      This directory contains at least the following files (if running only 
c      momentum only vel.out and zo.out are required):
c
c        vel.out     standard input/output velocity field file
c        t.out       standard input/output scalar field file
c        zo.out      ASCII surface roughness file (to be described later MATLAB Section)
c        t_surf.out  ASCII surface temperature file (to be described later MATLAB Section)
c       
c The following files are optional (depending on initialization options in dimen.h) 
c and may or may not be contained within
c
c       t.ini       standard name for a 'backed up initialization file' 
c                   (must be copied to t.out to use)
c       vel.ini     standard name for a 'backed up initialization file'
c                   (must be copied to vel.out to use)
c       t_surf.ini  standard name for a 'backed up initialization file'
c                   (must be copied to t_surf.out to use)
c       a1_e1.bin   Lagrangian SGS binary initialization file for a1_old-e1_old (momentum)
c       a2_e2.bin   Lagrangian SGS binary initialization file for a2_old-e2_old (momentum)
c       a4_e4.bin   Lagrangian SGS binary initialization file for a4_old-e4_old (scalars)
c       a8_e8.bin   Lagrangian SGS binary initialization file for a8_old-e8_old (scalars)
c
c 3) all of the output files are contained the directory output\.
c
c IMPORTANT NOTE: the output files will be appended unless inituA=0 and initsA=0
c (init params are described later) in which case the output files will
c be erased.  This means that if you don't want the output files to contain
c data form the last run you need to manually erase (or move) all of the files
c in the output folder.   
c
c ---------------------------
c Input parameters in dimen.h
c ---------------------------
c
c This code has the following input parameters grouped by function 
c (initialization, domain, etc.)
c
c INITIALIZATION:
c
c -resub     This number should match the number in the file resubmit
c            It tells the code how many times to run repeatedly.  Note that
c            it will only work if the code exits normally (ie not killed by wallclock limit). 
c            If resub is present and non zero it will override inituA and initsA.
c
c -inituA    If inituA=1 then the code will read vel.out,and all 
c            Lagrangian binary initialization files.  This is used to restart a
c            previous simulation from its last state. if initu is equal to 0 
c            the code will still read vel.out but it will NOT read the 
c            binary initialization files it will instead start by 
c            plane averaging for the first few timesteps. For Local/Plane average
c            inituA should be set to 0.
c
c -initsA    Same as for inituA but for scalars (reads binary initialization files).
c            A value of 1 is for a continuation and a value of 0 will read in 
c            a t.out file but NOT the binary initialization files this option 
c            with start the SGS model by plane averaging. For Local/Plane average
c            initsA should be set to 0. Note that this will also apply to q if
c            Q_FLAG=1.
c
c -ruler     Controls how often vel.out,t.out and t_surf.out are updated
c            
c SIMULATION DOMAIN, COMPUTATIONAL GRID AND SIMULATION TIME:
c
c -Nx,Ny,Nz  The number of grid points in the streamwise, spanwise and 
c            wall-normal directions respectively.
c
c -Nz2       Number of levels given to each processor where Nz2=Nz/nprocs+2
c            with nprocs=number of processors used in the simulation.
c
c -L_z       dimensioned depth of the domain
c
c -z_i       normalization height (Dimensionless depth  of domain=L_z/z_i 
c	    (recall Lx/z_i=2*PI; Ly/z_i=2*PI/l_r))
c
c -l_r       ratio of streawise to spanwise dimensions (l_r=Lx/Ly) this should
c            be an interger value.
c      
c -Nsteps    total number of timesteps
c
c -dt        non dimensional timestep (dt=dtr/(z_i/u_star), where dtr is 
c            the physical timestep in seconds).
c
c OUTPUTTING PARAMS:
c -aNx       Controls the dimension of the outputting of statistics. A value of
c	     one results in plane-averaged output (stats are only f(z)) and a 
c	     value of Nx results in spanwise averaged output (stats are f(x,z)).
c
c -c_count   frequency of statistical sampling (for example c_count=10 
c            samples stats every 10th timestep)
c
c -p_count   frequency of outputting statistic (for example p_count=100 averages
c	     c_count stats over 100 timesteps and then outputs them to files).
c
c SGS PARAMETERS:
c -cs_count  frequency that SGS coefficients are updated
c
c -model     model equals: 1->Smagorinsky; 2->Dynamic; 3->Scale dependent 
c            (possibly only 3 works)
c
c -averaging the averaging scheme for sgs optimization, 
c            averaging equals: 
c            0-> Plane average (Porte-Agel et al., JFM 2000; Porte-Agel, BLM, 2004)
c            1-> Lagrangian average (Stoll-Porte-Agel 2005, WRR), 
c            2-> Local average (Basu-Porte-Agel 2005, JAS)
c            3-> Wong-Lilly scale-dependent dynamic model (Wong and Lilly, Phys Fluids, 1994)
c
c -fgr       filter-to-grid ratio
c
c -tfr       test filter to grid filter ratio
c
c -mom_nodes computational level at which the dynamic momentum optimization is done: 
c            0 for w-nodes and 1 for uvp-nodes
c
c -scl_nodes computational level at which the dynamic scalar optimization is done: 
c            0 for w-nodes and 1 for uvp-nodes
c
c SCALAR PARAMETERS:
c -S_FLAG    S_FLAG=1 for including scalar computation set = 0 for NO scalars
c
c -Q_FLAG    Q_FLAG=1 for including a second scalar (e.g. water vapor)
c
c -pass_flag zero for active scalars 1 for passive
c
c -T_scale   Normalization temperature
c
c -Q_scale   Normalization factor for q (e.g., moisture content)
c
c -theta_0   buoyancy flux reference temperature
c
c SURFACE BOUNDARY CONDITIONS, MEAN FORCING AND DAMPENING LAYER:
c -press_cor zero for constant pressure forcing 1 for geostrophic wind
c
c -Ugeo_inf,Vgeo_inf  geostrophic wind components at domain top (for press_cor=1
c                     this value is constant throughout, isobaric ABLS)
c
c -Ugal      Gallilean transformation (should =0 for press_cor=0)
c
c -f_c       Coriolis parameter f_c=1.45*sin(latitude)
c
c -inversion inversion strength at domain top (lapse rate)
c
c -qinversion inversion strength at domain top for q
c
c -surf_flag zero for constant flux 1 for cooling rate
c
c -qsurf_flag same as above for q
c
c -s_flux    value of constant flux if surf_flag=0
c
c -q_flx     value of constant flux if qsurf_flag=0
c
c -c_coolrate cooling rate (K/hr) for surf_flag=1 
c             (see also coolrate.out defined in input file)
c
c -q_coolrate change in surface concentration of q with time for qsurf_flag=1
c
c -sponge    set sponge=1 for dampening layer starting at height z_d else none
c
c -rlx_time  relaxation time for dampening layer in seconds
c
c -Ri_flag   Richardson number criteria: if the pointwise gradient Richardson
c            number exceeds the pointwise SGS Prandtl number, the SGS terms are
c            switched off. USE THIS FLAG WITH CAUTION. set Ri_flag = 1 if you
c            wish to invoke this crtieria.
c
c 3D FIELD (FRAMES) OUTPUTTING:
c This outputs 3D fields of u,v,w and theta (if S_FLAG=1)
c
c -nframe    number of instantaneous 'frames' to output (a value of 0 will
c            skip all outputting of frames)
c -sframes   the number of timesteps at which to start saving frames   
c -framestep number of timesteps in between frames
c -frameH    how many levels in the vertical direction to include
c	
c --------------
c Job submission
c --------------
c
c Makeing the code:
c
c makefile    This is the file to compile the code. If a change to dimen.h (this file)
c             has been made or the code has not been compiled on a paticular machine
c             then use as:
c                   $ make clean
c                   $ make
c              if only individual *.f files have been altered use just:
c                   $ make
c              Note that the makefile must be altered for IBM AIX and SGI Altix
c              Machines.  The proper lines to comment/uncomment are listed in the makefile
c              where a comment is denoted by a # symbol.
c
c Submission to the Queue:
c
c jobsubSGI   $qsub jobsubSGI
c             a pbs batch submission script for running on the SGI Altix.
c             The max wallclock is 150 hrs and nodes have 16 procs (one with 256 procs).
c             The 256 proc node has 512 Gb of memory and the 8-16 proc nodes each have 32 Gb.
c             NOTE: make sure to check the start directory (startpath variable).
c
c jobsubIBM   llsubmit jobsubIBM
c             the batch job submission script for IBM systems contained in the main.
c             The max wallclock is 150 hrs and nodes have either 8 or 32 procs.
c             For a total of 344 procs with 680 Gb of memory.
c             NOTE: make sure to check the start path!!!
c
c jobsubBlade qsub jobsubBlade (production queue); qsub -q devel jobsubBlade (development queue)
c             a pbs batch submission script for running on the IBM Blade Center Linux
c             cluster.  The max wallclock is 24hrs for 1072 cores (processors) and 64 cores
c             have a 1 hour wallclock for development purposes.  Each node has 4 cores so the
c             parameter ppn in jobsubBlade should be <= 4.  Each node has 7Gb usable memory.
c             Note the options needed (see 1st line for jobsubBlade) to use the short queue.
c             Note: make sure the cd command has the right directory (full path).
c
c resubmit    a file containing a single ASCII number telling the code how many
c             times to resubmit. If resubmit does not exist it won't happen.
c
c Interactive running:
c
c interactive running uses the poe command (on AIX systems) as:
c
c            $ poe ./LES2 -tasks_per_node # -nodes # -rmpool 1
c
c where the # denotes a numeric value.
c NOTE:  The regatta has 1 interactive node with 8 procs and
c
c interactive running using the mpirun command (on SGI Altix):
c
c            $ module load fftw
c            $ mpirun -np # ./LES2
c
c where the # is the number of procs.
c NOTE:  the Altix has one 4 processor interactive node
c
c interactive running using the mpirun command (on IBM Blade):
c
c            $ mpirun -np # blade285 ./LES2
c
c Note # is the number of procs and the system has 4 interactive nodes each with 4 cores.
c For each processes (core) a intective node must be specified.  The interactive nodes
c are blade285, blade286, blade287, blade288 so the % above should be 5, 6, 7 or 8.  More than one node can
c be specified in the list, for example 4 procs total on 2 interactive nodes:
c
c            $ mpirun -np 4 blade287 blade287 blade288 blade288 ./LES2
c
c -------------------
c
c -------------------
c MATLAB AND THE CODE
c -------------------
c Although not the only way to work with the code, MATLAB is the preferred method
c for making initialization files and looking at code output.  For general help 
c with MATLAB see one of the Indian gurus (Sukanta or Venu) or an American 
c knock off (Rob or Matt).
c
c Two mfiles are included in the main directory (usually unless someone forgot
c to copy them).
c
c LESpatch.m     this mfile will create the ASCII delimited zo.out and t_surf.out
c                files needed for the BCs
c
c loadbin.m      function to load the binary output files (with the exception
c                of the frames).
c 
c Two mfiles to load the 3D frames data are also supplied as examples 
c (loadvel.m and loadtheta.m).
c
c Mfiles for plotting the output can be written by the user or obtained from a
c current LES user
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccc END CODE DESCRIPTION INCLUDE FILE FOLLOWS ccccccccccccccccccccc
c###
      include 'mpif.h'

      integer*4 Nx,Ny,Nz,Nx2,Ny2,aNx
      integer*4 nsteps,c_count,p_count,cs_count,t_count
      integer*4 frameH,nframe,sframe,framestep,nhrs
      integer*4 S_FLAG,Q_FLAG,pass_flag,qpass_flag,model
      integer*4 initu,inits,inituA,initsA,resub,nrsub      
      integer*4 press_cor,surf_flag,sponge,Ri_flag,ruler
      integer*4 mom_nodes,scl_nodes,averaging,qsurf_flag

	  integer*4 vel_instant,buffer,Nxbm,Nxbs,Nxbe,resubnum
      integer*4 T_instant,bufferT,NsT,Nxx	  
	  
      real*8 dx,dy,dz,delta,l_r,u2_surface
      real*8 Pi,Sc,u_star,vonk,Ugal,Vgal,g_hat,u_scale,nu,Pr
      real*8 theta_0,s_flux,qs_flux,t_scale,inversion
      real*8 Q_scale,cq_coolrate,qinversion,M_advec
      real*8 L_z,z_i,fc,Ugeo_inf,Vgeo_inf,f_c,z_d,rlx_time
      real*8 nnn,Co,dt,fgr,tfr,c_coolrate,T_advec,Q_advec
      real*8 idz,idy,idx,inxny,inx2ny2,dtl
	  
      real*8 zref

      integer*4 nzb, nz2, nprocs, me,nall,nps
      integer*4 ierr, status2(MPI_STATUS_SIZE)
      common /rune/nprocs, nzb,resubnum
      common /calc_dims/ dx,dy,dz,delta
      common /calc_norm/ idz,idy,idx,inxny,inx2ny2,dtl
      common /calc_pars/ g_hat,fc
      common /init_pars/ initu,inits,nrsub

cccc fixed parameters
      parameter (vonk=0.4d0)
      parameter (Sc=0.4,Co=0.1,nnn=2.)
      parameter (u_star=1.0d0,u_scale=1.0d0,u2_surface=1.0d0**2)
      parameter (Pi=3.14159265358979d0) 

cccc input parameters
      parameter (resub=0,inituA=0,initsA=0,ruler=48*20*5)
      parameter (Nx=48,Ny=48,Nz=48)
      parameter (nz2=5,aNx=1,l_r=1.0d0,nps=16)
      parameter (Nsteps=48*20*200)
	  
      parameter (c_count=10, p_count=48*20*20, t_count=25, cs_count=3)
      parameter (S_FLAG=0,Q_FLAG=0,pass_flag=1,qpass_flag=1)
      
      parameter (model=100)
      parameter (averaging=1)
      parameter (fgr=1.0d0)      
      parameter (tfr=2.0d0)
      parameter (mom_nodes=1,scl_nodes=1)

      parameter (z_i=1.0d0*Pi/Pi,l_z=1.0d0)
      parameter (nu=1.d-5/(z_i*u_star),Pr=0.7d0)
	  
      parameter (T_Scale = 1.d0, Q_Scale = 1.d0)
      parameter (theta_0 = 293.d0/T_Scale)
      parameter (s_flux = 0.03185d0, c_coolrate=0.0d0)
      parameter (qs_flux = 0.0d0, cq_coolrate=0.0d0)      
      parameter (surf_flag = 1, qsurf_flag = 1, nhrs=24)
      parameter (sponge = 0, z_d = 300.d0)
      parameter (rlx_time=1.0d0, Ugal=0.0d0/u_star, Vgal=0.0d0/u_star)
      parameter (Ri_flag=0)
	  
      parameter (press_cor=0, f_c=1.39d0)
      parameter (M_advec=0, T_advec=0, Q_advec=0)
      parameter (Ugeo_inf=8.0d0/u_star, Vgeo_inf=0.0d0/u_star)
      parameter (inversion=0.0d0*z_i/T_scale)
      parameter (qinversion=0.0d0*z_i/Q_scale)

      parameter (sframe=10000,nframe=0,framestep=80,frameH=48)
c.....calculated parameters......
      parameter (nx2=nx*3/2, ny2=ny*3/2)

      parameter (zref=0.1d0/z_i)
      parameter (dt=1.1450d-4*1.2d0/z_i)	  
	  
cc      parameter (dt=1.50d-4)
cc      parameter (dt=1.d0/(799.3d0*0.50d0/60.d0*z_i/u_star)/360)	    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer*4 Nxs,Nxe,Nys,Nye,Nzs,Nze
	  
      parameter (Nxs= 1,Nys= 1,Nzs= 1)
      parameter (Nxe=Nx,Nye=Ny,Nze=Nz)
	  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     wind turbine parameters
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer*4 num_turbine,turbine_model,nacelle_model,tower_model
      parameter (num_turbine=0)

      integer*4 wtx(num_turbine),wtz(num_turbine)
      real*8 wtR(num_turbine),Zhub(num_turbine),wtomega(num_turbine),
     +       wty(num_turbine),dir(num_turbine)   
 
      parameter (turbine_model=0,nacelle_model=0,tower_model=0)

      integer*4 N_ang,N_blade,N_blade_s
      parameter (N_ang=30,N_blade=40,N_blade_s=1)

      real*8 CTN(num_turbine),CTB(num_turbine),CTT(num_turbine)
      real*8 ALFA_CL_CD(3,181)
	  
      common /turbine_pars1/ wtx,wtz,ALFA_CL_CD
      common /turbine_pars2/ wtR,Zhub,wtomega,CTN,CTB,CTT,wty,dir
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     wind turbine parameters
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer*4 num_turbine_vawt,turbine_model_vawt,tower_model_vawt
      parameter (num_turbine_vawt=0)

      integer*4 wtx_vawt(num_turbine_vawt),wtz_vawt(num_turbine_vawt)
      real*8 wtR_vawt(num_turbine_vawt),Zhub_vawt(num_turbine_vawt),
     +       wtomega_vawt(num_turbine_vawt),
     +       wty_vawt(num_turbine_vawt),Dbar_vawt(num_turbine_vawt),
     +       wtH_vawt(num_turbine_vawt),dir_vawt(num_turbine_vawt),	 
     +       H1bar_vawt(num_turbine_vawt),H2bar_vawt(num_turbine_vawt),
     +       chord_vawt(num_turbine_vawt),TSR_vawt(num_turbine_vawt)	 
 
      parameter (turbine_model_vawt=0,tower_model_vawt=0)

      integer*4 N_ang_vawt,N_blade_vawt
      parameter (N_ang_vawt=102,N_blade_vawt=3)

      real*8 CTB_vawt(num_turbine_vawt),CTT_vawt(num_turbine_vawt)
      real*8 ALFA_CL_CD_vawt(3,361)
	  
      common /turbine_pars1/ wtx_vawt,wtz_vawt,ALFA_CL_CD_vawt
      common /turbine_pars2/ wtR_vawt,Zhub_vawt,wtomega_vawt,CTB_vawt,
     +            CTT_vawt,wty_vawt,Dbar_vawt,wtH_vawt,
     +            H1bar_vawt,H2bar_vawt,dir_vawt,chord_vawt,TSR_vawt 	  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  parameter (Nxbs=191,Nxbm=192,Nxbe=192,Nxx=192)
	  parameter (vel_instant=0,buffer=0)
	  parameter (T_instant=0,bufferT=0)

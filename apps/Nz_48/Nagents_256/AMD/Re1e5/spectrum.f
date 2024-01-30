      subroutine spectrum(uu,kk,spec,t,flag)
      implicit none
      include 'dimen.h'
      include 'fftw3.f'

      integer*4 i,j,k,ii,iii,im,kk,t,flag,kko 
      integer*8 plan_fspec

      real*8,dimension(Nx,Ny,Nz2) :: uu
      real*8,dimension(Nx):: vel_c
      real*8,dimension(Nx/2+1):: C_2,sumP,Pgm
      real*8,dimension(Nx/2+1,nz2):: Spec
      double complex, dimension(Nx/2+1) :: vel_hat
      save plan_fspec
      
      if(t.eq.p_count.and.flag.eq.1.and.kk.eq.2)then
         call dfftw_plan_dft_r2c_1d(plan_fspec,Nx,vel_c,vel_hat,
     +        FFTW_PATIENT)
      endif

      sumP=0.d0
      
      do j=1,Ny

         vel_c=uu(:,j,kk)
         call dfftw_execute_dft_r2c(plan_fspec,vel_c,vel_hat)
         
         do im=1,Nx/2+1
            C_2(im)=dreal(vel_hat(im))**2+dimag(vel_hat(im))**2
         end do
         
         do ii=1,Nx/2+1
            Pgm(ii)=1./(Nx**2.)*(2.d0*C_2(ii))
            Pgm(1)=1./(Nx**2.)*C_2(1)
            Pgm(Nx/2+1)=1./(Nx**2.)*C_2(Nx/2+1)
            sumP(ii)=sumP(ii)+Pgm(ii)
         end do
         
      end do

      IF(nprocs.eq.1)THEN
         kko=kk-1
      else
         kko=kk
      endif
      do iii=1,Nx/2+1
         spec(iii,kko)=sumP(iii)/(1.d0*Ny)
      end do			

      if(t.eq.nsteps.and.flag.eq.2.and.kk.eq.nzb+1)then
         call dfftw_destroy_plan(plan_fspec)
      endif

      return
      end



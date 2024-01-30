      Subroutine Dealias1(u1,u1_m,t,flag,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)

      implicit none
!      include 'dimen.h'
      include 'fftw3.f'
      integer*4 :: nx,ny,nzb,nx2,ny2,nz2,nsteps
      real*8 :: inxny
      integer*4 i,j,k,t,flag
      integer*8 plan_f,plan_b
      real*8, dimension(nx,ny,nz2) :: u1 
      real*8, dimension(nx2,ny2,nz2) :: u1_m
      real*8, dimension(nx,ny) :: u1_2d
      real*8, dimension(nx2,ny2) :: u1_m_2d
      double complex, dimension(nx/2+1,ny) :: u1_hat
      double complex, dimension(nx2/2+1,ny2) :: u1_m_hat
!      save plan_f,plan_b

      if(t.eq.1.and.flag.eq.1)then
         call dfftw_plan_dft_r2c_2d(plan_f,Nx,Ny,u1_2D,u1_hat,
     +        FFTW_PATIENT)
         call dfftw_plan_dft_c2r_2d(plan_b,Nx2,Ny2,u1_m_hat,u1_m_2d,
     +        FFTW_PATIENT)
      endif
       

      do k=2,Nzb+1

c ... perform forward FFT

         u1_2d=u1(:,:,k)
         call dfftw_execute_dft_r2c(plan_f,u1_2D,u1_hat)


         u1_hat=u1_hat*inxny

c ... padding with 3/2 rule

         u1_m_hat(:,:)=dcmplx(0.d0,0.d0)

         do j=1,Ny
            do i=1,Nx/2+1
               if(i.lt.nint(nx/2.+1.).and.j.lt.nint(ny/2.+1.))then
                  u1_m_hat(i,j)=u1_hat(i,j)
               elseif(i.lt.nint(nx/2.+1.).and.j.gt.nint(ny/2.+1.))then
                  u1_m_hat(i,j+Ny/2)=u1_hat(i,j)
               endif
            enddo
         enddo
       
c ... Back to physical space
         

         call dfftw_execute_dft_c2r(plan_b,u1_m_hat,u1_m_2d)
         
         u1_m(:,:,k)=u1_m_2d
         
      enddo

      if(t.eq.nsteps.and.flag.eq.2)then
         call dfftw_destroy_plan(plan_f)
         call dfftw_destroy_plan(plan_b)
      endif

      end

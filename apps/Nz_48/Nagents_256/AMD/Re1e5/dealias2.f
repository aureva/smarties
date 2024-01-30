      Subroutine Dealias2 (uu,uu_m,t,flag,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)

      implicit none
      integer*8 plan_f,plan_b
!      include 'dimen.h'
      include 'fftw3.f'

      integer*4 :: nx,ny,nzb,nx2,ny2,nz2,nsteps
      real*8 :: inxny,inx2ny2
      integer*4 i,j,k,t,flag
      integer*8 plan_ff,plan_bb
      real*8, dimension(nx,ny) :: uu_2D
      real*8, dimension(nx2,ny2) :: uu_m_2D
      real*8, dimension(nx,ny,nz2) :: uu
      real*8, dimension(nx2,ny2,nz2) :: uu_m
      double complex, dimension(nx/2+1,ny) :: uu_hat
      double complex, dimension(nx2/2+1,ny2) :: uu_m_hat
!      save plan_ff,plan_bb

 
      if(t.eq.1.and.flag.eq.1)then
         call dfftw_plan_dft_r2c_2d(plan_ff,Nx2,Ny2,uu_m_2D,uu_m_hat,
     +        FFTW_PATIENT)
         call dfftw_plan_dft_c2r_2d(plan_bb,Nx,Ny,uu_hat,uu_2D,
     +        FFTW_PATIENT)
      endif

      do k=2,Nzb+1

c.......perform forward FFT
         
         uu_m_2D=uu_m(:,:,k) 
         call dfftw_execute_dft_r2c(plan_ff,uu_m_2D,uu_m_hat)

         uu_m_hat=uu_m_hat*inx2ny2

         uu_hat(:,:)=dcmplx(0.d0,0.d0)

         do j=1,Ny
            do i=1,Nx/2+1
               if(i.lt.(nx/2+1).and.j.lt.(ny/2+1))then
                  uu_hat(i,j)=uu_m_hat(i,j)
               elseif(i.lt.(nx/2+1).and.j.gt.(ny/2+1))then
                  uu_hat(i,j)=uu_m_hat(i,j+Ny/2)
               end if
            end do
         end do

c.....Back to physical space

         call dfftw_execute_dft_c2r(plan_bb,uu_hat,uu_2D)
         uu(:,:,k)=uu_2D

      end do

      if(t.eq.nsteps.and.flag.eq.2)then
         call dfftw_destroy_plan(plan_ff)
         call dfftw_destroy_plan(plan_bb)
      endif

      end

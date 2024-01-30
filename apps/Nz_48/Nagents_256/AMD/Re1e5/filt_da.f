      Subroutine FILT_DA (F,DFDX,DFDY,t,flag,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,fgr,l_r)

      implicit none
      include 'fftw3.f'
      !include 'dimen.h'

      integer*4 :: nx,ny,nzb,nx2,ny2,nz2,nsteps
      real*8 :: inxny,fgr,l_r
      integer*4 i,j,k,ii,jj,t,flag
      integer*8 plan_f,plan_b1,plan_b2,plan_b3
      real*8, dimension(Nx,Ny,Nz2) :: F,DFDX,DFDY
      real*8, dimension(Nx,Ny) :: F_2D,DFDX_2D,DFDY_2D 
      double complex, dimension(Nx/2+1,Ny) :: F_hat,DFDX_hat,DFDY_hat
      save plan_f,plan_b1,plan_b2,plan_b3

 
      if(t.eq.1.and.flag.eq.1)then
         call dfftw_plan_dft_r2c_2d(plan_f,Nx,Ny,F_2D,F_hat,
     +        FFTW_PATIENT)
         call dfftw_plan_dft_c2r_2d(plan_b1,Nx,Ny,F_hat,F_2D,
     +        FFTW_PATIENT)
         call dfftw_plan_dft_c2r_2d(plan_b2,Nx,Ny,DFDX_hat,DFDX_2D,
     +        FFTW_PATIENT)
         call dfftw_plan_dft_c2r_2d(plan_b3,Nx,Ny,DFDY_hat,DFDY_2D,
     +        FFTW_PATIENT)
      endif

      do k=2,Nzb+1

c.......perform forward FFT
         F_2D=F(:,:,k)
         call dfftw_execute_dft_r2c(plan_f,F_2D,F_hat)

         F_hat=F_hat*inxny      

         do j=1,Ny

            jj=j-1
            if(jj.gt.nint(Ny/2.)) jj=jj-Ny  
            jj=jj*l_r

            do i=1,Nx/2+1
        
               ii=i-1

               if (ii.ge.nint(Nx/(2.0*fgr)))then
                  F_hat(i,j)=dcmplx(0.d0,0.d0)
               elseif(abs(jj).ge.nint(l_r*Ny/(2.0*fgr)))then
                  F_hat(i,j)=dcmplx(0.d0,0.d0)
               end if
               
               DFDX_hat(i,j)=dcmplx(dimag(F_hat(i,j))*(-1.d0),
     +                 dreal(F_hat(i,j)))*ii
               DFDY_hat(i,j)=dcmplx(dimag(F_hat(i,j))*(-1.d0),
     +                 dreal(F_hat(i,j)))*jj
        
            end do
         end do	     
 
c...Inverse transform

         call dfftw_execute_dft_c2r(plan_b1,F_hat,F_2D)
         call dfftw_execute_dft_c2r(plan_b2,DFDX_hat,DFDX_2D)
         call dfftw_execute_dft_c2r(plan_b3,DFDY_hat,DFDY_2D)
         
         F(:,:,k)=F_2D
         DFDX(:,:,k)=DFDX_2D
         DFDY(:,:,k)=DFDY_2D

      end do       

      if(t.eq.nsteps.and.flag.eq.2)then
         call dfftw_destroy_plan(plan_f)
         call dfftw_destroy_plan(plan_b1)
         call dfftw_destroy_plan(plan_b2)
         call dfftw_destroy_plan(plan_b3)
      endif

      return
      end

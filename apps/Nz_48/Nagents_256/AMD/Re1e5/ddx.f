      Subroutine DDX (DFDX,F,t,flag)              
      implicit none      
      include 'dimen.h'
      include 'fftw3.f'

      integer*4 i,j,k,ii,jj,t,flag
      integer*8 plan_fx,plan_bx

      real*8, dimension(Nx,Ny,Nz2) :: F,DFDX
      real*8, dimension(Nx,Ny) :: F_2D,DFDX_2D
      double complex, dimension(Nx/2+1,Ny) :: F_hat,DFDX_hat
      save plan_fx,plan_bx

 
      if(t.eq.1.and.flag.eq.1)then
         call dfftw_plan_dft_r2c_2d(plan_fx,Nx,Ny,F_2D,F_hat,
     +        FFTW_PATIENT)
         call dfftw_plan_dft_c2r_2d(plan_bx,Nx,Ny,DFDX_hat,DFDX_2D,
     +        FFTW_PATIENT)
      endif
     
      do k=2,Nzb+1

c.......perform forward FFT
         F_2D=F(:,:,k)
         call dfftw_execute_dft_r2c(plan_fx,F_2D,F_hat)
         F_hat=F_hat*inxny          
   
         do j=1,Ny
          
            jj=j-1
            if(jj.gt.nint(Ny/2.)) jj=jj-Ny  
            
            do i=1,Nx/2+1
        
               ii=i-1
               
               if ((ii.eq.nint(Nx/2.)).or.(abs(jj).ge.nint(Ny/2.)))then
                  DFDX_hat(i,j)=dcmplx(0.d0,0.d0)
               else
                  DFDX_hat(i,j)=dcmplx(dimag(F_hat(i,j))*(-1.d0),
     +                 dreal(F_hat(i,j)))*ii
               end if
                
            end do
         end do	     

c...  Inverse transform to get pseudospec. derivs.

         call dfftw_execute_dft_c2r(plan_bx,DFDX_hat,DFDX_2D)
         DFDX(:,:,k)=DFDX_2D

      end do

      if(t.eq.nsteps.and.flag.eq.2)then
         call dfftw_destroy_plan(plan_fx)
         call dfftw_destroy_plan(plan_bx)
      endif

      return
      end

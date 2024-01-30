      Subroutine FILTER_La (F_hat,F,t,flag)              
      implicit none      
      include 'dimen.h'
      include 'fftw3.f'

      integer*4 i,j,k,ii,jj,t,flag
      integer*8 plan_f,plan_bhat

      real*8,dimension(nx,ny,nz2) :: F,F_hat
      real*8,dimension(nx,ny) :: F_2D
      double complex, dimension(nx/2+1,ny) :: F_hat_ft
      real*8 d
      save plan_f,plan_bhat

      if(t.eq.1.and.flag.eq.1)then
         call dfftw_plan_dft_r2c_2d(plan_f,Nx,Ny,F_2D,F_hat_ft,
     +        FFTW_PATIENT)
         call dfftw_plan_dft_c2r_2d(plan_bhat,Nx,Ny,F_hat_ft,F_2D,
     +        FFTW_PATIENT)
      endif

      do k=2,nzb+1

c.......perform forward FFT         

         F_2D=F(:,:,k)
         call dfftw_execute_dft_r2c(plan_f,F_2D,F_hat_ft)

         F_hat_ft=F_hat_ft*inxny

         do j=1,Ny

            jj=j-1
            if(jj.gt.nint(Ny/2.)) jj=jj-Ny
            jj=jj*l_r

            do i=1,Nx/2+1

               ii=i-1

c%%%%%%%%%%%%%%%%%% Circular filter %%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$               if(ii.ge.nint(nx/(2.0*fgr*tfr*tfr)).and.
c$$$     +              abs(jj).ge.nint(l_r*Ny/(2.0*fgr*tfr*tfr)))then
c$$$                  ii=1000
c$$$                  jj=1000
c$$$               end if
c$$$               d=+sqrt((1.*ii/Nx)**2.+(1.*jj/Ny)**2.)
c$$$               if (d.ge.(1.0/(2.0*fgr*tfr*tfr))) then 
c$$$                  F_hat_ftd(i,j)=F_hat_ftd(i,j)*0.0
c$$$               end if
c$$$         
c$$$               if(ii.ge.nint(nx/(2.0*fgr*tfr)).and.
c$$$     +              abs(jj).ge.nint(l_r*Ny/(2.0*fgr*tfr)))then
c$$$                  ii=1000
c$$$                  jj=1000
c$$$               end if
c$$$               d=+sqrt((1.*ii/Nx)**2.+(1.*jj/Ny)**2.)
c$$$               if (d.ge.(1.0/(2.0*fgr*tfr))) then
c$$$                  F_hat_ft(i,j)=F_hat_ft(i,j)*0.0
c$$$               end if
c%%%%%%%%%%%%%%%%%% Square filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%
               if (ii.ge.nint(Nx/(2.0*fgr*tfr)))then
                  F_hat_ft(i,j)=dcmplx(0.d0,0.d0)
               elseif(abs(jj).ge.nint(l_r*Ny/(2.0*fgr*tfr)))then
                  F_hat_ft(i,j)=dcmplx(0.d0,0.d0)
               end if
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end do
         end do

c.... Inverse transform to get pseudospec. derivs.        

         call dfftw_execute_dft_c2r(plan_bhat,F_hat_ft,F_2D)

         F_hat(:,:,k)=F_2D

      end do

      if(t.eq.nsteps.and.flag.eq.2)then
         call dfftw_destroy_plan(plan_f)
         call dfftw_destroy_plan(plan_bhat)
      endif

      return

      end

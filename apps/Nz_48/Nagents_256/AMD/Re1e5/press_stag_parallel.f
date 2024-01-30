      Subroutine Press_stag(P, RHSx, RHSy, RHSz, RHSx_f, RHSy_f,
     +     RHSz_f,u,v,w,DFDX,DFDY,divtz,me,nall,t)

      implicit none
      include 'dimen.h'
      include 'fftw3.f'

      INTEGER*4 NPSX,NPSY,NI,NJ,MeI,MeJ,I_GLOBAL,J_GLOBAL,Lsize
      PARAMETER(NPSX=1, NPSY=nps)
      PARAMETER(NI=(Nx/2+1)/NPSX, NJ=Ny/NPSY)

      integer*4 i,j,k,ii,jj,size,t
      integer*8 plan_Hx,plan_Hy,plan_Hz,plan_tz,plan_bP,plan_bdpdx,
     +     plan_bdpdy

      real*8,dimension(nx,ny,nz2):: RHSx,RHSy,RHSz,RHSx_f,RHSy_f,RHSz_f,
     +     u,v,w,dfdx,dfdy,divtz,P
      real*8,dimension(nx,ny):: oned_x,oned_y,oned_z,tz_2D,F_2D,
     +     DFDX_2D,DFDY_2D
      real*8,dimension(nz+1):: RHS_col,a1,b1,c1,p_colr,p_coli

      double complex, dimension(Nx/2+1,Ny) :: H_x_2D,H_y_2D,H_z_2D,
     +     tz_hat,bottomw,topw,F_hat,DFDX_hat,DFDY_hat
      double complex, dimension(Nx/2+1,Ny,Nz2) :: H_x,H_y,H_z
      double complex, dimension(NI,NJ,Nz) :: H_2,H_3

      save plan_Hx,plan_Hy,plan_Hz,plan_tz,
     +     plan_bP,plan_bdpdx,plan_bdpdy

      if(t.eq.1)then
         call dfftw_plan_dft_r2c_2d(plan_Hx,Nx,Ny,oned_x,H_x_2D,
     +        FFTW_PATIENT)
         call dfftw_plan_dft_r2c_2d(plan_Hy,Nx,Ny,oned_y,H_y_2D,
     +        FFTW_PATIENT)
         call dfftw_plan_dft_r2c_2d(plan_Hz,Nx,Ny,oned_z,H_z_2D,
     +        FFTW_PATIENT)
         call dfftw_plan_dft_r2c_2d(plan_tz,Nx,Ny,tz_2D,tz_hat,
     +        FFTW_PATIENT)
         call dfftw_plan_dft_c2r_2d(plan_bP,Nx,Ny,F_hat,F_2D,
     +        FFTW_PATIENT)
         call dfftw_plan_dft_c2r_2d(plan_bdpdx,Nx,Ny,DFDX_hat,DFDX_2D,
     +        FFTW_PATIENT)
         call dfftw_plan_dft_c2r_2d(plan_bdpdy,Nx,Ny,DFDY_hat,DFDY_2D,
     +        FFTW_PATIENT)
      endif

      k = Nzb+2

      size=Nz+1

      call update3(RHSz,RHSz_f,w,me,nall)
 

      do k=2,Nzb+2

         do j=1,Ny
            do i=1,Nx
               oned_x(i,j)=RHSx(i,j,k)-(1.0/3.0)*RHSx_f(i,j,k)+
     +              u(i,j,k)*(2.0/3.0)/DT
               oned_y(i,j)=RHSy(i,j,k)-(1.0/3.0)*RHSy_f(i,j,k)+
     +              v(i,j,k)*(2.0/3.0)/DT
               oned_z(i,j)=RHSz(i,j,k)-(1.0/3.0)*RHSz_f(i,j,k)+
     +              w(i,j,k)*(2.0/3.0)/DT
            enddo
         enddo

c.....Forward transform the planes of RHS matrices
         call dfftw_execute_dft_r2c(plan_Hx,oned_x,H_x_2D)
         call dfftw_execute_dft_r2c(plan_Hy,oned_y,H_y_2D)
         call dfftw_execute_dft_r2c(plan_Hz,oned_z,H_z_2D)
       

         H_x(:,:,k-1)=H_x_2D*inxny
         H_y(:,:,k-1)=H_y_2D*inxny
         H_z(:,:,k-1)=H_z_2D*inxny

      end do

      IF (me==0) then
c         tz_2D=divtz(:,:,2)
c         call dfftw_execute(plan_tz)
c         bottomw=tz_hat*inxny
         bottomw(:,:)=H_z(:,:,1)
      END IF

      IF (me==nprocs-1) then
c         tz_2D=divtz(:,:,Nzb+1)
c         call dfftw_execute(plan_tz)
c         topw=tz_hat*inxny
         topw=H_z(:,:,nzb)
      END IF

      CALL MPI_BCAST(bottomw(1,1),(Nx/2+1)*ny,MPI_DOUBLE_COMPLEX,
     +               0,NALL,IERR)

      CALL MPI_BCAST(topw(1,1),(Nx/2+1)*ny,MPI_DOUBLE_COMPLEX,
     +               nprocs-1,NALL,IERR)

      do j=1,Ny
       jj=j-1
       if(jj.ge.(Ny/2)) jj=jj-Ny
       jj=jj*l_r
       do i=1,Nx/2+1
        ii=i-1
        do k=1,Nzb
        H_x(i,j,k)=ii*H_x(i,j,k)*dcmplx(0,1)+jj*H_Y(I,J,K)*dcmplx(0,1)
     +             +(H_z(i,j,k+1)-H_z(i,j,k))/DZ
        enddo
       enddo
      enddo

      Lsize=NI*NJ*NZB
      do j=1,NPSY
         do i=1,NPSX
            k=(j-1)*NPSX+i-1
            if ( k.eq.me ) then
               MeI = i
               MeJ = j
            end if
            H_3(:,:,k*NZB+1:(k+1)*NZB)=
     +      H_x((i-1)*NI+1:i*NI,(j-1)*NJ+1:j*NJ,1:NZB)
         end do
      end do

      CALL MPI_ALLTOALL(H_3(1,1,1),Lsize,MPI_DOUBLE_COMPLEX,
     +                  H_2(1,1,1),Lsize,MPI_DOUBLE_COMPLEX,NALL,IERR)

      if(me<(NPSX*NPSY))then
         do j=1,NJ

            j_global=j + (MeJ-1)*NJ
            jj=j_global-1
            if(jj.gt.(Ny/2)) jj=jj-Ny
            jj=jj*l_r

            do i=1,NI

               i_global=i + (MeI-1)*NI
               ii=i_global-1

c...  Zero wavenumber problem is non-unique for matrix solution
c.... So just integrate up from wall with arbitrary P=const starting point.

c..   Real PART .............................................
c...  Near Wall Nodes
                  a1(1)=0.d0
                  b1(1)=-1.d0
                  c1(1)=1.d0
                  RHS_col(1)=1.d0*DZ*dreal(bottomw(i_global,j_global))
                  if ((ii.eq.0).and.(jj.eq.0)) then
                    a1(1)=0.d0
                    b1(1)=1.d0
                    c1(1)=0.d0
                    RHS_col(1)=0.d0
                  end if
c...  Interior nodes
                  do k=2,Nz
                     RHS_col(k)=dreal(H_2(i,j,k-1))
                     a1(k)=1.d0/(DZ**2)
                     b1(k)=(-ii*ii-jj*jj-2.d0/(DZ**2))
                     c1(k)=1.d0/(DZ**2)
                  end do
c..   top nodes
                  RHS_col(Nz+1)=1.d0*dreal(topw(i_global,j_global))*DZ
                  a1(nz+1)= -1.d0
                  b1(Nz+1)=  1.d0
                  c1(nz+1)=  0.d0

                  call tridag (a1,b1,c1,RHS_col,P_colr,size)

c..   Imag PART .............................................
c...  Near Wall Nodes
                  a1(1)=  0.d0
                  b1(1)= -1.d0
                  c1(1)=  1.d0
                  RHS_col(1)=1.d0*DZ*dimag(bottomw(i_global,j_global))
                  if ((ii.eq.0).and.(jj.eq.0)) then
                    a1(1)=0.d0
                    b1(1)=1.d0
                    c1(1)=0.d0
                    RHS_col(1)=0.d0
                  end if
c...  Interior Nodes
                  do k=2,Nz
                     RHS_col(k)=dimag(H_2(i,j,k-1))
                     a1(k)=1.d0/(DZ**2)
                     b1(k)=(-ii*ii-jj*jj-2.d0/(DZ**2))
                     c1(k)=1.d0/(DZ**2)
                  end do
c...  Top nodes
                  RHS_col(Nz+1)=1.d0*dimag(topw(i_global,j_global))*DZ
                  c1(Nz+1)=  0.d0
                  b1(Nz+1)=  1.d0
                  a1(Nz+1)= -1.d0

                  call tridag (a1,b1,c1,RHS_col,P_coli,Size)

c...  Put Pressure amplitudes in Matrix
                  do k=1,Nz
                     H_2(i,j,k)=dcmplx(p_colr(k+1),p_coli(k+1))
                  end do

            end do
         end do
      end if

      CALL MPI_ALLTOALL(H_2(1,1,1),Lsize,MPI_DOUBLE_COMPLEX,
     +                  H_3(1,1,1),Lsize,MPI_DOUBLE_COMPLEX,NALL,IERR)

      do j=1,NPSY
         do i=1,NPSX
            k=(j-1)*NPSX+i-1
            H_z((i-1)*NI+1:i*NI,(j-1)*NJ+1:j*NJ,1:NZB)=
     +      H_3(:,:,k*NZB+1:(k+1)*NZB)
         end do
      end do

cc... Cut the Nyquist
      H_z(Nx/2+1,:,:) = dcmplx(0.d0,0.d0)
      H_z(:,Ny/2+1,:) = dcmplx(0.d0,0.d0)

c...  Now need to get P_hat (i.e. H_z) to physical P(i,j,k)

      do k=2,Nzb+1

         F_hat=H_z(:,:,k-1)

         do j=1,Ny

            jj=j-1
            if (jj.gt.(Ny/2)) jj=jj-Ny
            jj=jj*l_r

            do i=1,Nx/2+1

               ii=i-1

               DFDX_hat(i,j)=dcmplx(dimag(F_hat(i,j))*(-1.d0),
     +              dreal(F_hat(i,j)))*ii
               DFDY_hat(i,j)=dcmplx(dimag(F_hat(i,j))*(-1.d0),
     +              dreal(F_hat(i,j)))*jj

            end do

         end do

cc...... Inverse Transforms


         call dfftw_execute_dft_c2r(plan_bP,F_hat,F_2D)
         call dfftw_execute_dft_c2r(plan_bdpdx,DFDX_hat,DFDX_2D)
         call dfftw_execute_dft_c2r(plan_bdpdy,DFDY_hat,DFDY_2D)

         P(:,:,k)=F_2D
         DFDX(:,:,k)=DFDX_2D
         DFDY(:,:,k)=DFDY_2D

      end do

!      if(t.eq.nsteps)then
!         call dfftw_destroy_plan(plan_Hx)
!         call dfftw_destroy_plan(plan_tz)
!         call dfftw_destroy_plan(plan_bP)
!         call dfftw_destroy_plan(plan_bdpdx)
!         call dfftw_destroy_plan(plan_bdpdy)
!      endif

      return
      end

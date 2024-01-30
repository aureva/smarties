      Subroutine Scalar_RHS(s,u_m,v_m,w_m,qx,qy,qz,txz,tyz,dsdx,dsdy,
     +     dsdz,RHS,me,Surf_flux,tt,nall)

      implicit none
      include 'dimen.h'
      integer*4 i,j,k,tt

      real*8, dimension(Nx,Ny,Nz2):: dsdx,dsdy,dsdz,dtemp,RHS,s,txz,tyz,
     +     temp,qx,qy,qz

      real*8, dimension(Nx2,Ny2,Nz2)::
     +     dsdx_m,dsdy_m,dsdz_m,u_m,v_m,w_m,cc,temp_m,s_m

      real*8 S_Surf(Nx,Ny),surf_flux(Nx,Ny),l(Nz2)

      call dealias1(dsdx,dsdx_m,tt,0)
      call dealias1(dsdy,dsdy_m,tt,0)
      if(Q_flag.eq.1)then
         call dealias1(dsdz,dsdz_m,tt,0)
      else   
         call dealias1(dsdz,dsdz_m,tt,2)
      endif
      call update1_m(dsdz_m,me,nall)
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

      Do k=2,Nzb+1
         Do j=1,Ny2
            Do i=1,Nx2 
               
               cc(i,j,k) = u_m(i,j,k  )*dsdx_m(i,j,k  )+
     +              v_m(i,j,k  )*dsdy_m(i,j,k  )+ 
     +              (w_m(i,j,k  )*dsdz_m(i,j,k  )+
     +              w_m(i,j,k+1)*dsdz_m(i,j,k+1))/2.
               
            end do
         end do
      end do

C     ... TOP LIP

      if (me.eq.nprocs-1) then
         do j=1,Ny2
            do i=1,Nx2

               cc(i,j,Nzb+1)= u_m(i,j,Nzb+1)*dsdx_m(i,j,Nzb+1)+
     +              v_m(i,j,Nzb+1)*dsdy_m(i,j,Nzb+1)+
     +              w_m(i,j,Nzb+1)*dsdz_m(i,j,Nzb+1)

            end do
         end do
      endif

      if(Q_flag.eq.1)then
         call dealias2(RHS,cc,tt,0)
      else
         call dealias2(RHS,cc,tt,2)
      endif

c     ... Now building the SGS part of the RHS.
c     ... Note: Since we bring the Conective term to RHS its sign changes.
c     ... Below "Temp" is used for SGS flux; its divergence is added to RHS

C     ... XXXXXXXXXXXXXXXXXXXX

      Do k=2,Nzb+1
         Do j=1,Ny
            Do i=1,Nx
               temp(i,j,k)=qx(i,j,k)
            end do
         end do
      end do
     
      Call DDX(dtemp,temp,tt,1)        

      Do k=2,Nzb+1
         Do j=1,Ny
            Do i=1,Nx
                 RHS(i,j,k) = -RHS(i,j,k)-dtemp(i,j,k)
              end do
           end do
      end do

C     ... YYYYYYYYYYYYYYYYYYYY

      Do k=2,Nzb+1
           Do j=1,Ny
              Do i=1,Nx
                 temp(i,j,k)=qy(i,j,k)
              end do
           end do
      end do
        
      Call DDY (dtemp,temp,tt,1)     

c...  And Smagorinsky in interior
c.....Note dsdz on W nodes...	   

      Do k=2,Nzb+1
           Do j=1,Ny
              Do i=1,Nx
                 RHS(i,j,k)=RHS(i,j,k)-dtemp(i,j,k)
              end do
           end do
      end do

C     ... ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
      Do k=2,Nzb+2
           Do j=1,Ny
              Do i=1,Nx
                 temp(i,j,k)=qz(i,j,k)
              end do
           end do
      end do

c...  Use MO flux at wall

        if (me.eq.0) then
           Do j=1,Ny
                Do i=1,Nx
                   temp(i,j,2)= surf_flux(i,j)
                   qz(i,j,2)  = surf_flux(i,j)
                end do
           end do
        endif

c     ... The SGS_z flux is on the W nodes, 
c     ... but DDZ_W will put it back on UVP nodes

          call DDZ_w (dtemp, temp, me)

          Do k=2,Nzb+1
             Do j=1,Ny
                Do i=1,Nx
                   
                   RHS(i,j,k) = RHS(i,j,k)-dtemp(i,j,k)

                end do
             end do
          end do

          return

          end

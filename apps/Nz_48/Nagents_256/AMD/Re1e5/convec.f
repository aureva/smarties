      Subroutine CONVEC (Cx,Cy,Cz,u1_m,u2_m,u3_m,
     +     du1,du2,du3,du4,du5,du6,tt,me,nall,
     + nx,ny,nx2,ny2,nz2,nzb,nprocs,inxny,nsteps,inx2ny2,
     + plan_f,plan_b,
     + plan_ff,plan_bb)
      implicit none
!      include 'dimen.h'
      integer*4 i,j,k,tt
      Real*8,dimension (nx,ny,nz2):: du1,du2,du3,du4,du5,du6,
     +     cx,cy,cz
      real*8,dimension (nx2,ny2,nz2):: cc,u1_m,u2_m,u3_m,du1_m,du2_m,
     +     du3_m,du4_m,du5_m,du6_m
      real*8 arg1, arg2a, arg2b, arg2
      integer*4 :: nx,ny,nzb,nx2,ny2,nz2,nprocs,me,nall,nsteps
      real*8 :: inxny,inx2ny2 
      integer*8 :: plan_f,plan_b,
     + plan_ff,plan_bb   


      call dealias1(du1,du1_m,tt,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(du2,du2_m,tt,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(du3,du3_m,tt,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(du4,du4_m,tt,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(du5,du5_m,tt,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)
      call dealias1(du6,du6_m,tt,0 ,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,plan_f,plan_b)

      call update3_m(du2_m,du4_m,du5_m,me,nall)
      call update1_m(du6_m,me,nall)
      
c     ************************************************************
      Do k=2,nzb+1
         Do j=1,ny2
            Do i=1,nx2

               arg1=u2_m(i,j,k)*(du1_m(i,j,k)-du3_m(i,j,k))
               arg2a=u3_m(i,j,k+1)*(du2_m(i,j,k+1)-du5_m(i,j,k+1))
               arg2b=u3_m(i,j,k)*(du2_m(i,j,k)-du5_m(i,j,k)) 
               cc(i,j,k)=arg1+0.5d0*(arg2a+arg2b)
            End Do
         End Do
      End Do
      
      Do j=1,ny2
         Do i=1,nx2                 
            
            IF (me==0) then            
               k=2

               arg1=u2_m(i,j,k)*(du1_m(i,j,k)-du3_m(i,j,k))
               arg2a=u3_m(i,j,k+1)*(du2_m(i,j,k+1)-du5_m(i,j,k+1))
               arg2b=0d0 
               cc(i,j,k)=arg1+0.5d0*(arg2a+arg2b) 
            END IF
            
            IF (me==nprocs-1) then                    
               k=Nzb+1
               arg1=u2_m(i,j,k)*(du1_m(i,j,k)-du3_m(i,j,k))
               arg2a=u3_m(i,j,k)*(du2_m(i,j,k)-du5_m(i,j,k)) 
               cc(i,j,k)=arg1+arg2a                  
            END IF
            
         End Do
      End Do
     
      call dealias2(Cx,cc,tt,1,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)

c     **********************************************************
      
      Do k=2,nzb+1
         Do j=1,ny2
            Do i=1,nx2

               arg1=u1_m(i,j,k)*(du3_m(i,j,k)-du1_m(i,j,k))
               arg2a=u3_m(i,j,k+1)*(du4_m(i,j,k+1)-du6_m(i,j,k+1))
               arg2b=u3_m(i,j,k)*(du4_m(i,j,k)-du6_m(i,j,k)) 
               cc(i,j,k)=arg1+0.5d0*(arg2a+arg2b)                  
            End Do
         End Do
      End Do
      
      Do j=1,ny2
         Do i=1,nx2
            
            if (me==0) then
               k=2            
 
               arg1=u1_m(i,j,k)*(du3_m(i,j,k)-du1_m(i,j,k))
               arg2a=u3_m(i,j,k+1)*(du4_m(i,j,k+1)-du6_m(i,j,k+1))
               arg2b=0d0
               cc(i,j,k)=arg1+0.5d0*(arg2a+arg2b) 
            end if
            
            if (me==nprocs-1) then
               k=Nzb+1
               arg1=u1_m(i,j,k)*(du3_m(i,j,k)-du1_m(i,j,k))
               arg2a=u3_m(i,j,k)*(du4_m(i,j,k)-du6_m(i,j,k)) 
               cc(i,j,k)=arg1+arg2a                  
            end if
            
         End Do
      End Do
      
      call dealias2(Cy,cc,tt,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)

c     ******************************************************* 

      Do k=2,nzb+1
         Do j=1,ny2
            Do i=1,nx2
               arg1=(0.5d0*(u1_m(i,j,k)+u1_m(i,j,k-1)))*
     +              (du5_m(i,j,k)-du2_m(i,j,k))
               arg2=(0.5d0*(u2_m(i,j,k)+u2_m(i,j,k-1)))*
     +              (du6_m(i,j,k)-du4_m(i,j,k))
               cc(i,j,k)=arg1+arg2                  
            End Do
         End Do
      End Do
      
      do j=1,ny2
         do i=1,nx2
            
            if(me==0) then
               k=2 
               cc(i,j,k)=0d0
            end if
            
            if(me==nprocs-1)then
               k=nzb+1
               cc(i,j,k)=0d0
            end if
            
         end do
      end do
      
      call dealias2(Cz,cc,tt,0,
     + nx,ny,nz2,nx2,ny2,nzb,inxny,nsteps,inx2ny2,
     + plan_ff,plan_bb)
      
      Return
      End            
      

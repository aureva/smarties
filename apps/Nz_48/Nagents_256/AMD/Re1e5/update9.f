C     mpi update ghostlayers for 9 variables

C     ... TRICKS for MPI: SEND   .... TO   PROCESSOR
C     ...                 RECEIVE ... FROM PROCESSOR
C     ... For example, the following indicates the information
C     ... from the processor me==0 is sent to processor me+1!

      subroutine update9(dudx,dudy,dudz,dvdx,dvdy,dvdz,
     +     dwdx,dwdy,dwdz,me,nall)
      implicit none
      include 'dimen.h'
      real*8, dimension(nx,ny,nz2):: dudx,dudy,dudz,dvdx,dvdy,dvdz,
     +     dwdx,dwdy,dwdz

      IF(nprocs.gt.1)THEN

      IF (me==0) then

c     nzb+1 upper layer
      
         call MPI_SEND(dudx(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )
         call MPI_SEND(dudy(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )
         call MPI_SEND(dudz(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )   
         call MPI_SEND(dvdx(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )
         call MPI_SEND(dvdy(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )
         call MPI_SEND(dvdz(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )   
         call MPI_SEND(dwdx(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )
         call MPI_SEND(dwdy(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )
         call MPI_SEND(dwdz(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )

         call MPI_RECV(dudx(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )
         call MPI_RECV(dudy(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )
         call MPI_RECV(dudz(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )   
         call MPI_RECV(dvdx(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )
         call MPI_RECV(dvdy(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )
         call MPI_RECV(dvdz(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )   
         call MPI_RECV(dwdx(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )
         call MPI_RECV(dwdy(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )
         call MPI_RECV(dwdz(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )  
       
      ELSEIF(me<=nprocs-2) then
        
c     mpi RECV lower ghostlayer	

         call MPI_RECV(dudx(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )
         call MPI_RECV(dudy(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )
         call MPI_RECV(dudz(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )   
         call MPI_RECV(dvdx(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )
         call MPI_RECV(dvdy(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )
         call MPI_RECV(dvdz(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )   
         call MPI_RECV(dwdx(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )
         call MPI_RECV(dwdy(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )
         call MPI_RECV(dwdz(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr ) 

c     mpi send lower ghostlayer	

         call MPI_SEND(dudx(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )
         call MPI_SEND(dudy(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )
         call MPI_SEND(dudz(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )   
         call MPI_SEND(dvdx(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )
         call MPI_SEND(dvdy(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )
         call MPI_SEND(dvdz(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )   
         call MPI_SEND(dwdx(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )
         call MPI_SEND(dwdy(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )
         call MPI_SEND(dwdz(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )

c     mpi send upper layer

         call MPI_SEND(dudx(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )
         call MPI_SEND(dudy(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )
         call MPI_SEND(dudz(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )   
         call MPI_SEND(dvdx(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )
         call MPI_SEND(dvdy(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )
         call MPI_SEND(dvdz(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )   
         call MPI_SEND(dwdx(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )
         call MPI_SEND(dwdy(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )
         call MPI_SEND(dwdz(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr ) 

c     mpi RECV lower ghostlayer 	

         call MPI_RECV(dudx(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )
         call MPI_RECV(dudy(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )
         call MPI_RECV(dudz(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )   
         call MPI_RECV(dvdx(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )
         call MPI_RECV(dvdy(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )
         call MPI_RECV(dvdz(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )   
         call MPI_RECV(dwdx(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )
         call MPI_RECV(dwdy(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )
         call MPI_RECV(dwdz(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )

      ELSE
        
c     mpi	upper node

         call MPI_RECV(dudx(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )
         call MPI_RECV(dudy(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )
         call MPI_RECV(dudz(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )   
         call MPI_RECV(dvdx(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )
         call MPI_RECV(dvdy(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )
         call MPI_RECV(dvdz(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )   
         call MPI_RECV(dwdx(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )
         call MPI_RECV(dwdy(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )
         call MPI_RECV(dwdz(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr ) 
      
         call MPI_SEND(dudx(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )
         call MPI_SEND(dudy(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )
         call MPI_SEND(dudz(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )   
         call MPI_SEND(dvdx(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )
         call MPI_SEND(dvdy(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )
         call MPI_SEND(dvdz(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )
         call MPI_SEND(dwdx(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )
         call MPI_SEND(dwdy(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )
         call MPI_SEND(dwdz(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr ) 
     
      ENDIF

      ENDIF

      return
      
      end

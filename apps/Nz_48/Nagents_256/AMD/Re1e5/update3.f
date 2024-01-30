c     ... MPI UPDATE GHOSTLAYERS for any three variables

C     ... TRICKS for MPI: SEND   .... TO   PROCESSOR
C     ...                 RECEIVE ... FROM PROCESSOR
C     ... For example, the following indicates the information
C     ... from the processor me==0 is sent to processor me+1!

      subroutine update3(u,v,w,me,nall)
      implicit none
      include 'dimen.h'
      real*8, dimension(nx,ny,nz2):: u,v,w

      IF(nprocs.gt.1)THEN

      IF (me==0) then
         call MPI_SEND(u(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )   
         call MPI_SEND(v(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )   
         call MPI_SEND(w(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr ) 

         call MPI_RECV(u(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )   
         call MPI_RECV(v(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )   
         call MPI_RECV(w(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr ) 

      ELSEIF (me<=nprocs-2) then

c     mpi RECV lower ghostlayer	

         call MPI_RECV(u(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )   
         call MPI_RECV(v(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )   
         call MPI_RECV(w(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr ) 	

c     mpi send lower layer	

         call MPI_SEND(u(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )   
         call MPI_SEND(v(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )   
         call MPI_SEND(w(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr ) 

c     mpi send upper layer

         call MPI_SEND(u(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )   
         call MPI_SEND(v(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )   
         call MPI_SEND(w(1,1,nzb+1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr ) 

c     mpi recieve upper layers

         call MPI_RECV(u(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )   
         call MPI_RECV(v(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )   
         call MPI_RECV(w(1,1,nzb+2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr ) 

      ELSE

c     mpi	upper node

         call MPI_RECV(u(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )   
         call MPI_RECV(v(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )   
         call MPI_RECV(w(1,1,1),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr ) 

         call MPI_SEND(u(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )   
         call MPI_SEND(v(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )   
         call MPI_SEND(w(1,1,2),nx*ny, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )

      ENDIF

      ENDIF

      return
      
      end

cccc update for nx2 by ny2 cccc

      subroutine update3_m(u,v,w,me,nall)
      implicit none
      include 'dimen.h'
      real*8, dimension(nx2,ny2,nz2):: u,v,w
      real*8, dimension(nx2,ny2,nz2):: u_tmp,v_tmp,w_tmp
      integer*4 :: i,j,k
  
      u_tmp(:,:,:) = u(:,:,:)
      v_tmp(:,:,:) = v(:,:,:)
      w_tmp(:,:,:) = w(:,:,:)
      

      IF(nprocs.gt.1)THEN

      IF (me==0) then
         call MPI_SEND(u(1,1,nzb+1),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )   
         call MPI_SEND(v(1,1,nzb+1),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )   
         call MPI_SEND(w(1,1,nzb+1),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr ) 

         call MPI_RECV(u(1,1,nzb+2),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )   
         call MPI_RECV(v(1,1,nzb+2),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )   
         call MPI_RECV(w(1,1,nzb+2),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr ) 

      ELSEIF (me<=nprocs-2) then

c     mpi RECV lower ghostlayer	


         call MPI_RECV(u(1,1,1),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )   
         call MPI_RECV(v(1,1,1),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )   
         call MPI_RECV(w(1,1,1),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr ) 	

c     mpi send lower layer	


         call MPI_SEND(u(1,1,2),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )   
         call MPI_SEND(v(1,1,2),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )   
         call MPI_SEND(w(1,1,2),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr ) 


c     mpi send upper layer

         call MPI_SEND(u(1,1,nzb+1),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )   
         call MPI_SEND(v(1,1,nzb+1),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr )   
         call MPI_SEND(w(1,1,nzb+1),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me+1,me+1,nall,ierr ) 

c     mpi recieve upper layers

         call MPI_RECV(u(1,1,nzb+2),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )   
         call MPI_RECV(v(1,1,nzb+2),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr )   
         call MPI_RECV(w(1,1,nzb+2),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me+1,me,nall,status2,ierr ) 

      ELSE

c     mpi	upper node

         call MPI_RECV(u(1,1,1),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )   
         call MPI_RECV(v(1,1,1),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr )   
         call MPI_RECV(w(1,1,1),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me-1,me,nall,status2,ierr ) 

         call MPI_SEND(u(1,1,2),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )   
         call MPI_SEND(v(1,1,2),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )   
         call MPI_SEND(w(1,1,2),nx2*ny2, MPI_DOUBLE_PRECISION,
     +        me-1,me-1,nall,ierr )

      ENDIF
      
      ENDIF

      return
      
      end

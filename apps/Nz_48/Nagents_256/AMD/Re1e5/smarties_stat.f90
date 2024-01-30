!----------------------------------------------!
!     Module for interfacing with smarties     !
!----------------------------------------------!
Include 'smarties.f90'

Module smarties_stat

  ! Modules
  Use smarties
  ! prevent implicit typing
  Implicit None

Contains

  !----------------------------------------------!
  !                                              !
  !         communicate with smarties            !
  !                                              !
  !----------------------------------------------!
  Subroutine send_recv_state_action(u,v,smarties_comm,istep,num_actions,&
    state_size,nx,ny,nz2,dx,dy,dz,vonk,nu,Nsteps,nx_agents,ny_agents,&
    tauw,rewards_old,zo,M,z_i,me,neighbor,ct,c_count,RL_init,nall,&
    ierr,mpi_double_precision,RL_stop,nprocs,status2,mpi_status_size,train,&
    actions_2D,states_2D) &
    bind(c, name='send_recv_state_action')
    Use, intrinsic :: iso_c_binding

    Use smarties

      implicit none

      !include 'dimen.h'

      integer, parameter :: sp = kind(1.0)
      integer, parameter :: dp = kind(1.0d0)
      character(len=10) :: file_id
      real(dp) :: vonk,nu,dz,dx,dy,z_i
      integer(kind=4) :: nx,ny,nz2, nsteps,me,ct,c_count
      real(dp), dimension(nx,ny,nz2):: u,v
      type(c_ptr),    intent(in), value :: smarties_comm
      integer(kind=4) :: i, j, k, ii, jj, i1,i2,j1,j2
      integer(kind=4) :: nx_agents,ny_agents,RL_init
      Real(dp)    :: Uloc, dUloc,f11,f12,f11a,f12a
      Real(dp),dimension(2)    :: f11s,f12s
      Real(dp)    :: dPdx_ref, ustar, denom
      integer(c_int),intent(in),value    :: istep
      integer(kind=4) :: state_size, num_actions,nall
      real(c_double), dimension(num_actions), target ::&
                          upper_action_bound,lower_action_bound
      real   (c_double), dimension(num_actions), target :: action
      real   (c_double), dimension(state_size),  target :: state
      real   (c_double), dimension(nx_agents,ny_agents,2),target :: states
      real   (c_double), dimension(nx_agents,ny_agents),target :: actions
      real   (c_double), dimension(nx_agents,ny_agents),target :: rewards
      real*8, dimension(nx_agents,ny_agents) :: utau
      real*8, dimension(ny_agents) :: rewards_old
      real*8, dimension(nx_agents,ny_agents) :: actions_2D
      real*8, dimension(nx_agents,ny_agents,2) :: states_2D
      real   (c_double), dimension(nx,ny),target :: tauw
!      real*8, dimension(nx,ny,2),target :: states_2D
!      real*8, dimension(nx,ny),target :: actions_2D
      real(dp),dimension(nx,ny) :: zo,M
      real   (c_double) :: reward
      integer(kind=4) :: agent_id,ierr,mpi_double_precision,RL_stop
      integer*4, dimension(nx,ny,4) :: neighbor
      real :: rand
      integer*4 :: jmin,jmax,nprocs,proc
      real*8, dimension(ny_agents) :: buffer_action,buffer_reward
      real*8, dimension(2*ny_agents) :: buffer_state
      integer*4 :: mpi_status_size,train
      integer*4 :: status2(mpi_status_size)

      dPdx_ref = 1d0
      reward = 0d0

      Do i = 1,nx_agents
      Do j = 1,ny_agents

        ii = 1 + nx/nx_agents*(i-1)
        jj = 1 + ny/ny_agents*(j-1)


        denom=dlog((dz/2d0)/(zo(i,j)/z_i))
        if (istep .eq. 1) then
          call random_number(rand)
          if (train == 0) then
            rand = 0.5d0
          end if
          utau(i,j) = (dsqrt(u(ii,jj,2)**2+v(ii,jj,2)**2)*vonk/denom)
          tauw(ii,jj) = utau(i,j)**2*(1+(rand-0.5d0)/5d0)
        else
          utau(i,j) = dsqrt(dabs(tauw(ii,jj)))
        end if

        Uloc = 0.5d0*(u(ii,jj,2)+u(ii,jj,3))
        dUloc = (u(ii,jj,3) -&
                 u(ii,jj,2))/dz

        states(i,j,1) = (dUloc*dz/utau(i,j)-1d0/vonk)*dlog(dz*utau(i,j)/nu)
        states(i,j,2) = Uloc/utau(i,j) -&
                    dUloc*dz/utau(i,j)*dlog(dz*utau(i,j)/nu)

        rewards(i,j) = - abs(tauw(ii,jj)-dPdx_ref) / dPdx_ref

        state(1) = states(i,j,1)
        state(2) = states(i,j,2)
        
        states_2D(i,j,1) = states(i,j,1)
        states_2D(i,j,2) = states(i,j,2)
        DO k = 3,6
          state(k) = 0d0
        end do
        agent_id = (i-1) + nx_agents*(j-1)
        if (istep ==1) then
          call smarties_sendInitState(smarties_comm, c_loc(state) ,&
              state_size,agent_id)
          call smarties_recvAction(smarties_comm, c_loc(action),&
              num_actions, agent_id)
          
        else
          reward = 0d0
          call smarties_sendState (smarties_comm, c_loc(state) ,&
               state_size, reward, agent_id)
          call smarties_recvAction(smarties_comm, c_loc(action),&
               num_actions, agent_id)

        end if
        tauw(ii,jj) = tauw(ii,jj)*action(1)
        actions_2D(i,j) = action(1)

      end do
      End Do
      
      if (me ==0) then  
      DO i=1,nx
        tauw(i,ny) = tauw(i,1)
!        states_2D(i,ny,:) = states_2D(i,1,:)
!        actions_2D(i,ny) = actions_2D(i,1)
      end do
      DO j = 1,ny
        tauw(nx,j) = tauw(1,j)
!        states_2D(nx,j,:) = states_2D(1,j,:)
!        actions_2D(nx,j) = actions_2D(1,j)
      end do
 

      If (istep .le. nsteps) Then
        DO i = 1,nx
          DO j = 1,ny
            i1 = neighbor(i,j,1)
            i2 = neighbor(i,j,2)

            j1 = neighbor(i,j,3)
            j2 = neighbor(i,j,4)
 
!           Bilinear interpolation
            if (i2.gt.nx) then
              
              if (j2.gt.ny) then
                f11 = tauw(i1,j1) 
!                f11s(:) = states_2D(i1,j1,:) 
!                f11a = actions_2D(i1,j1) 

              else
                f11 = tauw(i1,j1)
                f12 = tauw(i1,j2)
!                f11s(:) = states_2D(i1,j1,:) 
!                f12s(:) = states_2D(i1,j2,:) 
!                f11a = actions_2D(i1,j1)
!                f12a = actions_2D(i1,j2)
              end if

            else 

              if (j2.gt.ny) then
                f11 = float(i2-i )/float(i2-i1)*tauw(i1,j1) &
                    + float(i -i1)/float(i2-i1)*tauw(i2,j1)
!                f11s(:) = float(i2-i )/float(i2-i1)*states_2D(i1,j1,:) &
!                        + float(i -i1)/float(i2-i1)*states_2D(i2,j1,:)
!                f11a = float(i2-i )/float(i2-i1)*actions_2D(i1,j1) &
!                     + float(i -i1)/float(i2-i1)*actions_2D(i2,j1)
              else
                f11 = float(i2-i )/float(i2-i1)*tauw(i1,j1) &
                    + float(i -i1)/float(i2-i1)*tauw(i2,j1)
                f12 = float(i2-i )/float(i2-i1)*tauw(i1,j2) &
                    + float(i -i1)/float(i2-i1)*tauw(i2,j2)
!                f11s(:) = float(i2-i )/float(i2-i1)*states_2D(i1,j1,:) &
!                        + float(i -i1)/float(i2-i1)*states_2D(i2,j1,:)
!                f12s(:) = float(i2-i )/float(i2-i1)*states_2D(i1,j2,:) &
!                        + float(i -i1)/float(i2-i1)*states_2D(i2,j2,:)
!                f11a = float(i2-i )/float(i2-i1)*actions_2D(i1,j1) &
!                     + float(i -i1)/float(i2-i1)*actions_2D(i2,j1)
!                f12a = float(i2-i )/float(i2-i1)*actions_2D(i1,j2) &
!                     + float(i -i1)/float(i2-i1)*actions_2D(i2,j2)
              end if

            end if

            if (j2.gt.ny) then
              tauw(i,j) = f11 
!              states_2D(i,j,:) = f11s(:) 
!              actions_2D(i,j) = f11a 
            else
              tauw(i,j) = float(j2-j )/float(j2-j1)*f11 &
                        + float(j -j1)/float(j2-j1)*f12
!              states_2D(i,j,:) = float(j2-j )/float(j2-j1)*f11s(:) &
!                               + float(j -j1)/float(j2-j1)*f12s(:)
!              actions_2D(i,j) = float(j2-j )/float(j2-j1)*f11a &
!                              + float(j -j1)/float(j2-j1)*f12a
            end if

          end do
        enddo
      end if
      end if
     
        


!      if (istep .lt. 5000 .AND. mod(istep,100)==0) then  
!      if (me==0) then
!      write(file_id, '(i0)') istep
!      open(unit=55,file = trim( 'output/tauw'// trim(adjustl(file_id)) &
!                // '.csv'),status="new",position="append")
!      do i=1,nx
!        do j=1,ny
!          write(55,*) i,',',j,',',tauw(i,j),',',states_2D(i,j,1)&
!                      ,',',states_2D(i,j,2),',',actions_2D(i,j)
!        end do
!      enddo
!      close(55)       
!      endif 
!      endif


 
      return



  End Subroutine send_recv_state_action
End Module smarties_stat

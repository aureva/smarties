      subroutine output_actions(states_2D,actions_2D,
     +                          t,nx_agents,ny_agents)

      implicit none
      integer*4 i,j,t,nx_agents,ny_agents
      
      real*8,dimension(nx_agents,ny_agents,2) :: states_2D
      real*8,dimension(nx_agents,ny_agents) :: actions_2D
      character(len=10) :: file_id
      
      write(file_id, '(i0)') t
      open(unit=55,file = trim( 'output/actions_map'// 
     +   trim(adjustl(file_id))
     +   // '.csv'),status="new",position="append")
      do i=1,nx_agents
        do j=1,ny_agents
          write(55,*) i,',',j,',',states_2D(i,j,1),',',
     +                 states_2D(i,j,2),',',actions_2D(i,j)
        end do
      enddo
      close(55)
 
      return    
      
      end



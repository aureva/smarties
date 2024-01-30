      subroutine output_tauw(tauw,t,nx,ny)

      implicit none
      integer*4 i,j,t,nx,ny
      
      real*8,dimension(nx,ny) :: tauw
      character(len=10) :: file_id
      
      write(file_id, '(i0)') t
      open(unit=55,file = trim( 'output/tauw'// trim(adjustl(file_id))
     +   // '.csv'),status="new",position="append")
      do i=1,nx
        do j=1,ny
          write(55,*) i,',',j,',',tauw(i,j)
        end do
      enddo
      close(55)
 
      return    
      
      end



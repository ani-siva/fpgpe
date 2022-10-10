module output
  contains
    subroutine write_den(fp,den)
      use cdata,only : nx,ny,x,y
      implicit none
      real(8),dimension(0:nx,0:ny),intent(in) :: den
      integer,intent(in) :: fp
      integer :: i,j

      do i = 0,nx ; do j = 0,ny
        write(fp,20) x(i),y(j),den(i,j)
      end do
        write(fp,*)
      end do
      20 format((2f12.6,g17.8e3))
    end subroutine write_den

end module output 

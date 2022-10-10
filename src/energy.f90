module energy
  contains
    subroutine norm(cp,zn)
      use cdata,only : dx,dy,nx,ny
      use omp_lib
      implicit none
      complex(8),dimension(0:nx,0:ny),intent(inout) :: cp
      real(8),intent(out) :: zn
      integer :: i,j
      real(8),dimension(0:nx,0:ny) :: tmp2d
      !$omp parallel do private(i,j)
      do i = 0,nx ; do j = 0,ny
        tmp2d(i,j) = cp(i,j)%re * cp(i,j)%re + cp(i,j)%im * cp(i,j)%im
      end do ; end do
      !$omp end parallel do

      zn = sqrt(integrate(tmp2d,dx,dy))

      !$omp parallel do private(i,j)
      do i = 0,nx ; do j = 0,ny
        cp(i,j) = cp(i,j)/zn
      end do ; end do
      !$omp end parallel do
    end subroutine norm

    function integrate(u2,dx,dy) result(res)
      use cdata,only : nx,ny
      implicit none
      real(8),dimension(0:,0:),intent(in) :: u2
      real(8),intent(in) :: dx,dy
      real(8) :: res
      real(8),dimension(:),allocatable :: tmp1d
      integer :: i
      allocate(tmp1d(0:nx))
      !$omp parallel do private(i)
      do i = 0,nx
        tmp1d(i) = simp(u2(i,0:),dy)
      end do
      !$omp end parallel do
      res = simp(tmp1d,dx)
      deallocate(tmp1d)
    end function integrate
    
    pure function simp(f,dx)
      implicit none
      real(8),dimension(0:),intent(in) :: f
      real(8),intent(in) :: dx
      real(8) :: simp
      real(8) :: f1,f2
      integer :: i,n
      n = size(f)-1
      f1 = f(1) + f(n-1)
      f2 = f(2)
      do i = 3,n-3,2
        f1 = f1 + f(i)
        f2 = f2 + f(i+1)
      end do
      simp = dx * (f(0) + 4.0d0 * f1 + 2.0d0 * f2 + f(n))/3.0d0
    end function simp

end module energy

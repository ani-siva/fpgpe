module derive
  use omp_lib
  contains
    subroutine calc_du(cp,dcp)
      use cdata,only : nx,ny,ci,g2d,v,lu,set_lhy,setup_dipole,ldg
      use initial,only : set_dipole
      implicit none

      complex(8),dimension(0:nx,0:ny),intent(inout) :: cp
      complex(8),dimension(0:nx,0:ny),intent(out) :: dcp
      real(8),dimension(:,:),allocatable :: dicp
      integer :: i,j
      real(8) :: ptmp,lhytmp
      call calc_lap(cp)
      !$omp parallel do private(i,j,ptmp)
      do i = 0,nx ; do j = 0,ny
        ptmp = cp(i,j)%re * cp(i,j)%re + cp(i,j)%im * cp(i,j)%im
        dcp(i,j) = -ci * ( -0.5 * lu(i,j) + v(i,j) * cp(i,j) + g2d * cp(i,j) * ptmp)
      end do ;  enddo
      !$omp end parallel do



      !---Modifications---!

      !--Set LHY correction--!
      if(set_lhy) then 
      if(setup_dipole) then
      !$omp parallel do private(i,j,ptmp,lhytmp)
        do i = 0,nx ; do j = 0,ny
        lhytmp = cp(i,j)%re * cp(i,j)%re * cp(i,j)%re + cp(i,j)%im * cp(i,j)%im * cp(i,j)%im
        dcp(i,j) = dcp(i,j) -ci * (-cp(i,j) * ldg * lhytmp)
      end do ;  end do
      !$omp end parallel do
      else
      !$omp parallel do private(i,j,ptmp,lhytmp)
        do i = 0,nx ; do j = 0,ny
        lhytmp = cp(i,j)%re * cp(i,j)%re * cp(i,j)%re + cp(i,j)%im * cp(i,j)%im * cp(i,j)%im
        dcp(i,j) = dcp(i,j) -ci * (-cp(i,j) * lhytmp)
      end do ;  end do
      !$omp end parallel do
      end if
      end if

      !--Dipole----!
      if(setup_dipole) then
        allocate(dicp(0:nx,0:ny))
        call set_dipole(cp,dicp)
        !$omp parallel do private(i,j)
        do i = 0,nx ; do j = 0,ny
        dcp(i,j) = dcp(i,j) - ci * (cp(i,j) * dicp(i,j))
        end do ; end do
        !$omp end parallel do
        deallocate(dicp)
      end if
      return
    end subroutine calc_du 

    subroutine calc_lap(cp)
      use cdata,only : nxx,nyy,nx,ny,lu,k2,fftw_execute_dft,planf,planb,s1,s2
      implicit none
      complex(8),dimension(0:nx,0:ny),intent(in) :: cp
      s1 = cp(0:nxx,0:nyy)
      call fftw_execute_dft(planf,s1,s2)
      s1 = -k2 * s2
      call fftw_execute_dft(planb,s1,s2)
      lu(0:nxx,0:nyy) = s2/(nx*ny)
      lu(nx,0:nyy) = lu(0,0:nyy)
      lu(0:nxx,ny) = lu(0:nxx,0)
      lu(nx,ny) = lu(0,0)
      return
    end subroutine calc_lap
end module derive

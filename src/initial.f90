module initial
  contains
    subroutine allocate_var ! allocate arrays and matrices
      use cdata
      implicit none

      allocate(x(0:nx))
      allocate(x2(0:nx))
      allocate(kx(0:nx))
      allocate(kx2(0:nx))
      allocate(y(0:ny))
      allocate(y2(0:ny))
      allocate(ky(0:ny))
      allocate(ky2(0:ny))

      allocate(v(0:nx,0:ny))
      allocate(cp(0:nx,0:ny))
      allocate(p2(0:nx,0:ny))

      allocate(lu(0:nx,0:ny))
      allocate(k2(0:nx,0:ny))
      
      allocate(du1(0:nx,0:ny))
      allocate(du2(0:nx,0:ny))
      allocate(du3(0:nx,0:ny))
      allocate(du4(0:nx,0:ny))
      allocate(cpv(0:nx,0:ny))
    end subroutine allocate_var
   
    subroutine free_var ! deallocate arrays and matrices
      use cdata
      implicit none
    
      if(allocated(x)) deallocate(x) 
      if(allocated(y)) deallocate(y) 
      if(allocated(x2)) deallocate(x2) 
      if(allocated(y2)) deallocate(y2) 
      if(allocated(kx)) deallocate(kx) 
      if(allocated(kx2)) deallocate(kx2) 
      if(allocated(ky)) deallocate(ky) 
      if(allocated(ky2)) deallocate(ky2) 

      if(allocated(v)) deallocate(v) 
      if(allocated(cp)) deallocate(cp) 
      if(allocated(p2)) deallocate(p2) 
      

      if(allocated(lu)) deallocate(lu) 
      if(allocated(k2)) deallocate(k2) 

      if(allocated(du1)) deallocate(du1) 
      if(allocated(du2)) deallocate(du2) 
      if(allocated(du3)) deallocate(du3) 
      if(allocated(du4)) deallocate(du4) 
      if(allocated(cpv)) deallocate(cpv) 
    end subroutine free_var

    subroutine init ! initialize space arrays and wave function
      use cdata
      implicit none
      integer :: i,j
      real(8),dimension((1+nx)*(1+ny)) :: randnum
      integer :: k
      integer,dimension(100) :: seed
      real(8) :: tmp,qr_al,theta
      do i = 0,nx ! initialize x
        x(i) = -lx/2.0d0 + (i) * dx
      end do
      do j = 0,ny ! initialize y
        y(j) = -ly/2.0d0 + (j) * dy
      end do

      do i = 1,nx2 ! initialize kx
        kx(i-1) = (i-1) * two_pi / (nx * dx)
        kx(nx2+i-1) = -(nx2-i+1) * two_pi / (nx * dx)
      end do
      do j = 1,ny2 ! initialize ky
        ky(j-1) = (j-1) * two_pi / (ny * dy)
        ky(ny2+j-1) = -(ny2-j+1) * two_pi / (ny * dy)
      end do

      x2 = x * x
      y2 = y * y
      kx2 = kx * kx
      ky2 = ky * ky
      
      do i = 0,nxx ; do j = 0,nyy
        k2(i,j) = kx2(i) + ky2(j)
      end do ; end do    


      select case(init_cond) ! provide initial conditon
        case('single_vortex') ! single vortex at center of harmonic trap with fixed vortex charge
          qr_al = sqrt(sqrt(gamma * nu)/pi)
          !$omp parallel do private(i,j,tmp,theta)
          do i = 0,nx ; do j = 0,ny
            tmp = 0.5d0 * (gamma * x2(i) + nu * y2(j))
            theta = atan2(y(j),x(i))
            cp(i,j) = qr_al * exp(-tmp) * (x(i) + ci*y(j)) * exp(ci * charge * theta)
          end do ; end do
          !$omp end parallel do
        case('random_phase') ! single vortex at center of harmonic trap with randomized phase
          seed = 13
          call random_seed(put=seed)
          call random_number(randnum)
          !$omp parallel do private(i,j,tmp) 
          do i = 0,nx ; do j = 0,ny
            k = (i+1) * (j+1)
            tmp = 0.5d0 * (gamma * x2(i) + nu * y2(j))
            cp(i,j) = qr_al * exp(-tmp) * (x(i) + ci*y(j)) * exp(2.0d0 * pi * ci * randnum(k))
          end do ; end do
          !$omp end parallel do
        end select
    end subroutine init
  
    subroutine calc_trap ! calculate potential
      use cdata
      implicit none
      integer :: i,j
      real(8) :: arg
      select case(potential)
      case('harmonic') ! harmonic potential
        !$omp parallel do private(i,j)
        do i = 0,nx ; do j = 0,ny
          v(i,j) = 0.5 * ( gamma2 * x2(i) + nu2 * y2(j) )
        end do ; end do
        !$omp end parallel do
      case('rotating') ! rotating potential
        arg = omega * t
        !$omp parallel do private(i,j)
        do i = 0,nx ; do j = 0,ny
          v(i,j) = 0.5 * ( gamma2 * x2(i) + nu2 * y2(j) + &
                  ecc *( x(i) * sin(arg) + y(j) * cos(arg))**2 - &
                  ecc * (x(i) * cos(arg) - y(j) * sin(arg))**2 )
        end do ; end do 
        !$omp end parallel do
      end select
    end subroutine calc_trap 
    

    subroutine set_dipole(cp,dicp) !setup dipole interactions
      use cdata,only : nx,ny,nxx,nyy,k2,two_pi,dkx,dky,rc,kx,cdd,p2,plandf,r1,s1,plandb,fftw_execute_dft_r2c,fftw_execute_dft_c2r
      use energy,only : integrate
      implicit none
      complex(8),dimension(0:nx,0:ny),intent(in) :: cp
      real(8),dimension(0:nx,0:ny),intent(out) :: dicp
      real(8),dimension(0:nxx,0:nyy) :: k,uddf
      real(8),dimension(0:nxx,0:nyy) :: delta,kappa,eta
      real(8) :: krc,udd,eps
      integer :: i,j

      k = abs(k2)  
      !$omp parallel do private(i,j,krc)
      do i = 0,nxx ; do j = 0,nyy
      if(k(i,j).ne.0.0) then
        delta(i,j) = 3 * (cos(rc * k(i,j)))/(rc * k(i,j))**2
        kappa(i,j) = 3 * (sin(rc * k(i,j)))/(rc * k(i,j))**3
        eta(i,j) = 3 * (kx(i)/k(i,j))**2
      else
        delta(i,j) = 3
        kappa(i,j) = 0
        eta(i,j) = 3 * kx(i)**2
      end if
        uddf(i,j) = 1 + (delta(i,j) - kappa(i,j)) * (eta(i,j) - 1)
      end do ; end do
      p2 = cp%re * cp%re + cp%im * cp%im
      r1 = p2(0:nxx,0:nyy) 
      call fftw_execute_dft_r2c(plandf,r1,s1)
      s1 = 0.33333 * cdd * uddf * s1 
      call fftw_execute_dft_c2r(plandb,s1,r1)
      dicp(0:nxx,0:nyy) = r1/(nx*ny)
      dicp(nx,0:nyy) = dicp(0,0:nyy)
      dicp(0:nxx,ny) = dicp(0:nxx,0)
      dicp(nx,ny) = dicp(0,0)
      return
   end subroutine set_dipole

   subroutine set_threads
     use cdata,only : threads
     use omp_lib
     implicit none
     if( threads/=0) then
      call omp_set_num_threads(threads)
     end if
     !$omp parallel
     !$omp master
       threads = omp_get_num_threads()
     !$omp end master
     !$omp end parallel  
   end subroutine set_threads

  subroutine set_fft
    use cdata
    implicit none
    planf = fftw_plan_dft_2d(nx,ny,s1,s2,fftw_forward,fftw_estimate)
    planb = fftw_plan_dft_2d(nx,ny,s1,s2,fftw_backward,fftw_estimate)
    plandf = fftw_plan_dft_r2c_2d(nx,ny,r1,s1,fftw_estimate)
    plandb = fftw_plan_dft_c2r_2d(nx,ny,s1,r1,fftw_estimate)
  end subroutine set_fft

  subroutine destroy_fft
    use cdata
    implicit none
    call fftw_destroy_plan(planf)
    call fftw_destroy_plan(planb)
    call fftw_destroy_plan(plandf)
    call fftw_destroy_plan(plandb)
  end subroutine destroy_fft
end module initial



module cdata
  use, intrinsic :: iso_c_binding
  implicit none
  include 'parameter.local'
  include 'fftw3.f03'
  real(8),dimension(:),allocatable :: x,x2,y,y2,kx,ky,kx2,ky2 ! space arrays
  real(8),dimension(:,:),allocatable :: v,p2,k2 ! potential, density, fourier space vector
  real(8) :: t ! time variable

  complex(8),dimension(:,:),allocatable :: cp,lu ! wavefunction and laplacian of wavefunction

  !---rk4_matrices---!
  complex(8),dimension(:,:),allocatable :: du1,du2,du3,du4,cpv

  !--fft_plans_and_matrices---!
  type(c_ptr) :: planf,planb,plandf,plandb
  complex(8),dimension(0:nxx,0:nyy) :: s1,s2
  real(8),dimension(0:nxx,0:nyy) :: r1,r2 

  !---energy_values---!
  real(8) :: zn ! Norm
end module 




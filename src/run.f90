program run
  use cdata
  use initial
  use output
  use energy
  use omp_lib
  use derive
  implicit none
  integer :: iter,i,s,step
  character(256) :: iter_file
  call allocate_var
  call set_threads
  call set_fft
  call init
  call calc_trap

  iter = nmax/snap
  step = 1
  t = 0
  !calculate initial density
  p2 = cp%re * cp%re + cp%im * cp%im
  open(1,file="density/initial_density.txt")
    call write_den(1,p2)
  close(1) 
  call norm(cp,zn)
  do s = 1,nmax ! rk4 time iteration begins
    if(time_dependent) call calc_trap !check whether potential is time-dependent
    
    call calc_du(cp,du1)
    cpv = cp + 0.5 * du1 * dt
    call calc_du(cpv,du2)
    cpv = cp + 0.5 * du2 * dt
    call calc_du(cpv,du3)
    cpv = cp + du3 * dt
    call calc_du(cpv,du4)

    cp = cp + (du1 + 2.0d0 * du2 + 2.0d0 * du3 + du4) * dt / 6.0d0 ! update wavefunction
    t = t+dt
    
    if(mod(s,iter) == 0) then !take snapshot of wavefunction
      write(iter_file,'(a,i0,a)') 'density/den-',step,'.txt'
      p2 = cp%re * cp%re + cp%im * cp%im
      open(2,file = iter_file)
        call write_den(2,p2)
      close(2)
      step = step + 1
      call norm(cp,zn)
    end if
  end do
  call free_var
  call destroy_fft
end program run

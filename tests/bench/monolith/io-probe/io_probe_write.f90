!> ORION's ASCII Tecplot writer on one block of n cells per side (the mesh and seven cell-centred variables, about
!  10 M values at n = 100), timed with the process's CPU times; one IPROBE line. Built four times (CMake): with and
!  without -assume buffered_io on this program, against ORION's Lib_Tecplot or a copy of it compiled with the flag.
!    io-probe-wr[-main|-mod|-both] <file.tec> <n>
program io_probe_write
#ifdef PROBE_TEC_COPY
  use Lib_Tecplot_bufio, only: tec_write_structured_multiblock
#else
  use Lib_Tecplot, only: tec_write_structured_multiblock
#endif
  use Lib_ORION_data, only: orion_data
  use, intrinsic :: iso_fortran_env, only: R8 => real64, I8 => int64
  use, intrinsic :: iso_c_binding, only: c_double
  implicit none
  interface
    subroutine probe_times(user, sys) bind(C, name='probe_times')
      import :: c_double
      real(c_double), intent(out) :: user, sys
    end subroutine probe_times
  end interface
  type(orion_data) :: orion
  character(len=4096) :: file, arg, exe
  integer :: n, i, j, k, s, err
  integer(I8) :: c0, c1, rate, nbytes, nval
  real(R8) :: u0, s0, u1, s1

  call get_command_argument(0, exe)
  call get_command_argument(1, file)
  call get_command_argument(2, arg)
  read(arg, *) n
  orion%tec%node = .false.; orion%tec%bc = .false.; orion%tec%format = 'ascii'
  orion%solutiontime = -10._R8
  allocate(orion%block(1))
  orion%block(1)%name = 'Block1'
  orion%block(1)%Ni = n; orion%block(1)%Nj = n; orion%block(1)%Nk = n
  allocate(orion%block(1)%mesh(1:3, 0:n, 0:n, 0:n), orion%block(1)%vars(1:7, 1:n, 1:n, 1:n))
  do k = 0, n; do j = 0, n; do i = 0, n
    orion%block(1)%mesh(:, i, j, k) = [real(i, R8), real(j, R8), real(k, R8)]*1.e-3_R8/3._R8
  enddo; enddo; enddo
  do k = 1, n; do j = 1, n; do i = 1, n; do s = 1, 7
    orion%block(1)%vars(s, i, j, k) = 10._R8**(s - 4)*(1._R8 + sin(real(i + 3*j + 7*k + s, R8)))
  enddo; enddo; enddo; enddo
  nval = size(orion%block(1)%mesh, kind=I8) + size(orion%block(1)%vars, kind=I8)

  call system_clock(c0, rate)
  call probe_times(u0, s0)
  err = tec_write_structured_multiblock(orion=orion, varnames='"a" "b" "c" "d" "e" "f" "g"', filename=trim(file))
  call system_clock(c1)
  call probe_times(u1, s1)
  inquire(file=trim(file), size=nbytes)
  i = index(exe, '/', back=.true.)
  write(*,'(A,A,A,I0,A,I0,A,I0,3(A,ES14.6E3),A,I0)') 'IPROBE step=buffered exe=', trim(exe(i+1:)), ' n=', n, &
    ' values=', nval, ' bytes=', nbytes, ' wall_s=', real(c1 - c0, R8)/real(rate, R8), ' user_s=', u1 - u0, &
    ' sys_s=', s1 - s0, ' err=', err
end program io_probe_write

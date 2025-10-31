!
! Include file to compute recycling and/or invariant subspace matrix U0
! U0 approximates the near-null space
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2025 Martin van Gijzen
!


if ( recycle ) then

! Initial solve to compute recycling subspace:
   call system_clock ( tb, clock_rate, clock_max )
   allocate( U0(size(b),s) )
   x0 = idrs_ritz( A, b, s, M1, tolerance = tol, maximum_iterations = maxit, variant = 1, &
        flag = flag, relres = relres, iterations = iter, U0=U0 )
   call system_clock ( te, clock_rate, clock_max )
   if ( my_proc == 1 ) then
      write(*,'(a)') 'Results for initial solve'
      write(*,'(a,f8.2,a)') 'Elapsed time           = ', real(te-tb)/real(clock_rate), 's.'
      write(*,'(a,i5)') 'Number of iterations   = ', iter
      write(*,'(a,e9.2)') 'Relative residual norm = ', relres
      if ( flag == 3 ) write(*,'(a)') 'Break down!'
      write(*,*)
   end if
end if

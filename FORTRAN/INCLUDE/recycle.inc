!
! Include file to compute recycling matrix U0
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2024 Martin van Gijzen
!

if ( recycle ) then
   call system_clock ( tb, clock_rate, clock_max )
   x0 = idrs_ritz( A, b, s, M, tol, maxit, 1, flag, relres, iter, x0, U0 )
   call system_clock ( te, clock_rate, clock_max )
   if ( this_image() == 1 ) then
      write(*,'(a)') 'Results for initial solve'
      write(*,'(a,f8.2,a)') 'Elapsed time           = ', real(te-tb)/real(clock_rate), 's.'
      write(*,'(a,i5)') 'Number of iterations   = ', iter
      write(*,'(a,e9.2)') 'Relative residual norm = ', relres
      if ( flag == 3 ) write(*,'(a)') 'Break down!'
      write(*,*)
   end if
end if

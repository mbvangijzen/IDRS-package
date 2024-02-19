!
! Include file to call the idrs-solvers and print some output
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2024 Martin van Gijzen
!

call system_clock ( tb, clock_rate, clock_max )
if ( use_idrs .or. use_bicgstab ) then
   x = idrs( A, b, s, M, tol, maxit, variant, flag, relres, iter, x0, U0, omega )
elseif ( use_qmridr ) then
   x = qmridr( A, b, s, M, tol, maxit, flag, relres, iter, in_s, in_tol, in_it, x0, omega )
else
   stop 'Error: no solution method specified!'
end if
call system_clock ( te, clock_rate, clock_max )

if ( this_image() == 1 ) then
   write(*,'(a,f8.2,a)') 'Elapsed time           = ', real(te-tb)/real(clock_rate), 's.'
   write(*,'(a,i5)') 'Number of iterations   = ', iter
   write(*,'(a,e9.2)') 'Relative residual norm = ', relres
   if ( flag > 0 ) then
      if ( flag == 1 ) write(*,'(a)') 'Maximum number of iterations reached!'
      if ( flag == 2 ) write(*,'(a)') 'Accuracy above prescribed tolerance!'
      if ( flag == 3 ) write(*,'(a)') 'Break down!'
   end if
   write(*,*)
end if

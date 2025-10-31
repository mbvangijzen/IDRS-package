!
! Include file to call the multishift idrs-solvers and print some output
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2025 Martin van Gijzen
!

!
! Compute shifts and parameters of preconditioned problem:
call ms_parameters( A, sigma, sigma_p )

call system_clock ( tb, clock_rate, clock_max )
if ( use_idrs .or. use_bicgstab .or. use_minsync ) then
   if ( plot_conv ) then
      x = msidrs( A, Pmul(A,b), sigma_p, s, M1, tol, maxit, variant, flag, resnrm, iter, omega, ms_resvec )
   else
      x = msidrs( A, Pmul(A,b), sigma_p, s, M1, tol, maxit, variant, flag, resnrm, iter, omega )
   end if
elseif ( use_qmridr ) then
   if ( plot_conv ) then
      x = msqmridr( A, Pmul(A,b), sigma_p, s, M1, tol, maxit, flag, resnrm, iter, in_s, in_tol, in_it, omega, ms_resvec )
   else
      x = msqmridr( A, Pmul(A,b), sigma_p, s, M1, tol, maxit, flag, resnrm, iter, in_s, in_tol, in_it, omega )
   end if
else
   error stop 'Error: no solution method specified!'
end if

call system_clock ( te, clock_rate, clock_max )

! Report performance measures:
if ( my_proc == 1 ) then
   write(*,'(a,f8.2,a)') 'Elapsed time           = ', real( te - tb )/real( clock_rate ), 's.'
   write(*,'(a,i5)') 'Number of iterations   = ', iter
   write(*,'(a,*(e9.2))') 'Relative residual norm = ', resnrm
   if ( flag > 0 ) then
      if ( flag == 1 ) write(*,'(a)') 'Maximum number of iterations reached!'
      if ( flag == 2 ) write(*,'(a)') 'Accuracy above prescribed tolerance!'
      if ( flag == 3 ) write(*,'(a)') 'Break down!'
   end if
end if

! Compute solution of unpreconditioned system:
x = scaleback( A, sigma, b, x )

! Store the solution:
if ( present(C) ) then
   solution = x(nn+1:neq,:)
else
   solution = x
end if

normb = norm(rhs)
if ( plot_conv )   ms_resvec = ms_resvec/normb

! Check the relative residuals of original systems
if ( present(M) .and. present(C) ) then
   do i_sigma = 1, Nsigma
      resnrm(i_sigma) = &
         norm( rhs - K*solution(:,i_sigma) - &
         sigma(i_sigma)*(C*solution(:,i_sigma))-sigma(i_sigma)**2*(M*solution(:,i_sigma)) )/normb
   end do
elseif ( present(C) ) then
   do i_sigma = 1, Nsigma
      resnrm(i_sigma) = &
         norm( rhs - K*solution(:,i_sigma) - &
         sigma(i_sigma)*(C*solution(:,i_sigma))-sigma(i_sigma)**2*(solution(:,i_sigma)) )/normb
   end do
elseif ( present(M) ) then
   do i_sigma = 1, Nsigma
      resnrm(i_sigma) = norm( rhs - K*solution(:,i_sigma) + sigma(i_sigma)*(M*solution(:,i_sigma)) )/normb
   end do
else
   do i_sigma = 1, size(sigma)
      resnrm(i_sigma) = norm( rhs - K*solution(:,i_sigma) + sigma(i_sigma)*solution(:,i_sigma) )/normb
   end do
end if
if ( my_proc == 1 ) write(*,'(a,*(e9.2))') 'Relative residual norms unpreconditioned systems = ', resnrm

! Plot convergence?
if ( plot_conv .and. my_proc == 1 ) then
   if ( preconditioner > 0 ) then
      write(title,'(a,a,a,a)') 'Convergence for ', trim(method), '. ', trim(pol_prec)
   else
      write(title,'(a,a)') 'Convergence for ', trim(method)
   end if
   open( 20, file = 'convergence.plt' )
   if ( gif ) then
      write(20,*) 'set terminal gif'
      write(20,*) 'set output "convergence.gif"' 
   elseif ( jpeg ) then
      write(20,*) 'set terminal jpeg'
      write(20,*) 'set output "convergence.jpeg"'
   end if
   write(20,*) 'set logscale y 10'
   write(20,*) 'set format y "%T"'
   write(20,*) 'set xlabel "Iteration"'
   write(20,*) 'set ylabel "log(||r||/||b||)"'
   write(20,*) 'set grid'
   write(20,*) 'set title "', trim(title), '"'
   write(20,*) 'plot for [i=1:*] "convergence.dat" using i with linespoints title "Shift ".i'
   close(20)
   open( 30, file = 'convergence.dat' )
   do i_conv = 1, iter+1
      write(30,*) ms_resvec(i_conv,:)
   end do
   close(30)
   call execute_command_line('gnuplot -p convergence.plt')
end if
if ( plot_conv ) deallocate( ms_resvec )

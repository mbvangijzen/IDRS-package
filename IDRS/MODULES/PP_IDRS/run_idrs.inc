!
! Include file to call the idrs-solvers and print some output
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2025 Martin van Gijzen
!

call system_clock ( tb, clock_rate, clock_max )
if ( use_idrs .or. use_bicgstab .or. use_minsync ) then
   if ( plot_conv ) then
      x = idrs( A, Pmul(A,b), s, M1, tol, maxit, variant, flag, relres, iter, x0, U0, omega, resvec=resvec )
   else
      x = idrs( A, Pmul(A,b), s, M1, tol, maxit, variant, flag, relres, iter, x0, U0, omega )
   end if
elseif ( use_qmridr ) then
   if ( plot_conv ) then
      x = qmridr( A, Pmul(A,b), s, M1, tol, maxit, flag, relres, iter, in_s, in_tol, in_it, x0, omega, resvec )
   else
      x = qmridr( A, Pmul(A,b), s, M1, tol, maxit, flag, relres, iter, in_s, in_tol, in_it, x0, omega )
   end if
else
   stop 'Error: no solution method specified!'
end if
call system_clock ( te, clock_rate, clock_max )

! Report performance measures:
if ( my_proc == 1 ) then
   write(*,'(a,f8.2,a)') 'Solution time          = ', real(te-tb)/real(clock_rate), 's.'
   write(*,'(a,i5)')     'Number of iterations   = ', iter
   write(*,'(a,e9.2)')   'Relative residual norm = ', relres
   if ( flag > 0 ) then
      if ( flag == 1 ) write(*,'(a)') 'Maximum number of iterations reached!'
      if ( flag == 2 ) write(*,'(a)') 'Accuracy above prescribed tolerance!'
      if ( flag == 3 ) write(*,'(a)') 'Break down!'
   end if
end if

! Save previous solution as initial starting vector: 
!(note that this is the solution of the (right-) preconditioned system)
x0 = x

! Compute solution of unpreconditioned system:
x = scaleback( A, b, x0 )

normb = norm(b)
if ( plot_conv ) resvec = resvec/normb

if ( .not. present(seed) ) then
   relres = norm( b - K*x )/normb
elseif ( present(M) .and. present(C) ) then
   relres = norm( b - K*x -seed*(C*x)-seed**2*(M*x) )/normb
elseif ( present(C) ) then
   relres = norm( b - K*x -seed*(C*x)-seed**2*(  x) )/normb
elseif ( present(M) ) then
   relres = norm( b - K*x + seed*(M*x) )/normb
else
   relres = norm( b - K*x + seed*x )/normb
end if
if ( my_proc == 1 ) write(*,'(a,e9.2,/)') 'Relative residual norm original system = ', relres

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
   write(20,*) 'unset key'
   write(20,*) 'set title "', trim(title), '"'
   write(20,*) 'plot "convergence.dat" with linespoints'
   close(20)

   open( 30, file = 'convergence.dat' )
   do i_conv = 1, iter+1
      write(30,*) resvec(i_conv)
   end do
   close(30)
   call execute_command_line('gnuplot -p convergence.plt')

end if

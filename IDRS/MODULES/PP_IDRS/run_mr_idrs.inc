!
! Include file to call the idrs-solvers and print some output
! Solves sequence of systems with different rhs-vectors
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2025 Martin van Gijzen
!

! Open file to store convergence curves:
if ( my_proc == 1 .and. plot_conv ) open( 10, file = 'convergence.dat' )

do irhs = 1, nrhs

   normb = norm(b_all(:,irhs))

! Deflated rhs-vector:
   b = Pmul(A,b_all(:,irhs))
!
   if ( my_proc == 1 ) write(*,'(a,i2)') 'Results for rhs ', irhs

   call system_clock ( tb, clock_rate, clock_max )
   if ( use_idrs .or. use_bicgstab .or. use_minsync ) then
      if ( plot_conv ) then
         x = idrs( A, b, s, M1, tol, maxit, variant, flag, relres, iter, x0, U0, omega, resvec=resvec )
      else
         x = idrs( A, b, s, M1, tol, maxit, variant, flag, relres, iter, x0, U0, omega )
      end if
   elseif ( use_qmridr ) then
      if ( plot_conv ) then
         x = qmridr( A, b, s, M1, tol, maxit, flag, relres, iter, in_s, in_tol, in_it, x0, omega, resvec )
      else
         x = qmridr( A, b, s, M1, tol, maxit, flag, relres, iter, in_s, in_tol, in_it, x0, omega )
      end if
   else
      stop 'Error: no solution method specified!'
   end if
   call system_clock ( te, clock_rate, clock_max )

! Report performance measures:
   if ( my_proc == 1 ) then
      write(*,'(a,f8.2,a)') 'Elapsed time           = ', real(te-tb)/real(clock_rate), 's.'
      write(*,'(a,i5)') 'Number of iterations   = ', iter
      write(*,'(a,e9.2)') 'Relative residual norm = ', relres
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
   x_all(:,irhs) = scaleback( A, b_all(:,irhs), x0 )

   if ( .not. present(seed) ) then 
      relres = norm( b_all(:,irhs) - K*x_all(:,irhs) )/normb
   elseif ( present(M) .and. present(C) ) then
      relres = norm( b_all(:,irhs) - K*x_all(:,irhs) -seed*(C*x_all(:,irhs))-seed**2*(M*x_all(:,irhs)) )/normb
   elseif ( present(C) ) then
      relres = norm( b_all(:,irhs) - K*x_all(:,irhs) -seed*(C*x_all(:,irhs))-seed**2*(  x_all(:,irhs)) )/normb
   elseif ( present(M) ) then
      relres = norm( b_all(:,irhs) - K*x_all(:,irhs) + seed*(M*x_all(:,irhs)) )/normb
   else
      relres = norm( b_all(:,irhs) - K*x_all(:,irhs) + seed*x_all(:,irhs) )/normb
   end if
   if ( my_proc == 1 ) write(*,'(a,e9.2,/)') 'Relative residual norm original system = ', relres

! Store convergence curves:
   if ( my_proc == 1 .and. plot_conv ) then
!     write( 10, '(a,i2.2)' )'RHS', irhs
      do i_conv = 1, iter+1
         write(10,*) resvec(i_conv)/normb
      end do
      write(10,*)
      write(10,*)
   end if
end do

! Close file to store convergence data:
if ( my_proc == 1 .and. plot_conv ) close(10)

! Make plot file for convergence plots:
if ( my_proc == 1 .and. plot_conv ) then

   if ( preconditioner > 0 ) then 
      write(title,'(a,a,a,a)') 'Convergence for ', trim(method), '. ', trim(pol_prec)
   else 
      write(title,'(a,a)') 'Convergence for ', trim(method)
   end if
   open( 20, file = 'convergence.plt' )
   if ( gif ) then
      write(20,'(a)') 'set terminal gif font "sans, 5"'
      write(20,'(a)') 'set output "convergence.gif'
   elseif ( jpeg ) then
      write(20,'(a)') 'set terminal jpeg font "sans, 5"'
      write(20,'(a)') 'set output "convergence.jpeg'
   else
      write(20,'(a)') 'set tics font "sans, 5"'
   end if
   write(20,'(a)') 'set logscale y 10'
   write(20,'(a)') 'set format y "%T"'
   write(20,'(a)') 'set xlabel "Iteration"'
   write(20,'(a)') 'set ylabel "log(||r||/||b||)"'
   write(20,'(a)') 'set grid'
   write(20,'(a)') 'set nokey'
   if ( mod(nrhs,3) == 0 ) then
      write(20,'(a,a,a,i2)') 'set multiplot title "', trim(title), '" layout 3,',nrhs/3
   elseif ( mod(nrhs,2) == 0 ) then
      write(20,'(a,a,a,i2)') 'set multiplot title "', trim(title), '" layout 2,',nrhs/2
   else 
      write(20,'(a,a,a)') 'set multiplot title "', trim(title), '"'
   end if
   do irhs = 1, nrhs
      write(20,'(a,i2,a)') 'set title "RHS', irhs,'"'
      write(20,'(a,i2,a)') 'plot "convergence.dat" index ',irhs-1,' with linespoints '
   end do
   close(20)
   call execute_command_line('gnuplot -p convergence.plt')
end if

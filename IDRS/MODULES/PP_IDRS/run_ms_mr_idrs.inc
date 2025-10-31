!
! Include file to call the multishift idrs-solvers and print some output
! This is the multiple rhs version
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2025 Martin van Gijzen
!

!
! Compute shifts and parameters of preconditioned problem:
call ms_parameters( A, sigma, sigma_p )

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
         x = msidrs( A, b, sigma_p, s, M1, tol, maxit, variant, flag, resnrm, iter, omega, ms_resvec )
      else
         x = msidrs( A, b, sigma_p, s, M1, tol, maxit, variant, flag, resnrm, iter, omega )
      end if
   elseif ( use_qmridr ) then
      if ( plot_conv ) then
         x = msqmridr( A, b, sigma_p, s, M1, tol, maxit, flag, resnrm, iter, in_s, in_tol, in_it, omega, ms_resvec )
      else
         x = msqmridr( A, b, sigma_p, s, M1, tol, maxit, flag, resnrm, iter, in_s, in_tol, in_it, omega )
      end if
   else
      error stop 'Error: no solution method specified!'
   end if

   call system_clock ( te, clock_rate, clock_max )

! Report performance measures:
   if ( my_proc == 1 ) then
      write(*,'(a,f8.2,a)')  'Elapsed time           = ', real( te - tb )/real( clock_rate ), 's.'
      write(*,'(a,i5)')      'Number of iterations   = ', iter
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
      solution(:,:,irhs) = x(nn+1:neq,:)
   else
      solution(:,:,irhs) = x
   end if

   if ( plot_conv )   ms_resvec = ms_resvec/normb

! Check the relative residuals of original systems
   if ( present(M) .and. present(C) ) then
      do i_sigma = 1, Nsigma
         resnrm(i_sigma) = &
            norm( rhs(:,irhs) - K*solution(:,i_sigma,irhs) - &
            sigma(i_sigma)*(C*solution(:,i_sigma,irhs))-sigma(i_sigma)**2*(M*solution(:,i_sigma,irhs)) )/normb
      end do
   elseif ( present(C) ) then
      do i_sigma = 1, Nsigma
         resnrm(i_sigma) = &
            norm( rhs(:,irhs) - K*solution(:,i_sigma,irhs) - &
            sigma(i_sigma)*(C*solution(:,i_sigma,irhs))-sigma(i_sigma)**2*(solution(:,i_sigma,irhs)) )/normb
      end do
   elseif ( present(M) ) then
      do i_sigma = 1, Nsigma
         resnrm(i_sigma) = norm( rhs(:,irhs) - K*solution(:,i_sigma,irhs) + sigma(i_sigma)*(M*solution(:,i_sigma,irhs)) )/normb
      end do
   else
      do i_sigma = 1, size(sigma)
         resnrm(i_sigma) = norm( rhs(:,irhs) - K*solution(:,i_sigma,irhs) + sigma(i_sigma)*solution(:,i_sigma,irhs) )/normb
      end do
   end if
   if ( my_proc == 1 ) then 
      write(*,'(a,*(e9.2))') 'Relative residual norms unpreconditioned systems = ', resnrm
      write(*,*)
   end if

! Store convergence curves:
   if ( my_proc == 1 .and. plot_conv ) then
      do i_conv = 1, iter+1
         write(10,*) ms_resvec(i_conv,:)/normb
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
   if ( mod(nrhs,3) == 0 ) then
      write(20,'(a,a,a,i2)') 'set multiplot title "', trim(title), '" layout 3,',nrhs/3
   elseif ( mod(nrhs,2) == 0 ) then
      write(20,'(a,a,a,i2)') 'set multiplot title "', trim(title), '" layout 2,',nrhs/2
   else
      write(20,'(a,a,a)') 'set multiplot title "', trim(title), '"'
   end if
   do irhs = 1, nrhs
      write(20,'(a,i2,a)')   'set title "RHS', irhs,'"'
      write(20,'(a,i2,a)') 'plot for [i=1:*] "convergence.dat" index ', irhs-1, ' using i with linespoints title "Shift ".i'
   end do
   close(20)
   call execute_command_line('gnuplot -p convergence.plt')
end if
if ( plot_conv ) deallocate( ms_resvec )

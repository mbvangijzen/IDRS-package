!
! Include file to determine real polynomial preconditioner or shift-and-invert preconditioner
! and/or user defined omega's
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2025 Martin van Gijzen
!

! Polynomial preconditioner:

if ( present(center) ) rcenter = center
if ( present(foci) )   rfoci   = foci

! Estimation of parameters, or are they given?
estimate_parameters = .true.
if ( preconditioner == 0 ) estimate_parameters = .false.
if ( preconditioner == 1 .and. present(center) ) estimate_parameters = .false.
if ( preconditioner == 2 .and. present(foci) ) estimate_parameters = .false.
initial_iterations = ( estimate_parameters .or. ( n_omega > 0 ) )

if ( initial_iterations ) then
   ini_iter = n_ritz + n_ritz/s
   allocate( omega(n_omega) )
   allocate( H(ini_iter+1,ini_iter) )
   x0 = idrs_ritz( A, b, s, M1, tol, ini_iter, 1, flag, relres, iter, omega = omega, H = H )
   n_ritz = (iter/(s+1))*s + mod(iter,s+1)
end if

if ( estimate_parameters ) then
   call real_preconditioner( A, preconditioner, degree, H = H(1:n_ritz,1:n_ritz) )
else
   call real_preconditioner( A, preconditioner, degree, rfoci, rcenter )
end if 

! Plot the Ritz values
if ( initial_iterations .and. plot_ritz .and. my_proc == 1 ) then

! Make plot file:
   write(title,'(a)') 'Ritz values'
   open( 20, file = 'ritz.plt' )
   if ( gif ) then
      write(20,*) 'set terminal gif '
      write(20,*) 'set output "ritz.gif"' 
   elseif ( jpeg ) then
      write(20,*) 'set terminal jpeg'
      write(20,*) 'set output "ritz.jpeg"'
   end if
   write(20,*) 'set xlabel "Real"'
   write(20,*) 'set ylabel "Imag"'
   write(20,*) 'set grid'
   write(20,*) 'unset key'
   write(20,*) 'set title "', trim(title), '"'
   write(20,*) 'plot "ritz.dat" with points'
   close(20)

! Make data file:
   allocate( ritzval(n_ritz), outlier(n_ritz) )
   open( 30, file = 'ritz.dat' )
   ritzval = QR_eig( H(1:n_ritz,1:n_ritz) )
   outlier = filter_ritzval( ritzval )
   do i_ritz = 1, n_ritz
      if ( .not. outlier(i_ritz) ) write(30,*) real(ritzval(i_ritz)), aimag(ritzval(i_ritz))
   end do
   close(30)

! Run:
   call execute_command_line('gnuplot -p ritz.plt')
end if


!
! Include file to call idrs-eigenvalue solver and print some output
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2026 Martin van Gijzen
!

call system_clock ( tb, clock_rate, clock_max )

allocate( resnrm(neig), resnrm_original(neig), U(neq,neig) )

lambda = idreig( A, b, s, neig, M1, tol, maxit, variant, U, flag, resnrm, iter, resvec=eig_resvec )

! Scale back the eigenvectors:
do ieig = 1, neig
   U(:,ieig) = U(:,ieig)/A
end do

call system_clock ( te, clock_rate, clock_max )

! Report performance measures:
if ( my_proc == 1 ) then
   write(*,'(a,f8.2,a)') 'Elapsed time           = ', real( te - tb )/real( clock_rate ), 's.'
   write(*,'(a,i5)')     'Number of iterations   = ', iter
   if ( flag == 1 ) write(*,'(a)') 'Maximum number of iterations reached!'
   if ( flag == 2 ) write(*,'(a)') 'Accuracy above prescribed tolerance!'
   if ( flag == 3 ) write(*,'(a)') 'Break down!'
end if

! Compute the eigenvalues of the original system via Rayleigh quotients:

! First remove the preconditioner from the system matrix:
! call remove_preconditioner( A )

! Now calculte rayleigh quotients to determine the eigenvalues:
! do ieig = 1, neig
!    lambda(ieig) = inprod( U(:,ieig), A*U(:,ieig) ) / inprod( U(:,ieig), U(:,ieig) )
! end do 

if ( present(M) .and. present(C) ) then
   do ieig = 1, neig
      lambda(ieig) = &
         ( -inprod(U(1:nn,ieig),C*U(1:nn,ieig)) - inprod(U(1:nn,ieig),K*U(nn+1:2*nn,ieig)) &
           +inprod(U(nn+1:2*nn,ieig),U(1:nn,ieig)) ) / &
         (  inprod(U(1:nn,ieig),M*U(1:nn,ieig)) + inprod(U(nn+1:2*nn,ieig),U(nn+1:2*nn,ieig)) ) 
   end do
elseif ( present(C) ) then
   do ieig = 1, neig
      lambda(ieig) = &
         ( -inprod(U(1:nn,ieig),C*U(1:nn,ieig)) - inprod(U(1:nn,ieig),K*U(nn+1:2*nn,ieig)) &
           +inprod(U(nn+1:2*nn,ieig),U(1:nn,ieig)) ) / inprod(U(:,ieig),U(:,ieig)) 
   end do
elseif ( present(M) ) then
   do ieig = 1, neig
      lambda(ieig) = &
         inprod(U(:,ieig),K*U(:,ieig)) / inprod(U(:,ieig),M*U(:,ieig))
   end do
else
   do ieig = 1, neig
      lambda(ieig) = &
         inprod(U(:,ieig),K*U(:,ieig))  / inprod(U(:,ieig),U(:,ieig))
   end do
end if

! Store the eigenvectors
if ( present(C) ) then
   do ieig = 1, neig
      eigvec(:,ieig) = U(nn+1:neq,ieig)
   end do
else
   do ieig = 1, neig
      eigvec(:,ieig) = U(:,ieig)
   end do
end if

! Normalise the eigenvectors
if ( present(M) ) then
   do ieig = 1, neig
      normu = sqrt( abs(inprod( eigvec(:,ieig), M*eigvec(:,ieig)) ) )
      eigvec(:,ieig) = eigvec(:,ieig)/normu
   end do
else
   do ieig = 1, neig
      normu = norm(eigvec(:,ieig))
      eigvec(:,ieig) = eigvec(:,ieig)/normu
   end do
end if

! Check the relative residuals of original systems
if ( present(M) .and. present(C) ) then
   do ieig = 1, neig
      resnrm_original(ieig) = &
         norm( K*eigvec(:,ieig) + &
         lambda(ieig)*(C*eigvec(:,ieig)) + lambda(ieig)**2*(M*eigvec(:,ieig)) )
   end do
elseif ( present(C) ) then
   do ieig = 1, neig
      resnrm_original(ieig) = &
         norm( K*eigvec(:,ieig) + &
         lambda(ieig)*(C*eigvec(:,ieig)) + lambda(ieig)**2*(  eigvec(:,ieig)) )
   end do
elseif ( present(M) ) then
   do ieig = 1, neig
      resnrm_original(ieig) = norm( K*eigvec(:,ieig) - lambda(ieig)*(M*eigvec(:,ieig)) )
   end do
else
   do ieig = 1, neig
      resnrm_original(ieig) = norm( K*eigvec(:,ieig) - lambda(ieig)*(  eigvec(:,ieig)) )
   end do
end if

if ( my_proc == 1 ) then
   write(*,*)
   write(*,'(a)') '      Eigenvalue       | Residual norm transformed system | Residual norm original system '
   write(*,'(a)') '------------------------------------------------------------------------------------------'
   do i = 1, neig
      ei = lambda(i)
      write(*,'(1x,a,e9.3,a,e9.3,a,1x,a,12x,e10.4,12x,a,12x,e10.4)') &
         '(', real(ei),',',aimag(ei), ')','|',resnrm(i),'|', resnrm_original(i)
   end do
   write(*,'(a)') '------------------------------------------------------------------------------------------'
   write(*,*)
end if

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
   write(20,*) 'set ylabel "log(||r||)"'
   write(20,*) 'set grid'
   write(20,*) 'set title "', trim(title), '"'
   write(20,*) 'plot for [i=1:*] "convergence.dat" using i with linespoints title "Eigenvalue ".i'
   close(20)
   open( 30, file = 'convergence.dat' )
   do i_conv = 1, iter
      write(30,*) eig_resvec(i_conv,:)
   end do
   close(30)
   call execute_command_line('gnuplot -p convergence.plt')
end if
if ( plot_conv ) deallocate( eig_resvec )

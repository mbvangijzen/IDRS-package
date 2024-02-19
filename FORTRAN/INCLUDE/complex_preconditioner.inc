!
! Include file to determine complex polynomial preconditioner and/or user-defined omega
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2024 Martin van Gijzen
!

if ( n_ritz > 0 .or. n_omega > 0 ) then
! Initial iterations
   allocate( omega(n_omega) )
   allocate( H(n_ritz+1,n_ritz) )
   ini_iter = n_ritz + ceiling(real(n_ritz)/real(s))
   x0 = idrs_ritz( A, b, s, M, tol, ini_iter, 1, flag, relres, iter, omega = omega, H = H )
end if

! Make the polynomial preconditioner:
if ( n_ritz > 0 ) then
   if ( allocated(D) ) then
      call complex_preconditioner( M, A%A, D, preconditioner, degree, H = H )
   else
      call complex_preconditioner( M, A%A, preconditioner=preconditioner, degree=degree, H = H )
   end if
else
   if ( allocated(D) ) then
      call complex_preconditioner( M, A%A, D, preconditioner, degree, cfoci, ccenter, cseed )
   else
      call complex_preconditioner( M, A%A, preconditioner=preconditioner, degree=degree, &
         foci=cfoci, center=ccenter, seed=cseed )
   end if
end if 


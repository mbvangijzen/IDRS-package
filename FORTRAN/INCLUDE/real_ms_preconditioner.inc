!
! Include file to determine complex preconditioner and/or user defined omega's
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2024 Martin van Gijzen
!

if ( n_ritz > 0 .or. n_omega > 0 ) then
   allocate( omega(n_omega) )
   allocate( H(n_ritz+1,n_ritz) )
   ini_iter = n_ritz + ceiling(real(n_ritz)/real(s))
   x0 = idrs_ritz( A, b, s, M, tol, ini_iter, 1, flag, relres, iter, omega = omega, H = H )
end if

! Make the polynomial preconditioner:
if ( n_ritz > 0 ) then
   if ( allocated(D) ) then
      call real_preconditioner( A, D=D, preconditioner=preconditioner, degree=degree, &
         H = H, sigma = sigma, sigma_p = sigma_p )
   else
      call real_preconditioner( A, preconditioner=preconditioner, degree=degree, &
         H = H,  sigma = sigma, sigma_p = sigma_p )
   end if
else
   if ( allocated(D) ) then
      call real_preconditioner( A, D=D, preconditioner=preconditioner, degree=degree, &
         foci=rfoci, center=rcenter, seed=rseed, sigma = sigma, sigma_p = sigma_p  )
   else
      call real_preconditioner( A, preconditioner=preconditioner, degree=degree, &       
         foci=rfoci, center=rcenter, seed=rseed, sigma = sigma, sigma_p = sigma_p )
   end if
end if


program ms_stencil
!
! This programme illustrates the use of the idrs-package on
! a convection-diffusion-reaction problem with constant coefficients
! and homogeneous dirichlet boundary conditions.
!
! The problem can be  solved with the different algorithms: IDRS, BiCGSTAB and QMRIDR,
! Neumann and Chebyshev polynomial preconditioners, and with and without recycling techniques.
! The problem is solved with the different idrs-algorithms: IDRS and QMRIDR 
! for the single shift (= reaction parameter) problem and MSIDRS and MSQMRIDR 
! for the multishift problem, which uses several reaction parameters.
!
! The matrix is not explicitly stored. The matrix-vector multiplication
! is done using a stencil operations.
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2024 Martin van Gijzen
!

   use interface_module          ! To read command line options
   use precision_module          ! Set real and complex precision
   use matrix_module             ! Always needed is a matrix module that defines the
                                 ! matrix-vector multiplication and preconditioning operation
   use user_module               ! User functions, user defined matrix-vector product
   use idrs_module               ! The module with the IDRS routines
   use ritz_module               ! Can be used for computing an initial search space 
                                 ! and user defined omegas

   implicit none 

   integer                       :: grid
   integer, parameter            :: default_grid = 3
!
   integer, parameter            :: Nsigma = 6              ! Number of sigma
   real(kind=rp)                 :: sigma(Nsigma)           ! Shifts
   real(kind=rp)                 :: dsigma                  ! Shift spacing 
   real(kind=rp)                 :: sigma_max               ! Maximum shift
   real, parameter               :: default_sigma_max = 5.  ! Maximum shift
!
   type(user_matrix)             :: A1, M1

   real(kind=rp), allocatable    :: x(:,:), x0(:), b(:)

   integer                       :: tb, te, clock_rate, clock_max  ! For system_clock
   integer                       :: my_proc, n_procs

   real(kind=rp)                 :: radius
   real(kind=rp)                 :: true_relres(Nsigma), normb
   integer                       :: i

! Declaration and initialization of IDRS parameters:
   include "../INCLUDE/initialize.inc"

   my_proc = this_image()
   n_procs = num_images()

! Which grid?
   grid = get_parameter( '-grid', default_grid )

! Make the matrices
   call cdr_matrix( A1, grid )
   call mass_matrix( A1, M1 )

! Output message, problem description:
   if ( this_image() == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      write(*,'(6x,a,i10,/)') &
         'Convection-diffusion-reaction problem, number of equations is ', n_procs*A1%Nn
   end if
!
! Output message, which methods:
   include "../INCLUDE/output_methods.inc"

! Allocate solution and right-hand-side vectors:
   allocate( x0(A1%Nn), b(A1%Nn), x(A1%Nn,Nsigma), resnrm(Nsigma) )

   A1%multishift = .true.

! Define the shifts:
   sigma_max = get_parameter( '-reaction', default_sigma_max )
   dsigma = 0.
   if ( Nsigma > 1 ) dsigma = sigma_max/(Nsigma-1)
   do i = 1, Nsigma
      sigma(i) = (i-1)*dsigma
   end do

! Compute right-hand-side
   call model_solution( x0, A1 )
   b = A1*x0
   normb = norm(b)

! Store the matrix in correct form
   A = real_matrix( A1 )

! Compute estimates for largest and smallest eigenvalues with Gershgorin's theorem:
   call gershgorin( M1, rcenter, radius )
   rfoci(1) = rcenter-radius
   rfoci(2) = rcenter+radius  
! Approximate inverse of Mass matrix by polynomial:
   call real_preconditioner( M, M1, preconditioner=preconditioner, degree=degree, &
         foci=rfoci, center=rcenter )

!
! Test the multishift algorithms
   call system_clock ( tb, clock_rate, clock_max )
   if ( use_idrs ) then
      variant = 1
      x = MSIDRS( A, b, sigma, s, M, tol, maxit, variant, flag, resnrm, iter )
   elseif ( use_bicgstab ) then
      variant = 2
      x = MSIDRS( A, b, sigma, s, M, tol, maxit, variant, flag, resnrm, iter )
   elseif ( use_qmridr ) then
      x = MSQMRIDR( A, b, sigma, s, M, tol, maxit, flag, resnrm, iter, in_s, in_tol, in_it )
   else
      error stop 'Error: no solution method specified!'
   end if
   call system_clock ( te, clock_rate, clock_max )

!
! Compute true residual norms
   do i = 1, Nsigma
      true_relres(i) = norm( b - A1*x(:,i) + sigma(i)*(M1*x(:,i)) )/normb
   end do

   if ( my_proc == 1 ) then
      write(*,'(a)') 'Results: '
      write(*,'(a,f8.2)') 'Elapsed time         = ', real ( te - tb ) / real ( clock_rate )
      write(*,'(a,i4)')   'Number of iterations = ', iter
      do i = 1, Nsigma
         write(*,'(a,f8.2,2x,a,e10.2)') 'Shift = ', sigma(i), '||r||/||b|| = ', resnrm(i)
         write(*,'(a,e10.2)') 'True relative residual norm   = ', true_relres(i)
      end do
      if ( flag > 0 ) then
         if ( flag == 1 ) write(*,*) 'Maximum number of iterations reached!'
         if ( flag == 2 ) write(*,*) 'Accuracy above prescribed tolerance!'
         if ( flag == 3 ) write(*,*) 'Break down!'
      end if
      write(*,'(a)') &
         '==============================================================================='
   end if

end program ms_stencil

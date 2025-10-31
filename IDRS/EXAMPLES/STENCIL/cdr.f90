program cdr

!
! This programme illustrates the use of the idrs-package on
! a convection-diffusion-reaction problem with constant coefficients
! and homogeneous dirichlet boundary conditions.
! The problem can be  solved with the different algorithms: MSIDRS, MSBiCGSTAB and MSQMRIDR,
! and Neumann and Chebyshev polynomial preconditioners.
!
! The matrices are not explicitly stored. The matrix-vector multiplication
! is done using a stencil operations.
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2024 Martin van Gijzen
!

   use interface_module          ! To read command line options
   use precision_module          ! Set real and complex precision
   use user_module               ! User functions, user defined matrix-vector product
   use pp_idrs_module            ! Preconditioner IDRS solver calls

   implicit none 

   integer                       :: grid
   integer, parameter            :: default_grid = 5
!
   integer, parameter            :: Nsigma = 6              ! Number of sigma
   real(kind=rp)                 :: sigma(Nsigma)           ! Shifts

   type(user_matrix)             :: A, M

   real(kind=rp), allocatable    :: x_ms(:,:), b(:), x(:)

   real(kind=rp)                 :: rcenter, rfoci(2)
   integer                       :: i

   logical                       :: multishift
   logical, parameter            :: default_multishift = .false.

   integer                       :: my_proc, n_procs

! Read command line
   call initialize()

   my_proc = this_image()
   n_procs = num_images()

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Definition of the problem. This is the part that has to be made by the user
! The routines used here should be put in the user_module, together with
! the definition of the user type user_matrix.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Which grid?
   grid = get_parameter( '-grid', default_grid )

! Multishift problem?
   multishift = get_parameter( '-multi', default_multishift )

! Make the matrix
   call cdr_matrix( A, grid )

! Output message, problem description:
   if ( my_proc == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      if ( multishift ) then
         write(*,'(a,/)') &
            'Multishift convection-diffusion-reaction problem.'
      else
         write(*,'(a,/)') &
            'Convection-diffusion-reaction problem.'
      end if
   end if

! Define the shifts:
   do i = 1, Nsigma
      sigma(i) = real(i-1,kind=rp)
   end do

! Allocate solution and right-hand-side vectors:
   allocate( b(A%Nn), x(A%Nn) )
   if ( multishift ) allocate( x_ms(A%Nn,Nsigma) )

! Compute right-hand-side
   call model_solution( x, A )
   b = A*x

! Parameters for polynomial preconditioner:
   rcenter = A%c
   rfoci(1) = 0.15_rp*A%c
   rfoci(2) = 1.85_rp*A%c

   if ( multishift ) then
      x_ms = pp_idrs( A, b, sigma )
   else
      if ( grid > 4 ) then
         x    = pp_idrs( A, b, foci=rfoci, center=rcenter )
      else
         x    = pp_idrs( A, b )
      end if
   end if

   if ( my_proc == 1 ) &
      write(*,'(a)') '========================================================================='

end program cdr

program stencil

!
! This programme illustrates the use of the idrs-module on
! a convection-diffusion-reaction problem with constant coefficients
! and homogeneous dirichlet boundary conditions.
! The problem can be  solved with the different algorithms: IDRS, BiCGSTAB and QMRIDR,
! and Neumann and Chebyshev polynomial preconditioners.
!
! The matrix is not explicitly stored. The matrix-vector multiplication
! is done using a stencil operations.
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2023 Martin van Gijzen
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
!
   integer                       :: grid                  
   integer, parameter            :: default_grid = 3      
!
   type(user_matrix)             :: A1

   real(kind=rp), allocatable    :: x(:), x0(:), b(:), U0(:,:), H(:,:), omega(:), D(:)

   real(kind=rp)                 :: radius

   integer                       :: tb, te, clock_rate, clock_max  ! For system_clock
   integer                       :: my_proc, n_procs ! For paralllel computing

! Declaration and initialization of IDRS parameters:
   include "../INCLUDE/initialize.inc"

   my_proc = this_image()
   n_procs = num_images()

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Definition of the problem. This is the part that has to be made by the user
! The routines used here should be put in the user_module, together with
! the definition of the user type user_matrix.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Which grid?
   grid = get_parameter( '-grid', default_grid )

   if ( my_proc == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      write(*,'(6x,a,/)') 'Convection-diffusion-reaction problem'
   end if

! Make the convection-diffusion-reaction matrix:
   call cdr_matrix( A1, grid )

! Allocate solution and right-hand-side vectors:
   allocate( x(A1%Nn), x0(A1%Nn), b(A1%Nn) )

! Compute right-hand-side
   call model_solution( x0, A1 )
   b = A1*x0
   x0 = 0.

! Store in the matrix type needed by IDRS
   A = real_matrix( A1 )

! Compute estimates for largest and smallest eigenvalues with Gershgorin's theorem:
   call gershgorin( A1, rcenter, radius )
   rfoci(1) = rcenter-radius   
   rfoci(2) = rcenter+radius   

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! User part is ready. The rest is only standard IDRS-test
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Inital solve to determine iteration parameters, set-up preconditioner:
   include "../INCLUDE/real_preconditioner.inc"

! Output message, which methods:
   include "../INCLUDE/output_methods.inc"

! Solve the system:
   include "../INCLUDE/run_idrs.inc"
!
   if ( this_image() == 1 ) &
      write(*,'(a)') '========================================================================='
!

end program stencil


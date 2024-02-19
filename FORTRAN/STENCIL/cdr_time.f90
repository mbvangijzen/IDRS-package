program stencil_time

!
! This programme illustrates the use of the idrs-module on
! a time dependent convection-diffusion-reaction problem with constant coefficients
! and homogeneous dirichlet boundary conditions.
! The problem can be  solved with the different algorithms: IDRS, BiCGSTAB and QMRIDR,
! Neumann and Chebyshev polynomial preconditioners, and with and without recycling techniques.
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
   real(kind=rp)                 :: reaction                ! Shift
   real, parameter               :: default_reaction = 5.   ! 
!
   type(user_matrix)             :: A1, At

   real(kind=rp), allocatable    :: x(:), x_mod(:), x0(:), b(:), U0(:,:), omega(:), rhs(:), D(:), H(:,:)

   integer, parameter            :: nt = 10
   real(kind=rp)                 :: dt = 1.

   integer                       :: tb, te, t1, t0, clock_rate, clock_max  ! For system_clock
   integer                       :: my_proc, n_procs

   real(kind=rp)                 :: radius
   integer                       :: it, total_iter


! Declaration and initialization of IDRS parameters:
   include "../INCLUDE/initialize.inc"

   my_proc = this_image()
   n_procs = num_images()

! Which grid?
   grid = get_parameter( '-grid', default_grid )

! Make the matrix
   call cdr_matrix( A1, grid )

   if ( my_proc == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      write(*,'(6x,a,i10,/)') &
         'Convection-diffusion-reaction problem, number of equations is ', n_procs*A1%Nn
   end if

! Make the time-integration matrix:
   At = A1
   At%c = 1. + dt*A1%c
   At%e = dt*A1%e
   At%w = dt*A1%w
   At%n = dt*A1%n
   At%s = dt*A1%s
   At%f = dt*A1%f
   At%b = dt*A1%b

   A = real_matrix( At )

! Allocate solution and right-hand-side vectors:
   allocate( x(A1%Nn), x_mod(A1%Nn), x0(A1%Nn), b(A1%Nn), rhs(A1%Nn), omega(n_omega) )
   x0 = 0.

! Compute right-hand-side
   call model_solution( x_mod, A1 )
   rhs = A1*x_mod
   b = x0+dt*rhs

! Compute estimates for largest and smallest eigenvalues with Gershgorin's theorem:
   call gershgorin( At, rcenter, radius )
   rfoci(1) = rcenter-radius   
   rfoci(2) = rcenter+radius   
   include "../INCLUDE/real_preconditioner.inc"

! Inital solve to determine to recycling subspace
   include "../INCLUDE/recycle.inc"

! Output message, which methods:
   include "../INCLUDE/output_methods.inc"

! Time integration
   total_iter = 0

! Inital solve to determine iteration parameters:
   call system_clock ( t0, clock_rate, clock_max )

   do it = 1, nt
      include "../INCLUDE/run_idrs.inc"
      x0 = x
      b = x0+dt*rhs
      total_iter = total_iter + iter
   end do
   call system_clock ( t1, clock_rate, clock_max )

   if ( my_proc == 1 ) then
      write(*,*)
      write(*,'(a,f8.2)')  'Total elapsed time            = ', real ( t1 - t0 ) / real ( clock_rate )
      write(*,'(a,i4)')    'Total number of iterations    = ', total_iter
      write(*,*)
      write(*,'(a)') &
         '==============================================================================='
   end if

end program stencil_time


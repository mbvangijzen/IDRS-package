program integer_random

!
! This programme illustrates the use of the idrs-module on an academic integer dense system. 
! The rhs-vector is such that the solution is equal to one.
!
! The matrix is row-partitioned for parallelisation.
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2023 Martin van Gijzen
!

   use precision_module ! Set the precisions for real and complex arithmatic

   use matrix_module    ! Always needed is a matrix module that defines the matrix type,
                        ! the matrix-vector multiplication and the preconditioning operation
                        ! Here we use the matrix module supplied with the package
                        ! The matrix module needs a user module that defines the 
                        ! user_matrix type and the multiplication with a user_matrix

   use user_module      ! This one has to be supplied by the user

   use idrs_module      ! The idrs solvers
   use ritz_module      ! To estimate iteration parameters

   use interface_module ! To read command line parameters

   implicit none 
  
   integer, parameter            :: default_nodes = 100
   integer                       :: N                ! Number of equations
   integer                       :: nrows            ! Number of rows per processor
   integer                       :: ncols            ! Number of columns

   real(kind=rp), allocatable    :: x(:), b(:), x0(:), omega(:), U0(:,:), H(:,:), D(:)

   integer                       :: tb, te, clock_rate, clock_max  ! For system_clock
   integer                       :: my_proc, n_procs

   real(kind=rp), allocatable    ::  R(:,:)
   integer, allocatable          :: IR(:,:)
   integer, parameter            :: maxint = 10

   integer, allocatable          :: seed(:) 
   integer                       :: n_seed

! Declaration and initialization of IDRS parameters:
   include "../INCLUDE/initialize.inc"

   my_proc = this_image()
   n_procs = num_images()

! Set the seed of the random-number generator to a fixed value for reproducability
   call random_seed(size=n_seed)
   allocate(seed(n_seed))
   seed = this_image()
   call random_seed(put=seed)

! How many nodes?
   N = get_parameter('-nodes', default_nodes)
   nrows = ceiling(real(N)/real(n_procs))
   ncols = nrows*n_procs

! Output message, problem description:
   if ( my_proc == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      write(*,'(5x,a,i8,/)') &
         'Randomly generated problem with integer matrix, size is ', ncols
   end if

! Generate the random matrix in dense format:
   allocate(R(nrows,ncols),IR(nrows,ncols))
   call random_number(R)
   IR = floor( R*maxint) + 1 

! Allocate solution and rhs-vectors
   allocate(x(nrows),b(nrows),x0(nrows))
   x0 = 0.

! Store the matrix in matrix-type (dense format)
   A = real_matrix( dense_format( imat = IR ) )

! Right-hand-side vector
   x = 1.
   b = A*x

! These estimates are based on "approximate" Gershgorin bounds:
   rcenter  = 5.5_rp
   rfoci(1) = 5.5_rp-(N-1)*5.5_rp
   rfoci(2) = 5.5_rp+(N-1)*5.5_rp
   rseed    = 5.5_rp-(N-1)*5.5_rp
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

end program integer_random



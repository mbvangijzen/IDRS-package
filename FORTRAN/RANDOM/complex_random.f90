program complex_random

!
! This programme illustrates the use of the idrs-module on an academic complex dense system. 
! The rhs-vector is such that the solution is equal to one.
!
! The matrix is row-partitioned for parallelisation.
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2024 Martin van Gijzen
!

   use precision_module ! Set the precisions for real and complex arithmatic

   use matrix_module    ! Always needed is a matrix module that defines the matrix type,
                        ! the matrix-vector multiplication and the preconditioning operation
                        ! Here we use the matrix module supplied with the package
                        ! The matrix module needs a user module that defines the 
                        ! user_matrix type and the multiplication with a user_matrix

   use user_module      ! This one has to be supplied by the user

   use idrs_module      ! The idrs solvers

   use interface_module ! To read command line parameters
   use ritz_module      ! To determine iteration parameters

   implicit none 
  
   integer, parameter            :: default_nodes = 100
   integer                       :: N                ! Number of equations
   integer                       :: nrows            ! Number of rows per processor
   integer                       :: ncols            ! Number of columns

   complex(kind=cp), allocatable :: x(:), b(:), x0(:), U0(:,:), omega(:), H(:,:), D(:)
   real(kind=rp)                 :: err
   integer                       :: i, j

   integer                       :: tb, te, clock_rate, clock_max  ! For system_clock
   integer                       :: my_proc, n_procs

   complex(kind=cp), allocatable :: R(:,:)
   real(kind=rp), allocatable    :: Rr(:,:), Ri(:,:)

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
         'Randomly generated problem with complex matrix, size is ', ncols
   end if

! Generate the random matrix in dense format:
   allocate(Rr(nrows,ncols),Ri(nrows,ncols),R(nrows,ncols))
   call random_number(Rr)
   call random_number(Ri)
   R = cmplx(Rr,Ri)

! Allocate solution and rhs-vectors
   allocate(x(nrows),b(nrows),x0(nrows))
   x0 = 0.

! Store the matrix in matrix-type (dense format)
   A = complex_matrix( dense_format( cmat = R ) )

! Right-hand-side vector
   x = 1.
   b = A*x

! These estimates are based on "approximate" Gershgorin bounds:
   ccenter  = (0.5_cp,0.5_cp)
   cfoci(1) = cmplx(0.5_cp-(N-1)*0.5_cp,0.5_cp-(N-1)*0.5_cp)
   cfoci(2) = cmplx(0.5_cp+(N-1)*0.5_cp,0.5_cp+(N-1)*0.5_cp)
   cseed    = cmplx(0.5_cp-(N-1)*0.5_cp,0.5_cp-(N-1)*0.5_cp)
! Inital solve to determine iteration parameters, set-up preconditioner:
   include "../INCLUDE/complex_preconditioner.inc"

! Output message, which methods:
   include "../INCLUDE/output_methods.inc"

! Solve the system:
   include "../INCLUDE/run_idrs.inc"
!
   if ( this_image() == 1 ) &
      write(*,'(a)') '========================================================================='
!

end program complex_random



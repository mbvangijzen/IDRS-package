program dense_integer

!
! This programme illustrates the use of the idrs-module on an academic integer dense system. 
! The rhs-vector is such that the solution is equal to one.
!
! The matrix is row-partitioned for parallelisation.
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2025 Martin van Gijzen
!

   use precision_module ! Set the precisions for real and complex arithmatic
   use matrix_module
   use user_module      ! This one has to be supplied by the user
   use pp_idrs_module   ! To solve the system

   implicit none 
  
   integer                       :: ncols = 2**12 
   integer                       :: nrows            ! Number of rows per processor

   type(user_matrix)             :: A, P
   real(kind=rp), allocatable    :: x(:), b(:)

   integer                       :: my_proc, n_procs

   real(kind=rp), allocatable    ::  R(:,:)
   integer, allocatable          :: IR(:,:)
   integer, parameter            :: maxint = 100

   integer, allocatable          :: i_seed(:) 
   integer                       :: n_seed

   my_proc = this_image()
   n_procs = num_images()

! Set the seed of the random-number generator to a fixed value for reproducability
   call random_seed(size=n_seed)
   allocate(i_seed(n_seed))
   i_seed = this_image()
   call random_seed(put=i_seed)

! How many nodes?
   nrows = ncols/n_procs

! Output message, problem description:
   if ( my_proc == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      write(*,'(a,/)') &
         'Randomly generated problem with integer matrix.'
   end if

! Generate the random matrix in dense format:
   allocate(R(nrows,ncols),IR(nrows,ncols))
   call random_number(R)
   IR = floor( R*maxint) + 1 

! Allocate solution and rhs-vectors
   allocate(x(nrows),b(nrows))

! Store the matrix in matrix-type (dense format)
   A = dense_format( imat = IR )

! Right-hand-side vector
   x = 1.
   b = A*x

   P = real_precon( 2, A )
   x = pp_idrs( A, b, P=P )
!
   if ( this_image() == 1 ) &
      write(*,'(a)') '========================================================================='
!

end program dense_integer


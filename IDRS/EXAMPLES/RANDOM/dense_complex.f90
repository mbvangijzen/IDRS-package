program dense_complex

!
! This programme illustrates the use of the idrs-module on an academic complex dense system. 
! The rhs-vector is such that the solution is equal to one.
!
! The matrix is row-partitioned for parallelisation.
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2025 Martin van Gijzen
!

   use precision_module ! Set the precisions for real and complex arithmatic
   use user_module      ! This one has to be supplied by the user
   use matrix_module
   use pp_idrs_module   ! To solve the system

   implicit none 
  
   integer                       :: ncols=2**12      ! Number of columns
   integer                       :: nrows            ! Number of rows per processor

   type(user_matrix)             :: A, P
   complex(kind=cp), allocatable :: x(:), b(:)

   integer                       :: my_proc, n_procs

   complex(kind=cp), allocatable :: R(:,:)
   real(kind=rp), allocatable    :: Rr(:,:), Ri(:,:)

   integer, allocatable          :: i_seed(:)
   integer                       :: n_seed

   my_proc = this_image()
   n_procs = num_images()

! Set the seed of the random-number generator to a fixed value for reproducability
   call random_seed(size=n_seed)
   allocate(i_seed(n_seed))
   i_seed = this_image()
   call random_seed(put=i_seed)

   nrows = ncols/n_procs

! Output message, problem description:
   if ( my_proc == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      write(*,'(a,/)') &
         'Randomly generated problem with complex matrix.'
   end if

! Generate the random matrix in dense format:
   allocate(Rr(nrows,ncols),Ri(nrows,ncols),R(nrows,ncols))
   call random_number(Rr)
   call random_number(Ri)
   R = cmplx(Rr,Ri)

! Allocate solution and rhs-vectors
   allocate(x(nrows),b(nrows))

! Store the matrix in matrix-type (dense format)
   A = dense_format( cmat = R )

! Right-hand-side vector
   x = 1.
   b = A*x

   P = complex_precon( 2, A )
   x = pp_idrs( A, b, P=P )
!
   if ( this_image() == 1 ) &
      write(*,'(a)') '========================================================================='
!

end program dense_complex


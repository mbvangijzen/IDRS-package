program band_complex

!
! This programme illustrates the use of the idrs-module on an academic complex random band system. 
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
  
   integer                       :: ncols= 2**16      ! Number of columns
   integer                       :: nrows             ! Number of rows per processor
   integer                       :: half_band = 2**8  ! Half bandwidth

   type(user_matrix)             :: A, P
   complex(kind=cp), allocatable :: x(:), b(:)

   integer                       :: my_proc, n_procs

   real(kind=rp), allocatable    :: Rr(:,:), Ri(:,:)
   complex(kind=cp), allocatable :: R(:,:)

   integer, allocatable          :: i_seed(:) 
   integer                       :: n_seed
   integer                       :: i, dia, offset

   my_proc = this_image()
   n_procs = num_images()

! Set the seed of the random-number generator to a fixed value for reproducability
   call random_seed(size=n_seed)
   allocate(i_seed(n_seed))
   i_seed = this_image()
   call random_seed(put=i_seed)

! How many local nodes?
   nrows = ncols/n_procs
   ncols = nrows*n_procs

! Output message, problem description:
   if ( my_proc == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      write(*,'(a,/)') &
         'Randomly generated problem with complex band matrix.'
   end if

! Generate the random matrix in band format:
   allocate(R(nrows,-half_band:half_band),Rr(nrows,-half_band:half_band),Ri(nrows,-half_band:half_band))
   call random_number(Rr)
   call random_number(Ri)
   R = cmplx(Rr,Ri,kind=cp)
   offset = (my_proc-1)*nrows
   do i = 1, nrows
      do dia = -half_band, half_band
         if ( (i+offset+dia) < 1     ) then 
            R(i,dia) = 0._cp
         elseif ( (i+offset+dia) > ncols ) then
            R(i,dia) = 0._cp
         end if
      end do
   end do

! Allocate solution and rhs-vectors
   allocate(x(nrows),b(nrows))

! Store as standard matrix:
   A = cds_format( nrows, ncols, half_band, cband=R )

! Right-hand-side vector
   x = 1.
   b = A*x

   P = complex_precon( 3, A )
   x = pp_idrs( A, b, P=P )
!
   if ( this_image() == 1 ) &
      write(*,'(a)') '========================================================================='
!

end program band_complex

program test_integer

!
! This programme illustrates the use of the idrs-module on an academic integer random sparse system.
! The matrix is row-partitioned for parallelisation.
! The storage format can be chosen from DENSE, CRS, COO, and CDS
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2026 Martin van Gijzen
!

   use interface_module          ! To read command line options
   use precision_module          ! Set real and complex precision
   use user_module               ! User functions, user defined matrix-vector product
   use idrs_driver_module        ! Simple call to IDRS

   implicit none 

   integer, parameter            :: default_nodes = 1024         ! Number of nodes
   integer                       :: N
   integer                       :: ncols                        ! Number of columns
   integer                       :: nrows                        ! Number of rows per processor

   integer, parameter            :: default_nsigma = 0           ! Number of shifts
   integer                       :: nsigma

   integer, parameter            :: default_nrhs = 0             ! Number of rhs
   integer                       :: nrhs

   type(user_matrix)             :: A, D
   real(kind=rp), allocatable    :: x(:), x_ms(:,:), x_mr(:,:), x_msmr(:,:,:), b(:), b_mr(:,:)
   integer, allocatable          :: diag(:,:)

   integer                       :: my_proc, n_procs

   real(kind=rp), allocatable    :: sigma(:)
   real(kind=rp)                 :: dsigma

! For random number generator:
   integer, allocatable          :: i_seed(:)
   integer                       :: n_seed
   real(kind=rp), allocatable    :: r(:)

! For generation of the matrix in CRS format 
   integer, allocatable          :: col_ind(:), row_ptr(:)
   integer, allocatable          :: ival(:)
   integer                       :: nnz_row

! For generation of the matrix in COO format
   integer                       :: nnz
   integer, allocatable          :: indx(:), jndx(:)

! For generation of the matrix in DENSE or CDS format
   integer                       :: half_band
   integer, allocatable          :: P(:,:)

! Other variables:
   real(kind=rp), parameter      :: seed = -1._rp
   real(kind=rp), parameter      :: shift = 2._rp*seed
   integer, parameter            :: maxint = 1000
   integer                       :: i, j, k, i_sigma, dia
   real(kind=rp)                 :: som

   call initialize()

   my_proc = this_image()
   n_procs = num_images()

! How many nodes?
   N = get_parameter( '-nodes', default_nodes )
   nrows = ceiling(real(N)/real(n_procs))
   ncols = nrows*n_procs

! Multishift analysis?
   nsigma = get_parameter( '-multi', default_nsigma )

! Multiple right-hand-sides?
   nrhs = get_parameter( '-rhs', default_nrhs )

! Set the seed of the random-number generator to a fixed value for reproducability
   call random_seed(size=n_seed)
   allocate(i_seed(n_seed))
   i_seed = this_image()
   call random_seed(put=i_seed)

! Output message, problem description:
   if ( my_proc == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      write(*,'(a)') 'Test if all four formats work correctly. A random test matrix is generated,'
      write(*,'(a)') 'subsequently converted to all four formats and then solved.'
      write(*,'(a)') 'The number of rhs-vectors and shifts can be specified as command line options.'
      write(*,*)
      write(*,'(a,i8,/)') &
         'Randomly generated problem with real matrix of size ', ncols
   end if

! First generate the random matrix in CRS-format:

! Number of nonzeroes per row and in total::
   nnz_row = int(sqrt(sqrt(real(ncols))))
   nnz     = nrows*nnz_row

! Allocate space:
   allocate( row_ptr(nrows+1), col_ind(nnz),ival(nnz),r(nnz) )
   call random_number(r)
   ival = int(r*maxint)

! Determine row pointers:
   row_ptr(1) = 1
   do i = 1, nrows
      row_ptr(i+1) = row_ptr(i) + nnz_row
   end do

! Determine column indices:
   call random_number(r)
   col_ind = int(r*ncols) + 1

! Generate the diagonal matrix with row sums:
   allocate(diag(nrows,0:0))
   do i = 1, nrows
      diag(i,0) = sum(ival(row_ptr(i):row_ptr(i+1)-1))
   end do
   D = cds_format( nrows, ncols, 0, iband=diag )

! Now convert to one of the other formats, store in user_matrix format, and solve the systems

! COO:
   if ( my_proc == 1 ) then
      write(*,'(a)') '-------------------------------------------------------------------------'
      write(*,'(a)') 'Format is COO'
      write(*,*) 
   end if

! Convert from CRS to COO format:
   allocate( indx(nnz), jndx(nnz) )
   k = 0
   do i = 1, nrows
      do j = row_ptr(i),row_ptr(i+1)-1
         k = k+1
         indx(k) = i
         jndx(k) = col_ind(k)
      end do
   end do

! Create the user matrix (real, coo format). 
   A = coo_format( nrows, ncols, nnz, indx, jndx, ival=ival )

! Solve the systems.
   include "real_solve.inc"

! DENSE:
   if ( my_proc == 1 ) then
      write(*,'(a)') '-------------------------------------------------------------------------'
      write(*,'(a)') 'Format is DENSE'
      write(*,*) 
   end if

! Convert from CRS to DENSE format:
   allocate( P(nrows,ncols) )
   P = 0._rp
   k = 0
   do i = 1, nrows
      do j = row_ptr(i),row_ptr(i+1)-1
         k = k+1
         P(i,col_ind(k)) = P(i,col_ind(k)) + ival(k)
      end do
   end do

! Create the user matrix (real, dense format). 
   A = dense_format( imat = P )
   deallocate(P)

! Solve the systems.
   include "real_solve.inc"

! CDS:
   if ( my_proc == 1 ) then
      write(*,'(a)') '-------------------------------------------------------------------------'
      write(*,'(a)') 'Format is CDS'
      write(*,*) 
   end if

! Convert from CRS to CDS format:
   half_band = ncols-1  ! this is very pessimistic
   allocate( P(nrows,-half_band:half_band))
   P = 0._rp
   k = 0
   do i = 1, nrows
      do j = row_ptr(i),row_ptr(i+1)-1
         k = k+1
         dia = col_ind(k) - (i+nrows*(my_proc-1))
         P(i,dia) = P(i,dia) + ival(k)
      end do
   end do

! Create the user matrix:
   A = cds_format( nrows, ncols, half_band, iband=P )

! Solve the systems.
   include "real_solve.inc"

! CRS:
   if ( my_proc == 1 ) then
      write(*,'(a)') '-------------------------------------------------------------------------'
      write(*,'(a)') 'Format is CRS'
      write(*,*) 
   end if

! Create the user matrix (real, dense format). 
   A = crs_format( nrows, ncols, row_ptr, col_ind, ival=ival ) 

! Solve the systems.
   include "real_solve.inc"

!
   if ( this_image() == 1 ) &
      write(*,'(a)') '========================================================================='

end program test_integer

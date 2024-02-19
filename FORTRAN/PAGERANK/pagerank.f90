program pagerank

!
! This programme illustrates the use of the idrs-module. A random PageRank problem
! is generated and solved with one of the different idrs-algorithms:
! IDRS and QMRIDR. The conventional value of 0.85 is used for the chance of 
! following an outlink.
!
! The PageRank matrix is stored in CRS format.
!
! The PageRank linear systems are quite easy to solve and therefore not very suited
! to illustrate the differences in performance of the variants: all methods
! converge in a number of iterations that is close to optimal (= the number of
! full GMRES iterations).
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

   use ritz_module      ! To estimate iteration parameters

   implicit none 

   integer, parameter            :: default_nodes = 10000
   integer                       :: N                ! Number of nodes
   real(kind=rp), parameter      :: p = 0.85         ! Chance of following outlink

   real(kind=rp), allocatable    :: x(:), x0(:), b(:), U0(:,:), omega(:), H(:,:), D(:)
   real(kind=rp)                 :: som
   integer                       :: i, j, k

   integer                       :: tb, te, clock_rate, clock_max  ! For system_clock
   integer                       :: my_proc, n_procs

   integer, parameter            :: max_links = 10
   integer                       :: nnz, nrows
   integer, allocatable          :: col_ind(:), row_ptr(:), nnz_row(:)
   real(kind=rp), allocatable    :: rval(:)
   real(kind=rp), allocatable    :: r(:), outlinks(:)

   include "../INCLUDE/initialize.inc"

   my_proc = this_image()
   n_procs = num_images()

! How many nodes?
   N = get_parameter('-nodes', default_nodes)
   nrows = ceiling(real(N)/real(n_procs))
   N = nrows*n_procs

! Output message, problem description:
   if ( my_proc == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      write(*,'(5x,a,i10,/)') &
         'Randomly generated PageRank problem, size is ', N
   end if

! Generate the Google matrix in CRS-format:

! Determine nonzeroes per row and allocate space:
   allocate(r(max(nrows,max_links+1)),nnz_row(nrows))
   call random_number(r)
   nnz_row = int(r(1:nrows)*(max_links) ) + 2
   nnz = sum( nnz_row )
   allocate(row_ptr(nrows+1),col_ind(nnz),rval(nnz))

! Determine row pointers:
   row_ptr(1) = 1
   do i = 1, nrows
      row_ptr(i+1) = row_ptr(i) + nnz_row(i)
   end do

! Determine column indices:
   do i = 1, nrows
      call random_number(r(2:nnz_row(i)))
      col_ind(row_ptr(i)) = i + (my_proc - 1)*nrows              ! Main diagonal element
      col_ind(row_ptr(i)+1:row_ptr(i+1)-1) = int(r(2:nnz_row(i))*N) + 1
   end do

! Determine number of outlinks per node (nonzeroes per column)
   allocate(outlinks(N))
   outlinks = 0.
   do i = 1, nnz
      outlinks(col_ind(i)) = outlinks(col_ind(i)) + 1._rp
   end do
   call co_sum( outlinks )
   outlinks = outlinks -1. ! correct for the ones on the main diagonal

! compute the elements of A:
   do i = 1, nrows
      rval(row_ptr(i)) = 1.
      do j = row_ptr(i)+1,row_ptr(i+1)-1
         rval(j) = -p/outlinks(col_ind(j))
      end do
   end do

! Allocate solution and rhs-vectors
   allocate(x(nrows),b(nrows),x0(nrows))
   x0 = 0.

! Right-hand-side vector
   b = 1. 

! make the user matrix (CRS format)
   A = real_matrix( crs_format( nrows, row_ptr, col_ind, rval ) )

! Polynomial preconditioner
   rfoci(1) = p
   rfoci(2) = 2.-p
   rcenter  = 1.
   include "../INCLUDE/real_preconditioner.inc"

! Output message, which methods:
   include "../INCLUDE/output_methods.inc"

! Solve the system:
   include "../INCLUDE/run_idrs.inc"

! Make the solution vectors stochastic:
   som = sum(x)
   call co_sum(som)
   x = x/som

!
   if ( this_image() == 1 ) &
      write(*,'(a)') '========================================================================='
!

end program pagerank



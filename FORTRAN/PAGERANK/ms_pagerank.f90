program ms_pagerank

!
! This programme illustrates the use of the idrs-module. A random PageRank problem 
! is generated and solved with the different multishift idrs-algorithms: 
! MSIDRS and MSQMRIDR. Here, only one shift is used, 1/0.85, which corresponds 
! to the conventional value of 0.85 for the chance of follwoing an outlink. 
! It is easy to simultaneously solve for different teleportation parameters by 
! adding more shifts.
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
  
   integer, parameter            :: default_nodes = 10000        ! Number of nodes
   real(kind=rp), parameter      :: p = 0.85                     ! Chance of following outlink 
   integer, parameter            :: N_sigma = 1                  ! Number of damping factor
   integer                       :: N
   real(kind=rp), parameter      :: shift(N_sigma) = [0.85]
   real(kind=rp)                 :: sigma(N_sigma), sigma_p(N_sigma)

   real(kind=rp), allocatable    :: b(:), x(:,:), x0(:), U0(:,:), omega(:), H(:,:), D(:)
   integer                       :: i

   integer, parameter            :: max_links = 10
   integer                       :: nnz, nrows, dangling
   integer, allocatable          :: col_ind(:), row_ptr(:), nnz_row(:)
   real(kind=rp), allocatable    :: r(:), outlinks(:)

   integer                       :: tb, te, clock_rate, clock_max  ! For system_clock
   integer                       :: my_proc, n_procs

   real(kind=rp)                 :: som

   include "../INCLUDE/initialize.inc"

   my_proc = this_image()
   n_procs = num_images()

! How many nodes?
   N = get_parameter( '-nodes', default_nodes )
   nrows = ceiling(real(N)/real(n_procs))
   N = nrows*n_procs

! Output message, problem description:
   if ( my_proc == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      write(*,'(a,i8,a,i2,/)') &
         'Randomly generated multishift PageRank problem, size is ', N, '. Number of shifts is ', N_sigma
   end if

! Generate the connectivity matrix in CRS-format:

! Determine nonzeroes per row and allocate space:
   allocate(r(nrows), nnz_row(nrows))
   call random_number(r)
   nnz_row = int(r(1:nrows)*max_links) + 1
   nnz = sum( nnz_row )
   allocate(row_ptr(nrows+1),col_ind(nnz))

! Determine row pointers:
   row_ptr(1) = 1
   do i = 1, nrows
      row_ptr(i+1) = row_ptr(i) + nnz_row(i)
   end do

! Determine column indices:
   do i = 1, nrows
      call random_number(r(1:nnz_row(i)))
      col_ind(row_ptr(i):row_ptr(i+1)-1) = int(r(1:nnz_row(i))*N) + 1
   end do

! Determine number of outlinks per node (nonzeroes per column)
   allocate(D(nrows))
   allocate(outlinks(N))
   outlinks = 0.
   do i = 1, nnz
      outlinks(col_ind(i)) = outlinks(col_ind(i)) + 1._rp
   end do
   call co_sum( outlinks )

! Put a 1 in outlinks if there are no outlinks (complete column is 0)
   dangling = 0
   do i = 1, N
      if ( outlinks(i) == 0 ) then
         outlinks(i) = 1._rp
         dangling = dangling + 1
      end if
   end do
   D = outlinks((my_proc-1)*nrows+1:my_proc*nrows)
   deallocate(outlinks)

! Rewrite eigenvalue problem as linear system:
! (pG/D + pevT + (1-p)ee^T)x = x
! (pG/Dx + ((1-p)+v^Tx) e = x
! (pG/D -I)x = beta e (beta is just scaling)
! (G/D - (1/p))x = beta e

! Create the user matrix (real, csr format). 
   A = real_matrix( crs_format( nrows, row_ptr, col_ind ) )
   
! Space for rhs and solution
   allocate( x(nrows,N_sigma), x0(nrows), b(nrows), resnrm(N_sigma) )
   x0 = 0.

! Right-hand-side vector
   b = 1. 

! Shifts
   do i = 1, N_sigma
      sigma(i) = 1./shift(i)
   end do

! Polynomial preconditioner
   rfoci(1) = -1.
   rfoci(2) =  1.
   rcenter  = 0.
   rseed    = -1./p

! Create the preconditioned system matrix:
   include "../INCLUDE/real_ms_preconditioner.inc"

! Output message, which methods:
   include "../INCLUDE/output_methods.inc"

! Solve the system:
   include "../INCLUDE/run_ms_idrs.inc"

!
   if ( this_image() == 1 ) &
      write(*,'(a)') '========================================================================='
!

end program ms_pagerank



program pagerank

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

   use interface_module          ! To read command line options
   use precision_module          ! Set real and complex precision
   use user_module               ! User functions, user defined matrix-vector product
   use pp_idrs_module            ! Preconditioner IDRS solver calls

   implicit none 
  
   integer, parameter            :: default_nodes = 1024         ! Number of nodes
   integer                       :: N

   integer, parameter            :: N_sigma = 2                  ! Number of damping factor
   real(kind=rp), parameter      :: shift(N_sigma) = [0.85_rp, 0.95_rp]
   real(kind=rp)                 :: sigma(N_sigma)

   type(user_matrix)             :: A, P

   real(kind=rp), allocatable    :: x(:,:), b(:)

   real(kind=rp)                 :: som
   integer                       :: i, j, i_sigma

   integer, parameter            :: max_links = 10
   integer                       :: nnz, nrows, ncols
   integer, allocatable          :: col_ind(:), row_ptr(:), nnz_row(:)
   real(kind=rp), allocatable    :: rval(:)
   real(kind=rp), allocatable    :: r(:), outlinks(:)

   logical                       :: multishift
   logical, parameter            :: default_multishift = .false.

   logical                       :: lu
   logical, parameter            :: default_lu = .false.

   real(kind=rp)                 :: rcenter, rfoci(2), rseed, rshift

   integer                       :: my_proc, n_procs

   call initialize()

   my_proc = this_image()
   n_procs = num_images()

! How many nodes?
   N = get_parameter( '-nodes', default_nodes )
   nrows = ceiling(real(N)/real(n_procs))
   ncols = nrows*n_procs

! Multishift problem?
   multishift = get_parameter( '-multi', default_multishift )

! LU decomposition?
   lu = get_parameter( '-lu', default_lu )

! Output message, problem description:
   if ( my_proc == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      if ( multishift ) then
         write(*,'(a,/)') 'Randomly generated multishift PageRank problem.'
      else
         write(*,'(a,/)') 'Randomly generated PageRank problem.'
      end if
   end if

! Generate the scaled connectivity matrix in CRS-format:

! Determine nonzeroes per row and allocate space:
   allocate(r(max(nrows,max_links+1)), nnz_row(nrows))
   call random_number(r)
   nnz_row = int(r(1:nrows)*max_links) + 1
   nnz = sum( nnz_row )
   allocate(row_ptr(nrows+1),col_ind(nnz),rval(nnz))

! Determine row pointers:
   row_ptr(1) = 1
   do i = 1, nrows
      row_ptr(i+1) = row_ptr(i) + nnz_row(i)
   end do

! Determine column indices:
   do i = 1, nrows
      call random_number(r(1:nnz_row(i)))
      col_ind(row_ptr(i):row_ptr(i+1)-1) = int(r(1:nnz_row(i))*ncols) + 1
   end do

! Determine number of outlinks per node (nonzeroes per column)
   allocate(outlinks(ncols))
   outlinks = 0.
   do i = 1, nnz
      outlinks(col_ind(i)) = outlinks(col_ind(i)) + 1._rp
   end do
   call co_sum( outlinks )

! Put a 1 in outlinks if there are no outlinks (complete column is 0)
   do i = 1, ncols
      if ( outlinks(i) == 0 ) then
         outlinks(i) = 1._rp
      end if
   end do

! Rewrite eigenvalue problem as linear system:
! (pG/D + pevT + (1-p)ee^T)x = x
! (pG/Dx + ((1-p)+v^Tx) e = x
! (pG/D -I)x = beta e (beta is just scaling)
! (G/D - (1/p))x = beta e

! compute the elements of A:
   do i = 1, nrows
      do j = row_ptr(i),row_ptr(i+1)-1
         rval(j) = 1._rp/outlinks(col_ind(j))
      end do
   end do

! Create the user matrix (real, csr format). 
   A = crs_format( nrows, ncols, row_ptr, col_ind, rval )
   
! Space for rhs and solution
   allocate( x(nrows,N_sigma), b(nrows) )

! Right-hand-side vector
   b = 1. 

! Shifts
   do i = 1, N_sigma
      sigma(i) = 1./shift(i)
   end do

! Polynomial preconditioner
   rfoci(1) = -1._rp+shift(N_sigma)
   rfoci(2) =  1._rp-shift(N_sigma)
   rcenter  = 0. 
   rseed    = sigma(N_sigma)
   rshift   = sum(sigma)/N_sigma

! Shift-and-invert preconditioner
   if ( lu ) P = real_precon( 2, A, shift = rshift )

! Solve the system:
   if ( multishift ) then
      if ( lu ) then
         x = pp_idrs( A, b, sigma, P=P, seed=rseed, shift=rshift )
      else
         x = pp_idrs( A, b, sigma, seed=rseed, center=rcenter, foci=rfoci )
      end if
   else
      if ( lu ) then
         do i_sigma = 1, N_sigma
            rseed = sigma(i_sigma)
            x(:,i_sigma) = pp_idrs( A, b, P=P, seed=rseed )
         end do
      else
         do i_sigma = 1, N_sigma
            rseed = sigma(i_sigma)
            x(:,i_sigma) = pp_idrs( A, b, seed=rseed, center=rcenter, foci=rfoci )
         end do
      end if
   end if

! Correct for scaling and make the solution vectors stochastic:
   do i_sigma = 1, N_sigma
      som = sum(x(:,i_sigma))
      call co_sum(som)
      x(:,i_sigma) = x(:,i_sigma)/som
   end do
!
   if ( this_image() == 1 ) &
      write(*,'(a)') '========================================================================='

end program pagerank



program pagerank

!
! This programme illustrates the use of the idrs-module. A random PageRank problem 
! is generated and solved with the different multishift idrs-algorithms: 
! MSIDRS and MSQMRIDR. Here, only one shift is used, 1/0.85, which corresponds 
! to the conventional value of 0.85 for the chance of follwing an outlink. 
! It is easy to simultaneously solve for different teleportation parameters by 
! adding more shifts.
!
! The PageRank matrix is stored in DENSE format. 
!
! The PageRank linear systems are quite easy to solve and therefore not very suited
! to illustrate the differences in performance of the variants: all methods
! converge in a number of iterations that is close to optimal (= the number of
! full GMRES iterations).
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2026 Martin van Gijzen
!

   use interface_module          ! To read command line options
   use precision_module          ! Set real and complex precision
   use mm_module                 ! To read the matrix-market matrices
   use user_module               ! User functions, user defined matrix-vector product
   use idrs_driver_module        ! Simple call to IDRS

   implicit none 

! Path to data files
   character(len=19), parameter    :: path = 'DATA/PAGERANK_DATA/'
   integer, parameter              :: default_graph = 1
   integer                         :: graph
   character(len=50)               :: graph_file, website_file
   character(len=160), allocatable :: website_name(:)
  
   integer, parameter              :: N_sigma = 4                  ! Number of damping factors
   real(kind=rp), parameter        :: shift(N_sigma) = [0.85_rp, 0.9_rp, 0.95_rp, 0.99_rp]
   real(kind=rp)                   :: sigma(N_sigma)

   type(user_matrix)               :: G, A

   real(kind=rp), allocatable      :: x(:,:)[:], b(:), K(:,:), eigvec(:,:)

   real(kind=rp)                   :: som
   integer                         :: i, j, i_sigma

   integer, parameter              :: top_pages = 5
   real(kind=rp), allocatable      :: pr(:)
   integer                         :: lc
   real(kind=rp)                   :: pr_top(top_pages,N_sigma)
   integer                         :: lc_top(top_pages,N_sigma)

   integer                         :: nrows, ncols
   logical(kind=rp), allocatable   :: mask(:)
   real(kind=rp), allocatable      :: outlinks(:)

   logical                         :: multishift
   logical, parameter              :: default_multishift = .false.

   logical                         :: eigen
   logical, parameter              :: default_eigen = .false.
   integer, parameter              :: neig = 1
   real(kind=rp)                   :: lambda(neig)

   real(kind=rp)                   :: rseed, rshift, p

   integer                         :: my_proc, n_procs, id

   call initialize()

   my_proc = this_image()
   n_procs = num_images()

   graph = get_parameter( '-graph', default_graph )

! Eigenvalue problem?
   eigen = get_parameter( '-eigen', default_eigen )

! Multishift problem?
   multishift = get_parameter( '-multi', default_multishift )

! Full filename to the graph:
   select case ( graph )
      case ( 1 ) 
         graph_file   = path//'homepage500.mtx'
         website_file = path//'homepage500.txt'
      case ( 2 ) 
         graph_file   = path//'homepage1000.mtx'
         website_file = path//'homepage1000.txt'
      case ( 3 ) 
         graph_file   = path//'homepage2000.mtx'
         website_file = path//'homepage2000.txt'
      case default
         error stop 'Graph number should be 1, 2, or 3 '
   end select

! Read the graph:
   G = mm_matrix( trim(graph_file), 'coo' )
   nrows = G%nrows 
   ncols = G%ncols 
   allocate( K(nrows,ncols), outlinks(ncols), mask(ncols) )
   K = real_dense( G )

! Output message, problem description:
   if ( my_proc == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      if ( eigen ) then
         write(*,'(a,a,/)') 'Eigenvalue PageRank problem ', graph
      elseif ( multishift ) then
         write(*,'(a,a,/)') 'Multishift PageRank problem ', graph
      else
         write(*,'(a,a,/)') 'Linear system PageRank problem ', graph
      end if
   end if

! Remove self_references:
   do i = 1, nrows
      if ( i+(my_proc-1)*nrows <= ncols ) then
         if ( K(i,i+(my_proc-1)*nrows) > 0._rp ) K(i,i+(my_proc-1)*nrows) = 0._rp
      end if
   end do

! Determine number of outlinks per node (nonzeroes per column).
   outlinks = 0.
   do j = 1, ncols
      outlinks(j) = sum(K(:,j))
   end do
   call co_sum( outlinks )

! Dangling nodes:
   do j = 1, ncols
      if ( outlinks(j) == 0 ) then
         K(:,j) = 1._rp
         outlinks(j) = real(ncols,kind=rp)
      end if
   end do

! Rewrite eigenvalue problem as linear system:
! (pG/D + pevT + (1-p)ee^T)x = x
! (pG/Dx + ((1-p)+v^Tx) e = x
! (pG/D -I)x = beta e (beta is just scaling)
! (G/D - (1/p))x = beta e

! compute the elements of A:
   do j = 1, ncols
      K(:,j) = K(:,j)/outlinks(j)
   end do
   
! Space for rhs and solution
   allocate( x(nrows,N_sigma)[*], b(nrows) )

! Right-hand-side vector
   b = 1._rp 

! Shifts
   do i = 1, N_sigma
      sigma(i) = 1._rp/shift(i)
   end do

! Solve the system:
   if ( eigen ) then
!
! Eigenvalue analysis
!
      allocate(eigvec(nrows,neig))
      rshift = 0.99_rp
      do i_sigma = 1, N_sigma
         p = shift(i_sigma)
         A = dense_format( rmat=(p*K+(1-p)/real(ncols,kind=rp)) )
         call idreig_driver( lambda, eigvec, A, b, shift = rshift  )
         x(:,i_sigma) = eigvec(:,1)
      end do

   elseif ( multishift ) then

! Create the user matrix (real, dense format). 
      A = dense_format( rmat=K )

      rseed = sigma(N_sigma)
      x = idrs_driver( A, b, sigma, seed=rseed )

   else

! Create the user matrix (real, dense format). 
      A = dense_format( rmat=K )
      do i_sigma = 1, N_sigma
         rseed = sigma(i_sigma)
         x(:,i_sigma) = idrs_driver( A, b, seed=rseed )
      end do
   end if

! Correct for scaling and make the solution vectors stochastic:
   do i_sigma = 1, N_sigma
      som = sum(x(:,i_sigma))
      call co_sum(som)
      x(:,i_sigma) = x(:,i_sigma)/som
   end do

! Find the top highest scores:
   allocate(pr(N_procs*nrows))
   do i_sigma = 1, N_sigma
      sync all
      do i = 1, N_procs
         pr(1+(i-1)*nrows:i*nrows) = x(:,i_sigma)[i]
      end do
      sync all
      mask = .true.
      do i = 1, top_pages
         lc = maxloc( pr(1:ncols), dim=1, mask=mask )
         lc_top(i,i_sigma) = lc
         pr_top(i,i_sigma) = pr(lc)
         mask(lc) = .false.
      end do
   end do

! Output the results
   if ( my_proc == 1 ) then
      allocate(website_name(ncols))
      open(newunit=id,file=website_file,status='old',action='read')
      do i = 1, ncols
         read(id,'(a)') website_name(i)
      end do

      do i_sigma = 1, N_sigma
         write(*,'(a)') '_________________________________________________________________________'
         write(*,'(a,f6.4)') 'Top ranked pages for damping parameter ', shift(i_sigma)
         write(*,*)
         write(*,'(2x,a,2x,a,2x,a)') 'Node', 'Pagerank','URL'
         do i = 1, top_pages
            j = lc_top(i,i_sigma)
            write(*,'(i6,3x,f6.4,3x,a)') j, pr_top(i,i_sigma), trim(website_name(j))
         end do
      end do
   end if

!
   if ( my_proc == 1 ) &
      write(*,'(a)') '========================================================================='

end program pagerank



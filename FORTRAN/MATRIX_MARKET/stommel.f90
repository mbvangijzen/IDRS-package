program stommel
!
! Wind-driven global ocean circulation model proposed by Stommel.
! The test problem is a second order partial differential equation reduced
! to two coupled second order equations.  The matrix is real and nonsymmetric. 
! Twelve rhs-vectors are supplied, corresponding to time-averaged wind fields 
! for the twelve months. The test problem is supplied on six different grids. 
! Selection of the grid ! is done by specifying as command line option
!    -grid k, with k a number from 1 to 6.
! The grid number k corresponds to the grid resolution in degrees.
!
! Reference:
!   M.B. van Gijzen, C.B. Vreugdenhil, and H. Oksuzoglu, 
!   The Finite Element Discretization for Stream-Function Problems on Multiply Connected Domains. 
!   Journal of Computational Physics, 140:30-46, 1998
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
   use mm_module        ! To read the matrix-market matrices

   implicit none 

! Path to data files
   character(len=22), parameter   :: path = '../../DATA/OCEAN_DATA/'

   integer                        :: t0, t1, tb, te, clock_rate, clock_max

! Test problems:
   integer                        :: grid
   integer, parameter             :: default_grid = 4

   character(len=10), parameter   :: month(12) = &
      ['January   ','February  ','March     ','April     ','May       ',&
       'June      ','July      ','August    ','September ','October   ',&
       'November  ','December  ']
   character(len=50)              :: matrix_file, rhs_file
  
! Matrices, number of equations etc.:
   type(user_matrix)              :: K, RHS
   real(kind=rp), allocatable     :: x_all(:,:), b_all(:,:)
   real(kind=rp), allocatable     :: x(:), b(:), x0(:), U0(:,:), omega(:), H(:,:)
   integer                        :: neq, irhs, nrhs

   real(kind=rp), allocatable     :: D(:)

! Declaration and initialization of IDRS parameters:
   include "../INCLUDE/initialize.inc"
  
! Which test problem?
   grid = get_parameter('-grid', default_grid )

! Output message, problem description:
   if ( this_image() == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      write(*,'(a,i2,/)') 'Test problem Stommel, grid = ', grid
   end if

   write(matrix_file,'(a,a,i1,a)') trim(path),'STOMMEL/stommel', grid,'.mtx'
   write(rhs_file,'(a,a,i1,a)') trim(path),'STOMMEL/stommel', grid,'_b.mtx'

! Read the system matrix:
   K = mm_matrix( trim(matrix_file) )

! Read the right-hand-side matrix:
   RHS = mm_matrix( rhs_file )

   neq = K%nrows
   nrhs = RHS%ncols

! Allocate space for the solution vectors and rhs-vector
   allocate( x(neq), b(neq), x0(neq) )
   x0 = 0.
   allocate( x_all(neq,nrhs), b_all(neq,nrhs) )
! Copy the rhs-vector b from the user-defined type RHS
   b_all = RHS%rmat
 
! Construct the matrix:
   A = real_matrix( K )

! Jacobi preconditioner:
   D = real_diagonal( K )

! Output message, which methods:
   include "../INCLUDE/output_methods.inc"

   call system_clock ( t0, clock_rate, clock_max )

! Inital solve to determine to set-up preconditioner:
   b = b_all(:,1)
   include "../INCLUDE/real_preconditioner.inc"

! Inital solve to determine to recycling subspace
   include "../INCLUDE/recycle.inc"

   do irhs = 1, nrhs
      b = b_all(:,irhs)
      if ( this_image() == 1 ) write(*,'(a,a)') 'Results for ', month(irhs)
! Solve the system
      include "../INCLUDE/run_idrs.inc"
      x0 = x
      x_all(:,irhs) = x
   end do
   call system_clock ( t1, clock_rate, clock_max )

   if ( this_image() == 1 ) then
      write(*,'(a,f8.2,a)') 'Total elapsed time     = ', real(t1-t0)/real(clock_rate), 's.'
      write(*,'(a)') '========================================================================='
   end if

end program stommel

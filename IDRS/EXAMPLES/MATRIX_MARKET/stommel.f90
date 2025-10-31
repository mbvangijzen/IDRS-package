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
   use user_module      ! This one has to be supplied by the user
   use interface_module ! To read command line parameters
   use mm_module        ! To read the matrix-market matrices
   use pp_idrs_module   ! To solve the system

   implicit none 

! Path to data files
   character(len=16), parameter   :: path = 'DATA/OCEAN_DATA/'

! Test problems:
   integer                        :: grid
   integer, parameter             :: default_grid = 4

! Matrices, number of equations etc.:
   type(user_matrix)              :: K, RHS, P, Q
   character(len=50)              :: matrix_file, rhs_file, map_file
   real(kind=rp), allocatable     :: x(:,:), b(:,:)
   integer                        :: nrows, ncols, nrhs, irhs

! For preconditioner:
   integer                        :: which_precon

   logical                        :: lu
   logical, parameter             :: default_lu = .false.

! Plotting:
   character(len=:), allocatable  :: plot
   character(len=4), parameter    :: default_plot = 'none'
   logical                        :: plot_solution
   logical                        :: gif, jpeg
   logical, parameter             :: default_gif  = .false.
   logical, parameter             :: default_jpeg = .false.
   real(kind=rp), allocatable     :: solution(:,:)

   character(len=10), parameter   :: month(12) = &
      ['January   ','February  ','March     ','April     ','May       ',&
        'June      ','July      ','August    ','September ','October   ',&
        'November  ','December  ']

   integer                        :: ngrid, igrid, nb, ne, nx, ny, i, j
   integer                        :: my_proc, num_procs

! Initialization for parameter routines (read command line)
   call initialize( )
   my_proc   = this_image()
   num_procs = num_images()
  
! Which test problem?
   grid = get_parameter('-grid', default_grid )

! Direct solution?
   lu = get_parameter( '-lu', default_lu )

! Output message, problem description:
   if ( my_proc == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      write(*,'(a,i2,/)') 'Test problem Stommel, grid = ', grid
   end if

   write(matrix_file,'(a,a,i1,a)') trim(path),'STOMMEL/stommel', grid,'.mtx'
   write(rhs_file,'(a,a,i1,a)') trim(path),'STOMMEL/stommel', grid,'_b.mtx'

! Read the system matrix:
   K = mm_matrix( trim(matrix_file), 'coo' )

! Read the right-hand-side matrix:
   RHS = mm_matrix( rhs_file )

   nrows = K%nrows
   ncols = K%ncols
   nrhs  = RHS%ncols

! Allocate space for the solution vectors and rhs-vector
   allocate( x(nrows,nrhs), b(nrows,nrhs) )

! Copy the rhs-vector b from the user-defined type RHS
   b = RHS%rmat

   which_precon = 1
   if ( lu ) then
! Shift-and-invert
      which_precon = 2
   end if
   P = real_precon( which_precon, K )

! Solve the system:
   x = pp_idrs( K, b, P=P )

! Make plot?
   include "ocean_plot.inc"
 
   if ( my_proc == 1 ) &
      write(*,'(a)') '========================================================================='

end program stommel

program wedge_neumann
!
! Sound propagation test problem at single frequencies.
! The test problem is a Helmholtz equation models sound propagation
! in the earth crust. The sound velocity has a wedge shape.
! The boundaries are reflection (Neumann condition). The
! matrix is real and symmetric. The test problem is supplied on six
! different grids of increasingly fine resolution. Selection of the grid
! is done by specifying as command line option
!    -grid k, with k a number form 1 to 6.
! The problem is solved for frequency equal to 2**(k-1)
!
! This test problem is among others used in:
!    M.B. van Gijzen, Y.A. Erlangga, and C. Vuik. 
!    Spectral analysis of the discrete Helmholtz operator preconditioned with a shifted Laplacian. 
!    SIAM Journal on Scientific Computing, 29:1942-1958, 2007.
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2024 Martin van Gijzen
!

   use interface_module
   use precision_module
   use matrix_module
   use user_module
   use mm_module
   use idrs_module
   use ritz_module

   implicit none 
!
   type(user_matrix)              :: K, RHS
   real(kind=rp), allocatable     :: x(:), b(:), x0(:), U0(:,:), omega(:), H(:,:)
   real(kind=rp)                  :: pi

   integer                        :: neq

   integer                        :: tb, te, clock_rate, clock_max

   real(kind=rp), allocatable     :: D(:)

   character(len=60)              :: matrix_file, rhs_file
   integer                        :: grid 
   integer, parameter             :: default_grid = 1

! Declaration and initialization of IDRS parameters:
   include "../INCLUDE/initialize.inc"

   pi = 4.*atan(1._rp)

! Which grid?
   grid = get_parameter( '-grid', default_grid )

! Output message, problem description:
   if ( this_image() == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      write(*,'(a,i2,/)') 'Single frequency wedge problem with Neumann condition, grid is ', grid
   end if

   write(matrix_file,'(a,i1,a)') '../../DATA/WEDGE_DATA/SINGLE_NEUMANN/wedge', grid,'.mtx'
   write(rhs_file,'(a,i1,a)')    '../../DATA/WEDGE_DATA/SINGLE_NEUMANN/wedge', grid,'_b.mtx'

! Read the system matrix:
   K = mm_matrix( matrix_file )

! Problem size
   neq = K%nrows

! Allocate solution and rhs vectors
   allocate( x(neq), b(neq), x0(neq) )
   x0 = 0.

! Read the rhs:
   RHS = mm_matrix( rhs_file )
   b = RHS%rmat(:,1)

! Store in the matrix structure:
   A = real_matrix( K )

! Jacobi preconditioner
   allocate( D(neq) )
   D = real_diagonal( K )

! Initial solve to determine iteration parameters, set-up preconditioner:
   include "../INCLUDE/real_preconditioner.inc"

! Output message, which methods:
   include "../INCLUDE/output_methods.inc"

! Solve the system:
   include "../INCLUDE/run_idrs.inc"
!
   if ( this_image() == 1 ) &
      write(*,'(a)') '========================================================================='
!

end program wedge_neumann

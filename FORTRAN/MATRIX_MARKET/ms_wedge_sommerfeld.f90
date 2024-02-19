program ms_wedge_sommerfeld
!
! Sound propagation test problem at multiple frequencies.
! The test problem is a Helmholtz equation models sound propagation
! in the earth crust. The sound velocity has a wedge shape.
! The boundaries are absorbing (Sommerfeld condition). The
! matrix is complex and nonsymmetric. The test problem is supplied on six
! different grids of increasingly fine resolution. Selection of the grid
! is done by specifying as command line option
!    -grid k, with k a number form 1 to 6.
! The problem is solved for k different frequencies equal to 2**(k-1)
!
! This test problem is among others used in:
!    M. Baumann and M.B. van Gijzen. 
!    Nested Krylov Methods for Shifted Linear Systems. 
!    SIAM Journal on Scientific Computing, 37(5):S90-S112, 2015.
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
   type(user_matrix)              :: K, M1, RHS
   complex(kind=cp), allocatable  :: x(:,:), b(:), x0(:), omega(:), H(:,:)
   real(kind=rp)                  :: pi
   integer                        :: i, j

   integer                        :: nn, neq

   integer                        :: tb, te, clock_rate, clock_max

   complex(kind=cp), allocatable  :: sigma(:), sigma_p(:)
   integer                        :: Nsigma

   real(kind=rp), allocatable     :: freq(:)
   complex(kind=cp), allocatable  :: D(:)
   real(kind=rp)                  :: Dmin

   character(len=60)              :: stiffness_file, rhs_file, mass_file
   integer                        :: grid 
   integer, parameter             :: default_grid = 1

! Declaration and initialization of IDRS parameters:
   include "../INCLUDE/initialize.inc"

   pi = 4._rp*atan(1._rp)

! Which grid?
   grid = get_parameter( '-grid', default_grid )

! Output message, problem description:
   if ( this_image() == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      write(*,'(a,i2,/)') 'Multiple frequency wedge problem with Sommerfeld condition, grid is ', grid
   end if

! Frequencies:
   allocate(freq(grid))
   freq(1) = 1._rp
   do i = 2, grid
      freq(i) = freq(i-1)*2.
   end do

! Set the shifts for the wedge problem:
   Nsigma = grid
   allocate( sigma(Nsigma), sigma_p(Nsigma), resnrm(Nsigma) )
   sigma = cmplx(0.,2.*pi*freq,kind=cp)

   write(stiffness_file,'(a,i1,a)') '../../DATA/WEDGE_DATA/MULTI_SOMMERFELD/wedge', grid,'_K.mtx'
   write(mass_file,'(a,i1,a)')      '../../DATA/WEDGE_DATA/MULTI_SOMMERFELD/wedge', grid,'_M.mtx'
   write(rhs_file,'(a,i1,a)')       '../../DATA/WEDGE_DATA/MULTI_SOMMERFELD/wedge', grid,'_b.mtx'

! Read the system matrix:
   K = mm_matrix( stiffness_file )

! Problem size
   neq = K%nrows

! Allocate solution and rhs vectors
   allocate( x(neq,grid), x0(neq), b(neq) )
   x0 = 0.

! Read the mass matrix:
   M1 = mm_matrix( mass_file )

! Read the rhs:
   RHS = mm_matrix( rhs_file )
   b = RHS%rmat(:,1)
 
! Store the system matrix in the matrix structure:
   A  = complex_matrix( K )

! Assemble mass matrix
   allocate( D(neq) )
   D = complex_diagonal( M1 )

! Make the polynomial preconditioner:
   Dmin = minval(real(D))
   call co_min(Dmin)
   ccenter = 0.
   cfoci(1) = -4._cp/Dmin
   cfoci(2) = +4._cp/Dmin
   cseed    = (0.,4._cp)/Dmin
! Make the polynomial preconditioner:
   include "../INCLUDE/complex_ms_preconditioner.inc"

! Output message, which methods:
   include "../INCLUDE/output_methods.inc"

! Solve the system:
   include "../INCLUDE/run_ms_idrs.inc"

!
   if ( this_image() == 1 ) &
      write(*,'(a)') '========================================================================='
!

end program ms_wedge_sommerfeld
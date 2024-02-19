program room

!
! This programme illustrates the use of the idrs-module on
! a convection-diffusion-reaction problem with constant coefficients
! and homogeneous dirichlet boundary conditions.
! The problem can be  solved with the different algorithms: IDRS, BiCGSTAB and QMRIDR,
! and Neumann and Chebyshev polynomial preconditioners.
!
! The matrix is not explicitly stored. The matrix-vector multiplication
! is done using a stencil operations.
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2023 Martin van Gijzen
!

   use interface_module          ! To read command line options
   use precision_module          ! Set real and complex precision
   use matrix_module             ! Always needed is a matrix module that defines the 
                                 ! matrix-vector multiplication and preconditioning operation
   use user_module               ! User functions, user defined matrix-vector product
   use idrs_module               ! The module with the IDRS routines
   use ritz_module               ! Can be used for computing an initial search space 
                                 ! and user defined omegas
   implicit none 

   integer                       :: grid
   integer, parameter            :: default_grid = 3
!
   type(user_matrix)             :: A1

   complex(kind=cp), allocatable :: x(:), x0(:), b(:), U0(:,:), omega(:), H(:,:), D(:)

   integer                       :: tb, te, clock_rate, clock_max  ! For system_clock
   integer                       :: my_proc, n_procs ! For parallel processing

   real(kind=rp)                 :: radius

   real(kind=rp), parameter      :: sound_velocity = 340.
   real(kind=rp)                 :: wavenumber, freq
   real(kind=rp), parameter      :: pi = 3.141592653589793
   logical, parameter            :: default_impedance = .false.
   logical                       :: impedance

! Declaration and initialization of IDRS parameters:
   include "../INCLUDE/initialize.inc"

! Initialization for parallel processing
   my_proc = this_image()
   n_procs = num_images()

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Definition of the problem. This is the part that has to be made by the user
! The routines used here should be put in the user_module, together with
! the definition of the user type user_matrix.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Which grid?
   grid = get_parameter( '-grid', default_grid )

! On grid = 1 the frequency = 100Hz
   freq = 2**(grid-1)*1.e2_rp
   wavenumber = 2._cp*pi*freq/sound_velocity

! Sommerfeld problem?
   impedance = get_parameter( '-impedance', default_impedance )

! Make the matrix
   call helmholtz_matrix( A1, grid, wavenumber, .false., impedance )

! Allocate solution and right-hand-side vectors:
   allocate( x(A1%Nn), x0(A1%Nn), b(A1%Nn) )
   x0 = 0.

! Right-hand-side
   call point_source(  b, A1, A1%Lx/2., A1%Ly/2., A1%Lz/2. )

! Compute estimates for largest and smallest eigenvalues with Gershgorin's theorem:
   call gershgorin( A1, ccenter, radius )
   cfoci(1) = ccenter-radius   
   cfoci(2) = ccenter+radius   
   cseed = (1.,1.)*wavenumber**2

   if ( my_proc == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      if ( impedance ) then
         write(*,'(a,i9,/)') 'Room problem with impedance boundary condition, number of equations is ', A1%Nn
      else
         write(*,'(a,i9,/)') 'Room problem with Sommerfeld boundary conditions, number of equations is ', A1%Nn
      end if
   end if

! Store in the matrix type needed by IDRS
   A = complex_matrix( A1 ) ! system matrix

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! User part is ready. The rest is only standard IDRS-test
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Initial solve to determine iteration parameters, set-up preconditioner:
   include "../INCLUDE/complex_preconditioner.inc"

! Output message, which methods:
   include "../INCLUDE/output_methods.inc"

! Solve the system:
   include "../INCLUDE/run_idrs.inc"
!
   if ( this_image() == 1 ) &
      write(*,'(a)') '========================================================================='

!

end program room


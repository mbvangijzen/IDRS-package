program ms_room

!
! This programme illustrates the use of the idrs-module on
! a Helmholtz problem with constant wave numbers, mutliple frequencies, 
! and homogeneous Sommerfeld and Neumann conditions
! The problem can be  solved with the different algorithms: MSIDRS, MSBiCGSTAB and MSQMRIDR,
! and Neumann and Chebyshev polynomial preconditioners.
!
! The matrix is not explicitly stored. The matrix-vector multiplication
! is done using a stencil operations.
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2024 Martin van Gijzen
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
   logical                       :: impedance
   logical, parameter            :: default_impedance = .false.
!
   integer                       :: Nsigma
   complex(kind=cp), allocatable :: sigma(:)           ! Shifts
   complex(kind=cp), allocatable :: sigma_p(:)         ! Shifts preconditioned problem

   real(kind=rp)                 :: freq                 
   real, parameter               :: default_frequency = 100.  
   real(kind=rp), parameter      :: sound_velocity = 340.  

   type(user_matrix)             :: A1

   complex(kind=cp), allocatable  :: x(:,:), b(:), x0(:), omega(:), H(:,:), D(:)

   integer                       :: tb, te, clock_rate, clock_max  ! For system_clock
   integer                       :: my_proc, n_procs

   real(kind=rp)                 :: radius
   real(kind=rp), parameter      :: pi = 3.141592653589793
   integer                       :: i

! Declaration and initialization of IDRS parameters:
   include "../INCLUDE/initialize.inc"

   my_proc = this_image()
   n_procs = num_images()

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Definition of the problem. This is the part that has to be made by the user
! The routines used here should be put in the user_module, together with
! the definition of the user type user_matrix.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Which grid?
   grid = get_parameter( '-grid', default_grid )

! Which grid?
   impedance = get_parameter( '-impedance', default_impedance )

! Make the Laplace matrix plus Sommerfeld boundary conditions
   call helmholtz_matrix( A1, grid, multishift = .true., impedance = impedance )

! Determine shifts:
   Nsigma = grid
   allocate( sigma(Nsigma), sigma_p(Nsigma), resnrm(Nsigma) )
   freq = get_parameter( '-frequency', default_frequency )
   do i = 1, Nsigma
      sigma(i) = (0.,1.)*(2*pi*freq/sound_velocity)
      freq = 2.*freq
   end do

! Allocate solution and right-hand-side vectors:
   allocate( x(2*A1%Nn,Nsigma), x0(2*A1%Nn), b(2*A1%Nn) )
   x0 = 0.

! Right-hand-side
   call point_source(b, A1, A1%Lx/2., A1%Ly/2., A1%Lz/2. )
   b = -b

! Compute estimates for largest and smallest eigenvalues with Gershgorin's theorem:
   call gershgorin( A1, ccenter, radius )
   cfoci(1) = ccenter-radius*(0.,1.)   
   cfoci(2) = ccenter+radius*(0.,1.)   
   cseed    = (1.,1.)*sigma(Nsigma)

   if ( my_proc == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      if ( impedance ) then
         write(*,'(a,i9,/)') 'Multishift room problem with impedance boundary condition, number of equations is ', A1%Nn
      else
         write(*,'(a,i9,/)') 'Multishift room problem with Sommerfeld boundary conditions, number of equations is ', A1%Nn
      end if
   end if

! Store in the matrix type needed by IDRS
! Note, here we will include the preconditioner in A, no M-matrix needed.
   A = complex_matrix( A1 ) ! system matrix

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! User part is ready. The rest is only standard IDRS-test
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

end program ms_room


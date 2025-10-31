program room

!
! This programme illustrates the use of the idrs-module on
! a Helmholtz problem with constant wave numbers, mutliple frequencies, 
! and homogeneous Sommerfeld and Neumann conditions.
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
   use user_module               ! User functions, user defined matrix-vector product
   use pp_idrs_module            ! Preconditioned IDRS solver calls

   implicit none 

   integer                       :: grid
   integer, parameter            :: default_grid = 1
!
   integer                       :: isigma, Nsigma
   complex(kind=cp), allocatable :: sigma(:)           ! Shifts

   real(kind=rp)                 :: freq                 
   real(kind=rp), parameter      :: sound_velocity = 340._rp

   integer                       :: neq

   type(user_matrix)             :: K, C

   complex(kind=cp), allocatable :: x(:,:,:), b(:,:)

   logical                        :: multishift
   logical, parameter             :: default_multishift = .false.

   logical                        :: deflation
   logical, parameter             :: default_deflation = .false.
   complex(kind=cp), allocatable  :: V(:,:)

   integer                       :: my_proc, n_procs

   real(kind=rp), parameter      :: pi = 3.141592653589793_cp
   complex(kind=cp)              :: cseed
   integer                       :: i


! Read command line
   call initialize()

   my_proc = this_image()
   n_procs = num_images()

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Definition of the problem. This is the part that has to be made by the user
! The routines used here should be put in the user_module, together with
! the definition of the user type user_matrix.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Which grid?
   grid = get_parameter( '-grid', default_grid )

! Multi frequency problem?
   multishift = get_parameter( '-multi', default_multishift )

! Deflation?
   deflation = get_parameter( '-deflation', default_deflation )

! Make the matrices: 
   call laplace_matrix( K, grid )
   call damping_matrix( K, C )

! Output message, problem description:
   if ( my_proc == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      write(*,'(a/)') 'Room problem with Sommerfeld boundary conditions.'
      if ( multishift ) then
         write(*,'(a)') 'Multiple frequency wedge problem.'
      else
         write(*,'(a)') 'Single frequency problem.'
      end if
   end if

! Determine shifts:
   Nsigma = grid
   allocate( sigma(Nsigma) )
   freq = 1.e2_rp
   do i = 1, Nsigma
      sigma(i) = (0.,1._rp)*(2._rp*pi*freq)/sound_velocity
      freq = 2._rp*freq
   end do

! Number of equations:
   neq = K%Nn

! Allocate solution and right-hand-side vectors:
   allocate( x(neq,Nsigma,n_procs), b(neq,n_procs) )

! Right-hand-side
   call point_source( b, K )

! Deflation:
   if ( deflation ) then
      if ( multishift ) then
         allocate(V(2*neq,1))
         V = 0._cp
         V(neq+1:2*neq,1) = 1._cp/sqrt(real(neq*n_procs,kind=cp))
      else
         allocate(V(neq,1))
         V(:,1) = 1._cp/sqrt(real(neq*n_procs,kind=cp))
      end if
   end if

   if ( multishift ) then
      cseed = sigma(Nsigma)
      if ( deflation ) then
         x = pp_idrs( K, b, sigma, C=C, seed=cseed, V = V )
      else
         x = pp_idrs( K, b, sigma, C=C, seed=cseed )
      end if
   else
      do isigma = 1, Nsigma
         cseed = sigma(isigma)
         if ( deflation ) then
            x(:,isigma,:) = pp_idrs( K, b, C=C, seed=cseed, V = V )
         else
            x(:,isigma,:) = pp_idrs( K, b, C=C, seed=cseed )
         end if
      end do
   end if

   if ( my_proc == 1 ) &
      write(*,'(a)') '========================================================================='

end program room


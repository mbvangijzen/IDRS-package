program wedge_sommerfeld
!
! Sound propagation test problem at multiple frequencies.
! The test problem is a Helmholtz equation models sound propagation 
! in the earth crust. The sound velocity has a wedge shape. 
! The boundaries are absorbing (Sommerfeld condition). The
! matrix is complex and general. The test problem is supplied on six
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
! Copyright:(c) 2025 Martin van Gijzen
!
!
   use interface_module
   use precision_module
   use user_module
   use mm_module
   use pp_idrs_module

   implicit none 
  
   type(user_matrix)              :: K, C, M, P, RHS
   complex(kind=cp), allocatable  :: x(:,:), b(:)
   real(kind=rp)                  :: pi

   integer                        :: nrows, ncols

   complex(kind=cp), allocatable  :: sigma(:)
   integer                        :: Nsigma, isigma

   real(kind=rp), allocatable     :: freq(:)

   complex(kind=cp)               :: cseed, cshift

   character(len=60)              :: stiffness_file, rhs_file, mass_file, damping_file
   integer                        :: grid 
   integer, parameter             :: default_grid = 1

   logical                        :: multishift
   logical, parameter             :: default_multishift = .false.

! For preconditioner:
   integer                        :: which_precon

   logical                        :: lu
   logical, parameter             :: default_lu = .false.

   logical                        :: deflation
   logical, parameter             :: default_deflation = .false.
   complex(kind=cp), allocatable  :: V(:,:)

! Plotting:
   character(len=:), allocatable  :: plot
   character(len=4), parameter    :: default_plot = 'none'
   logical                        :: plot_solution
   logical                        :: gif, jpeg
   logical, parameter             :: default_gif  = .false.
   logical, parameter             :: default_jpeg = .false.

   integer                        :: igrid, Ngrid, nb, ne, nx, ny, i, j
   integer                        :: my_proc, num_procs, bandwidth

! Initialization for parameter routines (read command line)
   call initialize( )
   my_proc   = this_image()
   num_procs = num_images()

   pi = 4._rp*atan(1._rp)

! Which grid?
   grid = get_parameter( '-grid', default_grid )

! Multishift problem?
   multishift = get_parameter( '-multi', default_multishift )

! Direct solution?
   lu = get_parameter( '-lu', default_lu )

! Deflation?
   deflation = get_parameter( '-deflation', default_deflation )

! Output message, problem description:
   if ( my_proc == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      if ( multishift ) then
         if ( lu ) then
            write(*,'(a)') 'Multiple frequency wedge problem with band Shift-and-Invert preconditioner.'
         else
            write(*,'(a)') 'Multiple frequency wedge problem with lumped mass matrix.'
         end if
      else
         if ( lu ) then
            write(*,'(a)') 'Single frequency wedge problem with band Shift-and-Invert preconditioner.'
         else
            write(*,'(a)') 'Single frequency wedge problem.'
         end if
      end if
      write(*,'(a,i2,/)') 'Sommerfeld boundary conditions. Grid is ',grid
   end if

! Frequencies:
   allocate(freq(grid))
   freq(1) = 1._rp
   do isigma = 2, grid
      freq(isigma) = freq(isigma-1)*2._rp
   end do
  
   write(stiffness_file,'(a,i1,a)') 'DATA/WEDGE_DATA/wedge', grid,'_K.mtx'
   write(mass_file,'(a,i1,a)')      'DATA/WEDGE_DATA/wedge', grid,'_M.mtx'
   write(damping_file,'(a,i1,a)')   'DATA/WEDGE_DATA/wedge', grid,'_C.mtx'
   write(rhs_file,'(a,i1,a)')       'DATA/WEDGE_DATA/wedge', grid,'_b.mtx'

! Read the stiffness matrix:
   K = mm_matrix( stiffness_file, 'crs' )

! Read the damping matrix:
   C = mm_matrix( damping_file )

! Read the mass matrix:
   M = mm_matrix( mass_file )

! Read the rhs:
   RHS = mm_matrix( rhs_file )

! Problem size
   nrows = K%nrows
   ncols = K%ncols
   Nsigma = grid

! Allocate solution and rhs vectors, and the shifts
   allocate( x(nrows,Nsigma), b(nrows),sigma(Nsigma) )

! Set the shifts and right-hand-side vector:
   b = RHS%rmat(:,1)
   sigma = (2._cp*pi*freq)*(0._cp,1._cp)

! Preconditioner:
   cshift = 2._cp*pi/(freq(Nsigma)+freq(1))*cmplx((freq(Nsigma)-freq(1))*sqrt(freq(1)*freq(Nsigma)),2._cp*freq(1)*freq(Nsigma))
   if ( lu ) then
! Shift-and-invert
      which_precon = 3
      P = complex_precon( which_precon, K, M, C, shift=cshift )
   elseif ( multishift ) then
! Lumped mass matrix
      which_precon = 1
      P = complex_precon( which_precon, M )
   end if

! Deflation:
   if ( deflation ) then
      if ( multishift ) then
         allocate(V(2*nrows,1))
         V = 0._cp
         V(nrows+1:2*nrows,1) = 1._cp/sqrt(real(ncols,kind=cp))
      else
         allocate(V(nrows,1))
         V(:,1) = 1._cp/sqrt(real(ncols,kind=cp))
      end if
   end if
 
   if ( .not. multishift ) then
      do isigma = 1, Nsigma
         cseed = sigma(isigma)
         if ( deflation ) then
            x(:,isigma) = pp_idrs( K, b, M, C, P, seed=cseed, V=V )
         else
            x(:,isigma) = pp_idrs( K, b, M, C, P, seed=cseed )
         end if
      end do
   else
      cseed = sigma(Nsigma)
      if ( lu .and. deflation) then
         x = pp_idrs( K, b, sigma, M, C, P, seed=cseed, shift=cshift, V=V )
      elseif ( lu ) then
         x = pp_idrs( K, b, sigma, M, C, P, seed=cseed, shift=cshift )
      elseif ( deflation ) then
         x = pp_idrs( K, b, sigma, M, C, P, seed=cseed, V=V )
      else
         x = pp_idrs( K, b, sigma, M, C, P, seed=cseed )
      end if
   end if

! Plot solution?
   include "wedge_plot.inc"

   if ( my_proc == 1 ) &
      write(*,'(a)') '========================================================================='
!
end program wedge_sommerfeld

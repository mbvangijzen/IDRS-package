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
! Copyright:(c) 2026 Martin van Gijzen
!
!
   use interface_module
   use precision_module
   use user_module
   use mm_module
   use idrs_driver_module

   implicit none 
  
   type(user_matrix)              :: K, C, M, RHS
   complex(kind=cp), allocatable  :: x(:,:), b(:)
   complex(kind=cp), allocatable  :: lambda(:)
   real(kind=rp), allocatable     :: eigfreq(:)
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

   logical                        :: eigen
   integer, parameter             :: default_neig = 0             ! Number of eigenvalues
   integer                        :: neig, ieig

! Plotting:
   character(len=:), allocatable  :: plot
   character(len=4), parameter    :: default_plot = 'none'
   logical                        :: plot_solution
   logical                        :: gif, jpeg
   logical, parameter             :: default_gif  = .false.
   logical, parameter             :: default_jpeg = .false.

! Points in each direction on coarsest grid:
   integer, parameter             :: px =  6
   integer, parameter             :: py = 10

   integer                        :: igrid, Ngrid, nb, ne, nx, ny, i, j
   integer                        :: my_proc, num_procs, n_solns

! Initialization for parameter routines (read command line)
   call initialize( )
   my_proc   = this_image()
   num_procs = num_images()

   pi = 4._rp*atan(1._rp)

! Which grid?
   grid = get_parameter( '-grid', default_grid )

! Eigenvalue problem?
   neig = get_parameter( '-eigen', default_neig )
   eigen = ( neig > 0 )

! Multishift problem?
   multishift = get_parameter( '-multi', default_multishift )

! Output message, problem description:
   if ( my_proc == 1 ) then
      write(*,*)
      write(*,'(a)') '========================================================================='
      if ( neig > 0 ) then
         write(*,'(a)') 'Eigenvalue wedge problem.'
      elseif ( multishift ) then
         write(*,'(a)') 'Multiple frequency wedge problem.'
      else
         write(*,'(a)') 'Single frequency wedge problem.'
      end if
      write(*,'(a,i2,/)') 'Sommerfeld boundary conditions. Grid is ',grid
   end if
  
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

! Number of shifts:
   Nsigma = grid

! Frequencies:
   allocate(freq(Nsigma),sigma(Nsigma))
   freq(1) = 1._rp
   do isigma = 2, grid
      freq(isigma) = freq(isigma-1)*2._rp
   end do

! Set the shifts:
   sigma = (0._cp,1._cp)*(2._cp*pi*freq)

! Eigenvalues and eigenfrequencies:
   allocate( lambda(neig), eigfreq(neig) )

! Number of solutions to be computed (neig for eigen, Nsigma for multishift)
   if ( eigen ) then
      n_solns = neig
   else
      n_solns = Nsigma
   end if

! Allocate solution and rhs vectors
   allocate( x(nrows,n_solns), b(nrows) )
   b = RHS%rmat(:,1)

   if ( eigen ) then
!
! Eigenvalue analysis
!

      cshift = sigma(Nsigma)

      call idreig_driver( lambda, x, K, b, M, C, shift = cshift  )
      eigfreq = aimag(lambda)/(2._cp*pi)

   elseif ( multishift ) then
!
! Multiple frequency problem
!

!     cshift = 2._cp*pi/(freq(Nsigma)+freq(1))*cmplx((freq(Nsigma)-freq(1))*sqrt(freq(1)*freq(Nsigma)),2._cp*freq(1)*freq(Nsigma))
      x = idrs_driver( K, b, sigma, M, C )

   else
!
! Sequence of shifted problems
!
      do isigma = 1, Nsigma
         cseed = sigma(isigma)
         x(:,isigma) = idrs_driver( K, b, M, C, seed=cseed )
      end do

   end if

! Plot solution?
   include "solns_plot.inc"

   if ( my_proc == 1 ) &
      write(*,'(a)') '========================================================================='
!
end program wedge_sommerfeld

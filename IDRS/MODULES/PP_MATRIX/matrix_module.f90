!
! Matrix module
!   This module provides routines for the matrix-vector multiplication
!   and preconditioning operation for use in the IDRS module. The preconditioners
!   implemented here are polynomial preconditioners and deflation preconditioners. 
!   These preconditioners can be combined with a user defined preconditioner.
!
!   The polynomial preconditioners can be used both for solving standard 
!   linear systems Ax = b and for sequences of shifted linear systems (A - sI) x = b. 
!
!   A user defined preconditioner can be used in combination with the polynomial 
!   preconditioner. In the standard case the scaled system
!        AP^-1 y = b, x = P^-1y is solved.
!   For shifted problems the problem
!        (AP^-1 - sI) y = b, x = P^-1y
!   is solved. Here P plays the role of a mass matrix, which
!   allows for solving problems of the form
!        (A - sP) x = b
!  
!   The only operations that have to be supplied by a user are the matrix-vector multiplication and 
!   preconditioning operation, in real and complex versions.
!   These operations are defined in the module user_module. Multiplication with matrices in the standard 
!   formats DENSE, CRS, COO, and CDS are supplied in the accompanying file user_module.f90.
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2025 Martin van Gijzen
!

module matrix_module

   use precision_module
   use dense_la_module
   use user_module

   implicit none

   type :: matrix

! User defined matrix:
      type(user_matrix)              :: K, M, C, P

! Seed shift:
      real(kind=rp)                  :: rseed = 0._rp
      complex(kind=cp)               :: cseed = (0._cp,0._cp)

! shift of SI preconditioner:
      real(kind=rp)                  :: rshift = 0._rp
      complex(kind=cp)               :: cshift = (0._cp,0._cp)

! Wich polynomial preconditioner?       
      integer                        :: preconditioner = 0

! Polynomial preconditioner:
      integer                        :: degree = 0
      real(kind=rp), allocatable     :: ralpha_p(:)
      real(kind=rp), allocatable     :: ralpha_s(:,:)
      complex(kind=cp), allocatable  :: calpha_p(:)
      complex(kind=cp), allocatable  :: calpha_s(:,:)

! For deflation preconditioner:
      integer                        :: n_deflation = 0
      real(kind=rp), allocatable     :: rV(:,:), rAV(:,:), rVtAV(:,:)
      complex(kind=cp), allocatable  :: cV(:,:), cAV(:,:), cVtAV(:,:)

! Do we have problem with mass matrix?
      logical                        :: mass_matrix = .false.

! Do we have problem with damping matrix?
      logical                        :: damping_matrix = .false.

! Do we have a multishift problem?
      logical                        :: multishift = .false.

! Do we have a block system?      
      logical                        :: block_system = .false.

! Do we have a shift-and-invert preconditioner?
      logical                        :: shift_invert = .false.

   end type matrix

! Overload * to define the matrix-vector multiplication using the matrix type

   INTERFACE OPERATOR(*)
      module procedure rmatvec, cmatvec
   END INTERFACE

   INTERFACE OPERATOR(/)
      module procedure rprecon, cprecon
   END INTERFACE

   interface CHEBYSHEV
      module procedure RCHEBYSHEV, CCHEBYSHEV
   end interface

   interface NEUMANN
      module procedure RNEUMANN, CNEUMANN
   end interface

   interface SCALEBACK
      module procedure CSCALEBACK, RSCALEBACK, MS_CSCALEBACK, MS_RSCALEBACK
   end interface

   interface COMPUTE_ALPHA
      module procedure CCOMPUTE_ALPHA, RCOMPUTE_ALPHA
   end interface

   interface POLPRE
      module procedure CPOLPRE, RPOLPRE
   end interface

   interface MS_PARAMETERS
      module procedure COMPLEX_MS_PARAMETERS, REAL_MS_PARAMETERS
   end interface

   interface ANALYSE_SPECTRUM
      module procedure CANALYSE_SPECTRUM, RANALYSE_SPECTRUM
   end interface

   interface Pmul
      module procedure rPmul, cPmul
   end interface

   interface make_deflation
      module procedure rmake_deflation, cmake_deflation
   end interface

   interface add_deflation
      module procedure radd_deflation, cadd_deflation, rms_add_deflation, cms_add_deflation
   end interface

contains

   function cmatvec( A, v ) result(w) 

      type(matrix), intent(in)          :: A
      complex(kind=cp), intent(in)      :: v(:)
      complex(kind=cp)                  :: w(size(v))

! Polynomial preconditioner:
      w = polpre( A, v )

! Multiplication with system matrix:
      w = complex_mv( A, w )

! Multiplication with deflation preconditioner:
      w = Pmul( A, w )

   end function cmatvec

   function rmatvec( A, v ) result(w)

      type(matrix), intent(in)          :: A
      real(kind=rp), intent(in)         :: v(:)
      real(kind=rp)                     :: w(size(v))

! Polynomial preconditioner:
      w = polpre( A, v )

! Multiplication with system matrix:
      w = real_mv( A, w )

! Multiplication with deflation preconditioner:
      w = Pmul( A, w )

   end function rmatvec

   function complex_mv( A, v ) result(w) 
!
! Multiplication with the (possibly shifted and scaled) system matrix
!

      type(matrix), intent(in)          :: A
      complex(kind=cp), intent(in)      :: v(:)
      complex(kind=cp)                  :: w(size(v))
      complex(kind=cp)                  :: t(size(v))
      integer                           :: neq, nn
!
      neq = size(v)
!
! Preconditioning:
      t = v/A
!
! Multiplication with the system matrix
!
! Block matrix?
      if ( A%block_system ) then
         nn = neq/2
         w(1:nn)     = -A%C*t(1:nn) -A%K*t(nn+1:neq)
         w(nn+1:neq) = t(1:nn)
         if ( A%cseed /= (0._cp,0._cp) ) then
            if ( A%mass_matrix ) then
               w(1:nn)      = w(1:nn)      - A%cseed*(A%M*t(1:nn))
               w(nn+1:2*nn) = w(nn+1:2*nn) - A%cseed*t(nn+1:2*nn)
            else
               w = w - A%cseed*t
            end if
         end if
      else
!
! Standard matrix
         if ( A%cseed == (0._cp,0._cp) ) then
            w = A%K*t
         elseif ( A%damping_matrix .and. A%mass_matrix ) then
            w = A%K*t + A%cseed*(A%C*t) + (A%cseed**2)*(A%M*t)
         elseif ( A%damping_matrix ) then
            w = A%K*t + A%cseed*(A%C*t) + (A%cseed**2)*t
         elseif ( A%mass_matrix ) then
            w = A%K*t - A%cseed*(A%M*t)
         else
            w = A%K*t - A%cseed*t
         end if
      end if

   end function complex_mv

   function real_mv( A, v ) result(w)
!
! Multiplication with the (possibly shifted and scaled) system matrix
!
      type(matrix), intent(in)          :: A
      real(kind=rp), intent(in)         :: v(:)
      real(kind=rp)                     :: w(size(v))
      real(kind=rp)                     :: t(size(v))
      integer                           :: neq, nn
!
      neq = size(v)
!
! Preconditioning:
      t = v/A
!
! Multiplication with the system matrix
!
! Block matrix?
      if ( A%block_system ) then
         nn = neq/2
         w(1:nn)     = -A%C*t(1:nn) -A%K*t(nn+1:neq)
         w(nn+1:neq) = t(1:nn)
         if ( A%rseed /= (0._rp,0._rp) ) then
            if ( A%mass_matrix ) then
               w(1:nn)      = w(1:nn)      - A%rseed*(A%M*t(1:nn))
               w(nn+1:2*nn) = w(nn+1:2*nn) - A%rseed*t(nn+1:2*nn)
            else
               w = w - A%rseed*t
            end if
         end if
      else
!
! Standard matrix
         if ( A%rseed == 0._rp ) then
            w = A%K*t
         elseif ( A%damping_matrix .and. A%mass_matrix ) then
            w = A%K*t + A%rseed*(A%C*t) + (A%rseed**2)*(A%M*t)
         elseif ( A%damping_matrix ) then
            w = A%K*t + A%rseed*(A%C*t) + (A%rseed**2)*t
         elseif ( A%mass_matrix ) then 
            w = A%K*t - A%rseed*(A%M*t)
         else
            w = A%K*t - A%rseed*t
         end if
      end if

   end function real_mv

   function rpolpre( A, v ) result(w)
!
! Evaluate polynomial preconditioner
!
      type(matrix), intent(in)              :: A
      real(kind=rp), intent(in)             :: v(:)
      real(kind=rp)                         :: w(size(v))
      integer                               :: i

      if ( A%degree == 0 ) then
         w = v
      else
         w = A%ralpha_p(A%degree)*v
         do i = A%degree-1, 0, -1
            w = real_mv(A, w)  + A%ralpha_p(i)*v
         end do
      end if

   end function rpolpre

   function cpolpre( A, v ) result(w)
!
! Evaluate polynomial preconditioner
!
      type(matrix), intent(in)              :: A
      complex(kind=cp), intent(in)          :: v(:)
      complex(kind=cp)                      :: w(size(v))
      integer                               :: i

      if ( A%degree == 0 ) then
         w = v
      else
         w = A%calpha_p(A%degree)*v
         do i = A%degree-1, 0, -1
            w = complex_mv(A, w)  + A%calpha_p(i)*v
         end do
      end if

   end function cpolpre

   function rprecon( v, A ) result(w)
!
! Solve system with the preconditioner
!
      type(matrix), intent(in)          :: A
      real(kind=rp), intent(in)         :: v(:)
      real(kind=rp)                     :: w(size(v))
      real(kind=rp)                     :: t(size(v))
      integer                           :: nn, neq

      neq = size(v)

      if ( A%block_system ) then
         nn = neq/2
         if ( A%shift_invert )then
!
! We have a block matrix of the form
!    A = [-C-sM   -K]
!        [  I    -sI]
!
! This matrix can be decomposed as
!  [-C-sM   -K] = [-C-sM I][I         O][I -sI]
!  [ I     -sI]   [I     O][O K+sC+s^2M][O  -I]
!
! Inversion gives:
! A^-1 = [I -sI][I              O][O    I]
!        [O  -I][0 (K+sC+s^2M)^-1][I C+sM]
!
            t(1:nn)     = v(nn+1:neq)
            if ( A%mass_matrix ) then
               t(nn+1:neq) = v(1:nn) + A%C*v(nn+1:neq) + A%rshift*(A%M*v(nn+1:neq))
            else
               t(nn+1:neq) = v(1:nn) + A%C*v(nn+1:neq) + A%rshift*v(nn+1:neq)
            end if
            t(nn+1:neq) = t(nn+1:neq)/A%P
            w(1:nn)     = t(1:nn) - A%rshift*t(nn+1:neq)
            w(nn+1:neq) = -t(nn+1:neq)
         else
            w(1:nn) = v(1:nn)/A%P
            w(nn+1:neq) = v(nn+1:neq)
         end if
      else
! Standard system
         w = v/A%P
      end if

   end function rprecon

   function cprecon( v, A ) result(w)
!
! Solve system with the preconditioner
!
      type(matrix), intent(in)          :: A
      complex(kind=cp), intent(in)      :: v(:)
      complex(kind=cp)                  :: w(size(v))
      complex(kind=cp)                  :: t(size(v))
      integer                           :: nn, neq

      neq = size(v)

      if ( A%block_system ) then
         nn = neq/2
         if ( A%shift_invert ) then
!
! We have a block matrix of the form
!    A = [-C-sM   -K]
!        [  I    -sI]
!
! This matrix can be decomposed as
!  [-C-sM   -K] = [-C-sM I][I         O][I -sI]
!  [ I     -sI]   [I     O][O K+sC+s^2M][O  -I]
!
! Inversion gives:
! A^-1 = [I -sI][I              O][O    I]
!        [O  -I][O (K+sC+s^2M)^-1][I C+sM]
!
            t(1:nn)     = v(nn+1:neq)
            if ( A%mass_matrix ) then
               t(nn+1:neq) = v(1:nn) + A%C*v(nn+1:neq) + A%cshift*(A%M*v(nn+1:neq))
            else
               t(nn+1:neq) = v(1:nn) + A%C*v(nn+1:neq) + A%cshift*v(nn+1:neq)
            end if
            t(nn+1:neq) = t(nn+1:neq)/A%P
            w(1:nn)     = t(1:nn) - A%cshift*t(nn+1:neq)
            w(nn+1:neq) = -t(nn+1:neq)
         else
            w(1:nn) = v(1:nn)/A%P
            w(nn+1:neq) = v(nn+1:neq)
         end if
      else
! Standard system
         w = v/A%P
      end if

   end function cprecon

   subroutine complex_ms_parameters( A, sigma, sigma_p )
!
! Compute the shifts of the preconditioned system
! Compute the scaleback polynomial in case of a polynomial preconditioner
! Note: the polynomial preconditioner itself must then already have been constructed
!
! Parameters: 
!   A: matrix, the parameters of the multishift preconditioner are added to A
!   sigma: complex, input. The shifts
!   sigma_p: complex, output. The shifts of the preconditioned system
!
! Derivatiopn shift-and-invert preconditioner:
!
! (A-sigmaM)P^-1 =  (A-seedM)(A-shiftM)^-1 -eta I
! (A-sigmaM)P^-1 = ((A-seedM)-eta(A-shiftM))(A-shiftM)^-1
! (A-sigmaM)P^-1 = (1-eta)(A - (seed-eta shift)/(1-eta)M)(A-shiftM)^-1
! Choose (seed-eta shift)/(1-eta) = sigma
! or     (seed-eta shift) = sigma(1-eta)
!        (seed-sigma)= eta shift - eta sigma
!        eta = (seed - sigma)/(shift - sigma )
!        eta = (sigma-seed)/(sigma-shift )
!
!    (1-eta) = (sigma-shift - sigma+seed)/((sigma-shift )
!    (1-eta) = (seed-shift)/(sigma-shift )
!
!       P^-1 = (seed-shift)/(sigma-shift )(A-shiftM)^-1
!      

!
      type(matrix), intent(inout)                 :: A
      complex(kind=cp), intent(in)                :: sigma(:)
      complex(kind=cp), intent(out)               :: sigma_p(:)
      complex(kind=cp)                            :: sigma_s(size(sigma))
      integer                                     :: i, i_sigma, n_sigma

      if ( A%shift_invert ) then
!
! Check that the shift is not equal to the seed
         if ( A%cshift == A%cseed ) &
            error stop 'Preconditioner shift should not be equal to seed, choose preconditioner shift and seed differently'
! Check that the shift is not equal to one of the shifts
         do i_sigma = 1, n_sigma
            if ( A%cshift == sigma(i_sigma) ) &
               error stop 'Preconditioner shift should not be equal to one of the shifts, choose different preconditioner shift'
         end do

!
! Shift-and-invert:
!
         sigma_s = ( sigma-A%cseed )/( sigma-A%cshift ) 
      else
!
! Only shift:
!
         sigma_s = sigma-A%cseed
      end if

!
! Compute parameters of polynomial preconditioner:
!
      if ( A%degree > 0 ) then
!
         if ( .not. allocated(A%calpha_p) ) error stop 'First call subroutine complex_preconditioner '
! Allocate space for the parameters for the scaleback polynomials
         n_sigma = size(sigma)
         if ( allocated(A%calpha_s) ) deallocate(A%calpha_s)
         allocate(A%calpha_s(0:A%degree,n_sigma))

! Compute the parameters of the scaleback polynomials
         A%calpha_s = 0.
         do i_sigma = 1,n_sigma
            A%calpha_s(A%degree,i_sigma) = A%calpha_p(A%degree)
            do i = A%degree-1,0,-1
               A%calpha_s(i,i_sigma) = A%calpha_p(i) + sigma_s(i_sigma)*A%calpha_s(i+1,i_sigma)
            end do
            sigma_p(i_sigma) = sigma_s(i_sigma)*A%calpha_s(0,i_sigma)
         end do
      else
         sigma_p = sigma_s
      end if

   end subroutine complex_ms_parameters

   subroutine real_ms_parameters( A, sigma, sigma_p )
!
! Compute the shifts of the preconditioned system
! Compute the scaleback polynomial in case of a polynomial preconditioner
! Note: the polynomial preconditioner itself must then already have been constructed
!
! Parameters: 
!   A: matrix, the parameters of the multishift preconditioner are added to A
!   sigma: real, input. The shifts
!   sigma_p: real, output. The shifts of the preconditioned system
!
      type(matrix), intent(inout)                 :: A
      real(kind=rp), intent(in)                   :: sigma(:)
      real(kind=rp), intent(out)                  :: sigma_p(:)
      real(kind=rp)                               :: sigma_s(size(sigma))
      integer                                     :: i, i_sigma, n_sigma

      if ( A%shift_invert ) then
!
! Check that the shift is not equal to the seed
         if ( A%rshift == A%rseed ) &
            error stop 'Preconditioner shift should not be equal to seed, choose preconditioner shift and seed differently'
! Check that the shift is not equal to one of the shifts
         do i_sigma = 1, n_sigma
            if ( A%rshift == sigma(i_sigma) ) &
               error stop 'Preconditioner shift should not be equal to one of the shifts, choose different preconditioner shift'
         end do
!
! Shift-and-invert:
!
         sigma_s = ( sigma-A%rseed )/( sigma-A%rshift )
      else
!
! Only shift:
!
         sigma_s = sigma-A%rseed
      end if

!
! Compute parameters of polynomial preconditioner:
!
      if ( A%degree > 0 ) then

!
! Polynomial preconditioner:
!
         if ( .not. allocated(A%ralpha_p) ) error stop 'First call subroutine real_preconditioner '
! Allocate space for the parameters for the shifted polynomials
         n_sigma = size(sigma)
         allocate(A%ralpha_s(0:A%degree,n_sigma))

! Compute the parameters of the shifted systems:
         A%ralpha_s = 0.
         do i_sigma = 1,n_sigma
            A%ralpha_s(A%degree,i_sigma) = A%ralpha_p(A%degree)
            do i = A%degree-1,0,-1
               A%ralpha_s(i,i_sigma) = A%ralpha_p(i) + sigma_s(i_sigma)*A%ralpha_s(i+1,i_sigma)
            end do
            sigma_p(i_sigma) = sigma_s(i_sigma)*A%ralpha_s(0,i_sigma)
         end do
      else
         sigma_p = sigma_s
      end if

   end subroutine real_ms_parameters

   subroutine real_preconditioner( A, preconditioner, degree, foci, center, H )
!
! Create the preconditioner A. Note that this routine can also be used to make a preconditioned
! system matrix. If the preconditioner is included in the system matrix, the final IDRS solution
! should be scaled back using the subroutine SCALEBACK.
!
! All the parameters are optional except the first are optional.
!    real_preconditioner( A )
! returns an identity matrix.
!
! Parameters:
!    A: matrix, at input system matrix, at output the preconditioner is added
!    preconditioner: integer, input.
!       preconditioner = 0: no preconditioner
!       preconditioner = 1: neumann
!       preconditioner = 2: chebyshev
!       preconditioner = 3: hessenberg
!    degree: integer, input. Degree of the polynomial preconditioner
!    foci: real, input. Foci of ellipse around spectrum.
!    center: real, input. Center of spectrum.
!    seed: real, input. Shift to the matrix.
!    H: real, hessenberg matrix, input. This matrix is required for the
!             hessenberg preconditioner and in case foci or center are not specified.
!
      type(matrix), intent(inout)                 :: A
      integer, intent(in), optional               :: preconditioner, degree
      real(kind=rp), intent(in), optional         :: foci(2), center
      real(kind=rp)                               :: f(2), c, shift
      real(kind=rp), intent(in), optional         :: H(:,:)

! Polynomial preconditioning?
      if ( present(preconditioner) ) then
         A%preconditioner = preconditioner
      else 
         A%preconditioner = 0
      end if

! Degree of the polynomial preconditioner?
      if ( present(degree) ) then
         A%degree = degree
      else 
         A%degree = 0
      end if

! Construct preconditioner:
      select case(A%preconditioner)
         case(1)
! Neumann preconditioner:

            if ( present(center) ) then
               c = center - A%rseed
            elseif ( present(H) ) then
               call analyse_spectrum( H, center = c )
            else
               c = 1._rp
            end if

            allocate(A%ralpha_p(0:A%degree))
            A%ralpha_p = neumann( A%degree, c )

         case(2)
! Chebyshev preconditioner:

            if ( present(foci) ) then
               f = foci - A%rseed
            elseif ( present(H) ) then
               call analyse_spectrum( H, foci = f )
            else
! These defaults are working well for the Stommel test problems
               f(1) = 0.25
               f(2) = 1.75
            end if

            allocate(A%ralpha_p(0:A%degree))
            A%ralpha_p = chebyshev( A%degree, f )

         case(3)

! Hessenberg preconditioner
            allocate(A%ralpha_p(0:A%degree))
            A%ralpha_p = compute_alpha( A%degree, H )

      end select

   end subroutine real_preconditioner

   subroutine complex_preconditioner( A, preconditioner, degree, foci, center, H )
!
! Create the preconditioner A. Note that this routine can also be used to make a preconditioned
! system matrix. If the preconditioner is included in the system matrix, the final IDRS solution
! should be scaled back using the subroutine SCALEBACK.
!
! All the parameters except the first are optional.
!    complex_preconditioner( A )
! returns an identity matrix.
!
! Parameters:
!    A: matrix
!    preconditioner: integer, input.
!       preconditioner = 0: no preconditioner
!       preconditioner = 1: neumann
!       preconditioner = 2: chebyshev
!       preconditioner = 3: hessenberg
!    degree: integer, input. Degree of the polynomial preconditioner
!    foci: complex, input. Foci of ellipse around spectrum.
!    center: complex, input. Center of spectrum.
!    H: complex, hessenberg matrix, input. This matrix is required for the
!             hessenberg preconditioner and in case foci or center are not specified.
!
      type(matrix), intent(inout)                 :: A
      integer, intent(in), optional               :: preconditioner, degree
      complex(kind=cp), intent(in), optional      :: foci(2), center
      complex(kind=cp)                            :: f(2), c, shift
      complex(kind=cp), intent(in), optional      :: H(:,:)

! Polynomial preconditioning?
      if ( present(preconditioner) ) then
         A%preconditioner = preconditioner
      else
         A%preconditioner = 0
      end if

! Degree of the polynomial preconditioner?
      if ( present(degree) ) then
         A%degree = degree
      else
         A%degree = 0
      end if

!
! Construct the preconditioner:
      select case(A%preconditioner)

         case(1)
! Neumann preconditioner:

            if ( present(center) ) then
               c = center - A%cseed
            elseif ( present(H) ) then
               call analyse_spectrum( H, center = c )
            else
               c = 1._cp
            end if

            allocate(A%calpha_p(0:A%degree))
            A%calpha_p = neumann( A%degree, c )

         case(2)
! Chebyshev preconditioner:

            if ( present(foci) ) then
               f = foci - A%cseed
            elseif ( present(H) ) then
               call analyse_spectrum( H, foci = f )
            else
! These defaults are god for diagonally dominant matrices with diagonal scaling
               f(1) = 0.25_cp
               f(2) = 1.75_cp
            end if

            allocate(A%calpha_p(0:A%degree))
            A%calpha_p = chebyshev( A%degree, f )

         case(3)
! Hessenberg preconditioner:

            allocate(A%calpha_p(0:A%degree))
            A%calpha_p = compute_alpha( A%degree, H )

      end select

   end subroutine complex_preconditioner

   function ms_cscaleback( A, sigma, b, y ) result(x)
!
! Scale back the solution
!
      type(matrix), intent(in)            :: A
      complex(kind=cp), intent(in)        :: y(:,:)
      complex(kind=cp), intent(in)        :: b(:)
      complex(kind=cp), intent(in)        :: sigma(:)
      complex(kind=cp)                    :: x(size(y,1),size(y,2))
      integer                             :: degree
      integer                             :: i, n_sigma, i_sigma

      degree = A%degree
      n_sigma = size(sigma)
      x = y

      if ( degree > 0 ) then
! Correct for polynomial preconditioner x = P_sigma(A)y (Chebyshev, Neumann or Hessenberg):
         do i_sigma = 1,n_sigma
            x(:,i_sigma) = A%calpha_s(degree,i_sigma)*y(:,i_sigma)
            do i = degree-1,0,-1
               x(:,i_sigma) = complex_mv( A, x(:,i_sigma) ) + A%calpha_s(i,i_sigma)*y(:,i_sigma)
            end do
         end do
      end if

! Correct for deflation:
      x = add_deflation( A, sigma, b, x )

      if ( A%shift_invert ) then
! Correct for SI-preconditioner:
!    (1-eta) = (seed-shift)/((sigma-seed)-shift )
         do i_sigma = 1,n_sigma
            x(:,i_sigma) = x(:,i_sigma)*((A%cseed - A%cshift)/( sigma(i_sigma)-A%cshift ))
         end do
      end if

! Correct for the preconditioner:
      do i_sigma = 1, n_sigma
         x(:,i_sigma) = x(:,i_sigma)/A
      end do

   end function ms_cscaleback

   function ms_rscaleback( A, sigma, b, y ) result(x)
!
! Scale back the solution
!
      type(matrix), intent(in)            :: A
      real(kind=rp), intent(in)           :: y(:,:)
      real(kind=rp), intent(in)           :: b(:)
      real(kind=rp), intent(in)           :: sigma(:)
      real(kind=rp)                       :: x(size(y,1),size(y,2))
      integer                             :: degree
      integer                             :: i, n_sigma, i_sigma

      degree = A%degree
      n_sigma = size(sigma)
      x = y

      if ( degree > 0 ) then
! Correct for polynomial preconditioner x = P_sigma(A)y (Chebyshev, Neumann or Hessenberg):
         do i_sigma = 1,n_sigma
            x(:,i_sigma) = A%ralpha_s(degree,i_sigma)*y(:,i_sigma)
            do i = degree-1,0,-1
               x(:,i_sigma) = real_mv( A, x(:,i_sigma) ) + A%ralpha_s(i,i_sigma)*y(:,i_sigma)
            end do
         end do
      end if

! Correct for deflation:
      x = add_deflation( A, sigma, b, x )

      if ( A%shift_invert ) then
! Correct for SI-preconditioner:
         do i_sigma = 1,n_sigma
            x(:,i_sigma) = x(:,i_sigma)*((A%rseed - A%rshift)/( sigma(i_sigma)-A%rshift ))
         end do
      end if

! Correct for the preconditioner:
      do i_sigma = 1, n_sigma
         x(:,i_sigma) = x(:,i_sigma)/A
      end do
      
   end function ms_rscaleback

   function cscaleback( A, b, y ) result(x)
!
! Scale back the solution
!
      type(matrix), intent(in)               :: A
      complex(kind=cp), intent(in)           :: y(:)
      complex(kind=cp), intent(in)           :: b(:)
      complex(kind=cp)                       :: x(size(y))

! Correct for the polynomial preconditioner:
      x = polpre( A, y )

! Correct for deflation preconditioner:
      x = add_deflation( A, b, x )

! Correct for user preconditioner:
      x = x/A

   end function cscaleback

   function rscaleback( A, b, y ) result(x)
!
! Scale back the solution
!
      type(matrix), intent(in)            :: A
      real(kind=rp), intent(in)           :: y(:)
      real(kind=rp), intent(in)           :: b(:)
      real(kind=rp)                       :: x(size(y))

! Correct for the polynomial preconditioner:
      x = polpre( A, y )

! Correct for deflation preconditioner:
      x = add_deflation( A, b, x )

! Correct for preconditioner:
      x = x/A

   end function rscaleback

   function rneumann( degree, center ) result(alpha)
!
! Generate Jacobi matrix for Neumann polynomial
!
! Recursion for the residual Neumann polynomial
!         q_0 = 1
!         q_i = q_i-1 - omega t q_i-1
!
      integer, intent(in)                      :: degree
      real(kind=rp)                            :: alpha(0:degree)
      real(kind=rp), intent(in)                :: center
      real(kind=rp)                            :: J(1:degree+2,1:degree+1)
      integer                                  :: i

! Generate the Jacobi matrix
      J = 0._rp
      do i = 1, degree+1
         J(i,i) =  center
         J(i+1,i) = -center
      end do

! Compute coefficients of polynomial preconditioner:
      alpha = rcompute_alpha( degree, J )

   end function rneumann

   function cneumann( degree, center ) result(alpha)
!
! Generate Jacobi matrix for Neumann polynomial
!
! Recursion for the residual Neumann polynomial
!         q_0 = 1
!         q_i = q_i-1 - omega t q_i-1
!
      integer, intent(in)                      :: degree
      complex(kind=cp)                         :: alpha(0:degree)
      complex(kind=cp), intent(in)             :: center
      complex(kind=cp)                         :: J(1:degree+2,1:degree+1)
      integer                                  :: i

! Generate the Jacobi matrix
      J = 0._cp
      do i = 1, degree+1
         J(i,i) =    center
         J(i+1,i) = -center
      end do

! Compute coefficients of polynomial preconditioner:
      alpha = ccompute_alpha( degree, J )

   end function cneumann

   function rchebyshev( degree, foci ) result(alpha) 
!
! Generate the Jacobi matrix for the shifted and scaled Chebyshev polynomials
!
! With
!      rho_k = 1/(2sigma - rho_k-1), rho_0 = 1/sigma
! the recursion for the residual polynomial is: 
!      q_k+1 = rho_k(2(sigma-t/delta)q_k - rho_k-1 q_k-1), q_1 = 1 - t/theta, q_0 = 1
! This can be written as
!      t q_k = (-delta sigma + delta rho_k-1/2)  q_k+1 + delta sigma q_k - delta rho_k-1/2 q_k-1
!      q_0 = 1, t q_0 = theta q_0 - theta q_1
!
      integer, intent(in)                      :: degree
      real(kind=rp)                            :: alpha(0:degree)
      real(kind=rp), intent(in)                :: foci(2)
      real(kind=rp)                            :: theta, delta, sigma, rho
      real(kind=rp)                            :: J(1:degree+2,1:degree+1)
      integer                                  :: i

      theta = (foci(1)+foci(2))/2.
      delta = (foci(1)-foci(2))/2.
      sigma = theta/delta

! Generate the Jacobi matrix
      J = 0._rp
      J(1,1) =  theta
      J(2,1) = -theta
      rho = 1._rp/sigma
      do i = 2, degree+1
         J(i-1,i) =              - delta*rho/2._rp
         J(i,i)   =  delta*sigma
         J(i+1,i) = -delta*sigma + delta*rho/2._rp
         rho = 1./(2._rp*sigma - rho)
      end do

! Compute the parameters of the polynomial preconditioner
      alpha = rcompute_alpha( degree, J )

   end function rchebyshev

   function cchebyshev( degree, foci ) result(alpha) 
!
! Generate the Jacobi matrix for the shifted and scaled Chebyshev polynomials
!
! With
!      rho_k = 1/(2sigma - rho_k-1), rho_0 = 1/sigma
! the recursion for the residual polynomial is: 
!      q_k+1 = rho_k(2(rho-t/delta)q_k - rho_k-1 q_k-1), q_1 = 1 - t/theta, q_0 = 1
! This can be written as
!      t q_k = (-delta sigma + delta rho_k-1/2)  q_k+1 + delta sigma q_k - delta rho_k-1/2 q_k-1
!      q_0 = 1, t q_0 = theta q_0 - theta q_1
!
      integer, intent(in)                      :: degree
      complex(kind=cp)                         :: alpha(0:degree)
      complex(kind=cp), intent(in)             :: foci(2)
      complex(kind=cp)                         :: theta, delta, sigma, rho
      complex(kind=cp)                         :: J(1:degree+2,1:degree+1)
      integer                                  :: i

      theta = (foci(1)+foci(2))/2.
      delta = (foci(1)-foci(2))/2.
      sigma = theta/delta

! Generate the Jacobi matrix
      J = 0._cp
      J(1,1) =  theta
      J(2,1) = -theta
      rho = 1._cp/sigma
      do i = 2, degree+1
         J(i-1,i) =              - delta*rho/2._cp
         J(i,i)   =  delta*sigma
         J(i+1,i) = -delta*sigma + delta*rho/2._cp
         rho = 1./(2._cp*sigma - rho)
      end do

! Compute the parameters of the polynomial preconditioner
      alpha = ccompute_alpha( degree, J )

   end function cchebyshev

   function rcompute_alpha( degree, H ) result(alpha)
!
! Compute the coeficients of the polynomial preconditioner
! from (partial) Hessenberg decomposition AV = VH of the matrix.
!

      integer, intent(in)                         :: degree
      real(kind=rp)                               :: alpha(0:degree)
      real(kind=rp)                               :: H(:,:)
      real(kind=rp)                               :: r(0:degree+1,0:degree+1);
      integer                                     :: i, j

      if ( degree == 0 ) then
         alpha(0) = 1._rp
      else

         if ( size(H,1) < degree+2 ) error stop "Hessenberg matrix too small for polynomial degree"
         if ( size(H,2) < degree+1 ) error stop "Hessenberg matrix too small for polynomial degree"

! Compute the parameters of the residual polynomial
         r = 0._rp
         r(0,0) = 1._rp
         do i = 1, degree+1
            r(i,1:degree+1) = r(i-1,0:degree) ! Polynomial becomes one degree higher r_i = t*r_i-1
            do j = 0, i-1
               r(i,:) = r(i,:) - H(j+1,i)*r(j,:)
            end do
            r(i,:) = r(i,:)/H(i+1,i)
         end do

! Now compute the parameters of the polynomial preconditioner
! Note that r = 1 - tp(t), therefore p corresponds to minus coefficient 1:degree+1 of r
         alpha = -r(degree+1,1:degree+1)

      end if

   end function rcompute_alpha

   function ccompute_alpha( degree, H ) result(alpha)
!
! Compute the coeficients of the polynomial preconditioner
! from (partial) Hessenberg decomposition AV = VH of the matrix.
!

      integer, intent(in)                         :: degree
      complex(kind=cp)                            :: alpha(0:degree)
      complex(kind=cp)                            :: H(:,:)
      complex(kind=cp)                            :: r(0:degree+1,0:degree+1);
      integer                                     :: i, j

      if ( degree == 0 ) then
         alpha(0) = 1._cp
      else

         if ( size(H,1) < degree+2 ) error stop "Hessenberg matrix too small for polynomial degree"
         if ( size(H,2) < degree+1 ) error stop "Hessenberg matrix too small for polynomial degree"

! Compute the parameters of the residual polynomial
         r = 0._cp
         r(0,0) = 1._cp
         do i = 1, degree+1
            r(i,1:degree+1) = r(i-1,0:degree) ! Polynomial becomes one degree higher r_i = t*r_i-1
            do j = 0, i-1
               r(i,:) = r(i,:) - H(j+1,i)*r(j,:)
            end do
            r(i,:) = r(i,:)/H(i+1,i)
         end do

! Now compute the parameters of the polynomial preconditioner
! Note that r = 1 - tp(t), therefore p corresponds to minus coefficient 1:degree+1 of r
         alpha = -r(degree+1,1:degree+1)

      end if

   end function ccompute_alpha

   function filter_ritzval( ritzval ) result(outlier)
!
! Determine outliers in the Ritz values
!    In:   ritzval, complex, size n      Ritz values
!    Out:  outlier, logical, size n      True of ritzvalue is outlier
!  
      complex(kind=cp), intent(in)      :: ritzval(:)
      logical                           :: outlier(size(ritzval))
      complex(kind=cp)                  :: center
      real(kind=rp)                     :: sigma
      integer                           :: n

      n = size(ritzval)
      outlier = .true.
      center = sum(ritzval)/n
      sigma = sqrt(sum(abs(ritzval-center)**2)/n)
      outlier = (abs(ritzval - center) > 1.5_rp*sigma)

   end function filter_ritzval

   subroutine ranalyse_spectrum( H, foci, center, shift )
!
! This subroutine computes:
!   - center of circle around the spectrum of H. The center is shifted
!     to make sure that the origin is not in the circle. The center
!     of the circle is used to construct the Neumann preconditioner
!   - shift needed to avoid that the circle contains the origin.
!     This shift s used by the shift-and-invert preconditioner
!   - foci of the ellipse around the spectrum of H. The foci are shifted
!     to make sure that the origin is not in the ellipse. The foci
!     of the ellipse are used to construct the Chebyshev preconditioner.
!
! Method:
! First compute the eigenvalues of H (upper square part).
! Then compute a rectangle around the eigenvalues.
! Determine its center and its size.
!
! Determine radius of the circle around the box. Shift the circle
! of the distance of the center of the box to the origin is smaller than
! the radius. Determine the required shift.
!
! Compute the lengths of the semi-axes a and b.
! The distance of the foci to 0 is then given by
!    c = sqrt(a^2-b^2) if a > b (a is semi major axis)
! or
!    c = sqrt(b^2-a^2) if b > a (b is semi major axis)
! The foci of the unshifted ellipse are center+c and center-c
! if a is the semi major axis, and center+ic and center-ic
! if b is the semi-major axis.
! Shift the ellipse if it contains the origin.
!
! Since this is the real version, we make sure that the output parameters are real.
!
      real(kind=rp), intent(in)            :: H(:,:)
      real(kind=rp), intent(out), optional :: foci(2), center, shift
      real(kind=rp)                        :: foci_ellipse(2), center_ellipse, &
                                              center_circle, shift_circle, center_box
      real(kind=rp)                        :: radius
      complex(kind=cp)                     :: ritzval(size(H,2))
      logical                              :: outlier(size(H,2))

      real(kind=rp)                        :: x0, y0, left, right, top, bottom
      real(kind=rp)                        :: a, b, c
      real(kind=rp)                        :: correction

      integer                              :: converged
      integer                              :: i, n

! Compute the eigenvalues:
      n = size(H,2)
      ritzval = QR_eig( H )
      outlier = filter_ritzval( ritzval )

! Compute the box:
      left       = minval( real(ritzval), .not. outlier)
      right      = maxval( real(ritzval), .not. outlier)
      bottom     = minval(aimag(ritzval), .not. outlier)
      top        = maxval(aimag(ritzval), .not. outlier)
      center_box = (right+left)/2._rp

! Compute the circle:
      radius = maxval( abs(ritzval-center_box), .not. outlier)
      center_circle = center_box
      if ( abs(center_box)/radius  < 1 ) then
         correction = abs(center_box)/radius
         center_circle = center_box/correction
      end if
      shift_circle = center_box - center_circle

! Compute the ellipse:

! Following formulae compute semi-axes of ellipse centered at 0 around a box
! with corners (pm x0,pm y0) with the constraint x0/y0 = (a/b)^2
!
      x0 = (right - left)/2.
      y0 = (top - bottom)/2.
      a = sqrt(x0*(x0+y0))
      b = sqrt(y0*(x0+y0))

! Move the ellipse such that it does not contain the origin
      center_ellipse = center_box
      if ( abs(center_box)/a < 1._rp ) then
         correction = abs(center_box)/a
         center_ellipse = center_box/correction
      end if

! Now determine the foci:
      if ( a > b ) then
! a is major axis
         c = sqrt(a**2 - b**2)
         foci_ellipse(1) = center_ellipse - c
         foci_ellipse(2) = center_ellipse + c
      else
! b is major axis, take real foci around center
         c = 0.15_rp*center_ellipse
         foci_ellipse(1) = center_ellipse - c
         foci_ellipse(2) = center_ellipse + c
      end if
      if ( present(foci) )   foci = foci_ellipse
      if ( present(center) ) center = center_circle
      if ( present(shift) )  shift = shift_circle

   end subroutine ranalyse_spectrum

   subroutine canalyse_spectrum( H, foci, center, shift )
!
! This subroutine computes:
!   - center of circle around the spectrum of H. The center is shifted 
!     to make sure that the origin is not in the circle. The center
!     of the circle is used to construct the Neumann preconditioner
!   - shift needed to avoid that the circle contains the origin.
!     This shift s used by the shift-and-invert preconditioner
!   - foci of the ellipse around the spectrum of H. The foci are shifted
!     to make sure that the origin is not in the ellipse. The foci
!     of the ellipse are used to construct the Chebyshev preconditioner.
!
! Method:
! First compute the eigenvalues of H (upper square part).
! Then compute a rectangle around the eigenvalues.
! Determine its center and its size. 
!
! Determine radius of the circle around the box. Shift the circle
! of the distance of the center of the box to the origin is smaller than
! the radius. Determine the required shift.
!
! Compute the lengths of the semi-axes a and b.
! The distance of the foci to 0 is then given by
!    c = sqrt(a^2-b^2) if a > b (a is semi major axis)
! or
!    c = sqrt(b^2-a^2) if b > a (b is semi major axis)
! The foci of the unshifted ellipse are center+c and center-c
! if a is the semi major axis, and center+ic and center-ic 
! if b is the semi-major axis.
! Shift the ellipse if it contains the origin.
!
      implicit none
      complex(kind=cp), intent(in)            :: H(:,:)
      complex(kind=cp), intent(out), optional :: foci(2), center, shift
      complex(kind=cp)                        :: foci_ellipse(2), center_ellipse, &
                                                 center_circle, shift_circle, center_box
      real(kind=rp)                           :: radius
      complex(kind=cp)                        :: ritzval(size(H,2))
      logical                                 :: outlier(size(H,2))

      real(kind=rp)                           :: x0, y0, left, right, top, bottom
      real(kind=rp)                           :: a, b, c
      real(kind=rp)                           :: correction
      integer                                 :: converged

! Compute the eigenvalues:
      ritzval = QR_eig( H )
      outlier = filter_ritzval( ritzval )

! Compute the box:
      left       = minval( real(ritzval), .not. outlier )
      right      = maxval( real(ritzval), .not. outlier)
      bottom     = minval(aimag(ritzval), .not. outlier)
      top        = maxval(aimag(ritzval), .not. outlier)
      center_box = cmplx(right+left,top+bottom,kind=cp)/2._cp

! Compute the circle:
      radius = maxval( abs(ritzval-center_box), .not. outlier)
      center_circle = center_box
      if ( abs(center_box)/radius  < 1 ) then
         correction = abs(center_box)/radius
         center_circle = center_box/correction
      end if
      shift_circle = center_box - center_circle

! Compute the ellipse:

! Following formulae compute semi-axes of ellipse centered at 0 around a box
! with corners (pm x0,pm y0) with the constraint x0/y0 = (a/b)^2
!
      x0 = (right - left)/2.
      y0 = (top - bottom)/2.
      a = sqrt(x0*(x0+y0))
      b = sqrt(y0*(x0+y0))

! Move the ellipse such that it does not contain the origin
      center_ellipse = center_box
      if ( (real(center_box)/a)**2 + (aimag(center_box)/b)**2 < 1._rp ) then
         correction = sqrt((real(center_box)/a)**2 + (aimag(center_box)/b)**2)
         center_ellipse = center_box/correction
      end if

! Now determine the foci:
      if ( a > b ) then
! a is major axis
         c = sqrt(a**2 - b**2)
         foci_ellipse(1) = center_ellipse - c
         foci_ellipse(2) = center_ellipse + c
      else
! b is major axis
         c = sqrt(b**2 - a**2)
         foci_ellipse(1) = center_ellipse - cmplx(0,c,kind=cp)
         foci_ellipse(2) = center_ellipse + cmplx(0,c,kind=cp)
      end if

      if ( present(foci) )   foci = foci_ellipse
      if ( present(center) ) center = center_circle
      if ( present(shift) )  shift = shift_circle

   end subroutine canalyse_spectrum

   subroutine rmake_deflation( A, V )
!
! Compute real deflation preconditioner P = I - AV(VtAV)^-1V^T
!
      type(matrix), intent(inout)       :: A
      real(kind=rp), intent(in)         :: V(:,:)
      integer                           :: nrows, n_deflation, i_deflation, j_deflation
!
      nrows = size(V,1)
      n_deflation = size(V,2)
      A%n_deflation = n_deflation

! Allocated the space for the deflation preconditioner:
      if ( allocated(A%rV) )    deallocate(A%rV)
      if ( allocated(A%rAV) )   deallocate(A%rAV)
      if ( allocated(A%rVtAV) ) deallocate(A%rVtAV)
      allocate(A%rV(nrows,n_deflation),  &
               A%rAV(nrows,n_deflation),  &
               A%rVtAV(n_deflation,n_deflation))

! Store deflation vectors:
      A%rV = V

! Compute AV
      do i_deflation=1,n_deflation
         A%rAV(:,i_deflation) = real_mv(A,V(:,i_deflation))
      end do

! Compute V^TAV:
      do j_deflation=1,n_deflation
         do i_deflation=1,n_deflation
            A%rVtAV(i_deflation,j_deflation) = dot_product(V(:,i_deflation),A%rAV(:,j_deflation))
         end do
      end do 
      call co_sum(A%rVtAV)
   end subroutine rmake_deflation

   subroutine cmake_deflation( A, V )
!
! Compute complex deflation preconditioner P = I - AV(VhAV)^-1V^h
!
      type(matrix), intent(inout)       :: A
      complex(kind=cp), intent(in)      :: V(:,:)
      integer                           :: nrows, n_deflation, i_deflation, j_deflation
!
      nrows       = size(V,1)
      n_deflation = size(V,2)
      A%n_deflation = n_deflation

! Allocated the space for the deflation preconditioner:
      if ( allocated(A%cV) )    deallocate(A%cV)
      if ( allocated(A%cAV) )   deallocate(A%cAV)
      if ( allocated(A%cVtAV) ) deallocate(A%cVtAV)
      allocate(A%cV(nrows,n_deflation),  &
               A%cAV(nrows,n_deflation),  &
               A%cVtAV(n_deflation,n_deflation))

! Store the deflation vectors:
      A%cV = V

! Compute AV
      do i_deflation=1,n_deflation
         A%cAV(:,i_deflation) = complex_mv(A,V(:,i_deflation))
      end do

! Compute V^TAV:
      do j_deflation=1,n_deflation 
         do i_deflation=1,n_deflation 
            A%cVtAV(i_deflation,j_deflation) = dot_product(V(:,i_deflation),A%cAV(:,j_deflation))
         end do
      end do 
      call co_sum(A%cVtAV)
   end subroutine cmake_deflation

   function cPmul( A, v ) result(w)
!
! Multiply with complex deflation preconditioner:
!
      type(matrix), intent(in)          :: A
      complex(kind=cp), intent(in)      :: v(:)
      complex(kind=cp)                  :: w(size(v))
      complex(kind=cp)                  :: s(A%n_deflation), t(A%n_deflation)
      integer                           :: i_deflation
      
      w = v
      if ( A%n_deflation > 0 ) then
         if ( A%multishift ) then
! Projector is: P = I - VV^t
            do i_deflation = 1, A%n_deflation
               t(i_deflation) = dot_product(A%cV(:,i_deflation),w)
            end do
            call co_sum(t)
            w = w - matmul(A%cV,t)
         else
! Projector is: P = I - AV(VtAV)^-1V^t 
            do i_deflation = 1, A%n_deflation
               s(i_deflation) = dot_product(A%cV(:,i_deflation),w)
            end do
            call co_sum(s)
            t = solve(A%cVtAV,s)
            w = w - matmul(A%cAV,t)
         end if
      end if
   end function cPmul

   function rPmul( A, v ) result(w)
!
! Multiply with real deflation preconditioner:
!
      type(matrix), intent(in)          :: A
      real(kind=rp), intent(in)         :: v(:)
      real(kind=rp)                     :: w(size(v))
      real(kind=rp)                     :: s(A%n_deflation), t(A%n_deflation)
      integer                           :: i_deflation
      
      w = v
      if ( A%n_deflation > 0 ) then
         if ( A%multishift ) then
! Projector is: P = I - VV^t
            do i_deflation = 1, A%n_deflation
               t(i_deflation) = dot_product(A%rV(:,i_deflation),w)
            end do
            call co_sum(t)
            w = w - matmul(A%rV,t)
        else
! Projector is: P = I - AV(VtAV)^-1V^t 
            do i_deflation = 1, A%n_deflation
               s(i_deflation) = dot_product(A%rV(:,i_deflation),w)
            end do
            call co_sum(s)
            t = solve(A%rVtAV,s)
            w = w - matmul(A%rV,t)
         end if
      end if
   end function rPmul

   function radd_deflation( A, b, y ) result(x)

! Final step of the deflation procedure, add the solution
! components in the deflation space and its complement together

      implicit none
      type(matrix), intent(in)     :: A
      real(kind=rp), intent(in)    :: b(:)
      real(kind=rp), intent(in)    :: y(:)
      real(kind=rp)                :: x(size(y))
      real(kind=rp)                :: x1(size(y))
      real(kind=rp)                :: x2(size(y))
      real(kind=rp)                :: Ay(size(y))
      real(kind=rp)                :: s(A%n_deflation), t(A%n_deflation)
      integer                      :: i_deflation, j_deflation

      if ( A%n_deflation == 0 ) then
         x = y
         return
      else 

! x1 = V (V^TAV)^-1 V^T b
         do i_deflation = 1, A%n_deflation
            s(i_deflation) = dot_product(A%rV(:,i_deflation),b)
         end do
         call co_sum(s)
         t = solve(A%rVtAV,s)
         x1 = matmul(A%rV,t)

! x2 = (I - V (V^TAV)^-1 V^TA)y
         Ay = real_mv( A, y)
         do i_deflation = 1, A%n_deflation
            s(i_deflation) = dot_product(A%rV(:,i_deflation),Ay)
         end do
         call co_sum(s)
         t = solve(A%rVtAV,s)
         x2 = y - matmul(A%rV,t)

! Now add the two components together
         x = x1 + x2
      end if

   end function radd_deflation

   function cadd_deflation( A, b, y ) result(x)

! Final step of the deflation procedure, add the solution
! components in the deflation space and its complement together

      implicit none
      type(matrix), intent(in)        :: A
      complex(kind=cp), intent(in)    :: b(:)
      complex(kind=cp), intent(in)    :: y(:)
      complex(kind=cp)                :: x(size(y))
      complex(kind=cp)                :: x1(size(y))
      complex(kind=cp)                :: x2(size(y))
      complex(kind=cp)                :: Ay(size(y))
      complex(kind=cp)                :: s(A%n_deflation), t(A%n_deflation)
      integer                         :: i_deflation

      if ( A%n_deflation == 0 ) then
         x = y
      else

! x1 = V (V^HAV)^-1 V^H b
         do i_deflation = 1, A%n_deflation
            s(i_deflation) = dot_product(A%cV(:,i_deflation),b)
         end do
         call co_sum(s)
         t = solve(A%cVtAV,s)
         x1 = matmul(A%cV,t)

! x2 = (I - V (V^HAV)^-1 V^HA)y
         Ay = complex_mv( A, y)
         do i_deflation = 1, A%n_deflation
            s(i_deflation) = dot_product(A%cV(:,i_deflation),Ay)
         end do
         call co_sum(s)
         t = solve(A%cVtAV,s)
         x2 = y - matmul(A%cV,t)

! Now add the two components together
         x = x1 + x2
      end if

   end function cadd_deflation

   function rms_add_deflation( A, sigma, b, y ) result(x)

! Final step of the deflation procedure, add the solution
! components in the deflation space and its complement together

      implicit none
      type(matrix), intent(in)     :: A
      real(kind=rp), intent(in)    :: b(:)
      real(kind=rp), intent(in)    :: sigma(:)
      real(kind=rp), intent(in)    :: y(:,:)
      real(kind=rp)                :: x(size(y,1),size(y,2))
      real(kind=rp)                :: x1(size(y,1))
      real(kind=rp)                :: x2(size(y,1))
      real(kind=rp)                :: sigma_s(size(sigma))
      real(kind=rp)                :: Ay(size(y,1))
      real(kind=rp)                :: s(A%n_deflation), t(A%n_deflation)
      real(kind=rp)                :: Hs(A%n_deflation,A%n_deflation)
      integer                      :: i_deflation, j_deflation
      integer                      :: n_sigma, i_sigma, i

      if ( A%n_deflation == 0 ) then
         x = y
      else
         n_sigma = size(sigma)
         if ( A%shift_invert ) then
            sigma_s = ( sigma-A%rseed )/( sigma-A%rshift )
         else
            sigma_s = sigma-A%rseed
         end if

! x1 = V (V^T(A-sI)V)^-1 V^T b
         do i_deflation = 1, A%n_deflation
            s(i_deflation) = dot_product(A%rV(:,i_deflation),b)
         end do
         call co_sum(s)
         do i_sigma = 1, n_sigma
            Hs = A%rVtAV
            do i = 1, A%n_deflation
               Hs(i,i) = Hs(i,i) - sigma_s(i_sigma)
            end do
            t = solve(Hs,s)
            x1 = matmul(A%rV,t)

! x2 = (I - V (V^T(A-sI)V)^-1 V^T(A-sI))y
            Ay = real_mv( A, y(:,i_sigma) ) - sigma_s(i_sigma)*y(:,i_sigma)
            do i_deflation = 1, A%n_deflation
               t(i_deflation) = dot_product(A%rV(:,i_deflation),Ay)
            end do
            call co_sum(t)
            t = solve(Hs,t)
            x2 = y(:,i_sigma) - matmul(A%rV,t)

! Now add the two components together
            x(:,i_sigma) = x1 + x2
         end do
      end if

   end function rms_add_deflation

   function cms_add_deflation( A, sigma, b, y ) result(x)

! Final step of the deflation procedure, add the solution
! components in the deflation space and its complement together

      implicit none
      type(matrix), intent(in)     :: A
      complex(kind=cp), intent(in) :: b(:)
      complex(kind=cp), intent(in) :: sigma(:)
      complex(kind=cp), intent(in) :: y(:,:)
      complex(kind=cp)             :: x(size(y,1),size(y,2))
      complex(kind=cp)             :: x1(size(y,1))
      complex(kind=cp)             :: x2(size(y,1))
      complex(kind=cp)             :: sigma_s(size(sigma))
      complex(kind=cp)             :: Ay(size(y,1))
      complex(kind=cp)             :: s(A%n_deflation), t(A%n_deflation)
      complex(kind=cp)             :: Hs(A%n_deflation,A%n_deflation)
      integer                      :: i_deflation, j_deflation
      integer                      :: n_sigma, i_sigma, i

      if ( A%n_deflation == 0 ) then
         x = y
      else
         n_sigma = size(sigma)
         if ( A%shift_invert ) then
            sigma_s = ( sigma-A%cseed )/( sigma-A%cshift )
         else
            sigma_s = sigma-A%cseed
         end if

! x1 = V (V^T(A-sI)V)^-1 V^T b
         do i_deflation = 1, A%n_deflation
            s(i_deflation) = dot_product(A%cV(:,i_deflation),b)
         end do
         call co_sum(s)
         do i_sigma = 1, n_sigma
            Hs = A%cVtAV
            do i = 1, A%n_deflation
               Hs(i,i) = Hs(i,i) - sigma_s(i_sigma)
            end do
            t = solve(Hs,s)
            x1 = matmul(A%cV,t)

! x2 = (I - V (V^T(A-sI)V)^-1 V^T(A-sI))y
            Ay = complex_mv( A, y(:,i_sigma) ) - sigma_s(i_sigma)*y(:,i_sigma)
            do i_deflation = 1, A%n_deflation
               t(i_deflation) = dot_product(A%cV(:,i_deflation),Ay)
            end do
            call co_sum(t)
            t = solve(Hs,t)
            x2 = y(:,i_sigma) - matmul(A%cV,t)

! Now add the two components together
            x(:,i_sigma) = x1 + x2
         end do
      end if

   end function cms_add_deflation

end module matrix_module

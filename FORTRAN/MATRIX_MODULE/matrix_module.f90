!
! Matrix module
!   This module provides routines for the matrix-vector multiplication
!   and preconditioning operation for use in the IDRS module. The preconditioners
!   implemented here are polynomial preconditioners combined with diagonal scaling.
!
!   The polynomial preconditioners can be used both for solving standard 
!   linear systems Ax = b and for sequences of shifted linear systems (A - sI) x = b. 
!
!   A diagonal preconditioner can be used in combination with the polynomial 
!   preconditioner. In the standard case the scaled system
!        AD^-1 y = b, x = D^-1y is solved.
!   For shifted problems the problem
!        (AD^-1 - sI) y = b, x = D^-1y
!   is solved. Here D plays the role of a (diagonal) mass matrix, which
!   allows for solving problems of the form
!        (A - sD) x = b
!
!   The only operation that has to be supplied by a user is the matrix-vector multiplication.
!   This operation is defined in the module user_module. Multiplication with matrices in the standard 
!   formats DENSE, CRS, and COO are supplied in the accompanying file user_module.f90.
!
!   The matrix structure for the system matrix is set with the call
!      M = REAL_MATRIX( A )
!   for a real matrix, or
!      M = COMPLEX_MATRIX( A ) 
!   for a complex matrix.
!   The parameters are:
!      A:        User defined matrix of type user_matrix. Input.
!      M:        Matrix of type matrix_type. Output. This matrix type is used by the
!                IDRS-routines.
!
!   The parameters for the polynomial preconditioner are set by one of the following calls:
!      CALL REAL_PRECONDITIONER( M, A, PRECONDITIONER, DEGREE, FOCI, CENTER, SEED, H, SIGMA, SIGMA_P )
!   for real matrices, and 
!      CALL COMPLEX_PRECONDITIONER( M, A, PRECONDITIONER, DEGREE, FOCI, CENTER, SEED, H, SIGMA, SIGMA_P )
!   for complex matrices.
!   See the descritption in the subroutines below for further explanation.
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2024 Martin van Gijzen
!

module matrix_module

   use precision_module
   use eigen_module
   use user_module

   implicit none

   type :: matrix

! Polynomial preconditioner:
      integer                        :: degree = 0
      real(kind=rp), allocatable     :: ralpha_p(:)
      real(kind=rp), allocatable     :: ralpha_s(:,:)
      complex(kind=cp), allocatable  :: calpha_p(:)
      complex(kind=cp), allocatable  :: calpha_s(:,:)

! Diagonal scaling:
      real(kind=rp), allocatable     :: rdiag(:)
      complex(kind=cp), allocatable  :: cdiag(:)

! User defined matrix:
      type(user_matrix)              :: A

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
      module procedure CSCALEBACK, RSCALEBACK
   end interface

   interface COMPUTE_ALPHA
      module procedure CCOMPUTE_ALPHA, RCOMPUTE_ALPHA
   end interface

   interface ADD_MSPRECONDITIONER
      module procedure ADD_COMPLEX_MSPRECONDITIONER, ADD_REAL_MSPRECONDITIONER
   end interface

   interface FOCI_SPECTRUM
      module procedure CFOCI_SPECTRUM, RFOCI_SPECTRUM
   end interface

   interface CENTER_SPECTRUM
      module procedure CCENTER_SPECTRUM, RCENTER_SPECTRUM
   end interface

contains

   function cmatvec( A, v ) result(w) 

      type(matrix), intent(in)          :: A
      complex(kind=cp), intent(in)      :: v(:)
      complex(kind=cp)                  :: w(size(v))

! Note: preconditioned matrix-vector operation if 
! A contains preconditioner. 
! If not, w = A*v, so just matrix-vector multiplication
      w = A%A*(v/A)

   end function cmatvec

   function rmatvec( A, v ) result(w)

      type(matrix), intent(in)          :: A
      real(kind=rp), intent(in)         :: v(:)
      real(kind=rp)                     :: w(size(v))

! Note: preconditioned matrix-vector operation if 
! A contains preconditioner. 
! If not, w = A*v, so just matrix-vector multiplication
      w = A%A*(v/A)

   end function rmatvec

   function rprecon( v, M ) result(w)

      type(matrix), intent(in)              :: M
      real(kind=rp), intent(in)             :: v(:)
      real(kind=rp)                         :: w(size(v))
      integer                               :: i

      if ( allocated( M%rdiag ) ) then
! Polynomial Preconditioner with diagonal scaling:
         if ( M%degree > 0 ) then 
            w = M%ralpha_p(M%degree)*v
            do i = M%degree-1, 0, -1
               w = M%A*(w/M%rdiag) + M%ralpha_p(i)*v
            end do
            w = w/M%rdiag
         else
! Only diagonal scaling:
            w = v/M%rdiag
         end if
      else
! Polynomial Preconditioner:
         if ( M%degree > 0 ) then 
            w = M%ralpha_p(M%degree)*v
            do i = M%degree-1, 0, -1
               w = M%A*w + M%ralpha_p(i)*v
            end do
         else
! No preconditioner:
            w = v
         end if
      end if

   end function rprecon

   function cprecon( v, M ) result(w)

      type(matrix), intent(in)              :: M
      complex(kind=cp), intent(in)          :: v(:)
      complex(kind=cp)                      :: w(size(v))
      integer                               :: i

      if ( allocated(M%cdiag) ) then
         if ( M%degree > 0 ) then 
! Polynomial Preconditioner with diagonal scaling:
            w = M%calpha_p(M%degree)*v  
            do i = M%degree-1, 0, -1
               w = M%A*(w/M%cdiag) + M%calpha_p(i)*v
            end do
            w = w/M%cdiag
         else
! Only diagonal scaling
            w = v/M%cdiag
         end if
      else
         if ( M%degree > 0 ) then 
! Polynomial Preconditioner:
            w = M%calpha_p(M%degree)*v  
            do i = M%degree-1, 0, -1
               w = M%A*w + M%calpha_p(i)*v
            end do
         else
! No Preconditioner:
            w = v
         end if
      end if

   end function cprecon

   function real_matrix( A ) result(M)
!
! Construct the matrix M
!
      type(matrix)                                :: M
      type(user_matrix), intent(in)               :: A

!
! Put user matrix in matrix type
      M%A = A
!
      M%degree = 0
      if ( allocated(M%ralpha_p) ) deallocate(M%ralpha_p)
      if ( allocated(M%ralpha_s) ) deallocate(M%ralpha_s)

   end function real_matrix

   function complex_matrix( A ) result(M)
!
! Construct a complex matrix M
!
      type(matrix)                                :: M
      type(user_matrix), intent(in)               :: A

!
! Put user matrix in matrix type
       M%A = A

!
      M%degree = 0
      if ( allocated(M%calpha_p) ) deallocate(M%calpha_p)
      if ( allocated(M%calpha_s) ) deallocate(M%calpha_s)

   end function complex_matrix

   subroutine add_complex_mspreconditioner( M, sigma, sigma_p )
!
! Add the multishift preconditioner.
! Note: the polynomial preconditioner must already have been constructed
!
! Parameters: 
!   M: matrix, the parameters of the multishift precondiationer are aded to M
!   sigma: complex, input. The shifts
!   sigma_p: complex, output. The shifts of the preconditioned system
!
      type(matrix), intent(inout)                 :: M
      complex(kind=cp), intent(in)                :: sigma(:)
      complex(kind=cp), intent(out)               :: sigma_p(:)
      integer                                     :: i, i_sigma, n_sigma

! Multishift polynomial preconditioner?
      if ( M%degree == 0 ) then
         sigma_p = sigma
      else
         if ( .not. allocated(M%calpha_p) ) stop 'First call subroutine add_preconditioner '
! Allocate space for the parameters for the shifted polynomials
         n_sigma = size(sigma)
         if ( allocated(M%calpha_s) ) deallocate(M%calpha_s)
         allocate(M%calpha_s(0:M%degree,n_sigma))

! Compute the parameters of the shifted systems:
         M%calpha_s = 0.
         do i_sigma = 1,n_sigma
            M%calpha_s(M%degree,i_sigma) = M%calpha_p(M%degree)
            do i = M%degree-1,0,-1
               M%calpha_s(i,i_sigma) = M%calpha_p(i) + sigma(i_sigma)*M%calpha_s(i+1,i_sigma)
            end do
            sigma_p(i_sigma) = sigma(i_sigma)*M%calpha_s(0,i_sigma)
         end do
      end if

   end subroutine add_complex_mspreconditioner

   subroutine add_real_mspreconditioner( M, sigma, sigma_p )
!
! Add the multishift preconditioner.
! Note: the polynomial preconditioner must already have been constructed
!
! Parameters: 
!   M: matrix, the parameters of the multishift preconditioner are added to M
!   sigma: complex, input. The shifts
!   sigma_p: complex, output. The shifts of the preconditioned system
!
      type(matrix), intent(inout)                 :: M
      real(kind=rp), intent(in)                   :: sigma(:)
      real(kind=rp), intent(out)                  :: sigma_p(:)
      integer                                     :: i, i_sigma, n_sigma

! Multishift polynomial preconditioner?
      if ( M%degree == 0 ) then
         sigma_p = sigma
      else
         if ( .not. allocated(M%ralpha_p) ) stop 'First call subroutine add_preconditioner '
! Allocate space for the parameters for the shifted polynomials
         n_sigma = size(sigma)
         if ( allocated(M%ralpha_s) ) deallocate(M%ralpha_s)
         allocate(M%ralpha_s(0:M%degree,n_sigma))

! Compute the parameters of the shifted systems:
         M%ralpha_s = 0.
         do i_sigma = 1,n_sigma
            M%ralpha_s(M%degree,i_sigma) = M%ralpha_p(M%degree)
            do i = M%degree-1,0,-1
               M%ralpha_s(i,i_sigma) = M%ralpha_p(i) + sigma(i_sigma)*M%ralpha_s(i+1,i_sigma)
            end do
            sigma_p(i_sigma) = sigma(i_sigma)*M%ralpha_s(0,i_sigma)
         end do
      end if

   end subroutine add_real_mspreconditioner

   subroutine real_preconditioner( M, A, &
      D, preconditioner, degree, foci, center, seed, H, sigma, sigma_p )
!
! Create the preconditioner M. Note that this routine can also be used to make a preconditioned
! system matrix. If the preconditioner is included in the system matrix, the final IDRS solution
! should be scaled back using the subroutine SCALEBACK.
!
! All the parameters are optional except the first are optional.
!    real_preconditioner( M )
! returns an identity matrix.
!
! Parameters:
!    M: matrix
!    A: user_matrix, input. The system matrix as defined by the user. 
!    D: real array, input. Main diagonal of the system matrix, or mass matrix for multishift problem.
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
!    sigma: real, input. Shifts of the multishift problem.
!    sigma_p: real, output. Shifts of the preconditioned multishift problem.
!
      type(matrix), intent(inout)                 :: M
      type(user_matrix), intent(in), optional     :: A
      real(kind=rp), optional                     :: D(:)
      integer, intent(in), optional               :: preconditioner, degree
      integer                                     :: pp, dp
      real(kind=rp), intent(in), optional         :: foci(2), center
      real(kind=rp)                               :: f(2), c
      real(kind=rp), intent(in), optional         :: seed
      real(kind=rp)                               :: seed_shift
      real(kind=rp), intent(in), optional         :: H(:,:)
      real(kind=rp), intent(in), optional         :: sigma(:)
      real(kind=rp), intent(out), optional        :: sigma_p(:)
      logical                                     :: multishift
      integer                                     :: Nsigma

! First add the user matrix:
      if ( present(A) ) M%A = A

! Diagonal scaling?
      if ( present(D) ) M%rdiag = D

! Polynomial preconditioning?
      if ( present(preconditioner) ) then
         pp = preconditioner
      else 
         pp = 0
      end if

! Degree of the polynomial preconditioner?
      if ( present(degree) ) then
         dp = degree
      else 
         dp = 0
      end if

      multishift = present(sigma) .and. present(sigma_p)
      if ( multishift ) Nsigma = size(sigma)

      if ( present(seed) ) then
         seed_shift = seed
      elseif ( multishift ) then
         seed_shift = sigma(Nsigma)
      else
         seed_shift = 0.
      end if

! Which polynomial preconditioner?
      select case(pp)
         case(1)
! Neumann preconditioner:

            if ( present(center) ) then
               c = center
            elseif ( present(H) ) then
               c = center_spectrum(H)
            else
               c = 1.
            end if

            M%degree = dp
            if ( allocated(M%ralpha_p) ) deallocate(M%ralpha_p)
            allocate(M%ralpha_p(0:dp))
            M%ralpha_p = neumann( dp, c, seed_shift )

         case(2)
! Chebyshev preconditioner:

            if ( present(foci) ) then
               f = foci
            elseif ( present(H) ) then
               f = foci_spectrum(H) 
            else
! These defaults are working well for the Stommel test problems
               f(1) = 0.25
               f(2) = 1.75
            end if

            M%degree = dp
            if ( allocated(M%ralpha_p) ) deallocate(M%ralpha_p)
            allocate(M%ralpha_p(0:dp))
            M%ralpha_p = chebyshev( dp, f, seed_shift )

         case(3)

! Hessenberg preconditioner
            M%degree = dp
            if ( allocated(M%ralpha_p) ) deallocate(M%ralpha_p)
            allocate(M%ralpha_p(0:dp))
            M%ralpha_p = compute_alpha( dp, H, seed_shift )

         case default 

! No polynomial preconditioner
            M%degree = 0
      end select

! Compute the shifts for the preconditioned multishift systems:
      if ( multishift ) call add_real_mspreconditioner( M, sigma, sigma_p )

   end subroutine real_preconditioner

   subroutine complex_preconditioner( M, A, &
      D, preconditioner, degree, foci, center, seed, H, sigma, sigma_p )
!
! Create the preconditioner M. Note that this routine can also be used to make a preconditioned
! system matrix. If the preconditioner is included in the system matrix, the final IDRS solution
! should be scaled back using the subroutine SCALEBACK.
!
! All the parameters are optional except the first are optional.
!    complex_preconditioner( M )
! returns an identity matrix.
!
! Parameters:
!    M: matrix
!    A: user_matrix, input. The system matrix as defined by the user. 
!    D: complex array, input. Main diagonal of the system matrix, or mass matrix for multishift problem.
!    preconditioner: integer, input.
!       preconditioner = 0: no preconditioner
!       preconditioner = 1: neumann
!       preconditioner = 2: chebyshev
!       preconditioner = 3: hessenberg
!    degree: integer, input. Degree of the polynomial preconditioner
!    foci: complex, input. Foci of ellipse around spectrum.
!    center: complex, input. Center of spectrum.
!    seed: complex, input. Shift to the matrix.
!    H: complex, hessenberg matrix, input. This matrix is required for the
!             hessenberg preconditioner and in case foci or center are not specified.
!    sigma: complex, input. Shifts of the multishift problem.
!    sigma_p: complex, output. Shifts of the preconditioned multishift problem.
!
      type(matrix), intent(inout)                 :: M
      type(user_matrix), optional                 :: A
      complex(kind=cp), optional                  :: D(:)
      integer, intent(in), optional               :: preconditioner, degree
      integer                                     :: pp, dp
      complex(kind=cp), intent(in), optional      :: foci(2), center
      complex(kind=cp)                            :: f(2), c
      complex(kind=cp), intent(in), optional      :: seed
      complex(kind=cp)                            :: seed_shift
      complex(kind=cp), intent(in), optional      :: H(:,:)
      complex(kind=cp), intent(in), optional      :: sigma(:)
      complex(kind=cp), intent(out), optional     :: sigma_p(:)
      logical                                     :: multishift
      integer                                     :: Nsigma

! Add the user matrix if present:
      if ( present(A) ) M%A = A

! Diagonal scaling?
      if ( present(D) ) M%cdiag = D

! Polynomial preconditioning?
      if ( present(preconditioner) ) then
         pp = preconditioner
      else
         pp = 0
      end if

! Degree of the polynomial preconditioner?
      if ( present(degree) ) then
         dp = degree
      else
         dp = 0
      end if

      multishift = present(sigma) .and. present(sigma_p) 
      if ( multishift ) Nsigma = size(sigma)

      if ( present(seed) ) then
         seed_shift = seed
      elseif ( multishift ) then
         seed_shift = (1._cp,-0.5_cp) 
      else
         seed_shift = 0.
      end if

!
! Which polynomial preconditioner?
      select case(pp)

         case(1)
! Neumann preconditioner:

            if ( present(center) ) then
               c = center
            elseif ( present(H) ) then
               c = center_spectrum(H) 
            else
               c = 1.
            end if

            M%degree = dp
            if ( allocated(M%calpha_p) ) deallocate(M%calpha_p)
            allocate(M%calpha_p(0:dp))
            M%calpha_p = neumann( dp, c, seed_shift )

         case(2)
! Chebyshev preconditioner:

            if ( present(foci) ) then
               f = foci
            elseif ( present(H) ) then
               f = foci_spectrum(H) 
            else
! These defaults are god for diagonally dominant matrices with diagonal scaling
               f(1) = 0.25
               f(2) = 1.75
            end if

            M%degree = dp
            if ( allocated(M%calpha_p) ) deallocate(M%calpha_p)
            allocate(M%calpha_p(0:dp))
            M%calpha_p = chebyshev( dp, f, seed_shift )

         case(3)
! Hessenberg preconditioner:

            M%degree = dp
            if ( allocated(M%calpha_p) ) deallocate(M%calpha_p)
            allocate(M%calpha_p(0:dp))
            M%calpha_p = compute_alpha( dp, H, seed_shift )

         case default 
! No polynomial preconditioner
            M%degree = 0
      end select

! Compute the shifts for the preconditioned multishift systems:
      if ( multishift ) call add_complex_mspreconditioner( M, sigma, sigma_p )

   end subroutine complex_preconditioner

   subroutine cscaleback( M, x ) 
!
! Scale back the solution
!
      type(matrix), intent(in)            :: M
      complex(kind=cp), intent(inout)     :: x(:,:)
      complex(kind=cp)                    :: y(size(x,1))
      integer                             :: degree
      integer                             :: i, n_sigma, i_sigma

      degree = M%degree
      n_sigma = size(x,2)

      if ( allocated(M%cdiag) ) then
         if ( degree > 0 .and. allocated(M%calpha_s) ) then
! Correct for polynomial preconditioner x = P_sigma(A)y:
            do i_sigma = 1,n_sigma
               y = x(:,i_sigma)
               x(:,i_sigma) = M%calpha_s(degree,i_sigma)*y
               do i = degree-1,0,-1
                  x(:,i_sigma) = M%A*(x(:,i_sigma)/M%cdiag) + M%calpha_s(i,i_sigma)*y
               end do
               x(:,i_sigma) = x(:,i_sigma)/M%cdiag
            end do
         elseif ( degree > 0 ) then
! Correct for Mass matrix x = P(A)y, with P(A) pol. approx. for inverse of mass matrix:
            do i_sigma = 1,n_sigma
               y = x(:,i_sigma)
               x(:,i_sigma) = M%calpha_p(M%degree)*y
               do i = M%degree-1, 0, -1
                  x(:,i_sigma) = M%A*(x(:,i_sigma)/M%cdiag) + M%calpha_p(i)*y
               end do
               x(:,i_sigma) = x(:,i_sigma)/M%cdiag
            end do
         else
            do i_sigma = 1, n_sigma
               x(:,i_sigma) = x(:,i_sigma)/M%cdiag
            end do
         end if
      else
         if ( degree > 0 .and. allocated(M%calpha_s) ) then
! Correct for polynomial preconditioner x = P_sigma(A)y:
            do i_sigma = 1,n_sigma
               y = x(:,i_sigma)
               x(:,i_sigma) = M%calpha_s(degree,i_sigma)*y
               do i = degree-1,0,-1
                  x(:,i_sigma) = M%A*x(:,i_sigma) + M%calpha_s(i,i_sigma)*y
               end do
            end do
         elseif ( degree > 0 ) then
! Correct for Mass matrix x = P(A)y, with P(A) pol. approx. for inverse of mass matrix:
            do i_sigma = 1,n_sigma
               y = x(:,i_sigma)
               x(:,i_sigma) = M%calpha_p(M%degree)*y
               do i = M%degree-1, 0, -1
                  x(:,i_sigma) = M%A*x(:,i_sigma) + M%calpha_p(i)*y
               end do
            end do
         end if
      end if

   end subroutine cscaleback

   subroutine rscaleback( M, x ) 
!
! Scale back the solution
!
      type(matrix), intent(in)            :: M
      real(kind=rp), intent(inout)        :: x(:,:)
      real(kind=rp)                       :: y(size(x,1))
      integer                             :: degree
      integer                             :: i, n_sigma, i_sigma

      n_sigma = size(x,2)
      degree = M%degree

      if ( allocated(M%rdiag) ) then
         if ( degree > 0 .and. allocated(M%ralpha_s) ) then
! Correct for polynomial preconditioner x = P_sigma(A)y:
            do i_sigma = 1,n_sigma
               y = x(:,i_sigma)
               x(:,i_sigma) = M%ralpha_s(degree,i_sigma)*y
               do i = degree-1,0,-1
                  x(:,i_sigma) = M%A*(x(:,i_sigma)/M%rdiag) + M%ralpha_s(i,i_sigma)*y
               end do
               x(:,i_sigma) = x(:,i_sigma)/M%rdiag
            end do
         elseif ( degree > 0 ) then
! Correct for Mass matrix x = P(A)y, with P(A) pol. approx. for inverse of mass matrix:
            do i_sigma = 1,n_sigma
               y = x(:,i_sigma)
               x(:,i_sigma) = M%ralpha_p(M%degree)*y
               do i = M%degree-1, 0, -1
                  x(:,i_sigma) = M%A*(x(:,i_sigma)/M%rdiag) + M%ralpha_p(i)*y
               end do
               x(:,i_sigma) = x(:,i_sigma)/M%rdiag
            end do
         else
            do i_sigma = 1, n_sigma
               x(:,i_sigma) = x(:,i_sigma)/M%rdiag
            end do
         end if
      else
         if ( degree > 0 .and. allocated(M%ralpha_s) ) then
! Correct for polynomial preconditioner x = P_sigma(A)y:
            do i_sigma = 1,n_sigma
               y = x(:,i_sigma)
               x(:,i_sigma) = M%ralpha_s(degree,i_sigma)*y
               do i = degree-1,0,-1
                  x(:,i_sigma) = M%A*x(:,i_sigma) + M%ralpha_s(i,i_sigma)*y
               end do
            end do
         elseif ( degree > 0 ) then
! Correct for Mass matrix x = P(A)y, with P(A) pol. approx. for inverse of mass matrix:
            do i_sigma = 1,n_sigma
               y = x(:,i_sigma)
               x(:,i_sigma) = M%ralpha_p(M%degree)*y
               do i = M%degree-1, 0, -1
                  x(:,i_sigma) = M%A*x(:,i_sigma) + M%ralpha_p(i)*y
               end do
            end do
         end if
      end if

   end subroutine rscaleback

   function rneumann( degree, center, seed ) result(alpha)
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
      real(kind=rp), intent(in), optional      :: seed
      real(kind=rp)                            :: shift
      real(kind=rp)                            :: J(1:degree+2,1:degree+1)
      integer                                  :: i

      shift = 0.
      if ( present(seed) ) shift = seed
      
! Generate the Jacobi matrix
      J = 0._rp
      do i = 1, degree+1
         J(i,i) =  center+shift
         J(i+1,i) = -center-shift
      end do

! Compute coefficients of polynomial preconditioner:
      alpha = rcompute_alpha( degree, J, shift )

   end function rneumann

   function cneumann( degree, center, seed ) result(alpha)
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
      complex(kind=cp), intent(in), optional   :: seed
      complex(kind=cp)                         :: shift
      complex(kind=cp)                         :: J(1:degree+2,1:degree+1)
      integer                                  :: i

      if ( present(seed) ) shift = seed

! Generate the Jacobi matrix
      J = 0._cp
      do i = 1, degree+1
         J(i,i) =  center+shift
         J(i+1,i) = -center - shift
      end do

! Compute coefficients of polynomial preconditioner:
      alpha = ccompute_alpha( degree, J, shift )

   end function cneumann

   function rchebyshev( degree, foci, seed ) result(alpha) 
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
      real(kind=rp)                            :: alpha(0:degree)
      real(kind=rp), intent(in)                :: foci(2)
      real(kind=rp), intent(in), optional      :: seed
      real(kind=rp)                            :: theta, delta, sigma, rho, shift
      real(kind=rp)                            :: J(1:degree+2,1:degree+1)
      integer                                  :: i

      shift = 0.
      if ( present(seed) ) shift = seed

      theta = (foci(1)+foci(2))/2. + shift
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
      alpha = rcompute_alpha( degree, J, shift )

   end function rchebyshev

   function cchebyshev( degree, foci, seed ) result(alpha) 
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
      complex(kind=cp), intent(in), optional   :: seed
      complex(kind=cp)                         :: theta, delta, sigma, rho, shift
      complex(kind=cp)                         :: J(1:degree+2,1:degree+1)
      integer                                  :: i

      shift = 0.
      if ( present(seed) ) shift = seed

      theta = (foci(1)+foci(2))/2. + shift
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
      alpha = ccompute_alpha( degree, J, shift )

   end function cchebyshev

   function rcompute_alpha( degree, H, shift ) result(alpha)
!
! Compute the coeficients of the polynomial preconditioner
! from (partial) Hessenberg decomposition AV = VH of the matrix.
!

      integer, intent(in)                         :: degree
      real(kind=rp)                               :: alpha(0:degree)
      real(kind=rp)                               :: H(:,:)
      real(kind=rp), intent(in)                   :: shift
      real(kind=rp)                               :: r(0:degree+1,0:degree+1);
      integer                                     :: i, j

      if ( degree == 0 ) then
         alpha(0) = 1._rp
      else

         if ( size(H,1) < degree+2 ) stop "Hessenberg matrix too small for polynomial degree"
         if ( size(H,2) < degree+1 ) stop "Hessenberg matrix too small for polynomial degree"

! Add shift to main diagonal of H
         do i = 1, degree+1
            H(i,i) = H(i,i) + shift
         end do

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

   function ccompute_alpha( degree, H, shift ) result(alpha)
!
! Compute the coeficients of the polynomial preconditioner
! from (partial) Hessenberg decomposition AV = VH of the matrix.
!

      integer, intent(in)                         :: degree
      complex(kind=cp)                            :: alpha(0:degree)
      complex(kind=cp)                            :: H(:,:)
      complex(kind=cp), intent(in)                :: shift
      complex(kind=cp)                            :: r(0:degree+1,0:degree+1);
      integer                                     :: i, j

      if ( degree == 0 ) then
         alpha(0) = 1._cp
      else

         if ( size(H,1) < degree+2 ) stop "Hessenberg matrix too small for polynomial degree"
         if ( size(H,2) < degree+1 ) stop "Hessenberg matrix too small for polynomial degree"

! Add shift to main diagonal of H
         do i = 1, degree+1
            H(i,i) = H(i,i) + shift
         end do

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

   function rfoci_spectrum( H ) result(rfoci)
!
! Compute foci of an ellipse around the eigenvalues of H.
! First compute the eigenvalues of H (upper square part).
! Then compute a rectangle around the eigenvalues, and its center. 
! Compute its center and the lengths of the semi-axes a and b.
! The distance of the foci to 0 is then given by
!    c = sqrt(a^2-b^2) if a > b (a is semi major axis)
! or
!    c = sqrt(b^2-a^2) if b > a (b is semi major axis)
! The foci of the unshifted ellipse are center+c and center-c
! if a is the semi major axis, and center+ic and center-ic 
! if b is the semi-major axis.
!
! Since this is the real version, we take the real part of the foci as output.
!
   real(kind=rp), intent(in)            :: H(:,:)
   complex(kind=cp)                     :: ritzval(size(H,2))
   complex(kind=cp)                     :: cfoci(2), center
   real(kind=rp)                        :: rfoci(2)

   real(kind=rp)                        :: x0, y0, left, right, top, bottom
   real(kind=rp)                        :: a, b, c

! Compute the eigenvalues:
   ritzval = QR_eig( H )

! Compute the box:
   left   = minval( real(ritzval))
   right  = maxval( real(ritzval))
   bottom = minval(aimag(ritzval))
   top    = maxval(aimag(ritzval))
   center = cmplx(right+left,top+bottom,kind=cp)/2.
   x0 = (right - left)/2.
   y0 = (top - bottom)/2.

! Following formulae compute semi-axes of ellipse centered at 0 around a box
! with corners (pm x0,pm y0) with the constraint x0/y0 = (a/b)^2
!
   a = sqrt(x0*(x0+y0))
   b = sqrt(y0*(x0+y0))

! Now determine the foci:
   if ( a > b ) then
! a is major axis
      c = sqrt(a**2 - b**2)
      cfoci(1) = center - c
      cfoci(2) = center + c
   else
! b is major axis
     c = sqrt(b**2 - a**2)
     cfoci(1) = center - cmplx(0,c,kind=cp)
     cfoci(2) = center + cmplx(0,c,kind=cp)
   end if

   rfoci = real(cfoci)

   end function rfoci_spectrum

   function cfoci_spectrum( H ) result(foci)
!
! Compute foci of an ellipse around the eigenvalues of H.
! First compute the eigenvalues of H (upper square part).
! Then compute a rectangle around the eigenvalues, and its center. 
! Compute its center and the lengths of the semi-axes a and b.
! The distance of the foci to 0 is then given by
!    c = sqrt(a^2-b^2) if a > b (a is semi major axis)
! or
!    c = sqrt(b^2-a^2) if b > a (b is semi major axis)
! The foci of the unshifted ellipse are center+c and center-c
! if a is the semi major axis, and center+ic and center-ic 
! if b is the semi-major axis.
!
   complex(kind=cp), intent(in)            :: H(:,:)
   complex(kind=cp)                        :: ritzval(size(H,2))
   complex(kind=cp)                        :: foci(2), center

   real(kind=rp)                           :: x0, y0, left, right, top, bottom
   real(kind=rp)                           :: a, b, c

! Compute the eigenvalues:
   ritzval = QR_eig( H )

! Compute the box:
   left   = minval( real(ritzval))
   right  = maxval( real(ritzval))
   bottom = minval(aimag(ritzval))
   top    = maxval(aimag(ritzval))
   center = cmplx(right+left,top+bottom,kind=cp)/2.
   x0 = (right - left)/2.
   y0 = (top - bottom)/2.

! Following formulae compute semi-axes of ellipse centered at 0 around a box
! with corners (pm x0,pm y0) with the constraint x0/y0 = (a/b)^2
!
   a = sqrt(x0*(x0+y0))
   b = sqrt(y0*(x0+y0))

! Now determine the foci:
   if ( a > b ) then
! a is major axis
      c = sqrt(a**2 - b**2)
      foci(1) = center - c
      foci(2) = center + c
   else
! b is major axis
     c = sqrt(b**2 - a**2)
     foci(1) = center - cmplx(0,c,kind=cp)
     foci(2) = center + cmplx(0,c,kind=cp)
   end if

   end function cfoci_spectrum

   function ccenter_spectrum( H ) result(center)
!
! Compute center of spectrum of H computed as center = trace(H)/N. 
!
   complex(kind=cp), intent(in)     :: H(:,:)
   complex(kind=cp)                 :: trace, center
   integer                          :: N, i

   N = size(H,2)
   trace = 0.
   do i = 1, N
      trace = trace + H(i,i)
   end do
   center = trace/N

   end function ccenter_spectrum

   function rcenter_spectrum( H ) result(center)
!
! Compute center of spectrum of H as center = trace(H)/N. 
!
   real(kind=rp), intent(in)        :: H(:,:)
   real(kind=rp)                    :: trace, center
   integer                          :: N, i

   N = size(H,2)
   trace = 0.
   do i = 1, N
      trace = trace + H(i,i)
   end do
   center = trace/N

   end function rcenter_spectrum

end module matrix_module

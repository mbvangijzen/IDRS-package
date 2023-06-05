!
! IDRS Induced Dimension Reduction module
!
!   IDR(s) is a family of iterative solution methods to solve (a sequency of)
!   nonsymmetric linear systems Ax = b, with A a matrix of size NxN, and x and b
!   vectors of size N. Using special choices of the parameters the classical Krylov 
!   subspace method BiCGSTAB and (flexible) GMRES can be recovered. Depending on the 
!   characteristics of the linear system, other IDR(s) variants can considerably outperform 
!   these methods.  In exact arithmatic, IDR(s) terminates in at most N+N/s iterations
!   
!   This module contains four IDR(s) variants:
!
!      IDRS: efficient implementation based on bi-orthogonalisation of the basis vectors. 
!            With special choice of parameters the method is mathematically equivalent to 
!            BiCGSTAB, and has the same storage and computational requirements.
!            For indefinite problems and/or strongly nonsymmetric problems IDRS 
!            with s > 1 is often significant faster than BiCGSTAB.
!            IDR(s) alllows for recycling of a previously computed subspace
!            to speed-up the solution of a sequence of linear systems.
!
!    QMRIDR: robust implementation based on orthogonalisation of the basis vectors.
!            As long as the number of iterations is smaller than s, the method is
!            mathematically equivalent with flexible GMRES, although slightly less efficient
!            in storage (needs more vectors) and vector operations per iteration. 
!            When the number of iterations is larger than s, QMRIDR can be viewed as 
!            a very effient restarted GMRES method, that terminates within N+N/s iterations
!            QMRIDR is a flexible method, it can handle non-constant preconditioners. 
!            IDRS can be used as an inner iterative method.
!
!    MSIDRS: Routine to solve sequence of shifted problem A-shift I = b for several values of 
!            shift. MSIDRS is an efficient implementation based on bi-orthogonalisation of the 
!            basis vectors. The resulting residuals are collinear.
!
!  MSQMRIDR: Routine to solve sequence of shifted problem A-shift I = b for several values of 
!            shift. MSQMRIDR is a robust implementation based on orthogonalisation of the 
!            basis vectors. MSQMRIDR is a flexible method. MSIDRS can be used as inner
!            iterative method.
!
!   The IDRS module requires user defined TYPE MATRIX, and TYPE PRECONDITIONER.
!   These types should be defined in a module MATRIX_MODULE. Tis module must contain functions 
!   XR = RMATVEC( A, YR ), XC = CMATVEC( A, YC ), XR = RPRECON( M, YR ) and XC = CPRECON( M, YR ),
!   for performing the matrix-vector product and preconditioning operation.
!   The vectors xr and yr must be of type real(kind=rp), and xc and yc must be of type real(kind=cp).
!   A must be of type matrix, and M of type preconditioner. 
!   The operator * must be overloaded to perform the operation A*y, and / to perform M/y.
!
!   The precision of the complex and real variables are defined by the 
!   parameters rp (real precision)and cp (complex precision). These parameters 
!   should also be set in the module matrix_module.
!
!
! This software is distributed under the MIT License:
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2023 Martin van Gijzen
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!

module idrs_module

   use matrix_module
   use preconditioner_module

   implicit none

   private
   public :: IDRS, QMRIDR, MSQMRIDR, MSIDRS, NORM

   interface IDRS
      module procedure CIDRS, RIDRS
   end interface

   interface QMRIDR
      module procedure CQMRIDR, RQMRIDR
   end interface

   interface MSQMRIDR
      module procedure CMSQMRIDR, RMSQMRIDR
   end interface

   interface MSIDRS
      module procedure CMSIDRS, RMSIDRS
   end interface

   interface P_DOT
      module procedure CP_DOT, RP_DOT
   end interface

   interface GROT
      module procedure CGROT, RGROT
   end interface

   interface MAKE_P
      module procedure RMAKE_P, CMAKE_P
   end interface

   interface SOLVE
      module procedure CSOLVE, RSOLVE
   end interface

   interface COMP_OM
      module procedure CCOMP_OM, RCOMP_OM
   end interface

   interface NORM
      module procedure CNORM, RNORM
   end interface

   interface INPROD
      module procedure CINPROD, RINPROD
   end interface

contains

   function CIDRS( A, b, M1, s, &                             ! required
                   tolerance, maximum_iterations, variant, &  ! optional input
                   flag, relres, iterations, &                ! optional output
                   x0, U0, omega, resvec, H )                 ! optional arrays
!
! CIDRS: complex version of IDRS
! Parameters:
!  A:                  type matrix, required, input                       - System matrix
!  b:                  complex(kind=cp), required, size is n, input       - Rhs of the system
!  M1:                 type preconditioner, required, input               - Preconditioner
!  s:                  integer > 0, required, input                       - Size of shadow space
!  tolerance:          real(kind=rp), optional, input, default 1e-8       - Termination criterion
!  maximum_iterations: integer, optional, input, default min(2*n,1000)    - Maximum number of iterations
!  variant:            integer, optional, input, default 1                - Variant = 1: idrs, 2: bicgstab
!  flag:               integer, optional, output                          - Convergence flags
!  relres:             real(kind=rp), optional, output                    - Relative residual norm
!  iterations:         integer, optional, output                          - Number of iterations to converge
!  x0:                 complex(kind=cp), array size n, optional, input    - Initial guess
!  U0:                 complex(kind=cp), array size n*s, optional, input  - Initial search space
!  omega:              complex(kind=cp), array, optinal, input            - User defined values for omega
!  resvec:             real(kind=rp), array size maxit+1, optional, output- Residual norm for every iteration
!  H:                  complex(kind=cp), array size (nritz+1)*nritz,      - Upper Hessenberger
!                                        optional, output                 - can be used to compute ritzvalues
!
   IMPLICIT NONE

! Required input parameters:
   type(matrix), intent(in)               :: A         
   complex(kind=cp), intent(in)           :: b(:)     
   type(preconditioner), intent(in)       :: M1      
   integer, intent(in)                    :: s      

! Solution:
   complex(kind=cp)                       :: CIDRS(size(b))
! Optional input parameters:
   real(kind=rp), optional, intent(in)    :: tolerance 
   integer, optional, intent(in)          :: maximum_iterations 
   integer, optional, intent(in)          :: variant     

! Optional output parameters:
   integer, optional, intent(out)         :: flag        
   real(kind=rp), optional, intent(out)   :: relres     
   integer, optional, intent(out)         :: iterations

! Optional input arrays:
   complex(kind=cp), optional, intent(in) :: x0(:)  
   complex(kind=cp), optional, intent(in) :: U0(:,:) 
   complex(kind=cp), optional, intent(in) :: omega(:) 

! Optional output arrays:
   real(kind=rp), optional, intent(out)   :: resvec(:) 
   complex(kind=cp), optional, intent(out):: H(:,:)   
    
! Local arrays:
   complex(kind=cp)                       :: x(size(b))
   complex(kind=cp), allocatable          :: PT(:,:) 
   complex(kind=cp)                       :: G(size(b),s)
   complex(kind=cp)                       :: U(size(b),s)
   complex(kind=cp)                       :: r(size(b)) 
   complex(kind=cp)                       :: v(size(b))   
   complex(kind=cp)                       :: t(size(b))  
   complex(kind=cp)                       :: M(s,s), f(s), mu(s)
   complex(kind=cp)                       :: alpha(s), beta(s), gamma(s)
   
   complex(kind=cp)                       :: om
   real(kind=rp)                          :: kappa
   complex(kind=cp)                       :: buffer(3)

   include 'idrs_body.f90'
   cidrs = x
   
   end function CIDRS

   function RIDRS( A, b, M1, s, &                            ! required
                   tolerance, maximum_iterations, variant, & ! optional input
                   flag, relres, iterations, &               ! optional output
                   x0, U0, omega, H, resvec )                ! optional arrays
!
! RIDRS: real version of IDRS
! Parameters:
!  A:                  type matrix, required, input                       - System matrix
!  b:                  real(kind=rp), required, size is n, input          - Rhs of the system
!  M1:                 type preconditioner, required, input               - Preconditioner
!  s:                  integer > 0, required, input                       - Size of shadow space
!  tolerance:          real(kind=rp), optional, input, default 1e-8       - Termination criterion
!  maximum_iterations: integer, optional, input, default min(2*n,1000)    - Maximum number of iterations
!  variant:            integer, optional, input, default 1                - Variant = 1: idrs, 2: bicgstab
!  flag:               integer, optional, output                          - Convergence flags
!  relres:             real(kind=rp), optional, output                    - Relative residual norm
!  iterations:         integer, optional, output                          - Number of iterations to converge
!  x0:                 real(kind=rp), array size n, optional, input       - Initial guess
!  U0:                 real(kind=rp), array size n*s, optional, input     - Initial search space
!  omega:              real(kind=rp), array, optinal, input               - User defined values for omega
!  resvec:             real(kind=rp), array size maxit+1, optional, output- Residual norm for every iteration
!  H:                  real(kind=rp), array size (nritz+1)*nritz,         - Upper Hessenberger
!                                        optional, output                 - can be used to compute ritzvalues
!

   IMPLICIT NONE

! Required input and output parameters:
   type(matrix), intent(in)               :: A             
   real(kind=rp), intent(in)              :: b(:)         
   type(preconditioner), intent(in)       :: M1          
   integer, intent(in)                    :: s          

! Solution
   real(kind=rp)                          :: RIDRS(size(b))

! Optional input parameters:
   real(kind=rp), optional, intent(in)    :: tolerance
   integer, optional, intent(in)          :: maximum_iterations
   integer, optional, intent(in)          :: variant

! Optional output parameters:
   integer, optional, intent(out)         :: flag
   real(kind=rp), optional, intent(out)   :: relres
   integer, optional, intent(out)         :: iterations 

! Optional input arrays:
   real(kind=rp), optional, intent(in)    :: x0(:)
   real(kind=rp), optional, intent(in)    :: U0(:,:)
   real(kind=rp), optional, intent(in)    :: omega(:)

! Optional output arrays
   real(kind=rp), optional, intent(out)   :: resvec(:)
   real(kind=cp), optional, intent(out)   :: H(:,:)  
    
! Local arrays:
   real(kind=rp)                          :: x(size(b))
   real(kind=rp), allocatable             :: PT(:,:) 
   real(kind=rp)                          :: G(size(b),s)
   real(kind=rp)                          :: U(size(b),s)
   real(kind=rp)                          :: r(size(b))            
   real(kind=rp)                          :: v(size(b))         
   real(kind=rp)                          :: t(size(b))        
   real(kind=rp)                          :: M(s,s), f(s), mu(s)
   real(kind=rp)                          :: alpha(s), beta(s), gamma(s) 
    
   real(kind=rp)                          :: om, kappa
   real(kind=rp)                          :: buffer(3)
    
   include 'idrs_body.f90'
   RIDRS = x
     
   end function RIDRS

   function CQMRIDR( A, b, M1, s, &                                ! required
                     tolerance, maximum_iterations, &              ! optional input
                     flag, relres, iterations, &                   ! optional output
                     inner_s, inner_tolerance, inner_iterations, & ! optional input                    
                     x0, omega, resvec )                           ! optional arrays
!
! CQMRIDR: complex version of Quasi-minimal residual IDRS
! Parameters:
!  A:                  type matrix, required, input                       - System matrix
!  b:                  complex(kind=cp), required, size is n, input       - Rhs of the system
!  M1:                 type preconditioner, required, input               - Preconditioner
!  s:                  integer > 0, required, input                       - Size of shadow space
!  tolerance:          real(kind=rp), optional, input, default 1e-8       - Termination criterion
!  maximum_iterations: integer, optional, input, default min(2*n,1000)    - Maximum number of iterations
!  flag:               integer, optional, output                          - Convergence flags
!  relres:             real(kind=rp), optional, output                    - Relative residual norm
!  iterations:         integer, optional, output                          - Number of iterations to converge
!  inner_s             integer, optional, input, default 4                - Parameter s for inner iterations
!  inner_tolerance     real(kind-rp), optional, input, default = 1e-1     - Tolerance for inner iterations
!  inner_iterations    integer, optional, input, default 0                - Maximum number of inner iterations
!  x0:                 complex(kind=cp), array size n, optional, input    - Initial guess
!  omega:              complex(kind=cp), array, optional, input           - User defined values for omega
!  resvec:             real(kind=rp), array size maxit+1, optional, output- Residual norm for every iteration
!

   IMPLICIT NONE

! Required input parameters:
   type(matrix), intent(in)               :: A       ! system matrix
   complex(kind=cp), intent(in)           :: b(:)    ! system rhs
   type(preconditioner), intent(in)       :: M1      ! preconditioner
   integer, intent(in)                    :: s       ! s parameter

! Solution:
   complex(kind=cp)                       :: CQMRIDR(size(b) )

! Optional input parameters:
   real(kind=rp), optional, intent(in)    :: tolerance 
   integer, optional, intent(in)          :: maximum_iterations 

   integer, optional, intent(in)          :: inner_s
   real(kind=rp), optional, intent(in)    :: inner_tolerance
   integer, optional, intent(in)          :: inner_iterations

! Optional output parameters:
   integer, optional, intent(out)         :: flag        
   real(kind=rp), optional, intent(out)   :: relres
   integer, optional, intent(out)         :: iterations

! Optional input arrays:
   complex(kind=cp), optional, intent(in) :: x0(:)  
   complex(kind=cp), optional, intent(in) :: omega(:) 

! Optional output arrays:
   real(kind=rp), optional, intent(out)   :: resvec(:) 
    
! Local arrays:
   complex(kind=cp)                       :: x(size(b))
   complex(kind=cp), allocatable          :: PT(:,:) 
   complex(kind=cp)                       :: G(size(b),s), gn(size(b))
   complex(kind=cp)                       :: W(size(b),s+1) 
   complex(kind=cp)                       :: wn(size(b))
   complex(kind=cp)                       :: v(size(b))   
   complex(kind=cp)                       :: vtilde(size(b))
   complex(kind=cp)                       :: M(s,s), f(s), mu(s)
   complex(kind=cp)                       :: alpha(s), beta(s), gamma(s)
   complex(kind=cp)                       :: c(s+2), h(s+2), r_sigma(s+3), tau
   complex(kind=cp)                       :: cs(s+2), sn(s+2)
   complex(kind=cp)                       :: phi_n, phi_n1
   complex(kind=cp)                       :: eta
   
   complex(kind=cp)                       :: om
   real(kind=rp)                          :: kappa
   complex(kind=cp)                       :: buffer(3)

   include 'qmridr_body.f90'
   cqmridr = x
   
   end function CQMRIDR

   function RQMRIDR( A, b, M1, s, &                             ! required
                     tolerance, maximum_iterations, &           ! optional input
                     flag, relres, iterations, &                ! optional output
                     inner_s, inner_tolerance, inner_iterations, & ! optional input
                     x0, omega, resvec )                        ! optional arrays
!
! RQMRIDR: real version of Quasi-minimal resdual IDRS
! Parameters:
!  A:                  type matrix, required, input                       - System matrix
!  b:                  real(kind=rp), required, size is n, input          - Rhs of the system
!  M1:                 type preconditioner, required, input               - Preconditioner
!  s:                  integer > 0, required, input                       - Size of shadow space
!  tolerance:          real(kind=rp), optional, input, default 1e-8       - Termination criterion
!  maximum_iterations: integer, optional, input, default min(2*n,1000)    - Maximum number of iterations
!  flag:               integer, optional, output                          - Convergence flags
!  relres:             real(kind=rp), optional, output                    - Relative residual norm
!  iterations:         integer, optional, output                          - Number of iterations to converge
!  inner_s             integer, optional, input, default 4                - Parameter s for inner iterations
!  inner_tolerance     real(kind-rp), optional, input, default = 1e-1     - Tolerance for inner iterations
!  inner_iterations    integer, optional, input, default 0                - Maximum number of inner iterations
!  x0:                 real(kind=rp), array size n, optional, input       - Initial guess
!  omega:              real(kind=rp), array, optional, input              - User defined values for omega
!  resvec:             real(kind=rp), array size maxit+1, optional, output- Residual norm for every iteration
!          

   IMPLICIT NONE

! Required input parameters:
   type(matrix), intent(in)               :: A       ! system matrix
   real(kind=rp), intent(in)              :: b(:)    ! system rhs
   type(preconditioner), intent(in)       :: M1      ! preconditioner
   integer, intent(in)                    :: s       ! s parameter

! Solution:
   real(kind=rp)                          :: RQMRIDR(size(b) )

! Optional input parameters:
   real(kind=rp), optional, intent(in)    :: tolerance 
   integer, optional, intent(in)          :: maximum_iterations 

   integer, optional, intent(in)          :: inner_s
   real(kind=rp), optional, intent(in)    :: inner_tolerance
   integer, optional, intent(in)          :: inner_iterations

! Optional output parameters:
   integer, optional, intent(out)         :: flag        
   real(kind=rp), optional, intent(out)   :: relres
   integer, optional, intent(out)         :: iterations

! Optional input arrays:
   real(kind=rp), optional, intent(in)    :: x0(:)  
   real(kind=rp), optional, intent(in)    :: omega(:) 

! Optional output arrays:
   real(kind=rp), optional, intent(out)   :: resvec(:) 
    
! Local arrays:
   real(kind=rp)                          :: x(size(b))
   real(kind=rp), allocatable             :: PT(:,:) 
   real(kind=rp)                          :: G(size(b),s), gn(size(b))
   real(kind=rp)                          :: W(size(b),s+1) 
   real(kind=rp)                          :: wn(size(b))
   real(kind=rp)                          :: v(size(b))   
   real(kind=rp)                          :: vtilde(size(b))   
   real(kind=rp)                          :: V_tilde(size(b))   
   real(kind=rp)                          :: M(s,s), f(s), mu(s)
   real(kind=rp)                          :: alpha(s), beta(s), gamma(s)
   real(kind=rp)                          :: c(s+2), h(s+2), r_sigma(s+3), tau
   real(kind=rp)                          :: cs(s+2), sn(s+2)
   real(kind=rp)                          :: phi_n, phi_n1
   real(kind=rp)                          :: eta
   
    
   real(kind=rp)                          :: om, kappa
   real(kind=rp)                          :: buffer(3)

   include 'qmridr_body.f90'
   rqmridr = x
   
   end function RQMRIDR

   function CMSQMRIDR( A, b, sigma, M1, s, &                                  ! required
                       tolerance, maximum_iterations,&                    ! optional input
                       flag, relres, iterations, &                        ! optional output
                       inner_s, inner_tolerance, inner_iterations, &      ! optional input
                       omega, resvec )                                    ! optional arrays
!
! CMSQMRIDR: complex version of Quasi-minimal residual multi-shift IDRS
! Parameters:
!  A:                  type matrix, required, input                       - System matrix
!  b:                  complex(kind=cp), required, size is n, input       - Rhs of the system
!  sigma:              complex(kind=cp), array size n_sigma, required, input - Shifts
!  M1:                 type preconditioner, required, input               - Shift matrix
!  s:                  integer > 0, required, input                       - Size of shadow space
!  tolerance:          real(kind=rp), optional, input, default 1e-8       - Termination criterion
!  maximum_iterations: integer, optional, input, default min(2*n,1000)    - Maximum number of iterations
!  flag:               integer, optional, output                          - Convergence flags
!  relres:             real(kind=rp), array size n_sigma, optional, output- Relative residual norms
!  iterations:         integer, optional, output                          - Number of iterations to converge
!  inner_s             integer, optional, input, default 4                - Parameter s for inner iterations
!  inner_tolerance     real(kind-rp), optional, input, default = 1e-1     - Tolerance for inner iterations
!  inner_iterations    integer, optional, input, default 0                - Maximum number of inner iterations
!  omega:              complex(kind=cp), array, optional, input           - User defined values for omega
!  resvec:             complex(kind=cp), array size maxit+1, optional, output- Residual norm for every iteration
!
   IMPLICIT NONE

! Required input parameters:
   type(matrix), intent(in)               :: A       ! system matrix
   complex(kind=cp), intent(in)           :: b(:)    ! system rhs
   complex(kind=cp), intent(in)           :: sigma(:)! shifts
   type(preconditioner), intent(in)       :: M1      ! shift matrix
   integer, intent(in)                    :: s       ! s parameter

! Solution:
   complex(kind=cp)                       :: CMSQMRIDR(size(b),size(sigma) )

! Optional input parameters:
   real(kind=rp), optional, intent(in)    :: tolerance 
   integer, optional, intent(in)          :: maximum_iterations 

   integer, optional, intent(in)          :: inner_s 
   real(kind=rp), optional, intent(in)    :: inner_tolerance 
   integer, optional, intent(in)          :: inner_iterations 

! Optional output parameters:
   integer, optional, intent(out)         :: flag        
   real(kind=rp), optional, intent(out)   :: relres(size(sigma))
   integer, optional, intent(out)         :: iterations

! Optional input arrays:
   complex(kind=cp), optional, intent(in) :: omega(:) 

! Optional output arrays:
   real(kind=rp), optional, intent(out)   :: resvec(:,:) 

! Local arrays:
   complex(kind=cp),    allocatable       :: PT(:,:) 
   complex(kind=cp)                       :: x(size(b),size(sigma))
   complex(kind=cp)                       :: G(size(b),s), gn(size(b))
   complex(kind=cp)                       :: W(size(b),s+1,size(sigma)) 
   complex(kind=cp)                       :: wn(size(b))
   complex(kind=cp)                       :: v(size(b))
   complex(kind=cp)                       :: V_tilde(size(b),size(sigma)+1)
   complex(kind=cp)                       :: vtilde(size(b))   
   complex(kind=cp)                       :: M(s,s), f(s), mu(s)
   complex(kind=cp)                       :: alpha(s), beta(s), gamma(s)
   complex(kind=cp)                       :: c(s+2), h(s+2), r_sigma(s+3), tau
   complex(kind=cp)                       :: cs(s+2,size(sigma)), sn(s+2,size(sigma))
   complex(kind=cp)                       :: phi_n(size(sigma)), phi_n1(size(sigma))
   complex(kind=cp)                       :: eta(size(sigma))
   complex(kind=cp)                       :: cf(size(sigma)+1)
   complex(kind=cp)                       :: sigma_p(size(sigma))
   complex(kind=cp)                       :: sigma_0(size(sigma)+1)
   
   complex(kind=cp)                       :: om
   real(kind=rp)                          :: kappa
   complex(kind=cp)                       :: buffer(3)

   include 'msqmridr_body.f90'
   cmsqmridr = x
   
   end function CMSQMRIDR

   function RMSQMRIDR( A, b, sigma, M1, s, &           ! required
                       tolerance, maximum_iterations,& ! optional input
                       flag, relres, iterations, &     ! optional output
                       inner_s, inner_tolerance, inner_iterations, &
                       omega, resvec )                 ! optional arrays
!
! RMSQMRIDR: real version of Quasi-minimal residual multi-shift IDRS
! Parameters:
!  A:                  type matrix, required, input                       - System matrix
!  b:                  real(kind=rp), required, size is n, input          - Rhs of the system
!  sigma:              real(kind=rp), array size n_sigma, required, input - Shifts
!  M1:                 type preconditioner, required, input               - Shift matrix
!  s:                  integer > 0, required, input                       - Size of shadow space
!  tolerance:          real(kind=rp), optional, input, default 1e-8       - Termination criterion
!  maximum_iterations: integer, optional, input, default min(2*n,1000)    - Maximum number of iterations
!  flag:               integer, optional, output                          - Convergence flags
!  relres:             real(kind=rp), array size n_sigma, optional, output- Relative residual norms
!  iterations:         integer, optional, output                          - Number of iterations to converge
!  inner_s             integer, optional, input, default 4                - Parameter s for inner iterations
!  inner_tolerance     real(kind-rp), optional, input, default = 1e-1     - Tolerance for inner iterations
!  inner_iterations    integer, optional, input, default 0                - Maximum number of inner iterations
!  omega:              real(kind=rp), array, optional, input              - User defined values for omega
!  resvec:             real(kind=rp), array size maxit+1, optional, output- Residual norm for every iteration
!
   IMPLICIT NONE

! Required input parameters:
   type(matrix), intent(in)               :: A       ! system matrix
   real(kind=rp), intent(in)              :: b(:)    ! system rhs
   real(kind=rp), intent(in)              :: sigma(:)! shifts
   type(preconditioner), intent(in)       :: M1      ! shift matrix
   integer, intent(in)                    :: s       ! s parameter
   integer                                :: inner   ! inner iterations

! Solution:
   real(kind=rp)                          :: RMSQMRIDR(size(b),size(sigma) )

! Optional input parameters:
   real(kind=rp), optional, intent(in)    :: tolerance 
   integer, optional, intent(in)          :: maximum_iterations 

   integer, optional, intent(in)          :: inner_s 
   real(kind=rp), optional, intent(in)    :: inner_tolerance 
   integer, optional, intent(in)          :: inner_iterations 

! Optional output parameters:
   integer, optional, intent(out)         :: flag        
   real(kind=rp), optional, intent(out)   :: relres(size(sigma))
   integer, optional, intent(out)         :: iterations

! Optional input arrays:
   real(kind=rp), optional, intent(in)    :: omega(:) 

! Optional output arrays:
   real(kind=rp), optional, intent(out)   :: resvec(:,:) 
    
! Local arrays:
   real(kind=rp), allocatable             :: PT(:,:) 
   real(kind=rp)                          :: x(size(b),size(sigma))
   real(kind=rp)                          :: G(size(b),s), gn(size(b))
   real(kind=rp)                          :: W(size(b),s+1,size(sigma)) 
   real(kind=rp)                          :: wn(size(b))
   real(kind=rp)                          :: v(size(b))   
   real(kind=rp)                          :: vtilde(size(b))   
   real(kind=rp)                          :: V_tilde(size(b),size(sigma)+1)
   real(kind=rp)                          :: M(s,s), f(s), mu(s)
   real(kind=rp)                          :: alpha(s), beta(s), gamma(s)
   real(kind=rp)                          :: c(s+2), h(s+2), r_sigma(s+3), tau
   real(kind=rp)                          :: cs(s+2,size(sigma)), sn(s+2,size(sigma))
   real(kind=rp)                          :: phi_n(size(sigma)), phi_n1(size(sigma))
   real(kind=rp)                          :: eta(size(sigma))
   real(kind=rp)                          :: cf(size(sigma)+1)
   real(kind=rp)                          :: sigma_p(size(sigma))
   real(kind=rp)                          :: sigma_0(size(sigma)+1)

   real(kind=rp)                          :: om, kappa
   real(kind=rp)                          :: buffer(3)

   include 'msqmridr_body.f90'
   rmsqmridr = x
   
   end function RMSQMRIDR

   function RMSIDRS( A, b, sigma, M1, s, &                     ! required
                     tolerance, maximum_iterations, variant, & ! optional input
                     flag, relres, iterations,      &          ! optional output
                     omega, resvec, colfac )                   ! optional arrays
!
! RMSIDRS: real version of multi-shift IDRS
! Parameters:
!  A:                  type matrix, required, input                       - System matrix
!  b:                  real(kind=rp), required, size is n, input          - Rhs of the system
!  sigma:              real(kind=rp), array size n_sigma, required, input - Shifts
!  s:                  integer > 0, required, input                       - Size of shadow space
!  M1:                 type preconditioner, required, input               - Shift matrix
!  tolerance:          real(kind=rp), optional, input, default 1e-8       - Termination criterion
!  maximum_iterations: integer, optional, input, default min(2*n,1000)    - Maximum number of iterations
!  variant:            integer, optional, input, default 1                - Variant = 1: msidrs, 2: msbicgstab
!  flag:               integer, optional, output                          - Convergence flags
!  relres:             real(kind=rp), array size n_sigma, optional, output- Relative residual norms
!  iterations:         integer, optional, output                          - Number of iterations to converge
!  omega:              real(kind=rp), array, optional, input              - User defined values for omega
!  resvec:             real(kind=rp), array size maxit+1, optional, output- Residual norm for every iteration
!  colfac:             real(kind=rp), array size n_sigma, optional, output- For internal use by MSQMRIDR
!

   IMPLICIT NONE

! Required input parameters:
   type(matrix), intent(in)                :: A       ! system matrix
   real(kind=rp), intent(in)               :: b(:)    ! system rhs
   integer, intent(in)                     :: s       ! s parameter
   type(preconditioner), intent(in)        :: M1      ! shift matrix
   real(kind=rp), intent(in)               :: sigma(:)! shifts

! Solution:
   real(kind=rp)                           :: RMSIDRS(size(b),size(sigma) )

! Optional input parameters:
   real(kind=rp), optional, intent(in)     :: tolerance 
   integer, optional, intent(in)           :: maximum_iterations 
   integer, optional, intent(in)           :: variant

! Optional output parameters:
   integer, optional, intent(out)          :: flag        
   real(kind=rp), optional, intent(out)    :: relres(size(sigma))
   integer, optional, intent(out)          :: iterations

! Optional input arrays:
   real(kind=rp), optional, intent(in)     :: omega(:) 

! Optional output arrays:
   real(kind=rp), optional, intent(out)    :: resvec(:,:) 
   real(kind=rp), optional, intent(out)    :: colfac(:)
    
! Local arrays:
   real(kind=rp), allocatable              :: PT(:,:) 
   real(kind=rp)                           :: x(size(b),size(sigma))
   real(kind=rp)                           :: G(size(b),s), gn(size(b))
   real(kind=rp)                           :: W(size(b),s,size(sigma)) 
   real(kind=rp)                           :: v(size(b)), v0(size(b))   
   real(kind=rp)                           :: t(size(b))
   real(kind=rp)                           :: M(s,s), f(s)
   real(kind=rp)                           :: alpha(s), gamma(s)
   real(kind=rp)                           :: h(s+3), l_sigma(s+2,size(sigma)), u_sigma(s+2)
   real(kind=rp)                           :: phi(size(sigma))
   real(kind=rp)                           :: om_sigma(size(sigma))

   real(kind=rp)                           :: om, kappa
   real(kind=rp)                           :: buffer(3)
 
   include 'msidrs_body.f90'
   rmsidrs = x
   
   end function RMSIDRS

   function CMSIDRS( A, b, sigma, M1, s, &                     ! required
                     tolerance, maximum_iterations, variant, & ! optional input
                     flag, relres, iterations, &               ! optional output
                     omega, resvec, colfac )                   ! optional arrays
!
! CMSIDRS: complex version of multi-shift IDRS
! Parameters:
!  A:                  type matrix, required, input                       - System matrix
!  b:                  real(kind=rp), required, size is n, input          - Rhs of the system
!  sigma:              complex(kind=cp), array size n_sigma, required, input - Shifts
!  M1                  type preconditioner, required, input               - System matrix
!  s:                  integer > 0, required, input                       - Size of shadow space
!  tolerance:          real(kind=rp), optional, input, default 1e-8       - Termination criterion
!  maximum_iterations: integer, optional, input, default min(2*n,1000)    - Maximum number of iterations
!  variant:            integer, optional, input, default 1                - Variant = 1: msidrs, 2: msbicgstab
!  flag:               integer, optional, output                          - Convergence flags
!  relres:             real(kind=rp), array size n_sigma, optional, output- Relative residual norms
!  iterations:         integer, optional, output                          - Number of iterations to converge
!  omega:              complex(kind=cp), array, optional, input           - User defined values for omega
!  resvec:             real(kind=rp), array size maxit+1, optional, output- Residual norm for every iteration
!  colfac:             complex(kind=cp), array size n_sigma, optional, output- For internal use by MSQMRIDR
!

   IMPLICIT NONE

! Required input parameters:
   type(matrix), intent(in)                :: A       ! system matrix
   complex(kind=cp), intent(in)            :: b(:)    ! system rhs
   complex(kind=cp), intent(in)            :: sigma(:)! shifts
   type(preconditioner), intent(in)        :: M1      ! shift matrix
   integer, intent(in)                     :: s       ! s parameter

! Solution:
   complex(kind=cp)                        :: CMSIDRS(size(b),size(sigma) )

! Optional input parameters:
   real(kind=rp), optional, intent(in)     :: tolerance 
   integer, optional, intent(in)           :: maximum_iterations 
   integer, optional, intent(in)           :: variant

! Optional output parameters:
   integer, optional, intent(out)          :: flag        
   real(kind=rp), optional, intent(out)    :: relres(size(sigma))
   integer, optional, intent(out)          :: iterations

! Optional input arrays:
   complex(kind=cp), optional, intent(in)  :: omega(:) 

! Optional output arrays:
   real(kind=rp), optional, intent(out)    :: resvec(:,:) 
   complex(kind=cp), optional, intent(out) :: colfac(:)
    
! Local arrays:
   complex(kind=cp), allocatable           :: PT(:,:) 
   complex(kind=cp)                        :: x(size(b),size(sigma))
   complex(kind=cp)                        :: G(size(b),s), gn(size(b))
   complex(kind=cp)                        :: W(size(b),s,size(sigma)) 
   complex(kind=cp)                        :: v(size(b)), v0(size(b))   
   complex(kind=cp)                        :: t(size(b))
   complex(kind=cp)                        :: M(s,s), f(s)
   complex(kind=cp)                        :: alpha(s), gamma(s)
   complex(kind=cp)                        :: h(s+3), l_sigma(s+2,size(sigma)), u_sigma(s+2)
   complex(kind=cp)                        :: phi(size(sigma))
   complex(kind=cp)                        :: om_sigma(size(sigma))

   complex(kind=cp)                        :: om
   real(kind=rp)                           :: kappa
   complex(kind=cp)                        :: buffer(3)
 
   include 'msidrs_body.f90'
   cmsidrs = x
   
   end function CMSIDRS

!
! The functions below are for use inside the IDRS routines
!

   function CSOLVE(A,b) result(x)
! 
! Function to solve linear system
!
   IMPLICIT NONE
      complex(kind=cp), intent(in)     :: A(:,:), b(:)
      complex(kind=cp)                 :: x(size(b))
      integer                          :: i, j, ind(1), k, n
      complex(kind=cp)                 :: p, t(size(b)+1)
      complex(kind=cp)                 :: M(size(b),size(b)+1)    ! Augmented matrix

      n = size(b)
      M(1:n,1:n) = A
      M(:,n+1) = b

      do j = 1, n-1
! Partial pivotting:
         ind = maxloc( abs(M(j:n,j)) )
         k = ind(1)+j-1
         t = M(k,:)
         M(k,:) = M(j,:)
         M(j,:) = t
         do i = j + 1, n
            p = M(i,j)/M(j,j)
            M(i,j+1:n+1) = M(i,j+1:n+1) - M(j,j+1:n+1)*p
         end do
      end do

      do i = n, 1, -1
         p = M(i,n+1)
         do j = i+1, n
            p = p - M(i,j) * M(j,n+1)
         end do
         M(i,n+1) = p/M(i,i)
      end do
      x = M(:,n+1)
   end function CSOLVE

   function RSOLVE(A,b) result(x)
! 
! Function to solve linear system
!
   IMPLICIT NONE
      real(kind=rp), intent(in)        :: A(:,:), b(:)
      real(kind=rp)                    :: x(size(b))
      integer                          :: i, j, ind(1), k, n
      real(kind=rp)                    :: p, t(size(b)+1)
      real(kind=rp)                    :: M(size(b),size(b)+1)    ! Augmented matrix

      n = size(b)
      M(1:n,1:n) = A
      M(:,n+1) = b

      do j = 1, n-1
! Partial pivotting:
         ind = maxloc( abs(M(j:n,j)) )
         k = ind(1)+j-1
         t = M(k,:)
         M(k,:) = M(j,:)
         M(j,:) = t
         do i = j + 1, n
            p = M(i,j)/M(j,j)
            M(i,j+1:n+1) = M(i,j+1:n+1) - M(j,j+1:n+1)*p
         end do
      end do

      do i = n, 1, -1
         p = M(i,n+1)
         do j = i+1, n
            p = p - M(i,j) * M(j,n+1)
         end do
         M(i,n+1) = p/M(i,i)
      end do
      x = M(:,n+1)
   end function RSOLVE

   subroutine rgrot(a,b,c,s,r);
!
! RGROT: construct and apply Givens rotation
! RGROT returns c,s,r such that
!        | c       s || x |   | r |
!        |           ||   | = |   |
!        |-s       c || y | = | 0 |
!
! RGROT is based on the BLAS routine RROTG
!
! Declarations:
      real(kind=rp)            :: a,b
      real(kind=rp)            :: c,s,r
      real(kind=rp)            :: scale, rho
      real(kind=rp)            :: alpha

   if ( abs(a) < epsilon(scale) ) then
      c = 0.d0
      s = 1.d0
      r = b
   else
      scale = abs(a) + abs(b)
      rho = scale*sqrt(abs(a/scale)**2 +  abs(b/scale)**2)
      alpha = a/abs(a)
      c = abs(a)/rho
      s = alpha*b/rho
      r = alpha*rho
   end if
!
   end subroutine RGROT

   subroutine cgrot(a,b,c,s,r);
!
! CGROT: construct and apply Givens rotation
! CGROT returns c,s,r such that
!        | c        s || x |   | r |
!        |            ||   | = |   |
!        |-conjg(s) c || y | = | 0 |
!
! CGROT is based on the BLAS routine CROTG
!
! Declarations:
      complex(kind=cp)         :: a,b
      complex(kind=cp)         :: c,s,r
      real(kind=rp)            :: scale, rho
      complex(kind=cp)         :: alpha

   if ( abs(a) < epsilon(scale) ) then
      c = 0.
      s = 1.
      r = b
   else
      scale = abs(a) + abs(b)
      rho = scale*sqrt(abs(a/scale)**2 +  abs(b/scale)**2)
      alpha = a/abs(a)
      c = abs(a)/rho
      s = alpha*conjg(b)/rho
      r = alpha*rho
   end if
!
   end subroutine CGROT
!
   function RCOMP_OM( t, r, kappa ) result(om)
!
! Calcultate new omega using "Minimal residual" or
! "Maintaining convergence" technique depending on kappa
!
   implicit none
   real(kind=rp), intent(in)         :: t(:), r(:)
   real(kind=rp), intent(in)         :: kappa
   real(kind=rp)                     :: om
   real(kind=rp)                     :: rho
   real(kind=rp)                     :: buffer(3)

   if ( kappa == 0. ) then
! Minimal residual (same as in Bi-CGSTAB):
      buffer(1) = dot_product(t,r)
      buffer(2) = dot_product(t,t)
      call co_sum(buffer)
      om = buffer(1)/buffer(2)
   else

! 'Maintaining the convergence':
      buffer(1) = dot_product(r,r)
      buffer(2) = dot_product(t,t)
      buffer(3) = dot_product(t,r)
      call co_sum( buffer )
      rho = abs(buffer(3))/sqrt(abs(buffer(1)*buffer(2)))
      om=buffer(3)/buffer(2)
      if ( rho < kappa ) then
         om = om*kappa/rho
      end if
   end if

   end function RCOMP_OM

   function CCOMP_OM( t, r, kappa ) result(om)
!
! Calcultate new omega using "Minimal residual" or
! "Maintaining convergence" technique depending on kappa
!
   implicit none
   complex(kind=cp), intent(in)      :: t(:), r(:)
   real(kind=rp), intent(in)         :: kappa
   complex(kind=cp)                  :: om
   real(kind=rp)                     :: rho
   complex(kind=cp)                  :: buffer(3)

   if ( kappa == 0. ) then
! Minimal residual (same as in Bi-CGSTAB):
      buffer(1) = dot_product(t,r)
      buffer(2) = dot_product(t,t)
      call co_sum(buffer)
      om = buffer(1)/buffer(2)
   else

! 'Maintaining the convergence':
      buffer(1) = dot_product(r,r)
      buffer(2) = dot_product(t,t)
      buffer(3) = dot_product(t,r)
      call co_sum( buffer )
      rho = abs(buffer(3))/sqrt(abs(buffer(1)*buffer(2)))
      om=buffer(3)/buffer(2)
      if ( rho < kappa ) then
         om = om*kappa/rho 
      end if
   end if

   end function CCOMP_OM

   function CNORM( w )

! Compute norm of a complex vector

      complex(kind=cp), intent(in)              :: w(:)
      real(kind=rp)                             :: normw, CNORM

      normw = sum( abs(w)**2 )
      call co_sum( normw )
      CNORM = sqrt( normw )

   end function CNORM

   function RNORM( w )

! Compute norm of a real vector

      real(kind=rp), intent(in)                 :: w(:)
      real(kind=rp)                             :: normw, RNORM

      normw = sum( w**2 )
      call co_sum( normw )
      RNORM = sqrt( normw )

   end function RNORM
   
   function CINPROD(u, v )

! Compute inner product of a complex vectors

      complex(kind=cp), intent(in)              :: u(:), v(:)
      real(kind=rp)                             :: uv, CINPROD

      uv = dot_product( u, v )
      call co_sum( uv )
      CINPROD = uv

   end function CINPROD

   function RINPROD( u, v )

! Compute inner product of a real vectors

      real(kind=rp), intent(in)                 :: u(:), v(:)
      real(kind=rp)                             :: uv, RINPROD

      uv = dot_product( u, v )
      call co_sum( uv )
      RINPROD = uv

   end function RINPROD
   
   subroutine RMAKE_P( PT, r, method )

! Generate real random matrix P (transposed)

      real(kind=rp)  :: PT(:,:) 
      real(kind=rp)  :: r(:) 
      integer        :: method
      real(kind=rp)  :: alpha(size(PT,1)), rho
      integer        :: s, j, k

      if ( method == 1 ) then
! IDRS: P has orthogonal s random vectors
         s = size(PT,1)
         call RANDOM_NUMBER(PT)
! Orthogonalise:
         do j = 1,s
            do k = 1,j-1
               alpha(k) = dot_product( PT(j,:),PT(k,:) )
            end do
            call co_sum( alpha )
            do k = 1,j-1
               PT(j,:) = PT(j,:) - alpha(k)*PT(k,:)
            end do
            rho = norm( PT(j,:) )
            PT(j,:) = PT(j,:)/rho
         end do
      elseif ( method == 2 ) then
! BiCGSTAB: shadow vector equal to initial residual
         PT(1,:) = r
      else 
! P consists of s=n_procs sparse vectors that are distributed over
! the n_procs processors, one vector p per processor.
         call RANDOM_NUMBER(PT)
         rho = norm( PT(1,:) )
         PT(1,:) = PT(1,:)/rho
      end if

   end subroutine RMAKE_P

   subroutine CMAKE_P( PT, r, method )

! Generate complex random matrix P (transposed)

      complex(kind=cp)  :: PT(:,:) 
      complex(kind=cp)  :: r(:) 
      integer           :: method
      real(kind=rp)     :: Pr(size(PT,1),size(PT,2)), Pi(size(PT,1),size(PT,2)) 
      complex(kind=cp)  :: alpha(size(PT,1))
      real(kind=rp)     :: rho
      integer           :: s, j, k

      if ( method == 1 ) then
! IDRS: P has orthogonal random numbers
         s = size(PT,1)
         call RANDOM_NUMBER(Pr)
         call RANDOM_NUMBER(Pi)
         PT = cmplx(Pr,Pi)
         do j = 1,s
            do k = 1,j-1
               alpha(k) = dot_product( PT(k,:),PT(j,:) )
            end do
            call co_sum( alpha )
            do k = 1,j-1
               PT(j,:) = PT(j,:) - alpha(k)*PT(k,:)
            end do
            rho = norm( PT(j,:) )
            PT(j,:) = PT(j,:)/rho
         end do
      elseif ( method == 2 ) then
! BiCGSTAB: shadow vector equal to initial residual
         PT(1,:) = conjg(r)
      else
! P consists of s=n_procss sparse vectors that are distributed over
! the n_procs processors, one vector p per processor.
         call RANDOM_NUMBER(Pr)
         call RANDOM_NUMBER(Pi)
         PT = cmplx(Pr,Pi)
         rho = norm( PT(1,:) )
         PT(1,:) = PT(1,:)/rho
      end if

   end subroutine CMAKE_P

   function CP_DOT(PT, w, s, method)

! P inner product of complex matrices

      integer                                   :: s, method
      complex(kind=cp), intent(in)              :: PT(:,:) 
      complex(kind=cp), intent(in)              :: w(:)
      complex(kind=cp)                          :: v(s), CP_DOT(s)
      complex(kind=cp), save                    :: vt(1)[*]
      integer                                   :: i

      if ( method < 3 ) then
! For the standard methods, P consists of dense vectors that
! are distributed over the processors
         v = matmul( PT, w )
         call co_sum( v )
      else
! P consists of s=n_procs sparse vectors that are distributed over
! the n_procs processors, one vector p per processor.
         vt = matmul( PT, w )
         sync all
         do i = 1,s
            v(i) = vt(1)[i]
         end do
      end if
      
      CP_DOT = v

   end function CP_DOT

   function RP_DOT( PT, w, s, method )

! P inner product of real matrices

      integer                              :: s, method
      real(kind=rp), intent(in)            :: PT(:,:) 
      real(kind=rp), intent(in)            :: w(:)
      real(kind=rp)                        :: v(s), RP_DOT(s)
      real(kind=rp), save                  :: vt(1)[*]
      integer                              :: i

      if ( method < 3 ) then
! For the standard methods, P consists of dense vectors that
! are distributed over the processors
         v = matmul( PT, w )
         call co_sum( v )
      else
! P consists of s=n_procs sparse vectors that are distributed over
! the n_procs processors, one vector p per processor.
         vt = matmul( PT, w )
         sync all
         do i = 1,s
            v(i) = vt(1)[i]
         end do
      end if
         
      RP_DOT = v

   end function RP_DOT 

end module idrs_module

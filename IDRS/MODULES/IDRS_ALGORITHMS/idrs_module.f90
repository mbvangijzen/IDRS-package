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
!   The IDRS module requires a user defined TYPE MATRIX.
!   This type should be defined in a module MATRIX_MODULE. 
!   The MATRIX_MODULE must furthermore contain functions 
!   XR = RMATVEC( A, YR ) and XC = CMATVEC( A, YC ) for
!   performing the matrix-vector product, and
!   XR = RPRECON( M, YR ) and XC = CPRECON( M, YR ) for performing the preconditioning operation.
!   The vectors XR and YR must be of type REAL(KIND=RP), and XC and YC must be of type COMPLEX(KIND=CP).
!   A and M must be of type matrix.
!   The operator * must be overloaded to perform the operation A*y, and / to perform y/M.
!
!   The precision of the complex and real variables are defined by the 
!   parameters RP (real precision)and CP (complex precision). These parameters 
!   should be set in the module precision_module.
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2025 Martin van Gijzen
!

module idrs_module

   use precision_module
   use dense_la_module
   use matrix_module

   implicit none

   private
   public :: IDRS, QMRIDR, MSQMRIDR, MSIDRS, IDREIG, IDR_HESSENBERG, COMPUTE_EIGVEC, COMPUTE_EIGVAL

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

   interface IDREIG
      module procedure CIDREIG, RIDREIG
   end interface

   interface IDR_HESSENBERG
      module procedure CIDR_HESSENBERG, RIDR_HESSENBERG
   end interface

   interface COMPUTE_EIGVEC
      module procedure CCOMPUTE_EIGVEC, RCOMPUTE_EIGVEC
   end interface

   interface COMPUTE_EIGVAL
      module procedure CCOMPUTE_EIGVAL, RCOMPUTE_EIGVAL
   end interface

   interface P_DOT
      module procedure CP_DOT, RP_DOT
   end interface

   interface MAKE_P
      module procedure RMAKE_P, CMAKE_P
   end interface

   interface COMP_OM
      module procedure CCOMP_OM, RCOMP_OM
   end interface

contains

   function CIDRS( A, b, s, &                                     ! required
                   M1, tolerance, maximum_iterations, variant, &  ! optional input
                   flag, relres, iterations, &                    ! optional output
                   x0, U0, omega, resvec, H ) result(x)           ! optional arrays
!
! CIDRS: complex version of IDRS
! Parameters:
!  A:                  type matrix, required, input                       - System matrix
!  b:                  complex(kind=cp), required, size is n, input       - Rhs of the system
!  s:                  integer > 0, required, input                       - Size of shadow space
!  M1:                 type matrix, optional, input                       - Preconditioner
!  tolerance:          real(kind=rp), optional, input, default 1e-6       - Termination criterion
!  maximum_iterations: integer, optional, input, default min(2*n,1000)    - Maximum number of iterations
!  variant:            integer, optional, input, default 1                - Variant = 1: idrs, 2: bicgstab
!  flag:               integer, optional, output                          - Convergence flags
!  relres:             real(kind=rp), optional, output                    - Relative residual norm
!  iterations:         integer, optional, output                          - Number of iterations to converge
!  x0:                 complex(kind=cp), array size n, optional, input    - Initial guess
!  U0:                 complex(kind=cp), array size n*s, optional, input  - Initial search space
!  omega:              complex(kind=cp), array, optional, input           - User defined values for omega
!  resvec:             real(kind=rp), array size maxit+1, optional, output- Residual norm for every iteration
!  H:                  complex(kind=cp), array size (nritz+1)*nritz,      - Upper Hessenberger
!                                        optional, output                 - can be used to compute ritzvalues
!
   IMPLICIT NONE

! Required input parameters:
   type(matrix), intent(in)               :: A         
   complex(kind=cp), intent(in)           :: b(:)     
   integer, intent(in)                    :: s      

! Solution:
   complex(kind=cp)                       :: x(size(b))

! Optional input parameters:
   type(matrix), optional, intent(in)     :: M1      
   real(kind=rp), optional, intent(in)    :: tolerance 
   integer, optional, intent(in)          :: maximum_iterations     
   integer, optional, intent(in)          :: variant     

! Optional output parameters:
   integer, optional, intent(out)         :: flag        
   real(kind=rp), optional, intent(out)   :: relres     
   integer, optional, intent(out)         :: iterations     

! Optional input arrays:
   complex(kind=cp), optional, intent(in) :: x0(:)  
   complex(kind=cp), optional, allocatable, &
                               intent(in) :: U0(:,:) 
   complex(kind=cp), optional, intent(in) :: omega(:) 

! Optional output arrays:
   real(kind=rp), optional, intent(out)   :: resvec(:) 
   complex(kind=cp), optional, intent(out):: H(:,:)   
    
! Local arrays:
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

   include 'idrs.inc'
   
   end function CIDRS

   function RIDRS( A, b, s, &                                    ! required
                   M1, tolerance, maximum_iterations, variant, & ! optional input
                   flag, relres, iterations, &                   ! optional output
                   x0, U0, omega, resvec, H ) result(x)          ! optional arrays
!
! RIDRS: real version of IDRS
! Parameters:
!  A:                  type matrix, required, input                       - System matrix
!  b:                  real(kind=rp), required, size is n, input          - Rhs of the system
!  s:                  integer > 0, required, input                       - Size of shadow space
!  M1:                 type matrix, optional, input                       - Preconditioner
!  tolerance:          real(kind=rp), optional, input, default 1e-6       - Termination criterion
!  maximum_iterations: integer, optional, input, default min(2*n,1000)    - Maximum number of iterations
!  variant:            integer, optional, input, default 1                - Variant = 1: idrs, 2: bicgstab
!  flag:               integer, optional, output                          - Convergence flags
!  relres:             real(kind=rp), optional, output                    - Relative residual norm
!  iterations:         integer, optional, output                          - Number of iterations to converge
!  x0:                 real(kind=rp), array size n, optional, input       - Initial guess
!  U0:                 real(kind=rp), array size n*s, optional, input     - Initial search space
!  omega:              real(kind=rp), array, optional, input              - User defined values for omega
!  resvec:             real(kind=rp), array size maxit+1, optional, output- Residual norm for every iteration
!  H:                  real(kind=rp), array size (nritz+1)*nritz,         - Upper Hessenberger
!                                        optional, output                 - can be used to compute ritzvalues
!

   IMPLICIT NONE

! Required input and output parameters:
   type(matrix), intent(in)               :: A             
   real(kind=rp), intent(in)              :: b(:)         
   integer, intent(in)                    :: s          

! Solution
   real(kind=rp)                          :: x(size(b))

! Optional input parameters:
   type(matrix), optional, intent(in)     :: M1          
   real(kind=rp), optional, intent(in)    :: tolerance
   integer, optional, intent(in)          :: maximum_iterations
   integer, optional, intent(in)          :: variant

! Optional output parameters:
   integer, optional, intent(out)         :: flag
   real(kind=rp), optional, intent(out)   :: relres
   integer, optional, intent(out)         :: iterations 

! Optional input arrays:
   real(kind=rp), optional, intent(in)    :: x0(:)
   real(kind=rp), optional, allocatable, &
                            intent(in)    :: U0(:,:)
   real(kind=rp), optional, intent(in)    :: omega(:)

! Optional output arrays
   real(kind=rp), optional, intent(out)   :: resvec(:)
   real(kind=cp), optional, intent(out)   :: H(:,:)  
    
! Local arrays:
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
    
   include 'idrs.inc'
     
   end function RIDRS

   function CQMRIDR( A, b, s, &                                    ! required
                     M1, tolerance, maximum_iterations, variant, & ! optional input
                     flag, relres, iterations, &                   ! optional output
                     inner_s, inner_tolerance, inner_iterations, & ! optional input                    
                     x0, omega, resvec ) result(x)                 ! optional arrays
!
! CQMRIDR: complex version of Quasi-minimal residual IDRS
! Parameters:
!  A:                  type matrix, required, input                       - System matrix
!  b:                  complex(kind=cp), required, size is n, input       - Rhs of the system
!  s:                  integer > 0, required, input                       - Size of shadow space
!  M1:                 type matrix, optional, input                       - Preconditioner
!  tolerance:          real(kind=rp), optional, input, default 1e-6       - Termination criterion
!  maximum_iterations: integer, optional, input, default min(2*n,1000)    - Maximum number of iterations
!  flag:               integer, optional, output                          - Convergence flags
!  relres:             real(kind=rp), optional, output                    - Relative residual norm
!  iterations:         integer, optional, output                          - Number of iterations to converge
!  variant:            integer, optional, input, default 1                - Variant = 1: qmridr, 2: qmrcgstab
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
   integer, intent(in)                    :: s       ! s parameter

! Solution:
   complex(kind=cp)                       :: x(size(b) )

! Optional input parameters:
   type(matrix), optional, intent(in)     :: M1      ! preconditioner
   real(kind=rp), optional, intent(in)    :: tolerance 
   integer, optional, intent(in)          :: maximum_iterations 
   integer, optional, intent(in)          :: variant     

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

   include 'qmridr.inc'
   
   end function CQMRIDR

   function RQMRIDR( A, b, s, &                                          ! required
                     M1, tolerance, maximum_iterations, variant, &       ! optional input
                     flag, relres, iterations, &                         ! optional output
                     inner_s, inner_tolerance, inner_iterations, &       ! optional input
                     x0, omega, resvec ) result(x)                       ! optional arrays
!
! RQMRIDR: real version of Quasi-minimal resdual IDRS
! Parameters:
!  A:                  type matrix, required, input                       - System matrix
!  b:                  real(kind=rp), required, size is n, input          - Rhs of the system
!  s:                  integer > 0, required, input                       - Size of shadow space
!  M1:                 type matrix, optional, input                       - Preconditioner
!  tolerance:          real(kind=rp), optional, input, default 1e-6       - Termination criterion
!  maximum_iterations: integer, optional, input, default min(2*n,1000)    - Maximum number of iterations
!  variant:            integer, optional, input, default 1                - Variant = 1: qmridr, 2: qmrcgstab
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
   integer, intent(in)                    :: s       ! s parameter

! Solution:
   real(kind=rp)                          :: x(size(b) )

! Optional input parameters:
   type(matrix), optional, intent(in)     :: M1      ! preconditioner
   real(kind=rp), optional, intent(in)    :: tolerance 
   integer, optional, intent(in)          :: maximum_iterations 
   integer, optional, intent(in)          :: variant     

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

   include 'qmridr.inc'
   
   end function RQMRIDR

   function CMSQMRIDR( A, b, sigma, s, &                                  ! required
                       M1, tolerance, maximum_iterations, variant, &      ! optional input
                       flag, relres, iterations, &                        ! optional output
                       inner_s, inner_tolerance, inner_iterations, &      ! optional input
                       omega, resvec ) result(x)                          ! optional arrays
!
! CMSQMRIDR: complex version of Quasi-minimal residual multi-shift IDRS
! Parameters:
!  A:                  type matrix, required, input                       - System matrix
!  b:                  complex(kind=cp), required, size is n, input       - Rhs of the system
!  sigma:              complex(kind=cp), array size n_sigma, required, input - Shifts
!  s:                  integer > 0, required, input                       - Size of shadow space
!  M1:                 type matrix, optional, input                       - Mass matrix
!  tolerance:          real(kind=rp), optional, input, default 1e-6       - Termination criterion
!  maximum_iterations: integer, optional, input, default min(2*n,1000)    - Maximum number of iterations
!  variant:            integer, optional, input, default 1                - Variant = 1: msqmrdrs, 2: msqmrcgstab
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
   integer, intent(in)                    :: s       ! s parameter

! Solution:
   complex(kind=cp)                       :: x(size(b),size(sigma) )

! Optional input parameters:
   type(matrix), optional, intent(in)     :: M1      ! mass matrix
   real(kind=rp), optional, intent(in)    :: tolerance 
   integer, optional, intent(in)          :: maximum_iterations 
   integer, optional, intent(in)          :: variant     

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

   include 'msqmridr.inc'
   
   end function CMSQMRIDR

   function RMSQMRIDR( A, b, sigma, s, &                             ! required
                       M1, tolerance, maximum_iterations, variant, & ! optional input
                       flag, relres, iterations, &                   ! optional output
                       inner_s, inner_tolerance, inner_iterations, & ! optional input
                       omega, resvec ) result(x)                     ! optional arrays
!
! RMSQMRIDR: real version of Quasi-minimal residual multi-shift IDRS
! Parameters:
!  A:                  type matrix, required, input                       - System matrix
!  b:                  real(kind=rp), required, size is n, input          - Rhs of the system
!  sigma:              real(kind=rp), array size n_sigma, required, input - Shifts
!  s:                  integer > 0, required, input                       - Size of shadow space
!  M1:                 type matrix, optional, input                       - Mass matrix
!  tolerance:          real(kind=rp), optional, input, default 1e-6       - Termination criterion
!  maximum_iterations: integer, optional, input, default min(2*n,1000)    - Maximum number of iterations
!  variant:            integer, optional, input, default 1                - Variant = 1: msqmrdrs, 2: msqmrcgstab
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
   integer, intent(in)                    :: s       ! s parameter

! Solution:
   real(kind=rp)                          :: x(size(b),size(sigma) )

! Optional input parameters:
   type(matrix), optional, intent(in)     :: M1      ! Mass matrix
   real(kind=rp), optional, intent(in)    :: tolerance 
   integer, optional, intent(in)          :: maximum_iterations 
   integer, optional, intent(in)          :: variant     

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
   integer                                :: inner   ! inner iterations

   include 'msqmridr.inc'
   
   end function RMSQMRIDR

   function RMSIDRS( A, b, sigma, s, &                             ! required
                     M1, tolerance, maximum_iterations, variant, & ! optional input
                     flag, relres, iterations, &                   ! optional output
                     omega, resvec, colfac ) result(x)             ! optional arrays
!
! RMSIDRS: real version of multi-shift IDRS
! Parameters:
!  A:                  type matrix, required, input                       - System matrix
!  b:                  real(kind=rp), required, size is n, input          - Rhs of the system
!  sigma:              real(kind=rp), array size n_sigma, required, input - Shifts
!  s:                  integer > 0, required, input                       - Size of shadow space
!  M1:                 type matrix, optional, input                       - Shift matrix
!  tolerance:          real(kind=rp), optional, input, default 1e-6       - Termination criterion
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
   real(kind=rp), intent(in)               :: sigma(:)! shifts

! Solution:
   real(kind=rp)                           :: x(size(b),size(sigma) )

! Optional input parameters:
   type(matrix), optional, intent(in)      :: M1      ! Mass matrix
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
 
   include 'msidrs.inc'
   
   end function RMSIDRS

   function CMSIDRS( A, b, sigma, s, &                             ! required
                     M1, tolerance, maximum_iterations, variant, & ! optional input
                     flag, relres, iterations, &                   ! optional output
                     omega, resvec, colfac ) result(x)             ! optional arrays
!
! CMSIDRS: complex version of multi-shift IDRS
! Parameters:
!  A:                  type matrix, required, input                       - System matrix
!  b:                  real(kind=rp), required, size is n, input          - Rhs of the system
!  sigma:              complex(kind=cp), array size n_sigma, required, input - Shifts
!  s:                  integer > 0, required, input                       - Size of shadow space
!  M1                  type matrix, optional, input                       - Mass matrix
!  tolerance:          real(kind=rp), optional, input, default 1e-6       - Termination criterion
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
   integer, intent(in)                     :: s       ! s parameter

! Solution:
   complex(kind=cp)                        :: x(size(b),size(sigma) )

! Optional input parameters:
   type(matrix), optional, intent(in)      :: M1      ! shift matrix
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
 
   include 'msidrs.inc'
   
   end function CMSIDRS

   function CIDREIG( A, b, s, neig, &                                        ! required input
                     M1, tolerance, maximum_iterations, variant, &           ! optional input
                     eigvec, flag, resnrm, iterations, omega, resvec ) &     ! optional output
                     result(eigval)                                          ! output
!
! CIDREIG: IDRS for eigenvalues, complex variant
! Parameters:
!  A:                  type matrix, required, input                         - System matrix
!  b:                  complex(kind=cp), array size N                       - starting vector
!  s:                  integer > 0, required, input                         - Size of shadow space
!  neig:               integer                                              - Number of eigenvalues
!  M1:                 type matrix, optional, input                         - Mass matrix
!  eigvec:             complex(kind=cp), array size N*neig, optional, output- Eigenvectors
!  flag:               integer, optional, output                            - Error flag
!  resnrm:             real(kind=rp), optional, output                      - Relative residual norms
!  iterations:         integer, optional, output                            - Number of iterations
!  omega:              complex(kind=cp), array, optional, input             - User defined shifts
!  resvec:             real(kind=rp), array size maximum_iterations*neig,   - Residual norms
!                      optional output                                        per iteration
!
! Output:    
!  eigval:             complex(kind=cp), array size(neig)
!

   IMPLICIT NONE

! Required input and output parameters:
   type(matrix), intent(in)                :: A             
   complex(kind=cp), intent(in)            :: b(:)
   integer, intent(in)                     :: s          
   integer, intent(in)                     :: neig

! Eigenvalues
   complex(kind=cp)                        :: eigval(neig)

! Optional input parameters:
   type(matrix), optional, intent(in)      :: M1          
   real(kind=rp), optional, intent(in)     :: tolerance 
   integer, intent(in), optional           :: maximum_iterations
   integer, optional, intent(in)           :: variant 
   complex(kind=cp), optional, intent(in)  :: omega(:) 

! Optional output parameters:
   complex(kind=cp), optional, intent(out) :: eigvec(size(b),neig)
   integer, optional, intent(out)          :: flag        
   real(kind=rp), optional, intent(out)    :: resnrm(neig) 
   integer, optional, intent(out)          :: iterations
   real(kind=rp), optional, intent(out)    :: resvec(:,:) 

! Local arrays:
   complex(kind=cp), allocatable           :: H(:,:)
   complex(kind=cp), allocatable           :: PT(:,:)
   complex(kind=cp)                        :: G(size(b),s)
   complex(kind=cp)                        :: r(size(b))
   complex(kind=cp)                        :: v(size(b))
   complex(kind=cp)                        :: t(size(b))
   complex(kind=cp)                        :: M(s,s), f(s), mu(s)
   complex(kind=cp)                        :: alpha(s), gamma(s)

   complex(kind=cp)                        :: shift
   complex(kind=cp)                        :: buffer(3)

   include "idreig.inc"
   
   end function CIDREIG

   function RIDREIG( A, b, s, neig, &                                        ! required input
                     M1, tolerance, maximum_iterations, variant, &           ! optional input
                     eigvec, flag, resnrm, iterations, omega, resvec ) &     ! optional output
                     result(eigval)                                          ! output
!
! RIDREIG: IDRS for eigenvalues, real variant
!
!    Note: THIS ALGORITHM ASSUMES THAT THE EIGENVALUES AND EIGENVALUES AER REAL!!!
!
! Parameters:
!  A:                  type matrix, required, input                         - System matrix
!  b:                  real(kind=rp), array size N                          - starting vector
!  s:                  integer > 0, required, input                         - Size of shadow space
!  neig:               integer, input                                       - Number of eigenvalues
!  M1:                 type matrix, optional, input                         - Mass matrix
!  eigvec:             real(kind=rp), array size N*neig, optional, output   - Eigenvectors
!  flag:               integer, optional, output                            - Error flag
!  resnrm:             real(kind=rp), optional, output                      - Relative residual norms
!  iterations:         integer, optional, output                            - Number of iterations
!  omega:              complex(kind=cp), array, optional, input             - User defined shifts
!  resvec:             real(kind=rp), array size maximum_iterations*neig,   - Residual norms
!                      optional output                                        per iteration
!
! Output:    
!  eigval:             real(kind=rp), array size(neig)
!

   IMPLICIT NONE

! Required input and output parameters:
   type(matrix), intent(in)                :: A             
   real(kind=rp), intent(in)               :: b(:)
   integer, intent(in)                     :: s 
   integer, intent(in)                     :: neig 

! Eigenvalues
   real(kind=rp)                           :: eigval(neig)

! Optional input parameters:
   type(matrix), optional, intent(in)      :: M1          
   real(kind=rp), optional, intent(in)     :: tolerance 
   integer, intent(in), optional           :: maximum_iterations
   integer, optional, intent(in)           :: variant
   real(kind=rp), optional, intent(in)     :: omega(:) 

! Optional output parameters:
   real(kind=rp), optional, intent(out)    :: eigvec(size(b),neig)
   integer, optional, intent(out)          :: flag        
   real(kind=rp), optional, intent(out)    :: resnrm(neig) 
   integer, optional, intent(out)          :: iterations
   real(kind=rp), optional, intent(out)    :: resvec(:,:) 

! Local arrays:
   real(kind=rp), allocatable              :: H(:,:)
   real(kind=rp), allocatable              :: PT(:,:)
   real(kind=rp)                           :: G(size(b),s)
   real(kind=rp)                           :: r(size(b))
   real(kind=rp)                           :: v(size(b))
   real(kind=rp)                           :: t(size(b))
   real(kind=rp)                           :: M(s,s), f(s), mu(s)
   real(kind=rp)                           :: alpha(s), gamma(s)

   real(kind=cp)                           :: shift
   real(kind=cp)                           :: buffer(3)

   include "idreig.inc"
   
   end function RIDREIG

   function CIDR_HESSENBERG( A, b, s, nritz, &                      ! required input
                      M1, variant, &                                ! optional input
                      flag, &                                       ! optional output
                      omega )           &                           ! optional array
                      result(H)                                     ! output
!
! IDR_HESSENBERG: compute IDRS Hessenberg matrix
! Parameters:
!  A:                  type matrix, required, input                         - System matrix
!  b:                  complex(kind=cp), array size(N)                      - starting vector
!  s:                  integer > 0, required, input                         - Size of shadow space
!  N:                  integer > 0, required, input                         - Local equations
!  nritz:              integer > 0, required, input                         - Number of ritzvalues
!  M1:                 type matrix, optional, input                         - Mass matrix
!  variant:            integer, optional, input                             - IDREIG or BiLanczos
!  flag:               integer, optional, output                            - Error flag
!  omega:              complex(kind=cp), array size(n_omega)                - omega parameters
!
! Output:    
!  H:                  complex(kind=cp), array size((nritz+1)xnritz)        - Hessenberg matrix
!

   IMPLICIT NONE

! Required input and output parameters:
   type(matrix), intent(in)                :: A             
   complex(kind=cp), intent(in)            :: b(:)
   integer, intent(in)                     :: s          
   integer, intent(in)                     :: nritz          

! Optional input parameters:
   type(matrix), optional, intent(in)      :: M1          

   integer, optional, intent(in)           :: variant

! Optional output parameters:
   integer, optional, intent(out)          :: flag        

! Optional input arrays:
   complex(kind=cp), optional, intent(in)  :: omega(:) 

! Output arrays:
   complex(kind=cp)                        :: H(nritz+1,nritz)

! Number of equations:                     :: N

! Local arrays:
   complex(kind=cp), allocatable           :: PT(:,:) 
   complex(kind=cp)                        :: G(size(b),s)
   complex(kind=cp)                        :: r(size(b))
   complex(kind=cp)                        :: v(size(b))         
   complex(kind=cp)                        :: t(size(b))        
   complex(kind=cp)                        :: M(s,s), f(s), mu(s)
   complex(kind=cp)                        :: alpha(s), gamma(s)
    
   complex(kind=cp)                        :: shift
   complex(kind=cp)                        :: buffer(3)

   include "idr_hessenberg.inc"

   end function CIDR_HESSENBERG

   function RIDR_HESSENBERG( A, b, s, nritz, &                      ! required input
                      M1, variant, &                                ! optional input
                      flag, &                                       ! optional output
                      omega ) &                                     ! optional array
                      result(H)                                     ! output
!
! IDR_HESSENBERG: compute IDRS Hessenberg matrix
! Parameters:
!  A:                  type matrix, required, input                         - System matrix
!  b:                  real(kind=rp), array size(N)                         - starting vector
!  s:                  integer > 0, required, input                         - Size of shadow space
!  nritz:              integer > 0, required, input                         - Number of ritzvalues
!  M1:                 type matrix, optional, input                         - Mass matrix
!  variant:            integer, optional, input                             - IDREIG or BiLanczos
!  flag:               integer, optional, output                            - Error flag
!  omega:              real(kind=rp), array size(n_omega)                   - omega parameters
!
! Output:    
!  ritzval:            complex(kind=cp), array size(nritz)
!

   IMPLICIT NONE

! Required input and output parameters:
   type(matrix), intent(in)                :: A             
   real(kind=rp), intent(in)               :: b(:)
   integer, intent(in)                     :: s          
   integer, intent(in)                     :: nritz          

! Optional input parameters:
   type(matrix), optional, intent(in)      :: M1          

   integer, optional, intent(in)           :: variant

! Optional output parameters:
   integer, optional, intent(out)          :: flag        

! Optional input arrays:
   real(kind=rp), optional, intent(in)     :: omega(:) 

! Output arrays:
   real(kind=rp)                           :: H(nritz+1,nritz)
    
! Number of equations:                     :: N

! Local arrays:
   real(kind=rp), allocatable              :: PT(:,:) 
   real(kind=rp)                           :: G(size(b),s)
   real(kind=rp)                           :: r(size(b))
   real(kind=rp)                           :: v(size(b))         
   real(kind=rp)                           :: t(size(b))        
   real(kind=rp)                           :: M(s,s), f(s), mu(s)
   real(kind=rp)                           :: alpha(s), gamma(s)
    
   real(kind=rp)                           :: shift
   real(kind=rp)                           :: buffer(3)

   include "idr_hessenberg.inc"

   end function RIDR_HESSENBERG

!
! The functions below are for use inside the IDRS routines
!

   function rcompute_eigval( H, neig, err_out ) result(eigval)
!
! Compute eigenvalues 
!    In:    H, real, size (nritz+1)*(nritz)       Hessenberg matrix from IDREIG
!           neig, integer                         Number of wanted eigenvalues
!    Out:   eigval, real, size neig               Wanted eigenvalues
!           err_out, real, optional, size(neig)   Relative shift in eigenvalue
!
      real(kind=rp), intent(in)            :: H(:,:)
      integer, intent(in)                  :: neig
      real(kind=rp)                        :: eigval(neig)
      real(kind=rp), optional              :: err_out(neig)

      complex(kind=cp)                     :: ritzval(size(H,2)), ritzval_old(size(H,2))
      real(kind=rp)                        :: err(size(H,2)), tol
      logical                              :: mask(size(H,2))
      integer                              :: nritz, i, j, k

! Number of Ritzvalues:
      nritz = size(H,2)

! Tolerance to determine if eigenvalues are the same
      tol = sqrt(epsilon(1._rp))

! Compute ritzvalues and check
      ritzval_old(1:nritz-1) = HQR( H(1:nritz-1,1:nritz-1) )
      ritzval                = HQR( H(1:nritz,1:nritz) )

! Determine the error (relative shift) in the ritzvalues
      do i = 1,nritz
         err(i) = minval( abs(ritzval_old(1:nritz-1)-ritzval(i))/abs(ritzval(i)) )
      end do

! Select all eigenvalues:
      mask = .true.

! Find the neig eigenvalues with the smallest errors:
      do i = 1, neig
         j = minloc( err, dim = 1, mask = mask )
         if ( j == 0 ) then
            exit
         else
            eigval(i) = real(ritzval(j),kind=rp)
            if ( present(err_out) ) err_out(i) = err(j)
         end if
! Unselect multiple eigenvalues:
         do k = 1, nritz
            if ( abs(ritzval(k) - ritzval(j))  < tol*abs(ritzval(j)) ) mask(k) = .false.
         end do
      end do

   end function rcompute_eigval

   function ccompute_eigval( H, neig, err_out ) result(eigval)
!
! Compute eigenvalues 
!    In:    H, complex, size (nritz+1)*(nritz)    Hessenberg matrix from IDREIG
!           neig, integer                         Number of wanted eigenvalues
!    Out:   eigval, complex, size neig            Wanted eigenvalues
!           err_out, real, optional, size(neig)   Relative shift in eigenvalue
!
      complex(kind=cp), intent(in)         :: H(:,:)
      integer, intent(in)                  :: neig
      complex(kind=cp)                     :: eigval(neig)
      real(kind=rp), optional              :: err_out(neig)

      complex(kind=cp)                     :: ritzval(size(H,2)), ritzval_old(size(H,2))
      real(kind=rp)                        :: err(size(H,2)), tol
      logical                              :: mask(size(H,2))
      integer                              :: nritz, ieig, i, j, k

! Number of Ritzvalues:
      nritz = size(H,2)

! Tolerance to determine if eigenvalues are the same
      tol = sqrt(epsilon(1._rp))

! Compute ritzvalues and check
      ritzval_old(1:nritz-1) = HQR( H(1:nritz-1,1:nritz-1) )
      ritzval                = HQR( H(1:nritz,1:nritz) )

! Determine the error (relative shift) in the ritzvalues
      do i = 1,nritz
         err(i) = minval( abs(ritzval_old(1:nritz-1)-ritzval(i))/abs(ritzval(i)) )
      end do

! Select all eigenvalues:
      mask = .true.

! Find the neig eigenvalues with the smallest errors:
      ieig = 0
      do i = 1, neig
         j = minloc( err, dim = 1, mask = mask )
         if ( j == 0 ) then
            exit
         else
            eigval(i) = ritzval(j)
            if ( present(err_out) ) err_out(i) = err(j)
         end if
! Unselect multiple eigenvalues:
         do k = 1, nritz
            if ( abs(ritzval(k) - ritzval(j))  < tol*abs(ritzval(j)) ) mask(k) = .false.
         end do
      end do
   end function ccompute_eigval

   function rcompute_eigvec( A, s, H, eigval, b, M, tolerance, resnrm, resvec ) result(eigvec)
!
! Compute eigenvectors
!    In:   A, type matrix                    System matrix
!          s, integer                        IDREIG parameter
!          H, real, size (nritz+1)*(nritz)   Hessenberg matrix from IDREIG
!          eigval, real, size neig           Real eigenvalues of Hessenberg matrix
!          b, real, size neq                 Starting vector Krylov space
!          M, type matrix, optional          Mass matrix
!          tolerance, real, optional         Tolerance              
!          resnrm, real, size neig, optional Residual norms (estimate)
!          resvec, real, size nritz*neig, 
!                  optional                  Convergence history
!    Out:  eigvec, real, size neq*neig       Real eigenvectors 
!
! Optional:
!    In:   M, type matrix                   Mass matrix
!
      type(matrix), intent(in)             :: A
      integer                              :: s
      real(kind=rp), intent(in)            :: H(:,:)
      real(kind=rp), intent(in)            :: eigval(:)
      real(kind=rp), intent(in)            :: b(:)
      type(matrix), intent(in), optional   :: M
      real(kind=rp), intent(in), optional  :: tolerance
      real(kind=rp), intent(out), optional :: resnrm(size(eigval))
      real(kind=rp), intent(out), optional :: resvec(:,:)
      real(kind=rp)                        :: eigvec(size(b),size(eigval))

! Local variables:
      logical                              :: converged(size(eigval))
      real(kind=rp)                        :: eigvec_H(size(H,2),size(eigval))
      real(kind=rp)                        :: resnrm_loc(size(eigval))
      real(kind=rp)                        :: W(size(b),s+1)
      real(kind=rp)                        :: v(size(b))
      real(kind=rp)                        :: normu, normv, tol
      integer                              :: ind(size(H,2))
      integer                              :: nritz, neig, neq, i, j, k, low
      logical                              :: mass_matrix

      nritz = size(H,2)
      neq   = size(b)
      neig  = size(eigval)
      mass_matrix = present(M)
      converged = .false.
      if ( present(tolerance) ) then
         tol = tolerance
      else
         tol = sqrt(epsilon(1._rp))
      end if

! Compute eigenvectors of H:
      do i = 1, neig
         eigvec_H(1:nritz,i) = SI_power( H(1:nritz,1:nritz), eigval(i) )
      end do

      eigvec = 0._rp
      W = 0._rp
      v = b/norm(b)
      ind = 0       ! Indicates which which basis vector corresponds to which column of W
      k = 0
      do i = 1,nritz
         k = k + 1
         if ( k > s+1 ) k = 1
         ind(i) = k
         W(:,k) = v
         do j = 1, neig
            if ( .not. converged(j) ) eigvec(:,j) = eigvec(:,j) + eigvec_H(i,j)*v
         end do
         low = max(1,i-s)
         if ( mass_matrix ) then
            v = (A*(v/M) - matmul(W(:,ind(low:i)),H(low:i,i))) /H(i+1,i)
         else
            v = (A*v - matmul(W(:,ind(low:i)),H(low:i,i))) /H(i+1,i)
         end if

!
! Has the eigenpair converged?
         normv = norm(v)
         do j = 1, neig
            if ( .not. converged(j) ) then
               normu = norm(eigvec(:,j))
               resnrm_loc(j) = normv*abs(H(i+1,i)*eigvec_H(i,j))/normu
               converged(j) = ( normv*abs(H(i+1,i)*eigvec_H(i,j))/normu < tol*abs(eigval(j)) )
            end if
         end do
!
! Save convergence history:
         if ( present(resvec) ) resvec(i,:) = resnrm_loc
      end do
!
! Normalise
      do j = 1,neig
         normu = norm(eigvec(:,j))
         eigvec(:,j) = eigvec(:,j)/normu
      end do
!
! Output the residual norms
      if ( present(resnrm) ) resnrm = resnrm_loc

   end function rcompute_eigvec

   function ccompute_eigvec( A, s, H, eigval, b, M, tolerance, resnrm, resvec ) result(eigvec)
!
! Compute initial search space
!    In:   A, type matrix                      System matrix
!          s, integer                          IDREIG parameter
!          H, complex, size (nritz+1)*(nritz)  Hessenberg matrix from IDREIG
!          eigval, complex, size nritz*neig    Complex eigenvalues of Hessenberg matrix
!          b, complex, size neq                Starting vector Krylov space
!          M, type matrix, optional            Mass matrix
!          tolerance, real, optional           Tolerance              
!          resnrm, real, size neig, optional   Residual norms (estimate)
!          resvec, real, size nritz*neig, 
!                  optional                  Convergence history
!    Out:  eigvec, complex, size neq*neig      Complex eigenvectors 
!
! Optional:
!    In:   M, type matrix                      Mass matrix
!
      implicit none

      type(matrix), intent(in)             :: A
      integer                              :: s
      complex(kind=cp), intent(in)         :: H(:,:)
      complex(kind=cp), intent(in)         :: eigval(:)
      complex(kind=cp), intent(in)         :: b(:)
      type(matrix), intent(in), optional   :: M
      real(kind=rp), intent(in), optional  :: tolerance
      real(kind=rp), intent(out), optional :: resnrm(size(eigval))
      real(kind=rp), intent(out), optional :: resvec(:,:)
      complex(kind=cp)                     :: eigvec(size(b),size(eigval))

! Local variables:
      logical                              :: converged(size(eigval))
      complex(kind=cp)                     :: eigvec_H(size(H,2),size(eigval))
      complex(kind=cp)                     :: W(size(b),s+1)
      real(kind=rp)                        :: resnrm_loc(size(eigval))
      complex(kind=cp)                     :: v(size(b))
      real(kind=rp)                        :: normu, normv, tol
      integer                              :: ind(size(H,2))
      integer                              :: nritz, neig, neq, i, j, k, low
      logical                              :: mass_matrix

      nritz = size(H,2)
      neq   = size(b)
      neig  = size(eigval)
      mass_matrix = present(M)
      converged = .false.
      if ( present(tolerance) ) then
         tol = tolerance
      else
         tol = sqrt(epsilon(1._cp))
      end if

! Compute eigenvectors of H:
      do i = 1, neig
         eigvec_H(1:nritz,i) = SI_power( H(1:nritz,1:nritz), eigval(i) )
      end do

      eigvec = 0._cp
      W = 0._cp
      v = b/norm(b)
      ind = 0       ! Indicates which which basis vector corresponds to which column of W
      k = 0

      do i = 1,nritz
         k = k + 1
         if ( k > s+1 ) k = 1
         ind(i) = k
         W(:,k) = v
         do j = 1, neig
            if ( .not. converged(j) ) eigvec(:,j) = eigvec(:,j) + eigvec_H(i,j)*v
         end do
         low = max(1,i-s)
         if ( mass_matrix ) then
            v = (A*(v/M) - matmul(W(:,ind(low:i)),H(low:i,i))) /H(i+1,i)
         else
            v = (A*v - matmul(W(:,ind(low:i)),H(low:i,i))) /H(i+1,i)
         end if

!
! Has the eigenpair converged?
         normv = norm(v)
         do j = 1, neig
            if ( .not. converged(j) ) then
               normu = norm(eigvec(:,j))
               resnrm_loc(j) = normv*abs(H(i+1,i)*eigvec_H(i,j))/normu
               converged(j) = ( normv*abs(H(i+1,i)*eigvec_H(i,j))/normu < tol*abs(eigval(j)) )
            end if
         end do
!
! Save convergence history:
         if ( present(resvec) ) resvec(i,:) = resnrm_loc
      end do
!
! Normalise the eigenvectors
      do j = 1,neig
         normu = norm(eigvec(:,j))
         eigvec(:,j) = eigvec(:,j)/normu
      end do

!
! Output the residual norms
      if ( present(resnrm) ) resnrm = resnrm_loc

   end function ccompute_eigvec

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

   subroutine RMAKE_P( PT, r, method )

! Generate real random matrix P (transposed)

      real(kind=rp)        :: PT(:,:) 
      real(kind=rp)        :: r(:) 
      integer              :: method
      real(kind=rp)        :: alpha(size(PT,1)), rho
      integer              :: s, j, k
      integer, allocatable :: seed(:)
      integer              :: n

! Set the seed of the random-number generator to a fixed value for reproducability
      call random_seed(size=n)
      allocate(seed(n))
      seed = this_image()
      call random_seed(put=seed)

      if ( method == 1 ) then
! IDRS: P has s random orthogonalised vectors
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

      complex(kind=cp)     :: PT(:,:) 
      complex(kind=cp)     :: r(:) 
      integer              :: method
      real(kind=rp)        :: Pr(size(PT,1),size(PT,2)), Pi(size(PT,1),size(PT,2)) 
      complex(kind=cp)     :: alpha(size(PT,1))
      real(kind=rp)        :: rho
      integer              :: s, j, k
      integer, allocatable :: seed(:) 
      integer              :: n

! Set the seed of the random-number generator to a fixed value for reproducability
      call random_seed(size=n)
      allocate(seed(n))
      seed = this_image()
      call random_seed(put=seed)

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
! P consists of s=n_procs sparse vectors that are distributed over
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
      integer                                   :: me

      if ( method < 3 ) then
! For the standard methods, P consists of dense vectors that
! are distributed over the processors
         v = matmul( PT, w )
         call co_sum( v )
      else
! P consists of s=n_procs sparse vectors that are distributed over
! the n_procs processors, one vector p per processor.
         v = 0._cp
         me = this_image()
         v(me) = dot_product( PT(1,:), w )
         call co_sum( v )
      end if
      
      CP_DOT = v

   end function CP_DOT

   function RP_DOT( PT, w, s, method )

! P inner product of real matrices

      integer                              :: s, method
      real(kind=rp), intent(in)            :: PT(:,:) 
      real(kind=rp), intent(in)            :: w(:)
      real(kind=rp)                        :: v(s), RP_DOT(s)
      integer                              :: me

      if ( method < 3 ) then
! For the standard methods, P consists of dense vectors that
! are distributed over the processors
         v = matmul( PT, w )
         call co_sum( v )
      else
! P consists of s=n_procs sparse vectors that are distributed over
! the n_procs processors, one vector p per processor.
         v = 0._rp
         me = this_image()
         v(me) = dot_product( PT(1,:), w )
         call co_sum( v )
      end if
         
      RP_DOT = v

   end function RP_DOT 

end module idrs_module

!
! The recycle-module contains subroutines to extract approximations 
! to the Schur vectors for IDRS. 
! These vectors can be used to create a recycling subspace, 
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2024 Martin van Gijzen
! 

module recycle_module

   use precision_module
   use dense_la_module
   use matrix_module
   use idrs_module

   implicit none

   interface compute_U0
      module procedure ccompute_U0, rcompute_U0
   end interface

   interface IDRS_SCHUR
      module procedure CIDRS_SCHUR, RIDRS_SCHUR
   end interface

contains

   function rcompute_U0( A, H, eigvec, b, M, resnrm ) result(U0)
!
! Compute initial search space
!    In:   A, type matrix                   System matrix
!          H, real, size (nritz+1)*(nritz)  Hessenberg matrix from IDRS
!          eigvec, real, size nritz*s       eigenvectors of Hessenberg matrix
!          b, real, size neq                Starting vector Krylov space
!    Out:  U0, real, size neq*s             Initial search space
! Optional:
!    In:   M, type matrix, optional         Mass matrix
!    Out:  resnrm, real, optional, size s   residual norms for ritz pairs
!
      type(matrix), intent(in)             :: A
      real(kind=rp), intent(in)            :: H(:,:)
      real(kind=rp), intent(in)            :: eigvec(:,:)
      real(kind=rp), intent(in)            :: b(:)
      type(matrix), intent(in), optional   :: M
      real(kind=rp), optional, intent(out) :: resnrm(:)
      real(kind=rp)                        :: U0(size(b),size(eigvec,2))
      real(kind=rp)                        :: W(size(b),size(eigvec,2)+1)
      real(kind=rp)                        :: v(size(b))
      real(kind=rp)                        :: nrm
      integer                              :: ind(size(H,2))
      integer                              :: nritz, s, neq, i, j, k, low

      nritz = size(H,2)
      neq   = size(b)
      s     = size(eigvec,2)

!
! Compute the basis vectors for the Krylov subspace and the Ritzvector:
      U0 = 0.
      W = 0.
      v = b
      ind = 0       ! Indicates which column of W corresponds to which basis vector
      k = 0
      do i = 1,nritz
         k = k + 1
         if ( k > s+1 ) k = 1
         ind(i) = k
         W(:,k) = v
         do j = 1, s
            U0(:,j) = U0(:,j)+ eigvec(i,j)*v
         end do
         low = max(1,i-s)
         if ( present(M) ) then
            v = (A*(v/M) - matmul(W(:,ind(low:i)),H(low:i,i))) /H(i+1,i)
         else
            v = (A*v - matmul(W(:,ind(low:i)),H(low:i,i))) /H(i+1,i)
         end if
      end do
!
! Normalise the Schur vectors, compute the residual norms
      if ( present(resnrm) ) then
         do j = 1,s
            nrm = norm(U0(:,j))
            U0(:,j) = U0(:,j)/nrm
            resnrm(j) = norm(v)*abs(H(nritz+1,nritz))*abs(eigvec(nritz,j))/nrm
         end do
      end if
   end function rcompute_U0

   function ccompute_U0( A, H, eigvec, b, M, resnrm ) result(U0)
!
! Compute inital search space U0
!    In:   A, type matrix                      System matrix
!          H, complex, size (nritz+1)*nritz    Hessenberg matrix from IDRS
!          eigvec, complex, size nritz*s       eigenvectors of Hessenberg matrix
!          b, complex, size neq                Starting vector Krylov space
!    Out:  U0, complex, size neq*s             Search space
! Optional:
!    In:   M, type matrix, optional            Mass matrix
!    Out:  resnrm, real, optional, size s      residual norms for ritz pairs
!
      type(matrix), intent(in)                :: A
      complex(kind=cp), intent(in)            :: H(:,:)
      complex(kind=cp), intent(in)            :: eigvec(:,:)
      complex(kind=cp), intent(in)            :: b(:)
      type(matrix), optional, intent(in)      :: M
      real(kind=rp), optional, intent(out)    :: resnrm(:)
      complex(kind=cp)                        :: U0(size(b),size(eigvec,2))
      complex(kind=cp)                        :: W(size(b),size(eigvec,2)+1)
      complex(kind=cp)                        :: v(size(b))
      real(kind=rp)                           :: nrm
      integer                                 :: ind(size(H,2))
      integer                                 :: nritz, s, neq, i, j, k, low

      nritz = size(H,2)
      neq   = size(b)
      s     = size(eigvec,2)

!
! Compute the basis vectors for the Krylov subspace and the Ritzvector:
      U0 = 0.
      W = 0.
      v = b
      ind = 0       ! Indicates which which basis vector basis vector corresponds to which column of W
      k = 0
      do i = 1,nritz
         k = k + 1
         if ( k > s+1 ) k = 1
         ind(i) = k
         W(:,k) = v
         do j = 1, s
            U0(:,j) = U0(:,j)+ eigvec(i,j)*v
         end do
         low = max(1,i-s)
         if ( present(M) ) then
            v = (A*(v/M) - matmul(W(:,ind(low:i)),H(low:i,i))) /H(i+1,i)
         else
            v = (A*v - matmul(W(:,ind(low:i)),H(low:i,i))) /H(i+1,i)
         end if
      end do
!
! Normalise the Schur vectors, compute the residual norms
      if ( present(resnrm) ) then
         do j = 1,s
            nrm = norm(U0(:,j))
            U0(:,j) = U0(:,j)/nrm
            resnrm(j) = norm(v)*abs(H(nritz+1,nritz))*abs(eigvec(nritz,j))/nrm
         end do
      end if

   end function ccompute_U0

   function CIDRS_SCHUR( A, b, s, &                                     ! required
                         M, tolerance, maximum_iterations, variant, &   ! optional input
                         flag, relres, iterations, &                    ! optional output
                         x0, U0, omega, resvec, H ) result(x)           ! optional arrays, output
!
! CIDRS_RECYCLE: complex version of IDRS_RECYCLE
! Parameters:
!  A:                  type matrix, required, input                       - System matrix
!  b:                  complex(kind=cp), required, size is n, input       - Rhs of the system
!  s:                  integer > 0, required, input                       - Size of shadow space
!  M:                  type matrix, optional, input                       - Preconditioner
!  tolerance:          real(kind=rp), optional, input, default 1e-6       - Termination criterion
!  maximum_iterations: integer, optional, input, default min(2*n,1000)    - Maximum number of iterations
!  variant:            integer, optional, input, default 1                - Variant = 1: idrs, 2: bicgstab
!  flag:               integer, optional, output                          - Convergence flags
!  relres:             real(kind=rp), optional, output                    - Relative residual norm
!  iterations:         integer, optional, output                          - Number of iterations to converge
!  x0:                 complex(kind=cp), array size n, optional, input    - Initial guess
!  U0:                 complex(kind=cp), array size n*s, optional, output - Approximate Schur vectors
!  omega:              complex(kind=cp), array, optional, input           - User omegas
!  resvec:             real(kind=rp), array size maxit+1, optional, output- Residual norm for every iteration
!  H:                  complex(kind=rp), array size (nritz+1)*nritz       - Upper Hessenberger
!
! This routines calls the standard IDRS method to solve a linear system. 
! IDRS_SCHUR, however, also returns approximations to the Schur vectors.
! 

   IMPLICIT NONE

! Required input parameters:
   type(matrix), intent(in)                             :: A
   complex(kind=cp), intent(in)                         :: b(:)
   integer, intent(in)                                  :: s

! Solution:
   complex(kind=cp)                                     :: x(size(b))

! Optional input parameters:
   type(matrix), optional, intent(in)                   :: M
   real(kind=rp), optional, intent(in)                  :: tolerance
   integer, optional, intent(in)                        :: maximum_iterations
   integer, optional, intent(in)                        :: variant

! Optional output parameters:
   integer, optional, intent(out)                       :: flag
   real(kind=rp), optional, intent(out)                 :: relres
   integer, optional, intent(out)                       :: iterations

! Optional input arrays:
   complex(kind=cp), optional, intent(in)               :: x0(:)

! Optional output arrays:
   complex(kind=cp), optional, allocatable, intent(inout) :: U0(:,:)
   complex(kind=cp), optional, intent(in)               :: omega(:)
   real(kind=rp), optional                              :: resvec(:)
   complex(kind=cp), optional, intent(out)              :: H(:,:)

   integer                                              :: nritz
   real(kind=rp)                                        :: relnrm, tol
   real(kind=rp), allocatable                           :: resval(:)
   real(kind=rp), allocatable                           :: resnrm(:)
   complex(kind=cp), allocatable                        :: ritzval(:)
   complex(kind=cp), allocatable                        :: x_start(:)
   complex(kind=cp), allocatable                        :: U(:,:), R(:,:)
   complex(kind=cp), allocatable                        :: Hf(:,:), vec(:,:), V(:,:)
   integer                                              :: n, maxit, converged, i, iter, conv, method, dimU0

   n = size(b)

   if ( present(tolerance) ) then
      if ( tolerance < 0 ) stop "Illegal value for parameter tolerance"
      tol = tolerance
   else
      tol = sqrt(epsilon(1._cp))
   endif

   maxit=min(2*n,1000)
   if ( present(maximum_iterations) ) maxit = maximum_iterations

   method = 1 ! idrs
   if ( present(variant) ) method = variant

   if ( present(x0) ) then
      x_start = x0
   else
      allocate(x_start(n))
      x_start = 0.
   end if

   allocate(resnrm(maxit+1), Hf(maxit+1,maxit))

   if ( present(M) ) then
      x = idrs( A, b, s, M, tol, maxit, method, conv, relnrm, iter, x_start, omega=omega, resvec = resnrm, H = Hf )
   else
      x = idrs( A, b, s, tolerance = tol, maximum_iterations = maxit, variant = method, flag = conv, &
                omega=omega, relres = relnrm, iterations = iter, x0 = x_start, resvec = resnrm, H = Hf )
   end if
   nritz = (iter/(s+1))*s + mod(iter,(s+1))

   if ( present(flag) ) flag = conv
   if ( present(relres) ) relres = relnrm
   if ( present(iterations) ) iterations = iter
   if ( present(resvec) ) resvec = resnrm

   if ( present(U0) .and. allocated(U0) ) then
      dimU0 = size(U0,2)
      if ( nritz < dimU0 ) stop "Recycle space is too large, reduce dimension of U0."

! Compute the Schur vectors of (upper part of) H:
      allocate( vec(nritz,dimU0), resval(dimU0), V(nritz,nritz) )
      ritzval = QR_eig( Hf(1:nritz,1:nritz), V )
      vec =  V(:,1:dimU0)

      allocate( U(size(U0,1),dimU0), R(dimU0,dimU0) )
      if ( present(M) ) then
         U = compute_U0( A, Hf(1:nritz+1,1:nritz), vec, b, M, resval )
         do i = 1, dimU0
            U(:,i) = U(:,i)/M
         end do
      else
         U = compute_U0( A, Hf(1:nritz+1,1:nritz), vec, b, resnrm = resval )
      end if
      U0 = orthogonalize( U ) 
   end if

   if ( present(H) ) then
      if ( size(H) > size(Hf) ) stop 'Dimensions of Hessenberg matrix too large'
      H = Hf(1:size(H,1),1:size(H,2))
   end if

   end function CIDRS_SCHUR

   function RIDRS_SCHUR( A, b, s, &                                     ! required
                         M, tolerance, maximum_iterations, variant, &   ! optional input
                         flag, relres, iterations, &                    ! optional output
                         x0, U0, omega, resvec, H ) result(x)           ! optional arrays, output
!
! RIDRS: real version of IDRS_SCHUR
! Parameters:
!  A:                  type matrix, required, input                       - System matrix
!  b:                  real(kind=rp), required, size is n, input          - Rhs of the system
!  s:                  integer > 0, required, input                       - Size of shadow space
!  M:                  type matrix, optional, input                       - Preconditioner
!  tolerance:          real(kind=rp), optional, input, default 1e-6       - Termination criterion
!  maximum_iterations: integer, optional, input, default min(2*n,1000)    - Maximum number of iterations
!  variant:            integer, optional, input, default 1                - Variant = 1: idrs, 2: bicgstab
!  flag:               integer, optional, output                          - Convergence flags
!  relres:             real(kind=rp), optional, output                    - Relative residual norm
!  iterations:         integer, optional, output                          - Number of iterations to converge
!  x0:                 real(kind=rp), array size n, optional, input       - Initial guess
!  U0:                 real(kind=rp), array size n*s, optional, output    - Initial search space
!  omega:              real(kind=rp), array, optional, input              - Omega's based on Ritzvalues
!  resvec:             real(kind=rp), array size maxit+1, optional, output- Residual norm for every iteration
!  H:                  real(kind=rp), array size (nritz+1)*nritz,         - Upper Hessenberger
!
! This routines calls the standard IDRS method to solve a linear system. 
! IDRS_SCHUR, however, also returns approximations to the Schur vectors.
! 

   IMPLICIT NONE

! Required input parameters:
   type(matrix), intent(in)                             :: A
   real(kind=rp), intent(in)                            :: b(:)
   integer, intent(in)                                  :: s

! Solution:
   real(kind=rp)                                        :: x(size(b))

! Optional input parameters:
   type(matrix), optional, intent(in)                   :: M
   real(kind=rp), optional, intent(in)                  :: tolerance
   integer, optional, intent(in)                        :: maximum_iterations
   integer, optional, intent(in)                        :: variant

! Optional output parameters:
   integer, optional, intent(out)                       :: flag
   real(kind=rp), optional, intent(out)                 :: relres
   integer, optional, intent(out)                       :: iterations

! Optional input arrays:
   real(kind=rp), optional, intent(in)                  :: x0(:)

! Optional output arrays:
   real(kind=rp), optional, allocatable, intent(inout)  :: U0(:,:)
   real(kind=rp), optional, intent(in)                  :: omega(:)
   real(kind=rp), optional                              :: resvec(:)
   real(kind=rp), optional, allocatable, intent(out)    :: H(:,:)

   integer                                              :: nritz, maxit
   real(kind=rp)                                        :: relnrm, tol
   real(kind=rp), allocatable                           :: resval(:)
   complex(kind=cp), allocatable                        :: ritzval(:), V(:,:)
   real(kind=rp), allocatable                           :: resnrm(:), x_start(:)
   real(kind=rp), allocatable                           :: U(:,:), R(:,:)
   real(kind=rp), allocatable                           :: Hf(:,:), vec(:,:)
   integer                                              :: n, converged, i, k, iter, conv, method, dimU0

   n = size(b)

   if ( present(tolerance) ) then
      if ( tolerance < 0 ) stop "Illegal value for parameter tolerance"
      tol = tolerance
   else
      tol = sqrt(epsilon(1._rp))
   endif

   maxit=min(2*n,1000)
   if ( present(maximum_iterations) ) maxit = maximum_iterations

   method = 1 ! idrs
   if ( present(variant) ) method = variant

   if ( present(x0) ) then
      x_start = x0
   else
      allocate(x_start(n))
      x_start = 0.
   end if

   allocate(resnrm(maxit+1), Hf(maxit+1,maxit))

   if ( present(M) ) then
      x = idrs( A, b, s, M, tol, maxit, method, conv, relnrm, iter, x_start, resvec = resnrm, H = Hf )
   else
      x = idrs( A, b, s, tolerance = tol, maximum_iterations = maxit, variant = method, flag = conv, &
                relres = relnrm, iterations = iter, x0 = x_start, resvec = resnrm, H = Hf )
   end if
   nritz = (iter/(s+1))*s + mod(iter,(s+1))

   if ( present(flag) ) flag = conv
   if ( present(relres) ) relres = relnrm
   if ( present(iterations) ) iterations = iter
   if ( present(resvec) ) resvec = resnrm

   if ( present(U0) .and. allocated(U0) ) then
      dimU0 = size(U0,2)
      if ( nritz < dimU0 ) stop "Recycle space is too large, reduce dimension of U0."

! Compute the Schur vectors of (upper part of) H:
      allocate( vec(nritz,dimU0), resval(dimU0), V(nritz,nritz) )
! Note: here the shift-and-invert power method is used to
! compute a real invariant subspace. This avoids extracting a real
! subspace from the complex V computed by QR_eig:
      vec = SI_block_power( Hf(1:nritz,1:nritz), dimU0 )

      allocate( U(size(U0,1),dimU0), R(dimU0,dimU0) )
      if ( present(M) ) then
         U = compute_U0( A, Hf(1:nritz+1,1:nritz), vec, b, M, resval )
         do i = 1, dimU0
            U(:,i) = U(:,i)/M
         end do
      else
         U = compute_U0( A, Hf(1:nritz+1,1:nritz), vec, b, resnrm = resval )
      end if
      U0 = orthogonalize( U ) 
   end if

   if ( present(H) ) then
      if ( size(H) > size(Hf) ) stop 'Dimensions of Hessenberg matrix too large'
      H = Hf(1:size(H,1),1:size(H,2))
   end if

   end function RIDRS_SCHUR

end module recycle_module

!
! The dense_la_module contains a number of linear algebra routines for small dense matrices.
! The module is provided to make the IDRS-package self-contained.
! 
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2025 Martin van Gijzen
! 

module dense_la_module

   use precision_module

   implicit none

   interface QR_factorization
      module procedure CQR_factorization, RQR_factorization
   end interface

   interface QR_solve
      module procedure CQR_solve, RQR_solve
   end interface

   interface QR_eig
      module procedure CQR_eig, RQR_eig
   end interface

   interface SI_power
      module procedure CSI_power, RSI_power
   end interface

   interface SOLVE
      module procedure CSOLVE, RSOLVE
   end interface

   interface GROT
      module procedure CGROT, RGROT 
   end interface 

contains

   subroutine CQR_factorization( V, Q, R, info )

!
! Compute QR decomposition of complex matrix V
!    In:   V, complex, size m*n
!    Out:  Q, complex, orthogonal, size m*n
!    Out:  R, complex, upper triangular, size n*n
!
      complex(kind=cp), intent(in)      :: V(:,:)
      complex(kind=cp), intent(out)     :: Q(size(V,1),size(V,2)), R(size(V,2),size(V,2))
      integer, optional, intent(out)    :: info
      logical                           :: out_info
      complex(kind=cp)                  :: t
      integer                           :: i, j, k, m, n
      integer, parameter                :: refine = 2

      m = size(V,1)
      n = size(V,2)
      out_info = present(info) 
      if ( out_info ) info = 0

      Q = V
      R = 0.
      do j = 1, n
         do k = 1, refine
            do i = 1, j-1
               t      = dot_product( Q(:,i), Q(:,j) )
               Q(:,j) = Q(:,j) - t*Q(:,i)
               R(i,j) = R(i,j) + t
            end do
         end do
         R(j,j) = R(j,j) + norm2( abs(Q(:,j)) )
         if ( abs(R(j,j)) > tiny(1._rp) ) then 
            Q(:,j) = Q(:,j)/R(j,j)
         else
            if ( out_info ) info = j
            return
         end if
      end do

   end subroutine CQR_factorization

   subroutine RQR_factorization( V, Q, R, info )
 
!
! Compute QR decomposition of real matrix V
!    In:   V, real, size m*n
!    Out:  Q, real, orthogonal, size m*n
!    Out:  R, real, upper triangular, size n*n
!
      real(kind=rp), intent(in)         :: V(:,:)
      real(kind=rp), intent(out)        :: Q(size(V,1),size(V,2)), R(size(V,2),size(V,2))
      integer, optional, intent(out)    :: info
      logical                           :: out_info
      real(kind=rp)                     :: t
      integer                           :: i, j, k, m, n
      integer, parameter                :: refine = 2

      m = size(V,1)
      n = size(V,2)
      out_info = present(info)
      if ( out_info ) info = 0

      Q = V
      R = 0.
      do j = 1, n
         do k = 1, refine
            do i = 1, j-1
               t      = dot_product( Q(:,i), Q(:,j) )
               Q(:,j) = Q(:,j) - t*Q(:,i)
               R(i,j) = R(i,j) + t
            end do
         end do
         R(j,j) = norm2( Q(:,j) )
         if ( R(j,j) > tiny(1._rp) ) then 
            Q(:,j) = Q(:,j)/R(j,j)
         else
            if ( out_info ) info = j
            return
         end if
      end do

   end subroutine RQR_factorization

   function CQR_solve( Q, R, B ) result(X)
!
! Compute solution of complex system AX = B using QR factorization of A
!
! In:
!    Q:    complex, size m*n
!    R:    complex, size n*n
!    B:    complex, size m*s
! Out:  
!    X:    complex, size n*s
!
      complex(kind=cp), intent(in)         :: B(:,:)
      complex(kind=cp), intent(in)         :: Q(:,:), R(:,:)
      complex(kind=cp)                     :: X(size(Q,2),size(B,2))
      integer                              :: i, j, m, n, s

      m = size(Q,1)
      n = size(Q,2)
      s = size(B,2)

! Compute Q^Tb
      X = matmul(transpose(conjg(Q)),B)

! Backsubstitution:
      do j = n, 1, -1
         X(j,:) = X(j,:)/R(j,j)
         do i = 1, j-1
            X(i,:) = X(i,:) - R(i,j)*X(j,:)
         end do
      end do

   end function CQR_solve

   function RQR_solve( Q, R, B ) result(X)
!
! Compute solution of complex system AX = B using QR factorization of A
!
! In:
!    Q:    real, size m*n
!    R:    real, size n*n
!    B:    real, size m*s
! Out:  
!    X:    real, size n*s
!
      real(kind=rp), intent(in)            :: B(:,:)
      real(kind=rp), intent(in)            :: Q(:,:), R(:,:)
      real(kind=rp)                        :: X(size(Q,2),size(B,2))
      integer                              :: i, j, m, n, s

      m = size(Q,1)
      n = size(Q,2)
      s = size(B,2)

! Compute Q^Tb
      X = matmul(transpose(Q),B)

! Backsubstitution:
      do j = n, 1, -1
         X(j,:) = X(j,:)/R(j,j)
         do i = 1, j-1
            X(i,:) = X(i,:) - R(i,j)*X(j,:)
         end do
      end do

   end function RQR_solve

   function CQR_eig( H, U, T, converged ) result(eig)
!
! Compute eigenvalues of H using QR algorithm
!    In:   H, complex, size n*n            Matrix
!    Out:  eig, complex, size n            Eigenvalues
!          U, complex, size n*n, optional  Schur vectors
!          T, complex, size n*n, optional  Upper triangular matrix in Schur decomposition
!          converged, integer, optional    number of converged eigenvalues
!
      complex(kind=cp), intent(in)            :: H(:,:)
      complex(kind=cp), intent(out), optional :: U(size(H,2),size(H,2))
      complex(kind=cp), intent(out), optional :: T(size(H,2),size(H,2))
      integer, optional                       :: converged
      complex(kind=cp)                        :: eig(size(H,2))
      complex(kind=cp)                        :: shift
      complex(kind=cp)                        :: V(size(H,2),size(H,2))
      complex(kind=cp), allocatable           :: Q(:,:), R(:,:), S(:,:)
      integer                                 :: i, iter, m, n
      integer                                 :: maxit
      real(kind=rp), parameter                :: tol = epsilon(1._rp)
      integer                                 :: info
      logical                                 :: compute_U, compute_T

      n = size(H,2)
      V = H(1:n,1:n)
      iter = 0
      maxit = 10*n
      m = n 

      compute_U = present(U)
      compute_T = present(T)
!
      if ( compute_U ) then 
         allocate(S(n,n))
         U = 0._cp
         do i = 1, n
            U(i,i) = 1._cp
         end do
      end if

      shift = V(n,n)
!
      do while ( m > 1 .and. iter < maxit ) 

!
! Apply shift:
         do i = 1, m
            V(i,i) = V(i,i) - shift
         end do

!
! QR-iteration step:
         allocate(Q(m,m), R(m,m))
         call QR_factorization( V(1:m,1:m), Q, R, info )
         if ( info > 0 ) exit
         V(1:m,1:m) = matmul( R, Q )
            
!
! Correct for shift:
         do i = 1, m
            V(i,i) = V(i,i) + shift
         end do

!
! For Schur vectors:
         if ( compute_U ) then
            S = 0._cp
            do i = 1, n
               S(i,i) = 1._cp
            end do
            S(1:m,1:m) = Q
            U = matmul( U, S )
         end if

!
! For Upper triangular matrix in the Schur decomposition:
         if ( compute_T ) then
            if ( m < n  ) V(1:m,m+1:n) = matmul( transpose(conjg(Q)), V(1:m,m+1:n) )
         end if

!
! Convergence of eigenvalue?
         if ( abs( V(m,m-1) ) < tol*abs(V(m,m) ) ) m = m - 1

!
! Prepare for new iteration:
         shift = V(m,m)
         iter = iter + 1
         deallocate(Q,R)

      end do

!
! Store the eigenvalues:
      do i = 1, n
         eig(i) = V(i,i)
      end do
      if ( present(converged) ) converged = n-m+1

      if ( compute_T ) T = V

   end function CQR_eig

   function RQR_eig( H, U, T, converged ) result(eig)
!
! Compute eigenvalues of H using QR algorithm
!    In:   H, real, size n*n               Matrix
!    Out:  eig, complex, size n            Eigenvalues
!          U, complex, size n*n, optional  Schur vectors
!          T, complex, size n*n, optional  Upper triangular matrix in Schur decomposition
!          converged, integer, optional    number of converged eigenvalues
!
      real(kind=rp), intent(in)               :: H(:,:)
      complex(kind=cp), intent(out), optional :: U(size(H,2),size(H,2))
      complex(kind=cp), intent(out), optional :: T(size(H,2),size(H,2))
      integer, optional                       :: converged
      complex(kind=cp)                        :: eig(size(H,2))
      complex(kind=cp)                        :: shift
      complex(kind=cp)                        :: V(size(H,2),size(H,2))
      complex(kind=cp), allocatable           :: Q(:,:), R(:,:), S(:,:)
      integer                                 :: i, iter, m, n
      integer                                 :: maxit
      real(kind=rp), parameter                :: tol = epsilon(1._rp)
      integer                                 :: info
      logical                                 :: compute_U, compute_T

      n = size(H,2)
      V = H(1:n,1:n)
      iter = 0
      maxit = 100*n
      m = n 

      compute_U = present(U)
      compute_T = present(T)
!
      if ( compute_U ) then 
         allocate(S(n,n))
         U = 0._cp
         do i = 1, n
            U(i,i) = 1._cp
         end do
      end if

      shift = 0._cp
      do i = 1, m
         if ( abs(V(i,i)) > abs(shift) ) shift = V(i,i) 
      end do
      shift = shift*cmplx(1._cp,1._cp)
!
      do while ( m > 1 .and. iter < maxit ) 

!
! Apply shift:
         do i = 1, m
            V(i,i) = V(i,i) - shift
         end do

!
! QR-iteration step:
         allocate(Q(m,m), R(m,m))
         call QR_factorization( V(1:m,1:m), Q, R, info )
         if ( info > 0 ) exit
         V(1:m,1:m) = matmul( R, Q )
            
!
! Correct for shift:
         do i = 1, m
            V(i,i) = V(i,i) + shift
         end do

!
! For Schur vectors:
         if ( compute_U ) then
            S = 0._cp
            do i = 1, n
               S(i,i) = 1._cp
            end do
            S(1:m,1:m) = Q
            U = matmul( U, S )
         end if

!
! For Upper triangular matrix in the Schur decomposition:
         if ( compute_T ) then
            if ( m < n  ) V(1:m,m+1:n) = matmul( transpose(conjg(Q)), V(1:m,m+1:n) )
         end if

!
! Convergence of eigenvalue?
         if ( abs( V(m,m-1) ) < tol*abs(V(m,m) ) ) m = m - 1

!
! Prepare for new iteration:
         shift = 0._cp
         do i = 1, m
            if ( abs(V(i,i)) > abs(shift) ) shift = V(i,i) 
         end do
         iter = iter + 1
         deallocate(Q,R)

      end do

!
! Store the eigenvalues:
      do i = 1, n
         eig(i) = V(i,i)
      end do

      if ( present(converged) ) converged = n-m+1
      if ( compute_T ) T = V

   end function RQR_eig

   function complex_invariant( U, eig, dim ) result(V)
!
! Compute complex invariant subspace corresponding to dim smallest eigenvalues
!    In:   U, complex, size n*n          Schur vectors
!          eig, complex, size n          Eigenvalues
!          dim, integer                  Dimension invariant subspace
!    Out:  V, complex, size n*dim        Invariant subspace
!
      complex(kind=cp), intent(in)            :: U(:,:), eig(:)
      integer, intent(in)                     :: dim
      complex(kind=cp)                        :: V(size(U,1),dim)
      integer                                 :: i, n

      V = U(:,1:dim)

   end function complex_invariant

   function real_invariant( U, eig, dim ) result(V)
!
! Compute complex invariant subspace corresponding to dim smallest eigenvalues
!    In:   U, complex, size n*n          Schur vectors
!          eig, complex, size n          Eigenvectors
!          dim, integer                  Dimension invariant subspace
!    Out:  V, real, size n*dim           Invariant subspace
!
      complex(kind=cp), intent(in)            :: U(:,:), eig(:)
      integer, intent(in)                     :: dim
      real(kind=rp)                           :: V(size(U,1),dim)
      integer                                 :: i, n

      n = size(U,1)
      i = 0
      do while ( i < dim )
         i = i + 1
         V(:,i) = real(U(:,i))

         if ( abs(aimag(eig(i))) > 1.e-8*abs(real(eig(i))) ) then
! It is a complex eigenvalue
            if ( i < dim ) then
               i = i + 1
               V(:,i) = aimag(U(:,i))
            else
               write(*,'(a)') &
                 'Warning: no space to store imaginary part of eigenvector in invariant subspace.'
            end if
         end if
      end do

   end function real_invariant

   function CSI_power( H, s, shift ) result(schurvec)
!
! Compute Schur vectors of H closest to shift using block shift-invert power method
!    In:   H, complex, size ?*n          Matrix (note: only H(n,n) is used)
!          s, integer                    Block size
!          shift, complex, optional      shift
!    Out:  schurvec, complex, size n     Schur vector
!
      complex(kind=cp), intent(in)            :: H(:,:)
      integer, intent(in)                     :: s
      complex(kind=cp), optional, intent(in)  :: shift
      complex(kind=cp)                        :: A(size(H,2),size(H,2)), Q(size(H,2),size(H,2)), R(size(H,2),size(H,2))
      complex(kind=cp)                        :: schurvec(size(H,2),s)
      real(kind=rp)                           :: re_schurvec(size(H,2),s), im_schurvec(size(H,2),s)
      complex(kind=cp)                        :: Qs(size(H,2),s), Rs(s,s)
      real(kind=rp)                           :: trace, trace0, trace_old, err
      complex(kind=cp)                        :: sh
      integer, parameter                      :: maxit = 1000
      real(kind=rp), parameter                :: tol = 1e-3
      integer                                 :: n, i, iter
      
      n = size(H,2)

      if ( present(shift) ) then
         sh = shift
      else
         sh = 0.
      end if

! Shift the matrix
      A = H(1:n,1:n)
      do i = 1,n
         A(i,i) = A(i,i) - sh
      end do

! Compute QR-factorization:
      call QR_factorization( A, Q, R )

! Initial Power iteration to determine convergence criterion:
!   Generate initial value for the Schur vectors
      call random_number( re_schurvec )
      call random_number( im_schurvec )
      schurvec = cmplx( re_schurvec, im_schurvec,kind=cp )
      
!   Orthonormalise the vectors
      call QR_factorization( schurvec, Qs, Rs )
!   Change in the sum of the lenghts of the vectors
!   will be used as termination criterion. The R(i,i) are the lengths.
!   Determine sum of lengths of the vectors:
      trace0 = 0.
      do i = 1,s
         trace0 = trace0 + abs(Rs(i,i))
      end do
      trace_old = trace0

      iter = 0
      err = 1.
! Inverse Power iterations
      do while ( err > tol .and. iter < maxit )
         iter = iter + 1
! Solve shifted system:
         schurvec = QR_solve( Q, R, Qs )
! Orthonormalise the vectors:
         call QR_factorization( schurvec, Qs, Rs )
         trace = 0.
         do i = 1,s
            trace = trace + abs(Rs(i,i))
         end do
         err = abs( trace - trace_old ) / trace0
         trace_old = trace
      end do
      
      schurvec = Qs

   end function CSI_power

   function RSI_power( H, s, shift ) result(schurvec)
!
! Compute Schur vectors of H closest to shift using block shift-invert power method
!    In:   H, real, size ?*n             Matrix (note: only H(n,n) is used)
!          s, integer                    Block size
!          shift, real, optional         shift
!    Out:  schurvec, real, size n*s      Schur vector
!
      real(kind=rp), intent(in)               :: H(:,:)
      integer, intent(in)                     :: s
      real(kind=rp), optional, intent(in)     :: shift
      real(kind=rp)                           :: A(size(H,2),size(H,2)), Q(size(H,2),size(H,2)), R(size(H,2),size(H,2))
      real(kind=rp)                           :: schurvec(size(H,2),s)
      real(kind=rp)                           :: Qs(size(H,2),s), Rs(s,s)
      real(kind=rp)                           :: trace, trace0, trace_old, err
      real(kind=rp)                           :: sh
      integer, parameter                      :: maxit = 1000
      real(kind=rp), parameter                :: tol = 1e-3
      integer                                 :: n, i, iter
      
      n = size(H,2)

      if ( present(shift) ) then
         sh = shift
      else
         sh = 0.
      end if

! Shift the matrix
      A = H(1:n,1:n)
      do i = 1,n
         A(i,i) = A(i,i) - sh
      end do

! Compute QR-factorization:
      call QR_factorization( A, Q, R )

! Initial Power iteration to determine convergence criterion:
!   Generate initial value for the Schur vectors
      call random_number( schurvec )
!   Orthonormalise the vectors
      call QR_factorization( schurvec, Qs, Rs )
!   Change in the sum of the lenghts of the vectors 
!   will be used as termination criterion. The R(i,i) are the lengths.
!   Determine sum of lengths of the vectors:
      trace0 = 0.
      do i = 1,s
         trace0 = trace0 + Rs(i,i)
      end do
      trace_old = trace0

      iter = 0
      err = 1.
! Inverse Power iterations
      do while ( err > tol .and. iter < maxit )
         iter = iter + 1
! Solve shifted system:
         schurvec = QR_solve( Q, R, Qs )
! Orthonormalise the vectors:
         call QR_factorization( schurvec, Qs, Rs )
! Determine error
         trace = 0.
         do i = 1,s
            trace = trace + Rs(i,i)
         end do
         err = abs( trace - trace_old ) / trace0
         trace_old = trace
      end do
      
! Output:
      schurvec = Qs

   end function RSI_power

   function CSOLVE(A,b) result(x)
!
! Function to solve linear system
!
   IMPLICIT NONE
      complex(kind=cp), intent(in)     :: A(:,:), b(:)
      complex(kind=cp)                 :: x(size(b))
      integer                          :: i, j, ind(1), k, n, iter
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
      real(kind=rp)            :: rho
      real(kind=rp)            :: alpha

   if ( abs(a) < epsilon(1._rp) ) then
      c = 0.d0
      s = 1.d0
      r = b
   else
      rho = hypot(abs(a),abs(b))
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
      real(kind=rp)            :: rho
      complex(kind=cp)         :: alpha

   if ( abs(a) < epsilon(1._cp) ) then
      c = 0.
      s = 1.
      r = b
   else
      rho = hypot(abs(a),abs(b))
      alpha = a/abs(a)
      c = abs(a)/rho
      s = alpha*conjg(b)/rho
      r = alpha*rho
   end if
!
   end subroutine CGROT

end module dense_la_module

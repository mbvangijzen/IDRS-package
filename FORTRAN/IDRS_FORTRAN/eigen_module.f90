!
! The eigen_module contains a number of eigenvalue routines for small dense matrices.
! The algorithm are of standard QR and power type. The module is provided to make 
! the IDRS-package self-contained.
! 
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2024 Martin van Gijzen
! 

module eigen_module

   use precision_module

   implicit none

   interface QR_factorization
      module procedure CQR_factorization, RQR_factorization
   end interface

   interface QR_solve
      module procedure CQR_solve, RQR_solve, CBQR_solve, RBQR_solve
   end interface

   interface QR_eig
      module procedure CQR_eig, RQR_eig
   end interface

   interface SI_power
      module procedure CSI_power, RSI_power
   end interface

contains

   subroutine CQR_factorization( V, Q, R )

!
! Compute QR decomposition of complex matrix V
!    In:   V, complex, size m*n
!    Out:  Q, complex, orthogonal, size m*n
!    Out:  R, complex, upper triangular, size n*n
!
      complex(kind=cp), intent(in)      :: V(:,:)
      complex(kind=cp), intent(out)     :: Q(size(V,1),size(V,2)), R(size(V,2),size(V,2))
      integer                           :: i, j, m, n

      m = size(V,1)
      n = size(V,2)

      Q = V
      R = 0.
      do j = 1, n
         do i = 1, j-1
            R(i,j) = dot_product( Q(:,i), Q(:,j) )
            Q(:,j) = Q(:,j) - R(i,j)*Q(:,i)
         end do
         R(j,j) = sqrt(dot_product( Q(:,j), Q(:,j) ))
         if ( abs(R(j,j)) > tiny(1._rp) ) then 
            Q(:,j) = Q(:,j)/R(j,j)
         else
            print *, m, n, j
            stop ' V has dependent columns '
         end if
      end do

   end subroutine CQR_factorization

   subroutine RQR_factorization( V, Q, R )
 
!
! Compute QR decomposition of real matrix V
!    In:   V, real, size m*n
!    Out:  Q, real, orthogonal, size m*n
!    Out:  R, real, upper triangular, size n*n
!
      real(kind=rp), intent(in)         :: V(:,:)
      real(kind=rp), intent(out)        :: Q(size(V,1),size(V,2)), R(size(V,2),size(V,2))
      integer                           :: i, j, m, n

      m = size(V,1)
      n = size(V,2)

      Q = V
      R = 0.
      do j = 1, n
         do i = 1, j-1
            R(i,j) = dot_product( Q(:,i), Q(:,j) )
            Q(:,j) = Q(:,j) - R(i,j)*Q(:,i)
         end do
         R(j,j) = sqrt(dot_product( Q(:,j), Q(:,j) ))
         if ( R(j,j) > tiny(1._rp) ) then 
            Q(:,j) = Q(:,j)/R(j,j)
         else
            stop ' V has dependent columns '
         end if
      end do

   end subroutine RQR_factorization

   function CQR_solve( Q, R, b ) result(x)
!
! Compute solution of complex system Ax = b using QR factorization
!
! In:
!    A:    complex, size m*n
!    b:    complex, size m
! Out:  
!    x:    complex, size n
!
      complex(kind=cp), intent(in)         :: b(:)
      complex(kind=cp), intent(in)         :: Q(:,:), R(:,:)
      complex(kind=cp)                     :: x(size(Q,2))
      integer                              :: i, j, m, n

      m = size(Q,1)
      n = size(Q,2)

! Compute Q^Tb
      x = matmul(transpose(conjg(Q)),b)

! Backsubstitution:
      do j = n, 1, -1
         x(j) = x(j)/R(j,j)
         do i = 1, j-1
            x(i) = x(i) - R(i,j)*x(j)
         end do
      end do

   end function CQR_solve

   function RQR_solve( Q, R, b ) result(x)
!
! Compute solution of real system Ax = b using QR factorization
!
! In:
!    A:    complex, size m*n
!    b:    complex, size m
! Out:  
!    x:    complex, size n
!
      real(kind=rp), intent(in)         :: b(:)
      real(kind=rp), intent(in)         :: Q(:,:), R(:,:)
      real(kind=rp)                     :: x(size(Q,2))
      integer                           :: i, j, m, n

      m = size(Q,1)
      n = size(Q,2)

! Compute Q^Tb
      x = matmul(transpose(Q),b)

! Backsubstitution:
      do j = n, 1, -1
         x(j) = x(j) /R(j,j)
         do i = 1, j-1
            x(i) = x(i) - R(i,j)*x(j)
         end do
      end do

   end function RQR_solve

   function CBQR_solve( Q, R, B ) result(X)
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

   end function CBQR_solve

   function RBQR_solve( Q, R, B ) result(X)
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

   end function RBQR_solve

   function CQR_eig( H, converged ) result(eig)
!
! Compute eigenvalues of H using QR algorithm
!    In:   H, complex, size n*n          Matrix
!    Out:  eig, complex, size n          (converged) eigenvalues
!          converged, integer, optional  number of converged eigenvalues
!
      complex(kind=cp), intent(in)      :: H(:,:)
      integer, intent(out), optional    :: converged
      complex(kind=cp)                  :: eig(size(H,2))
      complex(kind=cp)                  :: shift
      complex(kind=cp)                  :: V(size(H,2),size(H,2))
      complex(kind=cp), allocatable     :: Q(:,:), R(:,:)
      integer                           :: i, iter, m, n, ieig
      integer                           :: maxit
      real(kind=rp), parameter          :: tol = 1e-16

      n = size(H,2)
      V = H(1:n,1:n)
      eig = 0.
      iter = 0
      ieig = 0
      maxit = 10*n
      m = n - ieig
      shift = V(m,m)

      do while ( ieig < n-1 .and. iter < maxit ) 
         iter = iter + 1
         do i = 1, m
            V(i,i) = V(i,i) - shift
         end do
         allocate( Q(m,m), R(m,m) )
         call QR_factorization( V(1:m,1:m), Q, R )
         V(1:m,1:m) = matmul( R, Q )
         deallocate( Q, R )
         do i = 1, m
            V(i,i) = V(i,i) + shift
         end do
         if ( abs( V(m,m-1) ) < tol*abs(V(m,m) ) ) then
            eig(ieig+1) = V(m,m)
            ieig = ieig + 1
            m = m - 1
         end if
         shift = V(m,m)
      end do
      if ( ieig == n-1 ) then
         eig(n) = V(1,1)
         ieig = ieig + 1
      end if
      if ( present(converged) ) converged = ieig

   end function CQR_eig

   function RQR_eig( H, converged ) result(eig)
!
! Compute eigenvalues of H using QR algorithm
!    In:   H, real, size n*n             Matrix
!    Out:  eig, complex, size n          (converged) eigenvalues
!          converged, integer, optional  number of converged eigenvalues
!
      real(kind=rp), intent(in)         :: H(:,:)
      integer, intent(out), optional    :: converged
      complex(kind=cp)                  :: eig(size(H,2))
      complex(kind=cp)                  :: shift
      complex(kind=cp)                  :: V(size(H,2),size(H,2))
      complex(kind=cp), allocatable     :: Q(:,:), R(:,:)
      integer                           :: i, iter, m, n, ieig
      integer                           :: maxit
      real(kind=rp), parameter          :: tol = 1e-16

      n = size(H,2)
      V = H(1:n,1:n)
      eig = 0.
      iter = 0
      ieig = 0
      maxit = 10*n
      m = n - ieig
      shift = V(n,n)*cmplx(1,1,kind=cp)         ! To break real arithmatic

      do while ( ieig < n-1 .and. iter < maxit ) 
         iter = iter + 1
         do i = 1, m
            V(i,i) = V(i,i) - shift
         end do
         allocate( Q(m,m), R(m,m) )
         call QR_factorization( V(1:m,1:m), Q, R )
         V(1:m,1:m) = matmul( R, Q )
         deallocate( Q, R )
         do i = 1, m
            V(i,i) = V(i,i) + shift
         end do
         if ( abs( V(m,m-1) ) < tol*abs(V(m,m) ) ) then
            eig(ieig+1) = V(m,m)
            ieig = ieig + 1
            m = m - 1
         end if
         shift = V(m,m)
      end do
      if ( ieig == n-1 ) then
         eig(n) = V(1,1)
         ieig = ieig + 1
      end if
      if ( present(converged) ) converged = ieig

   end function RQR_eig

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
      integer, parameter                      :: maxit = 100
      real(kind=rp), parameter                :: tol = 1e-8
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
      integer, parameter                      :: maxit = 100
      real(kind=rp), parameter                :: tol = 1e-8
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

end module eigen_module

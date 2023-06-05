module preconditioner_module

   use matrix_module

   implicit none

   type preconditioner

! Polynomial preconditioner:
      integer                        :: degree
      real(kind=rp), allocatable     :: ralpha_p(:)
      real(kind=rp), allocatable     :: ralpha_s(:,:)
      complex(kind=cp), allocatable  :: calpha_p(:)
      complex(kind=cp), allocatable  :: calpha_s(:,:)
! Matrix:
      type(matrix)                   :: A

   end type preconditioner

   INTERFACE OPERATOR(/)
      module procedure rprecon, cprecon
   END INTERFACE

   interface MAKE_PRECON
      module procedure RMAKE_PRECON, CMAKE_PRECON
   end interface

   interface CHEBYSHEV
      module procedure RCHEBYSHEV, CCHEBYSHEV
   end interface

   interface NEUMANN
      module procedure RNEUMANN, CNEUMANN
   end interface

   interface MAKE_MSPRECON
      module procedure RMAKE_MSPRECON, CMAKE_MSPRECON
   end interface

   interface SCALEBACK_PRECON
      module procedure CSCALEBACK_PRECON, RSCALEBACK_PRECON
   end interface

contains

   function rprecon( v, M ) result(w)

      type(preconditioner), intent(in)      :: M
      real(kind=rp), intent(in)             :: v(:)
      real(kind=rp)                         :: w(size(v))
      integer                               :: i

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

   end function rprecon

   function cprecon( v, M ) result(w)

      type(preconditioner), intent(in)      :: M
      complex(kind=cp), intent(in)          :: v(:)
      complex(kind=cp)                      :: w(size(v))
      integer                               :: i

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

   end function cprecon

   subroutine cmake_precon( M, A, degree, H )
!
! Compute the coeficients of the polynomial preconditioner
! from (partial) Hessenberg decomposition AV = VH of the matrix.
!

      type(preconditioner), intent(out)           :: M
      integer, intent(in)                         :: degree
      type(matrix), intent(in)                    :: A
      complex(kind=cp), intent(in)                :: H(:,:)
      complex(kind=cp)                            :: alpha_p(0:degree)
      complex(kind=cp)                            :: r(0:degree+1,0:degree+1)
      integer                                     :: i, j, i_sigma, n_sigma

      if ( degree == 0 ) then
         M%degree = degree
      else
         if ( size(H,1) < degree+2 ) stop "Hessenberg matrix too small for polynomial degree"
         if ( size(H,2) < degree+1 ) stop "Hessenberg matrix too small for polynomial degree"

! Store the information
         M%degree = degree
         M%A = A
         if ( allocated(M%calpha_p) ) deallocate(M%calpha_p)
         allocate(M%calpha_p(0:degree))

! Compute the parameters of the residual polynomial
         r = 0.
         r(0,degree+1) = 1.
         do i = 1, degree+1
            r(i,0:degree) = r(i-1,1:degree+1)
            do j = 1, i
               r(i,:) = r(i,:) - H(j,i)*r(j-1,:)
            end do
            r(i,:) = r(i,:)/H(i+1,i)
         end do

! Now compute the parameters of the polynomial preconditioner
! Note that r = 1 - tp(t), therefore p corresponds to minus coefficient 1:degree+1 of r
         M%calpha_p = -r(degree+1,degree:0:-1)

      end if

   end subroutine cmake_precon

   subroutine rmake_precon( M, A, degree, H )
!
! Compute the coeficients of the polynomial preconditioner
! from (partial) Hessenberg decomposition AV = VH of the matrix.
!

      type(preconditioner), intent(out)           :: M
      integer, intent(in)                         :: degree
      type(matrix), intent(in)                    :: A
      real(kind=rp), intent(in)                   :: H(:,:)
      real(kind=rp)                               :: r(0:degree+1,0:degree+1);
      integer                                     :: i, j

      if ( degree == 0 ) then
         M%degree = degree
      else

         if ( size(H,1) < degree+2 ) stop "Hessenberg matrix too small for polynomial degree"
         if ( size(H,2) < degree+1 ) stop "Hessenberg matrix too small for polynomial degree"

! Store the information
         M%degree = degree
         M%A = A
         if ( allocated(M%calpha_p) ) deallocate(M%calpha_p)
         allocate(M%ralpha_p(0:degree))

! Compute the parameters of the residual polynomial
         r = 0.
!        r(0,degree+1) = 1.
         r(0,0) = 1.
         do i = 1, degree+1
            r(i,1:degree+1) = r(i-1,0:degree) ! Polynomial becomes one degree higher r_i = t*r_i-1
            do j = 1, i
               r(i,:) = r(i,:) - H(j,i)*r(j-1,:)
            end do
            r(i,:) = r(i,:)/H(i+1,i)
         end do

! Now compute the parameters of the polynomial preconditioner
! Note that r = 1 - tp(t), therefore p corresponds to minus coefficient 1:degree+1 of r
         M%ralpha_p = -r(degree+1,1:degree+1)
!        M%ralpha_p = -r(degree+1,degree:0:-1)
         print *, r(degree+1,:)
         print * 
         print *, M%ralpha_p

      end if

   end subroutine rmake_precon

   subroutine cmake_msprecon( M, sigma, sigma_p  )
!
! Compute the coeficients of the shifted polynomial preconditioner
!
      type(preconditioner), intent(inout)      :: M
      complex(kind=cp), intent(in)             :: sigma(:)
      complex(kind=cp), intent(out)            :: sigma_p(:)
      integer                                  :: i, i_sigma, n_sigma
      integer                                  :: degree

      degree = M%degree
      if ( degree == 0 ) then
         sigma_p = sigma
      else
! Allocate space for the parameters for the shifted polynomials
         n_sigma = size(sigma)
         if ( .not. allocated(M%ralpha_s) ) allocate(M%calpha_s(0:degree,n_sigma))

! Compute the parameters of the shifted systems:
         M%calpha_s = 0.
         do i_sigma = 1,n_sigma
            M%calpha_s(degree,i_sigma) = M%calpha_p(degree)
            do i = degree-1,0,-1
               M%calpha_s(i,i_sigma) = M%calpha_p(i) + sigma(i_sigma)*M%calpha_s(i+1,i_sigma);
            end do
            sigma_p(i_sigma) = sigma(i_sigma)*M%calpha_s(0,i_sigma);
         end do
      end if

   end subroutine cmake_msprecon

   subroutine rmake_msprecon( M, sigma, sigma_p  )
!
! Compute the coeficients of the shifted polynomial preconditioner
!
      type(preconditioner), intent(inout)      :: M
      integer                                  :: degree
      real(kind=rp), intent(in)                :: sigma(:)
      real(kind=rp), intent(out)               :: sigma_p(:)
      integer                                  :: i, i_sigma, n_sigma

      degree = M%degree
      if ( degree == 0 ) then
         sigma_p = sigma
      else
! Allocate space for the parameters for the shifted polynomials
         n_sigma = size(sigma)
         if ( .not. allocated(M%ralpha_s) ) allocate(M%ralpha_s(0:degree,n_sigma))

! Compute the parameters of the shifted systems:
         M%ralpha_s = 0.
         do i_sigma = 1,n_sigma
            M%ralpha_s(degree,i_sigma) = M%ralpha_p(degree)
            do i = degree-1,0,-1
               M%ralpha_s(i,i_sigma) = M%ralpha_p(i) + sigma(i_sigma)*M%ralpha_s(i+1,i_sigma);
            end do
            sigma_p(i_sigma) = sigma(i_sigma)*M%ralpha_s(0,i_sigma);
         end do
      end if

   end subroutine rmake_msprecon

   function cscaleback_precon( M, y ) result(x)
!
! Scale back the solution
!
      type(preconditioner), intent(in)    :: M
      complex(kind=cp), intent(in)        :: y(:,:)
      complex(kind=cp)                    :: x(size(y,1),size(y,2))
      integer                             :: degree
      integer                             :: i, n_sigma, i_sigma

      degree = M%degree
      n_sigma = size(y,2)

      if ( degree > 0 ) then
         x = 0.
         do i_sigma = 1,n_sigma
            x(:,i_sigma) = M%calpha_s(degree,i_sigma)*y(:,i_sigma)
            do i = degree-1,0,-1
               x(:,i_sigma) = M%A*x(:,i_sigma) + M%calpha_s(i,i_sigma)*y(:,i_sigma)
            end do
         end do
      else
         x = y
      end if

   end function cscaleback_precon

   function rscaleback_precon( M, y ) result(x)
!
! Scale back the solution
!
      type(preconditioner), intent(in)    :: M
      real(kind=rp), intent(in)           :: y(:,:)
      real(kind=rp)                       :: x(size(y,1),size(y,2))
      integer                             :: degree
      integer                             :: i, n_sigma, i_sigma

      n_sigma = size(y,2)
      degree = M%degree

      if ( degree > 0 ) then
         x = 0.
         do i_sigma = 1,n_sigma
            x(:,i_sigma) = M%ralpha_s(degree,i_sigma)*y(:,i_sigma)
            do i = degree-1,0,-1
               x(:,i_sigma) = M%A*x(:,i_sigma) + M%ralpha_s(i,i_sigma)*y(:,i_sigma)
            end do
         end do
      else
         x = y
      end if

   end function rscaleback_precon

   function rneumann( degree, omega ) result(alpha)
!
! Rercursion for the residual Neumann polynomial
!         q_0 = 1
!         q_i = q_i-1 - omega t q_i-1
! The preconditioning polynomial is then given by
!         p = (q - 1)/t
!
      real(kind=rp)                            :: alpha(0:degree)
      integer, intent(in)                      :: degree
      real(kind=rp), intent(in)                :: omega
      real(kind=rp)                            :: q(0:degree+1,0:degree+1)

      integer                                  :: i

! First compute the residual polynomial q:
      q = 0.
      q(0,0) = 1.
      do i = 1, degree+1
         q(:,i) = q(:,i-1)
         q(1:degree+1,i) = q(1:degree+1,i) - omega*q(0:degree,i-1)
      end do
! Now get the coefficients of the preconditioning polynomial:
      alpha(0:degree) = -q(1:degree+1,degree+1)

   end function rneumann

   function cneumann( degree, omega ) result(alpha)
!
! Rercursion for the residual Neumann polynomial
!         q_0 = 1
!         q_i = q_i-1 - omega t q_i-1
! The preconditioning polynomial is then given by
!         p = (q - 1)/t
!
      complex(kind=rp)                         :: alpha(0:degree)
      integer, intent(in)                      :: degree
      complex(kind=cp), intent(in)             :: omega
      complex(kind=cp)                         :: q(0:degree+1,0:degree+1)

      integer                                  :: i

! First compute the residual polynomial q:
      q = 0.
      q(0,0) = 1.
      do i = 1, degree+1
         q(:,i) = q(:,i-1)
         q(1:degree+1,i) = q(1:degree+1,i) - omega*q(0:degree,i-1)
      end do
! Now get the coefficients of the preconditioning polynomial:
      alpha(0:degree) = -q(1:degree+1,degree+1)

   end function cneumann

!  function cneumann( degree, omega ) result(H)
!
! Put the recursion parameters for Neumann polynomials in Hessenberg matrix
!         r_i+1 = r_i - omega Ar_i, so Ar_i = 1/omega r_i - 1/omega r_i+1
!
!     integer, intent(in)                :: degree
!     complex(kind=cp), intent(in)       :: omega
!     complex(kind=cp)                   :: H(degree+2,degree+1)

!     integer                            :: i

!     H = 0.
!     do i = 1,degree+1
!        H(i,i)   =  1./omega
!        H(i+1,i) = -1./omega
!     end do

!  end function cneumann

   function rchebyshev( degree, f1, f2 ) result(alpha) 
!
! Generate polynomial preconditioner based on shifted and scaled Chebyshev polynomials
!
      real(kind=rp)                            :: alpha(0:degree)
      integer, intent(in)                      :: degree
      real(kind=rp), intent(in)                :: f1, f2
      real(kind=rp)                            :: theta, delta, sigma, rho
      real                                     :: p(0:degree,-1:degree)

      integer                                  :: i

      theta = (f1+f2)/2.
      delta = (f2-f1)/2.
      sigma = theta/delta
      p     = 0.
      p(0,0)  = 1./theta

      do i = 1, degree
         p(1:i,i) = -2.*sigma/theta*p(0:i-1,i-1)
         p(:,i) = p(:,i) + 2.*sigma*p(:,i-1) - rho*p(:,i-2)
         p(0,i) = p(0,i) + 2./delta
         rho = 1./(2*sigma-rho)
         p(:,i) = p(:,i)*rho
      end do
      alpha = p(:,degree)

   end function rchebyshev

   function cchebyshev( degree, f1, f2 ) result(alpha)
!
! Generate polynomial preconditioner based on shifted and scaled Chebyshev polynomials
!
      complex(kind=cp)                         :: alpha(0:degree)
      integer, intent(in)                      :: degree
      complex(kind=cp), intent(in)             :: f1, f2
      complex(kind=cp)                         :: theta, delta, sigma, rho, rho_old
      complex(kind=cp)                         :: p(0:degree,-1:degree)

      integer                                  :: i

      theta = (f1+f2)/2.
      delta = (f2-f1)/2.
      sigma = theta/delta
      p     = 0.
      p(0,0)  = 1./theta

      do i = 1, degree
         p(1:i,i) = -2.*sigma/theta*p(0:i-1,i-1)
         p(:,i) = p(:,i) + 2.*sigma*p(:,i-1) - rho*p(:,i-2)
         p(0,i) = p(0,i) + 2./delta
         rho = 1./(2*sigma-rho)
         p(:,i) = p(:,i)*rho
      end do
      alpha = p(:,degree)

   end function cchebyshev

end module preconditioner_module

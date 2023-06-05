!
! The idrs module uses a user defines module to define matrix-vector product 
! and preconditioning operation. This module gives an example for a dense
! matrix-vector product, and no preconditioning (so the preconditioner is just the 
! identity matrix.
!
! The matrix module should always contain the following four functions:
!
!    rmatvec:       real matrix-vector product
!    function call: xr = rmatvec( A, yr )
!    data types:    xr, yr: real of kind rp
!                   A: matrix (user defined data structure)
!
!    cmatvec:       complex matrix-vector product
!    function call: xc = cmatvec( A, yc )
!    data types:    xc, yc: complex of kind cp
!                   A: matrix (user defined data structure)
!                  
! The functions rmatvex and cmatvec should be overloaded to "*" to define the 
! matrix-vector operation x = A*y, where x and y can be real or complex.
!     
! The functions rprecon and cprecon should be overloaded to "/" to define the 
! matrix-vector operation x = P^-1y, where x and y can be real or complex.
!
! The kind parameters rp and cp also need to be set in the matrix module.
!

module matrix_module

   implicit none

! Define the precision for real and complex arithmatic

   integer, parameter :: rp = kind(1.d0)               ! kind for real
   integer, parameter :: cp = kind(1.d0)               ! kind for complex

! Define the matrix type, as an example here for dense matrices

   type matrix
      real(kind=rp), allocatable     :: Ar(:,:)        ! Real matrix
      complex(kind=cp), allocatable  :: Ac(:,:)        ! Complex matrix
   end type matrix

! Overload * to define the matrix-vector multiplication using the matrix type

   INTERFACE OPERATOR(*)
      module procedure rmatvec, cmatvec
   END INTERFACE

   contains

   function rmatvec( A, v ) 

      type(matrix), intent(in)          :: A
      real(kind=rp), intent(in)         :: v(:)
      real(kind=rp)                     :: rmatvec(size(v,1))
      
      rmatvec = matmul( A%Ar, v )

   end function rmatvec

   function cmatvec( A, v ) 

      type(matrix), intent(in)          :: A
      complex(kind=cp), intent(in)      :: v(:)
      complex(kind=cp)                  :: cmatvec(size(v,1))
      
      cmatvec = matmul( A%Ac, v )

   end function cmatvec

end module matrix_module

!
! The idrs module uses a user defines module to define matrix-vector product 
! and preconditioning operation. This module gives an example for a dense
! matrix-vector product, and no preconditioning (so the preconditioner is just the 
! identity matrix.
!
! The matrix module should always contain the following functions:
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
      real(kind=rp)   :: c, e, w, n, s, b, f
      real(kind=rp)   :: hx, hy, hz
      integer         :: Nx, Ny, Nz, Nn
      integer         :: Sx, Sy, Sz
      integer         :: me(3)
   end type matrix

! Overload * to define the matrix-vector multiplication using the matrix type

   INTERFACE OPERATOR(*)
      module procedure rmatvec, cmatvec
   END INTERFACE

   contains

   function rmatvec( A, v ) result(w)

      type(matrix), intent(in)       :: A
      real(kind=rp), intent(in)      :: v(A%Nn) 
      real(kind=rp)                  :: w(A%Nn)
      real(kind=rp), allocatable,save:: v_grid(:,:,:)[:,:,:] 
      real(kind=rp)                  :: w_grid(1:A%Nx,1:A%Ny,1:A%Nz)
      integer                        :: i,j,k,l, me(3)
      
! Allocate co-array. Note that this performs a global synchronization.
      sync all
      if ( .not. allocated(v_grid) ) allocate(v_grid(0:A%Nx+1,0:A%Ny+1,0:A%Nz+1)[A%Sx,A%Sy,*])
      v_grid = 0._rp
      do concurrent ( i=1:A%Nx, j=1:A%Ny, k=1:A%Nz )
         l = i +(j-1)*A%nx +(k-1)*A%nx*A%Ny
         v_grid(i,j,k) = v(l)
      end do
      me = A%me
!
! Fill v_grid with information from neigbouring domains:
! Now put information for the adjacent nodes in the six neighbours.
! West:
!     if ( me(1) > 1 ) v_grid(0,1:A%Ny,1:A%Nz) = v_grid(A%Nx,1:A%Ny,1:A%Nz)[me(1)-1,me(2),me(3)]
      if ( me(1) > 1 ) v_grid(A%Nx+1,1:A%Ny,1:A%Nz)[me(1)-1,me(2),me(3)] = v_grid(1,1:A%Ny,1:A%Nz) 
! East:
!     if ( me(1) < A%Sx ) v_grid(A%Nx+1,1:A%Ny,1:A%Nz) = v_grid(1,1:A%Ny,1:A%Nz)[me(1)+1,me(2),me(3)]
      if ( me(1) < A%Sx ) v_grid(0,1:A%Ny,1:A%Nz)[me(1)+1,me(2),me(3)] = v_grid(A%Nx,1:A%Ny,1:A%Nz) 
! South:
!     if ( me(2) > 1 ) v_grid(1:A%Nx,0,1:A%Nz) = v_grid(1:A%Nx,A%Ny,1:A%Nz)[me(1),me(2)-1,me(3)]
      if ( me(2) > 1 ) v_grid(1:A%Nx,A%Ny+1,1:A%Nz)[me(1),me(2)-1,me(3)] = v_grid(1:A%Nx,1,1:A%Nz) 
! North:
!     if ( me(2) < A%Sy ) v_grid(1:A%Nx,A%Ny+1,1:A%Nz) = v_grid(1:A%Nx,1,1:A%Nz)[me(1),me(2)+1,me(3)]
      if ( me(2) < A%Sy ) v_grid(1:A%Nx,0,1:A%Nz)[me(1),me(2)+1,me(3)] = v_grid(1:A%Nx,A%Ny,1:A%Nz) 
! Front:
!     if ( me(3) > 1 ) v_grid(1:A%Nx,1:A%Ny,0) = v_grid(1:A%Nx,1:A%Ny,A%Nz)[me(1),me(2),me(3)-1]
      if ( me(3) > 1 ) v_grid(1:A%Nx,1:A%Ny,A%Nz+1)[me(1),me(2),me(3)-1] = v_grid(1:A%Nx,1:A%Ny,1) 
! Back:
!     if ( me(3) < A%Sz ) v_grid(1:A%Nx,1:A%Ny,A%Nz+1) = v_grid(1:A%Nx,1:A%Ny,1)[me(1),me(2),me(3)+1]
      if ( me(3) < A%Sz ) v_grid(1:A%Nx,1:A%Ny,0)[me(1),me(2),me(3)+1] = v_grid(1:A%Nx,1:A%Ny,A%Nz) 
      sync all

      do concurrent ( i=1:A%Nx, j=1:A%Ny, k=1:A%Nz )
         w_grid(i,j,k) =  A%c*v_grid(i,j,k) + &
                          A%w*v_grid(i-1,j,k) + A%e*v_grid(i+1,j,k) + &
                          A%s*v_grid(i,j-1,k) + A%n*v_grid(i,j+1,k) + &
                          A%f*v_grid(i,j,k-1) + A%b*v_grid(i,j,k+1) 
      end do
      w = reshape(w_grid,(/A%Nn/))

   end function rmatvec

   function cmatvec( A, v ) result(w) 

      type(matrix), intent(in)          :: A
      complex(kind=cp), intent(in)      :: v(:)
      complex(kind=cp)                  :: w(size(v,1))
      
      w = v

   end function cmatvec

end module matrix_module

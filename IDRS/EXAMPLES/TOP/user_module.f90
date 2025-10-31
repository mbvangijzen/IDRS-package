!
! This module contains problem dependent routines, to be defined by
! the user. The present module contains routines that help make
! a polynomial or Jacobi preconditioner
!
! In particular, it contains the following routines and functions:
!
!    ruser_mv: real matrix vector-multiplication
!    function call: w = ruser_mv( A, v )
!
!    cuser_mv: complex matrix-vector multiplication
!    function call: w = cuser_mv( A )
!    data types:    A: matrix (user defined data structure)
!                   v, w: complex arrays
!
!    The "*" operator is overloaded to allow w = A*v
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2024 Martin van Gijzen
!
!
module user_module

   use precision_module

   implicit none

! Define the user matrix type

   type user_matrix
      integer                    :: nelx, nely
      real(kind=rp)              :: Ke(8,8)
      real(kind=rp), allocatable :: rho(:,:)
      real(kind=rp)              :: penal
      logical, allocatable       :: fixed(:)
      real(kind=rp), allocatable :: D(:)
   end type user_matrix

! Overload * to define the matrix-vector multiplication using the matrix type

   INTERFACE OPERATOR(*)
      module procedure ruser_mv, cuser_mv
   END INTERFACE

   INTERFACE OPERATOR(/)
      module procedure ruser_precon, cuser_precon
   END INTERFACE

   contains

   function ruser_mv( A, v ) result(w)

      implicit none
      type(user_matrix), intent(in)  :: A
      real(kind=rp), intent(in)      :: v(:) 
      real(kind=rp)                  :: w(size(v,1))
      integer                        :: n1, n2, elx, ely
      real(kind=rp)                  :: ve(8), we(8)
      integer                        :: ind(8)

      w = 0.
      do ely = 1,A%nely
         do elx = 1,A%nelx
            n1 = (A%nely+1)*(elx-1)+ely
            n2 = (A%nely+1)* elx   +ely
! Indices of element dofs:
            ind = [2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n2+1, 2*n2+2, 2*n1+1,2*n1+2]
! Gather the element vector:
            ve = v(ind)
! Element matrix-vector multiplication:
            we = matmul( A%Ke, ve )
! Scatter, put the element vector in global vector
            w(ind) = w(ind) + (A%rho(ely,elx)**A%penal)*we
         end do
      end do
! Correct for boundary condition
      where ( A%fixed )
         w = 0.
      end where
      
   end function ruser_mv

   function cuser_mv( A, v ) result(w) 

      implicit none
      type(user_matrix), intent(in)  :: A
      complex(kind=cp), intent(in)   :: v(:)
      complex(kind=cp)               :: w(size(v))
      
      w = v

   end function cuser_mv

   function ruser_precon( v, P ) result(w) 

      implicit none
      type(user_matrix), intent(in)  :: P
      real(kind=rp), intent(in)      :: v(:)
      real(kind=rp)                  :: w(size(v))
      
      if ( allocated(P%D) ) then
         w = v/P%D
      else
         w = v
      end if

   end function ruser_precon

   function cuser_precon( v, P ) result(w) 

      implicit none
      type(user_matrix), intent(in)  :: P
      complex(kind=cp), intent(in)   :: v(:)
      complex(kind=cp)               :: w(size(v))
      
      w = v

   end function cuser_precon

   function real_diagonal( A ) result(D)

      implicit none
      type(user_matrix), intent(in)         :: A
      real(kind=rp)                         :: D(2*(A%nelx+1)*(A%nely+1))
      integer                               :: ind(8)
      integer                               :: elx, ely, n1, n2, i

      D = 0._rp
      do ely = 1,A%nely
         do elx = 1,A%nelx
            n1 = (A%nely+1)*(elx-1)+ely
            n2 = (A%nely+1)* elx   +ely
! Get the element diagonal:
            ind = [2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n2+1, 2*n2+2, 2*n1+1, 2*n1+2]
            do i = 1, 8
               D(ind(i)) = D(ind(i)) + (A%rho(ely,elx)**A%penal)*A%Ke(i,i)
            end do
         end do
      end do
   end function real_diagonal

   function complex_diagonal( A ) result(D)

      implicit none
      type(user_matrix), intent(in)         :: A
      complex(kind=cp)                      :: D(2*(A%nelx+1)*(A%nely+1))
      integer                               :: ind(8)
      integer                               :: elx, ely, n1, n2, i

      D = 0._cp
      do ely = 1,A%nely
         do elx = 1,A%nelx
            n1 = (A%nely+1)*(elx-1)+ely
            n2 = (A%nely+1)* elx   +ely
! Get the element diagonal:
            ind = [2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n2+1, 2*n2+2, 2*n1+1, 2*n1+2]
            do i = 1, 8
               D(ind(i)) = D(ind(i)) + (A%rho(ely,elx)**A%penal)*A%Ke(i,i)
            end do
         end do
      end do
   end function complex_diagonal

end module user_module

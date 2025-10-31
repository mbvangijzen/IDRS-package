module matrix_module

   use precision_module
   implicit none

! Define the matrix type, as an example here for Toeplitz matrices

   type matrix
      integer                        :: rind(3)          ! Diagonal indices
      real(kind=rp)                  :: rval(3)          ! Diagonal values
      integer                        :: cind(4)          ! Diagonal indices
      complex(kind=cp)               :: cval(4)          ! Diagonal values
   end type matrix

! Overload * to define the matrix-vector multiplication using the matrix type

   INTERFACE OPERATOR(*)
      module procedure rmatvec, cmatvec
   END INTERFACE

! Define the preconditioner type 

   type preconditioner
      integer                        :: precon
   end type preconditioner

! Overload / to define the preconditioner operation using the preconditioner type

   INTERFACE OPERATOR(/)
      module procedure rprecon, cprecon
   END INTERFACE

   contains

   function rmatvec( A, v ) result(w)

      type(matrix), intent(in)        :: A
      real(kind=rp), intent(in)       :: v(:)
      real(kind=rp)                   :: w(size(v))
      integer                         :: i,j,N, ndiag, jdiag
      real(kind=rp)                   :: val

      N     = size(v)
      ndiag = size(A%rind)
      w = 0.
      do j = 1, ndiag
         jdiag = A%rind(j)
         val   = A%rval(j)
         if ( jdiag < 0 ) then
            do i = 1,N+jdiag
               w(i-jdiag) = w(i-jdiag) + val*v(i)
            end do
         else
            do i = 1,N-jdiag
               w(i) = w(i) + val*v(i+jdiag)
            end do
         end if
      end do

   end function rmatvec

   function cmatvec( A, v ) result(w)

      type(matrix), intent(in)        :: A
      complex(kind=cp), intent(in)    :: v(:)
      complex(kind=cp)                :: w(size(v))
      integer                         :: i,j,N, ndiag, jdiag
      complex(kind=cp)                :: val

      N    = size(v,1)
      ndiag = size(A%cind)
      w = 0.
      do j = 1, ndiag
         jdiag = A%cind(j)
         val   = A%cval(j)
         if ( jdiag < 0 ) then
            do i = 1,N+jdiag
               w(i-jdiag) = w(i-jdiag) + val*v(i)
            end do
         else
            do i = 1,N-jdiag
               w(i) = w(i) + val*v(i+jdiag)
            end do
         end if
      end do

   end function cmatvec

   function rprecon( v, M )

      type(matrix), intent(in)              :: M
      real(kind=rp), intent(in)             :: v(:)
      real(kind=rp)                         :: rprecon(size(v))

      rprecon = v

   end function rprecon

   function cprecon( v, M )

      type(matrix), intent(in)              :: M
      complex(kind=cp), intent(in)          :: v(:)
      complex(kind=cp)                      :: cprecon(size(v))

      cprecon = v

   end function cprecon

end module matrix_module

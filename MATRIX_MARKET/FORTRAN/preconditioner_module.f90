!
! The idrs module uses a user defines module to define the
! preconditioning operation. The matrix should be given
! in matrix market format. The preconditioner that is 
! computed is is the main diagonal of the matrix
!
! The matrix module contains the following two functions:
!
!    rprecon:       real preconditioning operation
!    function call: xr = rprecon( M, yr )
!    data types:    xr, yr: real of kind rp
!                   M: preconditioner (user defined data structure)
!
!    cprecon:       complex matrix-vector product
!    function call: xc = cprecon( A, yc )
!    data types:    xc, yc: complex of kind cp
!                   M: preconditiponer (user defined data structure)
!
! The functions rprecon and cprecon are overloaded to "/" to define the
! matrix-vector operation x = M^-1y, where x and y can be real or complex.
!

module preconditioner_module

   use matrix_module
   implicit none

! Define the preconditioner type

   type preconditioner
      real(kind=rp),    allocatable  :: RM(:)
      complex(kind=cp), allocatable  :: CM(:)
      character(len=6)               :: which_preconditioner
   end type preconditioner

! Overload / to define the preconditioner operation using the preconditioner type

   INTERFACE OPERATOR(/)
      module procedure rprecon, cprecon
   END INTERFACE

   contains

   function rprecon( v, M )

      type(preconditioner), intent(in)      :: M
      real(kind=rp), intent(in)             :: v(:)
      real(kind=rp)                         :: rprecon(size(v)), w(size(v))
      integer                               :: i

      if ( M%which_preconditioner == 'none' ) then
         w = v
      elseif ( M%which_preconditioner == 'jacobi' ) then
         w = v/M%RM
      end if
      rprecon = w

   end function rprecon

   function cprecon( v, M )

      type(preconditioner), intent(in)      :: M
      complex(kind=cp), intent(in)          :: v(:)
      complex(kind=cp)                      :: cprecon(size(v)), w(size(v))
      integer                               :: i

      if ( M%which_preconditioner == 'none' ) then
         w = v
      else if ( M%which_preconditioner == 'jacobi' ) then
         w = v/M%CM
      end if
      cprecon = w

   end function cprecon

   subroutine make_preconditioner( A, M )

      integer                               :: prec
      type(matrix), intent(in)              :: A
      type(preconditioner), intent(out)     :: M
      integer                               :: neq, nnz, i, j, k, err

      neq = A%rows
      nnz = A%nnz
      if ( (A%rep == 'array') .or. (A%field == 'integer') .or. (A%field == 'pattern') ) then
         M%which_preconditioner = 'none' 
      else
         M%which_preconditioner = 'jacobi' 

         if  ( A%field == 'real' ) then
            allocate( M%RM(neq), stat = err ) 
            M%RM = 0.
            do k = 1,nnz
               i = A%indx(k)
               j = A%jndx(k)
               if ( i == j ) then
                  M%RM(i) = M%RM(i) + A%rval(k)
               end if
            end do
         elseif  ( A%field == 'complex' ) then
            allocate( M%CM(neq), stat = err ) 
            M%CM = 0.
            do k = 1,nnz
               i = A%indx(k)
               j = A%jndx(k)
               if ( i == j ) then
                  M%CM(i) = M%CM(i) + A%cval(k)
               end if
            end do
         end if

! Check if there are zero elements on the main diagonal:
!
         if ( A%field == 'complex' ) then
            do i = 1, neq
               if ( M%CM(i) == 0. ) then
                  M%which_preconditioner = 'none'
                  exit
               end if
            end do 
         else
            do i = 1, neq
               if ( M%RM(i) == 0. ) then
                  M%which_preconditioner = 'none'
                  exit
               end if
            end do 
         end if
                  
      end if

   end subroutine make_preconditioner

end module preconditioner_module

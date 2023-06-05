!
! The idrs module uses a user defines module to define matrix-vector product
! and preconditioning operation. This module implements a matrix-vector product
! for matrices in matrix-market format. 
!
! The matrix module contains the following two functions:
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
! The functions rmatvex and cmatvec are overloaded to "*" to define the
! matrix-vector operation x = A*y, where x and y can be real or complex.
!
! The kind parameters rp and cp also need to be set in the matrix module.
!

module matrix_module

   implicit none

! Define the precision for real and complex arithmatic

   integer, parameter :: rp = kind(1.d0)               ! kind for real
   integer, parameter :: cp = kind(1.d0)               ! kind for complex

! Define the matrix type, here for matrices in matrix-market format

   type matrix
      real(kind=rp),    allocatable  :: rval(:)
      integer,          allocatable  :: ival(:)
      complex(kind=cp), allocatable  :: cval(:)
      integer, allocatable           :: indx(:)
      integer, allocatable           :: jndx(:)
      real(kind=rp),    allocatable  :: rmat(:,:)
      integer,          allocatable  :: imat(:,:)
      complex(kind=cp), allocatable  :: cmat(:,:)
      integer                        :: rows, cols, nnz
      character(len=7)               :: field
      character(len=10)              :: rep
   end type matrix

! Overload * to define the matrix-vector multiplication using the matrix type

   INTERFACE OPERATOR(*)
      module procedure rmatvec, cmatvec, imatvec
   END INTERFACE

   contains

   subroutine read_mm_matrix( file, A )

      implicit none 
      integer, intent(in)           :: file
      type(matrix), intent(out)     :: A

      integer                       :: cols, rows, nnz, nnzmax 
      integer                       :: Allocate_status     
      integer, allocatable          :: ival(:), indx(:), jndx(:)
      double precision, allocatable :: rval(:)
      complex, allocatable          :: cval(:)
      character(len=10)             :: rep
      character(len=7)              :: field
      character(len=19)             :: symm
!
      integer                       :: i, j, k

!  
! Read header to determine which space is needed
      call mminfo(file,rep,field,symm,rows,cols,nnzmax)
!
! Store the information
      A%field = field
      A%rows = rows
      A%cols = cols
      A%rep = rep
!
! Allocate the space
      allocate( indx(nnzmax), jndx(nnzmax), stat = Allocate_status)
      allocate( rval(nnzmax), stat = Allocate_status)
      allocate( ival(nnzmax), stat = Allocate_status)
      allocate( cval(nnzmax), stat = Allocate_status)
      call mmread(file,rep,field,symm,rows,cols,nnz,nnzmax, indx,jndx,ival,rval,cval )
!
      if ( rep == 'coordinate' ) then
         if ( symm == 'general' ) then
            allocate( A%indx(nnz), A%jndx(nnz) )
            k = 1
            do i = 1, nnz
               A%indx(k) = indx(i)
               A%jndx(k) = jndx(i)
               k = k + 1
            end do
         else
            allocate( A%indx(2*nnz), A%jndx(2*nnz) )
            k = 1
            do i = 1, nnz
               A%indx(k) = indx(i)
               A%jndx(k) = jndx(i)
               k = k + 1
               if ( indx(i) .ne. jndx(i) ) then
                  A%indx(k) = jndx(i)
                  A%jndx(k) = indx(i)
                  k = k + 1
               end if
            end do
         end if
         A%nnz = k-1
!
! Store the matrix coefficients:
         if ( field == 'real' ) then
            allocate( A%rval(A%nnz) )
            if ( symm == 'general' ) then
               k = 1
               do i = 1, nnz
                  A%rval(k) = rval(i)
                  k = k + 1
               end do
            elseif ( symm == 'symmetric' ) then
               k = 1
               do i = 1, nnz
                  A%rval(k) = rval(i)
                  k = k + 1
                  if ( indx(i) .ne. jndx(i) ) then
                     A%rval(k) = rval(i)
                     k = k + 1
                  end if
               end do
            elseif ( symm == 'skew-symmetric' ) then
               k = 1
               do i = 1, nnz
                  A%rval(k) = rval(i)
                  k = k + 1
                  A%rval(k) =  -rval(i)
                  k = k + 1
               end do
            end if
         elseif ( field == 'integer' ) then
            allocate( A%ival(A%nnz) )
            if ( symm == 'general' ) then
               k = 1
               do i = 1, nnz
                  A%ival(k) = ival(i)
                  k = k + 1
               end do
            elseif ( symm == 'symmetric' ) then
               k = 1
               do i = 1, nnz
                  A%ival(k) = ival(i)
                  k = k + 1
                  if ( indx(i) .ne. jndx(i) ) then
                     A%ival(k) = ival(i)
                     k = k + 1
                  end if
               end do
            elseif ( symm == 'skew-symmetric' ) then
               k = 1
               do i = 1, nnz
                  A%ival(k) = ival(i)
                  k = k + 1
                  A%ival(k) = -ival(i)
                  k = k + 1
               end do
            end if
         elseif ( field == 'pattern' ) then
            allocate( A%ival(A%nnz) )
            if ( symm == 'general' ) then
               k = 1
               do i = 1, nnz
                  A%ival(k) = 1
                  k = k + 1
               end do
            elseif ( symm == 'symmetric' ) then
               k = 1
               do i = 1, nnz
                  A%ival(k) = 1
                  k = k + 1
                  if ( indx(i) .ne. jndx(i) ) then
                     A%ival(k) = 1
                     k = k + 1
                  end if
               end do
            elseif ( symm == 'skew-symmetric' ) then
               k = 1
               do i = 1, nnz
                  A%ival(k) = 1
                  k = k + 1
                  A%ival(k) = -1
                  k = k + 1
               end do
            end if

         elseif ( field == 'complex' ) then
            allocate( A%cval(A%nnz) )
            if ( symm == 'general' ) then
               k = 1
               do i = 1, nnz
                  A%cval(k) = cval(i)
                  k = k + 1
               end do
            elseif ( symm == 'symmetric' ) then
               k = 1
               do i = 1, nnz
                  A%cval(k) = cval(i)
                  k = k + 1
                  if ( indx(i) .ne. jndx(i) ) then
                     A%cval(k) = cval(i)
                     k = k + 1
                  end if
               end do
            elseif ( symm == 'skew-symmetric' ) then
               k = 1
               do i = 1, nnz
                  A%cval(k) = cval(i)
                  k = k + 1
                  A%cval(k) =  -cval(i)
                  k = k + 1
               end do
            elseif ( symm == 'hermitian' ) then
               k = 1
               do i = 1, nnz
                  A%cval(k) = cval(i)
                  k = k + 1
                  if ( indx(i) .ne. jndx(i) ) then
                     A%cval(k) = conjg(cval(i))
                     k = k + 1
                  end if
               end do

            end if
         end if
      else
! Dense matrix
         if ( field == 'real' ) then
            allocate(A%rmat(rows,cols))
            A%rmat = 0.
            if ( symm == 'general' ) then 
               k = 1
               do j = 1,cols
                  A%rmat(1:rows,j) = rval(k:k+rows-1)
                  k = k + rows
               end do
            elseif ( symm == 'symmetric' ) then
               k = 1
               do j = 1,cols
                  A%rmat(j:rows,j)   = rval(k:k+rows-j)
                  A%rmat(j,j+1:rows) = rval(k+1:k+rows-j)
                  k = k + rows-j
               end do
            elseif ( symm == 'skew-symmetric' ) then
               k = 1
               do j = 1,cols-1
                  A%rmat(j+1:rows,j) =  rval(k:k+rows-j-1)
                  A%rmat(j,j+1:rows) = -rval(k:k+rows-j-1)
                  k = k + rows-j-1
               end do
            end if
         elseif ( field == 'integer' ) then
            allocate(A%imat(rows,cols))
            A%imat = 0
            if ( symm == 'general' ) then 
               k = 1
               do j = 1,cols
                  A%imat(1:rows,j) = ival(k:k+rows-1)
                  k = k + rows
               end do
            elseif ( symm == 'symmetric' ) then
               k = 1
               do j = 1,cols
                  A%imat(j:rows,j)   = ival(k:k+rows-j)
                  A%imat(j,j+1:rows) = ival(k+1:k+rows-j)
                  k = k + rows-j
               end do
            elseif ( symm == 'skew-symmetric' ) then
               k = 1
               do j = 1,cols-1
                  A%imat(j+1:rows,j) =  ival(k:k+rows-j-1)
                  A%imat(j,j+1:rows) = -ival(k:k+rows-j-1)
                  k = k + rows-j-1
               end do
            end if
         elseif ( field == 'complex' ) then
            allocate(A%cmat(rows,cols))
            A%cmat = (0.,0.)
! General dense complex matrix:
            if ( symm == 'general' ) then 
               k = 1
               do j = 1,cols
                  A%cmat(1:rows,j) = cval(k:k+rows-1)
                  k = k + rows
               end do
            elseif ( symm == 'symmetric' ) then
               k = 1
               do j = 1,cols
                  A%cmat(j:rows,j)   = cval(k:k+rows-j)
                  A%cmat(j,j+1:rows) = cval(k+1:k+rows-j)
                  k = k + rows-j
               end do
            elseif ( symm == 'skew-symmetric' ) then
               k = 1
               do j = 1,cols-1
                  A%cmat(j+1:rows,j) =  cval(k:k+rows-j-1)
                  A%cmat(j,j+1:rows) = -cval(k:k+rows-j-1)
                  k = k + rows-j-1
               end do
            end if
            if ( symm == 'hermitian' ) then 
               k = 1
               do j = 1,cols
                  A%cmat(j:rows,j)   = cval(k:k+rows-j)
                  A%cmat(j,j+1:rows) = conjg(cval(k+1:k+rows-j))
                  k = k + rows-j
               end do
            end if
         end if
      end if
!
   end subroutine read_mm_matrix 

   function rmatvec( A, v ) result(w)

      type(matrix), intent(in)        :: A
      real(kind=rp), intent(in)       :: v(:)
      real(kind=rp)                   :: w(A%rows)
      integer                         :: nnz, j

      nnz = A%nnz

      if ( A%rep == 'array' ) then
!
! General real matrix, dense matrix-vector product:
         w = matmul( A%rmat, v )
!
      else
!
! General real matrix, sparse matrix-vector product:
         w = 0.
         do j = 1,nnz
            w(A%indx(j)) = w(A%indx(j)) + A%rval(j)*v(A%jndx(j))
         end do
!
      end if

   end function rmatvec

   function imatvec( A, v ) result(w)

      type(matrix), intent(in)        :: A
      integer, intent(in)             :: v(:)
      integer                         :: w(A%rows)
      integer                         :: nnz, j, neq

      nnz = A%nnz

      if ( A%rep == 'array' ) then
!
! Integer matrix, dense matrix-vector product:
         w = matmul( A%imat, v )
!
      else
!
! Integer matrix, sparse matrix-vector product:
         w = 0.
         do j = 1,nnz
            w(A%indx(j)) = w(A%indx(j)) + A%ival(j)*v(A%jndx(j))
         end do
!
      end if

   end function imatvec

   function cmatvec( A, v ) result(w)

      type(matrix), intent(in)        :: A
      complex(kind=cp), intent(in)    :: v(:)
      complex(kind=cp)                :: w(A%rows)
      integer                         :: nnz, j

      nnz = A%nnz

      if ( A%rep == 'array' ) then
!
! General complex matrix, dense matrix-vector product:
         w = matmul( A%cmat, v )
!
      else
!
! General complex matrix, sparse matrix-vector product:
         w = 0.
         do j = 1,nnz
            w(A%indx(j)) = w(A%indx(j)) + A%cval(j)*v(A%jndx(j))
         end do
!
      end if

   end function cmatvec

end module matrix_module

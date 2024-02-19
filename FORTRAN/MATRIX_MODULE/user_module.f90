!
! The user module defines the matrix-vector product operation. It has
! been designed to be used by the matrix module that is used
! by the IDRS-routines. The matrix-module adds to the user defined
! matrix-vector product preconditioning operations, based on polynomial
! preconditioners
!
! The user module should always contain the following two functions:
!
!    ruser_mv:      real matrix-vector product
!    function call: xr = ruser_mv( A, yr )
!    data types:    xr, yr: real arrays of kind rp
!                   A: user_matrix (user defined data structure)
!
!    cuser_mv:      complex matrix-vector product
!    function call: xc = cuser_mv( A, yc )
!    data types:    xc, yc: complex arrays of kind cp
!                   A: user_matrix (user defined data structure)
!                  
! The functions ruser_mv and cuser_mv should be overloaded to "*" to define the 
! matrix-vector operation x = A*y, where x and y can be real or complex.
!
! This specific user_module defines these operations for matrices that
! are given in one of the follwing three formats:
!  1) Dense
!  2) Compressed Row Storage (CRS)
!  3) Coordinate (COO)
! The matrix-vector product can be performed in parallel. In this case
! the matrix is assumed to be partitioned row-wise.
!
! Next to the user_mv functions, this module contains the following routines:
!
!    A = dense_format( rmat, imat, cmat )
!    Purpose: store dense matrix in user-type matrix format
!    Input: 
!       rmat:    Real two-dimensional allocatable array, size local rows x unknowns. 
!                If allocated contains the values local rows of A on the processor.
!       imat:    Integer two-dimensional allocatable array, size local rows x unknowns. 
!                If allocated contains the values local rows of A on the processor.
!       cmat:    Complex two-dimensional allocatable array, size local rows x unknowns. 
!                If allocated contains the values local rows of A on the processor.
!    Note that only one of the three arrays rmat, imat, cmat should be allocated.
!    Output: A. The matrix stored in the data-type user_matrix. 
!
!    A = crs_format( nrows, row_ptr, col_ind, rval, ival, cval )
!    Purpose: store CRS-matrix in user-type matrix format
!    Input:
!       nrow: integer, number of (local) rows. Note that the global number of rows is nprocs*nrows
!       row_ptr: Integer array, length nrows+1. Gives the (local) row pointers in col_ind.
!                Note that row_ptr should be given in *local* numbering, corresponding to the rows in the processor.
!       col_ind: Integer array, length number of nonzeros on processor. Contains(global) column indices.
!                The complete row is stored on a processor, therefore col_ind should be given in global numbering.
!       rval:    Real allocatable array, length number of nonzeros on processor. If allocated contains the
!                values of the nonzero elements of A on the processor.
!       ival:    Integer allocatable array, length number of nonzeros on processor. If allocated contains the
!                values of the nonzero elements of A on the processor.
!       cval:    Complex allocatable array, length number of nonzeros on processor. If allocated contains the
!                values of the nonzero elements of A on the processor.
!    Note that only one of the three arrays rval, ival, cval should be allocated.
!    Output: A. The matrix stored in the data-type user_matrix.
!
!    A = coo_format( nrows, nnz, indx, jndx, rval, ival, cval )
!    Purpose: store COO-matrix in user-type matrix format
!    Input:
!       nrows:   Integer, number of (local) rows. Note that the global number of rows is nprocs*nrows
!       nnz:     Integer, local number of nonzeros
!       indx:    Integer array, length nnz. Gives the (local) row indices
!       jndx:    Integer array, length number of nonzeros on processor. Contains(global) column indices.
!       rval:    Real allocatable array, length number of nonzeros on processor. If allocated contains the
!                values of the nonzero elements of A on the processor.
!       ival:    Integer allocatable array, length number of nonzeros on processor. If allocated contains the
!                values of the nonzero elements of A on the processor.
!       cval:    Complex allocatable array, length number of nonzeros on processor. If allocated contains the
!                values of the nonzero elements of A on the processor.
!    Note that only one of the three arrays rval, ival, cval should be allocated.
!    Output: A. The matrix stored in the data-type user_matrix.
!
!    D = real_diagonal( A )
!    Purpose: get the main diagonal of a real matrix
!    Input:
!       A:       Matrix of type user_matrix
!    Output: D, real array. Main diagonal
!
!    D = complex_diagonal( A )
!    Purpose: get the main diagonal of a complex matrix
!    Input:
!       A:       Matrix of type user_matrix
!    Output: D, complex array. Main diagonal
!     
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2024 Martin van Gijzen
!

module user_module

   use precision_module

   implicit none

! Define the matrix type, as an example here for dense matrices

   type user_matrix
      integer                        :: matrix_format
      integer                        :: nrows            ! Number of (local) rows
      integer                        :: ncols            ! Number of columns

! Dense format
      real(kind=rp), allocatable     :: rmat(:,:)        ! Real dense matrix
      integer, allocatable           :: imat(:,:)        ! Integer dense matrix
      complex(kind=cp), allocatable  :: cmat(:,:)        ! Complex dense matrix

! Nonzero coefficients for sparse formats:
      real(kind=rp), allocatable     :: rval(:)          ! Real sparse matrix
      integer, allocatable           :: ival(:)          ! Integer sparse matrix
      complex(kind=cp), allocatable  :: cval(:)          ! Complex sparse matrix

! CSR format:
      integer, allocatable           :: row_ptr(:)
      integer, allocatable           :: col_ind(:)

! COO format:
      integer                        :: nnz                   ! Number of (local) nonzeros
      integer, allocatable           :: indx(:), jndx(:)

! For parallelisation:
      integer                        :: P = 1                 ! Number of processors
      integer, allocatable           :: phi(:)                ! Processor array  (1D)

   end type user_matrix

! Overload * to define the matrix-vector multiplication using the user-matrix type

   INTERFACE OPERATOR(*)
      module procedure ruser_mv, cuser_mv
   END INTERFACE

   contains

   function ruser_mv( A, v ) result(w)
!
! Real matvec in DENSE, CRS or COO format
! 
      type(user_matrix), intent(in)     :: A
      real(kind=rp), intent(in)         :: v(:)
      real(kind=rp)                     :: w(size(v))
      
      select case ( A%matrix_format )
         case(1) ! DENSE
            w = rdense_mv( A, v )
         case(2) ! CRS
            w = rcrs_mv( A, v )
         case(3) ! COO
            w = rcoo_mv( A, v )
         case default 
            stop 'No matrix format specified.'
      end select

   end function ruser_mv

   function cuser_mv( A, v ) result(w)
!
! Complex matvec in DENSE, CRS, or COO format
!
      type(user_matrix), intent(in)     :: A
      complex(kind=cp), intent(in)      :: v(:)
      complex(kind=cp)                  :: w(size(v))

      select case ( A%matrix_format )
         case(1) ! DENSE
            w = cdense_mv( A, v )
         case(2) ! CRS
            w = ccrs_mv(A, v )
         case(3) ! COO
            w = ccoo_mv( A, v )
         case default 
            stop 'Illegal value for matrix format'
      end select

   end function cuser_mv

   function ccrs_mv( A, v ) result(w)
!
! Complex matvec in dense or in CRS format
! 
      type(user_matrix), intent(in)     :: A
      complex(kind=cp), intent(in)      :: v(:)
      complex(kind=cp)                  :: w(size(v))
      integer                           :: i, j, k
      complex(kind=cp), allocatable,save:: co_v(:)[:]
      integer                           :: n_procs
      logical                           :: sequential
      
      n_procs = num_images()
! Parallel or sequential?
      sequential = ( n_procs == 1 )
      
! Initialize
      w = 0.

      if ( sequential ) then
         if ( allocated(A%cval) ) then
! CRS, matrix is complex
            do i = 1, A%nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  w(i) = w(i) + A%cval(j)*v(A%col_ind(j))
               end do
            end do
         elseif ( allocated(A%rval) ) then
! CRS, matrix is real
            do i = 1, A%nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  w(i) = w(i) + A%rval(j)*v(A%col_ind(j))
               end do
            end do
         elseif ( allocated(A%ival) ) then
! CRS, matrix is integer
            do i = 1, A%nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  w(i) = w(i) + A%ival(j)*v(A%col_ind(j))
               end do
            end do
         else
! CRS, matrix is binary, only pattern given
            do i = 1, A%nrows
               do j = A%row_ptr(j), A%row_ptr(j+1)-1
                  w(i) = w(i) + v(A%col_ind(j))
               end do
            end do
         end if
      else
         if ( .not. allocated(co_v) ) allocate(co_v(A%nrows)[*])
         co_v = v
         sync all
         if ( allocated(A%cval) ) then
! CRS, matrix is complex
            do i = 1, A%nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  k = A%phi(j)
                  w(i) = w(i) + A%cval(j)*co_v(A%col_ind(j))[k]
               end do
            end do
         elseif ( allocated(A%rval) ) then
! CRS, matrix is real
            do i = 1, A%nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  k = A%phi(j)
                  w(i) = w(i) + A%rval(j)*co_v(A%col_ind(j))[k]
               end do
            end do
         elseif ( allocated(A%ival) ) then
! CRS, matrix is real
            do i = 1, A%nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  k = A%phi(j)
                  w(i) = w(i) + A%ival(j)*co_v(A%col_ind(j))[k]
               end do
            end do
         else
! CRS, matrix is binary, only pattern given
            do i = 1, A%nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  k = A%phi(j)
                  w(i) = w(i) + co_v(A%col_ind(j))[k]
               end do
            end do
         end if
         sync all
      end if

   end function ccrs_mv

   function rcrs_mv( A, v ) result(w)
!
! Matvec in dense or in CRS format
! 
      type(user_matrix), intent(in)     :: A
      real(kind=rp), intent(in)         :: v(:)
      real(kind=rp)                     :: w(size(v))
      integer                           :: i, j, k, l
      real(kind=rp), allocatable,save   :: co_v(:)[:]
      integer                           :: n_procs
      logical                           :: sequential
      
      n_procs = num_images()
! Parallel or sequential?
      sequential = ( n_procs == 1 )

! Initialize
      w = 0.

      if ( sequential ) then
         if ( allocated(A%rval) ) then
! CRS, real matrix
            do i = 1, A%nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  w(i) = w(i) + A%rval(j)*v(A%col_ind(j))
               end do
            end do
         elseif ( allocated(A%ival) ) then
! CRS, integer matrix
            do i = 1, A%nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  w(i) = w(i) + A%ival(j)*v(A%col_ind(j))
               end do
            end do
         else
! CRS, binary matrix, only pattern given
            do i = 1, A%nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  w(i) = w(i) + v(A%col_ind(j))
               end do
            end do
         end if
      else
         if ( .not. allocated(co_v) ) allocate(co_v(A%nrows)[*])
         co_v = v
         sync all
         if ( allocated(A%rval) ) then
! CRS, matrix is real
            do i = 1, A%nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  k = A%phi(j)
                  w(i) = w(i) + A%rval(j)*co_v(A%col_ind(j))[k]
               end do
            end do
         elseif ( allocated(A%ival) ) then
! CRS, matrix is integer
            do i = 1, A%nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  k = A%phi(j)
                  w(i) = w(i) + A%ival(j)*co_v(A%col_ind(j))[k]
               end do
            end do
         else
! CRS, matrix is binary, only pattern given
            do i = 1, A%nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  k = A%phi(j)
                  w(i) = w(i) + co_v(A%col_ind(j))[k]
               end do
            end do
         end if
         sync all
      end if

   end function rcrs_mv

   function ccoo_mv( A, v ) result(w)
!
! Complex matvec in in COO format
! 
      type(user_matrix), intent(in)     :: A
      complex(kind=cp), intent(in)      :: v(:)
      complex(kind=cp)                  :: w(size(v))
      integer                           :: i, j, k
      complex(kind=cp), allocatable,save:: co_v(:)[:]
      logical                           :: sequential
      
! Parallel or sequential?
      sequential = ( A%P == 1 )
      
! Initialize
      w = 0.

      if ( sequential ) then
         if ( allocated(A%cval) ) then
! COO, matrix is complex
            do i = 1, A%nnz
               w(A%indx(i)) = w(A%indx(i)) + A%cval(i)*v(A%jndx(i))
            end do
         elseif ( allocated(A%rval) ) then
! COO, matrix is real
            do i = 1, A%nnz
               w(A%indx(i)) = w(A%indx(i)) + A%rval(i)*v(A%jndx(i))
            end do
         elseif ( allocated(A%ival) ) then
! COO, matrix is integer
            do i = 1, A%nnz
               w(A%indx(i)) = w(A%indx(i)) + A%ival(i)*v(A%jndx(i))
            end do
         else
! COO, matrix is binary, only pattern given
            do i = 1, A%nnz
               w(A%indx(i)) = w(A%indx(i)) + v(A%jndx(i))
            end do
         end if
      else
         if ( .not. allocated(co_v) ) allocate(co_v(A%nrows)[*])
         co_v = v
         sync all
         if ( allocated(A%cval) ) then
! COO, matrix is complex
            do i = 1, A%nnz
               k = A%phi(i)
               w(A%indx(i)) = w(A%indx(i)) + A%cval(i)*co_v(A%jndx(i))[k]
            end do
         elseif ( allocated(A%rval) ) then
! COO, matrix is real
            do i = 1, A%nnz
               k = A%phi(i)
               w(A%indx(i)) = w(A%indx(i)) + A%rval(i)*co_v(A%jndx(i))[k]
            end do
         elseif ( allocated(A%ival) ) then
! COO, matrix is real
            do i = 1, A%nnz
               k = A%phi(i)
               w(A%indx(i)) = w(A%indx(i)) + A%ival(i)*co_v(A%jndx(i))[k]
            end do
         else
! COO, matrix is binary, only pattern given
            do i = 1, A%nnz
               k = A%phi(i)
               w(A%indx(i)) = w(A%indx(i)) + co_v(A%jndx(i))[k]
            end do
         end if
         sync all
      end if

   end function ccoo_mv

   function rcoo_mv( A, v ) result(w)
!
! Matvec in COO format
! 
      type(user_matrix), intent(in)     :: A
      real(kind=rp), intent(in)         :: v(:)
      real(kind=rp)                     :: w(size(v))
      integer                           :: i, j, k, l
      real(kind=rp), allocatable,save   :: co_v(:)[:]
      logical                           :: sequential
      
! Parallel or sequential?
      sequential = ( A%P == 1 )

! Initialize
      w = 0.

      if ( sequential ) then
         if ( allocated(A%rval) ) then
! COO, real
            do i = 1, A%nnz
               w(A%indx(i)) = w(A%indx(i)) + A%rval(i)*v(A%jndx(i))
            end do
         elseif ( allocated(A%ival) ) then
! COO, integer
            do i = 1, A%nnz
               w(A%indx(i)) = w(A%indx(i)) + A%ival(i)*v(A%jndx(i))
            end do
         else
! COO, binary matrix, only pattern given
            do j = 1, A%nnz
               w(A%indx(i)) = w(A%indx(i)) + v(A%jndx(i))
            end do
         end if
      else
         if ( .not. allocated(co_v) ) allocate(co_v(A%nrows)[*])
         co_v = v
         sync all
         if ( allocated(A%rval) ) then
! COO, matrix is real
            do i = 1, A%nnz
               k = A%phi(i)
               w(A%indx(i)) = w(A%indx(i)) + A%rval(i)*co_v(A%jndx(i))[k]
            end do
         elseif ( allocated(A%ival) ) then
! COO, matrix is integer
            do i = 1, A%nnz
               k = A%phi(i)
               w(A%indx(i)) = w(A%indx(i)) + A%ival(i)*co_v(A%jndx(i))[k]
            end do
         else
! COO, matrix is binary, only pattern given
            do i = 1, A%nnz
               k = A%phi(i)
               w(A%indx(i)) = w(A%indx(i)) + co_v(A%jndx(i))[k]
            end do
         end if
         sync all
      end if

   end function rcoo_mv

   function rdense_mv( A, v ) result(w)
!
! Matvec in dense format
! 
      type(user_matrix), intent(in)     :: A
      real(kind=rp), intent(in)         :: v(:)
      real(kind=rp)                     :: w(size(v))
      real(kind=rp), allocatable        :: t(:)
      integer                           :: i
      real(kind=rp), allocatable,save   :: co_v(:)[:]
      integer                           :: n_procs, nrows
      logical                           :: sequential
      
! Parallel or sequential?
      n_procs = A%P
      sequential = ( n_procs == 1 )
      nrows = A%nrows

      if ( sequential ) then
         if ( allocated( A%rmat ) ) then
            w = matmul( A%rmat, v )
         else
            w = matmul( A%imat, v )
         end if
      else
         if ( .not. allocated(co_v) ) allocate(co_v(nrows)[*])
         co_v = v
         w = 0.
         sync all
         if ( allocated( A%rmat) ) then
            do i = 1, n_procs
               w = w + matmul( A%rmat(:,(i-1)*nrows+1:i*nrows), co_v(:)[i] )
            end do
         else
            do i = 1, n_procs
               w = w + matmul( A%imat(:,(i-1)*nrows+1:i*nrows), co_v(:)[i] )
            end do
         end if
         sync all
      end if
   end function rdense_mv

   function cdense_mv( A, v ) result(w)
!
! Matvec in dense format
! 
      type(user_matrix), intent(in)     :: A
      complex(kind=cp), intent(in)      :: v(:)
      complex(kind=cp)                  :: w(size(v))
      integer                           :: i
      complex(kind=cp), allocatable,save:: co_v(:)[:]
      integer                           :: n_procs, nrows
      logical                           :: sequential
      
! Parallel or sequential?
      n_procs = A%P
      sequential = ( n_procs == 1 )
      nrows = size(A%cmat,1)

      if ( sequential ) then
         if ( allocated(A%cmat) ) then
            w = matmul( A%cmat, v )
         elseif ( allocated( A%rmat) ) then
            w = matmul( A%rmat, v )
         elseif ( allocated( A%imat) ) then
            w = matmul( A%imat, v )
         end if
      else
         if ( .not. allocated(co_v) ) allocate(co_v(nrows)[*])
         co_v = v
         w = 0.
         sync all
         if ( allocated( A%cmat) ) then
            do i = 1, n_procs
               w = w + matmul( A%cmat(:,(i-1)*nrows+1:i*nrows), co_v(:)[i] )
            end do
         elseif ( allocated( A%rmat) ) then
            do i = 1, n_procs
               w = w + matmul( A%rmat(:,(i-1)*nrows+1:i*nrows), co_v(:)[i] )
            end do
         else
            do i = 1, n_procs
               w = w + matmul( A%imat(:,(i-1)*nrows+1:i*nrows), co_v(:)[i] )
            end do
         end if
         sync all
      end if
   end function cdense_mv

   function  crs_format( nrows, row_ptr, col_ind, rval, ival, cval ) result(A)
!
! Make crs matrix
!
      type(user_matrix)                      :: A
      integer, optional                      :: nrows 
      integer, optional                      :: row_ptr(:), col_ind(:)
      real(kind=rp), optional, intent(in)    :: rval(:)
      integer, optional, intent(in)          :: ival(:)
      complex(kind=cp), optional, intent(in) :: cval(:)
      integer                                :: i, j

      A%matrix_format = 2
      A%P = num_images()
      A%ncols = A%P*nrows

      allocate(A%phi(size(col_ind)),A%col_ind(size(col_ind)))
! Make row-cyclic distribution
      do i = 1, nrows
         do j = row_ptr(i), row_ptr(i+1)-1
            A%phi(j) = (col_ind(j)-1)/nrows + 1
            A%col_ind(j) = col_ind(j) - (A%phi(j)-1)*nrows
         end do
      end do
! Store everything
      A%nrows = nrows
      A%row_ptr = row_ptr
      if ( present(rval) ) A%rval = rval
      if ( present(ival) ) A%rval = ival
      if ( present(cval) ) A%cval = cval

   end function crs_format

   function  coo_format( nrows, nnz, indx, jndx, rval, ival, cval ) result(A)
!
! Make coo matrix
!
      type(user_matrix)                      :: A
      integer, optional                      :: nrows, nnz 
      integer, optional                      :: indx(:), jndx(:)
      real(kind=rp), optional, intent(in)    :: rval(:)
      integer, optional, intent(in)          :: ival(:)
      complex(kind=cp), optional, intent(in) :: cval(:)
      integer                                :: i, j

      A%matrix_format = 3
      A%P = num_images()
      A%nrows = nrows
      A%ncols = A%P*nrows
      A%nnz   = nnz

      allocate(A%phi(nnz),A%jndx(nnz))
! Make row-cyclic distribution
      do i = 1, nnz
         A%phi(i)  = (jndx(i)-1)/nrows + 1
         A%jndx(i) =  jndx(i) - (A%phi(i)-1)*nrows
      end do
! Store everything
      A%indx  = indx(1:nnz)
      if ( present(rval) ) A%rval = rval(1:nnz)
      if ( present(ival) ) A%rval = ival(1:nnz)
      if ( present(cval) ) A%cval = cval(1:nnz)

   end function coo_format

   function dense_format( rmat, imat, cmat ) result(A)
!
! Make dense matrix
!
      type(user_matrix)                      :: A
      real(kind=rp), optional                :: rmat(:,:)
      integer, optional                      :: imat(:,:)
      complex(kind=cp), optional             :: cmat(:,:)

      A%matrix_format = 1
      A%P = num_images()
      if ( present(rmat) ) then
! Dense real matrix
         A%nrows = size(rmat,1)
         A%ncols = size(rmat,2)
         A%rmat = rmat
      elseif ( present(imat) ) then
! Dense real matrix
         A%nrows = size(imat,1)
         A%ncols = size(imat,2)
         A%imat = imat
      elseif ( present(cmat) ) then
! Dense complex matrix
         A%nrows = size(cmat,1)
         A%ncols = size(cmat,2)
         A%cmat = cmat
      end if

   end function dense_format

   function real_diagonal( A ) result(D)

      type(user_matrix), intent(in)         :: A
      real(kind=rp)                         :: D(A%nrows)
      integer                               :: nrows, nnz, i, j, k, me

      me = this_image()
      nrows = A%nrows

      D = 0.
      if ( A%matrix_format == 1 ) then
! Dense
         if ( allocated(A%rmat) ) then
! Matrix is real
            do i = 1, nrows
               j = (me-1) * nrows+i
               D(i) = A%rmat(i,j)
            end do
         elseif ( allocated(A%imat) ) then
            do i = 1, nrows
               j = (me-1) * nrows+i
               D(i) = A%imat(i,j)
            end do
         end if
      elseif ( A%matrix_format == 2 ) then
! CRS
         if ( allocated(A%rval) ) then
! Matrix is real
            do i = 1, nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  if ( ( i == A%col_ind(j) ) .and. ( A%phi(j) == me ) ) D(i) = D(i) + A%rval(j)
               end do
            end do
         elseif ( allocated(A%ival) ) then
! Matrix is integer
            do i = 1, nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  if ( ( i == A%col_ind(j) ) .and. ( A%phi(j) == me ) ) D(i) = D(i) + A%ival(j)
               end do
            end do
         else
! Matrix is binary
            do i = 1, nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  if ( ( i == A%col_ind(j) ) .and. ( A%phi(j) == me ) ) D(i) = D(i) + 1.
               end do
            end do
         end if
      elseif ( A%matrix_format == 3 ) then
! COO
         nnz   = A%nnz
         if ( allocated(A%rval) ) then
! Matrix is real
            do k = 1,nnz
               i = A%indx(k)
               j = A%jndx(k)
               if ( ( i == j ) .and. ( A%phi(k) == me ) ) D(i) = D(i) + A%rval(k)
            end do
         elseif  ( allocated(A%ival) ) then
! Matrix is integer
            do k = 1,nnz
               i = A%indx(k)
               j = A%jndx(k)
               if ( ( i == j ) .and. ( A%phi(k) == me ) ) D(i) = D(i) + A%ival(k)
            end do
         else
! Matrix is binary
            do k = 1,nnz
               i = A%indx(k)
               j = A%jndx(k)
               if ( ( i == j ) .and. ( A%phi(k) == me ) ) D(i) = D(i) + 1.
            end do
         end if
      end if
      where ( D == 0 )
         D = 1._rp
      end where

   end function real_diagonal

   function complex_diagonal( A ) result(D)

      type(user_matrix), intent(in)         :: A
      complex(kind=cp)                      :: D(A%nrows)
      integer                               :: nrows, nnz, i, j, k, err, me

      me = this_image()                    ! This processor
      nrows = A%nrows

      D = 0.
      if ( A%matrix_format == 1 ) then
! Dense
         if ( allocated(A%rmat) ) then
! Matrix is real
            do i = 1, nrows
               j = (me-1) * nrows+i
               D(i) = A%rmat(i,j)
            end do
         elseif ( allocated(A%cmat) ) then
! Matrix is complex
            do i = 1, nrows
               j = (me-1) * nrows+i
               D(i) = A%cmat(i,j)
            end do
         elseif ( allocated(A%imat) ) then
! Matrix is complex
            do i = 1, nrows
               j = (me-1) * nrows+i
               D(i) = A%imat(i,j)
            end do
         end if
      elseif ( A%matrix_format == 2 ) then
! CRS
         if ( allocated(A%rval) ) then
! Matrix is real
            do i = 1, nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  if ( ( i == A%col_ind(j) ) .and. ( A%phi(j) == me ) ) D(i) = D(i) + A%rval(j)
               end do
            end do
         elseif ( allocated(A%cval) ) then
! Matrix is complex
            do i = 1, nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  if ( ( i == A%col_ind(j) ) .and. ( A%phi(j) == me )  ) D(i) = D(i) + A%cval(j)
               end do
            end do
         elseif ( allocated(A%ival) ) then
! Matrix is integer
            do i = 1, nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  if ( ( i == A%col_ind(j) ) .and. ( A%phi(j) == me )  ) D(i) = D(i) + A%ival(j)
               end do
            end do

         else
! Matrix is binary
            do i = 1, nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  if ( ( i == A%col_ind(j) ) .and. ( A%phi(j) == me ) ) D(i) = D(i) + 1.
               end do
            end do
         end if

      elseif ( A%matrix_format == 3 ) then
! COO
         nnz   = A%nnz
         if ( allocated(A%rval) ) then
! Matrix is real
            do k = 1,nnz
               i = A%indx(k)
               j = A%jndx(k)
               if ( ( i == j ) .and. ( A%phi(k) == me ) ) D(i) = D(i) + A%rval(k)
            end do
         elseif  ( allocated(A%cval) ) then
! Matrix is complex
            do k = 1,nnz
               i = A%indx(k)
               j = A%jndx(k)
               if ( ( i == j ) .and. ( A%phi(k) == me ) ) D(i) = D(i) + A%cval(k)
            end do
         elseif  ( allocated(A%ival) ) then
! Matrix is integer
            do k = 1,nnz
               i = A%indx(k)
               j = A%jndx(k)
               if ( ( i == j ) .and. ( A%phi(k) == me ) ) D(i) = D(i) + A%ival(k)
            end do
         else
! Matrix is binary
            do k = 1,nnz
               i = A%indx(k)
               j = A%jndx(k)
               if ( ( i == j ) .and. ( A%phi(k) == me ) ) D(i) = D(i) + 1.
            end do
         end if
      end if
      where ( D == 0. ) 
         D = 1._cp
      end where

   end function complex_diagonal

end module user_module

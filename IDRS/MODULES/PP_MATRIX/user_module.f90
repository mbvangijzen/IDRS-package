!
! The user module defines the matrix-vector product operation. It has
! been designed to be used by the matrix module that is used
! by the IDRS-routines. The matrix-module adds to the user defined
! matrix-vector product preconditioning operations, based on polynomial
! preconditioners
!
! The user module should always contain the following four functions:
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
!    ruser_precon:  real preconditioning operation
!    function call: xr = ruser_precon( yr, A )
!    data types:    xr, yr: real arrays of kind rp
!                   A: user_matrix (user defined data structure)
!    
!    cuser_precon:  complex preconditioning operation
!    function call: xc = cuser_precon( yc, A )
!    data types:    xc, yc: complex arrays of kind cp
!                   A: user_matrix (user defined data structure)
!                  
! The functions ruser_mv and cuser_mv should be overloaded to "*" to define the 
! matrix-vector operation x = A*y, where x and y can be real or complex.
!                  
! The functions ruser_precon and cuser_precon should be overloaded to "/" to define the 
! preconditioning operation x = y/A, where x and y can be real or complex.
!
! This specific user_module defines these operations for matrices that
! are given in one of the follwing four formats:
!  1) Dense
!  2) Compressed Row Storage (CRS)
!  3) Coordinate (COO)
!  4) Compressed Diagonal Storage (CDS)
! The matrix-vector product can be performed in parallel. In this case
! the matrix is assumed to be partitioned row-wise.
!
! Next to the user_mv functions, this module contains the following routines:
!
!    A = dense_format( rmat, imat, cmat )
!    Purpose: store dense matrix in user-type matrix format
!    Input: 
!       rmat:    Real two-dimensional allocatable array, size local rows x unknowns. 
!                If allocated contains the local rows of A on the processor.
!       imat:    Integer two-dimensional allocatable array, size local rows x unknowns. 
!                If allocated contains the local rows of A on the processor.
!       cmat:    Complex two-dimensional allocatable array, size local rows x unknowns. 
!                If allocated contains the local rows of A on the processor.
!    Note that only one of the three arrays rmat, imat, cmat should be allocated.
!    Output: A. The matrix stored in the data-type user_matrix. 
!
!    A = crs_format( nrows, ncols, row_ptr, col_ind, rval, ival, cval )
!    Purpose: store CRS-matrix in user-type matrix format
!    Input:
!       nrows: integer, number of (local) rows. Note that the global number of rows is nprocs*nrows
!       ncols: integer, number of columns. 
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
!    A = coo_format( nrows, ncols, nnz, indx, jndx, rval, ival, cval )
!    Purpose: store COO-matrix in user-type matrix format
!    Input:
!       nrows:   Integer, number of (local) rows. Note that the global number of rows is nprocs*nrows
!       ncols:   Integer, number of columns. Note that this is the number of equations in the system
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
!    A = cds_format( nrows, ncols, half_band, rband, iband, cband ) 
!    Purpose: store CDS-matrix in user-type matrix format
!    Input: 
!       nrows:     Integer, number of (local) rows. Note that the global number of rows is nprocs*nrows
!       ncols:     Integer, number of columns. Note that this is the number of equations in the system
!       half_band: Integer, half band_width
!       rband:     Real allocatable array, length nrows*2x(half_band+1).
!                  If allocated contains the diagonals of A within the band.
!       iband:     Integer allocatable array, length nrows*2x(half_band+1).
!                  If allocated contains the diagonals of A within the band.
!       cband:     Complex allocatable array, length nrows*2x(half_band+1).
!                  If allocated contains the diagonals of A within the band.
!    Note that only one of the three arrays rband, iband, cband should be allocated.
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
! Copyright:(c) 2025 Martin van Gijzen
!

module user_module

   use precision_module

   implicit none

! Define the matrix type

   type user_matrix
      integer                        :: matrix_format
      integer                        :: which_precon = 0
      integer                        :: nrows            ! Number of (local) rows
      integer                        :: ncols            ! Number columns

! DENSE format
      real(kind=rp), allocatable     :: rmat(:,:)        ! Real dense matrix
      integer, allocatable           :: imat(:,:)        ! Integer dense matrix
      complex(kind=cp), allocatable  :: cmat(:,:)        ! Complex dense matrix

! Nonzero coefficients for sparse formats CRS and COO:
      real(kind=rp), allocatable     :: rval(:)          ! Real sparse matrix
      integer, allocatable           :: ival(:)          ! Integer sparse matrix
      complex(kind=cp), allocatable  :: cval(:)          ! Complex sparse matrix

! CRS format:
      integer, allocatable           :: row_ptr(:)
      integer, allocatable           :: col_ind(:)

! COO format:
      integer                        :: nnz               ! Number of (local) nonzeros
      integer, allocatable           :: indx(:), jndx(:)

! CDS format:
      real(kind=rp), allocatable     :: rband(:,:)        ! Real band matrix
      integer, allocatable           :: iband(:,:)        ! Integer band matrix
      complex(kind=cp), allocatable  :: cband(:,:)        ! Complex band matrix
      integer                        :: half_band = 0     ! half bandwidth

! For parallelisation:
      integer, allocatable           :: phi(:)            ! Processor array  (1D)

   end type user_matrix

! Overload * to define the matrix-vector multiplication using the user-matrix type

   INTERFACE OPERATOR(*)
      module procedure ruser_mv, cuser_mv
   END INTERFACE

   INTERFACE OPERATOR(/)
      module procedure ruser_precon, cuser_precon
   END INTERFACE

   interface LU_decomp
      module procedure CLU_decomp, RLU_decomp
   end interface

   interface LU_solve
      module procedure CLU_solve, RLU_solve
   end interface

   interface band_LU_decomp
      module procedure rband_LU_decomp, cband_LU_decomp
   end interface

   interface band_LU_solve
      module procedure cband_LU_solve, rband_LU_solve
   end interface

   contains

   function ruser_mv( A, v ) result(w)
!
! Real matvec in DENSE, CRS or COO format
! 
      type(user_matrix), intent(in)     :: A
      real(kind=rp), intent(in)         :: v(:)
      real(kind=rp)                     :: w(A%nrows)
      
      select case ( A%matrix_format )
         case(1) ! DENSE
            w = rdense_mv( A, v )
         case(2) ! CRS
            w = rcrs_mv( A, v )
         case(3) ! COO
            w = rcoo_mv( A, v )
         case(4) ! CDS
            w = rcds_mv( A, v )
         case default 
            stop 'No matrix format specified.'
      end select

   end function ruser_mv

   function cuser_mv( A, v ) result(w)
!
! Complex matvec in DENSE, CRS, COO, or CDS format
!
      type(user_matrix), intent(in)     :: A
      complex(kind=cp), intent(in)      :: v(:)
      complex(kind=cp)                  :: w(A%nrows)

      select case ( A%matrix_format )
         case(1) ! DENSE
            w = cdense_mv( A, v )
         case(2) ! CRS
            w = ccrs_mv(A, v )
         case(3) ! COO
            w = ccoo_mv( A, v )
         case(4) ! CDS
            w = ccds_mv( A, v )
         case default 
            stop 'Illegal value for matrix format'
      end select

   end function cuser_mv

   function ruser_precon( v, A ) result(w)
!
! Real preconditioner: diagonal scaling, dense LU or band LU
! 
      type(user_matrix), intent(in)     :: A
      real(kind=rp), intent(in)         :: v(:)
      real(kind=rp)                     :: w(size(v))
      
      select case ( A%which_precon )
         case(1) ! Diagonal scaling
            w = v/A%rband(:,0)
         case(2) ! Dense LU
            w = LU_solve( v, A%rmat )
         case(3) ! Band LU
            w = band_LU_solve( v, A%rband, A%half_band, A%ncols )
         case default 
            w = v
      end select

   end function ruser_precon

   function cuser_precon( v, A ) result(w)
!
! Complex preconditioner: diagonal scaling, dense LU or band LU
!
      type(user_matrix), intent(in)     :: A
      complex(kind=cp), intent(in)      :: v(:)
      complex(kind=cp)                  :: w(size(v))

      select case ( A%which_precon )
         case(1) ! Diagonal scaling
            w = v/A%cband(:,0)
         case(2) ! Dense LU
            w = LU_solve( v, A%cmat )
         case(3) ! Band LU
            w = band_LU_solve( v, A%cband, A%half_band, A%ncols )
         case default 
            w = v
      end select

   end function cuser_precon

   function ccrs_mv( A, v ) result(w)
!
! Complex matvec in CRS format
! 
      type(user_matrix), intent(in)     :: A
      complex(kind=cp), intent(in)      :: v(:)
      complex(kind=cp)                  :: w(A%nrows)
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
            do concurrent (i = 1:A%nrows)
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  w(i) = w(i) + A%cval(j)*v(A%col_ind(j))
               end do
            end do
         elseif ( allocated(A%rval) ) then
! CRS, matrix is real
            do concurrent (i = 1:A%nrows)
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  w(i) = w(i) + A%rval(j)*v(A%col_ind(j))
               end do
            end do
         elseif ( allocated(A%ival) ) then
! CRS, matrix is integer
            do concurrent (i = 1:A%nrows)
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  w(i) = w(i) + A%ival(j)*v(A%col_ind(j))
               end do
            end do
         else
! CRS, matrix is binary, only pattern given
            do concurrent (i = 1:A%nrows)
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
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
! Matvec in CRS format
! 
      type(user_matrix), intent(in)     :: A
      real(kind=rp), intent(in)         :: v(:)
      real(kind=rp)                     :: w(A%nrows)
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
            do concurrent (i = 1:A%nrows)
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  w(i) = w(i) + A%rval(j)*v(A%col_ind(j))
               end do
            end do
         elseif ( allocated(A%ival) ) then
! CRS, integer matrix
            do concurrent (i = 1:A%nrows)
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  w(i) = w(i) + A%ival(j)*v(A%col_ind(j))
               end do
            end do
         else
! CRS, binary matrix, only pattern given
            do concurrent (i = 1:A%nrows)
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
      complex(kind=cp)                  :: w(A%nrows)
      integer                           :: i, j, k
      complex(kind=cp), allocatable,save:: co_v(:)[:]
      logical                           :: sequential
      
! Parallel or sequential?
      sequential = ( num_images() == 1 )
      
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
      real(kind=rp)                     :: w(A%nrows)
      integer                           :: i, j, k, l
      real(kind=rp), allocatable,save   :: co_v(:)[:]
      logical                           :: sequential
      
! Parallel or sequential?
      sequential = ( num_images() == 1 )

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
            do i = 1, A%nnz
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
      real(kind=rp)                     :: w(A%nrows)
      real(kind=rp), allocatable        :: t(:)
      integer                           :: i
      real(kind=rp), allocatable,save   :: co_v(:)[:]
      integer                           :: n_procs, nrows
      logical                           :: sequential
      
! Parallel or sequential?
      n_procs = num_images()
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
      complex(kind=cp)                  :: w(A%nrows)
      integer                           :: i
      complex(kind=cp), allocatable,save:: co_v(:)[:]
      integer                           :: n_procs, nrows
      logical                           :: sequential
      
! Parallel or sequential?
      n_procs = num_images()
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

   function rcds_mv( A, v ) result(w)
!
! Matvec in cds format
! 
      type(user_matrix), intent(in)     :: A
      real(kind=rp), intent(in)         :: v(:)
      real(kind=rp)                     :: w(size(v))
      integer                           :: dia, i, j
      real(kind=rp), allocatable,save   :: co_v(:)[:]
      integer                           :: n_procs, my_proc, nrows, half_band, i_proc, rest
      integer                           :: low, up
      logical                           :: sequential

! Parallel or sequential?
      my_proc = this_image()
      n_procs = num_images()
      sequential = ( n_procs == 1 )
      nrows = A%nrows
      half_band = A%half_band

      if ( .not. allocated(co_v) ) allocate(co_v(-half_band:nrows-1+half_band)[*])
      co_v = 0._rp
      co_v(0:nrows-1) = v

      if ( .not. sequential ) then
! Get the data from the neighbours
         sync all
         do i_proc = 1, half_band/nrows
            if ( (my_proc-i_proc) >= 1       ) &
               co_v(-i_proc*nrows:(-i_proc+1)*nrows-1) = co_v(0:nrows-1)[my_proc-i_proc]
            if ( (my_proc+i_proc) <= n_procs ) &
               co_v( i_proc*nrows:( i_proc+1)*nrows-1) = co_v(0:nrows-1)[my_proc+i_proc]
         end do
         rest = mod(half_band,nrows)
         if ( rest > 0 ) then
            i_proc = half_band/nrows+1
            if ( (my_proc-i_proc) >= 1       ) then
               low = -(i_proc-1)*nrows-rest
               up  = -(i_proc-1)*nrows-1
               co_v(low:up) = co_v(nrows-rest:nrows-1)[my_proc-i_proc]
            end if
            if ( (my_proc+i_proc) <= n_procs ) then
               low = i_proc*nrows
               up  = i_proc*nrows+rest-1
               co_v(low:up) = co_v(0:rest-1)[my_proc+i_proc]
            end if
         end if
         sync all
      end if

! Perform the multiplication
      w = 0._rp
      if ( allocated( A%rband ) ) then
         do dia = -half_band,half_band
            do concurrent ( i = 1:nrows )
               w(i) = w(i) + A%rband(i,dia)*co_v(i-1+dia)
            end do
         end do
      elseif ( allocated( A%iband ) ) then
         do dia = -half_band,half_band
            do concurrent ( i = 1:nrows )
               w(i) = w(i) + A%iband(i,dia)*co_v(i-1+dia)
            end do
         end do
      else
         do dia = -half_band,half_band
            do concurrent ( i = 1:nrows )
               w(i) = w(i) + co_v(i-1+dia)
            end do
         end do
      end if
   end function rcds_mv

   function ccds_mv( A, v ) result(w)
!
! Matvec in cds format
! 
      type(user_matrix), intent(in)       :: A
      complex(kind=cp), intent(in)        :: v(:)
      complex(kind=cp)                    :: w(A%nrows)
      integer                             :: dia, i, j
      complex(kind=cp), allocatable,save  :: co_v(:)[:]
      integer                             :: n_procs, my_proc, nrows, half_band, i_proc, rest
      integer                             :: low, up
      logical                             :: sequential
      
! Parallel or sequential?
      my_proc = this_image()
      n_procs = num_images()
      sequential = ( n_procs == 1 )
      nrows = A%nrows
      half_band = A%half_band

      if ( .not. allocated(co_v) ) allocate(co_v(-half_band:nrows-1+half_band)[*])
      co_v = 0._cp
      co_v(0:nrows-1) = v

      if ( .not. sequential ) then
! Get the data from the neighbours
         sync all
         do i_proc = 1, half_band/nrows
            if ( (my_proc-i_proc) >= 1       ) &
               co_v(-i_proc*nrows:(-i_proc+1)*nrows-1) = co_v(0:nrows-1)[my_proc-i_proc]
            if ( (my_proc+i_proc) <= n_procs ) &
               co_v( i_proc*nrows:( i_proc+1)*nrows-1) = co_v(0:nrows-1)[my_proc+i_proc]
         end do
         rest = mod(half_band,nrows)
         if ( rest > 0 ) then
            i_proc = half_band/nrows+1
            if ( (my_proc-i_proc) >= 1       ) then
               low = -(i_proc-1)*nrows-rest
               up  = -(i_proc-1)*nrows-1
               co_v(low:up) = co_v(nrows-rest:nrows-1)[my_proc-i_proc]
            end if
            if ( (my_proc+i_proc) <= n_procs ) then
               low = i_proc*nrows
               up  = i_proc*nrows+rest-1
               co_v(low:up) = co_v(0:rest-1)[my_proc+i_proc]
            end if
         end if
         sync all
      end if

! Perform the multiplication
      w = 0._rp
      if ( allocated( A%rband ) ) then
         do dia = -half_band,half_band
            do concurrent ( i = 1:nrows )
               w(i) = w(i) + A%rband(i,dia)*co_v(i-1+dia)
            end do
         end do
      elseif ( allocated( A%iband ) ) then
         do dia = -half_band,half_band
            do concurrent ( i = 1:nrows )
               w(i) = w(i) + A%iband(i,dia)*co_v(i-1+dia)
            end do
         end do
      elseif ( allocated( A%cband ) ) then
         do dia = -half_band,half_band
            do concurrent ( i = 1:nrows )
               w(i) = w(i) + A%cband(i,dia)*co_v(i-1+dia)
            end do
         end do
      else
         do dia = -half_band,half_band
            do concurrent ( i = 1:nrows )
               w(i) = w(i) + co_v(i-1+dia)
            end do
         end do
      end if

   end function ccds_mv

   function crs_format( nrows, ncols, row_ptr, col_ind, rval, ival, cval ) result(A)
!
! Make crs matrix
!
      type(user_matrix)                      :: A
      integer                                :: nrows , ncols, loc_ncols
      integer                                :: row_ptr(:), col_ind(:)
      real(kind=rp), optional, intent(in)    :: rval(:)
      integer, optional, intent(in)          :: ival(:)
      complex(kind=cp), optional, intent(in) :: cval(:)
      integer                                :: i, j
      integer                                :: n_procs

      A%matrix_format = 2
      n_procs = num_images()

      allocate(A%phi(size(col_ind)),A%col_ind(size(col_ind)))
! Make row distribution
      loc_ncols = ceiling(real(ncols)/real(n_procs))
      do i = 1, nrows
         do j = row_ptr(i), row_ptr(i+1)-1
            A%phi(j) = (col_ind(j)-1)/loc_ncols + 1
            A%col_ind(j) = col_ind(j) - (A%phi(j)-1)*loc_ncols
         end do
      end do

! Store everything
      A%nrows = nrows
      A%ncols = ncols
      A%row_ptr = row_ptr
      if ( present(rval) ) A%rval = rval
      if ( present(ival) ) A%rval = ival
      if ( present(cval) ) A%cval = cval

   end function crs_format

   function coo_format( nrows, ncols, nnz, indx, jndx, rval, ival, cval ) result(A)
!
! Make coo matrix
!
      type(user_matrix)                      :: A
      integer                                :: nrows, ncols, nnz , loc_ncols
      integer                                :: indx(:), jndx(:)
      real(kind=rp), optional, intent(in)    :: rval(:)
      integer, optional, intent(in)          :: ival(:)
      complex(kind=cp), optional, intent(in) :: cval(:)
      integer                                :: i, j
      integer                                :: n_procs

      A%matrix_format = 3
      n_procs = num_images()

      allocate(A%phi(nnz),A%jndx(nnz))
! Make row distribution
      loc_ncols = ceiling(real(ncols)/real(n_procs))
      do i = 1, nnz
         A%phi(i)  = (jndx(i)-1)/loc_ncols + 1
         A%jndx(i) =  jndx(i) - (A%phi(i)-1)*loc_ncols
      end do

! Store everything
      A%nrows = nrows
      A%ncols = ncols
      A%nnz   = nnz
      A%indx  = indx(1:nnz)
      if ( present(rval) ) A%rval = rval(1:nnz)
      if ( present(ival) ) A%ival = ival(1:nnz)
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
      integer                                :: n_procs

      A%matrix_format = 1
      n_procs = num_images()

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

   function cds_format( nrows, ncols, half_band, rband, iband, cband ) result(A)
!
! Make dense matrix
!
      type(user_matrix)                      :: A
      integer                                :: nrows, ncols, half_band
      real(kind=rp), optional                :: rband(nrows,-half_band:half_band)
      integer, optional                      :: iband(nrows,-half_band:half_band)
      complex(kind=cp), optional             :: cband(nrows,-half_band:half_band)

      A%matrix_format = 4
      A%nrows = nrows
      A%ncols = ncols
      A%half_band = half_band

      if ( present(rband) ) then
! Real band matrix
         A%rband = rband
      elseif ( present(iband) ) then
! Integer band matrix
         A%iband = iband
      elseif ( present(cband) ) then
! Complex band matrix
         A%cband = cband
      end if

   end function cds_format

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
      elseif ( A%matrix_format == 4 ) then
! CDS
         if ( allocated(A%rband) ) then
            D = A%rband(:,0)
         elseif ( allocated(A%iband) ) then
            D = A%iband(:,0)
         else
            D = 1._rp
         end if
      end if
      where ( D == 0 )
         D = 1._rp
      end where

   end function real_diagonal

   function complex_diagonal( A ) result(D)

      type(user_matrix), intent(in)         :: A
      complex(kind=cp)                      :: D(A%nrows)
      integer                               :: nrows, nnz, i, j, k, me

      me = this_image()                    ! This processor
      nrows = A%nrows

      D = 0._cp
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
! Matrix is integer
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
      elseif ( A%matrix_format == 4 ) then
! CDS
         if ( allocated(A%rband) ) then
            D = A%rband(:,0)
         elseif ( allocated(A%iband) ) then
            D = A%iband(:,0)
         elseif ( allocated(A%cband) ) then
            D = A%cband(:,0)
         else
            D = 1._cp
         end if
      end if
      where ( D == 0._cp ) 
         D = 1._cp
      end where

   end function complex_diagonal

   function real_dense( A ) result(LU)

      type(user_matrix), intent(in)         :: A
      real(kind=rp)                         :: LU(A%nrows,A%ncols)
      integer                               :: nrows, nnz, i, j, k, me

      me = this_image()                    ! This processor
      nrows = A%nrows

      LU = 0._rp
      if ( A%matrix_format == 1 ) then
! Dense
         if ( allocated(A%rmat) ) then
            LU = A%rmat
         elseif ( allocated(A%imat) ) then
            LU = A%imat
         end if
      elseif ( A%matrix_format == 2 ) then
! CRS
         if ( allocated(A%rval) ) then
! Matrix is real
            do i = 1, nrows
               do k = A%row_ptr(i), A%row_ptr(i+1)-1
                  j = A%col_ind(k)+(A%phi(k)-1)*nrows
                  LU(i,j) = LU(i,j) + A%rval(k)
               end do
            end do
         elseif ( allocated(A%ival) ) then
! Matrix is integer
            do i = 1, nrows
               do k = A%row_ptr(i), A%row_ptr(i+1)-1
                  j = A%col_ind(k)+(A%phi(k)-1)*nrows
                  LU(i,j) = LU(i,j) + A%ival(k)
               end do
            end do
         else
! Matrix is binary
            do i = 1, nrows
               do k = A%row_ptr(i), A%row_ptr(i+1)-1
                  j = A%col_ind(k)+(A%phi(k)-1)*nrows
                  LU(i,j) = LU(i,j) + 1._rp
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
               j = A%jndx(k)+(A%phi(k)-1)*nrows
               LU(i,j) = LU(i,j) + A%rval(k)
            end do
         elseif  ( allocated(A%ival) ) then
! Matrix is integer
            do k = 1,nnz
               i = A%indx(k)
               j = A%jndx(k)+(A%phi(k)-1)*nrows
               LU(i,j) = LU(i,j) + A%ival(k)
            end do
         else
! Matrix is binary
            do k = 1,nnz
               i = A%indx(k)
               j = A%jndx(k)+(A%phi(k)-1)*nrows
               LU(i,j) = LU(i,j) + 1._rp
            end do
         end if
      end if

   end function real_dense

   function complex_dense( A ) result(LU)

      type(user_matrix), intent(in)         :: A
      complex(kind=cp)                      :: LU(A%nrows,A%ncols)
      integer                               :: nrows, nnz, i, j, k

      nrows = A%nrows

      LU = 0._cp
      if ( A%matrix_format == 1 ) then
! Dense
         if ( allocated(A%rmat) ) then
            LU = cmplx(A%rmat,kind=cp)
         elseif ( allocated(A%cmat) ) then
            LU = A%cmat
         elseif ( allocated(A%imat) ) then
            LU = cmplx(real(A%imat,kind=cp))
         end if
      elseif ( A%matrix_format == 2 ) then
! CRS
         if ( allocated(A%rval) ) then
! Matrix is real
            do i = 1, nrows
               do k = A%row_ptr(i), A%row_ptr(i+1)-1
                  j = A%col_ind(k)+(A%phi(k)-1)*nrows
                  LU(i,j) = LU(i,j) + A%rval(k)
               end do
            end do
         elseif ( allocated(A%cval) ) then
! Matrix is complex
            do i = 1, nrows
               do k = A%row_ptr(i), A%row_ptr(i+1)-1
                  j = A%col_ind(k)+(A%phi(k)-1)*nrows
                  LU(i,j) = LU(i,j) + A%cval(k)
               end do
            end do
         elseif ( allocated(A%ival) ) then
! Matrix is integer
            do i = 1, nrows
               do k = A%row_ptr(i), A%row_ptr(i+1)-1
                  j = A%col_ind(k)+(A%phi(k)-1)*nrows
                  LU(i,j) = LU(i,j) + A%ival(k)
               end do
            end do
         else
! Matrix is binary
            do i = 1, nrows
               do k = A%row_ptr(i), A%row_ptr(i+1)-1
                  j = A%col_ind(k)+(A%phi(k)-1)*nrows
                  LU(i,j) = LU(i,j) + 1._cp
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
               j = A%jndx(k)+(A%phi(k)-1)*nrows
               LU(i,j) = LU(i,j) + A%rval(k)
            end do
         elseif  ( allocated(A%cval) ) then
! Matrix is complex
            do k = 1,nnz
               i = A%indx(k)
               j = A%jndx(k)+(A%phi(k)-1)*nrows
               LU(i,j) = LU(i,j) + A%cval(k)
            end do
         elseif  ( allocated(A%ival) ) then
! Matrix is integer
            do k = 1,nnz
               i = A%indx(k)
               j = A%jndx(k)+(A%phi(k)-1)*nrows
               LU(i,j) = LU(i,j) + A%ival(k)
            end do
         else
! Matrix is binary
            do k = 1,nnz
               i = A%indx(k)
               j = A%jndx(k)+(A%phi(k)-1)*nrows
               LU(i,j) = LU(i,j) + 1._cp
            end do
         end if
      end if

   end function complex_dense

   function REAL_PRECON( which_precon, K, M, C, shift ) result(P)
!
! Computer preconditioning matrix:
!    which_precon = 0 : no preconditioner
!    which_precon = 1 : diagonal preconditioner
!    which_precon = 2 : (dense) shift-and-invert preconditioner
!    which_precon = 3 : (band) shift-and-invert preconditioner
!

      integer, intent(in)                              :: which_precon
      type(user_matrix), intent(in)                    :: K
      type(user_matrix), intent(in), optional          :: M
      type(user_matrix), intent(in), optional          :: C
      type(user_matrix)                                :: P
      real(kind=rp), optional                          :: shift

      integer                                          :: i, j, me, half_band, my_proc

      P%matrix_format = 0
      P%which_precon = which_precon
      P%nrows = K%nrows
      P%ncols = K%ncols

      my_proc = this_image()

      if ( which_precon == 1 ) then
         if ( my_proc == 1 ) write(*,'(a)') 'Preconditioner is diagonal scaling.'
!
! Diagonal scaling
!
         P%matrix_format = 4
         P%half_band = 0
         if ( allocated(P%rband) ) deallocate(P%rband)
         allocate(P%rband(P%nrows,-P%half_band:P%half_band))
         if ( present(shift) ) then
            if ( present(M) .and. present(C) ) then
               P%rband(:,0) = real_diagonal(K) &
                            + shift*real_diagonal(C) &
                            + (shift**2)*real_diagonal(M)
            elseif ( present(C) ) then
               P%rband(:,0) = real_diagonal(K) &
                            + shift*real_diagonal(C)
! Mass matrix is identity, add squared shift to main diagonal:
               P%rband(:,0) = P%rband(:,0) + shift**2
            elseif ( present(M) ) then
               P%rband(:,0) = real_diagonal(K) &
                            - shift*real_diagonal(M)
            else
               P%rband(:,0) = real_diagonal(K)
! Mass matrix is identity, substract shift from main diagonal:
               P%rband(:,0) = P%rband(:,0) - shift
            end if
         else
            P%rband(:,0) = real_diagonal( K )
         end if

      elseif ( which_precon == 2 ) then
         if ( my_proc == 1 ) write(*,'(a)') 'Preconditioner is dense LU.'
!
! Dense LU decomposition
!
         P%matrix_format = 1
         if ( allocated(P%rmat) ) deallocate(P%rmat)
         allocate(P%rmat(P%nrows,P%ncols))

         me = this_image()                    ! This processor

         if ( present(shift) ) then
            if ( present(M) .and. present(C) ) then
               P%rmat = real_dense(K) &
                  + shift*real_dense(C) &
                  + (shift**2)*real_dense(M)
            elseif ( present(C) ) then
               P%rmat = real_dense(K) &
                  + shift*real_dense(C)
! Mass matrix is identity, add squared shift to main diagonal:
               do i = 1, P%nrows
                  j = (me-1)*P%nrows+i
                  P%rmat(i,j) = P%rmat(i,j) + shift**2
               end do
            elseif ( present(M) ) then
               P%rmat = real_dense(K) &
                  - shift*real_dense(M)
            else
               P%rmat = real_dense(K)
! Mass matrix is identity, substract shift from main diagonal:
               do i = 1, P%nrows
                  j = (me-1)*P%nrows+i
                  P%rmat(i,j) = P%rmat(i,j) - shift
               end do
            end if
         else
            P%rmat = real_dense( K )
         end if
         call LU_decomp( P%rmat )

      elseif ( which_precon == 3 ) then
         if ( my_proc == 1 ) write(*,'(a)') 'Preconditioner is band LU.'
!
! Band LU decomposition
!
         P%matrix_format = 4
! Determine bandwidths of the matrices:
         half_band = get_band( K )
         if ( present(C) ) half_band = max( half_band, get_band( C ) )
         if ( present(M) ) half_band = max( half_band, get_band( M ) )

         if ( allocated(P%rband) ) deallocate(P%rband)
         allocate(P%rband(K%nrows,-half_band:half_band))
         P%half_band = half_band

         if ( present(shift) ) then
            if ( present(M) .and. present(C) ) then
               P%rband = real_band(K,half_band) &
                  + shift*real_band(C,half_band) &
                  + (shift**2)*real_band(M,half_band)
            elseif ( present(C) ) then
               P%rband = real_band(K,half_band) &
                  + shift*real_band(C,half_band)
! Mass matrix is identity, add squared shift to main diagonal:
               P%rband(:,0) = P%rband(:,0) + shift**2
            elseif ( present(M) ) then
               P%rband = real_band(K,half_band) &
                  - shift*real_band(M,half_band)
            else
               P%rband = real_band(K,half_band)
! Mass matrix is identity, substract shift from main diagonal:
               P%rband(:,0) = P%rband(:,0) - shift
            end if
         else
            P%rband = real_band( K,half_band )
         end if
         call band_LU_decomp( P%rband, P%ncols, half_band )

      end if

   end function REAL_PRECON

   function COMPLEX_PRECON( which_precon, K, M, C, shift ) result(P)
!
! Compute complex preconditioning matrix:
!    which_precon = 0 : no preconditioner
!    which_precon = 1 : diagonal preconditioner
!    which_precon = 2 : (dense) shift-and-invert preconditioner
!    which_precon = 3 : (band) shift-and-invert preconditioner
!

      integer, intent(in)                              :: which_precon
      type(user_matrix), intent(in)                    :: K
      type(user_matrix), intent(in), optional          :: M
      type(user_matrix), intent(in), optional          :: C
      type(user_matrix)                                :: P
      complex(kind=cp), optional                       :: shift

      integer                                          :: i, j, me, half_band, my_proc

      P%matrix_format = 0
      P%which_precon = which_precon
      P%nrows = K%nrows
      P%ncols = K%ncols
      my_proc = this_image()

      if ( which_precon == 1 ) then
         if ( my_proc == 1 ) write(*,'(a)') 'Preconditioner is diagonal scaling.'
!
! Diagonal scaling
!
         P%matrix_format = 4
         P%half_band = 0
         if ( allocated(P%cband) ) deallocate(P%cband)
         allocate(P%cband(P%nrows,-P%half_band:P%half_band))
         if ( present(shift) ) then
            if ( present(M) .and. present(C) ) then
               P%cband(:,0) = complex_diagonal(K) &
                            + shift*complex_diagonal(C) &
                            + (shift**2)*complex_diagonal(M)
            elseif ( present(C) ) then
               P%cband(:,0) = complex_diagonal(K) &
                            + shift*complex_diagonal(C)
! Mass matrix is identity, add squared shift to main diagonal:
               P%cband(:,0) = P%cband(:,0) + shift**2
            elseif ( present(M) ) then
               P%cband(:,0) = complex_diagonal(K) &
                            - shift*complex_diagonal(M)
            else
               P%cband(:,0) = complex_diagonal(K)
! Mass matrix is identity, substract shift from main diagonal:
               P%cband(:,0) = P%cband(:,0) - shift
            end if
         else
            P%cband(:,0) = complex_diagonal( K )
         end if

      elseif ( which_precon == 2 ) then
         if ( my_proc == 1 ) write(*,'(a)') 'Preconditioner is dense LU.'
!
! Dense LU decomposition
!
         P%matrix_format = 1
         if ( allocated(P%cmat) ) deallocate(P%cmat)
         allocate(P%cmat(P%nrows,P%ncols))

         me = this_image()                    ! This processor

         if ( present(shift) ) then
            if ( present(M) .and. present(C) ) then
               P%cmat = complex_dense(K) &
                  + shift*complex_dense(C) &
                  + (shift**2)*complex_dense(M)
            elseif ( present(C) ) then
               P%cmat = complex_dense(K) &
                  + shift*complex_dense(C)
! Mass matrix is identity, add squared shift to main diagonal:
               do i = 1, P%nrows
                  j = (me-1)*P%nrows+i
                  P%cmat(i,j) = P%cmat(i,j) + shift**2
               end do
            elseif ( present(M) ) then
               P%cmat = complex_dense(K) &
                  - shift*complex_dense(M)
            else
               P%cmat = complex_dense(K)
! Mass matrix is identity, substract shift from main diagonal:
               do i = 1, P%nrows
                  j = (me-1)*P%nrows+i
                  P%cmat(i,j) = P%cmat(i,j) - shift
               end do
            end if
         else
            P%cmat = complex_dense( K )
         end if
         call LU_decomp( P%cmat )

      elseif ( which_precon == 3 ) then
         if ( my_proc == 1 ) write(*,'(a)') 'Preconditioner is band LU.'
!
! Band LU decomposition
!
         P%matrix_format = 4
! Determine bandwidths of the matrices:
         half_band = get_band( K )
         if ( present(C) ) half_band = max( half_band, get_band( C ) )
         if ( present(M) ) half_band = max( half_band, get_band( M ) )

         if ( allocated(P%cband) ) deallocate(P%cband)
         allocate(P%cband(K%nrows,-half_band:half_band))
         P%half_band = half_band

         me = this_image()                    ! This processor

         if ( present(shift) ) then
            if ( present(M) .and. present(C) ) then
               P%cband = complex_band(K,half_band) &
                  + shift*complex_band(C,half_band) &
                  + (shift**2)*complex_band(M,half_band)
            elseif ( present(C) ) then
               P%cband = complex_band(K,half_band) &
                  + shift*complex_band(C,half_band)
! Mass matrix is identity, add squared shift to main diagonal:
               P%cband(:,0) = P%cband(:,0) + shift**2
            elseif ( present(M) ) then
               P%cband = complex_band(K,half_band) &
                  - shift*complex_band(M,half_band)
            else
               P%cband = complex_band(K,half_band)
! Mass matrix is identity, substract shift from main diagonal:
               P%cband(:,0) = P%cband(:,0) - shift
            end if
         else
            P%cband = complex_band( K,half_band )
         end if
         call band_LU_decomp( P%cband, P%ncols, half_band )

      end if

   end function COMPLEX_PRECON

   subroutine RLU_decomp( LU )
      
! Make in place LU decomposition of the dense matrix LU
! This is the real version.
!  
      
      real(kind=rp), intent(inout)       :: LU(:,:)
      real(kind=rp)                      :: pivot_row(size(LU,2))
      integer                            :: i, k, neq, nrows, offset
      integer                            :: nprocs, my_proc, owner_proc
! Timings:
      integer                            :: tb, te, clock_rate, clock_max

      nrows   = size(LU,1)        ! Number of local rows
      neq     = size(LU,2)        ! Number of equations
      nprocs  = num_images()      ! Number of processors
      my_proc = this_image()      ! This processor
      offset  = nrows*(my_proc-1) ! Last equation on previous processor

      call system_clock ( tb, clock_rate, clock_max )

      do k = 1, neq

! Find the owner of row k:
         owner_proc = (k-1)/nrows+1

! Get the pivot row:
         if ( my_proc == owner_proc ) pivot_row(k:neq) = LU(k-offset,k:neq)

! Send pivot row to all processors:
         call co_broadcast( pivot_row(k:neq), source_image = owner_proc )

! Stop if the pivot is zero:
         if ( pivot_row(k) == 0._rp ) error stop 'Pivot equal to zero in LU decomposition'

! Eliminate nonzero's below pivot:
         do concurrent ( i = 1:nrows )
            if ( offset+i > k .and. LU(i,k)/= 0._rp) then
               LU(i,k) = LU(i,k)/pivot_row(k)
               LU(i,k+1:neq) = LU(i,k+1:neq) - LU(i,k)*pivot_row(k+1:neq)
            end if
         end do
      end do

      call system_clock ( te, clock_rate, clock_max )
      if ( my_proc == 1 ) &
         write(*,'(a,f8.2,a,/)') 'Time for LU decomposition = ', real(te-tb)/real(clock_rate), 's.'

   end subroutine RLU_decomp

   subroutine CLU_decomp( LU )
!
! Make in place LU decomposition of the dense matrix LU
! This is the complex version.
!
      complex(kind=cp), intent(inout)       :: LU(:,:)
      complex(kind=cp)                      :: pivot_row(size(LU,2))
      integer                               :: i, k, neq, nrows, offset
      integer                               :: nprocs, my_proc, owner_proc
! Timings:
      integer                               :: tb, te, clock_rate, clock_max

      nrows   = size(LU,1)        ! Number of local rows
      neq     = size(LU,2)        ! Number of equations
      nprocs  = num_images()      ! Number of processors
      my_proc = this_image()      ! This processor
      offset  = nrows*(my_proc-1) ! Last equation on previous processor

      call system_clock ( tb, clock_rate, clock_max )

      do k = 1, neq

! Find the owner of row k:
         owner_proc = (k-1)/nrows+1

! Get the pivot row:
         if ( my_proc == owner_proc ) pivot_row(k:neq) = LU(k-offset,k:neq)

! Send pivot row to all processors:
         call co_broadcast( pivot_row(k:neq), source_image = owner_proc )

! Stop if the pivot is zero:
         if ( pivot_row(k) == 0._cp ) error stop 'Pivot equal to zero in LU decomposition'

! Eliminate nonzero's below pivot:
         do concurrent ( i = 1:nrows )
            if ( offset+i > k .and. LU(i,k)/= 0._cp) then
               LU(i,k) = LU(i,k)/pivot_row(k)
               LU(i,k+1:neq) = LU(i,k+1:neq) - LU(i,k)*pivot_row(k+1:neq)
            end if
         end do
      end do

      call system_clock ( te, clock_rate, clock_max )
      if ( my_proc == 1 ) &
         write(*,'(a,f8.2,a,/)') 'Time for LU decomposition = ', real(te-tb)/real(clock_rate), 's.'

   end subroutine CLU_decomp

   function RLU_solve( b, LU ) result(x)
!
! Solve system LUx = b.
! This is the real version.
!
      implicit none
      real(kind=rp), intent(in)             :: LU(:,:)
      real(kind=rp), intent(in)             :: b(:)
      real(kind=rp)                         :: x(size(b))

      real(kind=rp)                         :: x_pivot
      integer                               :: i, k, neq, nrows, offset
      integer                               :: i_proc, my_proc, owner_proc

      nrows   = size(LU,1)        ! Number of local rows
      neq     = size(LU,2)        ! Number of equations
      my_proc = this_image()      ! This processor
      offset  = nrows*(my_proc-1) ! Last equation on previous processor

      x = b

! Forwards substitution
      do k = 1, neq

! Find the owner of row k:
         owner_proc = (k-1)/nrows+1

         if ( my_proc == owner_proc ) x_pivot = x(k-offset)
         call co_broadcast(x_pivot, source_image=owner_proc )

         if ( my_proc == owner_proc ) then
            do i = k-offset+1, nrows
               x(i) = x(i) - x_pivot * LU(i,k)
            end do
         elseif ( my_proc > owner_proc ) then
            do i = 1, nrows
               x(i) = x(i) - x_pivot * LU(i,k)
            end do
         end if
      end do

! Back substitution
      do k = neq, 1, -1

! Find the owner of row k:
         owner_proc = (k-1)/nrows+1

         if ( my_proc == owner_proc ) then
            x(k-offset) = x(k-offset)/LU(k-offset,k)
            x_pivot = x(k-offset)
         end if
         call co_broadcast(x_pivot, source_image=owner_proc )

         if ( my_proc == owner_proc ) then
            do i = k-offset-1, 1, -1
               x(i) = x(i) - x_pivot * LU(i,k)
            end do
         elseif ( my_proc < owner_proc ) then
            do i = nrows, 1, -1
               x(i) = x(i) - x_pivot * LU(i,k)
            end do
         end if
      end do
   end function RLU_solve

   function CLU_solve( b, LU ) result(x)
!
! Solve system LUx = b.
! This is the real version.
!
      implicit none
      complex(kind=cp), intent(in)          :: LU(:,:)
      complex(kind=cp), intent(in)          :: b(:)
      complex(kind=cp)                      :: x(size(b))

      complex(kind=cp)                      :: x_pivot
      integer                               :: i, k, neq, nrows, offset
      integer                               :: my_proc, owner_proc

      nrows   = size(LU,1)        ! Number of local rows
      neq     = size(LU,2)        ! Number of equations
      my_proc = this_image()      ! This processor
      offset  = nrows*(my_proc-1) ! Last equation on previous processor
      x = b

! Forwards substitution
      do k = 1, neq

! Find the owner of row k:
         owner_proc = (k-1)/nrows+1

         if ( my_proc == owner_proc ) x_pivot = x(k-offset)
         call co_broadcast(x_pivot, source_image = owner_proc )

         if ( my_proc == owner_proc ) then
            do i = k-offset+1, nrows
               x(i) = x(i) - x_pivot * LU(i,k)
            end do
         elseif ( my_proc > owner_proc ) then
            do i = 1, nrows
               x(i) = x(i) - x_pivot * LU(i,k)
            end do
         end if
      end do

! Back substitution
      do k = neq, 1, -1

! Find the owner of row k:
         owner_proc = (k-1)/nrows+1

         if ( my_proc == owner_proc ) then
            x(k-offset) = x(k-offset)/LU(k-offset,k)
            x_pivot = x(k-offset)
         end if
         call co_broadcast(x_pivot, source_image=owner_proc )

         if ( my_proc == owner_proc ) then
            do i = k-offset-1, 1, -1
               x(i) = x(i) - x_pivot * LU(i,k)
            end do
         elseif ( my_proc < owner_proc ) then
            do i = nrows, 1, -1
               x(i) = x(i) - x_pivot * LU(i,k)
            end do
         end if

      end do
   end function CLU_solve

   function get_band( A ) result(half_band)
!
! Determine half bandwidth of matrix
!
      type(user_matrix), intent(in)         :: A
      integer                               :: half_band
      integer                               :: my_proc, row, col, dia, i, j, k

      my_proc = this_image()
      half_band = 0
 
      if ( A%matrix_format == 1 ) then
! DENSE format
         half_band = A%ncols
      elseif ( A%matrix_format == 2 ) then
! CRS format
         do i = 1, A%nrows
            do j = A%row_ptr(i), A%row_ptr(i+1)-1
               row = i+(my_proc-1)*A%nrows
               col = A%col_ind(j)+(A%phi(j)-1)*A%nrows
               dia = col-row
               half_band = max(half_band,abs(dia))
            end do
         end do
         call co_max(half_band)
      elseif ( A%matrix_format == 3 ) then
! COO format
         do k = 1, A%nnz
            row = A%indx(k)+(my_proc-1)*A%nrows
            col = A%jndx(k)+(A%phi(k)-1)*A%nrows
            dia = col-row
            half_band = max(half_band,abs(dia))
         end do
         call co_max(half_band)
      elseif ( A%matrix_format == 4 ) then
! CDS format
          half_band = A%half_band
      end if
    
   end function get_band

   function complex_band( A, half_band ) result(LU)
!
! Convert matrix from sparse format (CRS or COO) to band format (CDS)
!

      type(user_matrix), intent(in)         :: A
      integer, intent(in)                   :: half_band
      complex(kind=cp)                      :: LU(A%nrows,-half_band:half_band)
      integer                               :: i, j, k, dia, my_proc, row, col

      my_proc = this_image()
      LU = 0._cp
      if ( A%matrix_format == 2 ) then
! CRS format
         if ( allocated(A%rval) ) then
! Matrix is real
            do i = 1, A%nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  row = i+(my_proc-1)*A%nrows
                  col = A%col_ind(j)+(A%phi(j)-1)*A%nrows
                  dia = col-row
                  LU(i,dia) = LU(i,dia) + A%rval(j)
               end do
            end do
         elseif  ( allocated(A%cval) ) then
! Matrix is complex
            do i = 1, A%nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  row = i+(my_proc-1)*A%nrows
                  col = A%col_ind(j)+(A%phi(j)-1)*A%nrows
                  dia = col-row
                  LU(i,dia) = LU(i,dia) + A%cval(j)
               end do
            end do
         elseif  ( allocated(A%ival) ) then
! Matrix is integer
            do i = 1, A%nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  row = i+(my_proc-1)*A%nrows
                  col = A%col_ind(j)+(A%phi(j)-1)*A%nrows
                  dia = col-row
                  LU(i,dia) = LU(i,dia) + A%ival(j)
               end do
            end do
         else
! Matrix is binary
            do i = 1, A%nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  row = i+(my_proc-1)*A%nrows
                  col = A%col_ind(j)+(A%phi(j)-1)*A%nrows
                  dia = col-row
                  LU(i,dia) = LU(i,dia) + 1._cp
               end do
            end do
         end if
      elseif ( A%matrix_format == 3 ) then
! COO format
         if ( allocated(A%rval) ) then
! Matrix is real
            do k = 1,A%nnz
               i   = A%indx(k)
               row = A%indx(k)+(my_proc-1)*A%nrows
               col = A%jndx(k)+(A%phi(k)-1)*A%nrows
               dia = col-row
               LU(i,dia) = LU(i,dia) + A%rval(k)
            end do
         elseif  ( allocated(A%cval) ) then
! Matrix is complex
            do k = 1,A%nnz
               i   = A%indx(k)
               row = A%indx(k)+(my_proc-1)*A%nrows
               col = A%jndx(k)+(A%phi(k)-1)*A%nrows
               dia = col-row
               LU(i,dia) = LU(i,dia) + A%cval(k)
            end do
         elseif  ( allocated(A%ival) ) then
! Matrix is integer
            do k = 1,A%nnz
               i   = A%indx(k)
               row = A%indx(k)+(my_proc-1)*A%nrows
               col = A%jndx(k)+(A%phi(k)-1)*A%nrows
               dia = col-row
               LU(i,dia) = LU(i,dia) + A%ival(k)
            end do
         else
! Matrix is binary
            do k = 1,A%nnz
               i   = A%indx(k)
               row = A%indx(k)+(my_proc-1)*A%nrows
               col = A%jndx(k)+(A%phi(k)-1)*A%nrows
               dia = col-row
               LU(i,dia) = LU(i,dia) + 1._cp
            end do
         end if
      elseif ( A%matrix_format == 4 ) then
         if ( allocated(A%rband) ) then
            LU = A%rband
         elseif ( allocated(A%iband) ) then
            LU = A%iband
         elseif ( allocated(A%cband) ) then
            LU = A%cband
         else
            LU = 1._rp
         end if
      end if
   end function complex_band

   function real_band( A, half_band ) result(LU)
!
! Convert matrix from sparse format (CRS or COO) to band format (CDS)
!

      type(user_matrix), intent(in)         :: A
      integer, intent(in)                   :: half_band
      real(kind=rp)                         :: LU(A%nrows,-half_band:half_band)
      integer                               :: i, j, k, dia, my_proc, row, col

      my_proc = this_image()
      LU = 0._rp

      if ( A%matrix_format == 2 ) then
! CRS format
         if ( allocated(A%rval) ) then
! Matrix is real
            do i = 1, A%nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  row = i+(my_proc-1)*A%nrows
                  col = A%col_ind(j)+(A%phi(j)-1)*A%nrows
                  dia = col-row
                  LU(i,dia) = LU(i,dia) + A%rval(j)
               end do
            end do
         elseif  ( allocated(A%ival) ) then
! Matrix is integer
            do i = 1, A%nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  row = i+(my_proc-1)*A%nrows
                  col = A%col_ind(j)+(A%phi(j)-1)*A%nrows
                  dia = col-row
                  LU(i,dia) = LU(i,dia) + A%ival(j)
               end do
            end do
         else
! Matrix is binary
            do i = 1, A%nrows
               do j = A%row_ptr(i), A%row_ptr(i+1)-1
                  row = i+(my_proc-1)*A%nrows
                  col = A%col_ind(j)+(A%phi(j)-1)*A%nrows
                  dia = col-row
                  LU(i,dia) = LU(i,dia) + 1._rp
               end do
            end do
         end if
      elseif ( A%matrix_format == 3 ) then
! COO format
         if ( allocated(A%rval) ) then
! Matrix is real
            do k = 1,A%nnz
               i   = A%indx(k)
               row = A%indx(k)+(my_proc-1)*A%nrows
               col = A%jndx(k)+(A%phi(k)-1)*A%nrows
               dia = col-row
               LU(i,dia) = LU(i,dia) + A%rval(k)
            end do
         elseif  ( allocated(A%ival) ) then
! Matrix is integer
            do k = 1,A%nnz
               i   = A%indx(k)
               row = A%indx(k)+(my_proc-1)*A%nrows
               col = A%jndx(k)+(A%phi(k)-1)*A%nrows
               dia = col-row
               LU(i,dia) = LU(i,dia) + A%ival(k)
            end do
         else
! Matrix is binary
            do k = 1,A%nnz
               i   = A%indx(k)
               row = A%indx(k)+(my_proc-1)*A%nrows
               col = A%jndx(k)+(A%phi(k)-1)*A%nrows
               dia = col-row
               LU(i,dia) = LU(i,dia) + 1._cp
            end do
         end if
      elseif ( A%matrix_format == 4 ) then
         if ( allocated(A%rband) ) then
            LU = A%rband
         elseif ( allocated(A%iband) ) then
            LU = A%iband
         else
            LU = 1._rp
         end if
      end if
   end function real_band

   subroutine rband_LU_decomp( LU, neq, half_band )

      integer, intent(in)                :: neq, half_band
      real(kind=rp), intent(inout)       :: LU(:,-half_band:)
      real(kind=rp), allocatable         :: pivot_row(:)[:]

      integer                            :: nrows, start_row, dia, offset, i, j, k
      integer                            :: owner_proc, my_proc, n_procs
      integer                            :: last_row, owner_last
! Timings:
      integer                            :: tb, te, clock_rate, clock_max

      n_procs   = num_images()
      my_proc   = this_image()
      nrows     = size(LU,1)
      offset    = nrows*(my_proc-1) ! Last equation on previous processor

      call system_clock ( tb, clock_rate, clock_max )

      allocate(pivot_row(0:half_band)[*])

      do k = 1, neq-1

! Find the owner of row k:
         owner_proc = (k-1)/nrows+1

         pivot_row = 0._rp
         last_row = min(neq, k+half_band )
         owner_last = (last_row-1)/nrows+1

         if ( my_proc == owner_proc ) then
            pivot_row(0:half_band) = LU(k-offset, 0:half_band)
            if ( pivot_row(0) == 0._rp ) error stop 'Pivot equal to zero in LU decomposition'
         end if
         sync all
         if ( my_proc > owner_proc .and. my_proc <= owner_last ) &
            pivot_row(0:half_band) = pivot_row(:)[owner_proc]
         sync all

         if (my_proc >= owner_proc .and. my_proc <= owner_last) then
            do concurrent (i = 1:nrows)
               dia = k - (offset + i)
               if (offset + i > k .and. -dia <= half_band ) then
                  LU(i,dia) = LU(i,dia)/pivot_row(0)
                  LU(i,dia+1:dia+half_band) = LU(i,dia+1:dia+half_band) - LU(i,dia)*pivot_row(1:half_band)
               end if
            end do
         end if
      end do

      call system_clock ( te, clock_rate, clock_max )
      if ( my_proc == 1 ) &
         write(*,'(a,f8.2,a,/)') 'Time for LU decomposition = ', real(te-tb)/real(clock_rate), 's.'

   end subroutine rband_LU_decomp

   subroutine cband_LU_decomp( LU, neq, half_band )

      implicit none
      integer, intent(in)                :: neq, half_band
      complex(kind=cp), intent(inout)    :: LU(:,-half_band:)
      complex(kind=cp), allocatable      :: pivot_row(:)[:]

      integer                            :: nrows, start_row, dia, offset, i, j, k
      integer                            :: owner_proc, my_proc, n_procs
      integer                            :: global_row, last_row, owner_last
! Timings:
      integer                            :: tb, te, clock_rate, clock_max

      n_procs   = num_images()
      my_proc   = this_image()
      nrows     = size(LU,1)
      offset    = nrows*(my_proc-1) ! Last equation on previous processor

      call system_clock ( tb, clock_rate, clock_max )

      allocate(pivot_row(0:half_band)[*])

      do k = 1, neq-1

! Find the owner of row k:
         owner_proc = (k-1)/nrows+1

         pivot_row = 0._cp
         last_row = min(neq, k+half_band )
         owner_last = (last_row-1)/nrows+1

         if (my_proc == owner_proc ) then 
            pivot_row(0:half_band) = LU(k-offset, 0:half_band)
            if ( pivot_row(0) == 0._cp ) error stop 'Pivot equal to zero in LU decomposition'
         end if
         sync all
         if ( my_proc > owner_proc .and. my_proc <= owner_last ) &
            pivot_row(0:half_band) = pivot_row(:)[owner_proc]
         sync all

         if (my_proc >= owner_proc .and. my_proc <= owner_last) then
            do concurrent (i = 1:nrows)
               dia = k - (offset + i)
               if (offset+i > k .and. -dia <= half_band ) then
                  LU(i,dia) = LU(i,dia)/pivot_row(0)
                  LU(i,dia+1:dia+half_band) = LU(i,dia+1:dia+half_band) - LU(i,dia)*pivot_row(1:half_band)
               end if
            end do
         end if
      end do

      call system_clock ( te, clock_rate, clock_max )
      if ( my_proc == 1 ) &
         write(*,'(a,f8.2,a,/)') 'Time for LU decomposition = ', real(te-tb)/real(clock_rate), 's.'

   end subroutine cband_LU_decomp

   function rband_LU_solve( b, LU, half_band, neq ) result(x)
!
! Solve system LUx = b.
! This is the real version.
!
      implicit none
      integer                               :: half_band, neq
      real(kind=rp), intent(in)             :: LU(:,-half_band:)
      real(kind=rp), intent(in)             :: b(:)
      real(kind=rp)                         :: x(size(b))

      integer                               :: i, k, nrows, offset, dia, up, low
      integer                               :: my_proc, owner_proc, n_procs
      real(kind=rp)                         :: x_pivot


      nrows   = size(LU,1)        ! Number of local rows
      my_proc = this_image()      ! This processor
      n_procs = num_images()      ! Number of images
      offset  = nrows*(my_proc-1) ! Last equation on previous processor

      x = b

! Forwards substitution
      do k = 1, neq

! Find the owner of row k:
         owner_proc = (k-1)/nrows+1

         if ( my_proc == owner_proc ) x_pivot = x(k-offset)
         if ( n_procs > 1 ) call co_broadcast(x_pivot, source_image = owner_proc )

         if ( my_proc == owner_proc ) then
            up = min( k-offset+half_band,nrows)
            do i = k-offset+1, up
               dia = k - (offset + i)
               x(i) = x(i) - x_pivot * LU(i,dia)
            end do
         elseif ( my_proc > owner_proc ) then
            up = min( k-offset+half_band,nrows)
            do i = 1, up
               dia = k - (offset + i)
               x(i) = x(i) - x_pivot * LU(i,dia)
            end do
         end if

      end do

! Back substitution
      do k = neq, 1, -1

! Find the owner of row k:
         owner_proc = (k-1)/nrows+1

         if ( my_proc == owner_proc ) then
            x(k-offset) = x(k-offset)/LU(k-offset,0)
            x_pivot = x(k-offset)
         end if

         if ( n_procs > 1 ) call co_broadcast(x_pivot, source_image = owner_proc )

         if ( my_proc == owner_proc ) then
            low = max(k-offset-half_band,1)
            do i = k-offset-1, low, -1
               dia = k - (offset + i)
               x(i) = x(i) - x_pivot * LU(i,dia)
            end do
         elseif ( my_proc < owner_proc ) then
            low = max(k-offset-half_band,1)
            do i = nrows, low, -1
               dia = k - (offset + i)
               x(i) = x(i) - x_pivot * LU(i,dia)
            end do
         end if

      end do

   end function rband_LU_solve

   function cband_LU_solve( b, LU, half_band, neq ) result(x)
!
! Solve system LUx = b.
! This is the complex version.
!
      implicit none
      integer                               :: half_band, neq
      complex(kind=cp), intent(in)          :: LU(:,-half_band:)
      complex(kind=cp), intent(in)          :: b(:)
      complex(kind=cp)                      :: x(size(b))

      integer                               :: i, k, nrows, offset, dia, up, low
      integer                               :: my_proc, owner_proc, n_procs
      complex(kind=cp)                      :: x_pivot

      nrows   = size(LU,1)        ! Number of local rows
      my_proc = this_image()      ! This processor
      n_procs = num_images()      ! Number of images
      offset  = nrows*(my_proc-1) ! Last equation on previous processor

      x = b

! Forwards substitution
      do k = 1, neq

! Find the owner of row k:
         owner_proc = (k-1)/nrows+1

         if ( my_proc == owner_proc ) x_pivot = x(k-offset)
         if ( n_procs > 1 ) call co_broadcast(x_pivot, source_image = owner_proc )

         if ( my_proc == owner_proc ) then
            up = min( k-offset+half_band,nrows)
            do i = k-offset+1, up
               dia = k - (offset + i)
               x(i) = x(i) - x_pivot * LU(i,dia)
            end do
         elseif ( my_proc > owner_proc ) then
            up = min( k-offset+half_band,nrows)
            do i = 1, up
               dia = k - (offset + i)
               x(i) = x(i) - x_pivot * LU(i,dia)
            end do
         end if

      end do

! Back substitution
      do k = neq, 1, -1

! Find the owner of row k:
         owner_proc = (k-1)/nrows+1

         if ( my_proc == owner_proc ) then
            x(k-offset) = x(k-offset)/LU(k-offset,0)
            x_pivot = x(k-offset)
         end if

         if ( n_procs > 1 ) call co_broadcast(x_pivot, source_image = owner_proc )

         if ( my_proc == owner_proc ) then
            low = max(k-offset-half_band,1)
            do i = k-offset-1, low, -1
               dia  = k - (offset + i)
               x(i) = x(i) - x_pivot * LU(i,dia)
            end do
         elseif ( my_proc < owner_proc ) then
            low = max(k-offset-half_band,1)
            do i = nrows, low, -1
               dia  = k - (offset + i)
               x(i) = x(i) - x_pivot * LU(i,dia)
            end do
         end if

      end do

   end function cband_LU_solve

end module user_module

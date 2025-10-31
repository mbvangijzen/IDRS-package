module mm_module

!
! Module to read a matrix from an matrix-market file. 
! The matrix is changed to CRS format, row-partition for parallel processing, 
! and then store the result in the type user-matrix.
!
! Function to read a matrix (or vector): 
!   A = mm_matrix( mm_file )
!     mm_file: character(len=*), input,  Name (including path) of mm_file
!     A:       type(user_matrix) output, Matrix stored in (pre-defined) user_matrix type
!
! This module uses the modules:
!    user_module: defines the user_matrix type 
!    precision_module: defines the real and complex precision
!
! This software is distributed under the MIT License:
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2024 Martin van Gijzen
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 

   use precision_module
   use user_module
   implicit none

   contains

   function mm_matrix( mm_file, which_format ) result(A)

      implicit none 
      character(len=*)              :: mm_file
      type(user_matrix)             :: A
      
      character(len=3), optional    :: which_format
      character(len=3)              :: format

      integer                       :: cols, rows, nnz, nnzmax 
      integer                       :: Allocate_status     
      integer, allocatable          :: ival(:), indx(:), jndx(:)
      double precision, allocatable :: rval(:)
      complex, allocatable          :: cval(:)
      integer                       :: l_nnz, nrows, ncols
      integer, allocatable          :: row_ptr(:), col_ind(:)
      integer, allocatable          :: l_indx(:), l_jndx(:)
      real(kind=rp), allocatable    :: l_rval(:), l_rmat(:,:)
      integer, allocatable          :: l_ival(:), l_imat(:,:)
      complex(kind=cp), allocatable :: l_cval(:), l_cmat(:,:)
      character(len=10)             :: rep
      character(len=7)              :: field
      character(len=19)             :: symm
!
      integer                       :: my_proc, n_procs
      integer                       :: i, j, k, ierr

      if ( .not. present(which_format) ) then
         format = 'coo'
      else
         format = which_format
      end if
!
      open(unit = 10, file = mm_file, status = 'old', iostat = ierr )
      if ( ierr > 0 ) stop 'Error while opening matrix-market file '
      my_proc = this_image()
      n_procs = num_images()
!  
! Read header to determine which space is needed
      call mminfo(10,rep,field,symm,rows,cols,nnzmax)
!
! Allocate the space
      allocate( indx(nnzmax), jndx(nnzmax), stat = Allocate_status)
      allocate( rval(nnzmax), stat = Allocate_status)
      allocate( ival(nnzmax), stat = Allocate_status)
      allocate( cval(nnzmax), stat = Allocate_status)
      call mmread(10,rep,field,symm,rows,cols,nnz,nnzmax, indx,jndx,ival,rval,cval )
      close(10)
!
      nrows = ceiling(real(rows)/real(n_procs))    ! Local rows
      ncols = cols

      if ( rep == 'coordinate' ) then
         if ( format == 'coo' ) then
!
! Change to coo-format, only for equation on this processor
            call mm_coo( nrows, ncols, nnz, ival, indx, jndx, rval, cval, & 
                         my_proc, n_procs, symm, field, &
                         l_nnz, l_indx, l_jndx, l_rval, l_ival, l_cval )
!
! Store the matrix coefficients in user_matrix::
            if ( field == 'complex' ) then
               A = coo_format( nrows, ncols, l_nnz, l_indx, l_jndx, cval = l_cval )
            elseif ( field == 'integer' ) then
               A = coo_format( nrows, ncols, l_nnz, l_indx, l_jndx, ival = l_ival )
            elseif ( field == 'real' ) then
               A = coo_format( nrows, ncols, l_nnz, l_indx, l_jndx, rval = l_rval )
            else
               A = coo_format( nrows, ncols, l_nnz, l_indx, l_jndx )
            end if

         else
!
! Change to crs-format, only for equation on this processor
            call mm_crs( nrows, ncols, nnz, ival, indx, jndx, rval, cval, & 
                         my_proc, n_procs, symm, field, &
                         row_ptr, col_ind, l_rval, l_ival, l_cval )
!
! Store the matrix coefficients in user_matrix::
            if ( field == 'complex' ) then
               A = crs_format( nrows, ncols, row_ptr, col_ind, cval = l_cval )
            elseif ( field == 'integer' ) then
               A = crs_format( nrows, ncols, row_ptr, col_ind, ival = l_ival )
            elseif ( field == 'real' ) then
               A = crs_format( nrows, ncols, row_ptr, col_ind, rval = l_rval )
            else
               A = crs_format( nrows, ncols, row_ptr, col_ind )
            end if
         end if
      else
!
! Change to standard dense-format, only for equation on this processor
         call mm_dense( nrows, ncols, ival, rval, cval, & 
                        my_proc, n_procs, symm, field, rows, &
                        l_rmat, l_imat, l_cmat )
!
! Store the matrix coefficients in user_matrix::
         if ( field == 'complex' ) then
            A = dense_format( cmat = l_cmat )
         elseif ( field == 'integer' ) then
            A = dense_format( imat = l_imat )
         elseif ( field == 'real' ) then
            A = dense_format( rmat = l_rmat )
         end if
      end if

   end function mm_matrix

   subroutine mm_dense( nrows, ncols, ival, rval, cval, & 
                        my_proc, n_procs, symm, field, rows, &
                        l_rmat, l_imat, l_cmat )
!
! Subroutine to convert matrix-market dense format to
! standard dense format. This is done only for the nrows assigned to
! this processor. Variables starting with l_ are output and unique for the local
! processor.
!

      integer                                    :: ncols, nrows         ! local rows and columns
      integer, intent(in)                        :: ival(:)
      double precision, intent(in)               :: rval(:)
      complex, intent(in)                        :: cval(:)
      integer, intent(in)                        :: my_proc, n_procs
      character(len=7), intent(in)               :: field
      character(len=19), intent(in)              :: symm
      integer, intent(in)                        :: rows                 ! global rows
!
      real(kind=rp), allocatable, intent(out)    :: l_rmat(:,:)
      integer, allocatable, intent(out)          :: l_imat(:,:)
      complex(kind=cp), allocatable, intent(out) :: l_cmat(:,:)
!
      integer                                    :: i, j, k

      if ( field == 'real' ) then
         allocate(l_rmat(nrows,ncols))
         l_rmat = 0._rp
         if ( symm == 'general' ) then 
            k = 1
            do j = 1, ncols
               do i = 1, rows
                  if ( ( i > (my_proc-1)*nrows ) .and. ( i <= my_proc*nrows )  ) &
                     l_rmat(i-(my_proc-1)*nrows,j) = rval(k)
                  k = k + 1
               end do
            end do
         elseif ( symm == 'symmetric' ) then
            k = 1
            do j = 1, ncols
               do i = j, rows
                  if ( ( i > (my_proc-1)*nrows ) .and. ( i <= my_proc*nrows )  ) &
                     l_rmat(i-(my_proc-1)*nrows,j) = rval(k)
                  if ( ( j > (my_proc-1)*nrows ) .and. ( j <= my_proc*nrows )  ) &
                     l_rmat(j-(my_proc-1)*nrows,i) = rval(k)
                  k = k + 1
               end do
            end do
         elseif ( symm == 'skew-symmetric' ) then
            k = 1
            do j = 1, ncols
               do i = j+1, rows
                  if ( ( i > (my_proc-1)*nrows ) .and. ( i <= my_proc*nrows )  ) &
                     l_rmat(i-(my_proc-1)*nrows,j) =  rval(k)
                  if ( ( j > (my_proc-1)*nrows ) .and. ( j <= my_proc*nrows )  ) &
                     l_rmat(j-(my_proc-1)*nrows,i) = -rval(k)
                  k = k + 1
               end do
            end do
         end if
      elseif ( field == 'integer' ) then
         allocate(l_imat(nrows,ncols))
         l_imat = 0
         if ( symm == 'general' ) then 
            k = 1
            do j = 1, ncols
               do i = 1, rows
                  if ( ( i > (my_proc-1)*nrows ) .and. ( i <= my_proc*nrows )  ) &
                     l_imat(i-(my_proc-1)*nrows,j) = ival(k)
                  k = k + 1
               end do
            end do
         elseif ( symm == 'symmetric' ) then
            k = 1
            do j = 1, ncols
               do i = j, rows
                  if ( ( i > (my_proc-1)*nrows ) .and. ( i <= my_proc*nrows )  ) &
                     l_imat(i-(my_proc-1)*nrows,j) = ival(k)
                  if ( ( j > (my_proc-1)*nrows ) .and. ( j <= my_proc*nrows )  ) &
                     l_imat(j-(my_proc-1)*nrows,i) = ival(k)
                  k = k + 1
               end do
            end do
         elseif ( symm == 'skew-symmetric' ) then
            k = 1
            do j = 1, ncols
               do i = j+1, rows
                  if ( ( i > (my_proc-1)*nrows ) .and. ( i <= my_proc*nrows )  ) &
                     l_imat(i-(my_proc-1)*nrows,j) =  ival(k)
                  if ( ( j > (my_proc-1)*nrows ) .and. ( j <= my_proc*nrows )  ) &
                     l_imat(j-(my_proc-1)*nrows,i) = -ival(k)
                  k = k + 1
               end do
            end do
         end if
      elseif ( field == 'complex' ) then
         allocate(l_cmat(nrows,ncols))
         l_cmat = (0._cp,0._cp)
         if ( symm == 'general' ) then 
            k = 1
            do j = 1, ncols
               do i = 1, rows
                  if ( ( i > (my_proc-1)*nrows ) .and. ( i <= my_proc*nrows )  ) &
                     l_cmat(i-(my_proc-1)*nrows,j) = cval(k)
                  k = k + 1
               end do
            end do
         elseif ( symm == 'symmetric' ) then
            k = 1
            do j = 1, ncols
               do i = j, rows
                  if ( ( i > (my_proc-1)*nrows ) .and. ( i <= my_proc*nrows )  ) &
                     l_cmat(i-(my_proc-1)*nrows,j) = cval(k)
                  if ( ( j > (my_proc-1)*nrows ) .and. ( j <= my_proc*nrows )  ) &
                     l_cmat(j-(my_proc-1)*nrows,i) = cval(k)
                  k = k + 1
               end do
            end do
         elseif ( symm == 'skew-symmetric' ) then
            k = 1
            do j = 1, ncols
               do i = j+1, rows
                  if ( ( i > (my_proc-1)*nrows ) .and. ( i <= my_proc*nrows )  ) &
                     l_cmat(i-(my_proc-1)*nrows,j) =  cval(k)
                  if ( ( j > (my_proc-1)*nrows ) .and. ( j <= my_proc*nrows )  ) &
                     l_cmat(j-(my_proc-1)*nrows,i) = -cval(k)
                  k = k + 1
               end do
            end do
         elseif ( symm == 'hermitian' ) then 
            k = 1
            do j = 1,ncols
               do i = j, rows
                  if ( ( i > (my_proc-1)*nrows ) .and. ( i <= my_proc*nrows )  ) &
                     l_cmat(i-(my_proc-1)*nrows,j) =  cval(k)
                  if ( ( j > (my_proc-1)*nrows ) .and. ( j <= my_proc*nrows )  ) &
                     l_cmat(j-(my_proc-1)*nrows,i) = conjg(cval(k))
                  k = k + 1
               end do
            end do
         end if
      end if

   end subroutine  mm_dense 

   subroutine mm_crs( nrows, ncols, nnz, ival, indx, jndx, rval, cval, &
                      my_proc, n_procs, symm, field, &
                      row_ptr, col_ind, l_rval, l_ival, l_cval )
!
! Subroutine to convert matrix-market sparse format to
! crs format. This is done only for the nrows assigned to
! this processor. Variables starting with l_ are output and unique for the local
! processor.
!
      integer, intent(in)                        :: nrows, ncols, nnz
      integer, intent(in)                        :: ival(:), indx(:), jndx(:)
      double precision, intent(in)               :: rval(:)
      complex, intent(in)                        :: cval(:)
      integer, intent(in)                        :: my_proc, n_procs
      character(len=7), intent(in)               :: field
      character(len=19), intent(in)              :: symm
!
      integer, allocatable, intent(out)          :: row_ptr(:), col_ind(:)
      integer, allocatable                       :: col_ptr(:)
      real(kind=rp), allocatable, intent(out)    :: l_rval(:)
      integer, allocatable, intent(out)          :: l_ival(:)
      complex(kind=cp), allocatable, intent(out) :: l_cval(:)
      integer                                    :: l_nnz, l_row
!
      integer                                    :: i, j, k

      allocate(row_ptr(nrows+1),col_ptr(nrows))
      row_ptr = 0
      col_ptr = 0

      if ( symm == 'general' ) then
         k = 0
         do i = 1, nnz
            if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
               k = k + 1
               l_row = indx(i) - (my_proc-1)*nrows
               row_ptr(l_row+1) = row_ptr(l_row+1) + 1
            end if
         end do
      else
         k = 0
         do i = 1, nnz
            if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
               k = k + 1
               l_row = indx(i) - (my_proc-1)*nrows
               row_ptr(l_row+1) = row_ptr(l_row+1) + 1
            end if
            if ( indx(i) .ne. jndx(i) ) then
               if ( jndx(i) > ( my_proc-1)*nrows .and. jndx(i) <= my_proc*nrows ) then
                  k = k + 1
                  l_row = jndx(i) - (my_proc-1)*nrows
                  row_ptr(l_row+1) = row_ptr(l_row+1) + 1
               end if
            end if
         end do
      end if
      row_ptr(1) = 1
      do i = 2, nrows+1
         row_ptr(i) = row_ptr(i) + row_ptr(i-1)
      end do

      l_nnz = k
      allocate( col_ind(l_nnz) )
!
      if ( field == 'real' ) then
         allocate( l_rval(l_nnz) )
         if ( symm == 'general' ) then
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_row = indx(i) - (my_proc-1)*nrows
                  col_ind(row_ptr(l_row)+col_ptr(l_row)) = jndx(i)
                  l_rval(row_ptr(l_row)+col_ptr(l_row))  = rval(i)
                  col_ptr(l_row) = col_ptr(l_row) + 1
               end if
            end do
         elseif ( symm == 'symmetric' ) then
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_row = indx(i) - (my_proc-1)*nrows
                  col_ind(row_ptr(l_row)+col_ptr(l_row)) = jndx(i)
                  l_rval(row_ptr(l_row)+col_ptr(l_row))  = rval(i)
                  col_ptr(l_row) = col_ptr(l_row) + 1
               end if
               if ( indx(i) .ne. jndx(i) ) then
                  if ( jndx(i) > ( my_proc-1)*nrows .and. jndx(i) <= my_proc*nrows ) then
                     l_row = jndx(i) - (my_proc-1)*nrows
                     col_ind(row_ptr(l_row)+col_ptr(l_row)) = indx(i)
                     l_rval(row_ptr(l_row)+col_ptr(l_row))  = rval(i)
                     col_ptr(l_row) = col_ptr(l_row) + 1
                  end if
               end if
            end do
         elseif ( symm == 'skew-symmetric' ) then
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_row = indx(i) - (my_proc-1)*nrows
                  col_ind(row_ptr(l_row)+col_ptr(l_row)) = jndx(i)
                  l_rval(row_ptr(l_row)+col_ptr(l_row))  = rval(i)
                  col_ptr(l_row) = col_ptr(l_row) + 1
               end if
               if ( jndx(i) > ( my_proc-1)*nrows .and. jndx(i) <= my_proc*nrows ) then
                  l_row = jndx(i) - (my_proc-1)*nrows
                  col_ind(row_ptr(l_row)+col_ptr(l_row)) =  indx(i)
                  l_rval(row_ptr(l_row)+col_ptr(l_row))  = -rval(i)
                  col_ptr(l_row) = col_ptr(l_row) + 1
               end if
            end do
         end if
      elseif ( field == 'integer' ) then
         allocate( l_ival(l_nnz) )
         if ( symm == 'general' ) then
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_row = indx(i) - (my_proc-1)*nrows
                  col_ind(row_ptr(l_row)+col_ptr(l_row)) = jndx(i)
                  l_ival(row_ptr(l_row)+col_ptr(l_row))  = ival(i)
                  col_ptr(l_row) = col_ptr(l_row) + 1
               end if
            end do
         elseif ( symm == 'symmetric' ) then
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_row = indx(i) - (my_proc-1)*nrows
                  col_ind(row_ptr(l_row)+col_ptr(l_row)) = jndx(i)
                  l_ival(row_ptr(l_row)+col_ptr(l_row))  = ival(i)
                  col_ptr(l_row) = col_ptr(l_row) + 1
               end if
               if ( indx(i) .ne. jndx(i) ) then
                  if ( jndx(i) > ( my_proc-1)*nrows .and. jndx(i) <= my_proc*nrows ) then
                     l_row = jndx(i) - (my_proc-1)*nrows
                     col_ind(row_ptr(l_row)+col_ptr(l_row)) = indx(i)
                     l_ival(row_ptr(l_row)+col_ptr(l_row))  = ival(i)
                     col_ptr(l_row) = col_ptr(l_row) + 1
                  end if
               end if
            end do
         elseif ( symm == 'skew-symmetric' ) then
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_row = indx(i) - (my_proc-1)*nrows
                  col_ind(row_ptr(l_row)+col_ptr(l_row)) = jndx(i)
                  l_ival(row_ptr(l_row)+col_ptr(l_row))  = ival(i)
                  col_ptr(l_row) = col_ptr(l_row) + 1
               end if
               if ( jndx(i) > ( my_proc-1)*nrows .and. jndx(i) <= my_proc*nrows ) then
                  l_row = jndx(i) - (my_proc-1)*nrows
                  col_ind(row_ptr(l_row)+col_ptr(l_row)) =  indx(i)
                  l_ival(row_ptr(l_row)+col_ptr(l_row))  = -ival(i)
                  col_ptr(l_row) = col_ptr(l_row) + 1
               end if
            end do
         end if
      elseif ( field == 'pattern' ) then
         if ( symm == 'general' ) then
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_row = indx(i) - (my_proc-1)*nrows
                  col_ind(row_ptr(l_row)+col_ptr(l_row)) = jndx(i)
                  col_ptr(l_row) = col_ptr(l_row) + 1
               end if
            end do
         elseif ( symm == 'symmetric' ) then
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_row = indx(i) - (my_proc-1)*nrows
                  col_ind(row_ptr(l_row)+col_ptr(l_row)) = jndx(i)
                  col_ptr(l_row) = col_ptr(l_row) + 1
               end if
               if ( indx(i) .ne. jndx(i) ) then
                  if ( jndx(i) > ( my_proc-1)*nrows .and. jndx(i) <= my_proc*nrows ) then
                     l_row = jndx(i) - (my_proc-1)*nrows
                     col_ind(row_ptr(l_row)+col_ptr(l_row)) = indx(i)
                     col_ptr(l_row) = col_ptr(l_row) + 1
                  end if
               end if
            end do
         end if
      elseif ( field == 'complex' ) then
         allocate( l_cval(l_nnz) )
         if ( symm == 'general' ) then
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_row = indx(i) - (my_proc-1)*nrows
                  col_ind(row_ptr(l_row)+col_ptr(l_row)) = jndx(i)
                  l_cval(row_ptr(l_row)+col_ptr(l_row))  = cval(i)
                  col_ptr(l_row) = col_ptr(l_row) + 1
               end if
            end do
         elseif ( symm == 'symmetric' ) then
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_row = indx(i) - (my_proc-1)*nrows
                  col_ind(row_ptr(l_row)+col_ptr(l_row)) = jndx(i)
                  l_cval(row_ptr(l_row)+col_ptr(l_row))  = cval(i)
                  col_ptr(l_row) = col_ptr(l_row) + 1
               end if
               if ( indx(i) .ne. jndx(i) ) then
                  if ( jndx(i) > ( my_proc-1)*nrows .and. jndx(i) <= my_proc*nrows ) then
                     l_row = jndx(i) - (my_proc-1)*nrows
                     col_ind(row_ptr(l_row)+col_ptr(l_row)) = indx(i)
                     l_cval(row_ptr(l_row)+col_ptr(l_row))  = cval(i)
                     col_ptr(l_row) = col_ptr(l_row) + 1
                  end if
               end if
            end do
         elseif ( symm == 'skew-symmetric' ) then
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_row = indx(i) - (my_proc-1)*nrows
                  col_ind(row_ptr(l_row)+col_ptr(l_row)) = jndx(i)
                  l_cval(row_ptr(l_row)+col_ptr(l_row))  = cval(i)
                  col_ptr(l_row) = col_ptr(l_row) + 1
               end if
               if ( jndx(i) > ( my_proc-1)*nrows .and. jndx(i) <= my_proc*nrows ) then
                  l_row = jndx(i) - (my_proc-1)*nrows
                  col_ind(row_ptr(l_row)+col_ptr(l_row)) =  indx(i)
                  l_cval(row_ptr(l_row)+col_ptr(l_row))  = -cval(i)
                  col_ptr(l_row) = col_ptr(l_row) + 1
               end if
            end do
         elseif ( symm == 'hermitian' ) then
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_row = indx(i) - (my_proc-1)*nrows
                  col_ind(row_ptr(l_row)+col_ptr(l_row)) = jndx(i)
                  l_cval(row_ptr(l_row)+col_ptr(l_row))  = cval(i)
                  col_ptr(l_row) = col_ptr(l_row) + 1
               end if
               if ( indx(i) .ne. jndx(i) ) then
                  if ( jndx(i) > ( my_proc-1)*nrows .and. jndx(i) <= my_proc*nrows ) then
                     l_row = jndx(i) - (my_proc-1)*nrows
                     col_ind(row_ptr(l_row)+col_ptr(l_row)) = indx(i)
                     l_cval(row_ptr(l_row)+col_ptr(l_row))  = conjg(cval(i))
                     col_ptr(l_row) = col_ptr(l_row) + 1
                  end if
               end if
            end do
         end if
      end if
   end subroutine mm_crs

   subroutine mm_coo( nrows, ncols, nnz, ival, indx, jndx, rval, cval, &
                      my_proc, n_procs, symm, field, &
                      l_nnz, l_indx, l_jndx, l_rval, l_ival, l_cval )
!
! Subroutine to convert matrix-market sparse format to
! coordinate format. This is done only for the nrows assigned to
! this processor. Variables starting with l_ are output and unique for the local
! processor.
!
      integer, intent(in)                        :: nrows, ncols, nnz
      integer, intent(in)                        :: ival(:), indx(:), jndx(:)
      double precision, intent(in)               :: rval(:)
      complex, intent(in)                        :: cval(:)
      integer, intent(in)                        :: my_proc, n_procs
      character(len=7), intent(in)               :: field
      character(len=19), intent(in)              :: symm
!
      integer, intent(out)                       :: l_nnz
      integer, allocatable, intent(out)          :: l_indx(:), l_jndx(:)
      real(kind=rp), allocatable, intent(out)    :: l_rval(:)
      integer, allocatable, intent(out)          :: l_ival(:)
      complex(kind=cp), allocatable, intent(out) :: l_cval(:)
!
      integer                                    :: i, j, k

      if ( symm == 'general' ) then
         k = 1
         do i = 1, nnz
            if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) k = k + 1
         end do
      else
         k = 1
         do i = 1, nnz
            if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) k = k + 1
            if ( indx(i) .ne. jndx(i) ) then
               if ( jndx(i) > ( my_proc-1)*nrows .and. jndx(i) <= my_proc*nrows ) k = k + 1
            end if
         end do
      end if
      l_nnz = k-1
      allocate( l_indx(l_nnz), l_jndx(l_nnz) )

! Determine the numbering of the local equations:
      if ( symm == 'general' ) then
         k = 1
         do i = 1, nnz
            if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
               l_indx(k) = indx(i) - (my_proc-1)*nrows
               l_jndx(k) = jndx(i)
               k = k + 1
            end if
         end do
      else
         k = 1
         do i = 1, nnz
            if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
               l_indx(k) = indx(i) - (my_proc-1)*nrows
               l_jndx(k) = jndx(i)
               k = k + 1
            end if
            if ( indx(i) .ne. jndx(i) ) then
               if ( jndx(i) > ( my_proc-1)*nrows .and. jndx(i) <= my_proc*nrows ) then
                  l_indx(k) = jndx(i) - (my_proc-1)*nrows
                  l_jndx(k) = indx(i)
                  k = k + 1
               end if
            end if
         end do
      end if
!
      if ( field == 'real' ) then
         allocate( l_rval(l_nnz) )
         if ( symm == 'general' ) then
            k = 1
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_rval(k) = rval(i)
                  k = k + 1
               end if
            end do
         elseif ( symm == 'symmetric' ) then
            k = 1
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_rval(k) = rval(i)
                  k = k + 1
               end if
               if ( indx(i) .ne. jndx(i) ) then
                  if ( jndx(i) > ( my_proc-1)*nrows .and. jndx(i) <= my_proc*nrows ) then
                     l_rval(k) = rval(i)
                     k = k + 1
                  end if
               end if
            end do
         elseif ( symm == 'skew-symmetric' ) then
            k = 1
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_rval(k) = rval(i)
                  k = k + 1
               end if
               if ( jndx(i) > ( my_proc-1)*nrows .and. jndx(i) <= my_proc*nrows ) then
                  l_rval(k) =  -rval(i)
                  k = k + 1
               end if
            end do
         end if
      elseif ( field == 'integer' ) then
         allocate( l_rval(l_nnz) )
         if ( symm == 'general' ) then
            k = 1
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_ival(k) = ival(i)
                  k = k + 1
               end if
            end do
         elseif ( symm == 'symmetric' ) then
            k = 1
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_ival(k) = ival(i)
                  k = k + 1
               end if
               if ( indx(i) .ne. jndx(i) ) then
                  if ( jndx(i) > ( my_proc-1)*nrows .and. jndx(i) <= my_proc*nrows ) then
                     l_ival(k) = ival(i)
                     k = k + 1
                  end if
               end if
            end do
         elseif ( symm == 'skew-symmetric' ) then
            k = 1
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_ival(k) = ival(i)
                  k = k + 1
               end if
               if ( jndx(i) > ( my_proc-1)*nrows .and. jndx(i) <= my_proc*nrows ) then
                  l_ival(k) =  -ival(i)
                  k = k + 1
               end if
            end do
         end if
      elseif ( field == 'complex' ) then
         allocate( l_cval(l_nnz) )
         if ( symm == 'general' ) then
            k = 1
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_cval(k) = cval(i)
                  k = k + 1
               end if
            end do
         elseif ( symm == 'symmetric' ) then
            k = 1
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_cval(k) = cval(i)
                  k = k + 1
               end if
               if ( indx(i) .ne. jndx(i) ) then
                  if ( jndx(i) > ( my_proc-1)*nrows .and. jndx(i) <= my_proc*nrows ) then
                     l_cval(k) = cval(i)
                     k = k + 1
                  end if
               end if
            end do
         elseif ( symm == 'skew-symmetric' ) then
            k = 1
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_cval(k) = cval(i)
                  k = k + 1
               end if
               if ( jndx(i) > ( my_proc-1)*nrows .and. jndx(i) <= my_proc*nrows ) then
                  l_cval(k) =  -cval(i)
                  k = k + 1
               end if
            end do
         elseif ( symm == 'hermitian' ) then
            k = 1
            do i = 1, nnz
               if ( indx(i) > ( my_proc-1)*nrows .and. indx(i) <= my_proc*nrows ) then
                  l_cval(k) = cval(i)
                  k = k + 1
               end if
               if ( indx(i) .ne. jndx(i) ) then
                  if ( jndx(i) > ( my_proc-1)*nrows .and. jndx(i) <= my_proc*nrows ) then
                     l_cval(k) = conjg(cval(i))
                     k = k + 1
                  end if
               end if
            end do
         end if
      end if
   end subroutine mm_coo

end module mm_module

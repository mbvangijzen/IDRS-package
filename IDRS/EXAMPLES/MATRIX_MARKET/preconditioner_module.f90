!
! The preconditioner module contains routines to compute some standard
! preconditioners for use by IDRS. These precondioners should be 
! passed as argument in the IDRS_DRIVER routines.
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2026 Martin van Gijzen
!

module preconditioner_module

   use precision_module
   use user_module

   implicit none

   contains

   function jacobi( A ) result(P)
!
! Compute Jacobi preconditioner 
!
      type(user_matrix), intent(in)     :: A
      type(user_matrix)                 :: P

      P%matrix_format = 4
      P%half_band = 0
      P%nrows = A%nrows
      P%ncols = A%ncols

      if ( allocated(A%cval) .or. allocated(A%cmat) .or. allocated(A%cband) ) then 
! A is complex matrix 
         allocate(P%cband(P%nrows,-P%half_band:P%half_band))
         P%cband(:,0) = complex_diagonal(A)
      else
         allocate(P%rband(P%nrows,-P%half_band:P%half_band))
         P%rband(:,0) = real_diagonal(A)
      end if
   end function jacobi

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

end module preconditioner_module

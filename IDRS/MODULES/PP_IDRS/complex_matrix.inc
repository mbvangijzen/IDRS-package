!
! Include file to store the matrices in the matrix type 
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2025 Martin van Gijzen
!

!
! Put stiffness matrix in matrix type
   A%K = K

!
! Put mass matrix in matrix type (if present)
   if ( present(M) ) then
      A%mass_matrix = .true.
      A%M = M
   end if

!
! Put damping matrix in matrix type (if present)
   if ( present(C) ) then
      A%damping_matrix = .true.
      A%C = C
   end if

   if ( present(P) ) then
!
! Shift-and-invert?
      A%P = P
      if ( present(shift) ) shift_invert = .true.
      A%shift_invert = shift_invert
   end if

! Multishift problem?
   A%multishift = multishift

! Do we have a block system?
   A%block_system = ( multishift .and. present(C) )

! Shift to matrix:
   if ( present(seed) ) then
      A%cseed = seed
   else
      A%cseed = 0._cp
   end if

! Shift of SI preconditioner:
   if ( present(shift) ) then
      A%cshift = shift
   else
      A%cshift = 0._cp
   end if
!
! Remove preconditioner (if it was included)
   A%preconditioner = 0
   A%degree = 0
   if ( allocated(A%calpha_p) ) deallocate(A%calpha_p)
   if ( allocated(A%calpha_s) ) deallocate(A%calpha_s)


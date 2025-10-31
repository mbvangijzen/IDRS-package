!
! Include file to store the matrices in the matrix type 
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2025 Martin van Gijzen
!
! Put stiffness matrix in matrix type (if present)
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
      A%P = P
!
! Shift-and-invert?
      if ( present(shift) ) shift_invert = .true.
      A%shift_invert = shift_invert
   end if


! Multishift problem?
   A%multishift = multishift

! Do we have a block system?
   A%block_system = ( multishift .and. present(C) )

! Seed system
   if ( present(seed) ) then
      A%rseed = seed
   else
      A%rseed = 0._rp
   end if

! Shift to matrix of SI preconditioner:
   if ( present(shift) ) then
      A%rshift = shift
   else
      A%rshift = 0._rp
   end if
!
! Remove preconditioner (if it was included)
   A%preconditioner = 0
   A%degree = 0
   if ( allocated(A%ralpha_p) ) deallocate(A%ralpha_p)
   if ( allocated(A%ralpha_s) ) deallocate(A%ralpha_s)

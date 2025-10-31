program toeplitz

   write(*,'(a,/,a,/,a,/,a,/,a,/,a,/,a,/,a,/,a,/)') &
      'This programme illustrates IDR(s) using three different test problems:', &
      '   - a complex test problem with a system matrix that is Toeplitz,', &
      '   - a real 1D convection-diffusion problem, ', &
      'The three IDR(s) variants are tested:', & 
      '   - IDRS: the most efficient method, uses bi-orthogonalisation;', &
      '   - BiCGSTAB: IDR(1) with BiCGSTAB parameter choice;', &
      '   - QMRIDR: a robust method based on modified Gram-Schmidt orthogonalisation.'

    call complex_toeplitz
    call real_toeplitz

    write(*,'(a)') &
      '===================================================================='

end program toeplitz

subroutine complex_toeplitz

   use precision_module
   use idrs_module
   use matrix_module

   implicit none 
   integer, parameter            :: maxit=4000 
   integer                       :: flag            
   integer                       :: iter          
   integer                       :: s            
  
   integer, parameter            :: N = 200
   type(matrix)                  :: A, M
   complex(kind=cp)              :: x(N), b(N)
   real(kind=rp)                 :: resvec(maxit+1)
   real(kind=rp)                 :: tol, relres
   character(len=8)              :: method
   integer                       :: variant
   integer                       :: i

! Generate the test problems
   A%cind = (/-1, 0, 2, 3/)
   A%cval = (/(0.,3.6_cp), (4._cp,0.), (1._cp,0.), (0.7_cp,0.)/)

   b = (0.,1._cp)

! Solve accuracy test
   tol = 1e-12_rp
    write(*,'(a)') &
      '===================================================================='
   write(*,'(a,/,a,/)') &
      'Accuracy test problem, see TOMS paper (2011), example 2', &
      'The termination criterion is ||r||/||b|| < 1e-12. ' 

   method = 'IDRS    '
   variant = 1
   write(*,'(5x,a,/,5x,a)') &
      '   Variant   |    s    | Iterations |  ||r||/||b||  |  Flag  ',&
      '--------------------------------------------------------------'
   do i = 1,7
      s = 2**(i-1)
      x = idrs( A, b, s, M, tol, maxit, variant, flag, relres, iter )
      write(*,'(8x,a8, a,i3,a,i4,a,e10.3,a,i1)') method, &
         '  |  ', s,'    |   ', iter,'     |   ', relres,'   |    ',flag
   end do
   
   method = 'BICGSTAB'
   variant = 2
   s = 1
   x = idrs( A, b, s, M, tol, maxit, variant, flag, relres, iter, resvec=resvec )
   write(*,'(8x,a8, a,i3,a,i4,a,e10.3,a,i1)') method, &
      '  |  ', s,'    |   ', iter,'     |   ', relres,'   |    ',flag
!
! UNCOMMEND THE CODE BELOW TO COMPARE THE CONVERGENCE WITH MATLAB'S BICGSTAB
!  write(*,'(a)') 'The norms of the residuals of the first ten Bi-CGSTAB iterations are:'
!  do i = 1, 10
!     write(*,'(18x,a,i1,a,e21.15)') '||r_',i-1,'|| = ', resvec(i)
!  end do
!  write(*,*)

   method = 'QMRIDR  '
   do i = 1,7
      s = 2**(i-1)
      x = qmridr( A, b, s, M, tol, maxit, flag, relres, iter )
      write(*,'(8x,a8, a,i3,a,i4,a,e10.3,a,i1)') method, &
         '  |  ', s,'    |   ', iter,'     |   ', relres,'   |    ',flag
   end do
   write(*,'(5x,a)') &
      '--------------------------------------------------------------'
   write(*,'(18x,a)') 'Results for the complex Toeplitz problem'
   write(*,*)

end subroutine complex_toeplitz

subroutine real_toeplitz

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 1D convection diffusion test problem to test finite termination
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   use precision_module
   use matrix_module
   use idrs_module

   implicit none 
   integer                       :: flag            
   integer, parameter            :: maxit = 4000
   integer                       :: iter          
   integer                       :: s            
  
   real(kind=rp), parameter      :: Pe = 0.5
   integer,       parameter      :: N = 60

   type(matrix)                  :: A, M
   real(kind=rp)                 :: resvec(maxit+1)
   real(kind=rp)                 :: x(N), b(N)
   real(kind=rp)                 :: tol, relres
   character(len=8)              :: method
   integer                       :: variant
   integer                       :: i

! Generate the test problems
   A%rind = (/-1, 0, 1/)
   A%rval = (/-1._rp-Pe, 2._rp,  -1._rp+Pe/)

   b = 0.
   b(1) = 1._rp+Pe
   b(N) = 1._rp-Pe

! Solve 1D convection-diffusion test
   write(*,'(a)') &
      '===================================================================='
   tol = 1.e-8_rp
   write(*,'(a,/,a,/,a,/,a,/)') &
      '1D convection-diffusion problem, N=60 gridpoints.', & 
      'This problem tests the finite termination of IDR(s).', &
      'IDR(s) should terminate in at most N+N/s iterations.', &
      'The termination criterion is ||r||/||b|| < 1e-8. ' 

   write(*,'(5x,a,/,5x,a)') &
      '   Variant   |    s    | Iterations |  ||r||/||b||  |  Flag  ',&
      '--------------------------------------------------------------'
   method = 'IDRS    '
   variant = 1
   do i = 1, 4
      s = 2**(i-1)
      x = idrs( A, b, s, M, tol, maxit, variant, flag, relres, iter )
      write(*,'(8x,a8, a,i3,a,i4,a,e10.3,a,i1)') method, &
         '  |  ', s,'    |   ', iter,'     |   ', relres,'   |    ',flag
   end do

   method = 'BICGSTAB'
   s = 1
   variant = 2
   x = idrs( A, b, s, M, tol, maxit, variant, flag, relres, iter, resvec=resvec )
   write(*,'(8x,a8, a,i3,a,i4,a,e10.3,a,i1)') method, &
      '  |  ', s,'    |   ', iter,'     |   ', relres,'   |    ',flag

   method = 'QMRIDR  '
   do i = 1, 4
      s = 2**(i-1)
      x = qmridr( A, b, s, M, tol, maxit, flag, relres, iter, resvec=resvec )
      write(*,'(8x,a8, a,i3,a,i4,a,e10.3,a,i1)') method, &
         '  |  ', s,'    |   ', iter,'     |   ', relres,'   |    ',flag
   end do
   write(*,'(5x,a)') &
      '--------------------------------------------------------------'
   write(*,'(15x,a)') 'Results for the 1D convection-diffusion problem'
   write(*,*)
!
!  write(*,'(a)') 'The norms of the residuals of the last ten iterations are:'
!  do i = iter-10, iter
!     write(*,'(18x,a,i1,a,e21.15)') '||r_',i-1,'|| = ', resvec(i+1)
!  end do
!  write(*,*)

end subroutine real_toeplitz

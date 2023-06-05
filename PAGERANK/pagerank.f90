program pagerank

!
! This programme illustrates the use of the idrs-module. A random PageRank problem 
! is generated and solved with the different idrs-algorithms: IDRS and QMRIDR 
! for the single shift (= damping parameter) problem and MSIDRS and MSQMRIDR 
! for the multishift problem.
!
! For simplicity, the PageRank matrix is stored in dense form. From a computation
! point of view this is not the best way, but done here to focus only on the
! proper use of the idrs-module.
!
! The PageRank linear systems are quite easy to solve and therefore not very suited
! to illustrate the differences in performance of the variants: all methods
! converge in a number of iterations that is close to optimal (= the number of
! full GMRES iterations).
!
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2023 Martin van Gijzen
!

   use matrix_module  ! Always needed is a user module that defines the 
                      ! matrix-vector multiplication and the preconditioning operation
                      ! This module shoudl also define the parameters
                      ! rp (kind of real precision) and cp (kind of complex precision)
   use preconditioner_module
   use idrs_module

   implicit none 
   integer, parameter            :: maxit=4000 
   integer                       :: flag            
   integer                       :: iter          
   integer                       :: s            
  
   integer, parameter            :: N = 10000        ! Number of nodes
   integer, parameter            :: N_out = 3        ! Average number of outlinks per node
   integer, parameter            :: N_dangling = N/10! Number of dangling nodes (about 10%)
   real(kind=rp), parameter      :: p = 0.99         ! Teleportation chance
   integer, parameter            :: N_sigma = 4      ! Number of damping factor
   real(kind=rp), parameter      :: sigma(N_sigma) = [0.85, 0.9, 0.95, 0.99]
   real(kind=rp)                 :: shift(N_sigma)
   real(kind=rp)                 :: rdangling(N_dangling)
   integer                       :: dangling(N_dangling)
   real(kind=rp)                 :: G(N,N)
   type(matrix)                  :: A
   type(preconditioner)          :: M1
   real(kind=rp)                 :: x(N), b(N), x_ms(N,N_sigma)
   real(kind=rp)                 :: tol, relres, relres_ms(N_sigma)
   integer                       :: variant
   integer                       :: i

   real                          :: som, t, t0

! Generate a simple test problem: a PageRank matrix

! Generate random matrix with numbers between 0 and 1
   call RANDOM_SEED
   call RANDOM_NUMBER(G)

! Turn this matrix into a binary connectivity matrix with on-avarage N_out outlinks
   where ( G > real(N-N_out)/real(N) )
      G = 1.
   elsewhere 
      G = 0.
   end where

! Generate dangling nodes
   call RANDOM_NUMBER(rdangling)
   dangling = ceiling(N*rdangling)

! Put the dangling nodes in the matrix (colums are 1)
   do i = 1, N_dangling
      G(:,dangling(i)) = 0.
   end do

   do i = 1, N
! Remove self-references
      G(i,i) = 0.
! Scale the columns to sum 1
      som = sum(G(:,i))
      if ( som > 0 ) then
         G(:,i) = G(:,i)/som
      end if
   end do

! Generate the system matrix A = G-1/pI
   allocate(A%Ar(N,N)) 
   A%Ar = G
   do i = 1, N
      A%Ar(i,i) = A%Ar(i,i) - 1./p
   end do
   
! (pG + (1-p)ee^T)x = x
! (pGx + (1-p)e = x
! (pG -I)x = -(1-p)e
! (G - (1/p))x = (1-1/p) e

! Right-hand-side vector
   b = 1. 

! Solve accuracy test
   tol = 1d-9
   write(*,'(a)') &
      '===================================================================================='
   write(*,'(15x,a,i5,/)') &
      'Randomly generated PageRank problem, size is ', N

   write(*,'(5x,a,/,5x,a)') &
      '   Method      |    s    | Iterations |  ||r||/||b||  |  Flag  |    CPU time   ',&
      '-------------------------------------------------------------------------------'
   
   variant = 1
   do i = 1,4
      s = 2**(i-1)
      call cpu_time( t0 )
      x = idrs( A, b, M1, s, tol, maxit, variant, flag, relres, iter )
      call cpu_time( t )
      write(*,'(8x,a10, a,i3,a,i4,a,e9.3,a,i1,a,e9.3)') 'IDRS', &
         '  |  ', s,'    |   ', iter,'     |   ', relres,'   |    ',flag,'   |    ', t-t0
   end do
   
   variant = 2
   s = 1
   call cpu_time( t0 )
   x = idrs( A, b, M1, s, tol, maxit, variant, flag, relres, iter )
   call cpu_time( t )
   write(*,'(8x,a10, a,i3,a,i4,a,e9.3,a,i1,a,e9.3)') 'BICGSTAB', &
      '  |  ', s,'    |   ', iter,'     |   ', relres,'   |    ',flag,'   |    ', t-t0

   do i = 1,4
      s = 2**(i-1)
      call cpu_time( t0 )
      x = qmridr( A, b, M1, s, tol, maxit, flag, relres, iter )
      call cpu_time( t )
      write(*,'(8x,a10, a,i3,a,i4,a,e9.3,a,i1,a,e9.3)') 'QMRIDR', &
         '  |  ', s,'    |   ', iter,'     |   ', relres,'   |    ',flag,'   |    ', t-t0
   end do

! GMRES
   s = 100
   call cpu_time( t0 )
   x = qmridr( A, b, M1, s, tol, maxit, flag, relres, iter )
   call cpu_time( t )
   write(*,'(8x,a10, a,i3,a,i4,a,e9.3,a,i1,a,e9.3)') 'GMRES', &
         '  |  ', s,'    |   ', iter,'     |   ', relres,'   |    ',flag,'   |    ', t-t0
   write(*,'(5x,a)') &
      '-------------------------------------------------------------------------------'
   write(*,'(23x,a)') 'Results for the PageRank test problem'
   write(*,*)

   write(*,'(a)') &
      '====================================================================================='
   write(*,'(a,i5,a,i2,/)') &
      'Randomly generated multishift PageRank problem, size is ', N, '. Number of shifts is ', N_sigma
   do i = 1, N_sigma
      shift(i) = 1./sigma(i) - 1./p
   end do

   write(*,'(5x,a,/,5x,a)') &
      '   Method      |    s    | Iterations |  ||r||/||b||  |  Flag  |    CPU time   ',&
      '-------------------------------------------------------------------------------'

   variant = 1
   do i = 1,4
      s = 2**(i-1)
      call cpu_time( t0 )
      x_ms = msidrs( A, b, shift, M1, s, tol, maxit, variant, flag, relres_ms, iter )
      relres = maxval( relres_ms)
      call cpu_time( t )
      write(*,'(8x,a10, a,i3,a,i4,a,e9.3,a,i1,a,e9.3)') 'MSIDRS', &
         '  |  ', s,'    |   ', iter,'     |   ', relres,'   |    ',flag,'   |    ', t-t0
   end do
   
   variant = 2
   s = 1
   call cpu_time( t0 )
   x_ms = msidrs( A, b, shift, M1, s, tol, maxit, variant, flag, relres_ms, iter )
   relres = maxval( relres_ms)
   call cpu_time( t )
   write(*,'(8x,a10, a,i3,a,i4,a,e9.3,a,i1,a,e9.3)') 'MSBICGSTAB', &
      '  |  ', s,'    |   ', iter,'     |   ', relres,'   |    ',flag,'   |    ', t-t0

   do i = 1,4
      s = 2**(i-1)
      call cpu_time( t0 )
      x_ms = msqmridr( A, b, shift, M1, s, tol, maxit, flag, relres_ms, iter )
      relres = maxval( relres_ms)
      call cpu_time( t )
      write(*,'(8x,a10, a,i3,a,i4,a,e9.3,a,i1,a,e9.3)') 'MSQMRIDR', &
         '  |  ', s,'    |   ', iter,'     |   ', relres,'   |    ',flag,'   |    ', t-t0
   end do
   write(*,'(5x,a)') &
      '-------------------------------------------------------------------------------'
   write(*,'(18x,a)') 'Results for the multishift PageRank test problem'
   write(*,*)

end program pagerank



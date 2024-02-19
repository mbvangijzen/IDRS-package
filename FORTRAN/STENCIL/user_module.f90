!
! This module contains problem dependent routines, to be defined by
! the user. The present module contains routines that help make
! a polynomial or Jacobi preconditioner
!
! In particular, it contains the following routines and functions:
!
!    real_diagonal: assemble the main diagonal of a real matrix
!    function call: Dr = real_diagonal( A )
!    data types:    Dr: real dimension size neq of kind rp
!                   A: matrix (user defined data structure)
!
!    real_diagonal: assemble the main diagonal of a real matrix
!    function call: Dr = complex_diagonal( A )
!    data types:    Dr: complex dimension size neq of kind cp
!                   A: matrix (user defined data structure)
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2023 Martin van Gijzen
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
!%%%%%%
!
module user_module

   use interface_module
   use precision_module

   implicit none

! Define the user matrix type

   type user_matrix
      real(kind=rp)   :: c, e, w, n, s, b, f  ! Stencil
      real(kind=rp)   :: rk                   ! Real Robin parameter
      complex(kind=cp):: ck                   ! Complex Robin parameter
      real(kind=rp)   :: Lx, Ly, Lz           ! Size of the domain
      real(kind=rp)   :: hx, hy, hz           ! Grid sizes
      integer         :: Nx, Ny, Nz, Nn       ! Points per directions per subdomain
      integer         :: Sx, Sy, Sz           ! Number of subdomains in each direction
      integer         :: me(3)                ! Subdomain numbers of this processor
      integer         :: bc(6)                ! Type of boundary condition 
      logical         :: multishift           ! Multishift problem?
   end type user_matrix

! Overload * to define the matrix-vector multiplication using the matrix type

   INTERFACE OPERATOR(*)
      module procedure ruser_mv, cuser_mv
   END INTERFACE

   INTERFACE GERSHGORIN
      module procedure rgershgorin, cgershgorin
   END INTERFACE

   public :: GERSHGORIN, CDR_MATRIX, HELMHOLTZ_MATRIX, MASS_MATRIX, MODEL_SOLUTION

   contains

   function ruser_mv( A, v ) result(w)

      type(user_matrix), intent(in)  :: A
      real(kind=rp), intent(in)      :: v(A%Nn) 
      real(kind=rp)                  :: w(A%Nn)
      real(kind=rp), allocatable,save:: v_grid(:,:,:)[:,:,:] 
      real(kind=rp)                  :: w_grid(1:A%Nx,1:A%Ny,1:A%Nz)
      integer                        :: i,j,k,l, me(3), bc(6)
      
! Allocate co-array. Note that this performs a global synchronization.
      if ( .not. allocated(v_grid) ) allocate(v_grid(0:A%Nx+1,0:A%Ny+1,0:A%Nz+1)[A%Sx,A%Sy,*])
!
! map v to interior points of v_grid
      v_grid(1:A%Nx,1:A%Ny,1:A%Nz) = reshape( v, (/A%Nx, A%Ny, A%Nz/) )
      me = A%me
      bc = A%bc
!
! Fill the boundary points of v_grid
      sync all

! West:
      if ( bc(1) == 0 ) then
         v_grid(0,1:A%Ny,1:A%Nz) = v_grid(A%Nx,1:A%Ny,1:A%Nz)[me(1)-1,me(2),me(3)]   ! Internal boundary
      elseif ( bc(1) == 1 ) then
         v_grid(0,1:A%Ny,1:A%Nz) = 0.                                                ! Dirichlet boundary
      elseif ( bc(1) == 2 ) then
         v_grid(0,1:A%Ny,1:A%Nz) = v_grid(2,1:A%Ny,1:A%Nz)                           ! Neumann boundary
      elseif ( bc(1) == 3 ) then
         v_grid(0,1:A%Ny,1:A%Nz) = v_grid(2,1:A%Ny,1:A%Nz)                           ! Robin boundary
      end if

! East:
      if ( bc(2) == 0 ) then
         v_grid(A%Nx+1,1:A%Ny,1:A%Nz) = v_grid(1,1:A%Ny,1:A%Nz)[me(1)+1,me(2),me(3)] ! Internal boundary
      elseif ( bc(2) == 1 ) then
         v_grid(A%Nx+1,1:A%Ny,1:A%Nz) = 0.                                           ! Dirichlet boundary
      elseif ( bc(2) == 2 ) then
         v_grid(A%Nx+1,1:A%Ny,1:A%Nz) = v_grid(A%Nx-1,1:A%Ny,1:A%Nz)                 ! Neumann boundary
      elseif ( bc(2) == 3 ) then
         v_grid(A%Nx+1,1:A%Ny,1:A%Nz) = v_grid(A%Nx-1,1:A%Ny,1:A%Nz)                 ! Robin boundary
      end if

! South:
      if ( bc(3) == 0 ) then
         v_grid(1:A%Nx,0,1:A%Nz) = v_grid(1:A%Nx,A%Ny,1:A%Nz)[me(1),me(2)-1,me(3)]   ! Internal boundary
      elseif ( bc(3) == 1 ) then
         v_grid(1:A%Nx,0,1:A%Nz) = 0.                                                ! Dirichlet boundary
      elseif ( bc(3) == 2 ) then
         v_grid(1:A%Nx,0,1:A%Nz) = v_grid(1:A%Nx,2,1:A%Nz)                           ! Neumann boundary
      elseif ( bc(3) == 3 ) then
         v_grid(1:A%Nx,0,1:A%Nz) = v_grid(1:A%Nx,2,1:A%Nz)                           ! Robin boundary
      end if

! North:
      if ( bc(4) == 0 ) then
         v_grid(1:A%Nx,A%Ny+1,1:A%Nz) = v_grid(1:A%Nx,1,1:A%Nz)[me(1),me(2)+1,me(3)] ! Internal boundary
      elseif ( bc(4) == 1 ) then
         v_grid(1:A%Nx,A%Ny+1,1:A%Nz) = 0.                                           ! Dirichlet boundary
      elseif ( bc(4) == 2 ) then
         v_grid(1:A%Nx,A%Ny+1,1:A%Nz) = v_grid(1:A%Nx,A%Ny-1,1:A%Nz)                 ! Neumann boundary
      elseif ( bc(4) == 3 ) then
         v_grid(1:A%Nx,A%Ny+1,1:A%Nz) = v_grid(1:A%Nx,A%Ny-1,1:A%Nz)                 ! Robin boundary
      end if

! Front:
      if ( bc(5) == 0 ) then
         v_grid(1:A%Nx,1:A%Ny,0) = v_grid(1:A%Nx,1:A%Ny,A%Nz)[me(1),me(2),me(3)-1]   ! Internal boundary
      elseif ( bc(5) == 1 ) then
         v_grid(1:A%Nx,1:A%Ny,0) = 0.                                                ! Dirichlet boundary
      elseif ( bc(5) == 2 ) then
         v_grid(1:A%Nx,1:A%Ny,0) = v_grid(1:A%Nx,1:A%Ny,2)                           ! Neumann boundary
      elseif ( bc(5) == 3 ) then
         v_grid(1:A%Nx,1:A%Ny,0) = v_grid(1:A%Nx,1:A%Ny,2)                           ! Robin boundary
      end if

! Back:
      if ( bc(6) == 0 ) then
         v_grid(1:A%Nx,1:A%Ny,A%Nz+1) = v_grid(1:A%Nx,1:A%Ny,1)[me(1),me(2),me(3)+1] ! Internal boundary
      elseif ( bc(6) == 1 ) then
         v_grid(1:A%Nx,1:A%Ny,A%Nz+1) = 0.                                           ! Dirichlet boundary
      elseif ( bc(6) == 2 ) then
         v_grid(1:A%Nx,1:A%Ny,A%Nz+1) = v_grid(1:A%Nx,1:A%Ny,A%Nz-1)                 ! Neumann boundary
      elseif ( bc(6) == 3 ) then
         v_grid(1:A%Nx,1:A%Ny,A%Nz+1) = v_grid(1:A%Nx,1:A%Ny,A%Nz-1)                 ! Robin boundary
      end if
      sync all

! Perform the multiplication:
      do concurrent ( i=1:A%Nx, j=1:A%Ny, k=1:A%Nz )
         w_grid(i,j,k) =  A%c*v_grid(i,j,k) + &
                          A%w*v_grid(i-1,j,k) + A%e*v_grid(i+1,j,k) + &
                          A%s*v_grid(i,j-1,k) + A%n*v_grid(i,j+1,k) + &
                          A%f*v_grid(i,j,k-1) + A%b*v_grid(i,j,k+1) 
      end do

! Correct for Robin condition:
! West:
      if ( bc(1) == 3 ) then
         w_grid(1,1:A%Ny,1:A%Nz) = w_grid(1,1:A%Ny,1:A%Nz) + &
         v_grid(1,1:A%Ny,1:A%Nz)*(2./A%hx)*A%rk
      end if

! East:
      if ( bc(2) == 3 ) then
         w_grid(A%Nx,1:A%Ny,1:A%Nz) = w_grid(A%Nx,1:A%Ny,1:A%Nz) + &
         v_grid(A%Nx,1:A%Ny,1:A%Nz)*(2./A%hx)*A%rk
      end if

! South:
      if ( bc(3) == 3 ) then
         w_grid(1:A%Nx,1,1:A%Nz) = w_grid(1:A%Nx,1,1:A%Nz) + &
         v_grid(1:A%Nx,1,1:A%Nz)*(2./A%hy)*A%rk
      end if

! North:
      if ( bc(4) == 3 ) then
         w_grid(1:A%Nx,A%Ny,1:A%Nz) = w_grid(1:A%Nx,A%Ny,1:A%Nz) + &
         v_grid(1:A%Nx,A%Ny,1:A%Nz)*(2./A%hy)*A%rk
      end if

! Front:
      if ( bc(5) == 3 ) then
         w_grid(1:A%Nx,1:A%Ny,1) = w_grid(1:A%Nx,1:A%Ny,1) + &
         v_grid(1:A%Nx,1:A%Ny,1)*(2./A%hz)*A%rk
      end if

! Back:
      if ( bc(6) == 3 ) then
         w_grid(1:A%Nx,1:A%Ny,A%Nz) = w_grid(1:A%Nx,1:A%Ny,A%Nz) + &
         v_grid(1:A%Nx,1:A%Ny,A%Nz)*(2./A%hz)*A%rk
      end if

! Map w to vector
      w = reshape(w_grid,(/A%Nn/))

   end function ruser_mv

   function cuser_mv( A, v ) result(w) 

      type(user_matrix), intent(in)     :: A
      complex(kind=cp), intent(in)      :: v(:)
      complex(kind=cp)                  :: w(size(v,1))
      complex(kind=cp), allocatable,save:: v_grid(:,:,:)[:,:,:] 
      complex(kind=cp)                  :: w_grid(1:A%Nx,1:A%Ny,1:A%Nz)
      integer                           :: i,j,k,l, me(3), bc(6)
      
! Allocate co-array. Note that this performs a global synchronization.
      if ( .not. allocated(v_grid) ) allocate(v_grid(0:A%Nx+1,0:A%Ny+1,0:A%Nz+1)[A%Sx,A%Sy,*])

!
! Initialize w
      w = 0.

      if ( A%multishift ) then
! A is block matrix
!
!           A = [-C -K]
!               [ I O ]
!
! Multiplication with bottom row of A: 
         w(A%Nn+1:2*A%Nn) = v(1:A%Nn)                      

! Get - bottom half of v for the multiplication with the Laplacian:
         v_grid(1:A%Nx,1:A%Ny,1:A%Nz) = reshape( -v(A%Nn+1:2*A%Nn), (/A%Nx, A%Ny, A%Nz/) )
      else
         v_grid(1:A%Nx,1:A%Ny,1:A%Nz) = reshape( v, (/A%Nx, A%Ny, A%Nz/) )
      end if
      me = A%me
      bc = A%bc
!
! Get information from boundary nodes and neighbouring domains,
! put it in virtual nodes
      sync all

! West:
      if ( bc(1) == 0 ) then
         v_grid(0,1:A%Ny,1:A%Nz) = v_grid(A%Nx,1:A%Ny,1:A%Nz)[me(1)-1,me(2),me(3)]   ! Internal boundary
      elseif ( bc(1) == 1 ) then
         v_grid(0,1:A%Ny,1:A%Nz) = 0.                                                ! Dirichlet boundary
      elseif ( bc(1) == 2 ) then
         v_grid(0,1:A%Ny,1:A%Nz) = v_grid(2,1:A%Ny,1:A%Nz)                           ! Neumann boundary
      elseif ( bc(1) == 3 ) then
         v_grid(0,1:A%Ny,1:A%Nz) = v_grid(2,1:A%Ny,1:A%Nz)                           ! Sommerfeld boundary
      end if

! East:
      if ( bc(2) == 0 ) then
         v_grid(A%Nx+1,1:A%Ny,1:A%Nz) = v_grid(1,1:A%Ny,1:A%Nz)[me(1)+1,me(2),me(3)] ! Internal boundary
      elseif ( bc(2) == 1 ) then
         v_grid(A%Nx+1,1:A%Ny,1:A%Nz) = 0.                                           ! Dirichlet boundary
      elseif ( bc(2) == 2 ) then
         v_grid(A%Nx+1,1:A%Ny,1:A%Nz) = v_grid(A%Nx-1,1:A%Ny,1:A%Nz)                 ! Neumann boundary
      elseif ( bc(2) == 3 ) then
         v_grid(A%Nx+1,1:A%Ny,1:A%Nz) = v_grid(A%Nx-1,1:A%Ny,1:A%Nz)                 ! Sommerfeld boundary
      end if

! South:
      if ( bc(3) == 0 ) then
         v_grid(1:A%Nx,0,1:A%Nz) = v_grid(1:A%Nx,A%Ny,1:A%Nz)[me(1),me(2)-1,me(3)]   ! Internal boundary
      elseif ( bc(3) == 1 ) then
         v_grid(1:A%Nx,0,1:A%Nz) = 0.                                                ! Dirichlet boundary
      elseif ( bc(3) == 2 ) then
         v_grid(1:A%Nx,0,1:A%Nz) = v_grid(1:A%Nx,2,1:A%Nz)                           ! Neumann boundary
      elseif ( bc(3) == 3 ) then
         v_grid(1:A%Nx,0,1:A%Nz) = v_grid(1:A%Nx,2,1:A%Nz)                           ! Sommerfeld boundary
      end if

! North:
      if ( bc(4) == 0 ) then
         v_grid(1:A%Nx,A%Ny+1,1:A%Nz) = v_grid(1:A%Nx,1,1:A%Nz)[me(1),me(2)+1,me(3)] ! Internal boundary
      elseif ( bc(4) == 1 ) then
         v_grid(1:A%Nx,A%Ny+1,1:A%Nz) = 0.                                           ! Dirichlet boundary
      elseif ( bc(4) == 2 ) then
         v_grid(1:A%Nx,A%Ny+1,1:A%Nz) = v_grid(1:A%Nx,A%Ny-1,1:A%Nz)                 ! Neumann boundary
      elseif ( bc(4) == 3 ) then
         v_grid(1:A%Nx,A%Ny+1,1:A%Nz) = v_grid(1:A%Nx,A%Ny-1,1:A%Nz)                 ! Neumann boundary
      end if

! Front:
      if ( bc(5) == 0 ) then
         v_grid(1:A%Nx,1:A%Ny,0) = v_grid(1:A%Nx,1:A%Ny,A%Nz)[me(1),me(2),me(3)-1]   ! Internal boundary
      elseif ( bc(5) == 1 ) then
         v_grid(1:A%Nx,1:A%Ny,0) = 0.                                                ! Dirichlet boundary
      elseif ( bc(5) == 2 ) then
         v_grid(1:A%Nx,1:A%Ny,0) = v_grid(1:A%Nx,1:A%Ny,2)                           ! Neumann boundary
      elseif ( bc(5) == 3 ) then
         v_grid(1:A%Nx,1:A%Ny,0) = v_grid(1:A%Nx,1:A%Ny,2)                           ! Sommerfeld boundary
      end if

! Back:
      if ( bc(6) == 0 ) then
         v_grid(1:A%Nx,1:A%Ny,A%Nz+1) = v_grid(1:A%Nx,1:A%Ny,1)[me(1),me(2),me(3)+1] ! Internal boundary
      elseif ( bc(6) == 1 ) then
         v_grid(1:A%Nx,1:A%Ny,A%Nz+1) = 0.                                           ! Dirichlet boundary
      elseif ( bc(6) == 2 ) then
         v_grid(1:A%Nx,1:A%Ny,A%Nz+1) = v_grid(1:A%Nx,1:A%Ny,A%Nz-1)                 ! Neumann boundary
      elseif ( bc(6) == 3 ) then
         v_grid(1:A%Nx,1:A%Ny,A%Nz+1) = v_grid(1:A%Nx,1:A%Ny,A%Nz-1)                 ! Sommerfeld boundary
      end if
      sync all

! Perform the multiplication:
      do concurrent ( i=1:A%Nx, j=1:A%Ny, k=1:A%Nz )
         w_grid(i,j,k) =  A%c*v_grid(i,j,k) + &
                          A%w*v_grid(i-1,j,k) + A%e*v_grid(i+1,j,k) + &
                          A%s*v_grid(i,j-1,k) + A%n*v_grid(i,j+1,k) + &
                          A%f*v_grid(i,j,k-1) + A%b*v_grid(i,j,k+1) 
      end do

      if ( A%multishift ) then
! Update the upper part of the solution vector
         w(1:A%Nn) = w(1:A%Nn) + reshape(w_grid, (/A%Nn/) )
         w_grid = 0._cp

! Now get - the upper half of v for the multiplication with the damping matrix:
         v_grid(1:A%Nx,1:A%Ny,1:A%Nz) = reshape( -v(1:A%Nn), (/A%Nx, A%Ny, A%Nz/) )
      end if

! Correct for Sommerfeld condition:
      if ( bc(1) == 3 ) then
         w_grid(1,1:A%Ny,1:A%Nz) = w_grid(1,1:A%Ny,1:A%Nz) + &
         v_grid(1,1:A%Ny,1:A%Nz)*(2./A%hx)*A%ck
      end if

! East:
      if ( bc(2) == 3 ) then
         w_grid(A%Nx,1:A%Ny,1:A%Nz) = w_grid(A%Nx,1:A%Ny,1:A%Nz) + &
         v_grid(A%Nx,1:A%Ny,1:A%Nz)*(2./A%hx)*A%ck
      end if

! South:
      if ( bc(3) == 3 ) then
         w_grid(1:A%Nx,1,1:A%Nz) = w_grid(1:A%Nx,1,1:A%Nz) + &
         v_grid(1:A%Nx,1,1:A%Nz)*(2./A%hy)*A%ck
      end if

! North:
      if ( bc(4) == 3 ) then
         w_grid(1:A%Nx,A%Ny,1:A%Nz) = w_grid(1:A%Nx,A%Ny,1:A%Nz) + &
         v_grid(1:A%Nx,A%Ny,1:A%Nz)*(2./A%hy)*A%ck
      end if

! Front:
      if ( bc(5) == 3 ) then
         w_grid(1:A%Nx,1:A%Ny,1) = w_grid(1:A%Nx,1:A%Ny,1) + &
         v_grid(1:A%Nx,1:A%Ny,1)*(2./A%hz)*A%ck
      end if

! Back:
      if ( bc(6) == 3 ) then
         w_grid(1:A%Nx,1:A%Ny,A%Nz) = w_grid(1:A%Nx,1:A%Ny,A%Nz) + &
         v_grid(1:A%Nx,1:A%Ny,A%Nz)*(2./A%hz)*A%ck
      end if
      w(1:A%Nn) = w(1:A%Nn) + reshape(w_grid,(/A%Nn/))

   end function cuser_mv

   subroutine rgershgorin( A, center, radius )
!
! Compute center and radius of circel around spectrum
!
      type(user_matrix), intent(in) :: A
      real(kind=rp), intent(out)    :: center, radius
      center = A%c+A%rk
      radius = abs(A%n) + abs(A%s) + abs(A%e) + abs(A%w) + abs(A%f) +abs(A%b) + abs(A%rk)

   end subroutine rgershgorin

   subroutine cgershgorin( A, center, radius )
!
! Compute center and radius of circel around spectrum
!
      type(user_matrix), intent(in) :: A
      complex(kind=cp), intent(out) :: center
      real(kind=rp), intent(out)    :: radius

      if ( A%multishift ) then
         center = -A%ck
         radius = abs(A%n) + abs(A%s) + abs(A%e) + abs(A%w) + abs(A%f) +abs(A%b) + abs(A%c)
      else
         center = A%c+A%ck
         radius = abs(A%n) + abs(A%s) + abs(A%e) + abs(A%w) + abs(A%f) +abs(A%b)+abs(A%ck)
      end if

   end subroutine cgershgorin

   subroutine get_subdomains( A, grid )
! 
! Determine problem size, subdomain devision and number of unknowns per processor
! The subdomain division is based on the number of processors
! Not for every number of processors a subdomain division is given
! At this moment, the maximum number of subdomains is 64 and multiples of
! 1, 2, 3, 4, or 5 subdomains per direction are allowed
! Any desired division can be added to the list
! The information about the problem size and data decomposition is stored in the 
! matrix structure.
!

   type(user_matrix), intent(out)   :: A

      integer                       :: grid
      integer, parameter            :: Mx = 16             ! Minimum points in x-direction
      integer, parameter            :: My = 16             ! Minimum points in y-direction
      integer, parameter            :: Mz = 16             ! Minimum points in z-direction

      integer                       :: Sx, Sy, Sz, me(3)
      integer                       :: Nx, Ny, Nz, Nn
      integer, allocatable          :: Sd[:,:,:]
      integer                       :: n_procs

! Get the number of processors
      n_procs = num_images()

! Determine the number of subdomains per direction 
      select case ( n_procs)
         case(1)
            Sx = 1
            Sy = 1
            Sz = 1
         case(2)
            Sx = 2
            Sy = 1
            Sz = 1
         case(3)
            Sx = 3
            Sy = 1
            Sz = 1
         case(4)
            Sx = 2
            Sy = 2
            Sz = 1
         case(5)
            Sx = 5
            Sy = 1
            Sz = 1
         case(6)
            Sx = 3
            Sy = 2
            Sz = 1
         case(8)
            Sx = 2
            Sy = 2
            Sz = 2
         case(9)
            Sx = 3
            Sy = 3
            Sz = 1
         case(10)
            Sx = 5
            Sy = 2
            Sz = 1
         case(12)
            Sx = 3
            Sy = 2
            Sz = 2
         case(15)
            Sx = 5
            Sy = 3
            Sz = 1
         case(16)
            Sx = 4
            Sy = 2
            Sz = 2
         case(18)
            Sx = 3
            Sy = 3
            Sz = 2
         case(20)
            Sx = 5
            Sy = 2
            Sz = 2
         case(24)
            Sx = 4
            Sy = 3
            Sz = 2
         case(25)
            Sx = 5
            Sy = 5
            Sz = 1
         case(27)
            Sx = 3
            Sy = 3
            Sz = 3
         case(30)
            Sx = 5
            Sy = 3
            Sz = 2
         case(32)
            Sx = 4
            Sy = 4
            Sz = 2
         case(36)
            Sx = 4
            Sy = 3
            Sz = 3
         case(40)
            Sx = 5
            Sy = 4
            Sz = 2
         case(45)
            Sx = 5
            Sy = 3
            Sz = 3
         case(48)
            Sx = 4
            Sy = 4
            Sz = 3
         case(50)
            Sx = 5
            Sy = 5
            Sz = 2
         case(60)
            Sx = 5
            Sy = 4
            Sz = 3
         case(64)
            Sx = 4
            Sy = 4
            Sz = 4
         case(128)
            Sx = 8
            Sy = 4
            Sz = 4
         case(256)
            Sx = 8
            Sy = 8
            Sz = 4
         case(512)
            Sx = 8
            Sy = 8
            Sz = 8
         case(1024)
            Sx = 16
            Sy = 8
            Sz = 8
         case(2048)
            Sx = 16
            Sy = 16
            Sz = 8
         case(4096)
            Sx = 16
            Sy = 16
            Sz = 16
         case default
            stop 'Number of processors does not match subdomain division'
      end select

!
! Determine local subdomain coordinates:
      allocate(Sd[Sx,Sy,*])
      me = this_image(Sd)
   
!
! Determine local number of subdomain nodes:
      Nx = 2**(grid-1)*(Mx/Sx)                            ! Nodes per subdomain in x_direction
      Ny = 2**(grid-1)*(My/Sy)                            ! Nodes per subdomain in y_direction
      Nz = 2**(grid-1)*(Mz/Sz)                            ! Nodes per subdomain in z_direction
      Nn = Nx*Ny*Nz                                       ! Total number of nodes per subdomain

!
! Store the subdomain information in the user_matrix structure
      A%Nx = Nx
      A%Ny = Ny
      A%Nz = Nz
      A%Nn = Nn
      A%Sx = Sx
      A%Sy = Sy
      A%Sz = Sz
      A%me = me

   end subroutine get_subdomains

   subroutine cdr_matrix( A, grid ) 
!
! Make convection-diffusion matrix with Dirichlet conditions
!
      type(user_matrix)             :: A
      integer, intent(in)           :: grid

! Default problem, parameters can be modified here to create different problem
      real(kind=rp), parameter               :: eps = 5.e-3              ! Diffusion
!     real(kind=rp), parameter               :: eps = 1.e-1              ! Diffusion
      real(kind=rp), parameter               :: beta_x = 1.              ! Convection in x-direction
      real(kind=rp), parameter               :: beta_y = 1.              ! Convection in y-direction
      real(kind=rp), parameter               :: beta_z = 1.              ! Convection in z-direction
      real(kind=rp), parameter               :: r = -5.                  ! Reaction (note: here negative!)

! Determine the mapping of the subdomains onto the processors
      call get_subdomains( A, grid )

!
! Subdomain size
      A%Lx = 1._rp
      A%Ly = 1._rp
      A%Lz = 1._rp

! Grid sizes
      A%hx = A%Lx/(A%Nx*A%Sx+1)
      A%hy = A%Ly/(A%Ny*A%Sy+1)
      A%hz = A%Lz/(A%Nz*A%Sz+1)

! Store coefficients in matrix:
      A%c = 2.*eps/A%hx**2 + 2.*eps/A%hy**2 + 2.*eps/A%hz**2 + r 
      A%e = -eps/A%hx**2 + beta_x/(2.*A%hx)
      A%w = -eps/A%hx**2 - beta_x/(2.*A%hx)
      A%n = -eps/A%hy**2 + beta_y/(2.*A%hy)
      A%s = -eps/A%hy**2 - beta_y/(2.*A%hy)
      A%b = -eps/A%hz**2 + beta_z/(2.*A%hz)
      A%f = -eps/A%hz**2 - beta_z/(2.*A%hz)

! Boundary conditions: Dirichlet on the outside
      A%bc = 0  
      if ( A%me(1) == 1 )    A%bc(1) = 1 ! West boundary
      if ( A%me(1) == A%Sx ) A%bc(2) = 1 ! East boundary  
      if ( A%me(2) == 1 )    A%bc(3) = 1 ! South
      if ( A%me(2) == A%Sy ) A%bc(4) = 1 ! North
      if ( A%me(3) == 1 )    A%bc(5) = 1 ! Front
      if ( A%me(3) == A%Sz ) A%bc(6) = 1 ! Back

   end subroutine cdr_matrix

   subroutine helmholtz_matrix( A, grid, wavenumber, multishift, impedance ) 
!
! Make Helmholtz matrix with Neumann or Robin boundary conditions
!
! The default problem is a room with six damped walls (Sommerfeld conditions)
! If the input parameter impedance is .true. a different problem is defined:
! A room with five reflecting wall and one wall coated with a partly absorbing
! material (impedance condition). Note that this second problem is much harder
! to solve. Additional problems can be created by changing the parameter of the problem.
! 
! For the present problems we have that:
! - The number of points in each direction is 16*2**(grid-1)
! - The length of the domain is 4 m
! - The frequency is 100*2**(grid-1) Hz
! - The wavelenght is 100*2**(grid-1)/340 m
! So the number of points per wavelength is:
!   (16*2**(grid-1)/4)/(100*2**(grid-1)/340) = 13,6
!
      type(user_matrix)             :: A
      integer, intent(in)           :: grid
      real(kind=rp), optional       :: wavenumber
      logical, optional             :: multishift, impedance
      real(kind=rp)                 :: k
      logical                       :: ms, impedance_bc
      complex(kind=cp)              :: imp = (0.2,-1.5)

      if ( present(wavenumber) ) then
         k = wavenumber 
      else
         k = 1._rp
      end if
!
      if ( present(multishift) ) then
         ms = multishift 
      else
         ms = .false.
      end if
      A%multishift = ms
!
      if ( present(impedance) ) then
         impedance_bc = impedance
      else
         impedance_bc = .false.
      end if

! Determine the mapping of the subdomains onto the processors
      call get_subdomains( A, grid )

! Subdomain size
      A%Lx = 4._rp
      A%Ly = 4._rp
      A%Lz = 4._rp

! Grid sizes
      A%hx = A%Lx/(A%Nx*A%Sx-1)
      A%hy = A%Ly/(A%Ny*A%Sy-1)
      A%hz = A%Lz/(A%Nz*A%Sz-1)

! Store coefficients in matrix:
      A%c = 2./A%hx**2 + 2./A%hy**2 + 2./A%hz**2 
      A%e = -1./A%hx**2
      A%w = -1./A%hx**2 
      A%n = -1./A%hy**2
      A%s = -1./A%hy**2
      A%b = -1./A%hz**2
      A%f = -1./A%hz**2 

! Apply shift
      if ( .not. ms ) A%c = A%c - k**2

! Boundary conditions: 
      A%bc = 0  
      if ( impedance_bc ) then
! Five sides Neumann, and one impedance
         if ( A%me(1) == 1 )    A%bc(1) = 2 ! West boundary
         if ( A%me(1) == A%Sx ) A%bc(2) = 3 ! East boundary  
         if ( A%me(2) == 1 )    A%bc(3) = 2 ! South
         if ( A%me(2) == A%Sy ) A%bc(4) = 2 ! North
         if ( A%me(3) == 1 )    A%bc(5) = 2 ! Front
         if ( A%me(3) == A%Sz ) A%bc(6) = 2 ! Back
      else
! Sommerfeld on all sides
         if ( A%me(1) == 1 )    A%bc(1) = 3 ! West boundary
         if ( A%me(1) == A%Sx ) A%bc(2) = 3 ! East boundary  
         if ( A%me(2) == 1 )    A%bc(3) = 3 ! South
         if ( A%me(2) == A%Sy ) A%bc(4) = 3 ! North
         if ( A%me(3) == 1 )    A%bc(5) = 3 ! Front
         if ( A%me(3) == A%Sz ) A%bc(6) = 3 ! Back
      end if

! Sommerfeld condition
      if ( ms ) then
         A%ck = 1.
      else
         A%ck = (0._cp,1._cp)*k
      end if
      if ( impedance_bc ) A%ck = A%ck/imp

   end subroutine helmholtz_matrix

   subroutine mass_matrix( A, M )
!
! Make the mass matrix
!

      type(user_matrix), intent(in)                :: A
      type(user_matrix), intent(out)               :: M

      M = A

! Store coefficients in matrix:
      M%c = 12._rp/18._rp
      M%e =  1._rp/18._rp
      M%w =  1._rp/18._rp
      M%n =  1._rp/18._rp
      M%s =  1._rp/18._rp
      M%b =  1._rp/18._rp
      M%f =  1._rp/18._rp

   end subroutine mass_matrix

   subroutine point_source( b, A, x, y, z )
!
! Function b is point source in middle of domain
! Note that the number of gridpoints in each direction is even
! We therefore take the point right-up
!          delta(x_s) = Sum x_i,j,k/8 * h^-3
! in which x_i,j,k are the 8 gridpoints around the source
!
      type(user_matrix), intent(in)   :: A
      real(kind=rp)                   :: x, y, z
      complex(kind=cp), intent(out)   :: b(:)
      integer                         :: i_proc, j_proc, k_proc
      integer                         :: i, j, k, l

      b = 0.
!
! Convert coordinates to local indices and processor numbers
      i = int(x/A%hx)+1
      if ( i > A%Nx*A%Sx ) i = A%Nx*A%Sx
      i_proc = (i-1)/A%Nx+1
      i = i - A%Nx*(i_proc-1)

      j = int(y/A%hy)+1
      if ( j > A%Ny*A%Sy ) j = A%Ny*A%Sy
      j_proc = (j-1)/A%Ny+1
      j = j - A%Ny*(j_proc-1)

      k = int(z/A%hz)+1
      if ( k > A%Nz*A%Sz ) k = A%Nz*A%Sz
      k_proc = (k-1)/A%Nz+1
      k = k - A%Nz*(k_proc-1)
!
! To which subdomain do the coordinates correspond?
      if ( (i_proc == A%me(1)) .and. (j_proc == A%me(2)) .and. (k_proc == A%me(3)) ) then
         l = i+(j-1)*A%Nx+(k-1)*A%Nx*A%Ny
         b(l) = 1./(A%hx*A%hy*A%hz)
      end if
   end subroutine point_source

   subroutine model_solution( w, A ) 
!
! Create model solution w_model = xyz(1-x)(1-y)(1-z)
!
      type(user_matrix), intent(in) :: A
      real(kind=rp), intent(out)    :: w(A%Nn)
      real(kind=rp)                 :: w_model(A%Nx,A%Ny,A%Nz)
      integer                       :: Nx, Ny, Nz, Nn
      integer                       :: i, j, k

! Coordinates
      real(kind=rp)                 :: x, y, z

! Model solution
      do concurrent (i=1:A%Nx,j=1:A%Ny,k=1:A%Nz)
         x = (A%me(1)-1)*A%Nx*A%hx+i*A%hx
         y = (A%me(2)-1)*A%Ny*A%hy+j*A%hy
         z = (A%me(3)-1)*A%Nz*A%hz+k*A%hz
         w_model(i,j,k) = sqrt( x*y*z*(1.-x)*(1.-y)*(1.-z) )
      end do
      w = reshape(w_model,(/A%Nn/))

   end subroutine model_solution

end module user_module

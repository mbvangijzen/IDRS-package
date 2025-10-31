!
! The pp_idrs_module contains functions that take a minimum number of problem 
! dependend parameters as input and return the solution(s) of the linear system
! as output. They are designed for easy experimentation.
! When called, the functions 
!         -read parameters from the command line
!         -output which solution method is used
!         -store the matrix in the correct format 
!         -setup the polynomial preconditioner
!         -run the specified IDR(s) method
!         -report number of iterations and residual norms
!         -(optionally) plot convergence curves and/or ritzvalues
!
!    x = pp_idrs( K, b, sigma, M, C, P, center, foci, seed, shift )
! Input:
!      K: system matrix. User defined type user_matrix. Required.
!      b: right hand side vector(s). One or two dimensional complex or real array. Required.
!         Should be one dimensional array of length the number of equations on the processor for
!         single right-hand-side problems or multishift problems, or two-dimension for multiple 
!         rhs-problems. 
!  sigma: shifts of the multi-shift problem. Complex or real array of length the number of shifts.
!         Required for the multi-shift problem.
!      M: mass matrix or preconditioning matrix. User defined type user_matrix. Optional.
!      C: daming matrix. User defined type user_matrix. Optional.
!      P: Preconditioner. Diagonal matrix or LU factorisation.
!         Real or complex array. Optional
! center: center around the spectrum of the (scaled) system matrix. Complex or real scalar.
!         Optional.
!   foci: foci of ellips around the spectrum of the (scaled) system matrix.  
!         Complex or scaler array of dimension 2. Optional
!   seed: Seed shift to the system. Complex or real scalar. Optional.
!  shift: Shift to the preconditioner. Complex or real scalar. Optional.
!
! Note that all (non-matrix) parameter must be of the same type and kind, either complex or real.
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2025 Martin van Gijzen
!

module pp_idrs_module

   use precision_module ! Set the precisions for real and complex arithmatic

   use matrix_module    ! Always needed is a matrix module that defines the matrix type,
                        ! the matrix-vector multiplication and the preconditioning operation
                        ! Here we use the matrix module supplied with the package
                        ! The matrix module needs a user module that defines the
                        ! user_matrix type and the multiplication with a user_matrix

   use idrs_module      ! The idrs solvers

   use interface_module ! To read command line parameters

   use ritz_module      ! To estimate iteration parameters

   use dense_la_module  ! To calculate the Ritz values

   implicit none

   private
   public :: PP_IDRS

   interface PP_IDRS
      module procedure CPP_IDRS, RPP_IDRS, CPP_MR_IDRS, RPP_MR_IDRS, CPP_MS_IDRS, RPP_MS_IDRS, &
                       CPP_MS_MR_IDRS, RPP_MS_MR_IDRS
   end interface

contains

   function CPP_IDRS( K, b, M, C, P, center, foci, seed, shift, V ) result(x)
!
! PP_IDRS for complex single right-hand side system
! 

   type(user_matrix), intent(in)                          :: K
   complex(kind=cp), intent(in)                           :: b(:)
   type(user_matrix), intent(in), optional                :: M
   type(user_matrix), intent(in), optional                :: C
   type(user_matrix), intent(in), optional                :: P
   complex(kind=cp), intent(in), optional                 :: center
   complex(kind=cp), intent(in), optional                 :: foci(2)
   complex(kind=cp), intent(in), optional                 :: seed
   complex(kind=cp), intent(in), optional                 :: shift
   complex(kind=cp), intent(in), optional                 :: V(:,:)
   complex(kind=cp)                                       :: x(size(b))

   complex(kind=cp), allocatable                          :: x0(:), U0(:,:), omega(:), H(:,:)
   integer                                                :: neq

! Declaration and initialization of IDRS parameters:
   include "initialize.inc"

   neq  = size(b) 
   allocate( x0(neq) )
   x0 = 0.

! Save convergence history?
   if ( plot_conv ) allocate( resvec(maxit+1) )

! Output message, which methods:
   include "output_methods.inc"

! Store the matrices in the matrix structure:
   include "complex_matrix.inc"

! Initial solve to determine iteration parameters, set-up preconditioner:
   include "complex_preconditioner.inc"

! Compute deflation matrix:
   if ( present(V) ) call make_deflation( A, V )

! Solve the system:
   include "run_idrs.inc"

   end function CPP_IDRS

   function RPP_IDRS( K, b, M, C, P, center, foci, seed, shift, V ) result(x)
!
! PP_IDRS for real single right-hand side system
! 
   type(user_matrix), intent(in)                       :: K
   real(kind=rp), intent(in)                           :: b(:)
   type(user_matrix), intent(in), optional             :: M
   type(user_matrix), intent(in), optional             :: C
   type(user_matrix), intent(in), optional             :: P
   real(kind=rp), intent(in), optional                 :: center
   real(kind=rp), intent(in), optional                 :: foci(2)
   real(kind=rp), intent(in), optional                 :: seed
   real(kind=rp), intent(in), optional                 :: shift
   real(kind=rp), intent(in), optional                 :: V(:,:)
   real(kind=rp)                                       :: x(size(b))

   real(kind=rp), allocatable                          :: x0(:), U0(:,:), omega(:), H(:,:)
   integer                                             :: neq

! Declaration and initialization of IDRS parameters:
   include "initialize.inc"

   neq = size(b) 
   allocate( x0(neq) )
   x0 = 0.

! Save convergence history?
   if ( plot_conv ) allocate( resvec(maxit+1) )

! Output message, which methods:
   include "output_methods.inc"

! Store the matrices in the matrix structure:
   include "real_matrix.inc"

! Initial solve to determine iteration parameters, set-up preconditioner:
   include "real_preconditioner.inc"

! Compute deflation matrix:
   if ( present(V) ) call make_deflation( A, V )

! Solve the system:
   include "run_idrs.inc"

   end function RPP_IDRS

   function CPP_MR_IDRS( K, b_all, M, C, P, center, foci, seed, shift, V ) result(x_all)
!
! PP_IDRS for complex multiple right-hand side system
! 
   type(user_matrix), intent(in)                          :: K
   complex(kind=cp), intent(in)                           :: b_all(:,:)
   type(user_matrix), intent(in), optional                :: M
   type(user_matrix), intent(in), optional                :: C
   type(user_matrix), intent(in), optional                :: P
   complex(kind=cp), intent(in), optional                 :: center
   complex(kind=cp), intent(in), optional                 :: foci(2)
   complex(kind=cp), intent(in), optional                 :: seed
   complex(kind=cp), intent(in), optional                 :: shift
   complex(kind=cp), intent(in), optional                 :: V(:,:)
   complex(kind=cp)                                       :: x_all(size(b_all,1), size(b_all,2) )

   complex(kind=cp), allocatable                          :: x0(:), b(:), x(:), U0(:,:), omega(:), H(:,:)

   integer                                                :: irhs, neq, nrhs

! Declaration and initialization of IDRS parameters:
   include "initialize.inc"

   neq  = size(b_all,1)
   nrhs = size(b_all,2)
   allocate( x0(neq), x(neq), b(neq) )
   x0 = 0.

! Save convergence history?
   if ( plot_conv ) allocate( resvec(maxit+1) )

! Output message, which methods:
   include "output_methods.inc"

   include "complex_matrix.inc"

   call system_clock ( t0, clock_rate, clock_max )

! Initial solve to determine iteration parameters, set-up preconditioner:
   b = b_all(:,1)
   include "complex_preconditioner.inc"

! Compute deflation matrix:
   if ( present(V) ) call make_deflation( A, V )

! Compute recycle space:
   include "recycle.inc"

! Solve the system
   include "run_mr_idrs.inc"

   call system_clock( t1, clock_rate, clock_max )

   if ( my_proc == 1 ) &
      write(*,'(a,f8.2,a)') 'Total elapsed time     = ', real(t1-t0)/real(clock_rate), 's.'

   end function CPP_MR_IDRS

   function RPP_MR_IDRS( K, b_all, M, C, P, center, foci, seed, shift, V ) result(x_all)
!
! PP_IDRS for real multiple right-hand side system
! 
   type(user_matrix), intent(in)                       :: K
   real(kind=rp), intent(in)                           :: b_all(:,:)
   type(user_matrix), intent(in), optional             :: M
   type(user_matrix), intent(in), optional             :: C
   type(user_matrix), intent(in), optional             :: P
   real(kind=rp), intent(in), optional                 :: center
   real(kind=rp), intent(in), optional                 :: foci(2)
   real(kind=rp), intent(in), optional                 :: seed
   real(kind=rp), intent(in), optional                 :: shift
   real(kind=rp), intent(in), optional                 :: V(:,:)
   real(kind=rp)                                       :: x_all(size(b_all,1),size(b_all,2))

   real(kind=rp), allocatable                          :: x0(:), b(:), x(:), U0(:,:), omega(:), H(:,:)

   integer                                             :: irhs, neq, nrhs

! Declaration and initialization of IDRS parameters:
   include "initialize.inc"

   neq  = size(b_all,1)
   nrhs = size(b_all,2)
   allocate( x0(neq), x(neq), b(neq) )
   x0 = 0.

! Save convergence history?
   if ( plot_conv ) allocate( resvec(maxit+1) )

! Output message, which methods:
   include "output_methods.inc"

! Store the matrices in the matrix structure:
   include "real_matrix.inc"

! Compute deflation matrix:
   if ( present(V) ) call make_deflation( A, V )

   call system_clock ( t0, clock_rate, clock_max )

! Initial solve to determine iteration parameters, set-up preconditioner:
   b = b_all(:,1)
   include "real_preconditioner.inc"

! Compute recycle space:
   include "recycle.inc"

! Solve the systems
   include "run_mr_idrs.inc"

   call system_clock ( t1, clock_rate, clock_max )

   if ( my_proc == 1 ) &
      write(*,'(a,f8.2,a)') 'Total elapsed time     = ', real(t1-t0)/real(clock_rate), 's.'

   end function RPP_MR_IDRS

   function CPP_MS_IDRS( K, rhs, sigma, M, C, P, center, foci, seed, shift, V ) result(solution)
!
! PP_IDRS for complex multi-shift system
! 
   type(user_matrix), intent(in)                          :: K
   complex(kind=cp), intent(in)                           :: rhs(:)
   complex(kind=cp), intent(in)                           :: sigma(:)
   type(user_matrix), intent(in), optional                :: M
   type(user_matrix), intent(in), optional                :: C
   type(user_matrix), intent(in), optional                :: P
   complex(kind=cp), intent(in), optional                 :: center
   complex(kind=cp), intent(in), optional                 :: foci(2)
   complex(kind=cp), intent(in), optional                 :: seed
   complex(kind=cp), intent(in), optional                 :: shift
   complex(kind=rp), intent(in), optional                 :: V(:,:)
   complex(kind=cp)                                       :: solution(size(rhs),size(sigma))

   complex(kind=cp), allocatable                          :: x0(:), U0(:,:), omega(:), H(:,:), x(:,:), b(:)
   complex(kind=cp), allocatable                          :: sigma_p(:)
   integer                                                :: Nsigma

   integer                                                :: neq, nn, i_sigma

! Declaration and initialization of IDRS parameters:
   include "initialize.inc"
   multishift = .true.

   Nsigma = size(sigma)
   allocate( sigma_p(Nsigma), resnrm(Nsigma) )

! Save convergence history?
if ( plot_conv ) allocate( ms_resvec(maxit+1,Nsigma) )

   nn = size(rhs)
   if ( present(C) ) then
      neq  = 2*nn
   else
      neq = nn
   end if
   allocate( x0(neq), x(neq,Nsigma), b(neq) )
   x0 = 0.
   if ( present(C) ) then
      b(1:nn)     = -rhs
      b(nn+1:neq) = 0._cp
   else
      b = rhs
   end if

! Output message, which methods:
   include "output_methods.inc"

! Store the matrices in the matrix structure:
   include "complex_matrix.inc"

! Make the preconditioner:
   include "complex_preconditioner.inc"

! Compute deflation matrix:
   if ( present(V) ) call make_deflation( A, V )

! Solve the system:
   include "run_ms_idrs.inc"

   end function CPP_MS_IDRS

   function RPP_MS_IDRS( K, rhs, sigma, M, C, P, center, foci, seed, shift, V ) result(solution)
!
! PP_IDRS for real multi-shift system
! 
   type(user_matrix), intent(in)                       :: K
   real(kind=rp), intent(in)                           :: rhs(:)
   real(kind=rp), intent(in)                           :: sigma(:)
   type(user_matrix), intent(in), optional             :: M
   type(user_matrix), intent(in), optional             :: C
   type(user_matrix), intent(in), optional             :: P
   real(kind=rp), intent(in), optional                 :: center
   real(kind=rp), intent(in), optional                 :: foci(2)
   real(kind=rp), intent(in), optional                 :: seed
   real(kind=rp), intent(in), optional                 :: shift
   real(kind=rp), intent(in), optional                 :: V(:,:)
   real(kind=rp)                                       :: solution(size(rhs),size(sigma))

   real(kind=rp), allocatable                          :: x0(:), U0(:,:), omega(:), H(:,:), x(:,:), b(:)
   real(kind=rp), allocatable                          :: sigma_p(:)
   integer                                             :: Nsigma

   integer                                             :: neq, i_sigma, nn

! Declaration and initialization of IDRS parameters:
   include "initialize.inc"
   multishift = .true.

   Nsigma = size(sigma)
   allocate( sigma_p(Nsigma), resnrm(Nsigma) )

! Save convergence history?
if ( plot_conv ) allocate( ms_resvec(maxit+1,Nsigma) )

   nn = size(rhs)
   if ( present(C) ) then
      neq = 2*nn
   else
      neq = nn
   end if
   allocate( x0(neq), x(neq,Nsigma), b(neq) )
   x0 = 0.
   if ( present(C) ) then
      b(1:nn)     = -rhs
      b(nn+1:neq) = 0._rp
   else
      b = rhs
   end if

! Output message, which methods:
   include "output_methods.inc"

! Store the matrices in the matrix structure:
   include "real_matrix.inc"

! Make the preconditioner:
   include "real_preconditioner.inc"

! Compute deflation matrix:
   if ( present(V) ) call make_deflation( A, V )

! Solve the system:
   include "run_ms_idrs.inc"

   end function RPP_MS_IDRS

   function RPP_MS_MR_IDRS( K, rhs, sigma, M, C, P, center, foci, seed, shift, V ) result(solution)
!
! PP_IDRS for real multi-shift system
! 
   type(user_matrix), intent(in)                       :: K
   real(kind=rp), intent(in)                           :: rhs(:,:)
   real(kind=rp), intent(in)                           :: sigma(:)
   type(user_matrix), intent(in), optional             :: M
   type(user_matrix), intent(in), optional             :: C
   type(user_matrix), intent(in), optional             :: P
   real(kind=rp), intent(in), optional                 :: center
   real(kind=rp), intent(in), optional                 :: foci(2)
   real(kind=rp), intent(in), optional                 :: seed
   real(kind=rp), intent(in), optional                 :: shift
   real(kind=rp), intent(in), optional                 :: V(:,:)
   real(kind=rp)                                       :: solution(size(rhs,1),size(sigma),size(rhs,2))

   real(kind=rp), allocatable                          :: x0(:), U0(:,:), omega(:), H(:,:), x(:,:), b(:), b_all(:,:)
   real(kind=rp), allocatable                          :: sigma_p(:)
   integer                                             :: Nsigma

   integer                                             :: neq, i_sigma, nn, nrhs, irhs

! Declaration and initialization of IDRS parameters:
   include "initialize.inc"
   multishift = .true.

   Nsigma = size(sigma)
   allocate( sigma_p(Nsigma), resnrm(Nsigma) )

! Save convergence history?
if ( plot_conv ) allocate( ms_resvec(maxit+1,Nsigma) )

   nn   = size(rhs,1)
   nrhs = size(rhs,2)
   if ( present(C) ) then
      neq = 2*nn
   else
      neq = nn
   end if
   allocate( x0(neq), x(neq,Nsigma), b(neq), b_all(neq,nrhs) )
   x0 = 0.
   if ( present(C) ) then
      b_all(1:nn,:) = -rhs
      b_all(nn+1:neq,:) = 0._rp
   else
      b_all = rhs
   end if
   b = b_all(:,1)

! Output message, which methods:
   include "output_methods.inc"

! Store the matrices in the matrix structure:
   include "real_matrix.inc"

! Make the preconditioner:
   include "real_preconditioner.inc"

! Compute deflation matrix:
   if ( present(V) ) call make_deflation( A, V )

! Solve the system:
   include "run_ms_mr_idrs.inc"

   end function RPP_MS_MR_IDRS

   function CPP_MS_MR_IDRS( K, rhs, sigma, M, C, P, center, foci, seed, shift, V ) result(solution)
!
! PP_IDRS for real multi-shift system
! 
   type(user_matrix), intent(in)                       :: K
   complex(kind=cp), intent(in)                        :: rhs(:,:)
   complex(kind=cp), intent(in)                        :: sigma(:)
   type(user_matrix), intent(in), optional             :: M
   type(user_matrix), intent(in), optional             :: C
   type(user_matrix), intent(in), optional             :: P
   complex(kind=cp), intent(in), optional              :: center
   complex(kind=cp), intent(in), optional              :: foci(2)
   complex(kind=cp), intent(in), optional              :: seed
   complex(kind=cp), intent(in), optional              :: shift
   complex(kind=cp), intent(in), optional              :: V(:,:)
   complex(kind=cp)                                    :: solution(size(rhs,1),size(sigma),size(rhs,2))

   complex(kind=cp), allocatable                       :: x0(:), U0(:,:), omega(:), H(:,:), x(:,:), b(:), b_all(:,:)
   complex(kind=cp), allocatable                       :: sigma_p(:)
   integer                                             :: Nsigma

   integer                                             :: neq, i_sigma, nn, nrhs, irhs

! Declaration and initialization of IDRS parameters:
   include "initialize.inc"
   multishift = .true.

   Nsigma = size(sigma)
   allocate( sigma_p(Nsigma), resnrm(Nsigma) )

! Save convergence history?
if ( plot_conv ) allocate( ms_resvec(maxit+1,Nsigma) )

   nn   = size(rhs,1)
   nrhs = size(rhs,2)
   if ( present(C) ) then
      neq = 2*nn
   else
      neq = nn
   end if
   allocate( x0(neq), x(neq,Nsigma), b(neq), b_all(neq,nrhs) )
   x0 = 0.
   if ( present(C) ) then
      b_all(1:nn,:) = -rhs
      b_all(nn+1:neq,:) = 0._cp
   else
      b_all = rhs
   end if
   b = b_all(:,1)

! Output message, which methods:
   include "output_methods.inc"

! Store the matrices in the matrix structure:
   include "complex_matrix.inc"

! Make the preconditioner:
   include "complex_preconditioner.inc"

! Compute deflation matrix:
   if ( present(V) ) call make_deflation( A, V )

! Solve the system:
   include "run_ms_mr_idrs.inc"

   end function CPP_MS_MR_IDRS

end module pp_idrs_module

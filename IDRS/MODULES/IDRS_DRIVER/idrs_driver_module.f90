!
! The idrs_driver_module contains driver functions that take a minimum number of problem 
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
!    x = idrs_driver( K, b, sigma, M, C, P, center, foci, seed, shift )
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
! Copyright:(c) 2026 Martin van Gijzen
!

module idrs_driver_module

   use precision_module ! Set the precisions for real and complex arithmatic

   use matrix_module    ! Always needed is a matrix module that defines the matrix type,
                        ! the matrix-vector multiplication and the preconditioning operation
                        ! Here we use the matrix module supplied with the package
                        ! The matrix module needs a user module that defines the
                        ! user_matrix type and the multiplication with a user_matrix

   use idrs_module      ! The idrs solvers

   use interface_module ! To read command line parameters

   use recycle_module   ! To estimate Schur vectors for recycling

   use dense_la_module  ! To calculate the Ritz values

   implicit none

   private
   public :: IDRS_DRIVER, IDREIG_DRIVER

   interface IDRS_DRIVER
      module procedure CIDRS_DRIVER, RIDRS_DRIVER, CMR_IDRS_DRIVER, RMR_IDRS_DRIVER, CMS_IDRS_DRIVER, &
                       RMS_IDRS_DRIVER, CMS_MR_IDRS_DRIVER, RMS_MR_IDRS_DRIVER
   end interface

   interface IDREIG_DRIVER
      module procedure CIDREIG_DRIVER, RIDREIG_DRIVER
   end interface

contains

   function CIDRS_DRIVER( K, b, M, C, P, center, foci, seed, shift, V ) result(x)
!
! IDRS_DRIVER for complex single right-hand side system
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

   complex(kind=cp), allocatable                          :: x0(:), U(:,:), U0(:,:), H(:,:), omega(:)
   complex(kind=cp), allocatable                          :: eigval(:)
   integer                                                :: neq, i

! Declaration and initialization of IDRS parameters:
   include "initialize.inc"

   neq  = size(b) 
   allocate( x0(neq) )
   x0 = 0.

   if ( present(seed) )        cseed  = seed
   if ( .not. present(shift) ) cshift = cseed

! Save convergence history?
   if ( plot_conv ) allocate( resvec(maxit+1) )

! Output message, which methods:
   include "output_methods.inc"

! Store the matrices in the matrix structure:
   include "complex_system_matrix.inc"

! Solve the system:
   include "run_idrs.inc"

   end function CIDRS_DRIVER

   function RIDRS_DRIVER( K, b, M, C, P, center, foci, seed, shift, V ) result(x)
!
! IDRS_DRIVER for real single right-hand side system
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

   real(kind=rp), allocatable                          :: x0(:), U(:,:), U0(:,:), H(:,:), omega(:)
   real(kind=rp), allocatable                          :: eigval(:)
   integer                                             :: neq

! Declaration and initialization of IDRS parameters:
   include "initialize.inc"

   neq = size(b) 
   allocate( x0(neq) )
   x0 = 0.

   if ( present(seed) )        rseed  = seed
   if ( .not. present(shift) ) rshift = rseed

! Save convergence history?
   if ( plot_conv ) allocate( resvec(maxit+1) )

! Output message, which methods:
   include "output_methods.inc"

! Store the matrices in the matrix structure:
   include "real_system_matrix.inc"

! Solve the system:
   include "run_idrs.inc"

   end function RIDRS_DRIVER

   function CMR_IDRS_DRIVER( K, b_all, M, C, P, center, foci, seed, shift, V ) result(x_all)
!
! IDRS_DRIVER for complex multiple right-hand side system
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

   complex(kind=cp), allocatable                          :: x0(:), b(:), x(:), U(:,:), U0(:,:), H(:,:), omega(:)
   complex(kind=cp), allocatable                          :: eigval(:)

   integer                                                :: i, irhs, neq, nrhs

! Declaration and initialization of IDRS parameters:
   include "initialize.inc"

   neq  = size(b_all,1)
   nrhs = size(b_all,2)
   allocate( x0(neq), x(neq), b(neq) )
   x0 = 0.

   if ( present(seed) )        cseed  = seed
   if ( .not. present(shift) ) cshift = cseed

! Save convergence history?
   if ( plot_conv ) allocate( resvec(maxit+1) )

! Output message, which methods:
   include "output_methods.inc"

   call system_clock ( t0, clock_rate, clock_max )

! Initial solve to determine iteration parameters, set-up preconditioner:
   b = b_all(:,1)

! Compute deflation preconditioner:
   include "complex_system_matrix.inc"

! Compute recycle space:
   include "recycle.inc"

! Solve the system
   include "run_mr_idrs.inc"

   call system_clock( t1, clock_rate, clock_max )

   if ( my_proc == 1 ) &
      write(*,'(a,f8.2,a)') 'Total elapsed time     = ', real(t1-t0)/real(clock_rate), 's.'

   end function CMR_IDRS_DRIVER

   function RMR_IDRS_DRIVER( K, b_all, M, C, P, center, foci, seed, shift, V ) result(x_all)
!
! IDRS_DRIVER for real multiple right-hand side system
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

   real(kind=rp), allocatable                          :: x0(:), b(:), x(:), U(:,:), U0(:,:), H(:,:), omega(:)
   real(kind=rp), allocatable                          :: eigval(:)

   integer                                             :: i, irhs, neq, nrhs

! Declaration and initialization of IDRS parameters:
   include "initialize.inc"

   neq  = size(b_all,1)
   nrhs = size(b_all,2)
   allocate( x0(neq), x(neq), b(neq) )
   x0 = 0.

   if ( present(seed) )        rseed  = seed
   if ( .not. present(shift) ) rshift = rseed

! Save convergence history?
   if ( plot_conv ) allocate( resvec(maxit+1) )

! Output message, which methods:
   include "output_methods.inc"

   call system_clock ( t0, clock_rate, clock_max )

   b = b_all(:,1)

! Store the matrices in the matrix structure:
   include "real_system_matrix.inc"

! Compute recycle space:
   include "recycle.inc"

! Solve the systems
   include "run_mr_idrs.inc"

   call system_clock ( t1, clock_rate, clock_max )

   if ( my_proc == 1 ) &
      write(*,'(a,f8.2,a)') 'Total elapsed time     = ', real(t1-t0)/real(clock_rate), 's.'

   end function RMR_IDRS_DRIVER

   function CMS_IDRS_DRIVER( K, rhs, sigma, M, C, P, center, foci, seed, shift, V ) result(solution)
!
! IDRS_DRIVER for complex multi-shift system
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

   complex(kind=cp), allocatable                          :: x0(:), U(:,:), U0(:,:), H(:,:), x(:,:), b(:), omega(:)
   complex(kind=cp), allocatable                          :: sigma_p(:)
   complex(kind=cp), allocatable                          :: eigval(:)
   integer                                                :: Nsigma

   integer                                                :: neq, nn, i_sigma

! Declaration and initialization of IDRS parameters:
   include "initialize.inc"
   multishift = .true.

   Nsigma = size(sigma)

   if ( present(seed)  ) then
      cseed = seed
   else
      cseed  = sigma(Nsigma)
   end if
   if ( present(shift) ) then
      cshift = shift
   else
      cshift = sum(sigma)/Nsigma
   end if

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
   include "complex_system_matrix.inc"

! Solve the system:
   include "run_ms_idrs.inc"

   end function CMS_IDRS_DRIVER

   function RMS_IDRS_DRIVER( K, rhs, sigma, M, C, P, center, foci, seed, shift, V ) result(solution)
!
! IDRS_DRIVER for real multi-shift system
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

   real(kind=rp), allocatable                          :: x0(:), U(:,:), U0(:,:), H(:,:), x(:,:), b(:), omega(:)
   real(kind=rp), allocatable                          :: sigma_p(:)
   real(kind=rp), allocatable                          :: eigval(:)
   integer                                             :: Nsigma

   integer                                             :: neq, i_sigma, nn

! Declaration and initialization of IDRS parameters:
   include "initialize.inc"
   multishift = .true.

   Nsigma = size(sigma)

   if ( present(seed) ) then
      rseed = seed
   else
      rseed  = sigma(Nsigma)
   end if
   if ( present(shift) ) then
      rshift = shift
   else
      rshift = sum(sigma)/Nsigma
   end if

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
   include "real_system_matrix.inc"

! Solve the system:
   include "run_ms_idrs.inc"

   end function RMS_IDRS_DRIVER

   function CMS_MR_IDRS_DRIVER( K, rhs, sigma, M, C, P, center, foci, seed, shift, V ) result(solution)
!
! IDRS_DRIVER for real multi-shift system
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

   complex(kind=cp), allocatable                       :: x0(:), U(:,:), U0(:,:), H(:,:), x(:,:), b(:), b_all(:,:), omega(:)
   complex(kind=cp), allocatable                       :: sigma_p(:)
   complex(kind=cp), allocatable                       :: eigval(:)
   integer                                             :: Nsigma

   integer                                             :: neq, i_sigma, nn, nrhs, irhs

! Declaration and initialization of IDRS parameters:
   include "initialize.inc"
   multishift = .true.

   Nsigma = size(sigma)
   allocate( sigma_p(Nsigma), resnrm(Nsigma) )

   if ( present(seed)  ) then
      cseed = seed
   else
      cseed  = sigma(Nsigma)
   end if
   if ( present(shift) ) then
      cshift = shift
   else
      cshift = sum(sigma)/Nsigma
   end if

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
   include "complex_system_matrix.inc"

! Solve the system:
   include "run_ms_mr_idrs.inc"

   end function CMS_MR_IDRS_DRIVER

   function RMS_MR_IDRS_DRIVER( K, rhs, sigma, M, C, P, center, foci, seed, shift, V ) result(solution)
!
! IDRS_DRIVER for real multi-shift system
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

   real(kind=rp), allocatable                          :: x0(:), U(:,:), U0(:,:), H(:,:), x(:,:), b(:), b_all(:,:), omega(:)
   real(kind=rp), allocatable                          :: sigma_p(:)
   real(kind=rp), allocatable                          :: eigval(:)
   integer                                             :: Nsigma

   integer                                             :: neq, i_sigma, nn, nrhs, irhs

! Declaration and initialization of IDRS parameters:
   include "initialize.inc"
   multishift = .true.

   Nsigma = size(sigma)
   allocate( sigma_p(Nsigma), resnrm(Nsigma) )

   if ( present(seed)  ) then
      rseed = seed
   else
      rseed  = sigma(Nsigma)
   end if
   if ( present(shift) ) then
      rshift = shift
   else
      rshift = sum(sigma)/Nsigma
   end if

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
   include "real_system_matrix.inc"

! Solve the system:
   include "run_ms_mr_idrs.inc"

   end function RMS_MR_IDRS_DRIVER

   subroutine CIDREIG_DRIVER( lambda, eigvec, K, start_vec, M, C, P, center, foci, seed, shift, V )
!
! IDREIG_DRIVER for solving complex eigenvalue problem
! 

   complex(kind=cp), intent(out)                          :: lambda(:)
   complex(kind=cp), intent(out)                          :: eigvec(:,:)
   type(user_matrix), intent(in)                          :: K
   complex(kind=cp), intent(in)                           :: start_vec(:)
   type(user_matrix), intent(in), optional                :: M
   type(user_matrix), intent(in), optional                :: C
   type(user_matrix), intent(in), optional                :: P
   complex(kind=cp), intent(in), optional                 :: center
   complex(kind=cp), intent(in), optional                 :: foci(2)
   complex(kind=cp), intent(in), optional                 :: seed
   complex(kind=cp), intent(in), optional                 :: shift
   complex(kind=cp), intent(in), optional                 :: V(:,:)

   complex(kind=cp), allocatable                          :: U(:,:), b(:), H(:,:), omega(:)
   complex(kind=cp), allocatable                          :: eigval(:)
   integer                                                :: neq, nn
   integer                                                :: i, ieig
   complex(kind=cp)                                       :: ei
   real(kind=rp)                                          :: normu

! Declaration and initialization of IDRS parameters:
   include "initialize.inc"

   nn   = size(eigvec,1) 
   neig = size(eigvec,2)
   eigen = .true.

   if ( present(C) ) then
      neq  = 2*nn
   else
      neq = nn
   end if
   allocate( b(neq) )
   if ( present(C) ) then
      b(1:nn)     = -start_vec
      b(nn+1:neq) = 0._rp
   else
      b = start_vec
   end if

   if ( present(seed) ) cseed  = seed
   if ( present(shift) ) then
      cshift = shift
   else
      cshift = cseed
   end if

! Save convergence history?
   if ( plot_conv ) allocate( eig_resvec(maxit,neig) )

! Output message, which methods:
   include "output_methods.inc"

! Store the matrices in the matrix structure:
   include "complex_system_matrix.inc"

! Solve the system:
   include "run_idreig.inc"

   end subroutine CIDREIG_DRIVER

   subroutine RIDREIG_DRIVER( lambda, eigvec, K, start_vec, M, C, P, center, foci, seed, shift, V )
!
! IDREIG_DRIVER for solving real eigenvalue problem
! 

   real(kind=rp), intent(out)                             :: lambda(:)
   real(kind=rp), intent(out)                             :: eigvec(:,:)
   type(user_matrix), intent(in)                          :: K
   real(kind=rp), intent(in)                              :: start_vec(:)
   type(user_matrix), intent(in), optional                :: M
   type(user_matrix), intent(in), optional                :: C
   type(user_matrix), intent(in), optional                :: P
   real(kind=rp), intent(in), optional                    :: center
   real(kind=rp), intent(in), optional                    :: foci(2)
   real(kind=rp), intent(in), optional                    :: seed
   real(kind=rp), intent(in), optional                    :: shift
   real(kind=rp), intent(in), optional                    :: V(:,:)

   real(kind=rp), allocatable                             :: U(:,:), b(:), H(:,:), omega(:)
   real(kind=rp), allocatable                             :: eigval(:)
   integer                                                :: neq, nn
   integer                                                :: i, ieig
   complex(kind=cp)                                       :: ei
   real(kind=rp)                                          :: normu

! Declaration and initialization of IDRS parameters:
   include "initialize.inc"

   nn   = size(eigvec,1) 
   neig = size(eigvec,2)
   eigen = .true.

   if ( present(C) ) then
      neq  = 2*nn
   else
      neq = nn
   end if
   allocate( b(neq) )
   if ( present(C) ) then
      b(1:nn)     = -start_vec
      b(nn+1:neq) = 0._rp
   else
      b = start_vec
   end if

   if ( present(seed) )        rseed  = seed
   if ( present(shift) ) then
      rshift = shift
   else
      rshift = rseed
   end if

! Save convergence history?
   if ( plot_conv ) allocate( eig_resvec(maxit,neig) )

! Output message, which methods:
   include "output_methods.inc"

! Store the matrices in the matrix structure:
   include "real_system_matrix.inc"

! Solve the system:
   include "run_idreig.inc"

   end subroutine RIDREIG_DRIVER

end module idrs_driver_module

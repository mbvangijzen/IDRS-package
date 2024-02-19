!
! Include file to declare standard IDRS parameter and to read the command line options
!
! Allowed command line options are
!
!    Option                     Description                                  Default
! ======================================================================================
!  -idrs          Use idrs as solution method.                           Do not use idrs
!                 Should come with integer parameter s > 0                       s = 4
!  -bicgstab      Use bicgstab as solution method                        Do not use bicgstab
!  -qmridr        Use qmridr as solution method.                         Do not use qmridr
!                 Should come with integer parameter s > 0                       s = 4
!  -maxit         Maximum number of iterations                               maxit = 1000
!  -tol           Tolerance, stop if ||r||/||b|| < tol                         tol = 1e-6
!  -in_s          Parameter s for inner iterations, used if qmridr            in_s = 1
!                 is used as nested method, with idrs as inner method
!  -in_it         Maximum number of inner iterations (only for (ms)qmridr)   in_it = 0
!  -in_tol        Tolerance for inner iterations                            in_tol = 1e-1
!  -recycle       Create s-dimenional recylce space,                       No recycling
!  -omega         Number of user defined parameters omega                  n_omega = 0
!  -ritz          Number of ritz values                                     n_ritz = 0
!  -chebyshev     If specified a Chebyshev preconditioner should be used.   degree = 0
!                 Should come with the degree of the polynomial
!  -neumann       If specified a Neumann preconditioner should be used.     degree = 0
!                 Should come with the degree of the polynomial
!  -hessenberg    If specified a Hessenberg preconditioner should be used.  degree = 0
!                 Should come with the degree of the polynomial
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2024 Martin van Gijzen
!
!

! Method:
   logical                        :: use_idrs, use_qmridr, use_bicgstab

! Iteration parameters:
   integer                        :: s
   real(kind=rp)                  :: tol
   integer                        :: maxit
   logical                        :: recycle
   integer                        :: in_it
   integer                        :: in_s
   real(kind=rp)                  :: in_tol

! Preconditioner:
   logical                        :: use_chebyshev  
   logical                        :: use_neumann   
   logical                        :: use_hessenberg
   integer                        :: degree
   integer                        :: n_omega
   integer                        :: n_ritz

! Defaults for the method:
   logical, parameter             :: default_idrs       = .false.
   logical, parameter             :: default_qmridr     = .false.
   logical, parameter             :: default_bicgstab   = .false.

! Defaults for the iterations-parameters:
   integer, parameter             :: default_s = 4
   real, parameter                :: default_tol = 1e-6
   integer, parameter             :: default_maxit = 1000
   logical, parameter             :: default_recycle = .false.
   integer, parameter             :: default_in_it = 0
   integer, parameter             :: default_in_s = 1
   real, parameter                :: default_in_tol = 1e-1

! Defaults for the preconditioner:
   logical, parameter             :: default_chebyshev  = .false.
   logical, parameter             :: default_neumann    = .false.
   logical, parameter             :: default_hessenberg = .false.
   integer, parameter             :: default_degree = 0
   integer, parameter             :: default_ritz = 0
   integer, parameter             :: default_omega = 0

! System matrix and preconditioner (or Mass matrix):
   type(matrix)                   :: A, M

! Parameters to compute polynomial preconditioner
   real(kind=rp)                  :: rcenter = 1.
   real(kind=rp)                  :: rfoci(2) = [0.25, 1.75]
   real(kind=rp)                  :: rseed = 0.
   complex(kind=cp)               :: ccenter = 1.
   complex(kind=cp)               :: cfoci(2) = [0.25, 1.75]
   complex(kind=cp)               :: cseed = 0.

! Other declarations:
   integer                        :: preconditioner
   integer                        :: variant
   integer                        :: flag
   integer                        :: iter, ini_iter
   real(kind=rp)                  :: relres
   real(kind=rp), allocatable     :: resnrm(:)

! Initialization for parameter routines (read command line)
   call initialize( )

! Which method?
   use_idrs          = get_parameter('-idrs', default_idrs )
   use_bicgstab      = get_parameter('-bicgstab', default_bicgstab )
   use_qmridr        = get_parameter('-qmridr', default_qmridr )
   if ( use_idrs ) then
      s = get_parameter('-idrs', default_s )
   elseif ( use_bicgstab ) then
      s = 1
   elseif ( use_qmridr ) then
      s = get_parameter('-qmridr', default_s )
   end if
   if ( use_idrs )       variant = 1
   if ( use_bicgstab )   variant = 2
   if ( use_qmridr )     variant = 0

! Other iteration parameters:
   tol = get_parameter('-tol', default_tol )
   maxit = get_parameter('-maxit', default_maxit )
   recycle = get_parameter('-recycle', default_recycle )
   n_omega = get_parameter('-omega', default_omega )
   n_ritz = get_parameter('-ritz', default_ritz )
   in_s = get_parameter('-in_s', default_in_s )
   in_tol = get_parameter('-in_tol', default_in_tol )
   in_it = get_parameter('-in_it', default_in_it )

! Polynomial preconditioner:
   degree = 0
   use_chebyshev    = get_parameter('-chebyshev', default_chebyshev )
   use_neumann      = get_parameter('-neumann', default_neumann )
   use_hessenberg   = get_parameter('-hessenberg', default_hessenberg )
   if ( use_chebyshev ) then
      degree = get_parameter('-chebyshev', default_degree )
   elseif ( use_neumann ) then
      degree = get_parameter('-neumann', default_degree )
   elseif ( use_hessenberg ) then
      degree = get_parameter('-hessenberg', default_degree )
   end if
   if ( use_neumann )    preconditioner = 1
   if ( use_chebyshev )  preconditioner = 2
   if ( use_hessenberg ) preconditioner = 3
   if ( degree == 0 )    preconditioner = 0


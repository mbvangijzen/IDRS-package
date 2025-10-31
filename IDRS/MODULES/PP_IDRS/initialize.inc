!
! Include file to declare standard IDRS parameter and to read the command line options
!
! Allowed command line options are
!
!    Option                     Description                                  Default
! ======================================================================================
!  -idrs          Use idrs as solution method.                                 Use idrs
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
!  -recycle       Create s-dimenional recycle space,                       No recycling
!  -omega         Number of user defined parameters omega                  n_omega = 0
!  -chebyshev     If specified a Chebyshev preconditioner should be used.   degree = 0
!                 Should come with the degree of the polynomial
!  -neumann       If specified a Neumann preconditioner should be used.     degree = 0
!                 Should come with the degree of the polynomial
!  -hessenberg    If specified a Hessenberg preconditioner should be used.  degree = 0
!                 Should come with the degree of the polynomial
!  -ritz          Number of Ritz values to compute parameter of pp            ritz = 32
!  -plot          Plot convergence history and/or ritz values(needs gnuplot)  plot = none
!                 Options are: ritz, conv or all. Plots to screen.
!  -gif           Same as -plot, but writes to gif-file
!  -jpeg          Same as -plot, but writes to jpeg-file
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2024 Martin van Gijzen
!

! Method:
   logical                        :: use_idrs, use_qmridr, use_bicgstab, use_minsync

! Iteration parameters:
   integer                        :: s
   real(kind=rp)                  :: tol
   integer                        :: maxit
   logical                        :: recycle
   integer                        :: in_it
   integer                        :: in_s
   real(kind=rp)                  :: in_tol
   integer                        :: n_omega
   integer                        :: default_omega = 0

! Preconditioner:
   logical                        :: use_chebyshev  
   logical                        :: use_neumann   
   logical                        :: use_hessenberg
   integer                        :: degree

! Defaults for the method:
   logical, parameter             :: default_idrs       = .false.
   logical, parameter             :: default_qmridr     = .false.
   logical, parameter             :: default_bicgstab   = .false.
   logical, parameter             :: default_minsync    = .false.

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

! System matrix and preconditioner (or Mass matrix):
   type(matrix)                   :: A, M1

! Parameters to compute polynomial preconditioner
   real(kind=rp)                  :: rcenter    = 1._rp
   real(kind=rp)                  :: rfoci(2)   = (/0.15_rp, 1.85_rp/)
   complex(kind=cp)               :: ccenter    = 1._cp
   complex(kind=cp)               :: cfoci(2)   = (/0.15_rp, 1.85_rp/)

! Plotting:
   character(len=:), allocatable  :: plot
   character(len=4), parameter    :: default_plot = 'none'
   logical                        :: plot_ritz, plot_conv
   logical                        :: gif, jpeg
   logical, parameter             :: default_gif  = .false.
   logical, parameter             :: default_jpeg = .false.

! Timings:
   integer                        :: tb, te, clock_rate, clock_max, t0, t1

! Paralllel computing:
   integer                        :: my_proc, n_procs

! Other declarations:
   integer                        :: preconditioner
   integer                        :: variant
   integer                        :: flag
   integer                        :: iter, ini_iter
   integer                        :: i_conv
   integer                        :: i_ritz, n_ritz
   integer, parameter             :: default_ritz = 32
   complex(kind=cp), allocatable  :: ritzval(:)
   integer                        :: converged
   logical, allocatable           :: outlier(:)
   real(kind=rp)                  :: relres
   real(kind=rp), allocatable     :: resnrm(:), resvec(:), ms_resvec(:,:)
   real(kind=rp)                  :: normb, tol_def
   character(len=160)             :: title =''
   character(len=80)              :: method, pol_prec
   character(len=20)              :: filename
   logical                        :: initial_iterations = .false.
   logical                        :: estimate_parameters
   logical                        :: multishift = .false.
   logical                        :: shift_invert = .false.

! get the command line
   call initialize()

! Which processor?
   my_proc = this_image()
! How many processors?
   n_procs = num_images()

! Which method?
   use_idrs          = get_parameter('-idrs', default_idrs )
   use_bicgstab      = get_parameter('-bicgstab', default_bicgstab )
   use_qmridr        = get_parameter('-qmridr', default_qmridr )
   use_minsync       = get_parameter('-minsync', default_minsync )

   if ( use_idrs ) then
      s = get_parameter('-idrs', default_s )
   elseif ( use_bicgstab ) then
      s = 1
   elseif ( use_qmridr ) then
      s = get_parameter('-qmridr', default_s )
   elseif ( use_minsync ) then
      s = n_procs
   end if

! No solver specified? Use IDR(4)
   if ( (.not. use_idrs) .and. (.not. use_bicgstab) .and. (.not. use_qmridr ) .and. (.not. use_minsync ) ) then
      use_idrs = .true.
      s = 4
   end if

   if ( use_idrs )       variant = 1
   if ( use_bicgstab )   variant = 2
   if ( use_minsync )    variant = 3
   if ( use_qmridr )     variant = 0

! Other iteration parameters:
   tol = get_parameter('-tol', default_tol )
   maxit = get_parameter('-maxit', default_maxit )
   recycle = get_parameter('-recycle', default_recycle )
   n_omega = get_parameter('-omega', default_omega )
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

! Number of Ritz values?
   n_ritz = get_parameter('-ritz', default_ritz )

! Which plots?
   plot = get_parameter('-plot', default_plot )
   plot_ritz = ( plot == 'ritz' .or. plot == 'all' )
   plot_conv = ( plot == 'convergence' .or. plot == 'all' )

! Which gifs?
   gif = get_parameter('-gif', default_gif )
   if ( gif ) then
      plot = get_parameter('-gif', default_plot )
      plot_ritz = ( plot == 'ritz' .or. plot == 'all' )
      plot_conv = ( plot == 'convergence' .or. plot == 'all' )
   end if

! Which jpeg?
   jpeg = get_parameter('-jpeg', default_jpeg )
   if ( jpeg ) then
      plot = get_parameter('-jpeg', default_plot )
      plot_ritz = ( plot == 'ritz' .or. plot == 'all' )
      plot_conv = ( plot == 'convergence' .or. plot == 'all' )
   end if

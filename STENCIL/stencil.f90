program stencil

!
! This programme illustrates the use of the idrs-module on
! a convection-diffusion-reaction problem with constant coefficients
! and homogeneous dirichlet boundary conditions.
! The problem is solved with the different idrs-algorithms: IDRS and QMRIDR 
! for the single shift (= damping parameter) problem and MSIDRS and MSQMRIDR 
! for the multishift problem, which uses several reaction parameters.
!
! The matrix is not explicitly stored. The matrix-vector multiplication
! is done using a stencil operations.
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2023 Martin van Gijzen
!

   use matrix_module             ! Always needed is a user module that defines the 
                                 ! matrix-vector multiplication and the preconditioning operation
                                 ! This module should also define the parameters
                                 ! rp (kind of real precision) and cp (kind of complex precision)
   use preconditioner_module
   use idrs_module
   use general_module

   implicit none 
!
   integer, parameter            :: Nsigma = 6              ! Number of shifts
   real(kind=rp)                 :: sigma(Nsigma)           ! Shifts
   real(kind=rp)                 :: dsigma                  ! Shift spacing 
   real(kind=rp)                 :: sigma_max               ! Maximum shift
   real(kind=rp), parameter      :: default_sigma_max = -5. ! Maximum shift
   real(kind=rp)                 :: shifts(Nsigma)          ! Shifts of preconditioned system
!
! IDRS parameters
   integer, parameter            :: maxit=400 
   real(kind=rp), parameter      :: tol=1.e-6
   integer                       :: flag            
   integer                       :: iter          
   integer                       :: s            
   real(kind=rp)                 :: relres, relres_ms(Nsigma), error, norm_wmod
   integer                       :: variant

   type(matrix)                  :: A
   type(preconditioner)          :: M1
   real(kind=rp), allocatable    :: w(:), w_mod(:), b(:), y_ms(:,:), w_ms(:,:)

   integer, parameter            :: default_degree = 0
   integer                       :: degree
   real(kind=rp), allocatable    :: alpha(:)                       ! Coefficients of polynomial preconditioner

   real                          :: t, t0, t1                      ! For CPU time
   integer                       :: tb, te, clock_rate, clock_max  ! For system_clock
   integer                       :: my_proc, n_procs

   real(kind=rp)                 :: omega, r, c, f1, f2
   integer                       :: i

   my_proc = this_image()
   n_procs = num_images()

! Initialization for parameter routines (read command line)
   call initialize( )

! Make the matrix
   call define_matrix( A )

! Allocate solution and right-hand-side vectors:
   allocate( w(A%Nn), w_mod(A%Nn), b(A%Nn), y_ms(A%Nn,Nsigma), w_ms(A%Nn,Nsigma) )

! Compute right-hand-side
   call model_solution( w_mod, A )
   norm_wmod = norm( w_mod )
   b = A*w_mod

! Make the preconditioner
  degree = 0

  if ( degree == 0 ) then
     degree = get_parameter( '-chebyshev', default_degree )
     if ( degree > 0 ) then
        allocate(alpha(0:degree))
! Compute estimates for largest and smallest eigenvalueswith Gershgorin's theorem:
        c = A%c
        r = abs(A%n) + abs(A%s) + abs(A%e) + abs(A%w) + abs(A%f) +abs(A%b)
        f1 = c-r
        f2 = c+r
        alpha = chebyshev( degree, f1, f2 )
     end if
  end if

  if ( degree == 0 ) then
     degree = get_parameter( '-neumann', default_degree )
     if ( degree > 0 ) then
        allocate(alpha(0:degree))
! Compute relaxation parameter:
        omega = 1./A%c
        alpha = neumann( degree, omega )
     end if
  end if
 
! Store everything in the preconditioner structure:
  if ( degree == 0 )  then
     M1%degree = 0
  else
     allocate( M1%ralpha_p(0:degree) )
     M1%A = A
     M1%degree = degree
     M1%ralpha_p = alpha
  end if

! Solve accuracy test
   if ( my_proc == 1 ) then 
      write(*,'(a)') &
         '================================================================================================'
      write(*,'(19x,a,i8,/)') &
         'Convection-diffusion problem, size is ', n_procs*A%Nn

      write(*,'(1x,a,/,a)') &
         '   Method  |    s    | Iterations |  ||r||/||b||  |  Flag  | Relative error | Elapsed time',&
         '------------------------------------------------------------------------------------------------'
   end if
   
   variant = 1
   do i = 1,4
      s = 2**(i-1)
      call cpu_time( t0 )
      call system_clock ( tb, clock_rate, clock_max )
      w = idrs( A, b, M1, s, tol, maxit, variant, flag, relres, iter )
      call system_clock ( te, clock_rate, clock_max )
      t = real(te-tb)/real(clock_rate)
      call cpu_time( t1 )
      error = norm( w - w_mod )/norm_wmod
      if ( my_proc == 1 ) &
         write(*,'(a10, a,i3,a,i4,a,e9.3,a,i1,a,e9.3,a,e9.3)') 'IDRS', &
         '  |  ', s,'    |   ', iter,'     |   ', relres,'   |    ',flag,'   |   ', error, '    |   ', t
   end do
   
   variant = 2
   s = 1
   call cpu_time( t0 )
   call system_clock ( tb, clock_rate, clock_max )
   w = idrs( A, b, M1, s, tol, maxit, variant, flag, relres, iter )
   call system_clock ( te, clock_rate, clock_max )
   t = real(te-tb)/real(clock_rate)
   call cpu_time( t1 )
   error = norm( w - w_mod )/ norm_wmod
   if ( my_proc == 1 ) &
      write(*,'(a10, a,i3,a,i4,a,e9.3,a,i1,a,e9.3,a,e9.3)') 'BICGSTAB', &
         '  |  ', s,'    |   ', iter,'     |   ', relres,'   |    ',flag,'   |   ', error, '    |   ', t

   do i = 1,4
      s = 2**(i-1)
      call cpu_time( t0 )
      call system_clock ( tb, clock_rate, clock_max )
      w = qmridr( A, b, M1, s, tol, maxit, flag, relres, iter )
      call system_clock ( te, clock_rate, clock_max )
      t = real(te-tb)/real(clock_rate)
      call cpu_time( t1 )
      error = norm( w - w_mod )/norm_wmod
      if ( my_proc == 1 ) &
         write(*,'(a10, a,i3,a,i4,a,e9.3,a,i1,a,e9.3,a,e9.3)') 'QMRIDR', &
         '  |  ', s,'    |   ', iter,'     |   ', relres,'   |    ',flag,'   |   ', error, '    |   ', t
   end do

   if ( my_proc == 1 ) then
      write(*,'(a)') &
         '------------------------------------------------------------------------------------------------'
      write(*,'(19x,a)') 'Results for the Convection-Diffusion test problem'
      write(*,*)
      write(*,'(a)') &
         '================================================================================================'
      write(*,'(4x,a,i8,a,i2,/)') &
         'Convection-diffusion-reaction problem, size is ', n_procs*A%Nn, '. Number of shifts is ', Nsigma

      write(*,'(1x,a,/,a)') &
         '   Method  |    s    | Iterations |  ||r||/||b||  |  Flag  | Relative error | Elapsed time',&
         '------------------------------------------------------------------------------------------------'
   end if

! First define the shifts:
   sigma_max = get_parameter( '-reaction', default_sigma_max )
   dsigma = 0.
   if ( Nsigma > 1 ) dsigma = sigma_max/(Nsigma-1)
   do i = 1, Nsigma
      sigma(i) = (i-1)*dsigma
   end do

! Determine the preconditioned shifted systems:
   call make_msprecon( M1, sigma, shifts )

   variant = 1
   do i = 1,4
      s = 2**(i-1)
      call cpu_time( t0 )
      call system_clock ( tb, clock_rate, clock_max )
      y_ms = msidrs( A, b, shifts, M1, s, tol, maxit, variant, flag, relres_ms, iter )
      w_ms = scaleback_precon( M1, y_ms )
      call system_clock ( te, clock_rate, clock_max )
      t = real(te-tb)/real(clock_rate)
      call cpu_time( t1 )
      error = norm( w_mod - w_ms(:,1) )/norm_wmod
      relres = maxval( relres_ms)
      if ( my_proc == 1 ) &
         write(*,'(a10, a,i3,a,i4,a,e9.3,a,i1,a,e9.3,a,e9.3)') 'MSIDRS', &
         '  |  ', s,'    |   ', iter,'     |   ', relres,'   |    ',flag,'   |   ', error, '    |   ', t
   end do
   
   variant = 2
   s = 1
   call cpu_time( t0 )
   call system_clock ( tb, clock_rate, clock_max )
   y_ms = msidrs( A, b, shifts, M1, s, tol, maxit, variant, flag, relres_ms, iter )
   w_ms = scaleback_precon( M1, y_ms )
   call system_clock ( te, clock_rate, clock_max )
   t = real(te-tb)/real(clock_rate)
   call cpu_time( t1 )
   error = norm( w_mod - w_ms(:,1) )/norm_wmod
   relres = maxval( relres_ms)
   if ( my_proc == 1 ) &
      write(*,'(a10, a,i3,a,i4,a,e9.3,a,i1,a,e9.3,a,e9.3)') 'MSBICGSTAB', &
      '  |  ', s,'    |   ', iter,'     |   ', relres,'   |    ',flag,'   |   ', error, '    |   ', t

   do i = 1,4
      s = 2**(i-1)
      call cpu_time( t0 )
      call system_clock ( tb, clock_rate, clock_max )
      y_ms = msqmridr( A, b, shifts, M1, s, tol, maxit, flag, relres_ms, iter )
      w_ms = scaleback_precon( M1, y_ms )
      call system_clock ( te, clock_rate, clock_max )
      t = real(te-tb)/real(clock_rate)
      call cpu_time( t1 )
      error = norm( w_mod - w_ms(:,1) )/norm_wmod
      relres = maxval( relres_ms)
      if ( my_proc == 1 ) &
         write(*,'(a10, a,i3,a,i4,a,e9.3,a,i1,a,e9.3,a,e9.3)') 'MSQMRIDR', &
         '  |  ', s,'    |   ', iter,'     |   ', relres,'   |    ',flag,'   |   ', error, '    |   ', t
   end do
   if ( my_proc == 1 ) then
      write(*,'(a)') &
         '------------------------------------------------------------------------------------------------'
      write(*,'(18x,a)') 'Results for the Convection-diffusion-reaction problem'
      write(*,*)
   end if

end program stencil

subroutine get_subdomains( A )
! 
! Determine problem size, subdomain devision and number of unknowns per processor
! The subdomain division is based on the number of processors
! Not for every number of processors a subdomain division is given
! At this moment, the maximum number of subdomains is 48 and multiples of
! 1, 2, 3, 4, or 5 subdomains per direction are allowed
! Any desired division can be added to the list
! The information about the problem size and data decomposition is stored in the 
! matrix structure.
!
   use matrix_module
   use general_module

   type(matrix), intent(out)     :: A

   integer, parameter            :: Mn = 60             ! Minumum number of nodes per direction
   integer, parameter            :: default_refine = 1  ! Default problems size (60*1)^3 nodes
   integer                       :: refine

   integer                       :: Sx, Sy, Sz, me(3)
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
      case default
         stop 'Number of processors does not match subdomain division'
   end select

!
! Determine local subdomain coordinates:
   allocate(SD[Sx,Sy,*])
   me = this_image(SD)

!
! Determine local number of subdomain nodes:
   refine = get_parameter( '-refine', default_refine ) ! Refinement level
   Nx = (Mn*refine)/Sx                                 ! Nodes per subdomain in x_direction
   Ny = (Mn*refine)/Sy                                 ! Nodes per subdomain in y_direction
   Nz = (Mn*refine)/Sz                                 ! Nodes per subdomain in z_direction
   Nn = Nx*Ny*Nz                                       ! Total number of nodes per subdomain

!
! Store the subdomain information in the matrix structure
   A%Nx = Nx
   A%Ny = Ny
   A%Nz = Nz
   A%Nn = Nn
   A%Sx = Sx
   A%Sy = Sy
   A%Sz = Sz
   A%me = me

end subroutine get_subdomains

subroutine define_matrix( A ) 
!
! Make the matrix
!
   use matrix_module
   use general_module

   type(matrix)                  :: A

! Default problem, parameters can be modified through command line
   real(kind=rp), parameter      :: default_eps = 1.e-2              ! Diffusion 
   real(kind=rp), parameter      :: default_beta_x = 0./sqrt(5.)     ! Convection in x-direction
   real(kind=rp), parameter      :: default_beta_y = 1./sqrt(5.)     ! Convection in y-direction
   real(kind=rp), parameter      :: default_beta_z = 2./sqrt(5.)     ! Convection in z-direction
   real(kind=rp), parameter      :: default_r = -5.

! Actual parameters:
   real(kind=rp)                 :: eps, beta_x, beta_y, beta_z

! Coordinates
   real(kind=rp)                 :: x, y, z

! Determine the mapping of the subdomains onto the processors
   call get_subdomains( A )

! Grid sizes
   A%hx = 1._rp/(A%Nx*A%Sx+1)
   A%hy = 1._rp/(A%Ny*A%Sy+1)
   A%hz = 1._rp/(A%Nz*A%Sz+1)

! Get parameters:
   eps    = get_parameter( '-eps', default_eps )
   beta_x = get_parameter( '-beta_x', default_beta_x )
   beta_y = get_parameter( '-beta_y', default_beta_y )
   beta_z = get_parameter( '-beta_z', default_beta_z )
   r      = get_parameter( '-reaction', default_r )

! Store coefficients in matrix:
   A%c = 2.*eps/A%hx**2 + 2.*eps/A%hy**2 + 2.*eps/A%hz**2 + r
   A%e = -eps/A%hx**2 + beta_x/(2.*A%hx)
   A%w = -eps/A%hx**2 - beta_x/(2.*A%hx)
   A%n = -eps/A%hy**2 + beta_y/(2.*A%hy)
   A%s = -eps/A%hy**2 - beta_y/(2.*A%hy)
   A%b = -eps/A%hz**2 + beta_z/(2.*A%hz)
   A%f = -eps/A%hz**2 - beta_z/(2.*A%hz)

end subroutine define_matrix

subroutine model_solution( w, A ) 
!
! Create model solution w_model = xyz(1-x)(1-y)(1-z)
!
   use matrix_module
   use general_module

   type(matrix), intent(in)      :: A
   real(kind=rp), intent(out)    :: w(A%Nn)
   real(kind=rp)                 :: w_model(A%Nx,A%Ny,A%Nz)
   integer                       :: i, j, k

! Coordinates
   real(kind=rp)                 :: x, y, z

! Model solution
   do concurrent (i=1:A%Nx,j=1:A%Ny,k=1:A%Nz)
       x = (A%me(1)-1)*A%Nx*A%hx+i*A%hx
       y = (A%me(2)-1)*A%Ny*A%hy+j*A%hy
       z = (A%me(3)-1)*A%Nz*A%hz+k*A%hz
       w_model(i,j,k) = x*y*z*(1.-x)*(1.-y)*(1.-z)
   end do
   w = reshape(w_model,(/A%Nn/))

end subroutine model_solution


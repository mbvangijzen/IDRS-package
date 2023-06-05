program mm_msidrs

!
! This programme tests the correct behaviour of the idrs-routines msidrs and msqmridr
! to solve shifted systems (A-sI) x = b on matrix-market matrices
!
! Parameters are passed through the command line options as parameter name/value combinations:
!    ./idrs -name1 value1 -name2 value2 etc
!
! Allowed command line options are:
!
!    Option                     Description                                  Default
! ======================================================================================
!  -matrix        Filename of matrix-market matrix. The filename                 ""
!                 must include the path and the extension .mtx.
!                 This parameter is required. The program checks
!                 if a corresponding rhs-file ending on _b.mtx exists.
!                 The _b.mtx should contain only one rhs-vectors.
!                 If such a file does not exist b=1 is taken as rhs.
!                 The shifts must be specified in a corresponding 
!                 file ending on _s.mtx. This file must exist.
!  -msidrs        Use msidrs as solution method. Should come with
!                 integer parameter s > 1                                       s=4
!  -msbicgstab    Use msbicgstab as solution method                      Do not use bicgstab
!  -msqmridr      Use msqmridr as solution method. Should come with
!                 integer parameter s > 1                                       s=4
!  -maxit         Maximum number of iterations                              maxit=10000
!  -tol           Tolerance, stop if ||r||/||b|| < tol                        tol=1e-6
!  -in_s          Parameter s for inner iterations, used if msqmridr         s_in=4
!                 is used as nested method, with msidrs as inner method
!  -in_it         Maximum number of inner iterations (only for msqmridr)    in_it=0
!  -in_tol        Tolerance for inner iterations                           in_tol=1e-6
!  -save          Save solution vectors
!                 The solution vectors (one for every shift)
!                 are saved in a matching matrix market file 
!                 with extension _x.mtx
!  -conv          Save convergence history.                                 not saved
!                 The convergence history is saved in a matching
!                 matrix market file with extension _c.mtx
!
! Note that graphical user intefaces to specify the above parameters
! can be found in the MATLAB directory.
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2023 Martin van Gijzen
!

   use general_module
   use matrix_module
   use preconditioner_module
   use idrs_module

   implicit none 

! Timings
   real                           :: t0, t1

! Matrices, number of equations etc.:
   type(matrix)                   :: A, RHS, SHIFT
   type(preconditioner)           :: M1
   character(len=:), allocatable  :: matrix_name
   character(len=:), allocatable  :: matrix_file, default_matrix_file, rhs_file, &
                                     shift_file, solution_file, convergence_file
   logical                        :: rhs_exist, shift_exist
   real(kind=rp), allocatable     :: xr(:,:), br(:), shiftr(:)
   complex(kind=cp), allocatable  :: xc(:,:), bc(:), shiftc(:)
   integer                        :: neq, n_shift, ierr
   logical                        :: iscomplex

! Solution method:
   logical                        :: use_msidrs, default_msidrs, &
                                     use_msqmridr, default_msqmridr, &
                                     use_msbicgstab, default_msbicgstab
! IDRS parameters:
   integer                        :: maxit, default_maxit
   integer                        :: s, default_s
   real(kind=rp)                  :: tol, default_tol
   character(len=10)              :: variant
   integer                        :: method
   real(kind=rp)                  :: in_tol, default_in_tol
   integer                        :: in_it, default_in_it
   integer                        :: in_s, default_in_s
   real(kind=rp), allocatable     :: relres(:)
   real(kind=rp), allocatable     :: res(:,:)
   integer                        :: flag
   integer                        :: iter
   integer                        :: k

! Save solution?
   logical                        :: save_sol, default_save_sol
! Save convergence?
   logical                        :: save_conv, default_save_conv
   integer                        :: rows, cols, nnz
   integer, allocatable           :: ival(:), indx(:), jndx(:)
   double precision, allocatable  :: rval(:)
   complex, allocatable           :: cval(:)
   character(len=10)              :: rep
   character(len=7)               :: field
   character(len=19)              :: symm

! Initialization for parameter routines (read command line)
   call initialize( )

! Get name of the data-file and allocate space
   allocate(character(parameter_length) ::matrix_file)
   allocate(character(parameter_length) ::matrix_name)
   allocate(character(parameter_length) ::rhs_file)
   allocate(character(parameter_length) ::shift_file)

   default_matrix_file = ''
   matrix_file = trim(get_parameter( '-matrix', default_matrix_file ))
   k = index( matrix_file, ".mtx" )
   matrix_name = matrix_file(1:k-1)
   rhs_file    = matrix_name // '_b.mtx'
   shift_file  = matrix_name // '_s.mtx'

! Read the system matrix:
   open(unit = 10, file = matrix_file, status = 'old', iostat = ierr )
   if ( ierr > 0 ) stop 'Error while opening matrix file '
   call read_mm_matrix( 10, A )
   close(10)
   neq = A%rows

! Read the right-hand-side matrix:
   inquire( file=rhs_file, exist=rhs_exist )
   if ( rhs_exist ) then 
      open(unit = 20, file = rhs_file, status = 'old', iostat = ierr )
      call read_mm_matrix( 20, RHS )
      close(20)
   end if

   n_shift = 0
! Read the shifts:
   inquire( file=shift_file, exist=shift_exist )
   if ( .not. shift_exist ) stop 'File with shifts does not exist '
   open(unit = 30, file = shift_file, status = 'old', iostat = ierr )
   call read_mm_matrix( 30, SHIFT )
   close(30)
   n_shift = shift%cols

   iscomplex = ( A%field == 'complex' ) 
   if ( shift_exist ) then
      if ( shift%field == 'complex' ) iscomplex = .true.
   end if

! Allocate space for the solution vectors, shifts and rhs-vector
   if ( iscomplex ) then
      allocate( bc(neq), xc(neq,n_shift), shiftc(n_shift) )
   else
      allocate( br(neq), xr(neq,n_shift), shiftr(n_shift) )
   end if
   allocate( relres(n_shift) )
!
! Copy the rhs-vector b from the user-defined type RHS
   if ( rhs_exist ) then
      if ( iscomplex .and. RHS%field == 'complex' ) then
         bc = RHS%cmat(:,1)
      elseif ( iscomplex .and. RHS%field == 'real' ) then
         bc = RHS%rmat(:,1)
      else
         br = RHS%rmat(:,1)
      end if
   else
! Generate the rhs-vector equal to one
      if ( iscomplex ) then
         bc = (1.0,0.0)
      else
         br = 1.0
      end if
   end if
!
! Copy the copy the shifts from  user defined variable SHIFT
   if ( iscomplex .and. SHIFT%field == 'complex' ) then
      shiftc = SHIFT%cmat(1,:)
   elseif ( iscomplex .and. SHIFT%field == 'real' ) then
      shiftc = SHIFT%rmat(1,:)
   else
      shiftr = SHIFT%rmat(1,:)
   end if

! No preconditioning
   M1%which_preconditioner = 'none'
 
! Which method?
   default_msidrs      = .false.
   default_msqmridr    = .false.
   default_msbicgstab  = .false.

   use_msidrs          = get_parameter('-msidrs', default_msidrs )
   use_msqmridr        = get_parameter('-msqmridr', default_msqmridr )
   use_msbicgstab      = get_parameter('-msbicgstab', default_msbicgstab )
   if ( use_msidrs ) then
      variant = 'msidrs'
      method = 1
      default_s = 4
      s = get_parameter('-msidrs', default_s )
   elseif ( use_msbicgstab ) then
      variant = 'msbicgstab'
      method = 2
      s = 1
   elseif ( use_msqmridr ) then
      variant = 'msqmridr'
      default_s = 4
      s = get_parameter('-msqmridr', default_s )
   else
      stop 'Error: no solution method specified!'
   end if

! Read the parameters
   default_maxit = 10000
   maxit = get_parameter('-maxit', default_maxit )
   default_tol = 1e-6
   tol = get_parameter('-tol', default_tol )
   default_in_s = 4
   in_s = get_parameter('-in_s', default_in_s )
   default_in_it = 0
   in_it = get_parameter('-in_it', default_in_it )
   default_in_tol = 1e-1
   in_tol = get_parameter('-in_tol', default_in_tol )

! Allocate for vector with residual norms
   allocate( res(maxit+1,n_shift) )

   write(*,'(a,a)') 'Testproblem = ', matrix_file
   write(*,'(a,a,a,i3)') 'Solution method is ',variant,' with s = ', s

! Solve the system
   call cpu_time( t0 )
   if ( use_msidrs ) then
      if ( iscomplex ) then
         xc = msidrs( A, bc, shiftc, M1, s, tol, maxit, method, flag, relres, iter, resvec = res )
      else
         xr = msidrs( A, br, shiftr, M1, s, tol, maxit, method, flag, relres, iter, resvec = res )
      end if
   elseif ( use_msbicgstab ) then
      if ( iscomplex ) then
         xc = msidrs( A, bc, shiftc, M1, s, tol, maxit, method, flag, relres, iter, resvec = res )
      else
         xr = msidrs( A, br, shiftr, M1, s, tol, maxit, method, flag, relres, iter, resvec = res )
      end if
   elseif ( use_msqmridr ) then
      if ( A%field == 'complex' ) then
         xc = msqmridr( A, bc, shiftc, M1, s, tol, maxit, flag, relres, iter, &
                        in_s, in_tol, in_it, resvec = res )
      else
         xr = msqmridr( A, br, shiftr, M1, s, tol, maxit, flag, relres, iter, &
                        in_s, in_tol, in_it, resvec = res )
      end if
   end if

! Output: 
   write(*,*)
   write(*,'(a,i5)') 'Number of iterations    = ', iter
   if ( iscomplex ) then
      do k = 1, n_shift
         write(*,'(a,f7.2,a,f7.2,a,e11.4)') & 
            'Relative residual norm for shift (', real(shiftc(k)), ',', aimag(shiftc(k)), ') = ', relres(k)
      end do
   else
      do k = 1, n_shift
         write(*,'(a,f9.2,a,e11.4)') & 
            'Relative residual for shift ', shiftr(k), ' = ', relres(k)
      end do
   end if
   if ( flag > 0 ) then
      if ( flag == 1 ) write(*,'(a)') 'Maximum number of iterations reached!'
      if ( flag == 2 ) write(*,'(a)') 'Accuracy above prescribed tolerance!'
      if ( flag == 3 ) write(*,'(a)') 'Break down!'
   end if
   call cpu_time( t1 )

! Output:
   write(*,*)
   write(*,'(a,f7.2,a)') 'CPU time                = ', t1-t0,'s.'
   write(*,'(a)') '=================================================='

! Save solution?
   default_save_sol = .false.
   save_sol = get_parameter('-save', default_save_sol)
   if ( save_sol ) then
      if ( iscomplex ) then
         allocate( cval(neq*n_shift) )
!
! allocations unused variables:
         allocate(ival(1),rval(1),indx(1),jndx(1))
!
! copy solutions to complex array
         do k = 1, n_shift
            cval((k-1)*neq+1:k*neq) = xc(:,k)
         end do
         field = 'complex'
      else
         allocate( rval(neq*n_shift) )
!
! allocations unused variables:
         allocate(ival(1),cval(1),indx(1),jndx(1))
! copy solutions to real array
         do k = 1, n_shift
            rval((k-1)*neq+1:k*neq) = xr(:,k)
         end do
         field = 'real'
      end if
      symm = 'general'
      rep  = 'array'
      rows = neq
      cols = n_shift
      nnz  = rows*cols
      allocate(character(parameter_length) ::solution_file)
      solution_file = matrix_name // '_x.mtx'
      open(unit = 40, file = solution_file, status = 'unknown', iostat = ierr )
      call mmwrite( 40, rep, field, symm, rows, cols, nnz, indx, jndx,ival,rval,cval)
! Clean up
      deallocate(ival,rval,cval,indx,jndx)
   end if

! Save convergence?
   default_save_conv = .false.
   save_conv = get_parameter('-conv', default_save_conv)
   if ( save_conv ) then
      allocate( rval((iter+1)*n_shift) )
!
! allocations unused variables:
      allocate(ival(1),cval(1),indx(1),jndx(1))
!
! copy residual norms to rval
      do k = 1, n_shift
         rval((k-1)*(iter+1)+1:k*(iter+1)) = res(1:iter+1,k)
      end do
      symm = 'general'
      field = 'real'
      rep = 'array'
      rows = iter+1
      cols = n_shift
      nnz = (iter+1)*n_shift
      allocate(character(parameter_length) ::convergence_file)
      convergence_file = matrix_name // '_c.mtx'
      open(unit = 50, file = convergence_file, status = 'unknown', iostat = ierr )
      call mmwrite( 50, rep, field, symm, rows, cols, nnz, indx, jndx,ival,rval,cval)
! Clean up
      deallocate(ival,rval,cval,indx,jndx)
   end if

end program mm_msidrs

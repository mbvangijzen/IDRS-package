program mm_idrs

!
! This programme tests the correct behaviour of the idrs-routines idrs and qmridr
! on matrix-market matrices
!
! Parameters are passed through the command line options as parameter name/value combinations: 
!    ./idrs -name1 value1 -name2 value2 etc
! 
! Allowed command line options are
!
!    Option                     Description                                  Default
! ======================================================================================
!  -matrix        Filename of matrix-market matrix. The filename                 ""
!                 must include the path and the extension .mtx.             
!                 This parameter is required. The program checks 
!                 if a corresponding rhs-file ending on _b.mtx exists.
!                 The _b.mtx may contain multiple rhs-vectors.
!                 If such a file does not exist b=1 is taken as rhs. 
!  -precon        If specified a Jacobi preconditioner will be used.     No preconditioner
!  -idrs          Use idrs as solution method. Should come with
!                 integer parameter s > 1                                       s=4
!  -bicgstab      Use bicgstab as solution method                        Do not use bicgstab
!  -qmridr        Use qmridr as solution method. Should come with
!                 integer parameter s > 1                                       s=4
!  -maxit         Maximum number of iterations                              maxit=10000
!  -tol           Tolerance, stop if ||r||/||b|| < tol                        tol=1e-6
!  -in_s          Parameter s for inner iterations, used if qmridr           s_in=4
!                 is used as nested method, with idrs as inner method
!  -in_it         Maximum number of inner iterations (only for qmridr)      in_it=0
!  -in_tol        Tolerance for inner iterations                           in_tol=1e-6
!  -recycle       Recycle the previous s solution,                         No recycling
!                 for multiple rhs-problems
!  -save          Save solution vectors
!                 The solution vectors (one for every rhs-vector)
!                 are saved in a matching matrix market file 
!                 with extension _x.mtx
!  -conv          Save convergence history.                                 not saved
!                 The convergence history is saved in a matching
!                 matrix market file with extension _c.mtx
!  -hessenberg    Save hessenberg matrix of size (nritz+1)xnritz            nritz=30
!                 The eigenvalues of the upper nritzxnritz part
!                 are approximations to the eigenvalues of the
!                 preconditioned system matrix. Can be used for analysis.
!                 The hessenberg matrix is saved in a matching
!                 matrix market file with extension _h.mtx
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

! Hardware parameters
   real                           :: t0, t1

! Matrices, number of equations etc.:
   type(matrix)                   :: A, RHS
   type(preconditioner)           :: M1
   character(len=:), allocatable  :: matrix_name
   character(len=:), allocatable  :: matrix_file, default_matrix_file, rhs_file

   logical                        :: rhs_exist
   real(kind=rp), allocatable     :: xr(:,:), xr0(:), br(:)
   complex(kind=cp), allocatable  :: xc(:,:), xc0(:), bc(:)
   integer                        :: neq, i_rhs, n_rhs, ierr
   logical                        :: iscomplex

! Solution method:
   logical                        :: use_idrs, default_idrs, &
                                     use_qmridr, default_qmridr, &
                                     use_bicgstab, default_bicgstab, &
                                     use_precon, default_precon
! IDRS parameters:
   integer                        :: maxit, default_maxit, numit
   integer                        :: s, default_s
   real(kind=rp)                  :: tol, default_tol
   character(len=10)              :: variant
   integer                        :: method
   real(kind=rp)                  :: in_tol, default_in_tol
   integer                        :: in_it, default_in_it
   integer                        :: in_s, default_in_s
   real(kind=rp)                  :: relres
   real(kind=rp), allocatable     :: res(:,:)
   integer                        :: flag
   integer                        :: iter
   integer                        :: i,j,k

! Recycling:
   logical                        :: recycle, default_recycle
   real(kind=rp), allocatable     :: UR0(:,:)
   complex(kind=cp), allocatable  :: UC0(:,:)

! Save solution?
   logical                        :: save_sol, default_save_sol
   character(len=:), allocatable  :: solution_file

! Save convergence?
   logical                        :: save_conv, default_save_conv
   character(len=:), allocatable  :: convergence_file

! Save Hessenberg matrix?
   logical                        :: save_hessenberg
   integer                        :: nritz, default_nritz
   complex(kind=cp), allocatable  :: chessenberg(:,:)
   real(kind=rp), allocatable     :: rhessenberg(:,:)
   character(len=:), allocatable  :: hessenberg_file

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

   default_matrix_file = ''
   matrix_file = trim(get_parameter( '-matrix', default_matrix_file ))
   k = index( matrix_file, ".mtx" )
   matrix_name = matrix_file(1:k-1)
   rhs_file    = matrix_name // '_b.mtx'

! Read the system matrix:
   open(unit = 10, file = matrix_file, status = 'old', iostat = ierr )
   if ( ierr > 0 ) stop 'Error while opening matrix file '
   call read_mm_matrix( 10, A )
   close(10)
   neq = A%rows

! Are the rhs-vectors given?
   inquire( file=rhs_file, exist=rhs_exist )
   if ( rhs_exist ) then 
      open(unit = 20, file = rhs_file, status = 'old', iostat = ierr )
      call read_mm_matrix( 20, RHS )
      close(20)
      n_rhs = RHS%cols
   else
      n_rhs = 1
   end if

   iscomplex = ( A%field == 'complex' ) 

! Allocate space for the solution vectors, shifts and rhs-vector
   if ( iscomplex ) then
      allocate( xc(neq,n_rhs), xc0(neq), bc(neq) )
      xc0 = (0.,0)
   else
      allocate( xr(neq,n_rhs), xr0(neq), br(neq) )
      xr0 = 0.
   end if
!
! Construct the (diagonal) preconditioner:
   default_precon = .false.
   use_precon     = get_parameter('-precon', default_precon )
   if ( use_precon ) then 
      call make_preconditioner( A, M1 )
   else
      M1%which_preconditioner = 'none'
   end if

! Which method?
   default_idrs      = .false.
   default_qmridr    = .false.
   default_bicgstab  = .false.
   use_idrs          = get_parameter('-idrs', default_idrs )
   use_bicgstab      = get_parameter('-bicgstab', default_bicgstab )
   use_qmridr        = get_parameter('-qmridr', default_qmridr )
   if ( use_idrs ) then
      variant = 'idrs'
      method = 1
      default_s = 4
      s = get_parameter('-idrs', default_s )
   elseif ( use_bicgstab ) then
      variant = 'bicgstab'
      method = 2
      s = 1
   elseif ( use_qmridr ) then
      variant = 'qmridr'
      default_s = 4
      s = get_parameter('-qmridr', default_s )
   else
      stop 'Error: no solution method specified!'
   end if

! Read the iteration parameters
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
   allocate( res(maxit+1,n_rhs) )
 
! Recycling?
   default_recycle = .false.
   recycle         = get_parameter('-recycle', default_recycle )
   if ( recycle ) then
      if ( iscomplex ) then
         allocate( UC0(neq,s) )
      else
         allocate( UR0(neq,s) )
      end if
   end if

! Save the Hessenberg matrix (in IDRS)?
   default_nritz = 0
   nritz = get_parameter('-hessenberg', default_nritz )
   save_hessenberg = ( nritz > 0 )
   if ( iscomplex ) then
      allocate( chessenberg(nritz+1,nritz) )
   else
      allocate( rhessenberg(nritz+1,nritz) )
   end if

! Numit monitors the maximum number of iterations to solve the systems
   numit = 0

   write(*,'(a,a)') 'Testproblem = ', matrix_file
   write(*,'(a,a,a,i3)') 'Solution method is ',variant,' with s = ', s

! Solve the system
   call cpu_time( t0 )
   do i_rhs = 1,n_rhs
!
! Copy the rhs-vector b from the user-defined type RHS
      if ( rhs_exist ) then
         if ( iscomplex .and. RHS%field == 'complex' ) then
            bc = RHS%cmat(:,i_rhs)
         elseif ( iscomplex .and. RHS%field == 'real' ) then
            bc = RHS%rmat(:,i_rhs)
         else
            br = RHS%rmat(:,i_rhs)
         end if
      else
! Generate random rhs-vector 
         if ( iscomplex ) then
            bc = 1.
         else
            br = 1.
         end if
      end if
! Now solve the system
      if ( use_idrs .or. use_bicgstab ) then
         if ( recycle .and. ( i_rhs > s ) ) then
            if ( iscomplex ) then
               UC0 = xc(:,i_rhs-s:i_rhs-1)
               xc(:,i_rhs) = idrs( A, bc, M1, s, tol, maxit, method, flag, relres, iter, & 
                  xc0, UC0, resvec = res(:,i_rhs), H = chessenberg )
               xc0 = xc(:,i_rhs)
            else
               UR0 = xr(:,i_rhs-s:i_rhs-1)
               xr(:,i_rhs) = idrs( A, br, M1, s, tol, maxit, method, flag, relres, iter,  &
                  xr0, UR0, resvec = res(:,i_rhs), H = rhessenberg )
               xr0 = xr(:,i_rhs)
            end if
         else
            if ( iscomplex ) then
               xc(:,i_rhs) = idrs( A, bc, M1, s, tol, maxit, method, flag, relres, iter, xc0, & 
                  resvec = res(:,i_rhs), H = chessenberg )
               xc0 = xc(:,i_rhs)
            else
               xr(:,i_rhs) = idrs( A, br, M1, s, tol, maxit, method, flag, relres, iter, xr0, &
                  resvec = res(:,i_rhs), H = rhessenberg )
               xr0 = xr(:,i_rhs)
            end if
         end if
      elseif ( use_qmridr ) then
         if ( A%field == 'complex' ) then
            xc(:,i_rhs) = qmridr( A, bc, M1, s, tol, maxit, & 
                          flag, relres, iter, in_s, in_tol, in_it, &
                          xc0, resvec = res(:,i_rhs) )
            xc0 = xc(:,i_rhs)
         else
            xr(:,i_rhs) = qmridr( A, br, M1, s, tol, maxit, &
                          flag, relres, iter, in_s, in_tol, in_it, &
                          xr0, resvec = res(:,i_rhs) )
            xr0 = xr(:,i_rhs)
         end if
      end if
      numit = max( numit, iter )
      write(*,*)
      write(*,'(a,i2)') 'Results for system ', i_rhs
      write(*,'(a,i5)') 'Number of iterations    = ', iter
      write(*,'(a,e11.4)') 'Relative residual norm  = ', relres
      if ( flag > 0 ) then
         if ( flag == 1 ) write(*,'(a)') 'Maximum number of iterations reached!'
         if ( flag == 2 ) write(*,'(a)') 'Accuracy above prescribed tolerance!'
         if ( flag == 3 ) write(*,'(a)') 'Break down!'
      end if
   end do
   call cpu_time( t1 )

! Output: 
   write(*,*)
   write(*,'(a,f7.2,a)') 'CPU time                = ', t1-t0, 's.'
   write(*,'(a)') '=================================================='

! Save solution?
   default_save_sol = .false.
   save_sol = get_parameter('-save', default_save_sol)
   if ( save_sol ) then
      if ( iscomplex ) then
         allocate( cval(neq*n_rhs) )
!
! allocations unused variables:
         allocate(ival(1),rval(1),indx(1),jndx(1))
!
! copy solutions to complex array
         k = 0
         do i = 1, n_rhs
            k = k + 1
            cval((k-1)*neq+1:k*neq) = xc(:,i)
         end do
         field = 'complex'
      else
         allocate( rval(neq*n_rhs) )
!
! allocations unused variables:
         allocate(ival(1),cval(1),indx(1),jndx(1))
! copy solutions to real array
         k = 0
         do i = 1, n_rhs
            k = k + 1
            rval((k-1)*neq+1:k*neq) = xr(:,i)
         end do
         field = 'real'
      end if
      symm = 'general'
      rep  = 'array'
      rows = neq
      cols = n_rhs
      nnz  = rows*cols
      allocate(character(parameter_length) ::solution_file)
      solution_file = matrix_name // '_x.mtx'
      open(unit = 40, file = solution_file, status = 'unknown', iostat = ierr )
      call mmwrite( 40, rep, field, symm, rows, cols, nnz, indx, jndx,ival,rval,cval)
! Clean up
      deallocate(ival,rval,cval,indx,jndx)
   end if

   if ( save_hessenberg) then
      if ( iscomplex ) then
         allocate( cval((nritz+1)*nritz) )
!
! allocations unused variables:
         allocate(ival(1),rval(1),indx(1),jndx(1))
!
! copy Hessenberg matrix to complex array
         k = 0
         do i = 1, nritz
            k = k + 1
            cval((k-1)*(nritz+1)+1:k*(nritz+1)) = chessenberg(:,i)
         end do
         field = 'complex'
      else
         allocate( rval((nritz+1)*nritz) )
!
! allocations unused variables:
         allocate(ival(1),cval(1),indx(1),jndx(1))
! copy solutions to real array
         k = 0
         do i = 1, nritz
            k = k + 1
            rval((k-1)*(nritz+1)+1:k*(nritz+1)) = rhessenberg(:,i)
         end do
         field = 'real'
      end if
      symm = 'general'
      rep  = 'array'
      rows = nritz+1
      cols = nritz
      nnz  = rows*cols
      allocate(character(parameter_length) ::hessenberg_file)
      hessenberg_file = matrix_name // '_h.mtx'
      open(unit = 50, file = hessenberg_file, status = 'unknown', iostat = ierr )
      call mmwrite( 50, rep, field, symm, rows, cols, nnz, indx, jndx,ival,rval,cval)
! Clean up
      deallocate(ival,rval,cval,indx,jndx)
   end if

! Save convergence?
   default_save_conv = .false.
   save_conv = get_parameter('-conv', default_save_conv)
   if ( save_conv ) then
      allocate( rval( (numit+1)*n_rhs ) )
      rval = 0.d0
!
! allocations unused variables:
      allocate(ival(1),cval(1),indx(1),jndx(1))
!
! copy residual norms to rval
      do k = 1, n_rhs
         rval((k-1)*(numit+1)+1:k*(numit+1)) = res(1:numit+1,k)
      end do
      symm  = 'general'
      field = 'real'
      rep   = 'array'
      rows  = numit+1
      cols  = n_rhs
      nnz   = rows*cols
      allocate(character(parameter_length) ::convergence_file)
      convergence_file = matrix_name // '_c.mtx'
      open(unit = 60, file = convergence_file, status = 'unknown', iostat = ierr )
      call mmwrite( 60, rep, field, symm, rows, cols, nnz, indx, jndx,ival,rval,cval)
! Clean up
      deallocate(ival,rval,cval,indx,jndx)
   end if

end program mm_idrs

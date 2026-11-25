!
! Include file to store the user matrices in the matrix type 
! The system matrix and the preconditioner are stored in the
! matrix A. The preconditioner is composed of:
!    - shift-and-invert matrix
!    - polynomial preconditioner
!    - deflation matrix
! These operators can be combined.
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2026 Martin van Gijzen
!

!
!---------------------------------------------------------------------
! System matrix:
!---------------------------------------------------------------------
!

! Put stiffness matrix in matrix type
   A%K = K

!
! Put mass matrix in matrix type (if present)
   if ( present(M) ) then
      A%mass_matrix = .true.
      A%M = M
   end if

!
! Put damping matrix in matrix type (if present)
   if ( present(C) ) then
      A%damping_matrix = .true.
      A%C = C
   end if

! Multishift problem?
   A%multishift = multishift

! Do we have a block system?
   A%block_system = ( ( multishift .or. eigen ) .and. present(C) )

   A%cseed = cseed

!
!---------------------------------------------------------------------
! Shift-and-invert preconditioner (or decomposition  of mass matrix)):
!---------------------------------------------------------------------
!

   A%cshift = cshift

   if ( present(P) ) then
      A%P = P
      if ( present(shift) ) then
         shift_invert = .true.
      else
         shift_invert = .false.
      end if
   else
      call system_clock ( tb, clock_rate, clock_max )
      if ( use_lu ) then
! Shift-and-invert preconditioner
         shift_invert = .true.
         if ( present(C) .and. present(M) ) then
            A%P = complex_decomp( K=K, M=M, C=C, shift = cshift )
         elseif ( present(C) ) then
            A%P = complex_decomp( K=K, C=C, shift = cshift )
         elseif ( present(M) ) then
            A%P = complex_decomp( K=K, M=M, shift = cshift )
         else
            A%P = complex_decomp( K=K, shift = cshift )
         end if
         call system_clock ( te, clock_rate, clock_max )
         if ( my_proc == 1 ) &
            write(*,'(a,f8.2,a,/)') 'Time for LU decomposition = ', real(te-tb)/real(clock_rate), 's.'
      elseif ( ( multishift .or. eigen ) .and. present(M) ) then
! Decompose mass matrix 
         A%P = complex_decomp( M )
         shift_invert = .false.
         call system_clock ( te, clock_rate, clock_max )
         if ( my_proc == 1 ) &
            write(*,'(a,f8.2,a,/)') 'Time for decomposition of Mass matrix = ', real(te-tb)/real(clock_rate), 's.'
      else
! No SI preconditioner
         shift_invert = .false.
      end if
   end if
   A%shift_invert = shift_invert

!
!-------------------------------------------------------
! Polynomial preconditioner:
!-------------------------------------------------------
!

! Remove preconditioner (if it was included)
   A%preconditioner = 0
   A%degree = 0
   if ( allocated(A%ralpha_p) ) deallocate(A%ralpha_p)
   if ( allocated(A%ralpha_s) ) deallocate(A%ralpha_s)

   if ( present(center) ) ccenter = center
   if ( present(foci) )   cfoci   = foci

! Estimation of parameters, or are they given?
   estimate_parameters = .true.
   if ( preconditioner == 0 ) estimate_parameters = .false.
   if ( preconditioner == 1 .and. present(center) ) estimate_parameters = .false.
   if ( preconditioner == 2 .and. present(foci)   ) estimate_parameters = .false.

   if ( estimate_parameters ) then
      H = IDR_HESSENBERG( A, b, s, n_ritz, M1, variant, flag )
      call complex_preconditioner( A, preconditioner, degree, H = H(1:n_ritz,1:n_ritz) )
   else
      call complex_preconditioner( A, preconditioner, degree, cfoci, ccenter )
   end if 

!
!-------------------------------------------------------
! Deflation matrix:
!-------------------------------------------------------
!

   if ( present(V) ) then
      call make_deflation( A, V )
   elseif ( n_deflate > 0 ) then
      allocate(eigval(n_deflate), U(neq,n_deflate), resnrm(n_deflate) )

      call system_clock ( tb, clock_rate, clock_max )
      eigval = idreig( A, b, s, n_deflate, M1, tol, n_ritz, variant, U, flag, &
                       resnrm, iter )

      call system_clock ( te, clock_rate, clock_max )

! Report performance measures:
      if ( my_proc == 1 ) then
         write(*,'(a)') 'Results for IDREIG:'
         write(*,'(a,f8.2,a)') 'Elapsed time           = ', real( te - tb )/real( clock_rate ), 's.'
         write(*,'(a,i5)')     'Number of iterations   = ', iter
         if ( flag > 0 ) then
            if ( flag == 1 ) write(*,'(a)') 'Maximum number of iterations reached!'
            if ( flag == 2 ) write(*,'(a)') 'Accuracy above prescribed tolerance!'
            if ( flag == 3 ) write(*,'(a)') 'Break down!'
         end if
      end if

      if ( .not. flag == 3 ) then
         if ( my_proc == 1 ) then
            write(*,*)
            write(*,'(a)') '      Eigenvalue       |    Residual norm '
            write(*,'(a)') '----------------------------------------------'
            do i_deflate = 1, n_deflate
               write(*,'(1x,a,e9.3,a,e9.3,a,1x,a,6x,e10.4)') &
                '(', real(eigval(i_deflate)),',',aimag(eigval(i_deflate)), ')','|',resnrm(i_deflate)
            end do
            write(*,'(a)') '-----------------------------------------------'
            write(*,*)
         end if
         call make_deflation( A, U )
      end if
      deallocate( eigval, U, resnrm )
   end if

!
!-------------------------------------------------------
! Precomputed omega:
!-------------------------------------------------------
!

   if ( n_omega > 0 ) then
      allocate( omega(n_omega) )
      H = IDR_HESSENBERG( A, b, s, n_ritz, M1, variant, flag )
      omega = compute_omega( H, n_ritz, n_omega )
      if ( flag > 0 ) then
         write(*,'(a)') 'Error in pre-computation of omega, using default omega.'
         deallocate( omega )
      end if
   end if

!
!-------------------------------------------------------
! Plot the Ritz values
!-------------------------------------------------------
!

   if ( plot_ritz ) then
! Compute the ritz values:
      allocate( ritzval(n_ritz) )
      H = IDR_HESSENBERG( A, b, s, n_ritz, M1, variant, flag )
      ritzval = HQR( H(1:n_ritz,1:n_ritz) )

      if (  my_proc == 1 ) then

! Make plot file:
         write(title,'(a)') 'Ritz values'
         open( 20, file = 'ritz.plt' )
         if ( gif ) then
            write(20,*) 'set terminal gif '
            write(20,*) 'set output "ritz.gif"' 
         elseif ( jpeg ) then
            write(20,*) 'set terminal jpeg'
            write(20,*) 'set output "ritz.jpeg"'
         end if
         write(20,*) 'set xlabel "Real"'
         write(20,*) 'set ylabel "Imag"'
         write(20,*) 'set grid'
         write(20,*) 'unset key'
         write(20,*) 'set title "', trim(title), '"'
         write(20,*) 'plot "ritz.dat" with points'
         close(20)

! Make data file:
         open( 30, file = 'ritz.dat' )
         do i_ritz = 1, n_ritz
            write(30,*) real(ritzval(i_ritz)), aimag(ritzval(i_ritz))
         end do
         close(30)

! Run:
         call execute_command_line('gnuplot -p ritz.plt')
      end if
   end if


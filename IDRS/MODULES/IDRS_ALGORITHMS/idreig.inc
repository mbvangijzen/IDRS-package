   real(kind=rp)                           :: tol
   integer                                 :: N                 ! number of equations

   integer                                 :: info              ! convergence indicator
   logical                                 :: mass_matrix       ! mass matrix exists?
   logical                                 :: out_flag          ! store flag

   logical                                 :: user_shift        ! user defined shifts present
   integer                                 :: n_shift           ! number of user defined shifts

   logical                                 :: out_iterations    ! store number of iterations
   logical                                 :: out_resvec        ! store residual norms

   integer                                 :: ii                ! inner iterations index
   integer                                 :: jj                ! G-space index
   real(kind=rp)                           :: normr, normg
   integer                                 :: i,j,k,l           ! loop counters

   integer                                 :: method
   integer                                 :: maxit

   real(kind=rp)                           :: err(neig)

! Number of equations:
   N = size(b)

! Starting vector:
   normr = norm(b)
   r  = b/normr

! Check optional input parameters:
   mass_matrix = present(M1)

   if ( present(tolerance) ) then
      if ( tolerance < 0 ) error stop "Illegal value for parameter tolerance"
      tol = tolerance
   else
      tol = epsilon(1._rp)*1.e1
   endif

   if ( present(maximum_iterations) ) maxit = maximum_iterations

   method = 1                                 ! Random shadow vectors
   if ( present(variant) ) method = variant

! Initialize the output variables 
   out_flag       = present(flag)
   if ( out_flag )       flag = -1 
   out_iterations = present(iterations)

! Generate shadow space:
   if ( method == 1 ) then
! Standard IDREIG
      allocate( PT(s,n) )
   elseif ( method == 2 ) then
! BiLanczos
      allocate( PT(1,n) )
   elseif ( method == 3 ) then
! IDREIG, Sparse vectors, one per processor
      allocate( PT(1,n) )
   else
! BiLanczos, sparse vectors, one per processor
      allocate( PT(1,n) )
   endif
   call MAKE_P(PT, b, method)

   user_shift = present(omega)
   if ( user_shift ) then
      n_shift = size(omega)
      if ( n_shift < 1 ) user_shift = .false.
   end if

   allocate( H(maxit+1,maxit) )

! Initialize Hessenberg matrix:
   H = 0.

! Initialize local variables:
   M = 0.
   shift = 0.
   info = -1
   jj = 0
   ii = 0
   method = 1
    
! This concludes the initialisation phase

! Main iteration loop, build G-spaces:
    
   do while (  info < 0 )  ! start of iteration loop
     
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Generate s vectors in G_j
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

! New right-hand side for small system:
      f  = P_DOT( PT, r, s, method )
      mu = f

      do k=1, s

! Compute new v
         v = r 
         if ( jj > 0 ) then

! Solve small system (Note: M is lower triangular) and make v orthogonal to P:
            do i = k,s
               gamma(i) = mu(i)
               do j = k, i-1
                  gamma(i) = gamma(i) - M(i,j)*gamma(j)
               end do
               gamma(i) = gamma(i)/M(i,i)
               v = v - gamma(i)*G(:,i)
            end do
         end if

! Store old r in G-space:
         G(:,k)   = r
         M(k:s,k) = mu(k:s)

! Break down?
         if ( abs(M(k,k)) <= tiny(1._rp) ) then
            info = 3
            exit
         end if

! Compute new vector in G-space:
         if ( mass_matrix ) then
            r = v/M1
         else
            r = v
         end if
         t = A*r
         r = t - shift*v

! Bi-Orthogonalise the new basis vectors: 
         f = P_DOT( PT, r, s, method )
         mu = f
         do i = 1,k
            alpha(i) = f(i)
            do j = 1, i-1
               alpha(i) = alpha(i) - M(i,j)*alpha(j)
            end do
            alpha(i) = alpha(i)/M(i,i)
            r        = r - G(:,i)*alpha(i)
            mu(k:s)  = mu(k:s)  - M(k:s,i)*alpha(i)
         end do

         normg = norm(r);
! Break down?
         if ( normg <= tiny(1._rp) ) then
            info = 3
         end if
         r = r/normg;
         mu(k:s) = mu(k:s)/normg

!
! Store the information:
         ii = ii + 1;
         H(ii-k+1:ii,ii) = alpha(1:k)
         H(ii+1,ii) = normg
         H(ii,ii) = H(ii,ii) + shift
         if ( jj > 0 ) H(ii-s:ii-k,ii)=gamma(k:s)*normr

! Are we ready?
         if ( ii == maxit ) then
            info = 1
            exit
         end if

      end do ! Now we have computed s+1 vectors in G_j
      if ( info > 0 ) exit

!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Compute first residual in G_j+1
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Update G-space counter
      jj = jj + 1

! Compute first residual in G_j+1
! Note: r is already perpendicular to P so v = r
 
! Preconditioning:
      if ( mass_matrix ) then
         v = r/M1
      else
         v = r
      end if
      t = A*v

! Computation of a new shift
      if ( user_shift ) then
         i = mod(jj,n_shift)
         if ( i == 0 ) i = n_shift
         shift = omega(i)
      else
         shift = inprod( r, t )
      end if

! Update residual
      r = t - shift*r 
      normr = norm(r)
      r = r/normr

!
! Check for convergence 
      if ( ii > neig ) then 
         eigval = compute_eigval( H(1:ii,1:ii), neig, err )

! Convergence?
         if ( all( err < tol ) ) info = 0
      end if

   end do ! end of while loop

  if ( present(eigvec) ) then
! Compute the eigen vectors:
      if ( mass_matrix .and. present(resvec) ) then
         eigvec(:,1:neig) =  compute_eigvec( A, s, H(1:ii+1,1:ii), eigval, b, M1, tol, resvec=resvec )
      elseif ( mass_matrix ) then
         eigvec(:,1:neig) =  compute_eigvec( A, s, H(1:ii+1,1:ii), eigval, b, M1, tol )
      elseif ( present(resvec) ) then
         eigvec(:,1:neig) =  compute_eigvec( A, s, H(1:ii+1,1:ii), eigval, b, tolerance=tol, resvec=resvec )
      else
         eigvec(:,1:neig) =  compute_eigvec( A, s, H(1:ii+1,1:ii), eigval, b, tolerance=tol )
      end if

! Compute the residual norms:
      if ( present(resnrm) ) then
         do j = 1,neig
            if ( mass_matrix ) then
               resnrm(j) = norm( A*(eigvec(:,j)/M1) - eigval(j)*eigvec(:,j) )
            else
               resnrm(j) = norm( A*eigvec(:,j) - eigval(j)*eigvec(:,j) )
            end if
         end do
         if (( info == 0 ) .and. ( .not. all( resnrm < tol*abs(eigval) ) )) info = 2
      end if
   end if
  
   if ( out_flag )          flag = info
   if ( out_iterations )    iterations = ii


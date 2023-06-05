! Declarations:
   integer               :: n                 ! dimension of the system
   integer               :: maxit             ! maximum number of iterations
   integer               :: method            ! which IDR(s) variant?
   real(kind=rp)         :: tol               ! actual tolerance
   integer               :: info              ! convergence indicator
   logical               :: out_flag          ! output flag
   logical               :: out_relres        ! output relres
   logical               :: out_iterations    ! output number of iterations
   logical               :: out_colfac        ! output colinearity factors
   logical               :: user_omega        ! input user defined omega
   integer               :: n_omega           ! number of user defined omega's
   logical               :: out_resvec        ! output residual norms

   integer               :: iter              ! number of iterations
   integer               :: jj                ! G-space index
   real(kind=rp)         :: normb, tolb       ! for tolerance check
   real(kind=rp)         :: normg, normr(size(sigma))
   integer               :: k_sigma, n_sigma
   integer               :: i,j,k             ! loop counters

! Problem size:
   n    = size(b)
! Number of shifts
   n_sigma = size(sigma)
      
! Check optional input parameters:
   if ( present(tolerance) ) then
      if ( tolerance < 0 ) stop "Illegal value for parameter tolerance"
      tol = tolerance 
   else
      tol = 1e-6
   endif

   maxit=min(2*n,1000)
   if ( present(maximum_iterations) ) maxit = maximum_iterations 
   
   method = 1 ! msidrs   
   if ( present(variant) ) method = variant

! Initialize the output variables 
   out_flag       = present(flag)
   if ( out_flag )       flag = -1 
   out_relres     = present(relres)
   if ( out_relres)      relres = 1. 
   out_iterations = present(iterations)
   if ( out_iterations ) iterations = 0 
   out_colfac = present(colfac)

   x = 0.

   user_omega = present(omega)
   if ( user_omega ) then
      n_omega = size(omega)
   end if

! Check output arrays
   out_resvec     = present(resvec)
   if ( out_resvec ) then
      if ( maxit+1 > size(resvec) ) &
        stop "Length of vector with residual norms too small, should be maxit+1"
   end if

! Set tolerance
   normb = norm( b )
   tolb = tol * normb
   gn = b
   normg = normb
   normr = normg
   if ( out_resvec ) resvec(1,:)= normr

! check if the initial solution is not already a solution within the prescribed
! tolerance
   if (normg <= tolb) then      
      if ( out_iterations ) iterations = 0               
      if ( out_flag )       flag  = 0
      if ( out_relres )     relres = normr/normb
      return
   end if

! Define P and kappa (depending on the method)
   if ( method == 1 ) then
! Standard MSIDRS
      allocate( PT(s,n) )
      kappa = 0.7
   elseif ( method == 2 ) then
! MSBiCGSTAB
      allocate( PT(1,n) )
      kappa = 0.
   endif
   call MAKE_P(PT, gn, method)

! Initialize local variables:
   M = 0.
   G  = 0.
   W  = 0.
   v0 = 0.
   om = sqrt(2.)
   om_sigma = om/(1 - om*sigma)

! Last nonzero coefficients of projected rhs:
   phi = 1.

! Pivots of LU decomposition:
   l_sigma = 0.

   iter = 0
   info = -1
   jj = 0
    
! This concludes the initialisation phase

! Main iteration loop, build G-spaces:
    
   do while (  info < 0 )  ! start of iteration loop
     
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Generate s vectors in G_j
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

! New right-hand side for small system:
      f = P_DOT( PT, gn, s, method )

      do k=1, s
         iter = iter + 1

! Compute new v
         v = gn 
         if ( jj > 0 ) then

! Solve small system (Note: M is lower triangular) and make v orthogonal to P:
            do i = k,s
               gamma(i) = f(i)
               do j = k, i-1
                  gamma(i) = gamma(i) - M(i,j)*gamma(j)
               end do
               gamma(i) = gamma(i)/M(i,i)
               v = v - gamma(i)*G(:,i)
            end do

         end if

! Update M
         M(k:s,k) = f(k:s);
! Break down?
         if ( abs(M(k,k)) <= tiny(tol) ) then
            info = 3
            exit
         end if
! Store old gn
         G(:,k) = gn

!
! Compute new gn = (I - om A) v
         t = A*(v/M1)
         gn = v - om*t
        
! Bi-Orthogonalise the new basis vectors: 
         f = P_DOT( PT, gn, s, method )
         alpha  = 0.
         do i = 1,k
            alpha(i) = f(i)
            do j = 1, i-1
               alpha(i) = alpha(i) - M(i,j)*alpha(j)
            end do
            alpha(i) = alpha(i)/M(i,i)
            gn = gn - G(:,i)*alpha(i)
            f(k+1:s)  = f(k+1:s)  - M(k+1:s,i)*alpha(i)
         end do
         normg = norm( gn )

!
! Now form the next column of the Hessenberg matrix
         h = 0.
         if ( jj > 0 ) h(1:s-k+1) = -gamma(k:s)
         h(s+2-k+1:s+2) = h(s+2-k+1:s+2)-alpha(1:k)
         h(s+2) = h(s+2)+1.
         h(s+3) = -1.
         h = h/om;
         do k_sigma = 1, n_sigma
            u_sigma = h(1:s+2)
            if ( jj > 0 ) then
               u_sigma(1:s-k+1) = u_sigma(1:s-k+1) + sigma(k_sigma)*gamma(k:s)
            end if
            u_sigma(s+2) = u_sigma(s+2)-sigma(k_sigma)

            l_sigma(1:s+1,k_sigma) = l_sigma(2:s+2,k_sigma)
            do j = 2,s+2
               u_sigma(j) = u_sigma(j) - l_sigma(j-1,k_sigma)*u_sigma(j-1);
            end do
!
! Next pivot, update solution of projected system
            l_sigma(s+2,k_sigma) = -1/(om*u_sigma(s+2))

!
! New Update vector:
! Corresponds to Dimension Reduction step
            t = v - om_sigma(k_sigma)*v0*u_sigma(s+2-k)
! Corresponds to space j
            if ( k > 1 ) then
               do i = 1,k-1
                  t = t - W(:,i,k_sigma)*u_sigma(s+2-k+i)
               end do
            end if
! Corresponds to space j+1
            do i = k, s
               t = t - W(:,i,k_sigma)*u_sigma(i-k+1)
            end do
! Scale and store
            W(:,k,k_sigma) = t/u_sigma(s+2)
!
! Compute solution:
            x(:,k_sigma) = x(:,k_sigma) + phi(k_sigma)*W(:,k,k_sigma)
! Update rhs of projected system:
            phi(k_sigma) =  -l_sigma(s+2,k_sigma)*phi(k_sigma)
! Compute residual norms
            normr(k_sigma) = normg*abs(phi(k_sigma))
      
         end do
         if ( out_resvec ) resvec(iter+1,:) = normr
         if ( maxval(normr) < tolb ) then
            info = 0
            exit
         elseif ( iter == maxit ) then
            info = 1
            exit
         end if
      end do

!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Compute first residual in G_j+1
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if ( info >= 0 ) exit

! Update G-space counter
      jj = jj + 1
      iter = iter + 1

! Compute first residual in G_j+1
 
      v = gn
      t = A*(v/M1)

! Computation of a new omega
      if ( user_omega ) then
         i = mod(jj,n_omega)
         if ( i == 0 ) i = n_omega
         om = omega(i)
      else
         om = comp_om( t, v, kappa )
      end if
      if ( abs(om) <= epsilon(tol) ) then 
         info = 3
         exit
      end if 
      gn = v - om*t

      normg = norm( gn )
!
!
! It is not necessary to form a new colum of the Hessenberg matrix, it is simply
!        h_sigma(s+2) = 1/om
!        h_sigma(s+3) = -1/om
! For u_sigma we therefore get
!        u_sigma(s+2) = 1/om - sigma;
!
      do k_sigma = 1,n_sigma

! Compute pivots
         l_sigma(1:s+1,k_sigma) = l_sigma(2:s+2,k_sigma)
         l_sigma(s+2,k_sigma)   = -1/(1-om*sigma(k_sigma))

! Compute u_sigma(s+2), these are omega's for the shifted systems:
         om_sigma(k_sigma) = om/(1-om*sigma(k_sigma))
!
! Update vector. Note that this one is special:
! It is collinear for all the shifts, with collinearity factor om_sigma
         v0 = v
!
! Compute solution:
         x(:,k_sigma) = x(:,k_sigma) + (phi(k_sigma)*om_sigma(k_sigma))*v0
! Update rhs of projected system:
         phi(k_sigma) = -l_sigma(s+2,k_sigma)*phi(k_sigma)
! Compute residual norms
         normr(k_sigma) = normg*abs(phi(k_sigma))
      end do
      if ( out_resvec ) resvec(iter+1,:) = normr
      if ( maxval(normr) < tolb ) then
         info = 0
         exit
      elseif ( iter == maxit ) then
         info = 1
         exit
      end if

   end do ! end of while loop

   do k_sigma = 1,n_sigma
      gn = x(:,k_sigma)/M1
      t = b - A*gn + sigma(k_sigma)*x(:,k_sigma)
      normr(k_sigma) = norm( t )
   end do

   if ( info == 0 .and. maxval(normr) >= tolb ) info = 2
   if ( out_iterations ) iterations = iter
   if ( out_relres )     relres=normr/normb
   if ( out_flag )       flag = info
   if ( out_colfac) colfac = phi


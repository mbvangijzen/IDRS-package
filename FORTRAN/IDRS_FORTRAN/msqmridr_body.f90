! Declarations:
   integer               :: n                 ! dimension of the system
   integer               :: maxit             ! maximum number of iterations
   integer               :: in_s              ! parameter s for inner iterations
   integer               :: in_it             ! maximum number of inner iterations
   integer               :: method            ! which QMRIDR(s) variant?
   real(kind=rp)         :: tol               ! actual tolerance
   real(kind=rp)         :: in_tol            ! inner tolerance
   integer               :: info              ! convergence indicator
   logical               :: mass              ! Mass matrix exists?
   logical               :: out_flag          ! store flag
   logical               :: out_relres        ! store relres
   logical               :: out_iterations    ! store number of iterations
   logical               :: user_omega        ! user defined omega present
   integer               :: n_omega           ! number of user defined omega's
   logical               :: out_resvec        ! store residual norms

   logical               :: converged(size(sigma))
   integer               :: iter              ! number of iterations
   integer               :: jj                ! G-space index
   real(kind=rp)         :: normb, tolb       ! for tolerance check
   real(kind=rp)         :: normg, normr(size(sigma))
   integer               :: k_sigma, low, n_sigma
   integer               :: i,j,k,l           ! loop counters

! Problem size:
   n    = size(b)
! Number of shifts
   n_sigma = size(sigma)

! Mass matrix?
   mass = present(M1)

! Check optional input parameters:
   if ( present(tolerance) ) then
      if ( tolerance < 0 ) stop "Illegal value for parameter tolerance"
      tol = tolerance
   else
      tol = 1e-6
   end if

   maxit=min(2*n,1000)
   if ( present(maximum_iterations) ) maxit = maximum_iterations

   in_s=0
   if ( present(inner_s) ) in_s = inner_s
   if ( in_s < 0 ) stop "Illegal value for parameter inner_s"

   in_tol = 1.e-1
   if ( present(inner_tolerance) ) in_tol = inner_tolerance
   if ( in_tol < 0 ) stop "Illegal value for parameter inner_tolerance"

   in_it=0
   if ( present(inner_iterations) ) in_it = inner_iterations

   method = 1 ! msqmridr

! Initialize the output variables
   out_flag       = present(flag)
   if ( out_flag )  flag = -1
   out_relres     = present(relres)
   if ( out_relres) relres = 1.
   out_iterations = present(iterations)
   if ( out_iterations ) iterations = 0

   x = 0.

   user_omega = present(omega)
   if ( user_omega ) then
      n_omega = size(omega)
      if ( n_omega < 1 ) user_omega = .false.
   end if

! Check output arrays
   out_resvec     = present(resvec)
   if ( out_resvec ) then
      if ( maxit+1 > size(resvec,1) ) &
        stop "Length of vector with residual norms too small, should be maxit+1"
   end if

! compute initial residual, set absolute tolerance
   normb = norm( b )
   tolb = tol * normb
   gn = b
   normg = normb
   normr = normg
   if (out_resvec ) resvec(1,:) = normr

! check if the initial solution is not already a solution within the prescribed
! tolerance
   if (normg <= tolb) then
      if ( out_iterations ) iterations = 0
      if ( out_flag )       flag  = 0
      if ( out_relres )     relres = normg/normb
      return
   end if

! Define P and kappa (depending on the method)
   if ( method == 1 ) then
! Standard QMRIDR
      allocate( PT(s,n) )
      kappa = 0.7
   elseif ( method == 2 ) then
! QMRSTAB
      allocate( PT(1,n) )
      kappa = 0.
   endif
   call MAKE_P(PT, gn, method)

   converged = .false.
   M  = 0.
   G  = 0.
   gn = gn/normg
   W  = 0.
   wn = 0.
   mu = 0.
   om = sqrt(2.)

! Last two nonzero coefficients of projected rhs:
   phi_n = 0.
   phi_n1 = normg

   iter = 0
   cs = 0.
   sn = 0.

   jj = 0
   info = -1

   do while ( info < 0 )
   
      do k = 1,s+1
!
! Update counters:
         iter = iter + 1;
!
! First phase: make new vector in G-space:
         c = 0.d0
         c(s+1) = 1.d0
         mu = P_DOT( PT, gn, s, method )
         if ( iter > s ) then
!
! Construct v orthogonal to P
! Solve small system (Note: M is lower triangular) and make v orthogonal to P:
            gamma = solve( M, mu )
            v = gn
            do i = 1,s
               v = v - gamma(i)*G(:,i)
            end do
            c(1:s) = -gamma
         else  
!
! First s steps: Arnoldi
            v = gn;
         end if
         M(:,1:s-1) = M(:,2:s)
         M(:,s) = mu
         G(:,1:s-1) = G(:,2:s) 
         G(:,s) = gn 
!
! Inner iterations?
         if ( in_it == 0 ) then
            vtilde = v
         else
            sigma_0(1:n_sigma) = sigma
            sigma_0(n_sigma+1) = 0.
            if ( mass ) then
               V_tilde = MSIDRS( A, v,sigma_0, in_s, M1, in_tol, in_it, colfac = cf ) 
            else
               V_tilde = MSIDRS( A, v,sigma_0, in_s, &
                         tolerance = in_tol, maximum_iterations = in_it, colfac = cf ) 
            end if
            eta = cf(1:n_sigma)/cf(n_sigma+1)
            vtilde = V_tilde(:,n_sigma+1)
         end if
!
! Compute new vector in space G_j
         if ( mass ) then
            gn = A*(vtilde/M1)
         else
            gn = A*vtilde
         end if
!
! New G-space? 
         if ( k == s+1 ) then
! Computation of a new omega
            if ( user_omega ) then
               i = mod(jj,n_omega)
               if ( i == 0 ) i = n_omega
               om = omega(i)
            else
               om = comp_om( gn, v, kappa )
            end if
            if ( abs(om) <= epsilon(tol) ) then
               info = 3
               exit
            end if
            jj = jj+1
         end if
         gn = gn - v/om
!
! Orthogonalisation (modified Gram-Schmidt)
         h = c/om
         if ( k < s+1 ) then
            do i = s-k+1, s
               beta(i) = inprod(G(:,i),gn)
               gn      = gn - G(:,i)* beta(i)
               h(i+1) = h(i+1)+beta(i)
            end do
         end if
! Normalise
         normg = norm(gn)
         gn = gn/normg
         h(s+2) = normg
!
! Store the information
         do k_sigma = 1,n_sigma
            if ( .not. converged(k_sigma) ) then 
               r_sigma = 0.
               if ( in_it == 0 ) then
                  r_sigma(2:s+3) = h - sigma(k_sigma)*c
               else
                  r_sigma(2:s+3) = eta(k_sigma)*h - (eta(k_sigma)-1.)*c
               end if
!
               low = max(1,s+3-iter)
               do l = low,s+1
!
! Apply Givens rotation.
                  tau = r_sigma(l)
                  r_sigma(l)   =                       cs(l,k_sigma) *tau + sn(l,k_sigma)*r_sigma(l+1)
                  r_sigma(l+1) = -conjg(cmplx(sn(l,k_sigma),kind=cp))*tau + cs(l,k_sigma)*r_sigma(l+1)
               end do
!
! Form i-th rotation matrix.
               tau = r_sigma(s+2)
               call grot( tau, r_sigma(s+3), cs(s+2,k_sigma), sn(s+2,k_sigma), r_sigma(s+2) )
!
! Update projected right-hand side
               phi_n(k_sigma)    =                       cs(s+2,k_sigma) *phi_n1(k_sigma)
               phi_n1(k_sigma)   = -conjg(cmplx(sn(s+2,k_sigma),kind=cp))*phi_n1(k_sigma)
               cs(1:s+1,k_sigma) = cs(2:s+2,k_sigma)
               sn(1:s+1,k_sigma) = sn(2:s+2,k_sigma)
!
! Update vector:
               if ( abs(r_sigma(s+2)) < epsilon(tol) ) then
                  info = 3
                  exit 
               end if
               if ( in_it > 0 ) vtilde = V_tilde(:,k_sigma)
               wn = (vtilde - matmul(W(:,1:s+1,k_sigma),r_sigma(1:s+1)))/r_sigma(s+2)
               W(:,1:s,k_sigma) = W(:,2:s+1,k_sigma) 
               W(:,s+1,k_sigma) = wn 
!
! Compute solution:
               x(:,k_sigma) = x(:,k_sigma) + phi_n(k_sigma)*wn
               normr(k_sigma) = abs(phi_n1(k_sigma))*sqrt(real(jj+1))
               converged(k_sigma) = ( normr(k_sigma) < tolb)
            end if
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
   end do

! Compute the residual norms:
   do k_sigma = 1,n_sigma
      if ( mass ) then
         gn = b - A*(x(:,k_sigma)/M1) + sigma(k_sigma)*x(:,k_sigma)
      else
         gn = b - A*x(:,k_sigma) + sigma(k_sigma)*x(:,k_sigma)
      end if
      normr(k_sigma) = norm( gn )
   end do

! Scale back the solution
   if ( mass ) then
      do k_sigma = 1,n_sigma
         x(:,k_sigma) = x(:,k_sigma)/M1
      end do
   end if

   if ( info == 0 .and. maxval(normr) >= tolb ) info = 2
   if ( out_iterations ) iterations = iter
   if ( out_relres )     relres=normr/normb
   if ( out_flag )       flag = info


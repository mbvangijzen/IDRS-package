! Declarations:
   integer               :: n                 ! dimension of the system
   integer               :: maxit             ! maximum number of iterations
   integer               :: method            ! which QMRIDR(s) variant?
   real(kind=rp)         :: tol               ! actual tolerance
   integer               :: in_s              ! s of inner method
   real(kind=rp)         :: in_tol            ! inner tolerance
   integer               :: in_it             ! maximum number of inner iterations
   integer               :: info              ! convergence indicator
   logical               :: out_flag          ! store flag
   logical               :: out_relres        ! store relres
   logical               :: out_iterations    ! store number of iterations
   logical               :: inispace          ! initial search space
   logical               :: user_omega        ! user defined omega present
   integer               :: n_omega           ! number of user defined omega's
   logical               :: out_resvec        ! store residual norms

   integer               :: iter              ! number of iterations
   integer               :: ii                ! inner iterations index
   integer               :: jj                ! G-space index
   real(kind=rp)         :: normb, tolb       ! for tolerance check
   real(kind=rp)         :: normg, normr
   integer               :: low
   integer               :: i,j,k,l           ! loop counters

! Problem size:
   n    = size(b)

! Check optional input parameters:
   if ( present(tolerance) ) then
      if ( tolerance < 0 ) stop "Illegal value for parameter tolerance"
      tol = tolerance
   else
      tol = 1e-6
   end if

   maxit=min(2*n,1000)
   if ( present(maximum_iterations) ) maxit = maximum_iterations

   in_s = 4
   if ( present(inner_s) ) in_s = inner_s
   if ( in_s < 1 ) stop "Illegal value for parameter inner_s"

   in_tol = 1.e-1
   if ( present(inner_tolerance) ) in_tol = inner_tolerance
   if ( in_tol < 0 ) stop "Illegal value for parameter inner_tolerance"

   in_it=0
   if ( present(inner_iterations) ) in_it = inner_iterations

   method = 1 ! qmridr    

! Initialize the output variables
   out_flag       = present(flag)
   if ( out_flag )  flag = -1
   out_relres     = present(relres)
   if ( out_relres) relres = 1.
   out_iterations = present(iterations)
   if ( out_iterations ) iterations = 0

! Check optional input arrays:
   x = 0.
   if ( present(x0) ) x = x0

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

! compute initial residual, set absolute tolerance
   normb = norm( b )
   tolb = tol * normb
   gn = b - A*x
   normg = norm( gn )
   normr = normg
   if (out_resvec ) resvec(1) = normr;

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

   M  = 0.
   G  = 0.
   gn = gn/normg
   W  = 0.
   wn = 0.
   mu = 0.
   om = 1.

! Last two nonzero coefficients of projected rhs:
   phi_n = 0.
   phi_n1 = normg

   iter = 0
   cs = 0.d0
   sn = 0.d0

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
! Preconditioning:
         if ( in_it == 0 ) then
            vtilde = v/M1
         else
            vtilde = idrs( A, v, M1, in_s, in_tol, in_it )
         end if
!
! Compute new vector in space G_j
         gn = A*vtilde;
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
         h = c/om;
         if ( k < s+1 ) then
            do i = s-k+1, s
               beta(i) = inprod(G(:,i),gn) 
               gn      = gn - G(:,i)* beta(i)
               h(i+1) = h(i+1)+beta(i)
            end do
         end if
! Normalise
         normg = norm( gn )
         gn = gn/normg
         h(s+2) = normg
!
! Store the information
         r_sigma = 0.
         r_sigma(2:s+3) = h 
!
         low = max(1,s+3-iter)
         do l = low,s+1
!
! Apply Givens rotation.
            tau = r_sigma(l)
            r_sigma(l)   =                       cs(l) *tau + sn(l)*r_sigma(l+1)
            r_sigma(l+1) = -conjg(cmplx(sn(l),kind=cp))*tau + cs(l)*r_sigma(l+1)
         end do
!
! Form i-th rotation matrix.
         tau = r_sigma(s+2)
         call grot( tau, r_sigma(s+3), cs(s+2), sn(s+2), r_sigma(s+2) )
!
! Update projected right-hand side
         phi_n     =                       cs(s+2) *phi_n1
         phi_n1    = -conjg(cmplx(sn(s+2),kind=cp))*phi_n1
         cs(1:s+1) = cs(2:s+2)
         sn(1:s+1) = sn(2:s+2)
!
! Update vector:
         if ( abs(r_sigma(s+2)) < epsilon(tol) ) then
            info = 3
            exit 
         end if
         wn = (vtilde - matmul(W(:,1:s+1),r_sigma(1:s+1)))/r_sigma(s+2)
         W(:,1:s) = W(:,2:s+1) 
         W(:,s+1) = wn 
!
! Compute solution:
         x = x + phi_n*wn
         normr = abs(phi_n1)*sqrt(real(jj+1))
         if ( out_resvec ) resvec(iter+1) = normr
         if ( normr < tolb ) then
            info = 0
            exit
         elseif ( iter == maxit ) then
            info = 1
            exit
         end if
      end do
   end do

   gn = b - A*x
   normr = norm( gn )

   if ( info == 0 .and. normr >= tolb ) info = 2
   if ( out_iterations ) iterations = iter
   if ( out_relres )     relres=normr/normb
   if ( out_flag )       flag = info


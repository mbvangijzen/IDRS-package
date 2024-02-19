! Declarations:
   integer               :: n                 ! dimension of the system
   integer               :: maxit             ! maximum number of iterations
   integer               :: method            ! which IDR(s) variant?
   real(kind=rp)         :: tol               ! actual tolerance
   integer               :: info              ! convergence indicator
   logical               :: preconditioner    ! preconditioner matrix exists?
   logical               :: out_flag          ! store flag
   logical               :: out_relres        ! store relres
   logical               :: out_iterations    ! store number of iterations
   logical               :: inispace          ! initial search space
   logical               :: user_omega        ! user defined omega present
   integer               :: n_omega           ! number of user defined omega's
   logical               :: out_resvec        ! store residual norms
   logical               :: out_H             ! store iteration parameters in H
   integer               :: nritz             ! Number of wanted ritz values

   integer               :: iter              ! number of iterations
   integer               :: ii                ! inner iterations index
   integer               :: jj                ! G-space index
   real(kind=rp)         :: normb, normr, tolb! for tolerance check
   integer               :: i,j,k,l           ! loop counters

! Problem size:
   n    = size(b)
      
! Check optional input parameters:
   preconditioner = present(M1)

   if ( present(tolerance) ) then
      if ( tolerance < 0 ) stop "Illegal value for parameter tolerance"
      tol = tolerance 
   else
      tol = 1e-8
   endif

   maxit=min(2*n,1000)
   if ( present(maximum_iterations) ) maxit = maximum_iterations 
   
   method = 1 ! idrs   
   if ( present(variant) ) method = variant

! Initialize the output variables 
   out_flag       = present(flag)
   if ( out_flag )       flag = -1 
   out_relres     = present(relres)
   if ( out_relres)      relres = 1. 
   out_iterations = present(iterations)
   if ( out_iterations ) iterations = 0 

! Check optional input arrays:
   x = 0.
   if ( present(x0) ) x = x0
      
   U = 0.
   inispace =  present(U0)
   if ( inispace ) then
      if ( .not. allocated(U0) ) then
         inispace = .false.
      elseif ( size(U0) /= size(U) ) then
         inispace = .false.
      else
         U = U0
      end if
   end if

   user_omega = present(omega)
   if ( user_omega ) then
      n_omega = size(omega)
      if ( n_omega < 1 ) user_omega = .false.
   end if

! Check output arrays
   out_resvec     = present(resvec)
   if ( out_resvec ) then
      if ( maxit+1 > size(resvec) ) &
        stop "Length of vector with residual norms too small, should be maxit+1"
   end if

! Compute Hessenberg matrix?
   out_H = present(H)
   if ( out_H ) then
      nritz = size(H,1)-1
      if ( nritz <= 0 ) out_H = .false.
   end if
   if ( out_H ) then
      if ( size(H,2) /= nritz ) &
         stop "Second dimension of H incompatible, with first"
      H = 0.
   end if

! compute initial residual, set absolute tolerance
   normb = norm( b )
   tolb = tol * normb
   r = b - A*x
   normr = norm( r )
   if ( out_resvec ) resvec(1)= normr

! check if the initial solution is not already a solution within the prescribed
! tolerance
   if (normr <= tolb) then      
      if ( out_iterations ) iterations = 0               
      if ( out_flag )       flag  = 0
      if ( out_relres )     relres = normr/normb
      return
   end if

! Define P and kappa (depending on the method)
   if ( method == 1 ) then
! Standard IDRS 
      allocate( PT(s,n) )
      kappa = 0.7
   elseif ( method == 2 ) then
! BiCGSTAB 
      allocate( PT(1,n) )
      kappa = 0.
   elseif ( method == 3 ) then
! IDRS, Sparse vectors, one per processor
      allocate( PT(1,n) )
      kappa = 0.7
   else
! BiCGSTAB, sparse vectors, one per processor 
      allocate( PT(1,n) )
      kappa = 0.
   endif
   call MAKE_P(PT, r, method)

! Initialize local variables:
   M = 0.
   om = 1.
   iter = 0
   info = -1
   jj = 0
   ii = 0
    
! This concludes the initialisation phase

! Main iteration loop, build G-spaces:
    
   do while (  info < 0 )  ! start of iteration loop
     
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Generate s vectors in G_j
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

! New right-hand side for small system:
      f = P_DOT( PT, r, s, method )

      do k=1, s

! Update inner iteration counter
         ii = ii + 1

! Compute new v
         v = r 
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

! Compute new U(:,k)
            if ( preconditioner ) then
               t = om*(v/M1)
            else
               t = om*v
            end if
            do i = k,s
               t = t + gamma(i)*U(:,i)
            end do
            U(:,k) = t

! Compute Hessenberg matrix?
            if ( out_H .and. ii <= nritz ) &
               H(ii-s:ii-k,ii)   = -gamma(k:s)/beta(k:s)

         else if ( .not. inispace ) then

! Updates for the first s iterations (in G_0):
            if ( preconditioner ) then
               U(:,k) = v/M1
            else
               U(:,k) = v
            end if

         end if

! Compute new G(:,k), G(:,k) is in space G_j
         G(:,k) = A*U(:,k)
        
! Bi-Orthogonalise the new basis vectors: 
         mu = P_DOT( PT, G(:,k), s, method )
         do i = 1,k-1
            alpha(i) = mu(i)
            do j = 1, i-1
               alpha(i) = alpha(i) - M(i,j)*alpha(j)
            end do
            alpha(i) = alpha(i)/M(i,i)
            G(:,k) = G(:,k) - G(:,i)*alpha(i)
            U(:,k) = U(:,k) - U(:,i)*alpha(i)
            mu(k:s)  = mu(k:s)  - M(k:s,i)*alpha(i)
         end do
         M(k:s,k) = mu(k:s)

! Compute Hessenberg matrix?
         if ( out_H .and. ii <= nritz .and. k  > 1 ) &
            H(ii-k+1:ii-1,ii) =  alpha(1:k-1)/beta(1:k-1)

! Break down?
         if ( abs(M(k,k)) <= tiny(tol) ) then
            info = 3
            exit
         end if

! Make r orthogonal to p_i, i = 1..k, update solution and residual 
         beta(k) = f(k)/M(k,k)
         r = r - beta(k)*G(:,k)
         x = x + beta(k)*U(:,k)

! New f = P'*r (first k  components are zero)
         if ( k < s ) then
            f(k+1:s)   = f(k+1:s) - beta(k)*M(k+1:s,k)
         end if

! Compute Hessenberg matrix?
         if ( out_H .and. ii <= nritz ) then     
            H(ii,ii) = 1./beta(k)
            l = max(1,ii-s)
            H(l+1:ii+1,ii) = H(l+1:ii+1,ii) - H(l:ii,ii)
            H(l:ii+1,ii)   = H(l:ii+1,ii)/om
         end if

! Check for convergence
         normr = norm( r )
         iter = iter + 1
         if ( out_resvec ) resvec(iter + 1) = normr
         if ( normr < tolb ) then
            info = 0
            exit
         elseif ( iter == maxit ) then
            info = 1
            exit
         end if 

      end do ! Now we have computed s+1 vectors in G_j
      if ( info >= 0 ) exit

!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Compute first residual in G_j+1
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if ( info > 0 ) exit
! Update G-space counter
      jj = jj + 1

! Compute first residual in G_j+1
! Note: r is already perpendicular to P so v = r
 
! Preconditioning:
      if ( preconditioner ) then
         v = r/M1
      else
         v = r
      end if
      t = A*v

! Computation of a new omega
      if ( user_omega ) then
         i = mod(jj,n_omega)
         if ( i == 0 ) i = n_omega
         om = omega(i)
      else
         om = comp_om( t, r, kappa )
      end if
      if ( abs(om) <= epsilon(tol) ) then 
         info = 3
         exit
      end if 

! Update solution and residual
      r = r - om*t 
! Skip the update in case of BiCG
      x = x + om*v 

! Check for convergence
      normr = norm( r )
      iter = iter + 1
      if ( out_resvec ) resvec(iter + 1) = normr
      if ( normr < tolb ) then
         info = 0
      elseif ( iter == maxit ) then
         info = 1
      end if 

   end do ! end of while loop

! Set output parameters
   r = b - A*x
   normr = norm( r )

   if ( info == 0 .and. normr > tolb ) info = 2
   if ( out_iterations ) iterations = iter
   if ( out_relres )     relres=normr/normb
   if ( out_flag )       flag = info
   

   integer                                 :: N                 ! number of equations

   integer                                 :: info              ! convergence indicator
   logical                                 :: mass_matrix       ! mass matrix exists?
   logical                                 :: out_flag          ! store flag

   logical                                 :: user_shift        ! user defined shifts present
   integer                                 :: n_shift           ! number of user defined shifts

   integer                                 :: ii                ! inner iterations index
   integer                                 :: jj                ! G-space index
   real(kind=rp)                           :: normr, normg
   integer                                 :: i,j,k,l           ! loop counters

   integer                                 :: method

! Number of equations:
   N = size(b)

! Starting vector:
   normr = norm(b)
   r  = b/normr

! Check optional input parameters:
   mass_matrix = present(M1)

   method = 1                                 ! Random shadow vectors
   if ( present(variant) ) method = variant

! Initialize the output variables 
   out_flag       = present(flag)
   if ( out_flag )       flag = -1 

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
   endif
   call MAKE_P(PT, b, method)

   user_shift = present(omega)
   if ( user_shift ) then
      n_shift = size(omega)
      if ( n_shift < 1 ) user_shift = .false.
   end if

! Initialize Hessenberg matrix:
   H = 0.

! Initialize local variables:
   M = 0.
   shift = 0.
   info = 0
   jj = 0
   ii = 0
   method = 1
    
! This concludes the initialisation phase

! Main iteration loop, build G-spaces:
    
   do while (  ii < min( nritz, N) .and. info == 0 )  ! start of iteration loop
     
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
            info = 1
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
            info = 2
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
         if ( ii == nritz ) exit

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

   end do ! end of while loop

   if ( out_flag ) flag = info


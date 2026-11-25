! Include file to solve the systems

! Create (random) rhs-vectors and solve the correct problem:
   if ( nsigma > 0 .and. nrhs > 0 ) then
! Multishift-multiple rhs problem
      if ( .not. allocated(x_msmr) ) then
         allocate(x_msmr(nrows,nsigma,nrhs),sigma(nsigma),b_mr(nrows,nrhs))
         call random_number( b_mr )
      end if
      dsigma = 1._rp/real(nsigma,kind=rp)
      do i = 1, nsigma
         sigma(i) = seed -(i-1)*dsigma
      end do
      x_msmr = idrs_driver( A, b_mr, sigma, D, seed=seed, shift=shift )

   elseif ( nsigma > 0 ) then

! Multishift problem
      if ( .not. allocated(x_ms) ) then
         allocate(x_ms(nrows,nsigma),sigma(nsigma),b(nrows))
         call random_number( b )
      end if
      dsigma = 1._rp/real(nsigma,kind=rp)
      do i = 1, nsigma
         sigma(i) = seed -(i-1)*dsigma
      end do
      x_ms = idrs_driver( A, b, sigma, D, seed=seed, shift=shift )

   elseif ( nrhs > 0 ) then

! Multiple rhs problem
      if ( .not. allocated(x_mr) ) then
         allocate(x_mr(nrows,nrhs),b_mr(nrows,nrhs))
         call random_number( b_mr )
      end if
      x_mr = idrs_driver( A, b_mr, D, seed=seed )

   else

! Linear system
      if ( .not. allocated(x) ) then
         allocate(x(nrows),b(nrows))
         call random_number( b )
      end if
      x = idrs_driver( A, b, D, seed=seed )

   end if

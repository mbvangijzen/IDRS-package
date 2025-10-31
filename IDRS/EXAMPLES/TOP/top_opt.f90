program top_opt
!
! F90 translation of the matlab script
! A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE SIGMUND, JANUARY 2000 
! CODE MODIFIED FOR INCREASED SPEED, September 2002, BY OLE SIGMUND
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2025 Martin van Gijzen 
!
   use interface_module
   use precision_module
   use user_module
   use matrix_module
   use dense_la_module
   use pp_idrs_module

   implicit none

   integer                        :: grid
   integer, parameter             :: default_grid = 1
   integer                        :: nelx, nely
   real(kind=rp)                  :: volfrac
   real(kind=rp)                  :: penal
   real(kind=rp)                  :: rmin
   integer, parameter             :: default_load = 1
   integer                        :: load

   real(kind=rp), allocatable     :: rho(:,:)        ! Relative densities
   real(kind=rp), allocatable     :: rho_old(:,:)   
   real(kind=rp), allocatable     :: dc(:,:)         ! Relative densities
   real(kind=rp)                  :: rho_max
   integer                        :: N               ! Total number of DOFs
   integer                        :: nrhs            ! Number of load cases
   logical, allocatable           :: fixed(:)        ! Fixed degrees of freedom

! Plotting:
   character(len=:), allocatable  :: plot
   character(len=4), parameter    :: default_plot = 'none'
   logical                        :: plot_solution
   logical                        :: gif, jpeg
   logical, parameter             :: default_gif  = .false.
   logical, parameter             :: default_jpeg = .false.

   integer                        :: loop
   integer, parameter             :: max_loop = 1000

   real(kind=rp)                  :: change
   integer                        :: elx, ely, n1, n2

   real(kind=rp)                  :: xe(8), ye(8)
   integer                        :: ind(8)
   real(kind=rp)                  :: c

   integer                        :: i, j, irhs

! Element matrix:
   real(kind=rp), parameter       :: E  = 1._rp
   real(kind=rp), parameter       :: nu = 0.3_rp
   real(kind=rp), parameter       :: k(8) = E/(1.-nu**2) * & 
      [ 1._rp/2._rp-nu/ 6._rp,       1._rp/8._rp+ nu/8._rp, -1._rp/4._rp-nu/12._rp, &
       -1._rp/8._rp+3._rp*nu/8._rp, -1._rp/4._rp+nu/12._rp, -1._rp/8._rp- nu/8._rp, &
                         nu/ 6._rp,  1._rp/8._rp-3._rp*nu/8._rp]
   real(kind=rp), parameter       :: Ke(8,8) = reshape(([ &
        k(1), k(2), k(3), k(4), k(5), k(6), k(7), k(8), &
        k(2), k(1), k(8), k(7), k(6), k(5), k(4), k(3), &
        k(3), k(8), k(1), k(6), k(7), k(4), k(5), k(2), &
        k(4), k(7), k(6), k(1), k(8), k(3), k(2), k(5), &
        k(5), k(6), k(7), k(8), k(1), k(2), k(3), k(4), &
        k(6), k(5), k(4), k(3), k(2), k(1), k(8), k(7), &
        k(7), k(4), k(5), k(2), k(3), k(8), k(1), k(6), &
        k(8), k(3), k(2), k(5), k(4), k(7), k(6), k(1)]), [8,8])

   type(user_matrix)             :: A, P

   real(kind=rp), allocatable    :: x(:,:), b(:,:), D(:)

   complex(kind=cp)              :: eival(8)
   real(kind=rp)                 :: eimin, eimax, rcenter, rfoci(2)

! Read command line
   call initialize()

! Make plot?
   plot = get_parameter('-plot', default_plot ) 
   plot_solution = ( plot == 'solution' .or. plot == 'all' )

! Make gifs?
   gif = get_parameter('-gif', default_gif )
   if ( gif ) then
      plot = get_parameter('-gif', default_plot )
      plot_solution = ( plot == 'solution' .or. plot == 'all' )
   end if

! Make jpeg?
   jpeg = get_parameter('-jpeg', default_jpeg )
   if ( jpeg ) then
      plot = get_parameter('-jpeg', default_plot )
      plot_solution = ( plot == 'solution' .or. plot == 'all' )
   end if

! Which problem?
   grid = get_parameter('-grid', default_grid )
   load = get_parameter('-load', default_load )

   select case (load)
      case(1)
         nelx = grid*60
         nely = grid*20
         nrhs = 1
         volfrac = 0.5_rp
         penal   = 3.0_rp
         rmin    = 1.5_rp         
      case(2)
         nelx = grid*32
         nely = grid*20
         nrhs = 1
         volfrac = 0.4_rp
         penal   = 3.0_rp
         rmin    = 1.2_rp         
      case(3)
         nelx = grid*30
         nely = grid*30
         nrhs = 1
         volfrac = 0.4_rp
         penal   = 3.0_rp
         rmin    = 1.2_rp         
      case(4)
         nelx = grid*30
         nely = grid*30
         nrhs = 2
         volfrac = 0.4_rp
         penal   = 3.0_rp
         rmin    = 1.2_rp         
   end select

! Allocations and initializations:
   N = 2*(nelx+1)*(nely+1)
   allocate(rho(nely,nelx), rho_old(nely,nelx), dc(nely,nelx) )
   rho = volfrac
   allocate(x(N,nrhs),b(N,nrhs),fixed(N),D(N))

   select case (load)
      case(1)
! Force vector
         b = 0.
         b(2,1) = -1._rp

! Fixed degrees of freedom
        fixed = .false.
        fixed(1:2*(nely+1):2) = .true.
        fixed(2*(nelx+1)*(nely+1)) = .true.
     case (2)
! Force vector
         b = 0.
         b(2*(nelx+1)*(nely+1),1) = -1._rp

! Fixed degrees of freedom
        fixed = .false.
        fixed(1:2*(nely+1)) = .true.
     case (3)
! Force vector
         b = 0.
         b(2*(nelx+1)*(nely+1),1) = -1._rp
         b(2*(nelx)*(nely+1)+2,1) =  1._rp

! Fixed degrees of freedom
        fixed = .false.
        fixed(1:2*(nely+1)) = .true.
     case (4)
! Force vector
         b = 0.
         b(2*(nelx+1)*(nely+1),1) = -1._rp
         b(2*(nelx)*(nely+1)+2,2) =  1._rp

! Fixed degrees of freedom
        fixed = .false.
        fixed(1:2*(nely+1)) = .true.
   end select

! Estimate eigenvalue bounds
   eival = qr_eig( Ke )/k(1)
   eimin = minval( real(eival) )
   eimax = maxval( real(eival) )
   rcenter = (eimax+eimin)/2._rp
   rfoci(2) = eimax
   rfoci(1) = eimin+0.01*eimax

! Store data in user matrix:
   A%nelx  = nelx
   A%nely  = nely
   A%rho   = rho
   A%Ke    = Ke
   A%penal = penal
   A%fixed = fixed

   loop = 0 
   change = 1._rp

! START ITERATION
   do while ( change > 0.01 .and. loop < max_loop )
      loop = loop + 1

! Determine displacements:
      D   = real_diagonal( A )
      P%D = D
      x = pp_idrs( A, b, P = P, center = rcenter, foci = rfoci )

! OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
      c = 0.
      do ely = 1,nely
         do elx = 1,nelx
            n1 = (nely+1)*(elx-1)+ely
            n2 = (nely+1)* elx   +ely
! Get the element vector:
            ind = [2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n2+1, 2*n2+2, 2*n1+1, 2*n1+2]
            dc(ely,elx) = 0._rp
            do irhs = 1, nrhs
               xe = x(ind,irhs)
               ye = matmul(Ke,xe)
               c = c + (rho(ely,elx)**penal)*dot_product(xe,ye)
               dc(ely,elx) = dc(ely,elx)-penal*(rho(ely,elx)**(penal-1.))*dot_product(xe,ye)
            end do
         end do
      end do
 
! FILTERING OF SENSITIVITIES
      dc   = check(nelx,nely,rmin,rho,dc)

      rho_old = rho
      rho  = OC(nelx,nely,rho_old,volfrac,dc) 
 
! PRINT RESULTS
      change = maxval(abs(rho-rho_old))
      write(*,'(a,i4,a,f10.4,a,f6.3,a,f6.3)') ' It.:', loop, ' Obj.: ',c, &
                                             ' Vol.: ', sum(rho)/(nelx*nely), ' ch.: ', change  
! STORE rho IN MATRIX STRUCTURE FOR NEW ITERATION
      A%rho = rho

   end do

! Plot the result:
   if ( plot_solution ) then
      open(30,file='solution.dat')
      rho_max = maxval(rho)
      do j = 1, nely
         do i = 1, nelx
            write(30,'(i3,1x,i3,1x,f9.2)') i, nely-j+1, rho(j,i)/rho_max
         end do
         write(30,*)
      end do
      close(30)
      open(40,file='solution.plt')
      if ( gif ) then
         write(40,*) 'set terminal gif'
         write(40,*) 'set output "solution.gif"'
      elseif ( jpeg ) then
         write(40,*) 'set terminal jpeg'
         write(40,*) 'set output "solution.jpeg"'
      end if
      write(40,'(a)') 'unset key'
      write(40,'(a)') 'unset border'
      write(40,'(a)') 'unset colorbox'
      write(40,'(a)') 'unset tics'
      write(40,'(a,f5.2)') 'set size ratio ', real(nely)/real(nelx)
      write(40,'(a)') 'set palet gray negative'
      select case (load)
         case(1)
            write(40,'(a)') 'set title "MBB beam"'
         case(2)
            write(40,'(a)') 'set title "Cantilever beam, vertical load"'
         case(3)
            write(40,'(a)') 'set title "Cantilever beam, symmetric load"'
         case(4)
            write(40,'(a)') 'set title "Cantilever beam, dual load"'
      end select
      write(40,'(a)') 'plot "solution.dat" with image'
      close(40)
      call execute_command_line('gnuplot -p solution.plt')
   end if

contains

   function OC(nelx,nely,rho_old,volfrac,dc) result(rho) 
!%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      use precision_module

      implicit none
      integer             :: nelx, nely
      real(kind=rp)       :: l1, l2, lmid, move, volfrac
      real(kind=rp)       :: rho(nely,nelx), rho_old(nely,nelx), dc(nely,nelx)
      real(kind=rp)       :: t(nely,nelx)

      l1 = 0. 
      l2 = 1.e5_rp
      move = 0.2_rp
      do while ( l2-l1 > 1e-4 )
         lmid = 0.5_rp*(l2+l1)
         where ( rho_old+move < rho_old*sqrt(-dc/lmid) )
            t = rho_old+move
         elsewhere
            t = rho_old*sqrt(-dc/lmid)
         end where
         where ( 1._rp < t )
            t = 1._rp
         end where
         where ( rho_old-move > t )
            t = rho_old-move
         end where
         where ( 0.001 > t )
            t = 0.001
         end where
         rho = t

         if ( ( sum(rho) - volfrac*nelx*nely ) > 0. ) then
            l1 = lmid
         else
            l2 = lmid
         end if
      end do
   end function OC

   function check(nelx,nely,rmin,rho,dc) result(dcn)
!%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      use precision_module

      implicit none
      integer, intent(in)        :: nelx, nely
      real(kind=rp), intent(in)  :: dc(nely,nelx), rmin, rho(nely,nelx)
      real(kind=rp)              :: dcn(nely,nelx)
      real(kind=rp)              :: som, fac
      integer                    :: i, j, k, l

      dcn = 0.
      do i = 1,nelx
         do j = 1,nely
            som=0.0 
            do k = max(i-floor(rmin),1),min(i+floor(rmin),nelx)
               do l = max(j-floor(rmin),1),min(j+floor(rmin),nely)
                  fac = rmin-sqrt(real((i-k)**2+(j-l)**2,kind=rp))
                  som = som+max(0._rp,fac)
                  dcn(j,i) = dcn(j,i) + max(0._rp,fac)*rho(l,k)*dc(l,k)
               end do
            end do
            dcn(j,i) = dcn(j,i)/(rho(j,i)*som);
         end do
      end do
   end function check

end program top_opt

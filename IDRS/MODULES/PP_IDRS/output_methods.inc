!
! Include file to print the solver parameters
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2025 Martin van Gijzen
!

if ( my_proc == 1 ) then
   write(*,'(a,i5)') 'Number of processors   = ', n_procs
   write(*,'(a,i10,/)') 'Problem size           = ', n_procs*size(b,1)
   if ( use_idrs ) then
      write(method,'(a,i2,a)') 'IDR(',s,')'
      if ( recycle ) then
         write(*,'(a,a,a)') 'Solution method is ',trim(method), ' with recycling.'
      else
         write(*,'(a,a,a)') 'Solution method is ',trim(method),'.'
      end if
   elseif ( use_bicgstab ) then
      method = 'BiCGSTAB'
      if ( recycle ) then
         write(*,'(a,a,a)') 'Solution method is ',trim(method), ' with recycling.'
      else
         write(*,'(a,a,a)') 'Solution method is ',trim(method),'.'
      end if
   elseif ( use_qmridr ) then
      write(method,'(a,i2,a)') 'QMRIDR(',s,')'
      if ( in_it > 0 ) then
         write(*,'(a,a,a,i4,a)') &
            'Solution method is ', trim(method), ' with maximum ', in_it, ' inner iterations.'
      else
         write(*,'(a,a,a)') 'Solution method is ', trim(method),'.'
      end if
   elseif ( use_minsync ) then
      write(method,'(a,i2,a)') 'IDR(',s,')_minsync'
      if ( recycle ) then
         write(*,'(a,a,a)') 'Solution method is ',trim(method), ' with recycling.'
      else
         write(*,'(a,a,a)') 'Solution method is ',trim(method),'.'
      end if
   end if

   if ( use_neumann ) then
      if ( present(center) ) then
         write(pol_prec,'(a,i2,a)') &
            'Neumann Preconditioner of degree ', degree, ' with pre-computed parameter.'
         write(*,'(a)') pol_prec
      else
         write(pol_prec,'(a,i2,a)') &
            'Neumann Preconditioner of degree ', degree, ' with adaptively computed parameter.'
         write(*,'(a)') pol_prec
      end if
   elseif ( use_chebyshev ) then
      if ( present(foci) ) then
         write(pol_prec,'(a,i2,a)') &
            'Chebyshev Preconditioner of degree ', degree, ' with pre-computed parameters.'
         write(*,'(a)') pol_prec
      else
         write(pol_prec,'(a,i2,a)') &
            'Chebyshev Preconditioner of degree ', degree, ' with adaptively computed parameter.'
         write(*,'(a)') pol_prec
      end if
   elseif ( use_hessenberg ) then
      write(pol_prec,'(a,i2,a)') 'Hessenberg Preconditioner of degree ', degree, '.'
      write(*,'(a)') pol_prec
   else
      write(pol_prec,'(a)') 'No polynomial preconditioning.'
      write(*,'(a)') pol_prec
   end if

   if ( present(V) ) write(*,'(a,i2)') 'Deflation space of dimension ', size(V,2)

   write(*,*)

end if

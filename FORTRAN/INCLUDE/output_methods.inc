!
! Include file to print the solver parameters
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2024 Martin van Gijzen
!

if ( this_image() == 1 ) then
   write(*,'(a,i5)') 'Number of processors   = ', num_images()
   if ( use_idrs ) then
      if ( recycle ) then
         write(*,'(a,i2,a)') 'Solution method is IDR(s) with s = ', s, ' and with recycling. '
      else
         write(*,'(a,i2)') 'Solution method is IDR(s) with s = ', s
      end if
   end if
   if ( use_bicgstab ) write(*,'(a)') 'Solution method is BiCGSTAB '
   if ( use_qmridr ) then
      if ( in_it > 0 ) then
         write(*,'(a,i2,a,i4,a)') &
            'Solution method is QMRIDR(s) with s = ', s, ' and with maximum ', in_it, ' inner iterations. '
      else
         write(*,'(a,i2)') 'Solution method is QMRIDR(s) with s = ', s
      end if
   end if
   if ( use_neumann ) write(*,'(a,i2)') 'Preconditioner is Neumann of degree = ', degree
   if ( use_chebyshev ) write(*,'(a,i2)') 'Preconditioner is Chebyshev of degree = ', degree
   write(*,*)
end if

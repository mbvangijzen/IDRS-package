!
! Set the working precision. 
! Note that the parameters rp and cp are global to all routines and modules that use this module.
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2024 Martin van Gijzen
! 

module precision_module

! Define the precision for real and complex arithmatic
   use iso_fortran_env
   integer, parameter :: rp = real64                   ! kind for real
   integer, parameter :: cp = real64                   ! kind for complex

   public

end module precision_module

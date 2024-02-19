module interface_module

!
! Module to read the command line parameters
!
! The command line has to be read once at the beginning of the programme:
!    call get_command_line
! Parameters can then be read using a keyword parameter combination,
! You always have to specify a default value for the parameter of
! of the correct type.
!
! Examples:
!     grid      = get_parameter('-grid',3)
!     bicgstab  = get_parameter('-bicgstab', .false.)
!     tol       = get_parameter('-tol', 1.e-6)
!
! This software is distributed under the MIT License:
! http://www.opensource.org/licenses/mit-license.php
! Copyright:(c) 2024 Martin van Gijzen
!
!

   implicit none

! Global variables for command line parameters
   character(len=:), allocatable :: parameter_list
   integer                       :: parameter_length

   interface get_parameter
      module procedure get_real_parameter, get_integer_parameter, get_logical_parameter, get_character_parameter
   end interface

contains

subroutine initialize()

   implicit none

! Read the command line
   call get_command(length=parameter_length)
   allocate(character(parameter_length) :: parameter_list)
   call get_command(parameter_list)

end subroutine initialize

function get_real_parameter( parameter_name, default_value ) result(parameter_value)

   real                    :: parameter_value
   real                    :: default_value
   character(len = *)      :: parameter_name
   integer                 :: k, kb, ke, err_stat

! Set default_value:
   parameter_value = default_value

! Read the parameter (if it is present)
   k = index( parameter_list, parameter_name )
   if ( k > 0 ) then
! Read the value of parameter_name
      ke = k -1+scan( parameter_list(k:parameter_length)," " )
      kb = ke-1+verify( parameter_list(ke:parameter_length)," " )
      k  = scan( parameter_list(kb:parameter_length)," " )
      if ( k == 0 ) then
         ke = parameter_length
      else
         ke = kb-1+k
      end if
      read(parameter_list(kb:ke),*,iostat=err_stat) parameter_value
      if ( err_stat .ne. 0 ) then
          write(*,'(a,a,a)') 'Parameter for keyword ', parameter_name, ' missing.'
          error stop
      end if
   end if

end function get_real_parameter

function get_integer_parameter( parameter_name, default_value ) result(parameter_value)
   integer                 :: parameter_value
   integer                 :: default_value 
   character(len = *)      :: parameter_name
   integer                 :: k, kb, ke, err_stat

! Set default_value:
   parameter_value = default_value

! Read the parameter (if it is present)
   k = index( parameter_list, parameter_name )
   if ( k > 0 ) then
! Read the value of parameter_name
      ke = k -1+scan( parameter_list(k:parameter_length)," " )
      kb = ke-1+verify( parameter_list(ke:parameter_length)," " )
      k  = scan( parameter_list(kb:parameter_length)," " )
      if ( k == 0 ) then
         ke = parameter_length
      else
         ke = kb-1+k
      end if
      read(parameter_list(kb:ke),*,iostat=err_stat) parameter_value
      if ( err_stat .ne. 0 ) then
          write(*,'(a,a,a)') 'Parameter for keyword ', parameter_name, ' missing.'
          error stop
      end if
   end if

end function get_integer_parameter

function get_logical_parameter( parameter_name, default_value ) result(parameter_value)
   logical                 :: parameter_value
   logical                 :: default_value 
   character(len = *)      :: parameter_name
   integer                 :: k, kb, ke

! Set default_value:
   parameter_value = default_value

! Read the parameter (if it is present)
   k = index( parameter_list, parameter_name )
   if ( k > 0 ) parameter_value = .true.

end function get_logical_parameter

function get_character_parameter( parameter_name, default_value ) result(parameter_value)
   character(len = parameter_length) :: parameter_value
   character(len = *)                :: parameter_name, default_value
   integer                           :: k, kb, ke

! Set default_value:
   parameter_value = default_value

! Read the parameter (if it is present)
   k = index( parameter_list, parameter_name )

   if ( k > 0 ) then
! Read the value of parameter_name
      ke = k -1+scan( parameter_list(k:parameter_length)," " )
      kb = ke-1+verify( parameter_list(ke:parameter_length)," " )
      k  = scan( parameter_list(kb:parameter_length)," " )
      if ( k == 0 ) then
         ke = parameter_length
      else
         ke = kb-1+k
      end if
      parameter_value = parameter_list(kb:ke)
   end if

end function get_character_parameter

end module interface_module

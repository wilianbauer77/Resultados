!> Module containing modified routines
module leopardv3_mod
implicit none
integer :: maxit= 100 ! Save # of Muller steps as the original version

CONTAINS
   include 'muller.f90'
end module leopardv3_mod

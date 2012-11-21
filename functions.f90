!===============================================================================
!  functions module contains useful functions (not subroutines) for speeding
!  up the writing of calculations.  These are not geometrical, purely algebraic.
!
!  Andy Nowacki, University of Bristol
!  2011
!===============================================================================
module functions

   public kron_del
   
   contains
   
!===============================================================================
   function kron_del(i,j)
!-------------------------------------------------------------------------------
!  kron_del returns the Kronecker delta of two integers: 0 or 1.

	  integer,intent(in) :: i,j
	  integer :: kron_del
	  
	  if (i == j) kron_del = 1
	  if (i /= j) kron_del = 0
	  
	  return
	  
   end function kron_del
!-------------------------------------------------------------------------------

end module functions
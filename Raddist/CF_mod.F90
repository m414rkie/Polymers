module functions

contains

real function dist(x,y,z)
! Function checks the distance of two molecules
! in 3D space using simple distance formula
! xyz mol 1; abc mol2

implicit none
		real,intent(in) :: x, y, z

	dist = sqrt(x*x + y*y + z*z)

end function dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer function chain(beadNum)
! Function finds which chain a molecule
! belongs to. Each chain is assumed to be
! 40 molecules in length

implicit none
	real,intent(in) :: beadNum

	chain = ceiling(beadnum/40.0)

end function chain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module chain_Functions

contains

integer function chainEnds(beadNum)
! Function finds which half of a chain
! a molecule belongs to. Each chain is
! assumed to 40 molecules in length.

use functions, only : chain

implicit none
	integer,intent(in) :: beadNum
	integer	    	   :: chainNum

	chainNum = chain(real(beadNum))

	if (mod(float(beadNum),40.0) .gt. 20) then
		chainEnds = chainNum + 1
	else
		chainEnds = chainNum
	end if

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module

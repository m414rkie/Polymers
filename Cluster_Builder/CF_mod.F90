module chain_Functions
! Module contains a function that determines which chain a bead is on (chain),
! also a function that determine which half of a chain a bead is on (chainEnds)

contains

integer function chainEnds(beadNum)
! Function finds which half of a chain
! a bead belongs to. Each chain is
! assumed to 40 molecules in length.

implicit none
	integer,intent(in) :: beadNum
	integer	    	   :: chainNum

	chainNum = chain(real(beadNum))

	if (mod(float(beadNum),40.0) .gt. 20) then
		chainEnds = 2*chainNum
	else
		chainEnds = 2*chainNum - 1
	end if

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer function chain(beadNum)
! Function finds which chain a molecule
! belongs to. Each chain is assumed to be
! 40 molecules in length

implicit none
	real,intent(in) :: beadNum

	chain = ceiling(beadnum/40.0)

end function chain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module

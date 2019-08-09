module functions

contains

real function dist(x,y,z,a,b,c)
! Function checks the distance of two molecules
! in 3D space using simple distance formula
! xyz mol 1; abc mol2

implicit none
		real,intent(in) :: x, y, z, a, b, c

	dist = sqrt(((a-x)**2) + ((b-y)**2) + ((c-z)**2))

end function dist

end module

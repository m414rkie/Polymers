module functions

contains

real function dist(x,y,z)
! Function checks the distance of two molecules
! in 3D space using simple distance formula
! xyz mol 1; abc mol2

implicit none
		real,intent(in) :: x, y, z

	dist = sqrt((x**2) + (y**2) + (z**2))

end function dist

end module

! Program takes an output file from the LAMMPS software package and
! reads in the location and type of each molecule. From this, the radial
! distance function is determined from xyz data given in the input file.

! periodic boundaries rehandled

! Author: Jon Parsons
! Date: 2-1-19

program clusfinder

use functions

implicit none
	character*50			:: filename  ! Name of input file
	real				  		:: tstep ! Distance parameter, number of clusters found, cuurent timestep
	integer						:: numMols ! Number of molecules in system, ! Number of chains per molecule
	real			  			:: boxDim(3,2) ! Min and Max values of x,y,z respectively
	integer						:: numTsteps ! Number of time steps looked at
	integer						:: ioErr, j  ! System error variable, looping integer
	integer						:: r_num ! The discretized distances that a bead can be.
	real							:: r_max, dr ! Maximum distance bsad on box size
	real,allocatable	:: molData(:,:)    ! molnumber, moltype, x, y, z, cluster, molgroup
	real,allocatable	:: dist_arr(:) ! Holds the values for the radial distance functio.
	real						 	:: vol, xd, yd, zd ! Values for the box. Volume and largest axial distance
	character*1				:: type ! For user input, selects which bead type is being examined

! User input.
write(*,*) "Please choose which bond type to examine."
write(*,*) "Hydrogen (H), or Sulfide (S)"
read(*,*) type

! Parse
call chartoup(type,type)

100 write(*,*) "Please enter the name of the file with the data."
write(*,*) "Typical files will begin with the 'pict' prefix."
write(*,*) "If the file is not in this directory enter the full path."
read(*,*) filename

filename = trim(filename)

open(unit=15, file=filename, status="old", action="read", iostat=ioErr)

if (ioErr .ne. 0) then
	write(*,*) "File not found, please try again."
	goto 100
end if

! Initialize
numTsteps = 0

! Read in header data
read(15,*)
read(15,*) tstep
read(15,*)
read(15,*) numMols
read(15,*)
read(15,*) boxDim(1,1), boxDim(1,2)
read(15,*) boxDim(2,1), boxDim(2,2)
read(15,*) boxDim(3,1), boxDim(3,2)
read(15,*)

! Determine box characteristics
xd = boxDim(1,2) - boxDim(1,1)
yd = boxDim(2,2) - boxDim(2,1)
zd = boxDim(3,2) - boxDim(3,1)
vol = xd*yd*zd

r_max = xd*0.5

! Determine number of boxes for rad. distribution array. Assumes a cube
dr = 0.5
r_num = ceiling(r_max/dr) + 1
write(*,*) "Bin Size: ", dr, r_num, r_max

! Allocation
allocate(dist_arr(r_num), stat= ioErr)

if (ioErr .ne. 0) then
		write(*,*) "Failed to allocate distance array. Exiting"
		stop
end if

dist_arr = 0.0

! Reads in the data and calls the clustering subroutine until EOF
do
  ! Iterate timestep
	numTsteps = numTsteps + 1

	! Allocate data array
	allocate(molData(numMols,5), stat = ioErr)

	if (ioErr .ne. 0) then
		write(*,*) "Failed to allocated primary array. Exiting at timestep", tstep
		stop
	end if

	write(*,*) "Beginning time:", tstep

	! Initial values, to be over-ridden
	molData = 0.0

	! Read in molecule data
	fileread: do j = 1, numMols, 1

		read(15,*) molData(j,1), molData(j,2), molData(j,3), molData(j,4), molData(j,5)
							! bead #     , type        , x           , y           , z
	end do fileread

	! Call distribution subroutine
 	call rad_dist(moldata,numMols,5,dist_arr,r_num,vol,xd*0.5,yd*0.5,zd*0.5,type,dr)

	deallocate(molData)

	! Checks for EOF, if not then reads and discards header data for next step
	read(15,*,END=101)
	read(15,*) tstep
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)

end do

101 write(*,*) "Number of timesteps checked:", numTsteps
write(*,*) "End of input file reached. Goodbye"

! Final output
call dist_print(dist_arr,r_num,dr,numTsteps,type)
end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine rad_dist(arrin,dim1,dim2,arrout,r_num,vol,xbx,ybx,zbx,type,dr)
! Subroutine to determine the radial distribution function at each step.

use functions
use chain_Functions

implicit none
	integer,intent(in)		:: dim1, dim2  ! Dimensions of input array
	integer,intent(in)		:: r_num ! number of possible bins
	character*1,intent(in):: type !  bead type
	real,intent(in)				:: xbx, ybx, zbx ! Maximum axial distance beads can be from each other
	real,intent(in)				:: dr ! size of the bins
	real,intent(in)				:: vol ! total volume of box, time step
	real,intent(in)				:: arrin(dim1,dim2) ! Input array, holds bead location and type
																				! molnumber, moltype, x, y, z, cluster, molgroup
	real,intent(inout)		:: arrout(r_num) ! Holds the values for the radial distribution function
	real									:: arrout_temp(r_num)  ! holds values for current timestep
	integer								:: i, d, j, k  ! Looping integers
	integer								:: count_tot ! Holds number of beads in total
	real									:: xi, xj, yi, yj, zi, zj ! coordinates of beads 1, 2
	real									:: xd, yd, zd ! Axial distance of each pair of beads
	real									:: r, distance ! Distance being considered, distance from i-th particle
	real									:: shell ! Will hold volume of shell
	real,parameter				:: pi = acos(-1.0) ! Pi
	character*15					:: filename ! Holds filename. Changes based on which bead we are examining
	integer								:: b_type, b_type2 ! numbers associated with the bead type

if (type .eq. 'H') then
	b_type = 5
	b_type2 = 0
	filename = 'timrad_h.dat'
else if (type .eq. 'S') then
	filename = 'timrad_s.dat'
	b_type = 3
	b_type2 = 1
end if

write(*,*) xbx, r_num

! Initialize temporary array
arrout_temp = 0.0

! Initialize counting variables
count_tot = 0

! Bead being considered
outer_loop : do i = 1, dim1, 1

		! skip if not right bead type
		if (nint(arrin(i,2)) .ne. b_type) then
			if (nint(arrin(i,2)) .ne. b_type2) then
				cycle outer_loop
			end if
		end if

		! Add bead to total number of beads
		count_tot = count_tot + 1

		! Bead 1's (x,y,z) coord's
		xi = arrin(i,3)
		yi = arrin(i,4)
		zi = arrin(i,5)

			! Loop for second bead
			inner_loop : do j = i, dim1, 1

					! skip if not right bead type
					if (nint(arrin(j,2)) .ne. b_type) then
						if (nint(arrin(j,2)) .ne. b_type2) then
							cycle inner_loop
						end if
					end if

					if (j .eq. i) then
						cycle inner_loop
					end if

					! Bead 2's (x,y,z) coord's
					xj = arrin(j,3)
					yj = arrin(j,4)
					zj = arrin(j,5)
					! axial distances
					xd = abs(xj - xi)
					yd = abs(yj - yi)
					zd = abs(zj - zi)

					!! Periodic Boundary check, if it is closer to go through the
					! boundary wall this section does so.
					! X dimension check
					if (xd .ge. xbx) then
						xd = xd - xbx
					end if
					! Y dimension check
					if (yd .ge. ybx) then
						yd = yd - ybx
					end if
					! Z dimension check
					if (zd .ge. zbx) then
						zd = zd - zbx
					end if

					! Determine distance
					distance = dist(xd,yd,zd)

					! Place the beads in the appropriate bins
					k = ceiling(distance/dr)
					if (k .le. r_num*0.75) arrout_temp(k) = arrout_temp(k) + 1

			end do inner_loop

end do outer_loop

write(*,*) "Number of pairs counted: ", sum(arrout_temp)

! Divide each box by the volume of the shell it represents
do d = 1, r_num, 1
	r = dr*float(d)
	shell = 1.3333*pi*(((r+dr)**3) - (r**3))
	arrout_temp(d) = arrout_temp(d)/(shell*count_tot)
end do

! Normalize, divide by overall density of box (concerning only the beads we care about)
arrout_temp = arrout_temp*vol/count_tot

! User output
write(*,*) "total beads of desired type:", count_tot

! Add current array to total array. After all timesteps are checked this will be
! time-averaged and output
arrout = arrout + arrout_temp

! Print current distribution function to file
open(unit=20,file=trim(filename),status="unknown",position="append")

do d = 1, r_num, 1
	r = dr*float(d)
	write(20,*) r, arrout_temp(d)
end do

close(20)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dist_print(arrin,r_num,dr,t,btype)
! Subroutine to print time-averaged radial distribution function.

implicit none
	integer				:: r_num, t ! Size of array, number of timesteps
	real					:: arrin(r_num) ! Array to be time-averaged
	real					:: dr !  distance interval
	character*1		:: btype ! Holds bead of interest, for use in filename

	real					:: r ! Current distance
	integer				:: i ! Looping integer
	character*20	:: filename

	if (btype .eq. 'H') then
		filename = "rad_dist_h.dat"
	else
		filename = "rad_dist_s.dat"
	end if

! Time average, the array contains the sum of all previous distributions
arrin = arrin/float(t)

! Open file for printing
open(unit=18,file=trim(filename),status="replace",position="append")

! Printing loop
do i = 1, r_num, 1
	r = float(i)*dr

	write(18,*) r, arrin(i)

end do

close(18)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine chartoup(stringin,stringout)
! converts text input to upper case

implicit none
	character(*)					   :: stringin ! string to adjust
	character(len(stringin)) :: stringout ! capitalized string
	integer									 :: i, j ! looping integer, iachar indice value

do i = 1, len(stringin), 1 ! loop through string
	j = iachar(stringin(i:i)) ! get iachar indice
		! replace with uppercase version if needed
		if(j .ge. iachar("a") .and. j .le. iachar("z")) then
			stringout(i:i) = achar(iachar(stringin(i:i))-32)
		else
			stringout(i:i) = stringin(i:i)
		end if
end do

end subroutine

! Program takes the output from the Diffus.F90 code and computes the storage
! and loss moduli from the MSD data contained in diff_out.dat.

! Version 0.0
! - initial implementations:
!				. Data input

! Author: Jon Parsons
! Date 11-6-19

program moduli

implicit none
	character*50		 :: raw_in ! Input file, user entered

	real,allocatable :: mast_arr(:,:) ! Holds input data
	real,allocatable :: mod_arr(:,:) ! Holds moduli (dim2) at freq (dim1)
	real						 :: te ! temperature, user input
	real						 :: junk ! For discarding some unneeded data
	integer					 :: points ! Number of data points
	real						 :: prev, now

	integer					 :: i ! Looping integer
	integer					 :: ioErr ! Error handling variable

! Get file name
write(*,*) "Please enter the name of the input file:"
read(*,*) raw_in

write(*,*) "Please enter the temperature of the system"
read(*,*) te

! Determine number of data points
open(unit=15,file=trim(raw_in),status="old",action="read")

points = 0
read_loop: do

	read(15,*,END=101)
	points = points + 1
end do read_loop

101 close(15)

! Allocate arrays
allocate(mast_arr(points,2), stat=ioErr)

if (ioErr .ne. 0) then
	write(*,*) "Allocation of input array failed. Exiting"
	stop
end if
mast_arr = 0.0

allocate(mod_arr(points,3), stat=ioErr)

if (ioErr .ne. 0) then
	write(*,*) "Allocation of moduli array failed. Exiting"
	stop
end if
mod_arr = 0.0

! Collect data
prev = 0.0
open(unit=15,file=trim(raw_in),status="old",action="read")
data_in: do i = 1, points, 1
	read(15,*) mast_arr(i,1), now, junk
	mast_arr(i,2) = now - prev
	prev = now
end do data_in

call modulis(mast_arr,mod_arr,points,2,te)
call outputs(mod_arr,points,3)

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine modulis(ins,outs,dim1,dim2,temp)
	implicit none
		integer,intent(in)		:: dim1, dim2
		real,intent(in)				:: temp
		real,intent(out)			:: outs(dim1,dim2+1)
		real,intent(inout)		:: ins(dim1,dim2)

		real									:: k_b = 1.38064852! Boltzmann Constant
		real									:: pi = 3.1415926

		integer								:: i
		real									:: alpha, Gstar, G_prime, G2_prime
		real									:: gamm, time, radius

		real									:: f, x, y, z, w
		real									:: freq

! Logorithmic derivative
f(x,y,z,w) = (log(x-y))/((log(z-w)))

! Radius of particle
radius = 0.5

do i = 2, dim1-1, 1
	time = ins(i,1) ! Time interval
	freq = 2.0*pi/time ! Frequency

	! Factor alpha
	alpha = f(ins(i+1,2),ins(i-1,2),ins(i+1,1),ins(i,1))
	gamm = gamma(1.0+alpha)

	write(*,*) alpha, gamm
	! Get complex moduli
	Gstar = k_b*temp/(pi*radius*ins(i,2)*gamm)
	! Storage and loss factors
	G_prime = abs(Gstar)*cos(pi*alpha/2.0)
	G2_prime = abs(Gstar)*sin(pi*alpha/2.0)

	! Store for output
	outs(i,1) = freq
	outs(i,2) = G_prime
	outs(i,3) = G2_prime

end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine outputs(outs,dim1,dim2)
	implicit none
		integer,intent(in)		:: dim1, dim2
		real,intent(in)				:: outs(dim1,dim2)

		integer								:: i

open(unit=15,file="modulis.dat",status="replace",position="append")
write(15,*)"w    G'    G''"
do i = 1, dim1-1, 1
	write(15,*) outs(i,1), outs(i,2), outs(i,3)
end do

close(15)

end subroutine

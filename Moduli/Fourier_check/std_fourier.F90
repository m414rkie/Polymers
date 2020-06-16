! Program takes the output from the Diffus.F90 code and computes the storage
! and loss moduli from the MSD data contained in diff_out.dat.

! Version 3.0
! standard discrete fourier transform

! Author: Jon Parsons
! Date 6-6-2020

program moduli

implicit none
	character*50		 :: raw_in ! Input file, user entered

	real,allocatable :: mast_arr(:,:) ! Holds input data
	real,allocatable :: mod_arr(:,:) ! Holds moduli (dim2) at freq (dim1)
	real						 :: junk ! For discarding some unneeded data
	integer					 :: points ! Number of data points

	integer					 :: i ! Looping integer
	integer					 :: ioErr ! Error handling variable

! Get file name
write(*,*) "Please enter the name of the input file:"
read(*,*) raw_in

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
open(unit=15,file=trim(raw_in),status="old",action="read")
data_in: do i = 1, points, 1
	read(15,*) mast_arr(i,1), mast_arr(i,2), junk
end do data_in

call modulis(mast_arr,mod_arr,points,2)
call outputs(mod_arr,points,3)

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine modulis(ins,outs,dim1,dim2)
	implicit none
		integer,intent(in)	:: dim1, dim2
		real,intent(out)		:: outs(dim1-1,3)
		real,intent(in)			:: ins(dim1,dim2)

		real,parameter			:: pi = 3.14159 ! pi
		real								:: freq_i ! Current Frequency
		real								:: dt, df ! sample, frequency spacing
		real								:: arg ! weighted angle
		complex							:: ex_arg ! complex weighted angle

		integer							:: i, j
		complex							:: loc_sum ! Placeholder values
		complex							:: g_interim(dim1) ! holds values
		complex							:: freq_mod ! Holds complex moduli at a frequency

outs = 0.0
dt = ins(2,1) - ins(1,1)
df = (2.0/(dt))/(float(dim1)) ! scaling
write(*,*) "most freq:", 1.0/(2.0*dt)
! Loop over frequencies
freq_loop: do i = 1, dim1-1, 1

	! Set current Frequency
	freq_i = float(i)
	! Begin Fourier Portion
	loc_sum = (0.0,0.0)

	sum_loop: do j = 1, dim1-1, 1
		arg = -2.0*pi*float(j)*freq_i/float(dim1)
		ex_arg = cmplx(0.0,arg)
		loc_sum = loc_sum + cmplx(ins(j,2))*cexp(ex_arg)
	end do sum_loop

	g_interim(i) = 2.0*loc_sum/float(dim1)

	! adjust freq and save to out array
	outs(i,1) = df*float(i)!freq_i

end do freq_loop

! Assign individual values
sep_loop: do i = 1, dim1-1, 1
	outs(i,2) = real(g_interim(i))
	outs(i,3) = aimag(g_interim(i))
end do sep_loop

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine outputs(outs,dim1,dim2)
	implicit none
		integer,intent(in)		:: dim1, dim2
		real,intent(in)				:: outs(dim1,dim2)

		integer								:: i

open(unit=15,file="modulis.dat",status="replace",position="append")
do i = 1, dim1-1, 1
	write(15,*) outs(i,1), outs(i,2), outs(i,3)
end do

close(15)

end subroutine

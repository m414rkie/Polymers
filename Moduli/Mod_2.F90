! Program takes the output from the Diffus.F90 code and computes the storage
! and loss moduli from the MSD data contained in diff_out.dat.

! Version 2.0
! - initial implementations

! Author: Jon Parsons
! Date 11-6-19

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
		integer,intent(in)		:: dim1, dim2
		real,intent(out)			:: outs(dim1-1,3)
		real,intent(inout)		:: ins(dim1,dim2)

		complex								:: rt1 = (0.0,1.0) ! i
		real									:: om_i, om_f ! Initial and final frequencies of interest
		real									:: freq(dim1) ! Array holding frequencies
		real									:: freq_i ! Current Frequency
		real									:: d_freq ! Change in frequency value, used in fill
		real									:: slp_infty ! Slope at infinity

		integer								:: i, j
		complex								:: i_om_t_f, i_om_t_i ! Complex values
		complex								:: sum_in, loc_sum, loc_slp ! Placeholder values
		complex								:: g_interim(dim1)
		complex								:: freq_mod ! Holds complex moduli at a frequency

		real									:: f, x, y, z, w

! Local derivative function
f(x,y,z,w) = (x-y)/(z-w)

! Generate frequencies
om_i = 1E-4
om_f = 1E5
d_freq = (om_f - 1E-6)/(dim1-1)

freq_fill: do i = 1, dim1-1, 1
	freq(i) = om_i + float(i)*d_freq
end do freq_fill

! find slope at infinity
slp_infty = f(ins(dim1,2),ins(dim1-5,2),ins(dim1,1),ins(dim1-5,1))

! Loop over frequencies
freq_loop: do i = 1, dim1-1, 1

	! Initialize summation
	sum_in = (0.0,0.0)

	! Set current Frequency
	freq_i = freq(i)
	! Begin Fourier Portion
	sum_loop: do j = 2, dim1, 1
		loc_sum = (0.0,0.0) ! Internal summation of Fourier
		loc_slp = f(ins(j,2),ins(j-1,2),ins(j,1),ins(j-1,1)) ! Local slope
		i_om_t_i = cmplx(0.0,-freq_i*ins(j-1,1)) ! exponential arguement, t-1
		i_om_t_f = cmplx(0.0,-freq_i*ins(j,1)) ! exponential arguement, t
		loc_sum = cmplx(loc_slp)*(cexp(i_om_t_i) - cexp(i_om_t_f))
		sum_in = sum_in + loc_sum ! Add internal to overall summation
	end do sum_loop

	! final computation
	freq_mod = (cmplx(1.0,0.0) - cexp(-rt1*cmplx(freq_i,0.0))*cmplx(ins(1,1),0.0)/cmplx(ins(1,1),0.0)) + &
	 						cmplx(slp_infty,0.0)*cexp(-rt1*cmplx(freq_i*ins(dim1,1),0.0)) + sum_in

	! Isolate moduli
	freq_mod = -freq_mod/(freq_i*freq_i)

	! Scale to frequency
	freq_mod = cmplx(1.0,0.0)/(rt1*cmplx(freq_i,0.0)*freq_mod)

	g_interim(i) = freq_mod

	! save frequency to out array
	outs(i,1) = freq_i

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
write(15,*)"w    G'    G''"
do i = 1, dim1-1, 1
	write(15,*) outs(i,1), outs(i,2), outs(i,3)
end do

close(15)

end subroutine

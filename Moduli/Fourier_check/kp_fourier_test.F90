! Program takes the output from the Diffus.F90 code and computes the storage
! and loss moduli from data created by a fourier testing program

! Author: Jon Parsons
! Date: 6-6-2020

program moduli

implicit none
	character*50		 :: raw_in ! Input file, user entered

	real,allocatable :: mast_arr(:,:) ! Holds input data
	real,allocatable :: mod_arr(:,:) ! Holds moduli (dim2) at freq (dim1)
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
	read(15,*) mast_arr(i,1), mast_arr(i,2)
end do data_in

call modulis(mast_arr,mod_arr,points,2)
call outputs(mod_arr,points,3)

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine modulis(ins,outs,dim1,dim2)
	implicit none
		integer,intent(in)	:: dim1, dim2
		real,intent(out)		:: outs(dim1-1,3)
		real,intent(in)			:: ins(dim1,dim2) ! ins(:,1) == time ; ins(:,2) == MSD

		complex							:: rt1 = (0.0,1.0) ! i
		real								:: dt, df ! sample, frequency spacing
		real								:: freq_i ! Current frequency
		real								:: slp_infty ! Slope at infinity

		integer							:: i, j ! Looping integers
		complex							:: i_om_t_f, i_om_t_i ! Complex values
		complex							:: sum_in, loc_sum, loc_slp ! Placeholder values
		complex							:: g_interim(dim1) ! Holding variable for |G|
		complex							:: Diff_const ! Diffusion coefficient
		complex							:: freq_mod ! Holds complex moduli at a frequency

		real								:: f, x, y, z, w ! Function variables

! Local derivative function
f(x,y,z,w) = (x-y)/(z-w)

outs = 0.0
dt = ins(2,1) - ins(1,1)
df = 2.0*(1.0/dt)/float(dim1) ! scaling

! find slope at infinity -> Diffusion Coefficient D = slp/2t_N
slp_infty = f(ins(dim1,2),ins(dim1-(dim1/4),2),ins(dim1,1),ins(dim1-(dim1/4),1))
slp_infty = slp_infty/2.0!(2.0*ins(dim1,1))
Diff_const = cmplx(slp_infty,0.0)

! Loop over frequencies
freq_loop: do i = 1, dim1-1, 1

	! Set current Frequency
	freq_i = float(i)

	! Initialize summation
	sum_in = (0.0,0.0)
	freq_mod = (0.0,0.0)

	! Begin Fourier Portion
	sum_loop: do j = 2, dim1, 1
		loc_slp = f(ins(j,2),ins(j-1,2),ins(j,1),ins(j-1,1)) ! Local slope
		i_om_t_i = cmplx(0.0,-freq_i*ins(j-1,1)) ! exponential arguement, t-1
		i_om_t_f = cmplx(0.0,-freq_i*ins(j,1)) ! exponential arguement, t
		loc_sum = cmplx(loc_slp)*(cexp(i_om_t_i) - cexp(i_om_t_f)) ! Internal sum
		sum_in = sum_in + loc_sum ! Add internal to overall summation
	end do sum_loop

	! final computation
	freq_mod = (cmplx(1.0,0.0) - cexp(cmplx(0.0,-freq_i*ins(1,1))))*cmplx(ins(1,2),0.0)/cmplx(ins(1,1),0.0) + &
							(Diff_const)*cexp(cmplx(0.0,-freq_i*ins(dim1,1))) + sum_in

	!freq_mod = (Diff_const)*cexp(cmplx(0.0,-freq_i*ins(dim1,1))) + sum_in

	! Isolate moduli
	freq_mod = freq_mod/(-freq_i*freq_i)

	! Scale to frequency

	g_interim(i) = freq_mod

	! save frequency to out array
	outs(i,1) = freq_i*df

end do freq_loop

! Assign individual values
sep_loop: do i = 1, dim1-1, 1
	outs(i,2) = (2.0/float(dim1))*abs(g_interim(i))!real(g_interim(i))
	outs(i,3) = 0.0!aimag(g_interim(i))
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

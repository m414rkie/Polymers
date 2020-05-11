! Program takes an output file from the LAMMPS software package and
! reads in the location and type of each molecule. From this, this
! program determines the diffusion of each bead.

! This version does not track chains, only beads

! Version 5
! - Reads indefinite number of files, user will supply the number and names
! 	of the datafiles.
! - Iterates over the stride to output data as avg diffusion - dtime (dtau)
! - Flattens data before passing to diffusion now.
! - Each stride has consistent numbers of points

! Now with periodic boundary handling

! Author: Jon Parsons
! Date: 2-1-19

program diffusion

implicit none
	character*50,allocatable	:: file_arr(:)
	character*50							:: raw_in
	real				  						:: tstep, numMols_r
	integer										:: numMols
	integer										:: numTsteps, t_cur
	integer										:: num_files
	real											:: lj_t
	integer										:: ioErr
	integer										:: i, j
	integer										:: mol_num
	real											:: xmin, xmax, ymin, ymax, zmin, zmax
	real											:: xlen, ylen, zlen
	real											:: x, y, z
	real											:: junk
	real,allocatable					:: mast_arr(:,:,:) ! (molNum,x-y-z,time)
	real											:: time_prev
	integer										:: t_stride

write(*,*) "Please enter the number of LJ timesteps per data snapshot."
read(*,*) lj_t

write(*,*) "Please enter the stride for the timesteps"
read(*,*) t_stride

write(*,*) "Please enter the number of files."
read(*,*) num_files

allocate(file_arr(num_files), stat=ioErr)
if (ioErr .ne. 0) then
	write(*,*) "Allocation of file name array has failed. Exiting"
	stop
end if

file_in: do i = 1, num_files, 1
	write(*,*) "Please enter the name of file", i
	read(*,*) raw_in
	file_arr(i) = trim(raw_in)
end do file_in

! This section to determine the total number of timesteps and number of molecules
numTsteps = 0
t_cur = 0
time_prev = -1.0

file_iter_time: do i = 1, num_files, 1
			open(unit=15,file=file_arr(i),status="old",action="read")
			write(*,*) "Reading file ", file_arr(i)

! Reads in the data and calls the clustering subroutine
read_loop: do

	read(15,*,END=101)
	read(15,*) tstep
	read(15,*)
	read(15,*) numMols_r
	numMols = int(numMols_r)
	read(15,*)
	read(15,*) xmin, xmax
	read(15,*) ymin, ymax
	read(15,*) zmin, zmax
	read(15,*)

	if (tstep .ne. time_prev) then
		numTsteps = numTsteps + 1
		time_prev = tstep
	end if

	! Read in molecule data
	file_read: do j = 1, numMols, 1

		read(15,*)

	end do file_read

end do read_loop

101 close(15)

end do file_iter_time

write(*,*) "Number of timesteps found:", numTsteps

xlen = xmax - xmin
ylen = ymax - ymin
zlen = zmax - zmin

! Allocate master array
allocate(mast_arr(numMols,3,numTsteps), stat=ioErr)

if (ioErr .ne. 0) then
	write(*,*) "Allocation of master array has failed. Exiting"
	stop
end if

mast_arr = -1000.0

! Read in the data to master array

time_prev = -1.0


file_iter_main: do i = 1, num_files, 1
			open(unit=15,file=file_arr(i),status="old",action="read")
			write(*,*) "Taking data from file ", file_arr(i)

! Reads in the data
out_loop: do

	read(15,*,END=103)

	read(15,*) tstep
	! Avoid duplications of time data
	if (tstep .eq. time_prev) then

		read(15,*)
		read(15,*) numMols_r
		read(15,*); read(15,*);
		read(15,*); read(15,*); read(15,*);

		! Read in molecule data
		data_in: do j = 1, numMols, 1

			read(15,*)

		end do data_in

	else
		time_prev = tstep
		read(15,*)
		read(15,*) numMols_r
		numMols = int(numMols_r)
		read(15,*); read(15,*);
		read(15,*); read(15,*); read(15,*);

		t_cur = t_cur + 1

		! Read in molecule data
		dat_in: do j = 1, numMols, 1

			read(15,*) mol_num, junk, x, y, z
			mast_arr(mol_num,1,t_cur) = x
			mast_arr(mol_num,2,t_cur) = y
			mast_arr(mol_num,3,t_cur) = z

		end do dat_in

	end if

	time_prev = tstep

end do out_loop

103 close(15)

end do file_iter_main

! See that counted timesteps match
if (t_cur .ne. numTsteps) then
	write(*,*) "Number of timsteps found does not match number of steps of &
					&		data recorded."
	stop
end if

write(*,*) "Number of data elements not read: ", count(mast_arr .eq. -1000.0)

write(*,*) "Beginning Boundary handling"

call flat(numMols,3,numTsteps,mast_arr,xmin,ymin,zmin,xmax,ymax,zmax)

write(*,*) "Data Flattened. Beginning diffusion calculations."

call diffuse(numMols,3,numTsteps,mast_arr,t_stride,lj_t)

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine flat(dim1,dim2,dim3,arr_in,xmn,ymn,zmn,xmx,ymx,zmx)
! Subroutine to refactor data to fix periodic boundary crossovers.

implicit none
	integer, intent(in)		:: dim1, dim2, dim3 ! num. beads, xyz, num T steps
	real,intent(in)				:: xmn, ymn, zmn ! Box dimensions
	real,intent(in)				:: xmx, ymx, zmx
	real,intent(inout)		:: arr_in(dim1,dim2,dim3) ! Master array

	integer								:: i, j, k
	real									:: r1, r2, d
	real									:: xbx, ybx, zbx

! Determine box dimensions
xbx = xmx - xmn
ybx = ymx - ymn
zbx = zmx - zmn

bead_loop: do j = 1, dim1, 1

	time_loop: do i = 1, dim3-1, 1

	! X part
	! First and second positions
	r1 = arr_in(j,1,i)
	r2 = arr_in(j,1,i+1)
	d = r2 - r1
	! If displacement larger than 50% of box, make diffusion shortest path
	if (abs(d) .gt. 0.5*xbx) then
		! Displacement in positive direction
		if (d .gt. 0.0) then
			! Shift positions of this and all future values in x
			do k = i+1, dim3, 1
			  arr_in(j,1,k) = arr_in(j,1,k) - xbx
			end do
			! Displacement in negative direction
		else if (d .lt. 0.0) then
			do k = i+1, dim3, 1
			  arr_in(j,1,k) = arr_in(j,1,k) + xbx
			end do
		end if
	end if

	! Y part
	! First and second positions
	r1 = arr_in(j,2,i)
	r2 = arr_in(j,2,i+1)
	d = r2 - r1
	! If displacement larger than 50% of box, make diffusion shortest path
	if (abs(d) .gt. 0.5*ybx) then
		! Displacement in positive direction
		if (d .gt. 0.0) then
			! Shift positions of this and all future values in x
			do k = i+1, dim3, 1
			  arr_in(j,2,k) = arr_in(j,2,k) - ybx
			end do
			! Displacement in negative direction
		else if (d .lt. 0.0) then
			do k = i+1, dim3, 1
			  arr_in(j,2,k) = arr_in(j,2,k) + ybx
			end do
		end if
	end if

	! Z part
	! First and second positions
	r1 = arr_in(j,3,i)
	r2 = arr_in(j,3,i+1)
	d = r2 - r1
	! If displacement larger than 50% of box, make diffusion shortest path
	if (abs(d) .gt. 0.5*zbx) then
		! Displacement in positive direction
		if (d .gt. 0.0) then
			! Shift positions of this and all future values in x
			do k = i+1, dim3, 1
			  arr_in(j,3,k) = arr_in(j,3,k) - zbx
			end do
			! Displacement in negative direction
		else if (d .lt. 0.0) then
			do k = i+1, dim3, 1
			  arr_in(j,3,k) = arr_in(j,3,k) + zbx
			end do
		end if
	end if

	end do time_loop

end do bead_loop

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine diffuse(dim1,dim2,dim3,arr_in,f_step,lj_c)
! Subroutine to find displacement of molecules. Finds average as well
! Periodic boundary handling

implicit none
	integer,intent(in)		:: dim1, dim2, dim3 ! num. beads, xyz,num T steps
	integer,intent(in)		:: f_step ! number of times to skip
	real,intent(in)			  :: lj_c  ! number of LJ timesteps in each snapshot
	real,intent(in)				:: arr_in(dim1,dim2,dim3) ! master array

	integer								:: max_tau, max_stride, dtau ! Boundaries
	integer								:: stride ! Current value of dt
	integer								:: i, j, k! Looping integers
	real									:: disp_arr_ti(dim3) ! holds disp. for each mol
	real									:: st_dev
	integer								:: num_t_chk
	real									:: disp_avg_co ! avg for beads
	real									:: disp_avg_ti ! Average for time
	real									:: disp_tot ! Total displacement from t = 0
	real									:: disp, disp_avg
	real									:: x1, x2, y1, y2, z1, z2
	real									:: xd, yd, zd ! Axial distances

	character*12					:: filename

max_tau = lj_c*dim3
max_stride = dim3/2

write(*,*) "Largest Tau available:", max_tau
write(*,*) "Number of Strides checking:", max_stride
! Each timestep is 50 tau

filename = "diff_out.dat"

open(unit=16,file=filename,status="replace",position="append")

tau_loop: do k = 1, max_stride, 1

stride = k
dtau = lj_c*k

num_t_chk = 0
disp_avg_ti = 0.0
disp_avg = 0.0
disp_tot = 0.0

	! Iterate over time
	time_loop: do i = 1, dim3-max_stride, f_step

			num_t_chk = num_t_chk + 1
			disp_avg_co = 0.0

			! Iterate over each bead
			coord_loop: do j = 1, dim1, 1

				x1 = arr_in(j,1,i)
				x2 = arr_in(j,1,i+stride)
				y1 = arr_in(j,2,i)
				y2 = arr_in(j,2,i+stride)
				z1 = arr_in(j,3,i)
				z2 = arr_in(j,3,i+stride)

				! Axial distances traveled
				xd = (x2 - x1)
				yd = (y2 - y1)
				zd = (z2 - z1)

				! Determine
				call dist(xd,yd,zd,disp)

				disp_avg_co = disp_avg_co + disp**2

			end do coord_loop

			! Average by number of beads
			disp_avg_co = disp_avg_co/dim1
			! Add to displacement by time
			disp_avg_ti = disp_avg_ti + disp_avg_co
			! Fill array
			disp_arr_ti(i) = disp_avg_co

		end do time_loop
	! Average by number of times utilized
	disp_avg_ti = disp_avg_ti/float(num_t_chk)
	! Total displacement
	disp_tot = disp_tot + disp_avg_ti

	call std_dev(num_t_chk,disp_arr_ti(1:num_t_chk),st_dev)

	write(16,*) dtau, disp_tot, st_dev
	!write(16,*) dtau, disp_avg_ti, st_dev

end do tau_loop

close(16)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine std_dev(dim,arr_in,stdev)
! Subroutine accepts a 1D real array and returns the standard deviation of the
! values

implicit none
	integer,intent(in) :: dim
	real,intent(in)		 :: arr_in(dim)
	real,intent(out)	 :: stdev

	real							 :: mean
	real							 :: r_sum
	integer						 :: i

r_sum = 0.0

mean = sum(arr_in)/float(dim)

do i = 1, dim, 1
		r_sum = r_sum + (arr_in(i)-mean)**2
end do

stdev = r_sum/float(dim)

stdev = sqrt(stdev)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dist(x,y,z,d)
! Function checks the distance of two molecules
! in 3D space using simple distance formula
! xyz mol 1; abc mol2

implicit none
		real,intent(in)    :: x, y, z
		real,intent(inout) :: d

	d = 0.0

	d = sqrt((x**2) + (y**2) + (z**2))

end subroutine

! Program takes an output file from the LAMMPS software package and
! reads in the location and type of each molecule. From this, this
! program determines the diffusion of each bead.

! This version does not track chains, only beads

! Version 3
! - Reads indefinite number of files, user will supply the number and names
! 	of the datafiles.
! - Iterates over the stride to output data as avg diffusion - time (tau)

! Now with periodic boundary handling

! Author: Jon Parsons
! Date: 2-1-19

program diffusion

implicit none
	character*50,allocatable	:: file_arr(:)
	character*50							:: raw_in
	real				  						:: tstep
	integer										:: numMols
	integer										:: numTsteps, t_cur
	integer										:: num_files
	integer										:: ioErr
	integer										:: i, j
	integer										:: mol_num
	real											:: xmin, xmax, ymin, ymax, zmin, zmax
	real											:: xlen, ylen, zlen
	real											:: x, y, z
	real											:: junk
	real,allocatable					:: mast_arr(:,:,:) ! (molNum,x-y-z,time)

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

file_iter_time: do i = 1, num_files, 1
			open(unit=15,file=file_arr(i),status="old",action="read")
			write(*,*) "Reading file ", file_arr(i)

! Reads in the data and calls the clustering subroutine
read_loop: do

	read(15,*,END=101)
	read(15,*)
	read(15,*)
	read(15,*) numMols
	read(15,*)
	read(15,*) xmin, xmax
	read(15,*) ymin, ymax
	read(15,*) zmin, zmax
	read(15,*)

	numTsteps = numTsteps + 1

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


file_iter_main: do i = 1, num_files, 1
			open(unit=15,file=file_arr(i),status="old",action="read")
			write(*,*) "Taking data from file", file_arr(i)

! Reads in the data and calls the clustering subroutine
out_loop: do

	read(15,*,END=103)
	read(15,*) tstep
	read(15,*)
	read(15,*) numMols
	read(15,*); read(15,*);
	read(15,*); read(15,*); read(15,*);

	t_cur = t_cur + 1

	! Read in molecule data
	data_in: do j = 1, numMols, 1

		read(15,*) mol_num, junk, x, y, z
		mast_arr(mol_num,1,t_cur) = x
		mast_arr(mol_num,2,t_cur) = y
		mast_arr(mol_num,3,t_cur) = z

	end do data_in

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
write(*,*) "Data read. Beginning diffusion calculations."

call diffuse(numMols,3,numTsteps,mast_arr,xlen,ylen,zlen)

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine diffuse(dim1,dim2,dim3,arr_in,xbx,ybx,zbx)
! Subroutine to find displacement of molecules. Finds average as well
! Periodic boundary handling

use functions

implicit none
	integer,intent(in)		:: dim1, dim2, dim3 ! num. beads, xyz,num T steps
	real,intent(in)				:: xbx, ybx, zbx ! Box dimensions
	real,intent(in)				:: arr_in(dim1,dim2,dim3) ! master array

	integer								:: stride, tau ! length of time, LJ time steps
	integer								:: max_tau, max_stride ! Boundaries
	integer								:: stride_iter ! Amt. to update stride
	integer								:: i, j, k ! Looping integers
	real									:: disp_avg_co, disp_avg_ti ! avg for beads, time resp.
	real									:: disp, disp_avg, disp_max
	real									:: x1, x2, y1, y2, z1, z2
	real									:: xd, yd, zd ! Axial distances
	integer								:: max_mol ! Bead with most displacement

	character*12					:: filename

max_tau = 50*dim3
write(*,*) "Largest Tau available:", max_tau
! Each timestep is 50 tau
max_stride = dim3/50

stride_iter = max_stride/10

filename = "diff_out.dat"
open(unit=16,file=filename,status="unknown",position="append")

disp_avg_co = 0.0
disp_avg_ti = 0.0
disp_avg = 0.0

! Iterate through time differences
iter_stride: do k = 1, max_stride, stride_iter

	stride = k
	tau = 50*stride
	disp_avg_ti = 0.0
	write(*,*) "Calculating for", tau

	! Iterate over time
	time_loop: do i = 1, dim3, 1

			if ((i+stride) .gt. dim3) then
				exit time_loop
			end if

			disp_avg_co = 0.0
			disp_max = 0.0

			! Iterate over each bead
			coord_loop: do j = 1, dim1, 1

				x1 = arr_in(j,1,i)
				x2 = arr_in(j,1,i+stride)
				y1 = arr_in(j,2,i)
				y2 = arr_in(j,2,i+stride)
				z1 = arr_in(j,3,i)
				z2 = arr_in(j,3,i+stride)

				! Axial distances traveled
				xd = abs(x2 - x1)
				yd = abs(y2 - y1)
				zd = abs(z2 - z1)

				! Account for periodic boundaries
				! X
				if (xd .ge. xbx) then
					xd = xd - xbx
				end if
				! Y
				if (yd .ge. ybx) then
					yd = yd - ybx
				end if
				! Z
				if (zd .ge. zbx) then
					zd = zd - zbx
				end if
				! Determine
				disp = dist(xd,yd,zd)

				if (disp .gt. disp_max) then
					max_mol = j
					disp_max = disp
				end if

				disp_avg_co = disp_avg_co + disp

			end do coord_loop

			disp_avg_co = disp_avg_co/float(dim1)
			disp_avg_ti = disp_avg_ti + disp_avg_co

		end do time_loop

	disp_avg_ti = disp_avg_ti/float(dim3-stride)

	write(16,*) tau, disp_avg_ti**2, disp_max**2

end do iter_stride

close(16)

end subroutine

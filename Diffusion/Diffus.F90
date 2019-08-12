! Program takes an output file from the LAMMPS software package and
! reads in the location and type of each molecule. From this, this
! program determines the diffusion of each bead.

! This version does not track chains, only beads

! Second attempt. Reads files to determine number of timesteps total.
! Reads in all positions at all times, runs a trace through the data.

! Now with periodic boundary handling

! Author: Jon Parsons
! Date: 2-1-19

program clusfinder

implicit none
	character*50			:: file1, file2, file3
	real				  		:: tstep
	integer						:: numMols, numChains, file_track
	integer						:: numTsteps, t_cur
	integer						:: stride
	integer						:: ioErr
	integer						:: j
	integer						:: mol_num
	real							:: xmin, xmax, ymin, ymax, zmin, zmax
	real							:: xlen, ylen, zlen
	real							:: x, y, z
	real							:: junk
	real,allocatable	:: mast_arr(:,:,:) ! (molNum,x-y-z,time)

! Number of chains in the system
numChains = 2000

write(*,*) "Please enter the desired stride."
read(*,*) stride

! pict file input
write(*,*) "Please enter the name of the file one."
write(*,*) "If the file is not in this directory enter the full path."
read(*,*) file1

file1 = trim(file1)

write(*,*) "Please enter the name of the file two."
write(*,*) "If the file is not in this directory enter the full path."
read(*,*) file2

file2 = trim(file2)

write(*,*) "Please enter the name of the file three."
write(*,*) "If the file is not in this directory enter the full path."
read(*,*) file3

file3 = trim(file3)

! This section to determine the total number of timesteps and number of molecules
numTsteps = 0
t_cur = 0
file_track = 1

101 if (file_track .eq. 1) then
			open(unit=15,file=file1,status="old",action="read")
			write(*,*) "Reading file ", file1
			file_track = 2
		else if (file_track .eq. 2) then
			close(15)
			open(unit=15,file=file2,status="old",action="read")
			write(*,*) "Reading file ", file2
			file_track = 3
		else if (file_track .eq. 3) then
			close(15)
			open(unit=15,file=file3,status="old",action="read")
			write(*,*) "Reading file ", file3
			file_track = 4
		else if (file_track .eq. 4) then
			close(15)
			goto 102
		end if

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

102 write(*,*) "Number of timesteps found:", numTsteps

xlen = xmax - xmin
ylen = ymax - ymin
zlen = zmax - zmin

! Allocate master array
allocate(mast_arr(numMols,3,numTsteps), stat=ioErr)

if (ioErr .ne. 0) then
	write(*,*) "Allocation of master array has failed. Exiting"
	stop
end if

! Read in the data to master array

file_track = 1

103 if (file_track .eq. 1) then
			open(unit=15,file=file1,status="old",action="read")
			write(*,*) "Taking data from file 1"
			file_track = 2
		else if (file_track .eq. 2) then
			close(15)
			open(unit=15,file=file2,status="old",action="read")
			write(*,*) "Taking data from file 2"
			file_track = 3
		else if (file_track .eq. 3) then
			close(15)
			open(unit=15,file=file3,status="old",action="read")
			write(*,*) "Taking data from file 3"
			file_track = 4
		else if (file_track .eq. 4) then
			close(15)
			goto 104
		end if

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
	file_in: do j = 1, numMols, 1

		read(15,*) mol_num, junk, x, y, z
		mast_arr(mol_num,1,t_cur) = x
		mast_arr(mol_num,2,t_cur) = y
		mast_arr(mol_num,3,t_cur) = z

	end do file_in

end do out_loop

! See that counted timesteps match
104 if (t_cur .ne. numTsteps) then
			write(*,*) "Number of timsteps found does not match number of steps of &
							&		data recorded."
			stop
		end if

write(*,*) "Number of data elements not read: ", count(mast_arr .eq. 0)
write(*,*) "Data read. Beginning diffusion."

call diffuse(numMols,3,numTsteps,mast_arr,stride,xlen,ylen,zlen)

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine diffuse(dim1,dim2,dim3,arr_in,strd,xbx,ybx,zbx)
! Subroutine to find displacement of molecules. Finds average as well
! Periodic boundary handling

use functions

implicit none
	integer,intent(in)		:: dim1, dim2, dim3
	integer,intent(in)		:: strd
	real,intent(in)				:: xbx, ybx, zbx
	real,intent(in)				:: arr_in(dim1,dim2,dim3)

	integer								:: i, j
	integer								:: num_pts
	real									:: disp, disp_avg, disp_max
	real									:: disp_glob_avg, disp_arr(dim3)
	real									:: x1, x2, y1, y2, z1, z2
	real									:: xd, yd, zd
	integer								:: max_mol

	character*12					:: filename

filename = "diff_out.dat"
open(unit=16,file=filename,status="unknown",position="append")
write(16,*) "Step Avg_disp Max_disp Max_mol"


disp_arr = 0.0
num_pts = 0

time_loop: do i = 1, dim3, 1

		if ((i+strd) .gt. dim3) then
			exit time_loop
		end if

		num_pts = num_pts + 1

		disp_avg = 0.0
		disp_max = 0.0

		coord_loop: do j = 1, dim1, 1

			x1 = arr_in(j,1,i)
			x2 = arr_in(j,1,i+strd)
			y1 = arr_in(j,2,i)
			y2 = arr_in(j,2,i+strd)
			z1 = arr_in(j,3,i)
			z2 = arr_in(j,3,i+strd)

			! Axial distances traveled
			xd = (x2 - x1)
			yd = (y2 - y1)
			zd = (z2 - z1)

			! Account for periodic boundaries
			! X
			if (abs(xd) .ge. xbx) then
				if (xd .gt. 0) then
						xd = xd - xbx
				else
						xd = xd + xbx
				end if
			end if
			! Y
			if (abs(yd) .ge. ybx) then
				if (yd .gt. 0) then
					yd = yd - ybx
				else
					yd = yd + ybx
				end if
			end if
			! Z
			if (abs(zd) .ge. zbx) then
				if (zd .gt. 0) then
					zd = zd - zbx
				else
					zd = zd + zbx
				end if
			end if
			! Determine
			disp = dist(xd,yd,zd)

			if (disp .gt. disp_max) then
				max_mol = j
				disp_max = disp
			end if

			disp_avg = disp_avg + disp

		end do coord_loop

		disp_avg = disp_avg/float(dim1)

		write(16,*) i, disp_avg, disp_max, max_mol
		disp_arr(i) = disp_avg

end do time_loop

write(*,*) "Number of points checked:", num_pts
write(*,*) "Total number of timesteps:", dim3

disp_glob_avg = sum(disp_arr)/float(dim3)

write(16,*) "Total average displacement: ", disp_glob_avg

close(16)

end subroutine

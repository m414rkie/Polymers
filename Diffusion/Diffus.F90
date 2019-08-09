! Program takes an output file from the LAMMPS software package and
! reads in the location and type of each molecule. From this, this
! program determines the diffusion of each bead.

! This version does not track chains, only beads

! Author: Jon Parsons
! Date: 2-1-19

program clusfinder

implicit none
	character*50			:: filename
	real				  		:: tstep
	integer						:: numMols, numChains
	integer						:: numTsteps, t_loc, t_glob
	integer						:: ioErr
	integer						:: j
	real,allocatable	:: molData(:,:) ! molnumber, moltype, x, y, z, cluster
	real,allocatable	:: avg_pos(:,:,:) ! Holds average position of each mol. Mol number, x y z t, time step
	real,allocatable	:: cur_pos(:,:,:) ! Holds current positions of each mol. Mol number, x y z t, time step

! Number of chains in the system
numChains = 2000

! pict file input
100 write(*,*) "Please enter the name of the file with the data."
write(*,*) "If the file is not in this directory enter the full path."
read(*,*) filename

filename = trim(filename)

open(unit=15, file=filename, status="old", action="read", iostat=ioErr)

if (ioErr .ne. 0) then
	write(*,*) "File not found, please try again."
	goto 100
end if

numTsteps = 0
t_loc = 0
t_glob = 0

! Reads in the data and calls the clustering subroutine
out_loop: do

	read(15,*,END=101)
	read(15,*) tstep
	read(15,*)
	read(15,*) numMols
	read(15,*); read(15,*); read(15,*); read(15,*); read(15,*);

	if (numTsteps .eq. 0) then
		allocate(avg_pos(numMols,4,1000), stat = ioErr)

		if (ioErr .ne. 0) then
			write(*,*) "Failed to allocate average position array. Exiting at timestep", tstep
			stop
		end if

		avg_pos = 0.0

	end if

	allocate(cur_pos(numMols,4,10), stat = ioErr)

	if (ioErr .ne. 0) then
		write(*,*) "Failed to allocate temp. position array."
		stop
	end if

	cur_pos = 0.0

	numTsteps = numTsteps + 1
	t_loc = t_loc + 1

	allocate(molData(numMols,5), stat = ioErr)

	if (ioErr .ne. 0) then
		write(*,*) "Failed to allocated primary array. Exiting at timestep", tstep
		stop
	end if

	write(*,*) "Beginning time:", tstep

	! Initial values, to be over-ridden
	molData = -1.0

	! Read in molecule data
	fileread: do j = 1, numMols, 1

		read(15,*) molData(j,1), molData(j,2), molData(j,3), molData(j,4), molData(j,5)

	end do fileread

	pos_save: do j = 1, numMols, 1

		cur_pos(nint(molData(j,1)),1,t_loc) = molData(j,3) ! x
		cur_pos(nint(molData(j,1)),2,t_loc) = molData(j,4) ! y
		cur_pos(nint(molData(j,1)),3,t_loc) = molData(j,5) ! z
		cur_pos(nint(molData(j,1)),4,t_loc) = tstep ! time

	end do pos_save

	deallocate(molData)

	if (t_loc .eq. 10) then
		write(*,*) "Averaging"
		t_glob = t_glob + 1
		call cur_average(numMols,4,10,1000,t_glob,cur_pos,avg_pos)
		write(*,*) "Averaging Complete"
		t_loc = 0
	end if

	deallocate(cur_pos)

end do out_loop

101 write(*,*) "Number of timesteps checked:", numTsteps

write(*,*) "Finalizing output"
call disp_out(numMols,4,1000,avg_pos,t_glob)

write(*,*) "End of input file reached. Goodbye"

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cur_average(dim1,dim2_pos,dim3_loc,dim3_glob,t_cur,arr_cur,arr_glob)
! Subroutine takes ten posisions, averages them and saves to arr_glob

implicit none
	integer, intent(in)			:: dim1, dim2_pos, dim3_loc, dim3_glob
	integer, intent(in)			:: t_cur
	real, intent(inout)			:: arr_cur(dim1,dim2_pos,dim3_loc)
	real, intent(inout)			:: arr_glob(dim1,dim2_pos,dim3_glob)

integer			:: i, j


! Average positions
mol_loop : do i = 1, dim1, 1

	avg_loop : do j = 1, dim3_loc, 1

				arr_glob(i,1,t_cur) = arr_glob(i,1,t_cur) + arr_cur(i,1,j) ! x
				arr_glob(i,2,t_cur) = arr_glob(i,2,t_cur) + arr_cur(i,2,j) ! y
				arr_glob(i,3,t_cur) = arr_glob(i,3,t_cur) + arr_cur(i,3,j) ! z
				arr_glob(i,4,t_cur) = arr_glob(i,4,t_cur) + arr_cur(i,4,j) ! t

	end do avg_loop

	arr_glob(i,1,t_cur) = arr_glob(i,1,t_cur)/10.0 ! x
 	arr_glob(i,2,t_cur) = arr_glob(i,2,t_cur)/10.0 ! y
	arr_glob(i,3,t_cur) = arr_glob(i,3,t_cur)/10.0 ! z
	arr_glob(i,4,t_cur) = arr_glob(i,4,t_cur)/10.0 ! t

end do mol_loop

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine disp_out(dim1,dim2,dim3,arrin,t_count)
! Finalization subroutine. Finds displacement on average, as well as maximal
! displacement

use functions

implicit none
	integer, intent(in) 	:: dim1, dim2, dim3
	integer, intent(in)		:: t_count
	real, intent(in)			:: arrin(dim1,dim2,dim3) ! molnum, xyzt, timestep

	integer		:: i, j
	real			:: disp_avg, disp_max
	real			:: disp_avg_all
	real			:: disp_cur


disp_max = 0.0
disp_avg_all = 0.0

open(unit=16,file="displacement.dat",status="unknown",position="append")
open(unit=17,file="maxdisplacement.dat",status="unknown",position="append")

time_loop : do i = 1, t_count, 1

	disp_max = 0.0
	disp_avg = 0.0

	mol_loop : do j = 1, dim1, 1

		disp_cur = dist(arrin(j,1,i),arrin(j,2,i),arrin(j,3,i), &
										arrin(j,1,i+1),arrin(j,2,i+1),arrin(j,3,i+1))

		disp_avg = disp_avg + disp_cur

		if (disp_cur .gt. disp_max) then
			disp_max = disp_cur
		end if

	end do mol_loop

	disp_avg = disp_avg/float(dim1)
	disp_avg_all = disp_avg_all + disp_avg


	write(16,*) arrin(1,4,i), disp_avg
	write(17,*) arrin(1,4,i), disp_max

end do time_loop

disp_avg_all = disp_avg_all/float(t_count)

write(16,*) "Total_disp", disp_avg_all

close(16)
close(17)

end subroutine

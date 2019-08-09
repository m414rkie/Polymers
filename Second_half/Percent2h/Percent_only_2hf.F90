! Program takes an output file from the LAMMPS software package and
! reads in the linked beads. The program then builds a picture of
! the linked chains

! Author: Jon Parsons
! Date: 5-21-19

program clusfinder

use functions

implicit none
	character*50				:: filename  ! Name of input file
	real								:: tstep, junk ! Distance parameter, number of clusters found, cuurent timestep
	integer							:: numBonds, numChains ! Number of molecules in system, ! Number of chains per molecule
	integer							:: numTsteps ! Number of time steps looked at
	integer							:: ioErr, j, i  ! System error variable, looping integer
	real,allocatable		:: molData(:,:)    ! molnumber, moltype, x, y, z, cluster, molgroup
	real,allocatable	  :: chainTrack(:)   ! This array tracks chain interactions
	real								:: perc, perc_total ! Current Percentage, overall percentage
	integer							:: chain_a, chain_b ! Which chains are involved
	integer							:: t_count, t_half

! Number of chains in the system
numChains = 2000

! Inititalization
perc_total = 0

! Allocation
allocate(chainTrack(numChains), stat = ioErr)

if (ioErr .ne. 0) then
	write(*,*) "Failed to allocate the chain tracking array. Exiting."
	stop
end if

100 write(*,*) "Please enter the name of the file with the data."
write(*,*) "If the file is not in this directory enter the full path."
read(*,*) filename

filename = trim(filename)

open(unit=15, file=filename, status="old", action="read", iostat=ioErr)

if (ioErr .ne. 0) then
	write(*,*) "File not found, please try again."
	goto 100
end if

! First half
t_count = 0

! Read in header data
read(15,*)
read(15,*)
read(15,*)
read(15,*) numBonds
read(15,*); read(15,*)
read(15,*); read(15,*)
read(15,*)


! Reads in the data and calls the clustering subroutine until EOF
do

	! Iterate timestep
	t_count = t_count + 1

	! Read in molecule data
	do j = 1, numBonds, 1
		read(15,*)
	end do

	! Checks for EOF, if not then reads header data for next step
	read(15,*,END=102)
	read(15,*,END=102)
	read(15,*,END=102)
	read(15,*,END=102) numBonds
	read(15,*,END=102)
	read(15,*,END=102)
	read(15,*,END=102)
	read(15,*,END=102)
	read(15,*,END=102)

end do

102 rewind(unit=15)

write(*,*) "Number of timesteps:", t_count
t_half = t_count/2

! Second half data part
! Initialize timesteps
numTsteps = 0

! Read in header data
read(15,*)
read(15,*) tstep
read(15,*)
read(15,*) numBonds
read(15,*); read(15,*)
read(15,*); read(15,*)
read(15,*)


! Reads in the data and calls the clustering subroutine until EOF
do

	! Iterate timestep
	numTsteps = numTsteps + 1

	chainTrack = 0.0
	perc = 0.0

	! Allocate working array
	allocate(molData(numBonds,2), stat = ioErr)

	if (ioErr .ne. 0) then
		write(*,*) "Failed to allocated primary array. Exiting at timestep", tstep
		stop
	end if

	write(*,*) "Beginning time:", tstep

	! Initial values, to be over-ridden
	molData = 0.0
	chainTrack = 0

	! Read in molecule data
	fileread: do j = 1, numBonds, 1

		read(15,*) molData(j,1), molData(j,2), junk

	end do fileread

	if (numTsteps .gt. t_half) then
		do i = 1, numBonds, 1

			chain_a = chain(molData(i,1))
			chain_b = chain(molData(i,2))

			chainTrack(chain_a) = 1.0
			chainTrack(chain_b) = 1.0

		end do

		perc = sum(chainTrack)/real(numChains)
		perc_total = perc_total + perc

		open(unit=18,file="Percentages_2h.dat",status="unknown",position="append")
		write(18,*) tstep, perc
		close(18)
	end if

	deallocate(molData)


	! Checks for EOF, if not then reads header data for next step
	read(15,*,END=101)
	read(15,*) tstep
	read(15,*)
	read(15,*) numBonds
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)

end do

101 perc_total = perc_total/float(t_half)

open(unit=18,file="Percentages_2h.dat",status="unknown",position="append")
write(18,*) "Overall Percentage:", perc_total
close(18)

write(*,*) "Number of timesteps checked:", t_half
write(*,*) "End of input file reached. Goodbye"

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

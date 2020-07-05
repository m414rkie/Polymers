! Program takes an output file from the LAMMPS software package and
! reads in the linked beads from a bonds file. The program then builds a picture
! of the clusters of chains formed.

! Author: Jon Parsons
! Date: 5-21-19

program clusfinder

implicit none
	character*50				:: filename  ! Name of input file
	real								:: tstep, junk ! Distance parameter, number of clusters found, cuurent timestep
	real								:: perc ! Holds overall percentage of chains involved in a cluster
	integer							:: numBonds, numChains ! Number of molecules in system, ! Number of chains per molecule
	integer							:: numTsteps ! Number of time steps looked at
	integer							:: first,ClusFound	! Tracks if this is first time output subroutine is called
	integer							:: ioErr, j  ! System error variable, looping integer
	real,allocatable		:: molData(:,:)    ! molnumber, moltype, x, y, z, cluster, molgroup
	integer							:: statsArr(2000)	   ! Tracks number of clusters with (2,3,4...) chains involved
	integer,allocatable	:: chainTrack(:)   ! This array tracks chain interactions

! Number of chains in the system
numChains = 2000

! Inititalization
first = 0
perc = 0

! Allocation
allocate(chainTrack(numChains), stat = ioErr)

if (ioErr .ne. 0) then
	write(*,*) "Failed to allocate the chain tracking array. Exiting."
	stop
end if

100 write(*,*) "Please enter the name of the file with the data."
write(*,*) "Files will typically include the 'bonds' prefix"
write(*,*) "If the file is not in this directory enter the full path."
read(*,*) filename

filename = trim(filename)

open(unit=15, file=filename, status="old", action="read", iostat=ioErr)

if (ioErr .ne. 0) then
	write(*,*) "File not found, please try again."
	goto 100
end if

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
	statsArr = 0

	! Read in molecule data
	fileread: do j = 1, numBonds, 1

		read(15,*) molData(j,1), molData(j,2), junk

	end do fileread

	call ClusBuilder(molData,chainTrack,numBonds,2,numChains,ClusFound)
	call output(chainTrack,numChains,first,tstep,ClusFound,statsArr,perc)
	call statistics(statsArr,first,numTsteps,perc)

	! Set flag to indicate no longer first timestep
	first = 1

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

101 close(15)
first = 2 ! set flag to final output

call statistics(statsArr,first,numTsteps,perc)

write(*,*) "Number of timesteps checked:", numTsteps
write(*,*) "End of input file reached. Goodbye"

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ClusBuilder(beadIn,chainIn,dim1,dim2,chDim,maxClus)
! Builds the clusters for the current timestep

use functions
use chain_Functions

implicit none
	integer,intent(in)		:: dim1,dim2,chdim ! Dimensions of arrays
	real,intent(in)			  :: beadIn(dim1,dim2) ! Bead array
	integer,intent(inout)	:: chainIn(chDim)	! Track chains assigned to clusters
	integer,intent(out)		:: maxClus ! Number of unique clusters found

	integer					:: curClus ! Current working cluster
	integer					:: i ! Looping integer

! Initialize
maxClus = 1
curClus = 1

! Build the clusters
ClusLoop : do i = 1, dim1, 1

				! Set working cluster to largest cluster
				curClus = maxClus

				! If both beads are in same chain, skip
				if (chain(beadIn(i,1)) .eq. chain(beadIn(i,2))) then
					cycle ClusLoop
				end if

				! If both chains are not already involved in a cluster, assign them to new cluster
				if ((chainIn(chain(beadIn(i,1))) .eq. 0).and.(chainIn(chain(beadIn(i,2))) .eq. 0)) then
					chainIn(chain(beadIn(i,1))) = curClus
					chainIn(chain(beadIn(i,2))) = curClus
					maxClus = maxClus + 1 ! Update max number of clusters
					cycle ClusLoop
				end if

				! If first chain is in cluster, assign second chain to existing cluster
				if ((chainIn(chain(beadIn(i,1))) .ne. 0).and.(chainIn(chain(beadIn(i,2))) .eq. 0)) then
					curClus = chainIn(chain(beadIn(i,1)))
					chainIn(chain(beadIn(i,2))) = curClus
					cycle ClusLoop
				end if

				! If second chain is in cluster, assign first chain to existing cluster
				if ((chainIn(chain(beadIn(i,2))) .ne. 0).and.(chainIn(chain(beadIn(i,1))) .eq. 0)) then
					curClus = chainIn(chain(beadIn(i,2)))
					chainIn(chain(beadIn(i,1))) = curClus
					cycle ClusLoop
				end if

				! If both chains are in a cluster, assign all chains in second cluster to first cluster
				if ((chainIn(chain(beadIn(i,1))) .ne. 0).and.(chainIn(chain(beadIn(i,2))) .ne. 0)) then
					where (chainIn .eq. chainIn(chain(beadIn(i,2)))) chainIn = chainIn(chain(beadIn(i,1)))
				end if

end do ClusLoop

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine output(chainIn,chDim,flag,t,maxClus,numchains,perc_tot)
! Prints to output files

implicit none
	integer,intent(in)	  :: chDim ! Dimension, number of chains in system
	integer,intent(in)	  :: chainIn(chDim) ! Working array, row values are chains, column values are cluster value
	integer,intent(in)	  :: flag, maxClus ! Flag for first call, maximum number of clusters found
	real,intent(in)		    :: t ! Current timestep
	integer,intent(inout)	:: numchains(2000) ! Array for holding how many clusters have 2, 3, etc chains
	real,intent(inout)		:: perc_tot

	integer				:: i, j ! Looping integers
	integer				:: chainsIn, cluster ! temporary values for holding
	logical				:: chain_count(chDim) ! Will hold number of chains involved for counting
	logical				:: chain_no_clus(chDim) ! holds chains that are not in a cluster
	integer				:: in_count ! Collapsed number of chains involved
	real					:: perc_chains ! Percentage of chains
	integer				:: num_not_assoc ! number of chains not associated with a cluster

! Open files
open(unit=11,file="Clusters.dat",status="unknown",position="append")
open(unit=12,file="ClusSizes.dat",status="unknown",position="append")

! First call write headers
if (flag .eq. 0) then
	write(11,*) "TimeStep	NumCLusters Percentage"
end if

write(12,*) "Timestep: ", t
write(12,*) "Cluster	Chains"

cluster = 1

! Find the percentage of chains in a cluster
chain_count = (chainIn .ne. 0)
in_count = count(chain_count)
perc_chains = float(in_count)/float(chDim)
perc_tot = perc_tot + perc_chains
! find number of chains not in a cluster, assign to the output
chain_no_clus = (chainIn .eq. 0)
num_not_assoc = count(chain_no_clus)
numchains(1) = num_not_assoc
! Iterate through the clusters
ReadLoop: do i = 1, maxClus, 1

	! Initialize how many chains ber cluster
	chainsIn = 0

	! Find how many chains in current cluster
	ClusFndLoop: do j = 1, chDim, 1

					! Add one for each chain
					if (chainIn(j) .eq. i) then
						chainsIn = chainsIn + 1
					end if

	end do ClusFndLoop

	! If no chains, skip
	if (chainsIn .gt. 1) then
		! Write output
		write(12,*) cluster, chainsIn
		cluster = cluster + 1 ! Iterate cluster number for writing

		! If more than 14 chains in cluster, assign to 14+ box, inform user
		if (chainsIn .ge. 2000) then
			numchains(2000) = numchains(2000) + 1
		! All else write to appropriate bo
	else if ((chainsIn .ge. 2).and.(chainsIn .le. 2000)) then
			numchains(chainsIn) = numchains(chainsIn) + 1
		else
			write(*,*) "Cluster with chains outside of parameters found."
		end if

	end if

	! Write timestep and total number of clusters to file
	if (i .eq. maxClus) then
		write(11,*) t, (cluster-1), perc_chains
	end if

end do ReadLoop

close(11)
close(12)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine statistics(statsArr,flag,numsteps,perc_tot)
! Find some statistics about clusters

implicit none
	integer,intent(in)		:: statsArr(2000) ! Holds number of clusters of size 1, 2 etc
	integer,intent(in)		:: flag ! First time flag
	integer,intent(in)		:: numsteps ! Number of timesteps
	real									:: perc_tot ! Overall percentage of chains in a cluster

	integer					:: i ! Looping integer

	integer					:: n ! Will hold total number of clusters
	real, save 		  :: num, num_sq ! total chains and total chains squared
	real, save		  :: mean, std_dev, variance
	real, save			:: unbound_avg
	real, save		  :: numavg(2000) ! For printing

! first call initializations
if (flag .eq. 0) then
	num = 0
	num_sq = 0
	mean = 0
	std_dev = 0
	variance = 0
	numavg = 0
	unbound_avg = 0
end if
! Value updates
if (flag .le. 1) then
	n = 0
	num = 0.0
	num_sq = 0.0
	
	! Global stuff
	! x
	unbound_avg = unbound_avg + statsArr(1)
	numavg(1) = numavg(1) + statsArr(1)
	do i = 2, 2000, 1
		num = num + float((i)*statsArr(i))
		numavg(i) = numavg(i) + statsArr(i)
		n = n + statsArr(i)
	end do

	! x^2
	do i = 2, 2000, 1
		num_sq = num_sq + (float((i)*statsArr(i))**2)
	end do

	! mean
	mean = mean + num/float(n)

	! variance
	variance = variance + (num_sq - (num*num/float(n)))/float(n)

	write(*,*) "Summed mean", mean, n, variance

end if
! Final call only
if (flag .gt. 1) then

	write(*,*) "Statistics Output in file 'Averages.dat'"

	open(unit=13,file="Averages.dat",status="replace",position="append")

	std_dev = sqrt(variance)

	write(13,*) "Total Time Steps: ", numsteps
	write(13,*) "Average Size: ",(mean/float(numsteps))
	write(*,*) "Average Size: ", (mean/float(numsteps))
	write(13,*) "Std devation: ", (std_dev/float(numsteps))
	write(13,*) "Avg. Percentage: ", (perc_tot/float(numsteps))
	write(13,*) "Avg. Unbonded Chains: ", (unbound_avg/float(numsteps))
	! Output for box plot
	write(13,*) "Box Plot"
	do i = 1, 2000, 1
		write(13,*) i, numavg(i)/float(numsteps)
	end do
	close(13)

end if
end subroutine

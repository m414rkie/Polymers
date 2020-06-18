! Program takes an output file from the LAMMPS software package and
! reads in the location and type of each molecule. From this, this
! program attempts to determine which molecules are clustered together.

! This version does not track chains, only beads

! Author: Jon Parsons
! Date: 2-1-19

program clusfinder

implicit none
	character*50	:: filename
	real				  :: sigma, numClus, tstep
	integer				:: numMols, numChains
	integer				:: numTsteps
	integer				:: ioErr
	integer				:: j, first
	real,allocatable	:: molData(:,:)    ! molnumber, moltype, x, y, z, cluster

! Number of chains in the system
numChains = 2000

first = 0

! User input and initializations
write(*,*) "Enter sigma:"
read(*,*) sigma

100 write(*,*) "Please enter the name of the file with the data."
write(*,*) "Typical files will begin with the 'pict' prefix."
write(*,*) "If the file is not in this directory enter the full path."
read(*,*) filename

filename = trim(filename)

open(unit=15, file=filename, status="old", action="read", iostat=ioErr)

if (ioErr .ne. 0) then
	write(*,*) "File not found, please try again."
	goto 100
end if

numTsteps = 0

! Read in header data
read(15,*)
read(15,*) tstep
read(15,*)
read(15,*) numMols
read(15,*); read(15,*); read(15,*); read(15,*); read(15,*);

! Reads in the data and calls the clustering subroutine
do

	numTsteps = numTsteps + 1

	allocate(molData(numMols,6), stat = ioErr)

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

	call clusSort(molData,numMols,6,sigma,numClus,tstep)
	call hist_sort(molData,numMols,6,numClus,tstep)

	if (int(numClus) .ge. 1) then
		call output(molData,numMols,6,numClus,tstep,first)
	end if

	first = 1
!	call averages(numMols,6,molData,numClus,tstep)

	deallocate(molData)

	read(15,*,END=101)
	read(15,*) tstep
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)

end do

101 write(*,*) "Number of timesteps checked:", numTsteps
write(*,*) "End of input file reached. Goodbye"

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine clusSort(arrin,dim1,dim2,dParam,maxClus,t)
! Subroutine to determine which molecules of the appropriate type are
! joined together in a cluster.

use functions
use chain_Functions

implicit none
	integer						 :: dim1,dim2 ! Dimensions of input array
	real,intent(inout) :: arrin(dim1,dim2) ! Input array, holds output in last column
	real,intent(in)		 :: dParam,  ! Distance parameter

	integer					 	 :: i, j, match ! Looping and the match integer
	real							 :: whichClus, maxClus ! Determines the cluster number and highest cluster found
	real						   :: near
	character*12			 :: fileout

maxClus = 1.0

fileout = "hbonds.dat"

open(unit=19,file=trim(fileout),status="unknown",position="append")

write(19,*) "Timestep:", t

! Nested do loops iterate through the data. When conditions are met assigns the
! current molecules to a cluster
primLoop: do i = 1, dim1, 1

	match = 0

	! Outputs to show that it is still working, will suppress if found to be
	! excessively expensive
	if (i .eq. floor(0.1*real(dim1))) then
		write(*,*) "Clustering 10% complete"
	else if (i .eq. floor(0.2*real(dim1))) then
		write(*,*) "Clustering 20% complete"
	else if (i .eq. floor(0.2*real(dim1))) then
		write(*,*) "Clustering 20% complete"
	else if (i .eq. floor(0.3*real(dim1))) then
		write(*,*) "Clustering 30% complete"
	else if (i .eq. floor(0.4*real(dim1))) then
		write(*,*) "Clustering 40% complete"
	else if (i .eq. floor(0.5*real(dim1))) then
		write(*,*) "Clustering 50% complete"
	else if (i .eq. floor(0.6*real(dim1))) then
		write(*,*) "Clustering 60% complete"
	else if (i .eq. floor(0.7*real(dim1))) then
		write(*,*) "Clustering 70% complete"
	else if (i .eq. floor(0.8*real(dim1))) then
		write(*,*) "Clustering 80% complete"
	else if (i .eq. floor(0.9*real(dim1))) then
		write(*,*) "Clustering 90% complete"
	end if


	! Checks for correct molecule type ; 2 := H
	if (nint(arrin(i,2)) .ne. 2) then
		cycle primLoop
	end if

	! Checks if already in a cluster, if yes then sets that to the current cluster
	! If not, starts a new cluster
	if (nint(arrin(i,dim2)) .gt. 0) then
		whichClus = arrin(i,dim2)
	else
		whichClus = maxClus
	end if

	matchLoop: do j = dim1, 1, -1

		! Checks that each are on a different chain
		if (chain(arrin(i,1)) .eq. chain(arrin(j,1))) then
			cycle matchLoop
		end if

		! Checks for correct molecule type
		if (nint(arrin(j,2)) .ne. 3) then
			cycle matchLoop
		end if

		! Checks if already in a cluster
		if (nint(arrin(j,dim2)) .gt. 0) then
			cycle matchLoop
		end if

		! Checks the distance, final check
		near = dist(arrin(i,3),arrin(i,4),arrin(i,5),arrin(j,3),arrin(j,4),arrin(j,5))
		if (near .le. dParam) then
			arrin(j,dim2) = whichClus
			arrin(i,dim2) = whichClus
			match = 1
			write(19,*) arrin(j,1), arrin(i,1)
		end if

	end do matchLoop

	! If a new cluster is found updates the cluster count
	if ((match .ne. 0).and.(whichClus .eq. maxClus)) then
		maxClus = maxClus + 1.0
	end if

end do primLoop

close(19)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine output(arrin,dim1,dim2,Nclu,t,check)

use functions

implicit none
	integer,intent(in)				:: dim1, dim2
	real,intent(in)					:: arrin(dim1,dim2)
	real,intent(in)					:: t
	real,intent(in)					:: Nclu
	character*25					:: fileName, filename2
	logical							:: file_exist, file_exist2
	integer							:: check

	integer							:: beadNum
	integer							:: molsIn, twntyTimes
	integer							:: i, j, k
	integer							:: beadTimes(1:80002)
	logical							:: beadCheck(1:80002)

write(*,*) "Total number of Molecules:", dim1
write(*,*) "Timestep:", t
write(*,*) "Number of Clusters found:", int(Nclu)

fileName = "ClustersandBeads_hbond.dat"
fileName2 = "ClustTime_hbond.dat"

open(unit=17,file=trim(fileName),status="unknown",position="append")
open(unit=18,file=trim(fileName2),status="unknown",position="append")

write(17,*) "Timestep ", t
write(17,*) "Clusters ", int(Nclu)
write(17,*) "Cluster	Size"


if (check .eq. 0) then
		write(18,*) "Timestep	Number"
end if

beadTimes = 0
twntyTimes = 0

! Counts and outputs number of clusters found. Inner loop to determine which chains are involved.
clustercount: do i = 1, int(Nclu), 1

	molsIn = 0

	! Determines number of molecules in the i-th cluster
	molcount: do j = 1, dim1, 1
		if (int(arrin(j,dim2)) .eq. i) then
			molsIn = molsIn + 1
		end if
	end do molcount

	! Rejects cluster if more than 20 molecules found in the cluster
	if (molsIn .gt. 20) then
		twntyTimes = twntyTimes + 1
	end if

	if (molsIn .gt. 1) then

		write(17,*)"Cluster:", i, "With:", molsIn
		write(17,*) "Beads:"

		! Determines which chains are a part of the cluster
		chainfind: do k = 1, dim1, 1

			if (int(arrin(k,dim2)) .eq. i) then
				beadNum = k
				write(17,'(1i5," ")',ADVANCE="no") k
				beadTimes(beadNum) = beadTimes(beadNum) + 1
			end if

		end do chainfind

		write(17,*)

	end if

end do clustercount

write(18,*) t, int(Nclu)

beadCheck = (beadTimes .gt. 1)

write(*,*) "Number of times a bead is found more than twice:", count(beadCheck)

write(*,*) "Number of clusters with more than 20 groups:", twntyTimes

beadCheck = (beadTimes .gt. 0)

write(*,*) "Percentage of chains in a cluster:", float(count(beadCheck))/80002.0

close(17)
close(18)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine averages(dim1,dim2,arrin,nclu,t)

implicit none
	integer,intent(in)			:: dim1, dim2
	real,intent(in)				:: arrin(dim1,dim2)
	real,intent(in)				:: nclu
	real,intent(in)				:: t

	integer						:: i
	real						:: molcount

open(unit=19,file="averages_hbond.dat",status="unknown",position="append")

molcount = 0.0

do i = 1, dim1, 1

		if (int(arrin(1,dim2)) .ne. -1) then
			molcount = molcount + 1.0
		end if

end do

write(19,*) t, molcount/nclu

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine hist_sort(arrin,dim1,dim2,max_clus,t)

implicit none
	integer,intent(in)		:: dim1, dim1
	real,intent(in)				:: arrin(dim1,dim2)
	integer,intent(in)		:: max_clus
	real,intent(in)				:: t

	integer								:: i, j
	integer								:: num_h
	integer								:: hist(10)

hist = 0

do i = 1, max_clus, 1

	num_h = 0

	do j = 1, dim1, 1

		if (nint(arrin(j,dim2)) .eq. i) then
			num_h = num_h + 1
		end if

	end do

	hist(num_h) = hist(num_h) + 1

end do

open(unit=20,file="hbond_distributions.dat",status="unknown",position="append")
write(20,*) "Hbond distribution at time:", t
do i = 1, 10, 1

	write(20,*) i, hist(i)

end do

close(20)

end subroutine

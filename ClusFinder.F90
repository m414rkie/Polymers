! Program takes an output file from the LAMMPS software package and 
! reads in the location and type of each molecule. From this, this
! program attempts to determine which molecules are clustered together.

! Author: Jon Parsons
! Date: 2-1-19

! Version 2.0
! Includes subroutine to output a selection of chains for graphing.

program clusfinder

use functions

implicit none
	character*50		:: filename  ! Name of input file
	real				:: sigma, tstep ! Distance parameter, number of clusters found, cuurent timestep
	integer				:: numMols, numClus, numChains ! Number of molecules in system, ! Number of chains per molecule
	integer				:: numMG ! Number of bead groups
	real				:: boxDim(3,2) ! Min and Max values of x,y,z respectively
	integer				:: numTsteps ! Number of time steps looked at 
	integer				:: ioErr, j  ! System error variable, looping integer
	integer				:: first     ! Checks if this is the first timestep, for file handling purposes 
	real,allocatable	:: molData(:,:)    ! molnumber, moltype, x, y, z, cluster, molgroup
	integer,allocatable	:: chainTrack(:)   ! This array will track ends of chains 
										   ! and if they are in a cluster.

! Number of chains in the system
numChains = 2000

! Inititalization
first = 0

! Allocation
allocate(chainTrack(numChains), stat = ioErr)

if (ioErr .ne. 0) then
	write(*,*) "Failed to allocate the chain tracking array. Exiting."
	stop
end if

! User input and initializations
write(*,*) "Enter sigma:"
read(*,*) sigma

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

! Read in header data
read(15,*)
read(15,*) tstep
read(15,*) 
read(15,*) numMols
read(15,*)
read(15,*) boxDim(1,1), boxDim(1,2)
read(15,*) boxDim(2,1), boxDim(2,2)
read(15,*) boxDim(3,1), boxDim(3,2)
read(15,*)

! Reads in the data and calls the clustering subroutine until EOF
do
	
	numTsteps = numTsteps + 1
		
	allocate(molData(numMols,7), stat = ioErr)

	if (ioErr .ne. 0) then
		write(*,*) "Failed to allocated primary array. Exiting at timestep", tstep
		stop
	end if
			
	write(*,*) "Beginning time:", tstep

	! Initial values, to be over-ridden
	molData = 0.0
	chainTrack = -1
	
	! Read in molecule data
	fileread: do j = 1, numMols, 1
		
		read(15,*) molData(j,1), molData(j,2), molData(j,3), molData(j,4), molData(j,5)
				
	end do fileread
	call clusSort(molData,chainTrack,numMols,7,numChains,sigma,boxDim,numClus,numMG)
	call output(molData,chainTrack,numMols,7,numChains,numClus,tstep,first)
	call statspack(molData,chainTrack,numMols,7,numChains,numClus,numMG,tstep)
	
	first = 1
	
	deallocate(molData)

	! Checks for EOF, if not then reads header data for next step
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

subroutine clusSort(arrin,arrin2,dim1,dim2,dimch,dParam,boxArr,maxClus,molGroup)
! Subroutine to determine which molecules of the appropriate type are 
! joined together in a cluster.

use functions

implicit none
	integer					:: dim1,dim2 ! Dimensions of input array
	integer					:: dimch ! Dimension of chain vector (4000)
	real,intent(inout)		:: arrin(dim1,dim2) ! Input array, holds output in last column
	integer,intent(inout)	:: arrin2(dimch) ! Holds which chains are in use and at which end.
	real,intent(in)			:: dParam ! Distance parameter sigma
	real,intent(in)			:: boxArr(3,2) ! Holds max/min vals of box 
	
	integer					:: i, j, match ! Looping and the match integer
	integer					:: whichClus, maxClus ! Determines the cluster number and highest cluster found
	integer					:: molGroup, whichGroup	 ! Tracks number of beads attachec to each other directly
	real					:: near, maxdist ! Holds the distance of a bead from another bead. Maximum distance box allows 
	real					:: axisDist, boxLen ! distance of two molecule on a given axis, size of box for a given dimension

maxClus = 1
molGroup = 1
arrin2 = 0

maxDist = dist(boxArr(1,1),boxArr(1,2),boxArr(2,1),boxArr(2,2),boxArr(3,1),boxArr(3,2))

! Nested do loops iterate through the data. When conditions are met assigns the 
! current molecules to a cluster
primLoop: do i = 1, dim1, 1
	
	! Outputs to show progress
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
	
	
	! Checks for correct molecule type
	if (nint(arrin(i,2)) .ne. 3) then
		cycle primLoop
	end if
		
	! Checks if already in a cluster, if yes then sets that to the current cluster
	! If not, starts a new cluster
	if (arrin(i,dim2-1) .gt. 0) then
		whichClus = nint(arrin(i,dim2))
	end if
	
	matchLoop: do j = dim1, i, -1 
			match = 0
		
		! Checks that each are on a different chain
		if (chain(arrin(i,1)) .eq. chain(arrin(j,1))) then
			cycle matchLoop
		end if
				
		! Checks for correct molecule type 
		if (nint(arrin(j,2)) .ne. 3) then
			cycle matchLoop
		end if
		
		! If j-th mol is in a cluster, but i-th isnt, set to j's cluster
		if ((nint(arrin(i,dim2-1)).eq.0).and.(nint(arrin(j,dim2-1)).gt.0)) then
			whichClus = nint(arrin(j,dim2-1))
		! If both are in different clusters, cycle
		else if ((nint(arrin(i,dim2-1)).gt.0).and.(nint(arrin(j,dim2-1)).gt.0).and.(nint(arrin(i,dim2-1)).ne.nint(arrin(j,dim2-1)))) then
			cycle matchLoop
		! If bith are in the same cluster already, cycle
		else if ((nint(arrin(i,dim2-1)).ne.0).and.(nint(arrin(i,dim2-1)).eq.nint(arrin(j,dim2-1)))) then
			cycle matchLoop
		end if
		
		! If neither are in a cluster, start new cluster
		if ((nint(arrin(i,dim2-1)).eq.0).and.(nint(arrin(j,dim2-1)).eq.0)) then
			whichClus = maxClus
		end if
		
		! Check distance and set appropriately
		near = dist(arrin(i,3),arrin(i,4),arrin(i,5),arrin(j,3),arrin(j,4),arrin(j,5))
		if (near .le. dparam) then
						
			! If neither are in a cluster, start new cluster
			if ((nint(arrin(i,dim2-1)).eq.0).and.(nint(arrin(j,dim2-1)).eq.0)) then
				whichClus = maxClus
				maxClus = maxClus + 1
			end if
		
			call clusassign(arrin,arrin2,dim1,dim2,dimch,i,j,whichClus)
			
			! If neither are in a molGroup, start new group
			if ((nint(arrin(i,dim2)).eq.0).and.(nint(arrin(j,dim2)).eq.0)) then
				whichgroup = molGroup
				molGroup = molGroup + 1
				arrin(i,dim2) = float(whichgroup)
				arrin(j,dim2) = float(whichgroup)
			! If there is a molGroup in progress, attach j-th mol to it
			else if (nint(arrin(i,dim2)) .ne. 0) then
				whichgroup = nint(arrin(i,dim2))
				arrin(j,dim2) = float(whichgroup)
			end if
			
			cycle matchLoop
				
		end if 		
!! Periodic Boundary check !!
		! X dimension check
		axisDist = arrin(j,3) - arrin(i,3)
		boxLen = boxArr(1,2) - boxArr(1,1)
			! Below box boundaries
			if ((axisDist .gt. 0.75*boxLen).and.(arrin(i,3) .lt. 0)) then
			
				near = dist(arrin(i,3)+boxLen,arrin(i,4),arrin(i,5),arrin(j,3),arrin(j,4),arrin(j,5))
				if (near .le. dparam) then
								
					! If neither are in a cluster, start new cluster
					if ((nint(arrin(i,dim2-1)).eq.0).and.(nint(arrin(j,dim2-1)).eq.0)) then
						whichClus = maxClus
						maxClus = maxClus + 1
					end if
			
					call clusassign(arrin,arrin2,dim1,dim2,dimch,i,j,whichClus)
					
					cycle matchLoop
					
				end if
			! Above box boundaries
			else if	((axisDist .gt. 0.75*boxLen).and.(arrin(i,3) .gt. 0)) then
				
				near = dist(arrin(i,3)-boxLen,arrin(i,4),arrin(i,5),arrin(j,3),arrin(j,4),arrin(j,5))
				if (near .le. dparam ) then
													
					! If neither are in a cluster, start new cluster
					if ((nint(arrin(i,dim2-1)).eq.0).and.(nint(arrin(j,dim2-1)).eq.0)) then
						whichClus = maxClus
						maxClus = maxClus + 1
					end if
			
					call clusassign(arrin,arrin2,dim1,dim2,dimch,i,j,whichClus)
					
					cycle matchLoop
					
				end if
			end if
		! Y dimension check
		axisDist = arrin(j,4) - arrin(i,4)
		boxLen = boxArr(2,2) - boxArr(2,1)
			! Below box boundaries
			if ((axisDist .gt. 0.75*boxLen).and.(arrin(i,4) .lt. 0)) then
			
				near = dist(arrin(i,3),arrin(i,4)+boxlen,arrin(i,5),arrin(j,3),arrin(j,4),arrin(j,5))
				if (near .le. dparam) then
								
					! If neither are in a cluster, start new cluster
					if ((nint(arrin(i,dim2-1)).eq.0).and.(nint(arrin(j,dim2-1)).eq.0)) then
						whichClus = maxClus
						maxClus = maxClus + 1
					end if
			
					call clusassign(arrin,arrin2,dim1,dim2,dimch,i,j,whichClus)
					
					cycle matchLoop
					
				end if
			! Above box boundaries
			else if	((axisDist .gt. 0.75*boxLen).and.(arrin(i,4) .gt. 0)) then
				
				near = dist(arrin(i,3),arrin(i,4)-boxLen,arrin(i,5),arrin(j,3),arrin(j,4),arrin(j,5))
				if (near .le. dparam ) then
													
					! If neither are in a cluster, start new cluster
					if ((nint(arrin(i,dim2-1)).eq.0).and.(nint(arrin(j,dim2-1)).eq.0)) then
						whichClus = maxClus
						maxClus = maxClus + 1
					end if
			
					call clusassign(arrin,arrin2,dim1,dim2,dimch,i,j,whichClus)
					
					cycle matchLoop
				
				end if
			end if	
		! Z dimension check
		axisDist = arrin(j,5) - arrin(i,5)
		boxLen = boxArr(3,2) - boxArr(3,1)
			! Below box boundaries
			if ((axisDist .gt. 0.75*boxLen).and.(arrin(i,5) .lt. 0)) then
			
				near = dist(arrin(i,3),arrin(i,4),arrin(i,5)+boxLen,arrin(j,3),arrin(j,4),arrin(j,5))
				if (near .le. dparam) then
								
					! If neither are in a cluster, start new cluster
					if ((nint(arrin(i,dim2-1)).eq.0).and.(nint(arrin(j,dim2-1)).eq.0)) then
						whichClus = maxClus
						maxClus = maxClus + 1
					end if
			
					call clusassign(arrin,arrin2,dim1,dim2,dimch,i,j,whichClus)
					
					cycle matchLoop
					
				end if
			! Above box boundaries
			else if	((axisDist .gt. 0.75*boxLen).and.(arrin(i,5) .gt. 0)) then
				
				near = dist(arrin(i,3),arrin(i,4),arrin(i,5)-boxLen,arrin(j,3),arrin(j,4),arrin(j,5))
				if (near .le. dparam ) then
													
					! If neither are in a cluster, start new cluster
					if ((nint(arrin(i,dim2-1)).eq.0).and.(nint(arrin(j,dim2-1)).eq.0)) then
						whichClus = maxClus
						maxClus = maxClus + 1
					end if
			
					call clusassign(arrin,arrin2,dim1,dim2,dimch,i,j,whichClus)
					
					cycle matchLoop
				
				end if
				
			end if	
									
	end do matchLoop

end do primLoop

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine clusassign(arrin,arrin2,dim1,dim2,dimch,mol1,mol2,clus)
! Assigns the cluster number to involved molecules and chains

use functions

implicit none
	integer,intent(in)				:: dim1, dim2 ! Dimensions of input array
	integer,intent(in)				:: dimch ! Dimensions of chain-tracking array
	real,intent(inout)				:: arrin(dim1,dim2) ! Primary array. Holds bead location and cluster data.
	integer,intent(inout)			:: arrin2(dimch) ! Chain array. Holds which chains are in a cluster.
	integer,intent(in)				:: mol1, mol2,clus ! Which molecules are involved as well as cluster number
	
	integer							:: i ! looping integer
	integer							:: icount, jcount ! Count number of molecules found on each chain when assigning cluster

arrin(mol1,dim2-1) = float(clus)
arrin(mol2,dim2-1) = float(clus)
arrin2(chain(arrin(mol1,1))) = clus
arrin2(chain(arrin(mol2,1))) = clus
icount = 0
jcount = 0
! Run through and set both chains to correct cluster. 
chSet: do i = 1, dim1, 1
	if ((icount.lt.40).and.(chain(arrin(i,1)).eq.chain(arrin(mol1,1)))) then
		arrin(i,dim2-1) = float(clus)
		icount = icount + 1
	end if
	
	if ((jcount.lt.40).and.(chain(arrin(i,1)).eq.chain(arrin(mol2,1)))) then
		arrin(i,dim2-1) = float(clus)
		jcount = jcount + 1
	end if

	if ((icount.eq.40).and.(jcount.eq.40)) then
		exit chSet
	end if
end do chSet

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine output(arrin,arrin2,dim1,dim2,dimch,Nclu,t,check)
! Processes data for output.


use functions

implicit none
	integer,intent(in)				:: dim1, dim2 ! Dimensions of input array
	integer,intent(in)				:: dimch ! Dimensions of chain-tracking array
	real,intent(in)					:: arrin(dim1,dim2) ! Primary array. Holds bead location and cluster data.
	integer,intent(in)				:: arrin2(dimch) ! Chain array. Holds which chains are in a cluster.
	real,intent(in)					:: t ! Current timestep
	integer,intent(in)				:: Nclu ! nuber of cluster at this timestep
	character*25					:: fileName, filename2 ! Filename outputs
	integer							:: check ! If this is the first call to this subroutine this is 0

	integer							:: chnsIn ! Number of chains in a particular cluster and counts 
											
	integer							:: i, j, k ! Looping integers
	integer							:: chainTimes(dimch) ! Counts number of times a chain is found in a cluster
	logical							:: chainCheck(dimch) ! Logical version of the above for analysis reasons.

	integer							:: grphout	! determines if a graph output has been made for a cluster.
												
! Initialization
grphout = 0

! Output to screen
write(*,*) "Timestep:", t
write(*,*) "Number of Clusters found:", int(Nclu)

! Standard filenames
fileName = "ClustersandChains.dat"
fileName2 = "ClustTime.dat"

! File open statements
open(unit=17,file=trim(fileName),status="unknown",position="append")
open(unit=18,file=trim(fileName2),status="unknown",position="append")

! Write timestep information to file 1
write(17,*) "Timestep ", t
write(17,*) "Clusters ", Nclu
write(17,*) "Cluster	Size"

! Writes header data to file 2
if (check .lt. 1) then
		write(18,*) "Timestep	Number"
end if

! Write output to file 2
write(18,*) t, Nclu

! Initialization of variables
chainTimes = 0

! Counts and outputs number of clusters found. Inner loop to determine which chains are involved.
clustercount: do i = 1, Nclu, 1

	chnsIn = 0
	
	! Determines number of chains in the i-th cluster
	chcount: do j = 1, dimch, 1
		if (arrin2(j) .eq. i) then
			chnsIn = chnsIn + 1
		end if
	end do chcount

	If (chnsIn .gt. 1) then
		write(17,*)"Cluster:", i, "With ", chnsIn, "Chains"
		write(17,*) "Chains:"
	
		! Determines which chains are a part of the cluster
		chainfind: do k = 1, dimch, 1
		
			if (arrin2(k) .eq. i) then
				write(17,'(1i4," ")',ADVANCE="no") k
				chainTimes(k) = chainTimes(k) + 1
			end if
		
		end do chainfind

		write(17,*)

	end if

	! Logic to determine if a cluster gets graphed
	if ((chnsIn .ge. 3).and.(grphout .lt. 1)) then
		grphout = 1

		call chainGraph(arrin,dim1,dim2,arrin2,dimch,i,t)
					
	end if

end do clustercount

! Determine if any chains are found more than twice
chainCheck = (chainTimes .gt. 1)

! Output to screen
write(*,*) "Number of times a chain is found more than twice:", count(chainCheck)
write(18,*) "Number of times a chain is found more than twice:", count(chainCheck)

chainCheck = (chainTimes .gt. 0)

write(*,*) "Percentage of chains in a cluster:", float(count(chainCheck))/2000.0
write(18,*) "Percentage of chains in a cluster:", float(count(chainCheck))/2000.0

close(17)
close(18)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine chainGraph(arrin,dim1,dim2,chainArr,chdim,clus,t)
! Subroutine takes a set of chains and prints them for use in
! graphing. The set of chains is passed in by chainArr and 
! this subroutine finds the set of beads that belong to them
! from arrin.

use functions

implicit none
	integer,intent(in)		:: dim1, dim2, chdim ! Array dimension variables
	real,intent(in)			:: arrin(dim1,dim2) ! primary array
	integer,intent(in)		:: chainArr(chdim) ! Which chains we find
	integer,intent(in)		:: clus			   ! Holds which cluster we are finding
	real,intent(in)			:: t				! Holds the current timestep
	
	integer					:: j, k ! Looping integers
	character(30)			:: filename, chainName ! Holds filename and header info
	integer					:: beadCount, beadChain ! Counts number of beads, once it hits the number on a chain, exits 
													! interior loop. Holds number of beads on a chain.		

	real					:: OrdChain(40,3)		! holds the beads on a chain in order. 

! Determine if on chain then find where on the chain it goes. 
! Process:	chain(mol)
!			chain*40 is max value, on chain location from 40 - abs(mol - chain*40)


beadChain = 40

93 format (1f10.0,"Tcl",1i3,".dat")
94 format ('"Chain ', 1i4,'"')

write(filename,93) t, clus
	
open(unit=19,file=adjustl(trim(filename)),status="unknown",position="append")

chainLoop: do j = 1, chdim, 1

	if (chainArr(j) .eq. clus) then
		
		OrdChain = 0.0
	
		write(chainName,94) j
	
	!	write(*,*) "Chain", j, "In graph out."
	
		write(19,*) chainName

		beadCount = 0		
		
			beadLoop: do k = 1, dim1, 1 

				if (chain(arrin(k,1)) .eq. j) then
					OrdChain(nint(arrin(k,1)-(j-1)*40),1) = arrin(k,3)
					OrdChain(nint(arrin(k,1)-(j-1)*40),2) = arrin(k,4)
					OrdChain(nint(arrin(k,1)-(j-1)*40),3) = arrin(k,5)
					beadCount = beadCount + 1

				end if
			
				if (beadCount .eq. beadChain) then
					exit beadLoop
				end if
	
			end do beadLoop


!		printLoop: do k = 1, 40, 1
		
!			write(19,*) OrdChain(k,1), OrdChain(k,2), OrdChain(k,3)

!		end do printLoop

		firstfive: do k = 1, 5, 1
			write(19,*) OrdChain(k,1), OrdChain(k,2), OrdChain(k,3)
		end do firstfive
		lastfive: do k = 40, 35, -1	
			write(19,*) OrdChain(k,1), OrdChain(k,2), OrdChain(k,3)
		end do lastfive
		
		write(19,*) 
		write(19,*)

	end if
	
end do chainLoop

close(19)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine statspack(arrin1,arrin2,dim1,dim2,dimch,nClu,nMG,t)
! Subroutine tracks statistics on outputs. arrin1 is molData, arrin2 is chainData

use functions

implicit none
	integer,intent(in)	:: dim1,dim2,dimch ! Number of beads, number of attributes, number of chains
	real,intent(in)		:: arrin1(dim1,dim2) ! holds bead informations
	integer,intent(in)	:: arrin2(dimch)	! Holds chain information
	integer,intent(in)	:: nClu,nMG		! Number of clusters, number of bead groups
	real,intent(in)		:: t
	logical				:: arr1count(dim1),arr2count(dimch) ! Logical arrays to count number of beads/chains
																! in a group
	integer				:: i ! Looping integer
	integer				:: chainGrps(10), molGroups(6)  ! Holds counts for a histogram of outputs
	integer				:: counter	! Temp holding variable
	
open(unit=20,file="stats.dat",status="unknown",position="append")

chainGrps = 0
molGroups = 0

! Counts how many clusters have x-many chains
Cluscounter: do i = 1, nClu, 1

	arr2count = (arrin2 .eq. i)
	counter = count(arr2count)
	
	if ((counter .lt. 10).and.(counter .gt. 1)) then
		chainGrps(counter) = chainGrps(counter) + 1
	else if (counter .ge. 10 ) then
		chainGrps(10) = chainGrps(10) + 1
	end if
	
end do Cluscounter

! Counts how many bead groups have x-many beads
BeadCounter: do i = 1, nMG, 1
	
	arr1count = (arrin1(:,dim2) .eq. i)
	counter = count(arr1count)
	
	if ((counter .lt. 6).and.(counter .gt. 1)) then
		molGroups(counter) = molGroups(counter) + 1
	else if (counter .ge. 6) then
		molGroups(6) = molGroups(6) + 1
	end if 
	
end do BeadCounter

write(20,*) "Time", t
write(20,*) "Count of clusters with 2-10 chains"
write(20,*) "  2     3     4     5     6     7    8    9    10+"
do i = 1, 10, 1
	write(20,'(1i4," ")',ADVANCE="no") chainGrps(i)
end do
write(20,*)
write(20,*) "Count of bead-clusters with 2-6 beads"
write(20,*) "  2    3    4    5    6    7+"
do i = 1, 6, 1
	write(20,'(1i4," ")',ADVANCE="no") molGroups(i)
end do
write(20,*)
write(20,*)

close(20)

end subroutine
















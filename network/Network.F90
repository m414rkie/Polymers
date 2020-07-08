! Program takes an output file from the LAMMPS software package and
! reads in the linked beads. The program then builds a picture of
! the linked chains and determines, on average, how many chains are
! attached to an end group, given that the endgroup is connected.

! Author: Jon Parsons
! Date: 5-21-19

program clusfinder

implicit none
	character*50				:: filename  ! Name of input file
	real								:: tstep, junk ! current timestep
	integer							:: numBonds, numChains ! Number of molecules in system, Number of chains per molecule
	integer							:: numTsteps ! Number of time steps looked at
	integer							:: first	! Tracks if this is first time output subroutine is called
	integer							:: ioErr, j  ! System error variable, looping integer
	real,allocatable		:: molData(:,:)    ! molnumber, moltype, x, y, z, cluster, molgroup
	integer,allocatable :: chainBranch(:,:) ! This array is to be used to tracked how networked the system is
	 																				! dim1 is which chain and which end
																					! dim2 is the chains attached to that chain
	real								:: tavg_conn ! Time averaged number of connections

! Number of chains in the system
numChains = 2000

! Inititalizations
first = 0
tavg_conn = 0.0

! Allocation
allocate(chainBranch(2*numChains,10), stat = ioErr)

if (ioErr .ne. 0) then
	write(*,*) "Failed to allocate the chain tracking array. Exiting."
	stop
end if

! User input
100 write(*,*) "Please enter the name of the file with the data."
write(*,*) "A file will typically begin with the 'bonds' prefix."
write(*,*) "If the file is not in this directory enter the full path."
read(*,*) filename

filename = trim(filename)

open(unit=15, file=filename, status="old", action="read", iostat=ioErr)

if (ioErr .ne. 0) then
	write(*,*) "File not found, please try again."
	goto 100
end if

! Initialize number of timesteps
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

	numTsteps = numTsteps + 1

	allocate(molData(numBonds,2), stat = ioErr)

	if (ioErr .ne. 0) then
		write(*,*) "Failed to allocated primary array. Exiting at timestep", tstep
		stop
	end if

	! Initial values, to be over-ridden
	molData = 0.0
	chainBranch = 0

	! Read in molecule data
	fileread: do j = 1, numBonds, 1

		read(15,*) molData(j,1), molData(j,2), junk

	end do fileread

	call network(molData,chainBranch,numBonds,2,2*numChains,10,tavg_conn,numTsteps)
	first = 1

	deallocate(molData)

	! Checks for EOF, if not EOF then reads header data for next step
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

101 first = 2

! Average the number of connections by number of times steps and write to file
tavg_conn = tavg_conn/numTsteps

open(unit=16,file="network.dat",status="unknown",position="append")
write(16,*) "End of Time-steps"
write(16,*) "Average number of Connections (total):", tavg_conn
close(16)

write(*,*) "Number of timesteps checked:", numTsteps
write(*,*) "End of input file reached. Goodbye"

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine network(dat_in,net_arr,d_in1,d_in2,n_in1,n_in2,tavg,Tstep)

use functions
use chain_Functions

implicit none
	integer,intent(in)		:: d_in1, d_in2, n_in1, n_in2	! Dimensions of arrays
													! # bonds, 2, # endgroups, 10
	integer,intent(in)		:: Tstep ! current timestep
	real,intent(in)			 	:: dat_in(d_in1,d_in2) ! moldata and network data
	integer,intent(inout) :: net_arr(n_in1,n_in2) ! Holds which chains are attached to which endgroups
	real,intent(inout)		:: tavg

	integer		:: i, j ! Looping integers
	integer		:: chainend_a, chainend_b ! Chain endgroups of interest
	real			:: num_con_chains, num_connects ! Number of chains with clusters
	 						! and number of connections to a particular chain end
	real			:: avg_connects ! Average number of connections
	real			:: tot_connects, max_connects ! total connections, largest grouping
	real			:: max_interim ! temp counter
	real			:: sizes(n_in2)


! Loop through bonds passed in through dat_in
bond_loop: do i = 1, d_in1, 1

	! Determine chains of interest
	chainend_a = chainEnds(nint(dat_in(i,1)))
	chainend_b = chainEnds(nint(dat_in(i,2)))

	! If bonds are between beads on the same chain, skip
	if (chainend_a .eq. chainend_b) then
		cycle bond_loop
	end if

	! Connect chains (disallow same chain twice) in the net_arr
	! Loop for attaching chain end a to chain end b
	link_loop_b: do j = 1, d_in2, 1
		! If chains are already connected at that end, skip
		if (chainend_a .eq. net_arr(chainend_b,j)) then
			exit link_loop_b
		end if
		! If chains are not already connected, allow and move to the next
		if (net_arr(chainend_b,j) .eq. 0) then
			net_arr(chainend_b,j) = chainend_a
			exit link_loop_b
		end if
	end do link_loop_b
	! Loop for attaching chain end b to chain end a
	link_loop_a: do j = 1, d_in2, 1
		! If chains are already connected at that end, skip
		if (chainend_b .eq. net_arr(chainend_a,j)) then
			exit link_loop_a
		end if
		! If chains are not already connected, allow and move to the next.
		if (net_arr(chainend_a,j) .eq. 0) then
			net_arr(chainend_a,j) = chainend_b
			exit link_loop_a
		end if
	end do link_loop_a

end do bond_loop

! Count number of connections for outputs
! Initialize
max_interim = 0.0
num_con_chains = 0.0
num_connects = 0.0
max_connects = 0
check_loop: do i = 1, n_in1, 1
	! If the chain end group has no connections, skip
	if (net_arr(i,1) .eq. 0) then
			cycle check_loop
	else if (net_arr(i,1) .ne. 0) then
		num_con_chains = num_con_chains + 1.0 ! Add one to total number of ends involved
		max_interim = 0.0
		do j = 1, n_in2, 1
			if (net_arr(i,j) .ne. 0) then
				num_connects = num_connects + 1.0 ! Add one for each additional chain
																				! to the total number of connections
				max_interim = max_interim + 1.0 ! interim count of connections
			end if
		end do
		! assign +1 to bin with appropriate size
		sizes(nint(max_interim)) = sizes(nint(max_interim)) + 1
		if (max_interim .gt. max_connects) then ! update max connections
			max_connects = max_interim
		end if
	end if
end do check_loop

! Find average number of connections per chain. If a chain is not involved in a cluster,
! it is not counted towards this.
avg_connects = num_connects/num_con_chains
tavg = tavg + avg_connects

! File output
open(unit=16,file="network.dat",status="unknown",position="append")

! Output
write(16,*) "Timestep: ", Tstep
write(16,*) "Connected Chain Endgroups: ", num_con_chains
write(16,*) "Total Connections: ", num_connects
write(16,*) "Largest Group: ", max_connects
write(16,*) "Avg. Connections: ", avg_connects
write(16,*) "Time Averaged Size: ", tavg/float(Tstep)
write(16,*) "Size Breakdown"
write(16,*) "1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10"
write(16,*) sizes

close(16)

end subroutine

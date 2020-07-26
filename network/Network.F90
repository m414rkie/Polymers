! Program takes an output file from the LAMMPS software package and
! reads in the linked beads. The program then builds a picture of
! the linked chains and determines, on average, how many chains are
! attached to an end group, given that the endgroup is connected.
! Includes detection of 'looped' linkages. Includes multiple file entries.

! Author: Jon Parsons
! Date: 5-21-19

program clusfinder

implicit none
	character*50							:: raw_in  ! Name of input file
	character*50,allocatable	:: file_arr(:)
	real								:: tstep, junk ! current timestep
	integer							:: num_files ! number of input files
	integer							:: numChains ! Number of chains in system
	integer							:: time_int, bonds_read ! iteration variables
	integer							:: tot_time_steps ! Number of time steps looked at
	integer							:: numbonds ! bonds at current timestep
	integer							:: tot_bonds ! total number of bonds held in data
	integer							:: ioErr, AllErr  ! System error variables
	integer							:: i, j	! Looping integers
	real,allocatable		:: bond_data(:,:) ! Holds bond pairs at each tie
	real,allocatable		:: bond_time(:) ! Number of bonds at each timestep
	integer,allocatable :: chainBranch(:,:,:) ! This array is to be used to track
	 																				! how networked the system is
																					! dim1 is timestep
																					! dim2 is which chain and which end

! Number of chains in the system
numChains = 2000

! User input
write(*,*) "Please enter the number of files with data"
read(*,*) num_files

allocate(file_arr(num_files), stat=AllErr)
if (AllErr .ne. 0) then
	write(*,*) "File name array allocation failed"
	stop
end if

! Get file names
file_names: do i = 1, num_files, 1
	100 write(*,*) "Please enter the name of file", i
	write(*,*) "A file will typically begin with the 'bonds' prefix."
	write(*,*) "If the file is not in this directory enter the full path."
	read(*,*) raw_in

	file_arr(i) = trim(raw_in)

	open(unit=15, file=file_arr(i), status="old", action="read", iostat=ioErr)

	if (ioErr .ne. 0) then
		write(*,*) "File not found, please try again."
		goto 100
	end if

	close(15)

	file_arr(i) = trim(raw_in)

end do file_names

! get scale of data
tot_time_steps = 0
tot_bonds = 0
file_iter: do i = 1, num_files, 1

	open(unit=15, file=file_arr(i), status='old', action='read')
	write(*,*) "Getting Macro parameters from file", file_arr(i)

	! move through data on current file
	read_loop: do

		! Read data of current file
		read(15,*, END=101)
		read(15,*) tstep
		read(15,*)
		read(15,*) numBonds
		read(15,*); read(15,*); read(15,*); read(15,*);	read(15,*)

		tot_time_steps = tot_time_steps + 1
		tot_bonds = tot_bonds + numBonds

		pass_loop: do j = 1, numBonds, 1
			read(15,*) junk
		end do pass_loop

	end do read_loop

	101	close(15)

end do file_iter

! allocation statements
! goes time, bead 1, bead 2
allocate(bond_data(tot_bonds,3), stat=AllErr)
if (AllErr .ne. 0) then
	write(*,*) "Failed to allocate bond array"
	stop
end if
bond_data = 0

allocate(bond_time(tot_time_steps), stat=AllErr)
if (AllErr .ne. 0) then
	write(*,*) "Failed to allocate timestep array"
	stop
end if
bond_time = 0

allocate(chainBranch(tot_time_steps,2*numChains,5), stat = AllErr)
if (AllErr .ne. 0) then
	write(*,*) "Failed to allocate the chain tracking array. Exiting."
	stop
end if
chainBranch = 0

! Collect the data into the arrays
time_int = 0
bonds_read = 0
file_input: do i = 1, num_files, 1

	open(unit=15, file=file_arr(i), status='old', action='read')
	write(*,*) "Getting data from file", file_arr(i)

	data_input: do

		! Read data of current file
		read(15,*, END=102)
		read(15,*) tstep
		time_int = time_int + 1
		read(15,*)
		read(15,*) numBonds
		bond_time(time_int) = numBonds
		read(15,*); read(15,*); read(15,*); read(15,*);	read(15,*)

		in_loop: do j = 1, numBonds, 1
			bond_data(j+bonds_read,1) = time_int
			read(15,*)  bond_data(j+bonds_read,2), bond_data(j+bonds_read,3), junk
		end do in_loop
		bonds_read = bonds_read + numBonds

	end do data_input

	102	close(15)

end do file_input

write(*,*) "Building Network"

call network(tot_bonds,3,tot_time_steps,2*numChains,5,bond_data,bond_time, &
																																		chainBranch)

! Average the number of connections by number of times steps and write to file

write(*,*) "Number of timesteps checked:", tot_time_steps
write(*,*) "Results found in 'network.dat'"
write(*,*) "End of input file reached. Goodbye"

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine network(tot_bonds,bond_dim2,tot_time,twochains,chn_dim2,bonds, &
																													bonds_time,chain_arr)

use functions
use chain_Functions

implicit none
	integer,intent(in)		:: tot_bonds, bond_dim2, tot_time ! dims
	integer, intent(in)		:: twochains, chn_dim2 ! more dims
	real,intent(in)			 	:: bonds(tot_bonds,bond_dim2) ! holds all bonds
	real,intent(in)				:: bonds_time(tot_time) ! # of bonds at each timestep
	integer,intent(inout) :: chain_arr(tot_time,twochains,chn_dim2)
												! Holds which chains are attached to which endgroups
	integer		:: i, j, k, l ! Looping integers
	integer		:: chainend_a, chainend_b ! Chain endgroups of interest
	real			:: num_con_chains, num_connects ! Number of chains with clusters
	 						! and number of connections to a particular chain end
	real			:: avg_connects ! Average number of connections
	real			:: max_connects ! largest grouping
	real			:: max_interim ! temp counter
	real			:: tot_avg ! total average
	integer		:: bonds_counted
	real			:: sizes(chn_dim2), avg_sizes(chn_dim2) !

! initialize
tot_avg = 0.0
bonds_counted = 1
sizes = 0
avg_sizes = 0

! open output file
open(unit=16,file="network.dat",status="unknown",position="append")

time_loop: do l = 1, tot_time, 1
	! Loop through bonds passed in through dat_in
	bond_loop: do i = bonds_counted + 1, nint(bonds_counted + bonds_time(l)), 1

		! Determine chains of interest
		chainend_a = chainEnds(nint(bonds(i,2)))
		chainend_b = chainEnds(nint(bonds(i,3)))

		! If bonds are between beads on the same chain, skip
		if (chainend_a .eq. chainend_b) then
			cycle bond_loop
		end if

		! Connect chains (disallow same chain twice) in the net_arr
		! Loop for attaching chain end a to chain end b
		link_loop_b: do j = 1, twochains, 1
			! If chains are already connected at that end, skip
			if (chainend_a .eq. chain_arr(l,chainend_b,j)) then
				exit link_loop_b
			end if
			! If chains are not already connected, allow and move to the next
			if (chain_arr(l,chainend_b,j) .eq. 0) then
				chain_arr(l,chainend_b,j) = chainend_a
				exit link_loop_b
			end if
		end do link_loop_b
		! Loop for attaching chain end b to chain end a
		link_loop_a: do j = 1, chn_dim2, 1
			! If chains are already connected at that end, skip
			if (chainend_b .eq. chain_arr(l,chainend_a,j)) then
				exit link_loop_a
			end if
			! If chains are not already connected, allow and move to the next.
			if (chain_arr(l,chainend_a,j) .eq. 0) then
				chain_arr(l,chainend_a,j) = chainend_b
				exit link_loop_a
			end if
		end do link_loop_a

	end do bond_loop

	if (l .eq. 1) then
		bonds_counted = nint(bonds_time(l))
	else
		bonds_counted = bonds_counted + nint(bonds_time(l))
	end if

! Count number of connections for outputs
! Initialize
max_interim = 0.0
num_con_chains = 0.0
num_connects = 0.0
max_connects = 0
check_loop: do i = 1, twochains, 1 ! loops through the chain end groups
	! If the chain end group has no connections, skip
	if (chain_arr(l,i,1) .eq. 0) then
			cycle check_loop
	else if (chain_arr(l,i,1) .ne. 0) then
		num_con_chains = num_con_chains + 1.0 ! Add one to total number of ends involved
		max_interim = 0.0
		fit_loop:	do j = 1, chn_dim2, 1 ! add connections to th elist
			if (chain_arr(l,i,j) .ne. 0) then
				double_loop:	do k = 1, j-1, 1 ! skip over double linked chains.
					if ((chain_arr(l,i,j) .eq. chain_arr(l,i,k)) .and. (j .gt. 1)) then
						cycle fit_loop
					end if
				end do double_loop
				num_connects = num_connects + 1.0 ! Add one for each additional
					 											! chain to the total number of connections
				max_interim = max_interim + 1.0 ! interim count of connections
			else if (chain_arr(l,i,j) .eq. 0) then ! move on to the next when out of
				cycle fit_loop									! chains to attach
			end if
		end do fit_loop
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
tot_avg = tot_avg + avg_connects

! Output
write(16,*) "Timestep: ", l
write(16,*) "Connected Chain Endgroups: ", num_con_chains
write(16,*) "Total Connections: ", num_connects
write(16,*) "Largest Group: ", max_connects
write(16,*) "Avg. Connections: ", avg_connects
write(16,*) "Time Averaged Size: ", tot_avg/float(l)
write(16,*) "Size Breakdown"
write(16,*) "1 | 2 | 3 | 4 | 5 | "
write(16,*) sizes

avg_sizes = avg_sizes + sizes
sizes = 0

end do time_loop
avg_sizes = avg_sizes/tot_time
write(16,*) "Average distribution"
write(16,*) "1 | 2 | 3 | 4 | 5 | "
write(16,*) sizes

close(16)

end subroutine

! Program takes an output file from the LAMMPS software package and
! reads in the linked beads from a bonds file. The program then determines the
! number of clusters of size 1, 2, 3, 4, ...

! Author: Jon Parsons
! Date: 5-21-19

program clusfinder

use chain_Functions

implicit none
	character*50, allocatable	:: file_arr(:)  ! Holds names of input files
	character*50							:: raw_in ! Name of input file
	integer, allocatable	:: bond_time(:), bead_pairs(:,:) ! Number of bonds at
																! each timestep ; bead pairs at time (t,b1,b2)
	integer						:: num_files, tot_time_Steps ! number of input files; total
																									! timesteps found
	integer						:: tot_bonds, time_int ! total number of bonds, time tracker
	real							:: tstep, junk ! Distance parameter, number of clusters
																		! found, curent timestep
	real							:: percent
	integer						:: numBonds, numChains ! Number of molecules in system
	 																					! Number of chains per molecule
	integer						:: bonds_read ! keeps track of number of bonds read
	integer						:: lrg_track ! largest size cluster to track
	integer						:: AllErr, i, j  ! System error variable, looping integers
	integer						:: ioErr
	integer,allocatable	:: statsArr(:,:) ! Tracks number of clusters with
																	! (2,3,4...) chains involved at each timestep
	integer,allocatable	:: chainTrack(:,:), chain_ends(:,:) ! These arrays track
	                        ! chain interactions at each timestep

! Number of chains in the system, largest cluster to look for
numChains = 2000
lrg_track = 50
percent = 0.0

! get number of files to read
write(*,*) "Please enter the number of files with data"
read(*,*) num_files

allocate(file_arr(num_files), stat=AllErr)
if (AllErr .ne. 0) then
		write(*,*) "Failed to allocate file array. Exiting"
		stop
end if

! Get names of files
file_names: do i = 1, num_files, 1
	100 write(*,*) "Please enter the name of file", i
	write(*,*) "Files will typically include the 'bonds' prefix"
	write(*,*) "If the file is not in this directory enter the full path."
	read(*,*) raw_in

	file_arr(i) = trim(raw_in)

	open(unit=15, file=file_arr(i), status="old", action="read", iostat=ioErr)

	if (ioErr .ne. 0) then
		write(*,*) "File not found, please try again."
		goto 100
	end if
	close(15)

end do file_names

! initialize macro parameters
tot_time_steps = 0
tot_bonds = 0
! file stores number of bonds per timestep, for troubleshooting
open(unit=16, file='bonds_tstep.dat', status='replace', position='append')
! get total timesteps, total number of bonds
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
		write(16,*) tstep, numBonds
		read(15,*); read(15,*); read(15,*); read(15,*);	read(15,*)

		tot_time_steps = tot_time_steps + 1
		tot_bonds = tot_bonds + numBonds

		pass_loop: do j = 1, numBonds, 1
			read(15,*) junk
		end do pass_loop

	end do read_loop

	101	close(15)

end do file_iter
close(16)

write(*,*) "Total Timesteps:", tot_time_steps
write(*,*) "Total Bonds Found:", tot_bonds

! allocation statements
! which beads are bonded
allocate(bead_pairs(tot_bonds,3), stat = AllErr)
if (AllErr .ne. 0) then
	write(*,*) "Failed to allocate bonds array"
	stop
end if
bead_pairs = 0.0
! number of bonds at each timestep
allocate(bond_time(tot_time_steps), stat=AllErr)
if (AllErr .ne. 0) then
	write(*,*) "Failed to allocate time array"
	stop
end if
bond_time = 0
! Number of chains involved in a cluster at each timestep
allocate(chainTrack(tot_time_steps,numChains), stat = AllErr)
if (AllErr .ne. 0) then
	write(*,*) "Failed to allocate the chain tracking array. Exiting."
	stop
end if
chainTrack = 0
! Chain Ends
allocate(chain_ends(tot_time_steps,2*numChains), stat = AllErr)
if (AllErr .ne. 0) then
	write(*,*) "Failed to allocate the chain end tracking array"
	stop
end if
chain_ends = 0
! Number of clusters of each size
allocate(statsArr(tot_time_steps,lrg_track), stat = AllErr)
if (AllErr .ne. 0) then
	write(*,*) "Failed to allocate the statistics array"
	stop
end if
statsArr = 0.0

! get the data into the arrays
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

		if (time_int .eq. 1) then
			open(unit=23,file='pairs.dat', status='replace',position='append')
		end if

		in_loop: do j = 1, numBonds, 1
			bead_pairs(j+bonds_read,1) = time_int
			read(15,*)  bead_pairs(j+bonds_read,2), bead_pairs(j+bonds_read,3), junk
			if (time_int .eq. 1) then
				write(23,*) bead_pairs(j+bonds_read,2), bead_pairs(j+bonds_read,3)
			end if
		end do in_loop
		bonds_read = bonds_read + numBonds
		if (time_int .eq. 1) then
			close(23)
		end if
	end do data_input

	102	close(15)

end do file_input

write(*,*) "Building Clusters"
call ClusBuilder(tot_bonds,3,tot_time_steps,numChains,bead_pairs,bond_time, &
																									chainTrack,chain_ends,percent)
write(*,*) "Preparing Outputs"
call linloop(numChains,tot_time_steps,chainTrack,chain_ends)
call output(tot_time_steps,numChains,lrg_track,chainTrack,statsArr)
call statistics(tot_time_steps,lrg_track,statsArr,percent,numChains)

write(*,*) "Number of timesteps checked:", tot_time_steps

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ClusBuilder(num_bonds,dim2a1,num_tsteps,chain_amt,bond_arr, &
																					bond_time,chainIn,chain_ends,percent)
! Builds the clusters at each timestep.
! output array is chainIn.

use chain_Functions

implicit none
	integer,intent(in)		:: num_bonds, dim2a1 ! total bonds, dim2 of bonds array
	integer,intent(in)		:: num_tsteps, chain_amt ! total timesteps,
										! number of chains
	integer,intent(in)	  :: bond_arr(num_bonds,dim2a1) ! Array containing bonds
	integer,intent(in)		:: bond_time(num_tsteps) ! Array with number of bonds
																									! at each timestep
	integer,intent(inout)	:: chainIn(num_tsteps,chain_amt)
	integer,intent(inout) :: chain_ends(num_tsteps,2*chain_amt)
	real,intent(inout)		:: percent
	logical		:: chain_counter(chain_amt)
	integer		:: curClus, maxClus ! Current working cluster, number of
																			! clusters found in total
	integer		:: i, j ! Looping integers
	integer		:: bonds_counted, chn_count
	integer		:: clus_a, clus_b
	integer		:: chain_a, chain_b

! Initialize
bonds_counted = 0
chainIn = 0
chain_ends = 0
percent = 0.0

time_loop: do i = 1, num_tsteps, 1

	! Initialize
	maxClus = 1
	curClus = 1

	! Build the clusters
	ClusLoop : do j = bonds_counted + 1, bonds_counted+bond_time(i), 1

		! Ensure we are on the correct timestep
		if (bond_arr(j,1) .ne. i) then
			write(*,*) "Error in timestep tracking. Should be", i,"Is ", bond_arr(j,1)
			cycle ClusLoop
		end if

		! Set working cluster to largest cluster
		curClus = maxClus

		! get the chains involved
		chain_a = chain(float(bond_arr(j,2)))
		chain_b = chain(float(bond_arr(j,3)))
		chain_ends(i,chainEnds(bond_arr(j,2))) = 1
		chain_ends(i,chainEnds(bond_arr(j,3))) = 1

		! If both chains are not involved in a cluster, put them in new cluster
		if ((chainIn(i,chain_a) .eq. 0) .and. (chainIn(i,chain_b) .eq. 0)) then
			chainIn(i,chain_a) = curClus
			chainIn(i,chain_b) = curClus
			maxClus = maxClus + 1
			cycle ClusLoop
		end if

		! If first chain is in cluster, assign second chain to existing cluster
		if ((chainIn(i,chain_a) .ne. 0) .and. (chainIn(i,chain_b) .eq. 0)) then
			curClus = chainIn(i,chain_a)
			chainIn(i,chain_b) = curClus
			cycle ClusLoop
		end if

		! If second chain is in cluster, assign first chain to existing cluster
		if ((chainIn(i,chain_b) .ne. 0) .and. (chainIn(i,chain_a) .eq. 0)) then
			curClus = chainIn(i,chain_b)
			chainIn(i,chain_a) = curClus
			cycle ClusLoop
		end if

		! If both chains are in a cluster, assign all chains to first cluster
		if ((chainIn(i,chain_a) .ne. 0) .and. (chainIn(i,chain_b) .ne. 0) .and. &
					(chainIn(i,chain_a) .ne. chainIn(i,chain_b))) then
			clus_a = chainIn(i,chain_a)
			clus_b = chainIn(i,chain_b)
			where (chainIn(i,:) .eq. clus_b) chainIn(i,:) = clus_a
		end if

	end do ClusLoop

	if (i .eq. 1) then
		bonds_counted = bond_time(i)
	else
		bonds_counted = bonds_counted + bond_time(i)
	end if

	chain_counter = (chainIn(i,:) .ne. 0)
	chn_count = count(chain_counter)
	percent = percent + float(chn_count)/float(chain_amt)

end do time_loop

percent = percent/float(num_tsteps)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine linloop(num_chains,tot_time,chains,chain_ends)
! Determines if clusters are linear or loop on themselves
! output handled locally

implicit none
	integer,intent(in)	:: num_chains, tot_time ! describes size of system
	integer,intent(in)	:: chain_ends(tot_time,2*num_chains) ! holds which
																								! cluster a chain end belongs to
	integer,intent(in)	:: chains(tot_time,num_chains) ! holds which chains
																											! belong to which cluster
	integer							:: i, j, k ! looping integers
	character*50				:: file_out ! name of output file
	integer							:: looped, linear ! # looped an linear clusters at a time
	integer							:: tot_looped, tot_linear ! total # of looped, linear
	integer							:: tot_clusters ! local and total clusters
	integer							:: num_clus ! number of clusters at a timestep
	integer							:: unlooped_flag ! flag, if both ends not in a cluster
	integer							:: looped_flag ! flag, if both ends in a cluster
	real								:: tot_ratio, ratio ! ratio of looped to linear
	integer							:: chains_in ! number of chains in the cluster
	integer							:: loop_multiple ! loops with multiple chains
	integer							:: loop_single, loop_single_tot ! loops with a single chain

! initializations
file_out = "loops.dat"
tot_looped = 0
tot_linear = 0
tot_clusters = 0
tot_ratio = 0
loop_single_tot = 0
loop_multiple = 0
loop_single = 0

open(unit=17, file=trim(file_out), status='replace', position='append')

! iterate through time
time_loop: do i = 1, tot_time, 1

	num_clus = maxval(chains(i,:))
	looped = 0
	linear = 0
	loop_single = 0

	! iterate through clusters
	clus_loop: do j = 1, num_clus, 1
		unlooped_flag = 0
		looped_flag = 1
		chains_in = 0
		! iterate through the chains
		chain_loop: do k = 1, num_chains, 1

			! cycle if not cluster we want, else add a chain to the cluster
			if (chains(i,k) .ne. j) then
				cycle chain_loop
			else
				chains_in = chains_in + 1
			end if

			! statement to ensure nothing fishy happened
			if ((chain_ends(i,2*k) .eq. 0) .and. (chain_ends(i,2*k - 1) .eq. 0)) then
				write(*,*) "Neither endgroup in a cluster, but the chain is."
			end if

			if (chains(i,k) .ne. j) then
				write(*,*) "Chain from another cluster caught"
			end if

			! If both endgroups are not in a bond, can't be a loop
			if (chain_ends(i,2*k) .ne. chain_ends(i,2*k - 1)) then
				unlooped_flag = 1 ! not looped
			end if

		end do chain_loop

		if (chains_in .eq. 0) then
			cycle clus_loop
		end if

		if ((unlooped_flag .gt. 0).and.(chains_in .gt. 1)) then
			linear = linear + 1
		else
			if (chains_in .eq. 0) then
				write(*,*) "No chain cluster, Err"
			end if
			looped = looped + 1
			if (chains_in .gt. 1) then
				loop_multiple = loop_multiple + 1
			else if (chains_in .eq. 1) then
				loop_single = loop_single + 1
				loop_single_tot = loop_single_tot + 1
			end if
		end if

	end do clus_loop

	ratio = float(looped)/float(linear+looped)
	write(17,*) "Timestep:", i
	write(17,*) "Clusters:", looped+linear
	write(17,*) "Looped Clusters:", looped
	write(17,*) "Singly Looped Chains:", loop_single
	write(17,*) "Linear Clusters:", linear
	write(17,*) "Ratio Loop/Linear:", ratio
	tot_looped = tot_looped + looped
	loop_single_tot = loop_single_tot + loop_single
	tot_linear = tot_linear + linear
	tot_ratio = tot_ratio + ratio

end do time_loop

tot_clusters = tot_looped + tot_linear
write(17,*) "Total Timesteps:", tot_time
write(17,*) "Total Clusters:", tot_clusters
write(17,*) "Average Clusters per Timestep:", float(tot_clusters)/float(tot_time)
write(17,*) "Total Looped Clusters:", tot_looped
write(17,*) "Average Single Loops:", float(loop_single_tot)/float(tot_time)
write(17,*) "Total Linear Clusters:", tot_linear
write(17,*) "Average Ratio Loop/Linear:", tot_ratio/float(tot_time)
write(17,*) "Average Looped With Multiple Chains:", float(loop_multiple)/float(tot_time)
write(17,*) "Total Looped With Multiple Chains:", loop_multiple
write(17,*) "Ratio of Multiple Chain Loops to Total Looped:", float(loop_multiple)/float(loop_multiple+tot_linear)

close(17)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output(num_tsteps,numChains,lrg_clus,chainIn,clusArr)
! Determines size of clusters and prints to output files
! output is clusArr with a distribution of cluster sizes

implicit none
	integer,intent(in)	  :: num_tsteps ! Total timesteps we are examining
	integer,intent(in)		:: numChains ! Total number of chains in system
	integer,intent(in)		:: lrg_clus ! largest cluster possible
	integer,intent(in)	  :: chainIn(num_tsteps,numChains) ! holds which cluster
																			! chains are assigned to for each timestep
	integer,intent(inout)	:: clusArr(num_tsteps,lrg_clus) ! holds cluster sizes
																															! at each timestep
	integer				:: i, j, k ! Looping integers
	integer				:: clus_count ! temporary values for holding size of cluster
	integer				:: maxClus ! largest cluster in a timestep
	logical				:: chain_count(numChains) ! Holds number of chains involved at
																						! a given timestep
	logical				:: chain_no_clus(numChains) ! holds chains not in a cluster at
																						! a given timestep
	integer				:: in_count ! Collapsed number of chains involved
	real					:: perc_chains ! Percentage of chains
	integer				:: num_not_assoc ! number of chains not in a cluster
	integer				:: out_count

! Open files
open(unit=11,file="Clusters.dat",status="replace",position="append")
open(unit=12,file="ClusSizes.dat",status="replace",position="append")

!  write headers
write(11,*) "TimeStep	NumCLusters Percentage"

time_loop: do k = 1, num_tsteps, 1

	write(12,*) "Timestep: ", j
	write(12,*) "Cluster	Chains"

	! Find the percentage of chains in a cluster
	chain_count = (chainIn(k,:) .ne. 0)
	in_count = count(chain_count)
	perc_chains = float(in_count)/float(numChains)
	! find number of chains not in a cluster, assign to the output
	chain_no_clus = (chainIn(k,:) .eq. 0)
	num_not_assoc = count(chain_no_clus)
	!clusArr(k,1) = num_not_assoc
	out_count = 1
	! find maximum cluster
	maxClus = maxval(chainIn(k,:))
	! Iterate through the clusters
	ReadLoop: do i = 1, maxClus, 1

		! Initialize how many chains per cluster
		clus_count = 0

		! Find how many chains in current cluster
		ClusFndLoop: do j = 1, numChains, 1

			! Add one for each chain
			if (chainIn(k,j) .eq. i) then
				clus_count = clus_count + 1
			end if
			! Where are the singles coming from?


		end do ClusFndLoop

		! Write output
		if (clus_count .ne. 0) then
			write(12,*) out_count, clus_count
			out_count = out_count + 1
			! If more than lrg_clus chains in cluster, assign to max box
			if (clus_count .gt. lrg_clus) then
				clusArr(k,lrg_clus) = clusArr(k,lrg_clus) + 1
			! All else write to appropriate box
			else if ((clus_count .ge. 1).and.(clus_count .le. lrg_clus)) then
				clusArr(k,clus_count) = clusArr(k,clus_count) + 1
			end if
		end if

		! Write timestep and total number of clusters to file
		if (i .eq. maxClus) then
			write(11,*) k, out_count, perc_chains
		end if

	end do ReadLoop

end do time_loop

! close files
close(11)
close(12)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine statistics(num_tsteps,lrg_track,statsArr,percent,numChains)
! Find some statistics about clusters

implicit none
	integer,intent(in)	:: num_tsteps, lrg_track ! dimensions of working array
	integer,intent(in)	:: numChains
	integer,intent(in)	:: statsArr(num_tsteps,lrg_track) ! Holds cluster
	 																							! distributions at each timestep
	real,intent(in)		:: percent
	integer		:: i, j ! Looping integers
	integer		:: n ! Will hold total number of clusters, no singles
	integer		:: n_sngl ! total clusters with singles
	real			:: n_tot, n_tot_sngl
	real			:: sz_std_dev_tot, sz_var
	real			:: sz_std_dev_tot_sngl, sz_var_sngl
	real			:: loc_avg(num_tsteps) ! holds average cluster size, no singles
	real			:: loc_avg_sngl(num_tsteps) ! holds average cluster size, with singles
	real 		  :: num ! total chains in a cluster, no singles
	real			:: num_sngl ! total chains in a cluster, with singles
	real		  :: mean, std_dev, variance, mean_t ! no singles
	real			:: mean_sngl, std_dev_sngl, variance_sngl, mean_t_sngl ! w/ singles
	real			:: avg_sngls ! avg number of single chain clusters
	real		  :: numavg(lrg_track) ! For printing
	real			:: summed

num = 0
mean = 0
std_dev = 0
variance = 0
num_sngl = 0
mean_sngl = 0
std_dev_sngl = 0
variance_sngl = 0
numavg = 0
summed = 0.0
avg_sngls = 0.0
loc_avg = 0.0
loc_avg_sngl = 0.0
n_tot = 0.0
n_tot_sngl = 0.0
sz_std_dev_tot = 0.0
sz_std_dev_tot_sngl = 0.0
sz_var = 0.0
sz_var_sngl = 0.0

! Value updates
time_loop: do j = 1, num_tsteps, 1

	n = 0
	n_sngl = 0
	num = 0.0
	num_sngl = 0.0
	sz_var = 0.0
	sz_var_sngl = 0.0

	! Global stuff
	! x
	numavg(1) = numavg(1) + statsArr(j,1)
	n_sngl = n_sngl + statsArr(j,1)
	num_sngl = statsArr(j,1)
	do i = 2, lrg_track, 1
		num = num + float((i)*statsArr(j,i))
		num_sngl = num_sngl + float(i)*statsArr(j,i)
		numavg(i) = numavg(i) + statsArr(j,i)
		n = n + statsArr(j,i) ! get total number of clusters
		n_sngl = n_sngl + statsArr(j,i) ! total number of clusters w/ singles
	end do

	n_tot = n_tot + float(n)
	n_tot_sngl = n_tot_sngl + float(n_sngl)

	! no singles
	! mean number of chains per cluster, save to array
	mean = num/float(n)
	loc_avg(j) = mean

	! with singles
	! mean number of chains per cluster, save to array
	mean_sngl = num_sngl/float(n_sngl)
	loc_avg_sngl(j) = mean_sngl

	! std deviation of cluster size
	sz_var_sngl = float(statsArr(j,1)) * (1.0 - mean_sngl)**2
	do i = 2, lrg_track, 1
		sz_var = sz_var + float(statsArr(j,i)) * (float(i)-mean)**2
		sz_var_sngl = sz_var_sngl + float(statsArr(j,i)) * (float(i)-mean_sngl)**2
	end do

	sz_var = sz_var/num
	sz_var_sngl = sz_var_sngl/num_sngl

	sz_std_dev_tot = sz_std_dev_tot + sqrt(sz_var)
	sz_std_dev_tot_sngl = sz_std_dev_tot_sngl + sqrt(sz_var_sngl)

end do time_loop

! find std dev.
mean_t = sum(loc_avg)/float(num_tsteps)
mean_t_sngl = sum(loc_avg_sngl)/float(num_tsteps)

do j = 1, num_tsteps, 1
	! no singles
	variance = variance + (loc_avg(j) - mean_t)**2
	! singles
	variance_sngl = variance_sngl + (loc_avg_sngl(j) - mean_t_sngl)**2
end do
variance = variance/float(num_tsteps)
std_dev = sqrt(variance)

variance_sngl = variance_sngl/float(num_tsteps)
std_dev_sngl = sqrt(variance_sngl)

! Final output
write(*,*) "Statistics Output in file 'Averages.dat'"
open(unit=13,file="Averages.dat",status="replace",position="append")
write(13,*) "Box Plot section involves only chains involved in a bond."

write(13,*) "Total Time Steps: ", num_tsteps
write(13,*) "Average Cluster Size: ", mean_t
write(13,*) "Average Cluster Size with Singles: ", mean_t_sngl
write(13,*) "Std devation of Cluster Size (no singles): ", sz_std_dev_tot/float(num_tsteps)
write(13,*) "Std devation of Cluster Size (w/ singles): ", sz_std_dev_tot_sngl/float(num_tsteps)
write(13,*) "Std devation of Mean Over Time (no singles): ", std_dev
write(13,*) "Std devation of Mean Over Time (w/ singles): ", std_dev_sngl
write(13,*) "Avg. Number of Clusters (no singles): ", (n_tot)/float(num_tsteps)
write(13,*) "Avg. Number of Clusters (w/ singles): ", (n_tot_sngl)/float(num_tsteps)
write(13,*) "Avg. Percentage of Unbound Chains: ", (1.0-percent)
! Output for box plot
write(13,*) "Box Plot"
do i = 1, lrg_track, 1
	write(13,*) i, numavg(i)/float(num_tsteps)
	summed = summed + float(i)*numavg(i)/float(num_tsteps)
end do
close(13)
write(*,*) "Averages sum:",summed
end subroutine

! Program takes an output file from the LAMMPS software package and
! reads in the location and type of each molecule. From this, this
! program attempts to determine which beads of type 5 (H) are bonded in a
! pair-wise fashion.

! This version does not track chains, only beads

! Author: Jon Parsons
! Date: 2-1-19

program clusfinder

implicit none
	character*50, allocatable	:: file_arr(:)
	character*50							:: raw_in
	real				  :: sigma, tstep
	real					:: boxDim(3,2)
	real					:: xd, yd, zd, d_lim
	integer				:: numMols
	integer				:: num_files
	integer				:: bond_mx ! maximum number of bonds per micelle
	integer				:: tot_time_steps, num_step
	integer				:: ioErr, AllErr
	integer				:: j, k
	real,allocatable		:: molData(:,:,:) ! time,numMols,(bead,moltype,x,y,z)
	real,allocatable		:: bonds(:,:) ! time, numMols holds associations
	integer,allocatable	:: stats(:,:) ! time, numMols holds associations


! Number of beads in the system
numMols = 80002
! maximum size of a micelle we are considering
bond_mx = 100

! User input and initializations
write(*,*) "Enter sigma:"
read(*,*) sigma

write(*,*) "Please enter the number of files."
read(*,*) num_files

allocate(file_arr(num_files), stat=AllErr)
if (AllErr .ne. 0) then
	write(*,*) "Allocation of file name array has failed. Exiting."
	stop
end if
file_in: do k = 1, num_files, 1
	100	write(*,*) "Typical files will begin with the 'pict' prefix."
	write(*,*) "If the file is not in this directory enter the full path."
	write(*,*) "Enter the name of file", k
	read(*,*) raw_in

	open(unit=15, file=trim(raw_in), status="old", action="read", iostat=ioErr)

	if (ioErr .ne. 0) then
		write(*,*) "File not found, please try again."
		goto 100
	end if

	close(15)
	file_arr(k) = trim(raw_in)

end do file_in

! get macro parameters
tot_time_steps = 0
write(*,*) "Obtaining macro parameters"
! Reads in the macro data
file_macro_iter: do k = 1, num_files, 1
	open(unit=15, file=file_arr(k), status="old", action="read", iostat=ioErr)
	if (ioErr .ne. 0) then
		write(*,*) "Failed to open file #", k, "Exiting"
		stop
	end if
	read_loop: do

		read(15,*,END=101)
		read(15,*) tstep
		read(15,*)
		read(15,*) numMols
		read(15,*)
		read(15,*) boxDim(1,1), boxDim(1,2)
		read(15,*) boxDim(2,1), boxDim(2,2)
		read(15,*) boxDim(3,1), boxDim(3,2)
	  read(15,*)
		tot_time_steps = tot_time_steps + 1


		pass_loop: do j = 1, numMols, 1
			read(15,*)
		end do pass_loop

	end do read_loop
	101 close(15)
end do file_macro_iter

xd = boxDim(1,2) - boxDim(1,1)
yd = boxDim(2,2) - boxDim(2,1)
zd = boxDim(3,2) - boxDim(3,1)

! Assumes a cubic box
d_lim = 0.5*xd

write(*,*) "Number of timesteps found:", tot_time_steps
write(*,*) "Getting data"

! Allocation statements
allocate(molData(tot_time_steps,numMols,5), stat = AllErr)
if (AllErr .ne. 0) then
	write(*,*) "Failed to allocate primary array. Exiting"
	write(*,*) "Error status", AllErr
	stop
end if
molData = 0
allocate(bonds(tot_time_steps,numMols), stat = AllErr)
if (AllErr .ne. 0) then
	write(*,*) "Failed to allocate bonds array. Exiting"
	write(*,*) "Error status", AllErr
	stop
end if
bonds = 0
allocate(stats(tot_time_steps,bond_mx), stat = AllErr)
if (AllErr .ne. 0) then
	write(*,*) "Failed to allocate stats array. Exiting"
	write(*,*) "Error status", AllErr
	stop
end if
stats = 0

num_step = 0
! collect the data
file_data_in: do k = 1, num_files, 1
	open(unit=15, file=file_arr(k), status="old", action="read", iostat=ioErr)
	if (ioErr .ne. 0) then
		write(*,*) "Failed to open file #", k, "Exiting"
		stop
	end if
	data_collect: do

		read(15,*,END=102)
		read(15,*) tstep
		read(15,*)
		read(15,*) numMols
		read(15,*); read(15,*); read(15,*); read(15,*); read(15,*);
		num_step = num_step + 1

		input_loop: do j = 1, numMols, 1
			read(15,*) molData(num_step,j,1), molData(num_step,j,2) , &
				molData(num_step,j,3), molData(num_step,j,4), molData(num_step,j,5)
		end do input_loop

	end do data_collect
	102 close(15)
end do file_data_in

write(*,*) "Determining bonds"
call clusSort(tot_time_steps,numMols,5,molData,bonds,sigma,d_lim)
write(*,*) "Sorting output"
call output(tot_time_steps,numMols,bonds,stats)
write(*,*) "Finding Distribution"
call statistics(tot_time_steps,bond_mx,stats)

write(*,*) "End of data reached. Goodbye"

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine clusSort(tsteps,numMols,vals,datin,bonds_out,dParam,d_lim)
! Subroutine to determine which molecules of the appropriate type are
! joined together in a bond.

implicit none
	integer,intent(in)	:: tsteps,numMols,vals ! Dimensions of input array
	real,intent(in) 		:: datin(tsteps,numMols,vals) ! Input array
	real,intent(inout) 	:: bonds_out(tsteps,numMols) ! output array
	real,intent(in)			:: dParam, d_lim  ! Distance parameter

	integer			:: i, j, k, l ! Looping integers
	integer			:: b1, b2, b3 ! holding variables for bead number
	integer			:: up_frm, up_to ! update parameters
	integer			:: b_type = 5 ! type of bead we care about
	integer			:: bonds, bond_count ! number of bonds found
	integer			:: old_bond ! holds bond number we are changing if needed.
	integer			:: cur_bond, flg ! bond number being assigned, flg if bond exists
	real				:: x1,x2, y1,y2, z1,z2 ! xyz positions
	real				:: dx, dy, dz ! axial distances
	real				:: near ! distance between two beads
	real				:: r, a, b, c ! function variables


! simple 3d distance formula
r(a,b,c) = sqrt((a)**2 + (b)**2 + (c)**2)

! Nested do loops iterate through the data. When conditions are met assigns the
! current molecules to a cluster
time_Loop: do i = 1, tsteps, 1

	! Outputs to show that it is still working
	if (i .eq. floor(0.1*real(tsteps))) then
		write(*,*) "Clustering 10% complete"
	else if (i .eq. floor(0.2*real(tsteps))) then
		write(*,*) "Clustering 20% complete"
	else if (i .eq. floor(0.2*real(tsteps))) then
		write(*,*) "Clustering 20% complete"
	else if (i .eq. floor(0.3*real(tsteps))) then
		write(*,*) "Clustering 30% complete"
	else if (i .eq. floor(0.4*real(tsteps))) then
		write(*,*) "Clustering 40% complete"
	else if (i .eq. floor(0.5*real(tsteps))) then
		write(*,*) "Clustering 50% complete"
	else if (i .eq. floor(0.6*real(tsteps))) then
		write(*,*) "Clustering 60% complete"
	else if (i .eq. floor(0.7*real(tsteps))) then
		write(*,*) "Clustering 70% complete"
	else if (i .eq. floor(0.8*real(tsteps))) then
		write(*,*) "Clustering 80% complete"
	else if (i .eq. floor(0.9*real(tsteps))) then
		write(*,*) "Clustering 90% complete"
	end if

	bonds = 0
	bond_count = 1
	bond_loop: do j = 1, numMols, 1
		bonds = 0
		flg = 0
		b1 = nint(datin(i,j,1))
		! cycle if not the correct bead type
		if (nint(datin(i,j,2)) .ne. b_type) then
			cycle bond_loop
		end if
		! assign position variables for first bead
		x1 = datin(i,j,3)
		y1 = datin(i,j,4)
		z1 = datin(i,j,5)

		! if already in a bond group get that bond group
		if (bonds_out(i,b1) .ne. 0) then
			cur_bond = int(bonds_out(i,b1))
			flg = 1
		else
			flg = 0
		end if

		up_frm = 1 ! reset update loop counter

		attachment_loop: do k = 1, numMols, 1

			b2 = nint(datin(i,k,1))

			! cycle if same bead or incorrect type
			if ((b1 .eq. b2) .or. (nint(datin(i,k,2)) .ne. b_type)) then
				cycle attachment_loop
			end if

			! assign position variables for second bead
			x2 = datin(i,k,3)
			y2 = datin(i,k,4)
			z2 = datin(i,k,5)
			! Get distances
			dx = abs(x2-x1)
			dy = abs(y2-y1)
			dz = abs(z2-z1)
			! Handle periodic boundary conditions
			if (dx .gt. d_lim) then
				dx = dx - d_lim
			end if
			if (dx .gt. d_lim) then
				dy = dy - d_lim
			end if
			if (dz .gt. d_lim) then
				dz = dz - d_lim
			end if

			near = r(dx,dy,dz)

			if (near .gt. dParam) then
				cycle attachment_loop
			end if

			if ((bonds_out(i,b2) .ne. 0) .and. (flg .eq. 0)) then
			! second bead has a bond and first does not
				cur_bond = int(bonds_out(i,b2))
			else if ((bonds_out(i,b2) .ne. 0) .and. (flg .eq. 1)) then
			! both beads have a bond, set all beads w/ second bond to first bond
				old_bond = int(bonds_out(i,b2))
				up_to = k
				update_loop: do l = up_frm, up_to, 1
					b3 = nint(datin(i,l,1))
					if (bonds_out(i,b3) .eq. old_bond) then
						bonds_out(i,b3) = cur_bond
					end if
				end do update_loop
				up_frm = up_to
			else if ((bonds_out(i,b2) .eq. 0) .and. (flg .eq. 0)) then
			! neither has a bond.
				cur_bond = bond_count
			end if

			! final assignment to a micelle
			bonds_out(i,b1) = cur_bond
			bonds_out(i,b2) = cur_bond

			bonds = bonds + 1

		end do attachment_loop

		if (bonds .gt. 0) then
			bond_count = bond_count + 1
		end if

	end do bond_loop

end do time_loop

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output(tsteps,numMols,bonds,micelle_bins)
! subroutine to output the results

implicit none
	integer,intent(in)		:: tsteps, numMols
	real,intent(in) 			:: bonds(tsteps,numMols)
	integer,intent(inout)	:: micelle_bins(tsteps,100) ! Holds micelle sizes


	character*25				:: out_file , stats_file! filename for output
	integer							:: i, j, k ! looping integers
	integer							:: max_bonds ! number of bonds in the timestep
	integer							:: count ! number of beads in the cluster

out_file = "h_bonds.dat"
stats_file = "micelle_stats.dat"

open(unit=19,file=trim(out_file),status="replace",position="append")

! iterate through time
time_loop: do i = 1, tsteps, 1

	if (i .eq. floor(0.1*real(tsteps))) then
		write(*,*) "Output 10% complete"
	else if (i .eq. floor(0.2*real(tsteps))) then
		write(*,*) "Output 20% complete"
	else if (i .eq. floor(0.2*real(tsteps))) then
		write(*,*) "Output 20% complete"
	else if (i .eq. floor(0.3*real(tsteps))) then
		write(*,*) "Output 30% complete"
	else if (i .eq. floor(0.4*real(tsteps))) then
		write(*,*) "Output 40% complete"
	else if (i .eq. floor(0.5*real(tsteps))) then
		write(*,*) "Output 50% complete"
	else if (i .eq. floor(0.6*real(tsteps))) then
		write(*,*) "Output 60% complete"
	else if (i .eq. floor(0.7*real(tsteps))) then
		write(*,*) "Output 70% complete"
	else if (i .eq. floor(0.8*real(tsteps))) then
		write(*,*) "Output 80% complete"
	else if (i .eq. floor(0.9*real(tsteps))) then
		write(*,*) "Output 90% complete"
	end if

	write(19,*)
	write(19,*) "Timestep: ", i
	max_bonds = nint(maxval(bonds(i,:)))
	write(19,*) "Bonds: ", max_bonds
	write(19,*) "Beads ... Num"
	! iterate through the clusters
	cluster_loop: do j = 1, max_bonds, 1
		count = 0

		bead_loop: do k = 1, numMols, 1
			if (bonds(i,k) .eq. j) then
				write(19,'(1i6, " ")',ADVANCE='no') k
				count = count + 1
			end if
		end do bead_loop

		if (count .eq. 0) then
			cycle cluster_loop
		end if

		if (count .gt. 0) then
			write(19,*) count
		end if

		if (count .ge. 100) then
			micelle_bins(i,100) =	micelle_bins(i,100) + 1
		else
			micelle_bins(i,count) =	micelle_bins(i,count) + 1
		end if

	end do cluster_loop

end do time_loop

close(19)

write(*,*) "Micelle beads found in ", out_file

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine statistics(num_tsteps,lrg_track,statsArr)
! Find some statistics about clusters

implicit none
	integer,intent(in)	:: num_tsteps, lrg_track ! dimensions of working array
	integer,intent(in)	:: statsArr(num_tsteps,lrg_track) ! Holds cluster
	 																							! distributions at each timestep
	integer		:: i, j ! Looping integers
	integer		:: n ! Will hold total number of clusters, no singles
	integer		:: n_sngl ! total clusters with singles
	real			:: n_tot
	real			:: sz_std_dev_tot, sz_var
	real			:: loc_avg(num_tsteps) ! holds average cluster size, no singles
	real 		  :: num ! total chains in a cluster, no singles
	real		  :: mean, std_dev, variance, mean_t ! no singles
	real		  :: numavg(lrg_track) ! For printing
	real			:: summed

num = 0
mean = 0
std_dev = 0
variance = 0
numavg = 0
summed = 0.0
loc_avg = 0.0
n_tot = 0.0
sz_std_dev_tot = 0.0
sz_var = 0.0

! Value updates
time_loop: do j = 1, num_tsteps, 1

	n = 0
	n_sngl = 0
	num = 0.0
	sz_var = 0.0

	! Global stuff
	! x
	numavg(1) = numavg(1) + statsArr(j,1)
	do i = 2, lrg_track, 1
		num = num + float((i)*statsArr(j,i))
		numavg(i) = numavg(i) + statsArr(j,i)
		n = n + statsArr(j,i) ! get total number of clusters
	end do

	n_tot = n_tot + float(n)

	! no singles
	! mean number of chains per cluster, save to array
	mean = num/float(n)
	loc_avg(j) = mean

	! std deviation of cluster size
	do i = 2, lrg_track, 1
		sz_var = sz_var + float(statsArr(j,i)) * (float(i)-mean)**2
	end do

	sz_var = sz_var/num

	sz_std_dev_tot = sz_std_dev_tot + sqrt(sz_var)

end do time_loop

! find std dev.
mean_t = sum(loc_avg)/float(num_tsteps)

do j = 1, num_tsteps, 1
	variance = variance + (loc_avg(j) - mean_t)**2
end do
variance = variance/float(num_tsteps)
std_dev = sqrt(variance)

! Final output
write(*,*) "Statistics Output in file 'hh_Averages.dat'"
open(unit=13,file="hh_Averages.dat",status="replace",position="append")
write(13,*) "Box Plot section involves only chains involved in a bond."
write(13,*) "Total Time Steps: ", num_tsteps
write(13,*) "Average Micelle Size: ", mean_t
write(13,*) "Std devation of Cluster Size: ", sz_std_dev_tot/float(num_tsteps)
write(13,*) "Std devation of Mean Over Time: ", std_dev
write(13,*) "Avg. Number of Clusters: ", (n_tot)/float(num_tsteps)
! Output for box plot
write(13,*) "Box Plot"
do i = 1, lrg_track, 1
	write(13,*) i, numavg(i)/float(num_tsteps)
	summed = summed + float(i)*numavg(i)/float(num_tsteps)
end do
close(13)
write(*,*) "Averages sum:",summed
end subroutine

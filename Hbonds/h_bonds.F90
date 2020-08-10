! Program takes an output file from the LAMMPS software package and
! reads in the location and type of each molecule. From this, this
! program attempts to determine which beads of type 2 (H) are bonded in a
! pair-wise fashion.

! This version does not track chains, only beads

! Author: Jon Parsons
! Date: 2-1-19

program clusfinder

implicit none
	character*50	:: filename
	real				  :: sigma, tstep
	integer				:: numMols
	integer				:: tot_time_steps, num_step
	integer				:: ioErr, AllErr
	integer				:: j, first
	real,allocatable	:: molData(:,:,:) ! time,numMols,(bead,moltype,x,y,z)
	real,allocatable	:: bonds(:,:) ! time, numMols holds associations


! Number of beads in the system
numMols = 80000

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

! get macro parameters
tot_time_steps = 0
write(*,*) "Obtaining macro parameters"
! Reads in the data
read_loop: do

	read(15,*,END=101)
	read(15,*) tstep
	read(15,*)
	read(15,*) numMols
	read(15,*); read(15,*); read(15,*); read(15,*); read(15,*);
	tot_time_steps = tot_time_steps + 1


	pass_loop: do j = 1, numMols, 1
		read(15,*)
	end do pass_loop

end do read_loop
101 close(15)

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

! collect the data
open(unit=15, file=filename, status="old", action="read", iostat=ioErr)
num_step = 0
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

write(*,*) "Determining bonds"
call clusSort(tot_time_steps,numMols,5,molData,bonds,sigma)
write(*,*) "Sorting output"
call output(tot_time_steps,numMols,bonds)

write(*,*) "End of data reached. Goodbye"

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine clusSort(tsteps,numMols,vals,datin,bonds_out,dParam)
! Subroutine to determine which molecules of the appropriate type are
! joined together in a bond.

implicit none
	integer,intent(in)	:: tsteps,numMols,vals ! Dimensions of input array
	real,intent(in) 		:: datin(tsteps,numMols,vals) ! Input array
	real,intent(inout) 	:: bonds_out(tsteps,numMols) ! output array
	real,intent(in)			:: dParam  ! Distance parameter

	integer					:: i, j, k ! Looping integers
	integer					:: bonds, bond_count ! number of bonds found
	integer					:: cur_bond, flg ! bond number being assigned
	real						:: x1,x2, y1,y2, z1,z2 ! xyz positions
	real						:: near ! distance between two beads
	real						:: r, a,b, c,d, e,f ! function variables

! simple 3d distance formula
r(a,b,c,d,e,f) = sqrt((b-a)**2 + (d-c)**2 + (f-e)**2)

! Nested do loops iterate through the data. When conditions are met assigns the
! current molecules to a cluster
time_Loop: do i = 1, tsteps, 1

	! Outputs to show that it is still working, will suppress if found to be
	! excessively expensive
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
		! cycle if not the correct bead type
		if (nint(datin(i,j,2)) .ne. 2) then
			cycle bond_loop
		end if
		! assign position variables for first bead
		x1 = datin(i,j,3)
		y1 = datin(i,j,4)
		z1 = datin(i,j,5)

		! if already in a bond group get that bond group
		if (bonds_out(i,nint(datin(i,j,1))) .ne. 0) then
			cur_bond = int(bonds_out(i,nint(datin(i,j,1))))
			flg = 1
		end if

		attachment_loop: do k = 1, numMols, 1
			! cycle if same bead or incorrect type
			if ((nint(datin(i,k,1)).eq.nint(datin(i,j,1))).or. &
																							(nint(datin(i,k,2)) .ne. 2)) then
				cycle attachment_loop
			end if

			if ((bonds_out(i,nint(datin(i,k,1))) .ne. 0) .and. (flg .eq. 0)) then
				cur_bond = int(bonds_out(i,nint(datin(i,k,1))))
			else
				cur_bond = bond_count
			end if

			! assign position variables for second bead
			x2 = datin(i,k,3)
			y2 = datin(i,k,4)
			z2 = datin(i,k,5)

			near = r(x1,x2,y1,y2,z1,z2)

			if (near .le. dParam) then
				bonds = bonds + 1
				bonds_out(i,nint(datin(i,j,1))) = cur_bond
				bonds_out(i,nint(datin(i,k,1))) = cur_bond
			end if

		end do attachment_loop

		if (bonds .gt. 0) then
			bond_count = bond_count + 1
		end if

	end do bond_loop


end do time_loop

close(19)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output(tsteps,numMols,bonds)
! subroutine to output the results

implicit none
	integer,intent(in)	:: tsteps, numMols
	real,intent(in) 		:: bonds(tsteps,numMols)

	character*25				:: out_file ! filename for output
	integer							:: i, j, k ! looping integers
	integer							:: max_bonds ! number of bonds in the timestep
	integer							:: count ! number of beads in the cluster

out_file = "h_bonds.dat"

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
		if (count .gt. 0) then
			write(19,*) count
		end if

	end do cluster_loop

end do time_loop

close(19)

write(*,*) "Outputs found in ", out_file

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

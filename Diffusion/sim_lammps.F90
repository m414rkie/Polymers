! Program to output a simulated LAMMPS output file.
! In particular the file will move a particle in a known path for use in
! testing diffusion-finding codes.

! Author: Jon Parsons
! Date: 9-21-19

program lammp_test

implicit none
	character*50			:: file_out, file_diag, file_msd, file_exp

	real							:: f, r, m, d, t, a
	real							:: x, y, z
	real							:: junk
	integer						:: num_beads
	real,allocatable	:: bead_pos_li(:,:,:), bead_pos_di(:,:,:)
	real,allocatable	:: bead_pos_ex(:,:,:)
	real							:: xmax, xmin, ymin, ymax, zmax, zmin
	real							:: dt, disp, ran, walk, ran2
	integer						:: t_step
	integer						:: numsteps
	integer						:: eqn_flag

	integer						:: i, j

xmax = 10.0
xmin = -10.0
ymax = 10.0
ymin = -10.0
zmax = 10.0
zmin = -10.0

junk = 100.0
t_step = 1

x = 0.0
y = 0.0
z = 0.0

num_beads = 12000

write(*,*) "Number of timesteps to simulate:"
read(*,*) numsteps
write(*,*) "Type of Equation:"

file_msd = "sim_rw_3d.dat"

allocate(bead_pos_li(num_beads,3,numsteps))
bead_pos_li = 0.0
allocate(bead_pos_di(num_beads,3,numsteps))
bead_pos_di = 0.0
allocate(bead_pos_ex(num_beads,3,numsteps))
bead_pos_di = 0.0


! Simulation part - random walk
! 3D

x = 0.0
y = 0.0
z = 0.0

dt = 1.0/float(numsteps)

disp = 1.0

time_m: do i = 2, numsteps, 1

	bead_m: do j = 1, num_beads, 1

			call random_number(ran)

			if (ran .lt. 1.0/2.0) then
				walk = disp
			else
				walk = -disp
			end if

			x = bead_pos_ex(j,1,i-1)
			x = x + walk
			if (x .gt. xmax) then
				x = xmin + disp
			end if
			if (x .lt. xmin) then
				x = xmax - disp
			end if
			bead_pos_ex(j,1,i) = x

			call random_number(ran)

			if (ran .lt. 1.0/2.0) then
				walk = disp
			else
				walk = -disp
			end if

			y = bead_pos_ex(j,2,i-1)
			y = y + walk
			if (y .gt. ymax) then
				y = ymin + disp
			end if
			if (y .lt. ymin) then
				y = ymax - disp
			end if
			bead_pos_ex(j,2,i) = y

			call random_number(ran)

			if (ran .lt. 1.0/2.0) then
				walk = disp
			else
				walk = -disp
			end if

			z = bead_pos_ex(j,3,i-1)
			z = z + walk
			if (z .gt. zmax) then
				z = zmin + disp
			end if
			if (z .lt. zmin) then
				z = zmax - disp
			end if
			bead_pos_ex(j,3,i) = z

	end do bead_m

end do time_m

! output
open(unit=17,file=file_msd,status="replace",position="append")

time: do i = 1, numsteps, 1

	! Header things
	write(17,*) junk
	write(17,*) t_step
	write(17,*) junk
	write(17,*) num_beads
	write(17,*) junk
	write(17,*) xmin, xmax
	write(17,*) ymin, ymax
	write(17,*) zmin, zmax
	write(17,*) junk

	bead: do j = 1, num_beads, 1
			write(17,*) j, 1, bead_pos_ex(j,1,i), bead_pos_ex(j,2,i), bead_pos_ex(j,3,i)
	end do bead

	t_step = t_step + 1

end do time

close(17)

end program

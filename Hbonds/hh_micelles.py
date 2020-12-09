#! usr/bin/python3

# A code to determine the average size of the hydrophilic associations of a
# mucin network. Input is a pict file from the LAMMPS software suite.

# Jon Parsons
# 11-4-2020

import matplotlib.pyplot as plt
import numpy as np
import os.path

################################################################################
###################### Classes #################################################
################################################################################
class bead:
# holds information for beads
    def __init__(self,number,type,x,y,z,micelle):
        self.number = number # number assigned to bead by lammps
        self.type = type # bead sub-type
        self.x = x # x coordinate
        self.y = y # y coordinate
        self.z = z # z coordinage
        self.micelle = micelle # if in a micelle, holds micelle value, else 0

################################################################################
class micelle:
# Holds beads which are associated in a particular micelle
    def __init__(self,beads):
         self.beads = []

################################################################################
###################### Functions ###############################################
################################################################################
# Function to get data out of a pict file. Assigns data of appropriate type to
# an array using the bead object. File structure is known, function takes
# advantage of this.
def get_data(f_name):
# INPUTS
# f_name : name of file data is in - char
# RETURNS
# data_list : 2d list. d1 is timesteps, d2 is bead data including coord and type
# t_steps : number of timesteps in the file - int

    # Initializations
    t_steps = 0; # timesteps
    data_list = []
    line_count = 0
    bd_count = 0

    print("\nOpening file {}".format(f_name))
    with open(f_name,'r') as d:
        f_dat = d.read().splitlines()

    for line in f_dat:
        line_count += 1
        d_ln = line.split()
        if line_count == 2: # line contains timestep information
            t_steps = t_steps + 1
        if line_count == 10:
            beads = [] # initialize bead list
        if line_count >= 10:
            if int(d_ln[1]) == 3: # only need information from this type
                nw_bead = bead(
                    int(d_ln[0]),
                    int(d_ln[1]),
                    float(d_ln[2]),
                    float(d_ln[3]),
                    float(d_ln[4]),
                    int(0))
                beads.append(nw_bead)
                bd_count += 1
        if line_count == 80011: # known
            line_count = 0
            data_list.append(beads)

    # inform user of number of timesteps found.
    print("\nTotal number of Timesteps: {}\n".format(t_steps))
    print("Total number of Beads Found: {}\n".format(bd_count))
    return data_list, t_steps

################################################################################
# Function that gives the distance between two data points in cartesian space
def distance(x1,y1,z1,x2,y2,z2):
# INPUTS
# all values are floats containing the cartesian coordinates of a bead in 3d
# RETURNS
# d : distance between the two beads - float

    d = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

    return d

################################################################################
# Function that updates bead micelles when combining micelles together
def combine(bead_list,keep,lose):
# INPUTS
# bead_list : list of the bead data, single timestep
# keep : micelle value to replace the lose value with
# lose : micelle value that is being replaced
# RETURNS
# updated bead_list

    for i in range(len(bead_list)):
        if bead_list[i].micelle == lose:
            bead_list[i].micelle = keep
        if bead_list[i].micelle > lose:
            bead_list[i].micelle -= 1

    return  bead_list

################################################################################
# Function handles assigning beads to micelles and when needed combines micells
def micelle_assign(b1,b2,bead_list,micelle_list):
# INPUTS
# b1 : number of bead 1 - int
# b2 : number of bead 2 - int
# bead_list : bead data, only single timestep
# micelle_lst : list of micelle data at current timestep
# RETURNS
# returns modified micelle list, bead_list

    m1 = bead_list[b1].micelle # get current micelle associations
    m2 = bead_list[b2].micelle

    if m1 == 0 and m2 == 0:
    # Neither bead is in a micelle yet
    # make a new one
        m_cnt = len(micelle_list) # new micelle count
        nw_micelle = micelle([])
        nw_micelle.beads.append(bead_list[b1].number) # assign beads to
        nw_micelle.beads.append(bead_list[b2].number) # new micelle
        bead_list[b1].micelle = m_cnt # update bead micelle associations
        bead_list[b2].micelle = m_cnt
        micelle_list.append(nw_micelle)
        return micelle_list, bead_list

    if m1 == 0 and m2 != 0:
    # Second bead is in a micelle
        # assign to old one
        bead_list[b1].micelle = m2
        micelle_list[m2-1].beads.append(bead_list[b1].number)
        return micelle_list, bead_list

    if m1 != 0 and m2 == 0:
    # First bead is in a micelle
        # assign b2 to old one
        bead_list[b2].micelle = m1
        micelle_list[m1-1].beads.append(bead_list[b2].number)
        return micelle_list, bead_list

    if m1 != 0 and m2 != 0 and m1 != m2:
    # Both already have an assigned micelle, which is not the same micelle
        # combine the two, keepig micelle from m1
        m_keep = m1
        m_lose = m2
        micelle_list[m_keep-1].beads.extend(
            micelle_list[m_lose-1].beads
        )
        del micelle_list[m_lose-1]
        bead_list = combine(bead_list,m_keep,m_lose)
        return micelle_list, bead_list

    return micelle_list, bead_list

################################################################################
# Function that creates micelles based on distances between beads
def micelle_maker(bead_list,timesteps,rng):
# INPUTS
# bead_list : 2d list of bead data, 1d is time.
# timesteps : number of timesteps we are checking - int
# rng - distance at which we consider the beads bonded - float
# RETURNS
# micelles_time : list of micelle data - micelles
# data_list : modified bead data

    # initializations
    micelles_time = []
    # Number of beads we are interested in determined in the input subroutine,
    # get that number for use here.
    # Should be same number of beads at every timestep barring a disaster.
    bead_num = len(bead_list[0][:])
    # Loop through time
    for t in range(timesteps):
        micelles = []
        # Loop through beads
        # b1 loop is bead we are comparing distance to
        for b1 in range(bead_num):
            # b2 loop is beads we are a comparing to b1
            for b2 in range(b1,bead_num):
                if b1 != b2:
                    bd1 = bead_list[t][b1] # obtain and find distances between
                    bd2 = bead_list[t][b2] # beads of interest
                    dist = distance(bd1.x,bd1.y,bd1.z,bd2.x,bd2.y,bd2.z)
                    if dist < rng: # if below criteria, determine micelle
                        micelles, bead_list[t][:] = micelle_assign(
                                                        b1,
                                                        b2,
                                                        bead_list[t][:],
                                                        micelles
                                                    )

        micelles_time.append(micelles)

        micelles_found = len(micelles)
        print("Timestep : {}".format(t))
        print("\nNumber of Micelles found: {}\n".format(micelles_found))

    return micelles_time, bead_list

################################################################################
# Function to output number of micelles at each timestep
def output_num(data_list,f_name_count,num_tsteps):
# INPUTS
# data_list : contains micelle data at each timestep
# f_name_count : name of output file - char
# num_tstepes : number of timesteps in file - int
# RETURNS
# datafile with number of micelles found at each timestep
# micelle_num : number of micelles at each timestep

    micelle_num = []

    with open(f_name_count,'w') as f:
        f.write("T-step \t Micelle Num\n")
        for t in range(num_tsteps):
            count = len(data_list[t][:])
            f.write("{}, {}\n".format(t,count))
            micelle_num.append(count)

    return micelle_num

################################################################################
# Function to find the average micelle size at each timestep and outputs to
# a file
def avg_size(data_list,num_tsteps):
# INPUTS
# data_list : holds micelle data
# num_tsteps : number of timesteps in the data - int
# RETURNS
# t_avg : a list of the average micelle size at each timestep

    t_avg = []
    with open("avg_micelle_sz.dat",'w') as f:
        f.write("T-step\tAvg_sz\n")
        for t in range(num_tsteps):
            num_micelles = len(data_list[t][:])
            avg = 0
            for obj in data_list[t][:]:
                bead_count = len(obj.beads)
                avg += bead_count
            avg /= num_micelles
            f.write("{}, {}\n".format(t,avg))
            t_avg.append(avg)

    return t_avg

################################################################################
# Function to plot data
def plot(f_name,title,x_ax,y_ax,x_data,y_data):
# INPUTS
# f_name : output file name - char
# title : plot title - char
# x_ax : x axis title - char
# y_ax : y axis title - char
# x_data : data for independant values - list
# y_data : data for dependant values - list
# OUTPUTS
# A png file of the plot

    # x range
    x_l = min(x_data)*0.9
    x_h = max(x_data)*1.1
    # y range
    y_l = min(y_data)*0.9
    y_h = max(y_data)*1.1

    plt.plot(x_data,y_data,'o')
    plt.xlim((x_l,x_h))
    plt.ylim((y_l,y_h))
    plt.ylabel(y_ax)
    plt.xlabel(x_ax)
    plt.title(title)
    plt.savefig(f_name,bbox_inches='tight')
    plt.clf()

    return

################################################################################
########### Main ###############################################################
################################################################################
# calls the various subroutines and handles user input

# output file names
f_mpt = "micelles_per_timestep.dat"
f_apt = "avg_micelles_per_timestep.dat"

print("Welcome to Micelle Maker\n")

# User input is checked
while True: # Get interaction range
    try:
        rng = float(input("Please enter the interaction range: \n") or '0')
    except ValueError:
        print("Input not recognized.\n")
        continue
    if rng == 0:
        print("Range entry must be non-zero, please try again.\n")
        continue
    else:
        break

while True: # get name of file with data
    data_file = input("Please enter the name of the file containing the data\n")
    if os.path.isfile(data_file):
        break
    else:
        print("File not found, please try again.\n")
        continue

print("Gathering data ... \n")
main_data, num_tsteps = get_data(data_file)

print("Building Micelles ... \n")
micelles, main_data = micelle_maker(main_data,num_tsteps,rng)

print("Building Averages ... \n")
avgs_tstep = avg_size(micelles,num_tsteps)

print("Outputting data ... \n")
# output and plot data
plot_name = "avg_micelle_time.png"
plot_title = "Micelle Size as a Function of Time"
x_axis = "Time (LJ Timesteps)"
y_axis = "Avg. Micelle Size"
x_data = np.linspace(0,num_tsteps,num=num_tsteps).tolist()
plot(plot_name,plot_title,x_axis,y_axis,x_data,avgs_tstep)

micelle_num = output_num(micelles,f_mpt,num_tsteps)

plot_name = "num_micelle_time.png"
plot_title = "Micelle Number as a Function of Time"
x_axis = "Time (LJ Timesteps)"
y_axis = "Number of Micelles"
plot(plot_name,plot_title,x_axis,y_axis,x_data,micelle_num)

print("All done. Goodbye.")

#! usr/bin/python3

# A code to determine the average size of the hydrophilic associations of a
# mucin network. Input is a pict file from the LAMMPS software suite.

# Jon Parsons
# 11-4-2020

import matplotlib.pyplot as plt
import numpy as np


################################################################################
###################### Classes #################################################
################################################################################
class bead:
# holds information for beads
    def __init__(number,type,x,y,z,micelle):
        self.number = number
        self.type = type
        self.x = x
        self.y = y
        self.z = z
        self.micelle = micelle

################################################################################
class micelle:
# holds which beads are in a particular micelle
    def __init__(beads=None):
        if beads is None:
             beads = []


################################################################################
###################### Functions ###############################################
################################################################################
# Function to get data out of a pict file. Assigns data of appropriate type to
# an array using the bead object
def get_data(f_name):
    t_steps = 0; # timesteps
    data_list = []


    print("\nOpening file {}".format(f_name))
    with open(f_name,'r') as d:
        f_dat = d.read().splitlines()

    line_count = 0
    bd_count = 0
    for line in f_dat:
        line_count = line_count + 1
        d_ln = line.split()
        if line_count == 2:
            t_steps = t_steps + 1
            print("\nObtaining data from timestep {}".format(t_step))
        if line_count == 10:
            beads = []
            for x in range(0,80002):
                if d_ln[1] == 5:
                    beads.append(bead(d_ln[0],d_ln[1],d_ln[2],d_ln[3],d_ln[4],0))
        data_list.append(beads)
        line_count = 0

    print("\nTotal number of Timesteps: {} \n".format(t_steps))

    return data_list, t_steps

################################################################################
# Function that gives the distance between two data points in cartesian space
def distance(x1,y1,z1,x2,y2,z2):

    d = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

    return d

################################################################################
# Function that creates micelles based on distances between beads
def micelle_maker(data_list,timesteps,rng):
    # len - distance at which we consider them bonded

    micelles = []
    m_cnt = 0
    # Should be same number of beads at every timestep barring a disaster
    bead_num = data_list.shape[1]
    for t in range(timesteps):

        for b1 in range(bead_num):
            for b2 in range(bead_num):
                if b1 != b2:
                    dist = distance(b1.x,b1.y,b1.z,b2.x,b2.y,b2.z)
                    if dist < rng:
                        if b1.micelle == 0  and b2.micelle == 0:
                        # make a new one
                            m_cnt = m_cnt + 1
                            b1.micelle = m_cnt
                            b2.micelle = m_cnt
                            micelles.append(micelle.beads.extend((b1.number,b2.number)))
                        if b1.micelle == 0 and b2.micelle != 0:
                            # assign to old one
                            m_old = b2.micelle
                            micelles[m_old].beads.append(b1.number)
                        if b1.micelle != 0 and b2.micelle == 0:
                            # assign to old one
                            m_old = b1.micelle
                            micelles[m_old].beads.append(b2.number)
                        if b1.micelle != 0 and b2.micelle != 0:
                            # combine the two
                            m_keep = b1.micelle
                            m_lose = b2.micelle
                            micelles[m_keep].beads.extend(micelles[m_lose].beads)
                            del micelles[m_lose]

    micelles_found = len(micelles)
    print("\nNumber of Micelles found: {}\n".format(micelles_found))
    print("\nExample Micelle:\n")
    print(micelles[1].beads)

    return micelles

################################################################################
# Function to output our data
def output(data_list,f_name_count,f_name_each):

    with open(f_name_count,'w') as f:
        f.write("T-step \t Micelle Num")
        for item in

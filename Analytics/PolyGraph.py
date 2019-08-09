#!/usr/bin/python3

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os

# Plotting for various outputs of the polymers set.

# Jon Parsons
# 8-1-19

################################################################################

## Function for rolling average

def roll_avg(x, n):
    sum = np.cumsum(np.insert(x,0,0))
    return (sum[n:] - sum[:-n])/float(n)

## find which runs to analise

runs = []
dirs = []

exit_input = ''

while True:
    _=os.system('clear')
    print("Runs to check: \n")
    print(*runs, sep=' ')
    print("\n Enter the run number to check. Leave blank when all desired", \
     "runs are entered. \n")

    num = input("Run Number: ")
    if num == exit_input:
        break
    if num not in runs:
     runs.append(num)

for i in runs:
    directory = "{}data".format(i)
    dirs.append(directory)

## directory management

dir_master = os.getcwd()

time = []
perc = []
per_rsum = []

## Percentages ##
for i in dirs:

    dir_cur = i

    os.chdir(i)

    # Read in data
    percfile = open('Percentages.dat','r')
    perdat = percfile.read().splitlines()
    percfile.close()

    perdat.pop(-1)

    for line in perdat:
        xy = line.split()
        time.append(float(xy[0]))
        perc.append(float(xy[1]))
    os.chdir(dir_master)


sz = len(perc) - 1
for i in range(1,sz):
    x1 = perc[i]
    x2 = perc[i+1]
    avg = (x1 + x2)/2
    per_rsum.append(avg)

# Plot stuff
plt.plot(time,perc)
plt.plot(time[1:sz],per_rsum)
plt.legend(['Percentage','Average'])
plt.xlabel('Time')
plt.ylabel('Percentage of Chains')
plt.title('Percentage of Chains \n in a Cluster')
plt.savefig('Percentage.png')

plt.clf()

os.chdir(dir_master)

################################################################################

##  Networkness ##

network = []
netsum = []

for i in dirs:

    dir_cur = i

    os.chdir(i)

    # Read in data
    netfile = open('network.dat','r')
    netdat = netfile.read().splitlines()
    netfile.close()

    # Convert to graphable form
    ## Time array read in from percentages

    netdat.pop(-1)
    netdat.pop(-1)

    for line in netdat:
        xy = line.split()
        network.append(float(xy[2]))
    os.chdir(dir_master)

netsum.append(0)
sz = len(network) - 1
for i in range(1,sz):
    x1 = netsum[-1] + network[i]
    x2 = network[i+1]
    avg = (x1 + x2)/i
    netsum.append(avg)

netsum.pop(0)
# Plot stuff
plt.plot(time,network)
plt.plot(time[1:sz],netsum)
plt.legend(['Networkness','Average'])
plt.xlabel('Time')
plt.ylabel('Networkness')
plt.title('Network Ratio')

plt.savefig('Networkness.png')

plt.clf()

os.chdir(dir_master)

################################################################################

## Hbonds over time

hbnd = []
time = []
hbnd_rsum = []

for i in dirs:

    dir_cur = i

    os.chdir(i)

    # Read in data
    hbondfile = open('ClustTime_hbond.dat','r')
    hbnddat = hbondfile.read().splitlines()
    hbondfile.close()

    # Convert to graphable form
    hbnddat.pop(0)

    for line in hbnddat:
        xy = line.split()
        time.append(float(xy[0]))
        hbnd.append(float(xy[1]))

    os.chdir(dir_master)


sz = len(hbnd) - 1
for i in range(1,sz):
    x1 = hbnd[i]
    x2 = hbnd[i+1]
    avg = (x1 + x2)/2
    hbnd_rsum.append(avg)

# Plot stuff
plt.plot(time,hbnd)
plt.plot(time[1:sz],hbnd_rsum)
plt.legend(['H-bonds','Average'])
plt.xlabel('Time')
plt.ylabel('Number of H-Bonds')
plt.title('Number of Hydrogen Bonds over Time')

plt.savefig('Hbonds.png')

plt.clf()

os.chdir(dir_master)

################################################################################

## Hbond Histogram
for i in dirs:

    dir_cur = i

    os.chdir(i)

    # Read in data
    hhistfile = open('hist.dat','r')
    hhistdat = hhistfile.read().splitlines()
    hhistfile.close()

    # Process data

    del hhistdat[::11]

    hhist = []

    for line in hhistdat:
        xy = line.split()
        hhist.append(float(xy[1]))

    box = [2,3,4,5,6,7,8,9,10]
    two = []
    three = []
    four = []
    five = []
    six = []
    seven = []
    eight = []
    nine = []
    ten = []
    avg = []
    stdev = []

    two = hhist[2::10]
    three = hhist[3::10]
    four = hhist[4::10]
    five = hhist[5::10]
    six = hhist[6::10]
    seven = hhist[7::10]
    eight = hhist[8::10]
    nine = hhist[9::10]
    ten = hhist[10::10]

    avg.append(np.mean(two))
    avg.append(np.mean(three))
    avg.append(np.mean(four))
    avg.append(np.mean(five))
    avg.append(np.mean(six))
    avg.append(np.mean(seven))
    avg.append(np.mean(eight))
    avg.append(np.mean(nine))
    avg.append(np.mean(ten))

    stdev.append(np.std(two))
    stdev.append(np.std(three))
    stdev.append(np.std(four))
    stdev.append(np.std(five))
    stdev.append(np.std(six))
    stdev.append(np.std(seven))
    stdev.append(np.std(eight))
    stdev.append(np.std(nine))
    stdev.append(np.std(ten))

    N = 9
    width = 0.2
    ind = np.arange(N)

    os.chdir(dir_master)


plt.bar(ind,avg,width,yerr=stdev)

plt.ylabel('Number of Clusters')
plt.xticks(ind,box)
plt.xlabel('Number of Beads')
plt.title('Average Number of Hbond Clusters of Size:')

plt.savefig('hbond_hist.png')

plt.clf()

os.chdir(dir_master)

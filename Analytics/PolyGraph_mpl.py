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

    two = hhist[1::10]
    three = hhist[2::10]
    four = hhist[3::10]
    five = hhist[4::10]
    six = hhist[5::10]
    seven = hhist[6::10]
    eight = hhist[7::10]
    nine = hhist[8::10]
    ten = hhist[9::10]

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

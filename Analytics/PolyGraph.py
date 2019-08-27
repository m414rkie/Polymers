#!/usr/bin/python3
import os
import math
import csv
import subprocess
import statistics

# Plotting for various outputs of the polymers set.
# Outputs the data in x-y format with no headers


# Jon Parsons
# 8-1-19

################################################################################

# Function for average
def ave(lst):
    tot = 0
    for i in lst:
        tot += i

    ave = tot/len(lst)
    return ave

# Function for std. deviation
def stdv(lst,ave):
    sum = 0
    for i in lst:
        sum += (i - ave)*(i - ave)

    sum = sum/len(lst)

    stddev = math.sqrt(sum)
    return stddev

# Function for running average
def run_avg(lst):
    run_avg = []
    N = 1
    sum = 0
    for i in lst:
        sum += i
        ave = sum/N
        run_avg.append(ave)
        N += 1

    return run_avg


################################################################################

## find which runs to analyse

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

################################################################################

time = [0]
perc = []

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
        time.append(time[-1]+1000)
        perc.append(float(xy[1]))

    os.chdir(dir_master)

time.pop(0)

perc_ave = run_avg(perc)

perc_out = "perc_grace.dat"
out = open(perc_out,'w')

writer = csv.writer(out, delimiter='\t')
writer.writerows(zip(time,perc,perc_ave))

out.close()

################################################################################

##  Networkness ##

network = []

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

net_ave = run_avg(network)

net_out = "net_grace.dat"
out = open(net_out,'w')

writer = csv.writer(out, delimiter='\t')
writer.writerows(zip(time,network,net_ave))

out.close()


################################################################################

## Hbonds over time

hbnd = []

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
        hbnd.append(float(xy[1]))

    os.chdir(dir_master)


hbnd_out = "hbnd_grace.dat"
out = open(hbnd_out,'w')

hbnd_ave = run_avg(hbnd)


writer = csv.writer(out, delimiter='\t')
writer.writerows(zip(time,hbnd,hbnd_ave))

out.close()

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

    avg.append(ave(two))
    avg.append(ave(three))
    avg.append(ave(four))
    avg.append(ave(five))
    avg.append(ave(six))
    avg.append(ave(seven))
    avg.append(ave(eight))
    avg.append(ave(nine))
    avg.append(ave(ten))

#    stdev.append(stdv(two,avg[0]))
#    stdev.append(stdv(three,avg[1]))
#    stdev.append(stdv(four,avg[2]))
#    stdev.append(stdv(five,avg[3]))
#    stdev.append(stdv(six,avg[4]))
#    stdev.append(stdv(seven,avg[5]))
#    stdev.append(stdv(eight,avg[6]))
#    stdev.append(stdv(nine,avg[7]))
#    stdev.append(stdv(ten,avg[8]))
    stdev.append(statistics.pstdev(two))
    stdev.append(statistics.pstdev(three))
    stdev.append(statistics.pstdev(four))
    stdev.append(statistics.pstdev(five))
    stdev.append(statistics.pstdev(six))
    stdev.append(statistics.pstdev(seven))
    stdev.append(statistics.pstdev(eight))
    stdev.append(statistics.pstdev(nine))
    stdev.append(statistics.pstdev(ten))


    os.chdir(dir_master)

boxes = list(range(2,10))

hist_out = "hist_grace.dat"
out = open(hist_out,'w')

writer = csv.writer(out, delimiter='\t')
writer.writerows(zip(boxes,avg,perc))

out.close()

subprocess.call("mkdir -p ~/GraceData/", shell=True, executable='/bin/bash')
subprocess.call("mv *.dat GraceData/", shell=True, executable='/bin/bash')

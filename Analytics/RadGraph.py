#!/usr/bin/python3
import os
import csv
import subprocess

# Plotting for various outputs of the polymers set.
# Outputs the data in x-y format with no headers


# Jon Parsons
# 8-1-19

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

r = []
s = []
h = []
sh = []

N = 0

# Gather Data
for i in dirs:

    dir_cur = i
    os.chdir(i)

    # Read data
    s_file = open('rad_dist_h.dat','r')
    s_dat = s_file.read().splitlines()
    s_file.close()

    h_file = open('rad_dist_s.dat','r')
    h_dat = h_file.read().splitlines()
    h_file.close()

    if i == dirs[0]:
        for line in s_dat:
            xy = line.split()
            r.append(float(xy[0]))
            s.append(float(xy[1]))

        for line in h_dat:
            xy = line.split()
            h.append(float(xy[1]))

    N = 0

    if i != dirs[0]:
        for line in s_dat:
            xy = line.split()
            s[N] += float(xy[1])
            N += 1

        N = 0

        for line in h_dat:
            xy = line.split()
            h[N] += float(xy[1])
            N += 1

    os.chdir(dir_master)

for i,v in enumerate(s, start=0):
    s[i] = s[i]/3
    h[i] = h[i]/3


for i,v in enumerate(h, start = 0):
    sh.append(s[i]+h[i])

rad_out = "rad_grace.dat"
out = open(rad_out,'w')

writer = csv.writer(out, delimiter='\t')
writer.writerows(zip(r,s,h,sh))

out.close()

## Find Maxima

max = 0
sz = 0
for i in range(4,25):
    if max <= sh[i]:
        max = sh[i]
        sz = r[i]

print("Maxima at: {}".format(int(round(sz))))

subprocess.call("mkdir -p ~/GraceData/", shell=True, executable='/bin/bash')
subprocess.call("mv *.dat GRaceData/", shell=True, executable='/bin/bash')

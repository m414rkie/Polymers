#!/usr/bin/python3

import subprocess
import os

# Handles the various programs for the polymers set.

# Jon Parsons
# 6-3-19

################################################################################

## Gather which runs to check ##

runs = []

exit_input = ' '

while True:
    _=os.system('clear')
    print("Current runs to check: \n")
    print(*runs, sep=' ')
    print("\n Enter the run number to check. Leave blank when all desired", \
     "runs are entered. \n")

    # Validate input
    try:
        num = int(input("Run Number: ") or '0')
    except ValueError:
        print("\n Integers only please. \n")
        continue
    # Exit if no number entered
    if num == 0:
        break
    if num not in runs:
     runs.append(num)

data_filep = "datafilep.txt"
data_fileb = "datafileb.txt"

# Run Programs
for i in runs:

    # Find datafiles for current run
    bondfile = "pict.a{}".format(i)
    linkfile = "bonds.a{}".format(i)

    in_file_pict = open(data_filep,'w')
    in_file_bond = open(data_fileb,'w')

    in_file_pict.write(bondfile)
    in_file_bond.write(linkfile)

    in_file_pict.close()
    in_file_bond.close()

    subprocess.call("./LinkClus.x < {}".format(data_fileb), shell=True, executable='/bin/bash')
    subprocess.call("./network.x < {}".format(data_fileb), shell=True, executable='/bin/bash')
    subprocess.call("./raddist.x < {}".format(data_filep), shell=True, executable='/bin/bash')

    subprocess.call("mkdir -p ~/{}data/".format(i), shell=True, executable='/bin/bash')
    subprocess.call("mv *.dat {}data/".format(i), shell=True, executable='/bin/bash')

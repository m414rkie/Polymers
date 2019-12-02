#!/usr/bin/python3

from scipy.interpolate import griddata
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import math

# This code takes known output files and plots the data within.

# 1.1 - Replaces the repetitive graphing code with a function
# 1.2 - Polynomial fit of error implemented

# Python code for HW2, COMP670

# Jon Parsons
# 10-1-19

################################################################################

# File names known
err_file = "diff_out.dat"

dt = []
msd = []
# linear constant displacement
x = []
y = []

# Read in the error file
with open(err_file,'r') as f:
    errdat = f.read().splitlines()

# Convert data to usable form
for line in errdat:
    dat = line.split()
    dt.append(float(dat[0]))
    msd.append((float(dat[1])))

t = []
for val in dt:
    t.append(val)


fit = np.polyfit(t,msd,1)

# Define and plot values for the l2 norm
plt.loglog(t, msd,label="y = {0:.2f}T^2 + {0:.2f}".format(fit[0],fit[1]))
plt.legend(loc='upper left')
plt.ylabel('MSD')
plt.xlabel('T')
plt.title('MSD Random Walk')
plt.savefig('rw_2.png',bbox_inches='tight')

plt.clf()

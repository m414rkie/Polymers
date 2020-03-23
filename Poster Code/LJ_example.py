#!usr/bin/python3

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# Program to output a graph of the Lennard Jones potential curve as a simple
# demonstration

################################################################################
# Function containing the LJ equation
def lj(r):
    sigma = 1.0

    V = 4.0*20*(((sigma/r)**12) - ((sigma/r)**6))
    return V

################################################################################
# Main

Radius = np.arange(0.9,5.0,0.01)
zerox = np.arange(0.0,5.0,0.01)
pot = []
zeroes = []

for i in zerox:
    zeroes.append(0)

for i in Radius:
    r_pot = lj(i)
    pot.append(r_pot)

plt.plot(Radius,pot,0)
plt.plot(zerox,zeroes,0,lw=0.5)
plt.xlim((0.5,3.0))
plt.rc('axes',titlesize=22)
plt.rc('axes',labelsize=22)
plt.rc('figure',titlesize=22)
plt.ylabel('Potential Energy', size=20)
plt.xlabel('Radial Distance', size=20)
plt.title('Lennard-Jones Potential Curve')
plt.savefig('LJ_ex.png')

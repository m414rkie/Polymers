#!usr/bin/python3

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# Program to output a graph of the fene interaction curve as a simple
# demonstration

################################################################################
# Function containing the LJ equation
def lj(r):
    r0 = 1.5
    k = 30

    V = -0.5*k*r0*r0*np.log(1.0-((r/r0)**2))
    return V

################################################################################
# Main

Radius = np.arange(0.0,1.4,0.01)
zerox = np.arange(0.0,1.4,0.01)
pot = []
zeroes = []

for i in zerox:
    zeroes.append(0)

for i in Radius:
    r_pot = lj(i)
    pot.append(r_pot)

plt.plot(Radius,pot,0)
plt.plot(zerox,zeroes,0,lw=0.5)
plt.xlim((0.0,1.5))
plt.rc('axes',titlesize=22)
plt.rc('axes',labelsize=22)
plt.rc('figure',titlesize=22)
plt.ylabel('Potential Energy', size=20)
plt.xlabel('Radial Distance', size=20)
plt.title('FENE Potential Curve')
plt.savefig('fene_ex.png')

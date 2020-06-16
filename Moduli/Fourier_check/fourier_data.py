#! usr/bin/python3

# Short program to output data with a known and simple Fourier transform
# Jon Parsons
# 6-6-2020

from matplotlib import pyplot as plt
import numpy as np
import csv
from scipy.fft import fft
import cmath

def outputs(x):
    # x - spatial position

    # pi
    pi = np.pi
    # component frequencies, f = [50,80] hz
    f1 = np.sin(50*pi*x)
    f2 = 0.5*np.sin(80*pi*x)
    f3 = 2*np.cos(300*pi*x)

    out = f1 + f2 + f3

    return out

################################################################################
## Main ##

N = 600
T = 1.0/800.0

# spatial values
x = np.linspace(0.0,N*T,N)
# initialize y vector
y = np.empty(N)
z = np.empty(N)



# file name
f_out = "fourier_test.dat"

# fill values
for i, val in enumerate(x):
    y[i] = outputs(val)

# output
with open(f_out,'w') as fl:
    writer = csv.writer(fl, delimiter = '\t')
    writer.writerows(zip(x,y,z))

# quick graph of outputs

plt.plot(x,y)
plt.title("Fourier Testing Data \n Frequencies of 50, 80, and 300 Hz")
plt.savefig("four_test.jpg")
plt.clf()

yf = fft(y)
xf = np.linspace(0.0,1/(T),N//2)

plt.plot(xf,2/N * np.abs(yf[0:N//2]))
plt.title("Fourier Transform \n Frequencies of 50, 80, and 300 Hz")
plt.savefig("four_trnsfm.jpg")
plt.clf()

# manual implementation bit
n = len(y)
ym = []
for k in range(n):  # For each output element, iterates through the freq's
	s = complex(0)
	for t in range(n):  # For each input element, iterates through data points
		angle = 2j * cmath.pi * t * k / n
		s += y[t] * cmath.exp(-angle)
	ym.append(s)

plt.plot(xf,2/N * np.abs(ym[0:N//2]))
plt.title("Fourier Transform \n manual implementation")
plt.savefig("four_man.jpg")
plt.clf()

## End ##

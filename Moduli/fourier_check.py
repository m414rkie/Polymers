#! usr/bin/python3

# Short program to output data with a known and simple Fourier transform
# Jon Parsons
# 6-6-2020

from matplotlib import pyplot as plt
import numpy as np
import csv
from scipy.fft import fft
import cmath

################################################################################
## Main ##

datfile = open('diff_32j.dat')
dat = datfile.read().splitlines()
datfile.close()

x = []
y = []

for line in dat:
    ls = line.split()
    x.append(float(ls[0]))
    y.append(float(ls[1]))

N = len(x)
T = x[1] - x[0]

print(N, T)

# file name
f_out = "fourier_check.dat"

# quick graph of outputs

yf = fft(y)
xf = np.linspace(0.0,1/(T),N//2)
plt.plot(xf,2/N * np.abs(yf[0:N//2].real))
plt.plot(xf,2/N * np.abs(yf[0:N//2].imag))
plt.ylim(-50,200)
plt.xlim(0.0,0.5)
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
plt.ylim(-50,200)
plt.xlim(0.0,0.5)
plt.title("Fourier Transform \n manual implementation")
plt.savefig("four_man.jpg")
plt.clf()

## End ##

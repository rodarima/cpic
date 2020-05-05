import numpy as np
from math import *
import matplotlib.pyplot as plt

L=8
N=64
e0 = 1.0
nperiods = 10

E0 = np.genfromtxt('E.csv', delimiter=',', skip_header=1)[:,0]

q = np.array([-1, -1])
qx = np.array([1, 7])

E = np.zeros(N)
r = np.arange(N) * L / N

for i in range(len(q)):
	E += q[i] / (e0 * 2.0 * pi * (r - qx[i]))

for p in range(1, nperiods):
	for i in range(len(q)):
		E += q[i] / (e0 * 2.0 * pi * (r - qx[i] + p * L))
	for i in range(len(q)):
		E += q[i] / (e0 * 2.0 * pi * (r - qx[i] - p * L))

print(E[7])
print(E[8])
print(E[9])
print(E0[7])
print(E0[8])
print(E0[9])

#plt.plot(r, E)
#plt.plot(r, E0)
#plt.grid()
plt.show()

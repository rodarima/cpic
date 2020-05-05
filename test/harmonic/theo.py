from math import *
import numpy as np
import matplotlib.pyplot as plt

E = np.genfromtxt('harm.E0x')
rraw = np.genfromtxt('harm.r0x')
r = 2 - rraw

L = 8
q = -1
e0 = 1
nc = 516.9 #cycles per period
f = 1/nc
nt = r.shape[0]
t = np.arange(nt)

rr = np.cos(t*f*2*pi)
dr = r - rr

dE = np.diff(E)
EE = - max(E) * np.cos((t-1)*f*2*pi)


def dist(i, alpha):
	return - (i + alpha) * L

def getEi(i, x):
	alpha = (L - 2.0 * x) / L
	if np.any(alpha < 0) or np.any(alpha > 1.0):
		print("ERROR alpha = {}".format(alpha))
		exit(1)
	d = dist(i, alpha)
	#print(i,alpha,d)
	# Coulomb law for 2 dimensions, see equation (1) in this paper:
	# https://arxiv.org/pdf/0912.0225.pdf
	return q / (e0 * 2 * pi * d)

def getE(x, maxit):
	E0 = getEi(0, x)
	E = E0
	for i in range(1, maxit):
		E += getEi(i, x) + getEi(-i, x)
	return E, E0

EEE, EEE0 = getE(rraw, 500)

plt.rcParams['axes.grid'] = True

plt.margins()
plt.subplots_adjust(hspace=0.3)
plt.subplot(3, 2, 1)
plt.title('Experimental position $x(t)$')
plt.plot(r, color='blue')
plt.subplot(3, 2, 3)
plt.title('Theoretical position $x_0(t) = \cos(\omega_p t)$')
plt.plot(rr, color='green')
plt.subplot(3, 2, 5)
plt.title('Position error $x(t) - x_0(t)$')
plt.plot(dr, color='red')

plt.subplot(3, 2, 2)
plt.title('Experimental electric field $E(t)$')
plt.plot(E, color='blue')
plt.subplot(3, 2, 4)
plt.title('Theoretical electric field $E_0(t)$')
plt.plot(EEE, color='green')
plt.subplot(3, 2, 6)
plt.title('Electric field error $E(t) - E_0(t)$')
plt.plot(E - EEE, color='red')
#plt.subplot(4, 2, 8)
#plt.title('Experimental acceleration $d^2x/dt^2$')
#plt.plot(np.diff(np.diff(r)), color='brown')

plt.show()

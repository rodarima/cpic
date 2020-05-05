from math import *
import numpy as np
import matplotlib.pyplot as plt

E = np.genfromtxt('harm.E0x')
r = 2 - np.genfromtxt('harm.r0x')

nc = 516.9 #cycles per period
f = 1/nc
nt = r.shape[0]
t = np.arange(nt)

rr = np.cos(t*f*2*pi)
dr = r - rr

dE = np.diff(E)
EE = - max(E) * np.cos((t-1)*f*2*pi)

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
plt.title('Theoretical electric field $E_0(t) = - A \cos(\omega_p t)$')
plt.plot(EE, color='green')
plt.subplot(3, 2, 6)
plt.title('Electric field error $E(t) - E_0(t)$')
plt.plot(E - EE, color='red')
#plt.subplot(4, 2, 8)
#plt.title('Experimental acceleration $d^2x/dt^2$')
#plt.plot(np.diff(np.diff(r)), color='brown')

plt.show()

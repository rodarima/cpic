import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(n, a, b):
	return a * n * np.log(n) + b

data = np.genfromtxt('time.csv', delimiter=' ')

x = data[:,0]#.reshape((-1, 1))
y = data[:,3]#.reshape((-1, 1))
std = data[:,4].reshape((-1, 1))

#print(x)

#print(y)

popt, pcov = curve_fit(func, x, y)
print(popt)

#plt.plot(x, y, 'b-', label='data')
#plt.plot(x, func(x, *popt), 'r-', label='fit: a=%e, b=%e' % tuple(popt))
#plt.legend()
#plt.show()


x = x.reshape((-1, 1))
y = y.reshape((-1, 1))

yy = func(x, *popt)

print(yy)

residuals = y - yy 

reg = np.concatenate([x, yy, residuals, std], axis=1)


np.savetxt("csv/regression.csv", reg, delimiter=' ')

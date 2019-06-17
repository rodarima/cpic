#!/usr/bin/env python

import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


data = np.genfromtxt('working.log', delimiter='\t')


t = data[:,0].reshape(-1, 1)
percent = data[:,2].reshape(-1, 1)
#print(t)
#print(percent)

reg = LinearRegression().fit(t, percent)
reg.score(t, percent)

#print(reg.coef_)

#print(reg.intercept_)

t0 = np.min(t)
tf = 1561326600
tm = 1560376200
tmm = 1560894600

tt = np.linspace(t0, tf, 100).reshape(-1, 1)

pf = reg.predict(tt)

plt.plot(tt, pf, 'g')
plt.scatter(t, percent, s=2, color='k')
plt.axvline(x=tm, color='k', linewidth=1.0)
plt.axvline(x=tmm, color='brown', linewidth=1.0)
plt.axvline(x=tf, color='r', linewidth=3.0)
plt.xlim(np.min(t), tf+3600)
plt.ylim(0, 100)
plt.grid(True)
#plt.show()
plt.savefig("prediction.png")

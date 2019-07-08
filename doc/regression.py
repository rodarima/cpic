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
tm = 1560376200
tmm = 1560894600
tf = 1561424400 #1561413000
tF = 1561672200 #1561585800

tt = np.linspace(t0, tF, 100).reshape(-1, 1)

pf = reg.predict(tt)

plt.plot(tt, pf, 'g')
plt.axvline(x=tm, color='k', linewidth=1.0)
plt.axvline(x=tmm, color='brown', linewidth=1.0)
plt.axvline(x=tf, color='orangered', linewidth=1.0)
plt.axvline(x=tF, color='r', linewidth=3.0)
plt.scatter(t, percent, s=2, color='k', zorder=100)
plt.xlim(np.min(t), tF+3600*24)
plt.ylim(0, 100)
plt.grid(True)
#plt.show()
plt.savefig("prediction.png")

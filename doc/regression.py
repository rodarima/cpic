import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt


data = np.genfromtxt('working.log', delimiter='\t')


t = data[:,0].reshape(-1, 1)
percent = data[:,2].reshape(-1, 1)
print(t)
print(percent)

reg = LinearRegression().fit(t, percent)
reg.score(t, percent)

print(reg.coef_)

print(reg.intercept_)

t0 = np.min(t)
tf = 1560894600

tt = np.linspace(t0, tf, 100).reshape(-1, 1)

pf = reg.predict(tt)

plt.plot(t, percent, 'bo')
plt.plot(tt, pf, 'r')
plt.xlim(np.min(t), tf)
plt.ylim(0, 100)
plt.show()

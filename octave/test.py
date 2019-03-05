import numpy as np
import matplotlib
matplotlib.use('tkAgg')
import matplotlib.pyplot as plt


N = 50
L = 2*N
q = -1;
dx = L/N;
periodic = 0;

ones = np.zeros(N) + 1

A = np.diag(ones[:N-1], k=-1) -3 * np.diag(ones) + np.diag(ones[:N-1], k=+1)
b = np.zeros(N)

b[int(N/2)] = 1

x = np.linalg.solve(A, b)

print(x)

plt.plot(x)
plt.show()

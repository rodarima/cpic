import numpy as np
from sklearn.linear_model import LinearRegression

data = np.genfromtxt('time.csv', delimiter=' ')

x = data[:,0].reshape((-1, 1))
y = data[:,3].reshape((-1, 1))
std = data[:,4].reshape((-1, 1))

model = LinearRegression()
model.fit(x, y)

r_sq = model.score(x, y)

print('coefficient of determination:', r_sq)

print('intercept:', model.intercept_)

print('slope:', model.coef_)

yy = model.coef_ * x + model.intercept_

residuals = yy - y

reg = np.concatenate([x, yy, residuals, std], axis=1)


np.savetxt("csv/regression.csv", reg, delimiter=' ')

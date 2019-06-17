import numpy as np
import fileinput

iteration = []
iteration_timer = []

for line in fileinput.input():
	parts = line.split(' ')

	if len(parts) != 4:
		continue

	i = int(parts[1])
	t = float(parts[3])

	if i == 0:
		continue

	iteration.append(i)
	iteration_timer.append(t)

i = np.array(iteration)
t = np.array(iteration_timer)

t_mean = np.mean(t)
t_std = np.std(t)
n = len(iteration)

print("For {} iterations: mean={:e} dev={:e}".format(n, t_mean, t_std))

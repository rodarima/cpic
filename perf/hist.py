#!/usr/bin/env python

import numpy as np
import sys

def get_field(line, name):
	for c in line.split():
		w = c.split('=')
		if w[0] == name:
			return w[1]

	return None


times = []
fn = sys.argv[1]

with open(fn, 'r') as f:
	for line in f:
		parts = line.split()
		if parts[0] != 'stats': continue

		last = float(get_field(line, "last"))
		if last == 0: continue
		times.append(last)

times = np.array(times)
print(times)

import matplotlib.pyplot as plt

# An "interface" to matplotlib.axes.Axes.hist() method
n, bins, patches = plt.hist(x=times, bins=80, alpha=0.5, histtype='bar', ec='black')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Simulation time per iteration')
#plt.text(23, 45, r'$\mu=15, b=3$')
maxfreq = n.max()
# Set a clean upper y-axis limit.
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)

plt.show()

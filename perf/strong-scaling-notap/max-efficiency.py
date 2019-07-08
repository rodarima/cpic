#stats iter=99 last=3.046556e-01 mean=3.424577e-01 std=8.783869e-02 sem=8.828121e-03 rsem=2.577872e-02 mem=6715500 solver=1.370335e-01

ts = 1.370335e-01
t = 3.424577e-01
C = 48
P = 16

tp = t - ts
alpha = ts / (tp*C + ts)

print("alpha = %e" % alpha)

for i in range(5):
	p = 2 ** i

	smax = 1 / (alpha + (1-alpha)/p)

	print("%d %e %e" % (p, smax, smax/p))

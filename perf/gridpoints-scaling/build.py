import sh
import math

BASE_CONF = 'base.conf'
BASE_JOB = 'base.job'

CP_list = [1,16,32]
C_list = [32*i for i in [1,2,4,8,16,32]]
Nc = 64
ng_list = [1024*i for i in [1,2,4,8,16]]

for C in C_list:
	for CP in CP_list:
		for ng in ng_list:

			Ng = ng ** 2

			P = int(C / CP)
			N = int(math.ceil(C / 32))

			print("Conf:\tN=%d\tCP=%d\tP=%d\tC=%d\tng=%d" %
					(N, CP, P, C, ng))

			name = 'C%d-CP%d-ng%d' % (C, CP, ng)
			conf = 'conf/' + name
			job = 'jobs/' + name

			# Set the config
			sh.sed('s/@NCHUNKS@/{}/g;s/@NPOINTS@/{}/g'.format(Nc, ng),
					BASE_CONF, _out=conf)

			# Set the jobs
			sh.sed('s/@NAME@/{}/g;s/@NCPUS@/{}/g;s/@NNODES@/{}/g;s/@NPROCS@/{}/g'.format(
					name, CP, N, P), BASE_JOB, _out=job)

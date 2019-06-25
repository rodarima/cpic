import sh
import math

BASE_CONF = 'base.conf'
BASE_JOB = 'base.job'

P_list = [1,2,4,8,16,32]
CP_list = [1,16,48]
Nc = 64
ng_list = [1024*i for i in [1,2,4,8,16]]

for P in P_list:
	for CP in CP_list:
		for ng in ng_list:

			Ng = ng ** 2

			C = P*CP
			N = int(math.ceil(C / 48))

			print("Conf:\tN=%d\tCP=%d\tP=%d\tC=%d\tng=%d" %
					(N, CP, P, C, ng))

			name = 'P%d-CP%d-ng%d' % (P, CP, ng)
			conf = 'conf/' + name
			job = 'jobs/' + name

			# Set the config
			sh.sed('s/@NCHUNKS@/{}/g;s/@NPOINTS@/{}/g'.format(Nc, ng),
					BASE_CONF, _out=conf)

			# Set the jobs
			sh.sed('s/@NAME@/{}/g;s/@NCPUS@/{}/g;s/@NNODES@/{}/g;s/@NPROCS@/{}/g'.format(
					name, CP, N, P), BASE_JOB, _out=job)

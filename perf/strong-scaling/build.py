import sh
import math

BASE_CONF = 'base.conf'
BASE_JOB = 'base.job'

P_list = [1,2,4,8,16,32]
CP_list = [1,2,4,8,16,32,48]
Nc_list = [32,64,128]

for P in P_list:
	for CP in CP_list:
		for Nc in Nc_list:

			C = P*CP
			N = int(math.ceil(C / 48))

			print("Conf:\tN=%d\tCP=%d\tP=%d\tC=%d\tNc=%d" %
					(N, CP, P, C, Nc))

			name = 'P%d-CP%d-Nc%d' % (P, CP, Nc)
			conf = 'conf/' + name
			job = 'jobs/' + name

			# Set the config
			sh.sed('s/@NCHUNKS@/{}/g'.format(Nc),
					BASE_CONF, _out=conf)

			# Set the jobs
			sh.sed('s/@NAME@/{}/g;s/@NCPUS@/{}/g;s/@NNODES@/{}/g;s/@NPROCS@/{}/g'.format(
					name, CP, N, P), BASE_JOB, _out=job)

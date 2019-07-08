import sh

BASE_CONF = 'base.conf'
BASE_JOB = 'base.job'

P_list = [1,2,4,8,16,32]
NCPUS = 48
Nc = 128
Nw = 32

for P in P_list:

	name = 'P{}'.format(P)
	conf = 'conf/' + name
	job = 'jobs/' + name

	print("Building conf for P=%d" % P)

	# Set the number of plasma chunks
	sh.sed('s/@NCHUNKS@/{}/g;s/@NWORKERS@/{}/g'.format(Nc, Nw),
			BASE_CONF, _out=conf)

	# Set the number of processes
	sh.sed('s/@NAME@/{}/g;s/@NCPUS@/{}/g;s/@NPROCS@/{}/g'.format(
			name, NCPUS, P), BASE_JOB, _out=job)

import sh

#BASE_CONF = 'base.conf'
BASE_JOB = 'base.job'

P_list = [1,2,4,8,16,32,64]
NCPUS = 1
Nc = 128

for P in P_list:
	print("Building conf for P={}".format(P))

	name = 'P{}'.format(P, Nc)
	conf = 'conf/' + name
	job = 'jobs/' + name

	# Set the number of plasma chunks
	#sh.sed('s/@NCHUNKS@/{}/g'.format(Nc),
	#		BASE_CONF, _out=conf)

	# Set the number of processes
	sh.sed('s/@NAME@/{}/g;s/@NCPUS@/{}/g;s/@NPROCS@/{}/g'.format(
			name, NCPUS, P), BASE_JOB, _out=job)

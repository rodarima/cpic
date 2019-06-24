import sh

BASE_CONF = 'base.conf'
BASE_JOB = 'base.job'

Nc_list = [32,64,128,256,512]
P_list = [1,2,4,8,16]
NCPUS = 48

for Nc in Nc_list:
	for P in P_list:
		print("Building conf for P={} Nc={}".format(P, Nc))

		name = 'P{}-Nc{}'.format(P, Nc)
		conf = 'conf/' + name
		job = 'jobs/' + name

		# Set the number of plasma chunks
		sh.sed('s/@NCHUNKS@/{}/g'.format(Nc),
				BASE_CONF, _out=conf)

		# Set the number of processes
		sh.sed('s/@NAME@/{}/g;s/@NCPUS@/{}/g;s/@NPROCS@/{}/g'.format(
				name, NCPUS, P), BASE_JOB, _out=job)

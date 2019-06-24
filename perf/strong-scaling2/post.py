import sh
from io import StringIO


BASE_CONF = 'base.conf'
BASE_JOB = 'base.job'

P_list = [1,2,4,8,16,32,64]
NCPUS = 48
Nc = 128

def get_field(line, name):
	for c in line.split():
		w = c.split('=')
		if w[0] == name:
			return w[1]

	return None


out = 'csv/time.csv'
t0 = 0
with open(out, 'w+') as f:

	f.write("P\tmean\trel-err\tspeedup\tefficiency\n")

	for P in P_list:

		name = 'P{}'.format(P)
		out = 'out/' + name

		try:
			#sh.tail('-1', sh.grep('^stats', out), _out=buf)
			buf = StringIO()
			sh.tail(sh.grep('^stats', out), '-1', _out=buf)
			line = buf.getvalue()
		except:
			line = ""

		if line == "": continue

		mean = float(get_field(line, "mean"))
		sem = float(get_field(line, "sem"))
		rel_error = sem * 1.96 / mean
		if P == 1:
			speedup = 1
			t0 = mean
		else:
			speedup = t0 / mean

		efficiency = speedup / P

		f.write("%d\t%e\t%e\t%e\t%e\n" %
				(P, mean, rel_error, speedup, efficiency))

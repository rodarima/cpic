import sh
from io import StringIO


BASE_CONF = 'base.conf'
BASE_JOB = 'base.job'

Nc_list = [32,64,128,256,512]
P_list = [1,2,4,8,16]
NCPUS = 48

def get_field(line, name):
	for c in line.split():
		w = c.split('=')
		if w[0] == name:
			return w[1]

	return None

for Nc in Nc_list:

	out = 'csv/Nc%d' % Nc
	with open(out, 'w+') as f:

		f.write("P\tNc\tmean\trel-err\n")

		for P in P_list:

			name = 'P{}-Nc{}'.format(P, Nc)
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

			f.write("%d\t%d\t%e\t%e\n" % (P, Nc, mean, rel_error))

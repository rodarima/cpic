import sh
import math
from io import StringIO

BASE_CONF = 'base.conf'
BASE_JOB = 'base.job'

C_list = [32*i for i in [1,2,4,8,16,32]]
CP_list = [1,16,32]
Nc = 64
ng_list = [1024*i for i in [1,2,4,8,16]]


def get_field(line, name):
	for c in line.split():
		w = c.split('=')
		if w[0] == name:
			return w[1]

	return None

for CP in CP_list:
	for ng in ng_list:

		out = 'csv/CP%d-ng%d' % (CP, ng)

		t0 = 0
		with open(out, 'w+') as f:

			f.write("C\tmean\trel-err\tspeedup\tefficiency\n")

			for C in C_list:

				P = int(C/CP)
				N = int(math.ceil(C / 32))

				name = 'C%d-CP%d-ng%d' % (C, CP, ng)
				out = 'out/' + name
				print('Reading %s' % out)

				try:
					#sh.tail('-1', sh.grep('^stats', out), _out=buf)
					buf = StringIO()
					sh.tail(sh.grep('^stats', out), '-1', _out=buf)
					line = buf.getvalue()
				except:
					line = ""

				if line == "":
					continue

				mean = float(get_field(line, "mean"))
				sem = float(get_field(line, "sem"))

				if mean == 0.0:
					print('Cannot get a proper mean')
					continue

				rel_error = sem * 1.96 / mean

				if C == 32:
					speedup = 1
					t0 = mean
				else:
					speedup = t0 / mean

				efficiency = speedup / (C/32)

				f.write("%d\t%e\t%e\t%e\t%e\n" %
						(C, mean, rel_error, speedup, efficiency))

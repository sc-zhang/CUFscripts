#!/usr/bin/env python
import sys


def gen_data(in_beds, in_blast, th, out_bed):
	blast_db = {}
	with open(in_blast, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			chrn = data[1]
			rsp = int(data[6])
			rep = int(data[7])
			if rsp > rep:
				rsp, rep = rep, rsp
			qsp = int(data[8])
			qep = int(data[9])
			if qsp > qep:
				qsp, qep = qep, qsp
			iden = float(data[2])
			if chrn not in blast_db:
				blast_db[chrn] = {}
			key = "%d-%d"%(qsp, qep)
			blast_db[chrn][key] = [rsp, rep, iden]
	
	with open(out_bed, 'w') as fout:
		for bed in in_beds.split(','):
			type = bed.split('.')[-2]
			with open(bed, 'r') as fin:
				for line in fin:
					data = line.strip().split()
					chrn = data[0]
					sp = int(data[1])
					ep = int(data[2])
					if sp > ep:
						sp, ep = ep, sp
					if chrn[:3].lower() != 'chr':
						continue
					key = "%d-%d"%(sp, ep)
					if key not in blast_db[chrn]:
						print("%s\t%d\t%d\t"%(chrn, sp, ep))
						continue
					rsp, rep, iden = blast_db[chrn][key]
					if iden < th:
						continue
					fout.write("%s\t%d\t%d\t%d\t%d\t%f\t%s\n"%(chrn, sp, ep, rsp, rep, iden, type))


if __name__ == "__main__":
	if len(sys.argv) < 5:
		print("Usage: python "+sys.argv[0]+" <in_beds> <in_blast> <threshold> <out_bed>")
	else:
		in_beds, in_blast, th, out_bed = sys.argv[1:]
		gen_data(in_beds, in_blast, float(th), out_bed)

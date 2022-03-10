#!/usr/bin/env python
import sys


def filter_blast_result(in_blast, out_bed):
	blast_db = {}
	with open(in_blast, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			chrn = data[1]
			sp = int(data[8])
			ep = int(data[9])
			if sp <= ep:
				dir = "+"
			else:
				sp, ep = ep, sp
				dir = "-"
			if chrn not in blast_db:
				blast_db[chrn] = {}
			is_skip = False
			pop_list = []
			for key in blast_db[chrn]:
				dsp, dep, ddir =key.split("|")
				dsp = int(dsp)
				dep = int(dep)
				if sp < dsp and ep > dep:
					pop_list.append(key)
				elif sp >dsp and ep < dep:
					is_skip = True
					break
				elif sp == dsp and ep == dep:
					if dir == "-" or ddir == "+":
						is_skip = True
						break
					else:
						pop_list.append(key)
			for key in pop_list:
				blast_db[chrn].pop(key)
			if is_skip == True:
				continue
			key = "%d|%d|%s"%(sp, ep, dir)
			blast_db[chrn][key] = ""
	bed_db = {}
	for chrn in blast_db:
		bed_db[chrn] = []
		for key in blast_db[chrn]:
			sp, ep, dir = key.split('|')
			sp = int(sp)
			ep = int(ep)
			if dir == '-':
				sp, ep = ep, sp
			bed_db[chrn].append([sp, ep])
	with open(out_bed, 'w') as fout:
		for chrn in sorted(bed_db):
			for data in sorted(bed_db[chrn]):
				fout.write("%s\t%d\t%d\n"%(chrn, data[0], data[1]))


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python "+sys.argv[0]+" <in_blast> <out_bed>")
	else:
		in_blast, out_bed = sys.argv[1:]
		filter_blast_result(in_blast, out_bed)


#!/usr/bin/env python
import sys


def stat_snp_indel(in_vcf, out_stat):
	snp_db = {}
	ins_db = {}
	del_db = {}
	with open(in_vcf, 'r') as fin:
		for line in fin:
			if line[0] == '#':
				continue
			data = line.strip().split()
			chrn = data[0]
			ref_len = len(data[3])
			alt_len = len(data[4])
			filter = data[6].lower()
			if filter == 'lowqual':
				continue
			if chrn[:3] == 'tig':
				chrn = 'contig'
			if chrn not in snp_db:
				snp_db[chrn] = 0
			if chrn not in ins_db:
				ins_db[chrn] = 0
			if chrn not in del_db:
				del_db[chrn] = 0
			if ref_len > alt_len:
				del_db[chrn] += 1
			elif ref_len < alt_len:
				ins_db[chrn] += 1
			else:
				snp_db[chrn] += 1
		
	with open(out_stat, 'w') as fout:
		fout.write("##########SNP#############\n")
		for chrn in sorted(snp_db):
			fout.write("%s\t%d\n"%(chrn, snp_db[chrn]))
		fout.write("##########INS#############\n")
		for chrn in sorted(ins_db):
			fout.write("%s\t%d\n"%(chrn, ins_db[chrn]))
		fout.write("##########DEL#############\n")
		for chrn in sorted(del_db):
			fout.write("%s\t%d\n"%(chrn, del_db[chrn]))


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python "+sys.argv[0]+" <in_vcf> <out_stat>")
	else:
		in_vcf, out_stat = sys.argv[1:]
		stat_snp_indel(in_vcf, out_stat)

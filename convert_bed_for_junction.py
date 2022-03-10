#!/usr/bin/env python
import sys


def convert_bed(in_bed, j_size, out_bed):
	with open(in_bed, 'r') as fin:
		with open(out_bed, 'w') as fout:
			for line in fin:
				data = line.strip().split()
				sp, ep = list(map(int, data[1: 3]))
				ssp = sp-j_size
				if ssp <= 0:
					ssp = 1
				sep = sp+j_size
				esp = ep-j_size
				if esp <= 0:
					esp = 1
				eep = ep+j_size
				fout.write("%s\t%d\t%d\t%d\t%d\n"%(data[0], ssp, sep, sp, ep))
				fout.write("%s\t%d\t%d\t%d\t%d\n"%(data[0], esp, eep, sp, ep))


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python "+sys.argv[0]+" <in_bed> <junction_size> <out_bed>")
	else:
		in_bed, j_size, out_bed = sys.argv[1:]
		j_size = int(int(j_size)/2)
		convert_bed(in_bed, j_size, out_bed)

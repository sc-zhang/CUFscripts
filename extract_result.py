#!/usr/bin/env python
import sys


def extract_res(in_res, th):
	uniq_fn = '.'.join(in_res.split('.')[:-1])+'.uniq.bed'
	comm_fn = '.'.join(in_res.split('.')[:-1])+'.common.bed'
	uniq_list = []
	comm_list = []
	with open(in_res, 'r') as fin:
		tmp_list = []
		for line in fin:
			data = line.strip().split()
			tmp_list.append([data[0], data[3], data[4], int(data[5])])
			if len(tmp_list) == 2:
				if tmp_list[0][-1] < th or tmp_list[1][-1] < th:
					uniq_list.append('\t'.join(tmp_list[0][:3]))
				else:
					comm_list.append('\t'.join(tmp_list[0][:3]))
				tmp_list = []
	with open(uniq_fn, 'w') as fout:
		fout.write("\n".join(uniq_list))
	
	with open(comm_fn, 'w') as fout:
		fout.write("\n".join(comm_list))


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python "+sys.argv[0]+" <in_result> <threshold>")
	else:
		in_res, th = sys.argv[1:]
		th = int(th)
		extract_res(in_res, th)

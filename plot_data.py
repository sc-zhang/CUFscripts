#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys


def plot_data(in_data, chr_len, out_pic):
	chr_len_db = {}
	chr_list = []
	print("Getting chr length")
	with open(chr_len, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			if data[0][:3].lower() != 'chr':
				continue
			chr_list.append(data[0])
			chr_len_db[data[0]] = int(data[1])
	base_x = []
	base = 0
	i = 0
	chrn_idx = {}
	for chrn in chr_list:
		base_x.append(base)
		base += chr_len_db[chrn]
		chrn_idx[chrn] = i
		i += 1
	base_x.append(base)
	
	x_ticks = []
	for i in range(0, len(base_x)-1):
		x_ticks.append((base_x[i]+base_x[i+1])/2)
	
	color_db = {'common98': 'red', 'common90': 'blue', 'uniq98': 'orange', 'uniq90': 'cyan'}
	label_db = {'common98': 'Common 100%-98%', 'common90': 'Common 98%-90%', 'uniq98': 'Uniq 100%-98%', 'uniq90': 'Uniq 98%-90%'}
	
	print("Reading data")
	data_db = {}
	max_x = base_x[-1]
	max_y = 0
	with open(in_data, 'r') as fin:
		for line in fin:
			chrn, qsp, qep, rsp, rep, iden, qt = line.strip().split()
			base = base_x[chrn_idx[chrn]]
			if float(iden) >= 98:
				type = qt+'98'
			else:
				type = qt+'90'
			qsp = int(qsp)
			qep = int(qep)
			rsp = int(rsp)
			rep = int(rep)
			
			w = qep-qsp
			h = rep-rsp
			if rep > max_y:
				max_y = rep
			if type not in data_db:
				data_db[type] = []
			data_db[type].append([qsp+base, rsp, w, h])
	
	print("Plotting")
	plt.figure(figsize=(15, 5), dpi=300)
	plt.xlim(0, max_x)
	plt.ylim(0, max_y)
	y_ticks = []
	for i in range(0, max_y, 60000):
		y_ticks.append(i)
	y_labels = []
	for t in y_ticks:
		y_labels.append("%d0k"%(t//10000))
	plt.yticks(y_ticks[1:], y_labels[1:], fontsize=20)
	plt.xticks(x_ticks, chr_list, fontsize=20)
	for i in base_x[1: -1]:
		plt.plot([i, i], [0, max_y], color='grey', alpha=0.5, lw=0.5)

	for type in data_db:
		for x, y, w, h in data_db[type]:
			if h < 500:
				continue
			plt.gca().add_patch(plt.Rectangle((x, y), w, h, fill=True, facecolor=color_db[type], edgecolor=color_db[type], lw=1))
	
	for type in color_db:
		plt.scatter(0, 0, color=color_db[type], label=label_db[type], marker='s', s=0.1, alpha=1)

	plt.legend(fontsize=20, markerscale=25, loc=[1.01, 0.5], frameon=False)
	plt.savefig(out_pic, filetype=out_pic.split('.')[-1], bbox_inches='tight')
	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python "+sys.argv[0]+" <in_data> <chr_len> <out_pic>")
	else:
		in_data, chr_len, out_pic = sys.argv[1:]
		plot_data(in_data, chr_len, out_pic)

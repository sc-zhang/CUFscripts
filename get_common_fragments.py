#!/usr/bin/env python
import sys
import re


def read_gff(in_gff, rq):
	gff_db = {}
	gene_db = {}
	with open(in_gff, 'r') as fin:
		for line in fin:
			if line.strip() == '':
				continue
			data = line.strip().split()
			if data[2] != 'gene':
				continue
			chrn = data[0]
			if chrn[:3].lower() == 'chr':
				chrn = int(re.findall('[Cc]hr([0-9]*)', chrn)[0])
			sp = int(data[3])
			ep = int(data[4])
			id = data[8].split(';')[0].split('=')[1]
			gene_db[id] = [chrn, sp, ep]
			if chrn not in gff_db:
				gff_db[chrn] = []
			gff_db[chrn].append([sp, ep, id])
	if rq.lower() == 'r':
		return gff_db
	else:
		return gene_db


def read_bed(in_bed):
	bed_db = {}
	with open(in_bed, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			chrn = data[0]
			if chrn[:3].lower() == 'chr':
				chrn = int(re.findall('[Cc]hr([0-9]*)', chrn)[0])
			sp = int(data[1])
			ep = int(data[2])
			if sp > ep:
				sp, ep = ep, sp
			if chrn not in bed_db:
				bed_db[chrn] = []
			bed_db[chrn].append([sp, ep])
	return bed_db


def read_col(in_col):
	col_db = {}
	with open(in_col, 'r') as fin:
		for line in fin:
			if line[0] == '#':
				continue
			gene1, gene2 = line.strip().split()[:2]
			if gene1 not in col_db:
				col_db[gene1] = gene2
	return col_db


def get_gene_id_pos(gene_list, pos, lr):
	s = 0
	e = len(gene_list)-1
	while s<=e:
		mid = int((s+e)/2)
		if gene_list[mid][0] < pos:
			s = mid+1
		elif gene_list[mid][0] > pos:
			e = mid-1
		else:
			return mid
	if gene_list[e][1] >= pos:
		return e
	if lr.lower() == 'l':
		return e
	else:
		return s


def get_common_fragments(r_gff, q_gff, r_bed, q_bed, coll, out_reg):
	r_gff_db = read_gff(r_gff, 'r')
	r_bed_db = read_bed(r_bed)
	q_gff_db = read_gff(q_gff, 'q')
	q_bed_db = read_bed(q_bed)
	col_db = read_col(coll)
	
	pos_db = []
	for chrn in r_bed_db:
		if chrn not in r_gff_db:
			continue
		for r_pos in r_bed_db[chrn]:
			spos = get_gene_id_pos(r_gff_db[chrn], r_pos[0], 'l')
			epos = get_gene_id_pos(r_gff_db[chrn], r_pos[1], 'r')
			sp = spos-2
			ep = epos+3
			if sp < 0:
				sp = 0
			if ep >= len(r_gff_db[chrn]):
				ep = len(r_gff_db[chrn])-1
			bef_genes = []
			aft_genes = []
			for i in range(sp, spos+1):
				gene = r_gff_db[chrn][i][2]
				if gene not in col_db:
					continue
				gene2 = col_db[gene]
				if gene2 not in q_gff_db:
					continue
				if 'G' in gene2:
					res = re.findall(r'.*\.([0-9]*)G([0-9]*)', gene2)
					cn = int(res[0][0])
					gn = res[0][1]
				else:
					cn = chrn
					gn = re.findall(r'[A-Za-z]*([0-9]*)', gene2)[0]
				bef_genes.append([cn, int(gn), gene2])
			for i in range(epos, ep):
				gene = r_gff_db[chrn][i][2]
				if gene not in col_db:
					continue
				gene2 = col_db[gene]
				if 'G' in gene2:
					res = re.findall(r'.*\.([0-9]*)G([0-9]*)', gene2)
					cn = int(res[0][0])
					gn = res[0][1]
				else:
					cn = chrn
					gn = re.findall(r'[A-Za-z]*([0-9]*)', gene2)[0]
				aft_genes.append([cn, int(gn), gene2])
			#if len(bef_genes) != 3 or len(aft_genes) != 3:
			#	continue
			for i in range(0, 3-len(bef_genes)):
				bef_genes.append(['', -1, ''])
			for i in range(0, 3-len(aft_genes)):
				aft_genes.append(['', -1, ''])
			bef_genes = sorted(bef_genes)
			aft_genes = sorted(aft_genes)
			left = []
			right = []
			l, m, r = bef_genes
			if m[0] == r[0]:
				if r[1] - m[1] <= 50:
					left = r
			elif l[0] == r[0]:
				if r[1] - l[0] <= 50:
					left = r
			elif l[0] == m[0]:
				if m[1] - l[1] <= 50:
					left = m
			
			l, m, r = aft_genes
			if m[0] == r[0]:
				if r[1] - m[1] <= 50:
					right = m
			elif l[0] == r[0]:
				if r[1] - l[0] <= 50:
					right = l
			elif l[0] == m[0]:
				if m[1] - l[1] <= 50:
					right = l
			if left == [] or right == [] or left[0] == '' or right[0] == '':
				continue
			l = q_gff_db[left[2]][2]
			r = q_gff_db[right[2]][1]
			pos_db.append([[chrn, r_pos[0], r_pos[1]], [left[0], l, r]])
	
	with open(out_reg, 'w') as fout:
		for ref, qry in pos_db:
			if qry[0] not in q_bed_db:
				continue
			qs, qe = qry[1:]
			is_match = False
			for s, e in q_bed_db[qry[0]]:
				if min(qe, e)-max(qs, s) > 0:
					is_match = True
					break
			if is_match:
				if len(str(ref[0])) < 3:
					fout.write("Chr%s\t%d\t%d\t"%(ref[0], ref[1], ref[2]))
				else:
					fout.write("%s\t%d\t%d\t"%(ref[0], ref[1], ref[2]))
				if len(str(qry[0])) < 3:
					fout.write("Chr%d\t%d\t%d\n"%(qry[0], s, e))
				else:
					fout.write("%s\t%d\t%d\n"%(qry[0], s, e))
				

if __name__ == "__main__":
	if len(sys.argv) < 7:
		print("Usage: python "+sys.argv[0]+" <ref_gff> <query_gff> <ref_bed> <query_bed> <collinearity> <out_regions>")
	else:
		r_gff, q_gff, r_bed, q_bed, coll, out_reg = sys.argv[1:]
		get_common_fragments(r_gff, q_gff, r_bed, q_bed, coll, out_reg)


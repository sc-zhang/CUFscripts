# PapayaSVscripts
## Step 1
Use following scripts to remove overlap regions and generate non-overlaping bed file and junction bed file.
```
filter_blast_result.py <in_blast> <out_bed>
# The bed file contain 3 columns
# Chr Startpos Endpos

convert_bed_for_junction.py <in_bed> <junction_size> <out_bed>
# Extending regions with junction size
```

## Step 2
1. Mapping reads to reference genome with "minimap2".
2. calculate reads coverage with junction bed file using "bedtools coverage".

## Step 3
Filter coverage results.
```
extract_result.py <in_result> <threshold>
```

## Step 4
Plot figure with following scripts.
```
gen_data.py <in_beds> <in_blast> <threshold> <out_bed>
# in_beds are bed files generate with step 3, and can be split by comma like: 1.bed,2.bed

plot_data_V8.py <in_data> <chr_len> <out_pic>
```
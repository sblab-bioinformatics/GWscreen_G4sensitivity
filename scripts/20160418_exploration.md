

This script is about the HT1080 cell line (Pool 8)






### edgeR tutorial

Trying http://bioinf.wehi.edu.au/shRNAseq/pooledScreenAnalysis.pdf
Dai2014
http://bioinf.wehi.edu.au/shRNAseq/




### Data preparation

In sblab-srv001:

```bash
cd /media/solexa3/ProcessedRuns/160418_NS500222_0156_H2F33BGXY/fastq
ls *.fq.gz | wc -l # 24
rsync --progress *.fq.gz martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq
```

In uk-cri-lcst01:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq
zcat SLX-11624.H2F33BGXY.s_1.r_1.lostreads.fq.gz | grep -A 3 ':GGTAGC' | grep -v '\-\-' | gzip > SLX-11624.RPI24.H2F33BGXY.s_1.r_1.fq.gz
bsub -R "rusage[mem=8192]" -q interactive -Is bash # to go to interactive node
```




### Prepare hairpin, sample and fastq files for processAmplicons function in edgeR


```python


# hairpin file - 22nt and 21nt

with open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3.csv", "rb") as f:
	data = f.readlines()

len(data) # 113003

of=open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_22nt.txt", "w") # containing all hairpins with 22nt
of2=open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_22nt_pool8.txt", "w") # containing pool#8 hairpins with 22nt
of3=open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_21nt.txt", "w") # containing all hairpins with 21nt
of4=open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_21nt_pool8.txt", "w") # containing pool#8 hairpins with 21nt


of.write("ID\tSequences\tGene_name\tGene_id\tPool\n")
of2.write("ID\tSequences\tGene_name\tGene_id\n")
of3.write("ID\tSequences\tGene_name\tGene_id\tPool\n")
of4.write("ID\tSequences\tGene_name\tGene_id\n")


for hp in data[1:]:
	fields=hp.split(",")
	id=fields[0]
	genename=fields[1]
	geneid=fields[2]
	antisense_22=fields[6][59:81]
	antisense_21=fields[6][60:81]
	pool="".join(fields[8].split())
	of.write("%s\t%s\t%s\t%s\t%s\n" % (id, antisense_22, genename, geneid, pool))
	of3.write("%s\t%s\t%s\t%s\t%s\n" % (id, antisense_21, genename, geneid, pool))
	if pool=="Pool#8":
		of2.write("%s\t%s\t%s\t%s\n" % (id, antisense_22, genename, geneid))
		of4.write("%s\t%s\t%s\t%s\n" % (id, antisense_21, genename, geneid))


of.close()
of2.close()
of3.close()
of4.close()



# hairpin files - 22nt and 21nt - no duplicates in shRNA sequences and ids and no hashes in shRNA ids and pools

with open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_22nt.txt", "rb") as f:
	data = f.readlines()

len(data[1:]) # 113002

data_dict={}

for line in data[1:]:
	fields = line.split()
	id = fields[0]
	if id == "#N/A":
		id = "ULTRA-NA"
	seq = fields[1]
	genename = fields[2]
	geneid = fields[3]
	pool = fields[4].replace("#","")
	if seq in data_dict:
		data_dict[seq].append((id, genename, geneid, pool))
	else:
		data_dict[seq] = [(id, genename, geneid, pool)]


len(data_dict.keys()) # 109180

for seq in data_dict:
	if len(data_dict[seq]) > 20:
		print seq

data_dict["TCCACTTCATCACTGAACTGCT"] # this is just an example of a sequence associated to multiple genes



with open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_21nt.txt", "rb") as f:
	data_21 = f.readlines()

len(data_21[1:]) # 113002

data_dict_21={}

for line in data_21[1:]:
	fields = line.split()
	id = fields[0]
	if id == "#N/A":
		id = "ULTRA-NA"
	seq = fields[1]
	genename = fields[2]
	geneid = fields[3]
	pool = fields[4].replace("#","")
	if seq in data_dict_21:
		data_dict_21[seq].append((id, genename, geneid, pool))
	else:
		data_dict_21[seq] = [(id, genename, geneid, pool)]


len(data_dict_21.keys()) # 108465

for seq in data_dict_21:
	if len(data_dict_21[seq]) > 20:
		print seq

data_dict_21["CAGGACCATCTTCACACTCAC"] # this is just an example of a sequence associated to multiple genes


of5=open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_22nt_nodup.txt", "w") # containing all hairpins with 22nt with no duplicate sequences, no duplicate ids and no #s
of6=open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_22nt_pool8_nodup.txt", "w") # containing pool8 hairpins with 22nt with no duplicate sequences, no duplicate ids and no #s

of5.write("ID\tSequences\tshRNA_name\tGene_name\tGene_id\tPool\n")
of6.write("ID\tSequences\tshRNA_name\tGene_name\tGene_id\n")

count = 0
count_pool8 = 0
for seq in data_dict:
	count += 1
	fields = data_dict[seq][0]
	id = "hp%s" % str(count) # to create unique hairpin ids
	shrnaname = fields[0]
	genename = fields[1]
	geneid = fields[2]
	pool = fields[3]
	of5.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (id, seq, shrnaname, genename, geneid, pool))
	if pool == "Pool8":
		count_pool8 += 1
		id_pool8 = "hp%s" % str(count_pool8)
		of6.write("%s\t%s\t%s\t%s\t%s\n" % (id_pool8, seq, shrnaname, genename, geneid))


of5.close()
of6.close()



of7=open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_21nt_nodup.txt", "w") # containing all hairpins with 21nt with no duplicate sequences, no duplicate ids and no #s
of8=open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_21nt_pool8_nodup.txt", "w") # containing pool8 hairpins with 21nt with no duplicate sequences, no duplicate ids and no #s

of7.write("ID\tSequences\tshRNA_name\tGene_name\tGene_id\tPool\n")
of8.write("ID\tSequences\tshRNA_name\tGene_name\tGene_id\n")

count = 0
count_pool8 = 0
for seq in data_dict_21:
	count += 1
	fields = data_dict_21[seq][0]
	id = "hp%s" % str(count) # to create unique hairpin ids
	shrnaname = fields[0]
	genename = fields[1]
	geneid = fields[2]
	pool = fields[3]
	of7.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (id, seq, shrnaname, genename, geneid, pool))
	if pool == "Pool8":
		count_pool8 += 1
		id_pool8 = "hp%s" % str(count_pool8)
		of8.write("%s\t%s\t%s\t%s\t%s\n" % (id_pool8, seq, shrnaname, genename, geneid))


of7.close()
of8.close()



# sample files - created manually :
# /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Samples.txt
# /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Samples_t0.txt




# fastq file


# The idea here is to copy the index from the header to the beginning of the sequence in each fastq file. Then merge all 24 modified fastq files.
# The fourth quality line also needs modification
# Need to read file in chunks of four lines
# Check using biopython for shiqing's duplex tags

import os
import gzip

d='/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq'

for file in os.listdir(d):
	if "lostreads" not in file:
		print file
		fqf = gzip.open("%s/%s" % (d, file), "rb")
		fq = fqf.read()
		fqf.close()
		newchunks=[]
		for chunk in fq.split('@')[1:]:
			lines=chunk.split('\n')
			index=lines[0].split()[1].split(':')[3]
			newseq=index+lines[1]
			newqua='AAAAAA'+lines[3]
			newchunk='\n'.join([lines[0], newseq, lines[2], newqua])
			newchunks.append(newchunk)
		of = open("%s/%s" % (d, ".".join(file.split(".")[:-2])+".mod.fq"), "w")
		of.write('@'+'\n@'.join(newchunks)+'\n')
		of.close()

```

Merging *.mod.fq files and deleting individual files.


```bash
cat *.mod.fq > SLX-11624.H2F33BGXY.s_1.r_1.mod.merged.fq # merging all 24 modified fastq files
cat SLX-11624.H2F33BGXY.s_1.r_1.mod.merged.fq | paste -d " "  - - - - | wc -l # 192032042
cat SLX-11624.RPI01.H2F33BGXY.s_1.r_1.mod.fq SLX-11624.RPI02.H2F33BGXY.s_1.r_1.mod.fq SLX-11624.RPI03.H2F33BGXY.s_1.r_1.mod.fq > SLX-11624.RPI01-03.H2F33BGXY.s_1.r_1.mod.fq # merging only 3 t0 modified fastq files
cat SLX-11624.RPI01-03.H2F33BGXY.s_1.r_1.mod.fq | paste -d " "  - - - - | wc -l # 24637385
#rm *.mod.fq
```



### Explore representation of shRNAs from Pool #8 in t0

```python



hpfile = "Hairpins_21nt_pool8.txt"

d='/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420'

hpf = open("%s/%s" % (d, hpfile), "r")
hps = hpf.readlines()
hpf.close()

len(hps[1:]) # 9600

shrna_pool8_ids = []
shrna_pool8_seqs = []
shrna_pool8_genes = []

for hp in hps[1:]:
	fields = hp.split()
	id = fields[0]
	seq = fields[1]
	gene = fields[2]
	shrna_pool8_ids.append(id)
	shrna_pool8_seqs.append(seq)
	shrna_pool8_genes.append(gene)


len(set(shrna_pool8_ids)) # 9600
len(set(shrna_pool8_seqs)) # 9415
len(set(shrna_pool8_genes)) # 7823

from collections import Counter

shrna_pool8_seqs_cnt = Counter()
for seq in shrna_pool8_seqs:
	shrna_pool8_seqs_cnt[seq] += 1

len(shrna_pool8_seqs_cnt) # 9415
shrna_pool8_seqs_cnt.most_common(10) # these are shRNAs targeting different versions of the same gene, by "grep 'CTCTATGAGGATGTGGTCGTT' Hairpins_21nt_pool8.txt" we can find out the genes that are been targeted

shrna_pool8_genes_cnt = Counter()
for gene in shrna_pool8_genes:
	shrna_pool8_genes_cnt[gene] += 1

len(shrna_pool8_genes_cnt) # 7823
shrna_pool8_genes_cnt.most_common(10) # these are genes targeted by different shRNA sequences, by doing e.g. "grep 'RAB2A' Hairpins_21nt_pool8.txt" we can find out the shRNAs used to target each gene
# BRCA2 is another example of a gene targeted by two shRNAs





import gzip

file_rep1 = "SLX-11624.RPI01.H2F33BGXY.s_1.r_1.fq.gz"

fqf_rep1 = gzip.open("%s/fastq/%s" % (d, file_rep1), "rb")
fq_rep1 = fqf_rep1.read()
fqf_rep1.close()

seqs_rep1 = []

for chunk in fq_rep1.split('@')[1:]:
	lines=chunk.split('\n')
	seq = lines[1][1:22]
	seqs_rep1.append(seq)

len(seqs_rep1) # 9716031

seqs_rep1_cnt = Counter()
for seq in seqs_rep1:
	seqs_rep1_cnt[seq] += 1

len(seqs_rep1_cnt) # 267400


intersection_rep1_pool8 = set(seqs_rep1_cnt.keys()).intersection(set(shrna_pool8_seqs_cnt.keys()))
len(intersection_rep1_pool8) # 9385
100*float(9385)/9415 # 99.7%

count_rep1_pool8 = 0
for seq in intersection_rep1_pool8:
	count_rep1_pool8 = count_rep1_pool8 + seqs_rep1_cnt[seq]

count_rep1_pool8 # 7399686
100*float(7399686)/9716031 # 76.1%



file_rep2 = "SLX-11624.RPI02.H2F33BGXY.s_1.r_1.fq.gz"

fqf_rep2 = gzip.open("%s/fastq/%s" % (d, file_rep2), "rb")
fq_rep2 = fqf_rep2.read()
fqf_rep2.close()

seqs_rep2 = []

for chunk in fq_rep2.split('@')[1:]:
	lines=chunk.split('\n')
	seq = lines[1][1:22]
	seqs_rep2.append(seq)


len(seqs_rep2) # 8506386

seqs_rep2_cnt = Counter()
for seq in seqs_rep2:
	seqs_rep2_cnt[seq] += 1


len(seqs_rep2_cnt) # 247622


intersection_rep2_pool8 = set(seqs_rep2_cnt.keys()).intersection(set(shrna_pool8_seqs_cnt.keys()))
len(intersection_rep2_pool8) # 9385
100*float(9385)/9415 # 99.7%

count_rep2_pool8 = 0
for seq in intersection_rep2_pool8:
	count_rep2_pool8 = count_rep2_pool8 + seqs_rep2_cnt[seq]

count_rep2_pool8 # 6263934
100*float(6263934)/8506386 # 73.6%




file_rep3 = "SLX-11624.RPI03.H2F33BGXY.s_1.r_1.fq.gz"

fqf_rep3 = gzip.open("%s/fastq/%s" % (d, file_rep3), "rb")
fq_rep3 = fqf_rep3.read()
fqf_rep3.close()

seqs_rep3 = []

for chunk in fq_rep3.split('@')[1:]:
	lines=chunk.split('\n')
	seq = lines[1][1:22]
	seqs_rep3.append(seq)


len(seqs_rep3) # 6414968

seqs_rep3_cnt = Counter()
for seq in seqs_rep3:
	seqs_rep3_cnt[seq] += 1


len(seqs_rep3_cnt) # 228130


intersection_rep3_pool8 = set(seqs_rep3_cnt.keys()).intersection(set(shrna_pool8_seqs_cnt.keys()))
len(intersection_rep3_pool8) # 9382
100*float(9382)/9415 # 99.6

count_rep3_pool8 = 0
for seq in intersection_rep3_pool8:
	count_rep3_pool8 = count_rep3_pool8 + seqs_rep3_cnt[seq]

count_rep3_pool8 # 4812815
100*float(4812815)/6414968 # 75.0%


```




### Continuation ... edgeR tutorial



```R
library(edgeR)

# Read in sample & hairpin information
samples <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Samples.txt", sep = "\t", header = TRUE)
samples
hairpins <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_22nt.txt", sep = "\t", header = TRUE, comment.char = "")
hairpins[1:10,]

# Process raw sequences from fastq file
z = processAmplicons("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq/SLX-11624.H2F33BGXY.s_1.r_1.mod.merged.fq", barcodefile = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Samples.txt", hairpinfile = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_22nt.txt", barcodeStart=1, barcodeEnd=6, hairpinStart = 7, hairpinEnd = 28, allowMismatch=TRUE, barcodeMismatchBase=0, hairpinMismatchBase=2, verbose = TRUE)
# Error in processAmplicons("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq/SLX-11624.H2F33BGXY.s_1.r_1.mod.merged.fq",  : There are duplicate hairpin sequences.
# Error in processAmplicons("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq/SLX-11624.H2F33BGXY.s_1.r_1.mod.merged.fq",  : There are duplicate hairpin IDs.

# Read file with no duplicate hairpin sequences
hairpins_nodup <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_22nt_nodup.txt", sep = "\t", header = TRUE, comment.char = "")
hairpins_nodup[1:10,]

# Process raw sequences from fastq file
z = processAmplicons("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq/SLX-11624.H2F33BGXY.s_1.r_1.mod.merged.fq", barcodefile = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Samples.txt", hairpinfile = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_22nt_nodup.txt", barcodeStart=1, barcodeEnd=6, hairpinStart = 7, hairpinEnd = 28, allowMismatch=TRUE, barcodeMismatchBase=0, hairpinMismatchBase=2, verbose = TRUE)
# Error in scan(file, what, nmax, sep, dec, quote, skip, nlines, na.strings,  : line 25570 did not have 6 elements


# This is because line 25570 in Hairpins_22nt_nodup.txt is:
hairpins_nodup[25570,]
#            ID              Sequences shRNA_name Gene_name Gene_id   Pool
# 25570 hp25570 TGAGGCAAGAGCGTCCGCATGC       #N/A    09-Mar   92979 Pool#8
# Possibly the hairpinfile function in processAmplicons is not using comment.char = "", 
# I decided to rewrite Hairpins_22nt_nodup.txt with no hashes in shRNA_name or Pool:

# Read file with no #s
hairpins_nodup <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_22nt_nodup.txt", sep = "\t", header = TRUE, comment.char = "")
hairpins_nodup[1:10,]
hairpins_nodup[25570,]

# Process raw sequences from fastq file - 22nt
z <- processAmplicons("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq/SLX-11624.H2F33BGXY.s_1.r_1.mod.merged.fq", barcodefile = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Samples.txt", hairpinfile = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_22nt_nodup.txt", barcodeStart=1, barcodeEnd=6, hairpinStart = 7, hairpinEnd = 28, allowMismatch=TRUE, barcodeMismatchBase=0, hairpinMismatchBase=2, verbose = TRUE)
# It works but it takes long to process

# Process raw sequences from fastq file - 21nt
z <- processAmplicons("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq/SLX-11624.H2F33BGXY.s_1.r_1.mod.merged.fq", barcodefile = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Samples.txt", hairpinfile = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_21nt_nodup.txt", barcodeStart=1, barcodeEnd=6, hairpinStart = 8, hairpinEnd = 28, allowMismatch=FALSE, verbose = TRUE)

# -- Number of Barcodes : 24
# -- Number of Hairpins : 108465
#Processing reads in /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq/SLX-11624.H2F33BGXY.s_1.r_1.mod.merged.fq.
# -- Processing 10 million reads
# -- Processing 20 million reads
# -- Processing 30 million reads
# -- Processing 40 million reads
# -- Processing 50 million reads
# -- Processing 60 million reads
# -- Processing 70 million reads
# -- Processing 80 million reads
# -- Processing 90 million reads
# -- Processing 100 million reads
# -- Processing 110 million reads
# -- Processing 120 million reads
# -- Processing 130 million reads
# -- Processing 140 million reads
# -- Processing 150 million reads
# -- Processing 160 million reads
# -- Processing 170 million reads
# -- Processing 180 million reads
# -- Processing 190 million reads
# -- Processing 200 million reads
#Number of reads in file /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq/SLX-11624.H2F33BGXY.s_1.r_1.mod.merged.fq : 192032042
#
#The input run parameters are: 
# -- Barcode: start position 1	 end position 6	 length 6
# -- Hairpin: start position 8	 end position 28	 length 21
# -- Hairpin sequences need to match at specified positions. 
# -- Mismatch in barcode/hairpin sequences not allowed. 
#
#Total number of read is 192032042 
#There are 187600304 reads (97.6922 percent) with barcode matches
#There are 147475762 reads (76.7975 percent) with hairpin matches
#There are 144300154 reads (75.1438 percent) with both barcode and hairpin matches





# Read in t0 samples & pool8 hairpin information - 22nt
library(edgeR)
samples_t0 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Samples_t0.txt", sep = "\t", header = TRUE)
samples_t0

# Read file with no duplicate hairpin sequences and ids and no #s
hairpins_nodup <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_22nt_pool8_nodup.txt", sep = "\t", header = TRUE, comment.char = "")
hairpins_nodup[1:10,]

# Process raw sequences from fastq file
z <- processAmplicons("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq/SLX-11624.RPI01-03.H2F33BGXY.s_1.r_1.mod.fq", barcodefile = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Samples_t0.txt", hairpinfile = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_22nt_pool8_nodup.txt", barcodeStart=1, barcodeEnd=6, hairpinStart = 7, hairpinEnd = 28, allowMismatch=TRUE, barcodeMismatchBase=0, hairpinMismatchBase=2, verbose = TRUE)
# It works but it also takes long to process.
# In the end it finishes:

# -- Number of Barcodes : 3
# -- Number of Hairpins : 9416
#Processing reads in /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq/SLX-11624.RPI01-03.H2F33BGXY.s_1.r_1.mod.fq.
# -- Processing 10 million reads
# -- Processing 20 million reads
# -- Processing 30 million reads
#Number of reads in file /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq/SLX-11624.RPI01-03.H2F33BGXY.s_1.r_1.mod.fq : 24637385
#
#The input run parameters are: 
# -- Barcode: start position 1	 end position 6	 length 6
# -- Hairpin: start position 7	 end position 28	 length 22
# -- Hairpin sequences need to match at specified positions. 
# -- Allow sequence mismatch, <= 0 base in barcode sequence and <= 2 base in hairpin sequence. 
#
#Total number of read is 24637385 
#There are 24084720 reads (97.7568 percent) with barcode matches
#There are 20841363 reads (84.5924 percent) with hairpin matches
#There are 20403790 reads (82.8164 percent) with both barcode and hairpin matches

#  Note that a very high proportion of the reads (82.8%) match to expected combinations from our screen, which is an indication that  the  sequencing  for  this  screen  has  gone  well.


# Read in t0 samples & pool8 hairpin information - 21nt
library(edgeR)
samples_t0 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Samples_t0.txt", sep = "\t", header = TRUE)
samples_t0

# Read file with no duplicate hairpin sequences and ids and no #s
hairpins_nodup <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_21nt_pool8_nodup.txt", sep = "\t", header = TRUE, comment.char = "")
hairpins_nodup[1:10,]

# Process raw sequences from fastq file - no mismatch
z <- processAmplicons("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq/SLX-11624.RPI01-03.H2F33BGXY.s_1.r_1.mod.fq", barcodefile = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Samples_t0.txt", hairpinfile = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_21nt_pool8_nodup.txt", barcodeStart=1, barcodeEnd=6, hairpinStart = 8, hairpinEnd = 28, allowMismatch=FALSE, verbose = TRUE)

# -- Number of Barcodes : 3
# -- Number of Hairpins : 9333
#Processing reads in /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq/SLX-11624.RPI01-03.H2F33BGXY.s_1.r_1.mod.fq.
# -- Processing 10 million reads
# -- Processing 20 million reads
# -- Processing 30 million reads
#Number of reads in file /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq/SLX-11624.RPI01-03.H2F33BGXY.s_1.r_1.mod.fq : 24637385
#
#The input run parameters are: 
# -- Barcode: start position 1	 end position 6	 length 6
# -- Hairpin: start position 8	 end position 28	 length 21
# -- Hairpin sequences need to match at specified positions. 
# -- Mismatch in barcode/hairpin sequences not allowed. 
#
#Total number of read is 24637385 
#There are 24084720 reads (97.7568 percent) with barcode matches
#There are 18322445 reads (74.3685 percent) with hairpin matches
#There are 17966903 reads (72.9254 percent) with both barcode and hairpin matches

# Process raw sequences from fastq file - 1 mismatch
z <- processAmplicons("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq/SLX-11624.RPI01-03.H2F33BGXY.s_1.r_1.mod.fq", barcodefile = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Samples_t0.txt", hairpinfile = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hairpins_21nt_pool8_nodup.txt", barcodeStart=1, barcodeEnd=6, hairpinStart = 8, hairpinEnd = 28, allowMismatch=TRUE, barcodeMismatchBase=0, hairpinMismatchBase=1, verbose = TRUE)

# -- Number of Barcodes : 3
# -- Number of Hairpins : 9333
#Processing reads in /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq/SLX-11624.RPI01-03.H2F33BGXY.s_1.r_1.mod.fq.
# -- Processing 10 million reads
# -- Processing 20 million reads
# -- Processing 30 million reads
#Number of reads in file /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq/SLX-11624.RPI01-03.H2F33BGXY.s_1.r_1.mod.fq : 24637385
#
#The input run parameters are: 
# -- Barcode: start position 1	 end position 6	 length 6
# -- Hairpin: start position 8	 end position 28	 length 21
# -- Hairpin sequences need to match at specified positions. 
# -- Allow sequence mismatch, <= 0 base in barcode sequence and <= 1 base in hairpin sequence. 
#
#Total number of read is 24637385 
#There are 24084720 reads (97.7568 percent) with barcode matches
#There are 20617143 reads (83.6824 percent) with hairpin matches
#There are 20185708 reads (81.9312 percent) with both barcode and hairpin matches


# Clearly, if we allow for 1 mismatch in the 21nt, it still takes really long








# Continue with the data loaded of all samples, all pools, 21nt, 0 mismatches 
# Note that a high proportion of the reads (75.1%) match to expected combinations from our screen, which is an indication that  the  sequencing  for  this  screen  has  gone  ok.

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160427_readcounts.pdf")
par(mar=c(7, 6, 4, 2) + 0.1, mgp = c(4, 1, 0))
barplot(colSums(z$counts), las=2, ylab="Reads with both barcode and 21nt hairpin matches", ylim=c(0,10000000), names.arg = z$samples$group)
dev.off()


# Plot pool8 hairpin totals across all samples
plateinfo <- z$genes$Pool
selhp <- plateinfo == "Pool8"
barplot(rowSums(z$counts[selhp,]), las=2, main="Counts per hairpin")
a <- data.table(z$counts[selhp,])
a[,total := rowSums(.SD)]
a[,id := rownames(z$counts[selhp,])]
a[total > 2000000]
a[total > 500000]
# "hp23696" particularly abundant in 20_HT1080_DMSO_12_41_pd_rep2
data.table(z$genes)[ID == "hp23696" | ID == "hp54058"]
# Where it accounts for 56% of the reads:
colSums(z$counts)


#        ID             Sequences    shRNA_name Gene_name Gene_id  Pool
#1: hp23696 TTATTTCCATGGCTCATCGCC ULTRA-3306733   ZC3H12D  340152 Pool8
#2: hp54058 TTAATTGCTCCATCACAGTTG ULTRA-3316623     ITGAM    3684 Pool8
# hp23696 ZC3H12D http://www.uniprot.org/uniprot/A2A288
# hp54058 ITGAM http://www.uniprot.org/uniprot/P11215 


# Look at CPMs of BRCA2
# grep "BRCA2" data/20160420/Hairpins_21nt_nodup.txt
# hp12474	TTTCCAACTGGATCTGAGCTT	ULTRA-3373304	BRCA2	675	Pool8
# hp42275	AAGTTATGAGAATTTCTACTG	ULTRA-3373301	BRCA2	675	Pool8
# hp49418	AAATTATCACTTAAGAGCTTA	ULTRA-3373302	BRCA2	675	Pool4

z$counts[c("hp12474", "hp42275", "hp49418"),]
colSums(z$counts)
cpm(z$counts, log = FALSE)[c("hp12474", "hp42275"),]

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160427_cpm_BRCA2.pdf")
par(mar=c(7, 4, 4, 2) + 0.1)
barplot(cpm(z$counts, log = FALSE)[c("hp12474", "hp42275"),], las=2, ylab="Counts per million", names.arg = z$samples$group, beside = TRUE, legend.text = TRUE, main = "BRCA2", ylim = c(0, 600))
dev.off()

# Look at CPMs of PARP1
# grep "PARP1" data/20160420/Hairpins_21nt_nodup.txt | grep "Pool8"
# hp60533	AGAGAGTGAGCAGCTCGTCGG	ULTRA-3406334	PARP10	84875	Pool8
# hp61947	ATGTTAGCTGAAGGATCAGGG	ULTRA-3216281	PARP1	142	Pool8
# hp66301	TATCTTGAGTGAAGATGTGCT	ULTRA-3342668	PARP14	54625	Pool8
# hp87927	ATTTCAAACAGGAAGTCCGGT	ULTRA-3344146	PARP16	54956	Pool8
# hp101012	AAAGAGGAGGCACTCTGCTGG	ULTRA-3365352	PARP12	64761	Pool8

barplot(cpm(z$counts, log = FALSE)["hp61947",], las=2, ylab="Counts per million", names.arg = z$samples$group)

# Look at CPMs of HNRNPA1
# grep "HNRNPA1" Hairpins_21nt_nodup.txt | grep "Pool8"
# hp13947	TCCCAAAATCATTGTAGCTTC	ULTRA-3302579	HNRNPA1	3178	Pool8

barplot(cpm(z$counts, log = FALSE)["hp13947",], las=2, ylab="Counts per million", names.arg = z$samples$group)


# Make a MDS plot to visualise relationships between replicate samples
seltf1r <- plateinfo=="Pool8"
seltf1c <- z$samples$ID != "19_HT1080_DMSO_12_41_pd_rep1"
z1 <- z[seltf1r, seltf1c]
plotMDS(z1, labels = z1$samples$group, col = c(rep(1,3), rep(2,9), rep(3,6), rep(4,5))) # coloured according to time point
legend("topright", legend=c("t0", "t1", "t2", "t3"), col=1:4, pch=15)
plotMDS(z1, labels = z1$samples$group, col = c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(2,3), rep(4,3), rep(2,2), rep(4,3))) # coloured according to drug used
legend("topright", legend=c("t0", "DMSO", "PhenDC3", "PDS"), col=1:4, pch=15)


# Begin differential representation analysis

# "t2_DMSO" vs. "t2_PDS" : comparing DMSO and PDS treatments at t2
# Estimate dispersions
seltf1r <- plateinfo == "Pool8"
seltf1c <- z$samples$group == "t2_DMSO" | z$samples$group == "t2_PDS"
z2 <- z[seltf1r, seltf1c]
z2 <- estimateDisp(z2)
sqrt(z2$common.dispersion) # 59% - the variation between replicates samples is large. This means that we might not be able to detect differences

# Assess differential representation between DMSO and PDS treatments at t2
# using classic exact testing methodology in edgeR
de.t2.DMSOvsPDS <- exactTest(z2, pair=c("t2_DMSO", "t2_PDS")) # Robinson2008

# Show top ranked hairpins
topTags(de.t2.DMSOvsPDS)
# There are no hairpins with FDR < 0.05
# Despite not being significant, just see LENEP (hp31273) to check:
barplot(cpm(z$counts, log = FALSE)["hp31273",], las=2, ylab="Counts per million", names.arg = z$samples$group)
barplot(cpm(z2$counts, log = FALSE)["hp31273",], las=2, ylab="Counts per million", names.arg = z2$samples$group)

# Plot logFC versus logCPM
plotSmear(de.t2.DMSOvsPDS, ylim = c(-10, 10))
abline(h = c(-5, 0, 5), col = c("dodgerblue", "yellow", "dodgerblue"), lty=2)

# Make a MDS plot
plotMDS(z2, labels = z2$samples$group, col = c(rep(1,3), rep(2,3)))
legend("topright", legend=c("DMSO", "PDS"), col=1:2, pch=15)

# Compare t0 against the rest of time points and conditions as suggested by Darcie
pools <- z$genes$Pool
seltf1r <- pools == "Pool8"
treatments <- levels(z$samples$group)[-1] # everything except t0
for (t in treatments){
	print(c("t0", t))
	seltf1c <- z$samples$group == "t0" | z$samples$group == t
	z2 <- estimateDisp(z[seltf1r, seltf1c])
	disp <- sqrt(z2$common.dispersion)
	de <- exactTest(z2, pair=c("t0", t))
	top <- topTags(de, n=Inf, p.value = 0.05)
	write.table(top, file = sprintf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160429_t0vs%s_pool8.txt", t), sep = "\t", row.names = FALSE)
	pdf(sprintf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160429_t0vs%s_pool8_mds.pdf", t))
	plotMDS(z2, labels = z2$samples$group, col = c(rep(1,3), rep(2,3)))
	legend("topright", legend=c("t0", t), col=1:2, pch=15)
	dev.off()
	thresh <- 0.05
	top2 <- topTags(de, n = Inf)
	top3 <- top2$table[top2$table$FDR < thresh, 1]
	if (length(top3) > 0) {
		pdf(sprintf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160429_t0vs%s_pool8_spear.pdf", t))
		plotSmear(de, de.tags = top3, ylim = c(-10, 10))
		abline(h = c(-5, 0, 5), col = c("dodgerblue", "yellow", "dodgerblue"), lty=2)
		dev.off()
	}
	else {
		pdf(sprintf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160429_t0vs%s_pool8_spear.pdf", t))
		plotSmear(de, de.tags = NULL, ylim = c(-10, 10))
		abline(h = c(-5, 0, 5), col = c("dodgerblue", "yellow", "dodgerblue"), lty=2)
		dev.off()
	}	
}

# Looks like there is an overall loss of shRNA upon treatment with PDS, see t2 and t3.
# I arranged all the above spear plots in three slides (20160429.odp), one for each time point: first slide t1 (DMSO, PDS, PhenDC3), second slide t2 (DMSO, PDS) and third slide t3 (DMSO, PDS)


```


Need to continue with the tutorial



Dario suggested to use:
- fastqc - to have a quick look at the quality of the data
- cutadapt - to remove illumina adapters
- bwa mem - to align
- edgeR/DESeq to work with the count matrices


Barbara has suggested to align to the genome to see if there are any off-targets of the shRNAs.









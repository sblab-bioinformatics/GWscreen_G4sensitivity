
This script aims to process the data of the mini screen performed by Darcie. It is writen on the basis of: [20160627_processing.md](20160627_processing.md), [20160902_processing.md](20160902_processing.md), [20160902_processing.md](20160902_processing.md), [20161108_processing.md](20161108_processing.md) and [20170105_processing.md](20170105_processing.md).


This is an email from Darcie on 20170330 summarising the priorities for this part:

I think this is what we decided for the analysis:

FDR 0.05; no FC threshold (I'll apply this myself) for -

1) comparing PDS tf v t0, PhenDC3 tf v t0, DMSO tf v t0. Tables for PDS
only; PhenDC3 only and the overlap, for the cell lines separately. Both
sensitisers and resistant

2) comparing PDS tf v t0, PhenDC3 tf v t0, DMSO tf v t0 and SA100128 tf
vs t0. Tables for PDS only; PhenDC3 only; SA28 only and the different
overlaps. Both sensitisers and resistant.

3) I will then look at the overlap between the two cell lines myself

4) At a lower priority, just a comparison between SA100128 and DMSO in
the two cell lines (again if time is an issue, I think I can do this
using the tables from (2)

I will aim to get the sequencing data to you (if all goes to plan) the
week beginning 10 th April. I think number 1) is the biggest priority in
terms of the paper, but number 2) also needs doing, to help with
Santosh's experiments.


There will be two sequencing runs (in principle). Both will contain A375 and HT1080 as follows:

- A375 (15 libraries x 2 runs)
  - t0 (x3)
  - t15 DMSO (x3)
  - t15 PDS (x3)
  - t15 PhenDC3 (x3)
  - t15 SA (x3)

- HT1080 (15 libraries x 2 runs)
  - t0 (x3)
  - t15 DMSO (x3)
  - t15 PDS (x3)
  - t15 PhenDC3 (x3)
  - t15 SA (x3)



## Preparing alignment reference

Copy excel sheet sent by Darcie on 20170330 email and also the sheet number 1, which is the relevant here. In `uk-cri-lcst01`:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data
mkdir -p 20170404/reference
cd 20170404/reference
rsync -arvuP martin03@nm149s012925:/Users/martin03/Desktop/20170120_Miniscreen_directory.xls .
rsync -arvuP martin03@nm149s012925:/Users/martin03/Desktop/20170120_Miniscreen_directory.csv .

```

Need to figure out which sequences I do need to extract from this file:

```python
with open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/reference/20170120_Miniscreen_directory.csv", "rb") as f:
	data = f.readlines()

len(data) # 8021 lines


#Playing with the data a bit to have an idea how best to generate the reference fasta file
gene_names=[]
gene_ids=[]
seqs=[]
cshl_ids=[]

for hp in data[1:]:
  fields=hp.split(",")
  gene_names.append(fields[0])
  gene_ids.append(fields[1])
  seqs.append(fields[4][59:81].upper())
  cshl_ids.append(fields[6])

len(set(gene_names)) # 1365 gene names, including RFP, EGFP, REN and OR*
len(set(gene_ids)) # 1360 gene ids
len(set(seqs)) # 7885 unique shRNA seqs
len(set(cshl_ids)) # 7988 gene ids

from collections import Counter

seqs_cnt = Counter()
for seq in seqs:
  seqs_cnt[seq] += 1

seqs_cnt.most_common(10)
#[('TTCTAGGAGAGGTTGCGCCTGC', 6), ('TCTGAATCCTGGACTCCGGGAG', 6), ('AAACCAAATCTGGACCCTGGGC', 6), ('GAAACCAGATCTGAATCCTGGA', 6), ('TAATCTGGACCCTGGGCTCCGG', 6), ('TCATTCTGAAACCAAATCTGGA', 6), ('TACAGAGAAAGGTCTGCGTCAT', 4), ('TAGGGTTTCTCTCCAGTGTGAA', 4), ('TATCAAAGGCATGAAGATCTAA', 4), ('TTTGATAAGGCCATCAGGTTTC', 4)]

#grep "TTCTAGGAGAGGTTGCGCCTGC" 20170120_Miniscreen_directory.csv
#DUX2,26583,hg19,10,TGCTGTTGACAGTGAGCGACAGGCGCAACCTCTCCTAGAATAGTGAAGCCACAGATGTATTCTAGGAGAGGTTGCGCCTGCTGCCTACTGCCTCGGA,204-0161_G10,DUX2__sensor__1,ULTRA-3269141,10,3,6
#DUX4L2,728410,hg19,10,TGCTGTTGACAGTGAGCGACAGGCGCAACCTCTCCTAGAATAGTGAAGCCACAGATGTATTCTAGGAGAGGTTGCGCCTGCTGCCTACTGCCTCGGA,204-0161_C10,DUX4L2__sensor__1,ULTRA-3385534,10,3,6
#DUX4L3,653548,hg19,10,TGCTGTTGACAGTGAGCGACAGGCGCAACCTCTCCTAGAATAGTGAAGCCACAGATGTATTCTAGGAGAGGTTGCGCCTGCTGCCTACTGCCTCGGA,204-0161_A10,DUX4L3__sensor__1,ULTRA-3367260,10,3,6
#DUX4L4,441056,hg19,10,TGCTGTTGACAGTGAGCGACAGGCGCAACCTCTCCTAGAATAGTGAAGCCACAGATGTATTCTAGGAGAGGTTGCGCCTGCTGCCTACTGCCTCGGA,204-0161_B10,DUX4L4__sensor__1,ULTRA-3326166,10,3,6
#DUX4L5,653545,hg19,10,TGCTGTTGACAGTGAGCGACAGGCGCAACCTCTCCTAGAATAGTGAAGCCACAGATGTATTCTAGGAGAGGTTGCGCCTGCTGCCTACTGCCTCGGA,204-0161_E10,DUX4L5__sensor__1,ULTRA-3367254,10,3,6
#DUX4L7,653543,hg19,10,TGCTGTTGACAGTGAGCGACAGGCGCAACCTCTCCTAGAATAGTGAAGCCACAGATGTATTCTAGGAGAGGTTGCGCCTGCTGCCTACTGCCTCGGA,204-0093_G12,DUX4L7__sensor__1,ULTRA-3367242,10,3,6

# They all belong to the same family - they are perhaps products of gene duplication



# Produce table of common shRNA sequences that appear repeatedly in several lanes

seqs_dict = {}

for hp in data[1:]:
  fields=hp.split(",")
  seq = fields[4][59:81].upper()
  cshl_id = fields[6]
  ultra_id = fields[7]
  i = "__".join([cshl_id, ultra_id])
  if seq in seqs_dict:
    seqs_dict[seq].append(i)
  else:
    seqs_dict[seq] = [i]

len(seqs_dict) # 7885, ok, as above

ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/reference/20170120_Miniscreen_directory_commonseqs.txt", "w")

for seq in seqs_cnt.most_common(100):
  if seq[1] > 1:
    ofile.write("%s\t%s\t%s\n" % (seq[0], int(len(seqs_dict[seq[0]])), ",".join(seqs_dict[seq[0]])))

ofile.close()


# See cases like KDM1A__shERWOOD__1U__Sensor__2__ULTRA-3482295, they same id appears multiple times with the same sequence, these are duplicates!

cshl_ids_cnt = Counter()
for id in cshl_ids:
  cshl_ids_cnt[id] += 1

cshl_ids_cnt.most_common(32)
#[('RPA3__sherwood_1u__5', 2), ('OR2A42__1__1', 2), ('OR2A1__3__3', 2), ('NSD1__sherwood_1u__5', 2), ('SUB1__sensor__2', 2), ('C9orf128__sensor__2', 2), ('OR2A42__3__3', 2), ('RECQL4__sensor__2', 2), ('RECQL4__sensor__1', 2), ('AHCYL1__sherwood__4__4', 2), ('SMARCD3__sherwood_1u__5', 2), ('PCNA__sherwood__2__2', 2), ('OR2A1__1__1', 2), ('BRD2__sherwood_1u__8', 2), ('KDM1A__shERWOOD__1U__Sensor__2', 2), ('KDM1A__shERWOOD__1U__Sensor__1', 2), ('SIN3A__sherwood_1u__8', 2), ('IMP5__sensor__2', 2), ('NAP1L1__sherwood__5__4', 2), ('RECQL4__sherwood_1u__2', 2), ('PRDM4__sherwood_1u__2', 2), ('MRPL23__sherwood__4__4', 2), ('RBM3__sherwood_1u__2__2', 2), ('ADH5__sherwood_1u__4', 2), ('SUB1__sensor__1', 2), ('RNF207__sherwood__2__2', 2), ('FGFR3__sherwood_1u__8', 2), ('MBD2__sherwood_1u__8', 2), ('BCL2L14__sherwood__1__1', 2), ('SUB1__6__4', 2), ('TOMM20__sherwood__3__1', 2), ('INO80B__sherwood__3__3', 2)]

# Sometimes the same CSHLID is associated to different sequences, we need to use both CSHLID and UltraId as identifiers for the shRNA sequence. Like with OR2A42__1__1 and SUB1__sensor__2, if CSHLID+UltraId is present more than once, then the lines are duplicated.

id_seq_dict = {}
id_score_dict = {}

for hp in data[1:]:
  fields=hp.split(",")
  score = fields[3]
  if score == "N/A": # for some reason some of the scores in the file are just this "Need score"
    score = "NA"
  seq = fields[4][59:81].upper()
  cshl_id = fields[6]
  ultra_id = fields[7]
  i = "__".join([cshl_id, ultra_id])
  id_seq_dict[i] = seq
  id_score_dict[i] = score

len(id_seq_dict) # 8006, there were some duplicates
len(id_score_dict) # 8006, fine

#Produce reference fasta file
ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/reference/20170120_Miniscreen_directory.fa", "w")

for hp in id_seq_dict:
  ofile.write(">%s\n%s\n" % (hp, id_seq_dict[hp]))

ofile.close()

#Produce reference scores file
ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/reference/20170120_Miniscreen_directory_scores.txt", "w")

for hp in id_score_dict:
  ofile.write("%s\t%s\n" % (hp, id_score_dict[hp]))

ofile.close()

```

Checking things are as expected:

```bash
grep "RPA3__sherwood_1u__5" 20170120_Miniscreen_directory.fa -A 1
#>RPA3__sherwood_1u__5__ULTRA-3359186
#TAAACTGGACATAAGATGTACA
#--
#>RPA3__sherwood_1u__5__ULTRA-3430608
#TCATGGATAATTTTCACAGCTT

grep "OR2A42__1__1" 20170120_Miniscreen_directory.fa -A 1
#>OR2A42__1__1__ULTRA-3322547
#TACAGAGAAAGGTCTGCGTCAT

grep "TTCTAGGAGAGGTTGCGCCTGC" 20170120_Miniscreen_directory.fa -B 1
#>DUX4L7__sensor__1__ULTRA-3367242
#TTCTAGGAGAGGTTGCGCCTGC
#--
#>DUX4L4__sensor__1__ULTRA-3326166
#TTCTAGGAGAGGTTGCGCCTGC
#--
#>DUX2__sensor__1__ULTRA-3269141
#TTCTAGGAGAGGTTGCGCCTGC
#--
#>DUX4L5__sensor__1__ULTRA-3367254
#TTCTAGGAGAGGTTGCGCCTGC
#--
#>DUX4L2__sensor__1__ULTRA-3385534
#TTCTAGGAGAGGTTGCGCCTGC
#--
#>DUX4L3__sensor__1__ULTRA-3367260
#TTCTAGGAGAGGTTGCGCCTGC

```

What does bowtie2 do when it finds the same reference sequence multiple times (but with different ids)? Does it randomly align to any of the two?




## Location of fastq files

Run id: 170405_NS500222_0289_HJTY3BGX2    SLX11932

It is in [Basepace](https://basespace.illumina.com/run/23807804/SLX11932_DM) too.

First, create receiving folders in `uk-cri-lcst01`:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404
mkdir fastq
```

Second, copy them from `sblab-srv001` to `uk-cri-lcst01`,

```bash
rsync --progress /media/staging/170405_NS500222_0289_HJTY3BGX2/fastq/*.fq.gz martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/fastq
rsync --progress /media/staging/170405_NS500222_0289_HJTY3BGX2/fastq/SLX-11932.HJTY3BGX2.s_1.barcodesummary.txt martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/fastq
rsync --progress /media/staging/170405_NS500222_0289_HJTY3BGX2/fastq/SLX-11932.HJTY3BGX2.s_1.contents.csv martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/fastq

```

Looking at SLX-11932.HJTY3BGX2.s_1.barcodesummary.txt and comparing them with the previous barcodesummary files. Here we have:

* 294280855 reads. Lost: 6.68%


## Quality check

In uk-cri-lcst01,

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/fastq
mkdir ../fastqc

for fq in SLX-11932.RPI*.fq.gz
do
bname=${fq%.HJTY3BGX2.s_1.r_1.fq.gz}
bsub -J $bname -oo ../fastqc/$bname.log -R "rusage[mem=4096]" "/home/martin03/sw/fastqc/fastqc --noextract --nogroup -q -o ../fastqc $fq && \
rsync ../fastqc/*.html martin03@nm149s012925:/home/cri.camres.org/martin03/Desktop/fastqc"
done

```

In `nm149s012925`, having a look at the reports:

```bash
cd /home/cri.camres.org/martin03/Desktop/fastqc
firefox *.html &
```

There is drop in quality of base 1 too. Rest of things look fine.


## Trimming

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/fastq
mkdir ../fastq_trim

for fq in SLX-11932.RPI*;
do
  bsub -R "rusage[mem=4096]" "zcat $fq | fastx_trimmer -Q33 -l22 | pigz > ../fastq_trim/${fq%.fq.gz}.trim.fq.gz"
done

```


## Alignment using `bowtie2`

### Generate fasta file containing hairpins

The file was already generated above. It is here:

`/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/reference/20170120_Miniscreen_directory.fa`

### Generate index for hairpins fasta file using bowtie2-build

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/reference
bowtie2-build 20170120_Miniscreen_directory.fa 20170120_Miniscreen_directory
```


### Align

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/fastq_trim
IND="/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/reference/20170120_Miniscreen_directory"
mkdir ../sam/

for fq in SLX-11932.RPI*;
do
bsub -J $fq -oo ../sam/${fq%.trim.fq.gz}.log -R "rusage[mem=4096]" "zcat $fq | bowtie2 --norc --no-hd -x $IND -U - -S ../sam/${fq%.trim.fq.gz}.sam"
done

cd ../sam


for log in SLX-11932.RPI*.log;
do
echo $log
grep "overall alignment rate" $log
done

# 81-93% of reads with at least one reported alignment. Lower than in previous rounds


for sam in *.sam
do
  echo $sam, `cat $sam | awk '$3 !="*"' | wc -l`
done

#SLX-11932.RPI01.HJTY3BGX2.s_1.r_1.sam, 6885755
#SLX-11932.RPI02.HJTY3BGX2.s_1.r_1.sam, 6790473
#SLX-11932.RPI03.HJTY3BGX2.s_1.r_1.sam, 7436086
#SLX-11932.RPI04.HJTY3BGX2.s_1.r_1.sam, 9568823
#SLX-11932.RPI05.HJTY3BGX2.s_1.r_1.sam, 8374428
#SLX-11932.RPI06.HJTY3BGX2.s_1.r_1.sam, 8013880
#SLX-11932.RPI07.HJTY3BGX2.s_1.r_1.sam, 7945567
#SLX-11932.RPI08.HJTY3BGX2.s_1.r_1.sam, 8147781
#SLX-11932.RPI09.HJTY3BGX2.s_1.r_1.sam, 9299761
#SLX-11932.RPI10.HJTY3BGX2.s_1.r_1.sam, 9481987
#SLX-11932.RPI11.HJTY3BGX2.s_1.r_1.sam, 8109324
#SLX-11932.RPI12.HJTY3BGX2.s_1.r_1.sam, 7753721
#SLX-11932.RPI13.HJTY3BGX2.s_1.r_1.sam, 7026772
#SLX-11932.RPI14.HJTY3BGX2.s_1.r_1.sam, 7186912
#SLX-11932.RPI15.HJTY3BGX2.s_1.r_1.sam, 7310312
#SLX-11932.RPI16.HJTY3BGX2.s_1.r_1.sam, 7662640
#SLX-11932.RPI17.HJTY3BGX2.s_1.r_1.sam, 8229337
#SLX-11932.RPI18.HJTY3BGX2.s_1.r_1.sam, 7021856
#SLX-11932.RPI19.HJTY3BGX2.s_1.r_1.sam, 7676106
#SLX-11932.RPI20.HJTY3BGX2.s_1.r_1.sam, 8248709
#SLX-11932.RPI21.HJTY3BGX2.s_1.r_1.sam, 8384807
#SLX-11932.RPI22.HJTY3BGX2.s_1.r_1.sam, 8226260
#SLX-11932.RPI23.HJTY3BGX2.s_1.r_1.sam, 9079780
#SLX-11932.RPI24.HJTY3BGX2.s_1.r_1.sam, 9036228
#SLX-11932.RPI25.HJTY3BGX2.s_1.r_1.sam, 10135631
#SLX-11932.RPI26.HJTY3BGX2.s_1.r_1.sam, 7436887
#SLX-11932.RPI27.HJTY3BGX2.s_1.r_1.sam, 7762096
#SLX-11932.RPI28.HJTY3BGX2.s_1.r_1.sam, 7509401
#SLX-11932.RPI29.HJTY3BGX2.s_1.r_1.sam, 8515472
#SLX-11932.RPI30.HJTY3BGX2.s_1.r_1.sam, 7994584


cat SLX-11932.RPI01.HJTY3BGX2.s_1.r_1.sam | cut -f2 | sort | uniq -c
#6885755 0
#1606552 4

```


## Build table of counts

Generating counts:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/sam
mkdir ../counts

for sam in SLX-11932.RPI*.sam;
do
bsub -J $sam -R "rusage[mem=4096]" "cat $sam | cut -f3 | sort | uniq -c > ../counts/${sam%.sam}.counts"
done

#e.g. DHX36 counts in t0 libraries
grep "DHX36" ../counts/SLX-11932.RPI0{1..6}.HJTY3BGX2.s_1.r_1.counts
#../counts/SLX-11932.RPI01.HJTY3BGX2.s_1.r_1.counts:    344 DHX36__sensor__1__ULTRA-3226137
#../counts/SLX-11932.RPI01.HJTY3BGX2.s_1.r_1.counts:    379 DHX36__sensor__2__ULTRA-3226132
#../counts/SLX-11932.RPI01.HJTY3BGX2.s_1.r_1.counts:   2788 DHX36__sherwood_1u__3__ULTRA-3226135
#../counts/SLX-11932.RPI01.HJTY3BGX2.s_1.r_1.counts:   1714 DHX36__sherwood_1u__4__ULTRA-3226134
#../counts/SLX-11932.RPI01.HJTY3BGX2.s_1.r_1.counts:    579 DHX36__sherwood_1u__5__ULTRA-3424948
#../counts/SLX-11932.RPI01.HJTY3BGX2.s_1.r_1.counts:    775 DHX36__sherwood_1u__7__ULTRA-3424950
#../counts/SLX-11932.RPI01.HJTY3BGX2.s_1.r_1.counts:   1081 DHX36__sherwood_1u__8__ULTRA-3424951
#../counts/SLX-11932.RPI02.HJTY3BGX2.s_1.r_1.counts:    432 DHX36__sensor__1__ULTRA-3226137
#../counts/SLX-11932.RPI02.HJTY3BGX2.s_1.r_1.counts:    414 DHX36__sensor__2__ULTRA-3226132
#../counts/SLX-11932.RPI02.HJTY3BGX2.s_1.r_1.counts:   2992 DHX36__sherwood_1u__3__ULTRA-3226135
#../counts/SLX-11932.RPI02.HJTY3BGX2.s_1.r_1.counts:   2135 DHX36__sherwood_1u__4__ULTRA-3226134
#../counts/SLX-11932.RPI02.HJTY3BGX2.s_1.r_1.counts:    816 DHX36__sherwood_1u__5__ULTRA-3424948
#../counts/SLX-11932.RPI02.HJTY3BGX2.s_1.r_1.counts:    771 DHX36__sherwood_1u__7__ULTRA-3424950
#../counts/SLX-11932.RPI02.HJTY3BGX2.s_1.r_1.counts:    865 DHX36__sherwood_1u__8__ULTRA-3424951
#../counts/SLX-11932.RPI03.HJTY3BGX2.s_1.r_1.counts:    491 DHX36__sensor__1__ULTRA-3226137
#../counts/SLX-11932.RPI03.HJTY3BGX2.s_1.r_1.counts:    617 DHX36__sensor__2__ULTRA-3226132
#../counts/SLX-11932.RPI03.HJTY3BGX2.s_1.r_1.counts:   3437 DHX36__sherwood_1u__3__ULTRA-3226135
#../counts/SLX-11932.RPI03.HJTY3BGX2.s_1.r_1.counts:   2261 DHX36__sherwood_1u__4__ULTRA-3226134
#../counts/SLX-11932.RPI03.HJTY3BGX2.s_1.r_1.counts:    680 DHX36__sherwood_1u__5__ULTRA-3424948
#../counts/SLX-11932.RPI03.HJTY3BGX2.s_1.r_1.counts:    870 DHX36__sherwood_1u__7__ULTRA-3424950
#../counts/SLX-11932.RPI03.HJTY3BGX2.s_1.r_1.counts:   1119 DHX36__sherwood_1u__8__ULTRA-3424951
#../counts/SLX-11932.RPI04.HJTY3BGX2.s_1.r_1.counts:    419 DHX36__sensor__1__ULTRA-3226137
#../counts/SLX-11932.RPI04.HJTY3BGX2.s_1.r_1.counts:    699 DHX36__sensor__2__ULTRA-3226132
#../counts/SLX-11932.RPI04.HJTY3BGX2.s_1.r_1.counts:   4958 DHX36__sherwood_1u__3__ULTRA-3226135
#../counts/SLX-11932.RPI04.HJTY3BGX2.s_1.r_1.counts:   4028 DHX36__sherwood_1u__4__ULTRA-3226134
#../counts/SLX-11932.RPI04.HJTY3BGX2.s_1.r_1.counts:   1841 DHX36__sherwood_1u__5__ULTRA-3424948
#../counts/SLX-11932.RPI04.HJTY3BGX2.s_1.r_1.counts:   1124 DHX36__sherwood_1u__7__ULTRA-3424950
#../counts/SLX-11932.RPI04.HJTY3BGX2.s_1.r_1.counts:   1515 DHX36__sherwood_1u__8__ULTRA-3424951
#../counts/SLX-11932.RPI05.HJTY3BGX2.s_1.r_1.counts:    468 DHX36__sensor__1__ULTRA-3226137
#../counts/SLX-11932.RPI05.HJTY3BGX2.s_1.r_1.counts:    756 DHX36__sensor__2__ULTRA-3226132
#../counts/SLX-11932.RPI05.HJTY3BGX2.s_1.r_1.counts:   4726 DHX36__sherwood_1u__3__ULTRA-3226135
#../counts/SLX-11932.RPI05.HJTY3BGX2.s_1.r_1.counts:   3516 DHX36__sherwood_1u__4__ULTRA-3226134
#../counts/SLX-11932.RPI05.HJTY3BGX2.s_1.r_1.counts:    877 DHX36__sherwood_1u__5__ULTRA-3424948
#../counts/SLX-11932.RPI05.HJTY3BGX2.s_1.r_1.counts:   1201 DHX36__sherwood_1u__7__ULTRA-3424950
#../counts/SLX-11932.RPI05.HJTY3BGX2.s_1.r_1.counts:   1118 DHX36__sherwood_1u__8__ULTRA-3424951
#../counts/SLX-11932.RPI06.HJTY3BGX2.s_1.r_1.counts:    403 DHX36__sensor__1__ULTRA-3226137
#../counts/SLX-11932.RPI06.HJTY3BGX2.s_1.r_1.counts:    550 DHX36__sensor__2__ULTRA-3226132
#../counts/SLX-11932.RPI06.HJTY3BGX2.s_1.r_1.counts:   4781 DHX36__sherwood_1u__3__ULTRA-3226135
#../counts/SLX-11932.RPI06.HJTY3BGX2.s_1.r_1.counts:   3378 DHX36__sherwood_1u__4__ULTRA-3226134
#../counts/SLX-11932.RPI06.HJTY3BGX2.s_1.r_1.counts:    677 DHX36__sherwood_1u__5__ULTRA-3424948
#../counts/SLX-11932.RPI06.HJTY3BGX2.s_1.r_1.counts:   1176 DHX36__sherwood_1u__7__ULTRA-3424950
#../counts/SLX-11932.RPI06.HJTY3BGX2.s_1.r_1.counts:   1153 DHX36__sherwood_1u__8__ULTRA-3424951


grep "RFP" ../counts/SLX-11932.RPI*.HJTY3BGX2.s_1.r_1.counts # some are there

grep "EGFP" ../counts/SLX-11932.RPI*.HJTY3BGX2.s_1.r_1.counts # none are there

```

Use python to parse the `.counts` files and put together a table:

```python
import os
import csv
from Bio import SeqIO

libraries_dict = {}

with open('/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/fastq/SLX-11932.HJTY3BGX2.s_1.contents.csv', 'rb') as csvfile:
  reader = csv.reader(csvfile)
  header = reader.next()
  for row in reader:
    libraries_dict["%s.%s.HJTY3BGX2.s_1.r_1" % (row[0], row[1])] = row[3]


hp_lib = {}
hps = set()

for f in os.listdir("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/counts"):
  ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/counts/%s" % f, "r")
  ilines = ifile.readlines()
  ifile.close()
  print f.replace(".counts", "")
  for l in ilines[1:]:
    fields = l.split()
    c = fields[0]
    hp = fields[1]
    hp_lib[(hp, f.replace(".counts", ""))] = c
    hps.add(hp)


len(hps) # 7919, number of unique hairpins aligned in all libraries, 100 * (7919/8006) = 98.9% of the total found in 20170120_Miniscreen_directory.fa.

all_hps = set()
for record in SeqIO.parse("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/reference/20170120_Miniscreen_directory.fa", "fasta"):
  all_hps.add(record.id)


len(all_hps) # 8006 - 7919 = 87 are missing

all_hps.issuperset(hps) # True
all_hps.difference(hps)


ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/20170407_notfound.txt", "w")

for hp in sorted(all_hps.difference(hps)):
  ofile.write("%s\n" % hp)

ofile.close()


ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/20170407_counts.txt", "w")

ofile.write("key\t%s\n" % ("\t".join([libraries_dict[lib] for lib in sorted(libraries_dict)])))

for hp in sorted(hps):
  ofile.write("%s" % hp)
  for lib in sorted(libraries_dict):
    if (hp,lib) in hp_lib:
      ofile.write("\t%s" % hp_lib[(hp, lib)])
    else:
      ofile.write("\t0")
  ofile.write("\n")

ofile.close()

```




## Share relevant count table, fasta file and other relevant files in the cri-public

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404
rsync -arvuP reference/20170120_Miniscreen_directory.fa reference/20170120_Miniscreen_directory_commonseqs.txt 20170407_notfound.txt 20170407_counts.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170410/tables/

```






# Moving `.fastq` files to archive and remove unnecessary files

Moving `.fastq` files to `martin03@nas-srv001:/archive/Groups/SBLab/fs05/`. In `uk-cri-lcst01`,

```bash
rsync -arvuP /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404/fastq/* martin03@nas-srv001:/archive/Groups/SBLab/fs05/martin03/repository/fastq/
```

Remove unnecessary files. In `uk-cri-lcst01`,

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170404
rm -rf fastq fastqc fastq_trim sam counts

```

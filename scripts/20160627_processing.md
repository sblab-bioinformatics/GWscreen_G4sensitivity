# Processing of fastq files to count table (Pools 1 and 2)


## Location of fastq files

The fastq files from the sequencing of the A375 cell line are available in sblab-srv001 here:

```bash
/media/solexa3/ProcessedRuns/160626_NS500222_0187_H557MBGXY/fastq # They appear here first and then they go to:
/media/staging/160626_NS500222_0187_H557MBGXY/fastq
```

There are 24 sample files and a lostreads file corresponding to pools 1 and 2

- t0: 6 (3 pool1 and 3 pool2)
- t15: 18 (9 pool1 (3 DMSO, 3 PDS and 3 PhenDC3) and 9 pool2 (3 DMSO, 3 PDS and 3 PhenDC3))

In uk-cri-lcst01,

```bash
mkdir /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160627
mkdir /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160627/fastq
rsync --progress martin03@nm149s012925:/home/cri.camres.org/martin03/Desktop/20160627_SLX11618_Pool1_and_2_Seq_Samples.xlsx /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160627
```

In sblab-srv001,

```bash
rsync --progress /media/staging/160626_NS500222_0187_H557MBGXY/fastq/*.fq.gz martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160627/fastq
rsync --progress /media/staging/160626_NS500222_0187_H557MBGXY/fastq/SLX-11618.H557MBGXY.s_1.barcodesummary.txt martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160627/fastq
rsync --progress /media/staging/160626_NS500222_0187_H557MBGXY/fastq/SLX-11618.H557MBGXY.s_1.contents.csv martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160627/fastq
```


## Quality check

In uk-cri-lcst01,

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160627/fastq
mkdir ../fastqc

bsub -J "fastqc" -oo ../fastqc/SLX-11618.H557MBGXY.log -R "rusage[mem=4096]" "/home/martin03/sw/fastqc/fastqc --noextract -q -o ../fastqc SLX-11618.RPI* &&
rsync ../fastqc/*.html martin03@nm149s012925:/home/cri.camres.org/martin03/Desktop
"
```

In nm149s012925, having a look at the reports:

```bash
cd /home/cri.camres.org/martin03/Desktop
firefox *.html &
```

The drop in quality of base 1 is not as significant as in SLX-11622


## Trimming

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160627/fastq
mkdir ../fastq_trim

for fq in SLX-11618.RPI*;
do
  bsub -R "rusage[mem=4096]" "zcat $fq | fastx_trimmer -Q33 -l22 | pigz > ../fastq_trim/${fq%.fq.gz}.trim.fq.gz"
done
```


## Alignment using `bowtie2`

### Generate fasta file containing hairpins

The file was already generated before. It is here:

`/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3.fa`

It was generated for SLX-11622 before.


### Generate index for hairpins fasta file using bowtie2-build

The indexes were already generated before and are available here:

`/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/*.bt2`


### Align

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160627/fastq_trim
IND="/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3"
mkdir ../sam/

for fq in SLX-11618.RPI*;
do
bsub -J $fq -oo ../sam/${fq%.trim.fq.gz}.log -R "rusage[mem=4096]" "zcat $fq | bowtie2 --norc --no-hd -x $IND -U - -S ../sam/${fq%.trim.fq.gz}.sam"
done

cd ../sam

for log in SLX-11618.RPI*.log;
do
echo $log
grep "overall alignment rate" $log
done

```

82-92% of reads with at least one reported alignment


### Build table of counts

Generating counts:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160627/sam
mkdir ../counts

for sam in SLX-11618.RPI*.sam;
do
bsub -J $sam -R "rusage[mem=4096]" "cat $sam | cut -f3 | sort | uniq -c > ../counts/${sam%.sam}.counts"
done

```

Now in an interactive session, use python to parse the `.counts` files and put together a table:

```bash
bsub -R "rusage[mem=8192]" -q interactive -Is bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160627/counts
```

```python
import os
import csv

libraries_dict = {}

with open('/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160627/fastq/SLX-11618.H557MBGXY.s_1.contents.csv', 'rb') as csvfile:
  reader = csv.reader(csvfile)
  header = reader.next()
  for row in reader:
    libraries_dict["%s.%s.H557MBGXY.s_1.r_1" % (row[0], row[1])] = row[3]


hp_lib = {}
hps = set()

for f in os.listdir("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160627/counts"):
  ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160627/counts/%s" % f, "r")
  ilines = ifile.readlines()
  ifile.close()
  print f.replace(".counts", "")
  for l in ilines[1:]:
    fields = l.split()
    c = fields[0]
    hp = fields[1]
    hp_lib[(hp, f.replace(".counts", ""))] = c
    hps.add(hp)


len(hps) # 57454, number of unique hairpins aligned in all libraries, 100 * (57454.0/ 113002) = 51% of the total found in Hs_WG_pool_data_P003_with_12Nic_3.csv (see 113002 above)


ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160627/20160627_counts.txt", "w")

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



# Processing of fastq files to count table (Pools 5 and 10)


## Location of fastq files

The fastq files from the sequencing of the A375 cell line are available in sblab-srv001 here:

```bash
/media/solexa3/160629_NS500222_0189_H3WYFBGXY/fastq # They appear here first, then they probably go to the ProcessedRuns/ folder, and finally they probably go to /media/staging:
/media/staging/160629_NS500222_0189_H3WYFBGXY/fastq
```

There are 24 sample files and a lostreads file corresponding to pools 5 and 10

- t0: 6 (3 pool5 and 3 pool10)
- t15: 18 (9 pool5 (3 DMSO, 3 PDS and 3 PhenDC3) and 9 pool10 (3 DMSO, 3 PDS and 3 PhenDC3))

In uk-cri-lcst01,

```bash
mkdir /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160629
mkdir /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160629/fastq
```

In sblab-srv001,

```bash
rsync --progress /media/staging/160629_NS500222_0189_H3WYFBGXY/fastq/*.fq.gz martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160629/fastq
rsync --progress /media/staging/160629_NS500222_0189_H3WYFBGXY/fastq/SLX-12004.H3WYFBGXY.s_1.barcodesummary.txt martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160629/fastq
rsync --progress /media/staging/160629_NS500222_0189_H3WYFBGXY/fastq/SLX-12004.H3WYFBGXY.s_1.contents.csv martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160629/fastq
```

After comparing with `../../20160627/fastq/SLX-11618.H557MBGXY.s_1.barcodesummary.txt`, there are more reads overall (196240664 to 249078388) but also higher proportion of lost reads (9.31% to 11.69%)


## Quality check

In uk-cri-lcst01,

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160629/fastq
mkdir ../fastqc

bsub -J "fastqc" -oo ../fastqc/SLX-12004.H3WYFBGXY.log -R "rusage[mem=4096]" "/home/martin03/sw/fastqc/fastqc --noextract -q -o ../fastqc SLX-12004.RPI* &&
rsync ../fastqc/*.html martin03@nm149s012925:/home/cri.camres.org/martin03/Desktop
"
```

In nm149s012925, having a look at the reports:

```bash
cd /home/cri.camres.org/martin03/Desktop
firefox *.html &
```

The drop in quality of base 1 is again significant here compared to `../../20160627/fastq/SLX-11618.*.fq.gz`


## Trimming

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160629/fastq
mkdir ../fastq_trim

for fq in SLX-12004.RPI*;
do
  bsub -R "rusage[mem=4096]" "zcat $fq | fastx_trimmer -Q33 -l22 | pigz > ../fastq_trim/${fq%.fq.gz}.trim.fq.gz"
done
```


## Alignment using `bowtie2`

### Generate fasta file containing hairpins

The file was already generated before. It is here:

`/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3.fa`


### Generate index for hairpins fasta file using bowtie2-build

The indexes were already generated before and are available here:

`/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/*.bt2`


### Align

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160629/fastq_trim
IND="/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3"
mkdir ../sam/

for fq in SLX-12004.RPI*;
do
bsub -J $fq -oo ../sam/${fq%.trim.fq.gz}.log -R "rusage[mem=4096]" "zcat $fq | bowtie2 --norc --no-hd -x $IND -U - -S ../sam/${fq%.trim.fq.gz}.sam"
done

cd ../sam

for log in SLX-12004.RPI*.log;
do
echo $log
grep "overall alignment rate" $log
done

```

80-89% of reads with at least one reported alignment. Slightly lower than in `../../20160629/sam/SLX-12004.*.sam`


### Build table of counts

Generating counts:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160629/sam
mkdir ../counts

for sam in SLX-12004.RPI*.sam;
do
bsub -J $sam -R "rusage[mem=4096]" "cat $sam | cut -f3 | sort | uniq -c > ../counts/${sam%.sam}.counts"
done

```

Now in an interactive session, use python to parse the `.counts` files and put together a table:

```bash
bsub -R "rusage[mem=8192]" -q interactive -Is bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160629/counts
```

```python
import os
import csv

libraries_dict = {}

with open('/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160629/fastq/SLX-12004.H3WYFBGXY.s_1.contents.csv', 'rb') as csvfile:
  reader = csv.reader(csvfile)
  header = reader.next()
  for row in reader:
    libraries_dict["%s.%s.H3WYFBGXY.s_1.r_1" % (row[0], row[1])] = row[3]


hp_lib = {}
hps = set()

for f in os.listdir("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160629/counts"):
  ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160629/counts/%s" % f, "r")
  ilines = ifile.readlines()
  ifile.close()
  print f.replace(".counts", "")
  for l in ilines[1:]:
    fields = l.split()
    c = fields[0]
    hp = fields[1]
    hp_lib[(hp, f.replace(".counts", ""))] = c
    hps.add(hp)


len(hps) # 49570, number of unique hairpins aligned in all libraries, 100 * (49570.0/ 113002) = 44% of the total found in Hs_WG_pool_data_P003_with_12Nic_3.csv (see 113002 above)


ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160629/20160629_counts.txt", "w")

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


## Moving `.fastq` files to archive and remove unnecessary files (Pools 1, 2, 5 and 10)

Moving `.fastq` files to archive. `/archive/Groups/SBLab/fs04/` is no longer available, was set to read-only (we can still copy from here). `/archive/Groups/SBLab/fs05/` is the new archive. In `nas-srv001`.

```bash
cd /archive/Groups/SBLab/fs05
mkdir -p martin03/repository/fastq # make subdirectories recursively
cp /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160627/fastq/* /archive/Groups/SBLab/fs05/martin03/repository/fastq/
cp /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160629/fastq/* /archive/Groups/SBLab/fs05/martin03/repository/fastq/
```

Remove unnecessary files. In `uk-cri-lcst01`.

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160627
rm -rf fastq fastqc fastq_trim sam counts
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160629
rm -rf fastq fastqc fastq_trim sam counts

```

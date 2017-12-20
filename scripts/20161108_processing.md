
This script is writen on the basis of what has already been developed before in [20160627_processing.md](20160627_processing.md), [20160902_processing.md](20160902_processing.md) and [20161011_processing.md](20161011_processing.md).

There will be two sequencing runs: one by Darcie and a second one by Katie. These data are reruns so this it complement (in some cases substitute) data that was already generated before.


# Processing of fastq files to count table (Pools 1 (x12 libraries), 2 (x12), 6 (x2) 10 (x2) and 12 (x1)) (Darcie)

Pools 1 and 2 have been resequenced completely. Pool 6 (t0 replicates 2 and 3). Pool 10 (PDS replicate 3 and PhenDC3 replicate 1). Pool 12 (DMSO replicate 2).


## Location of fastq files

ids: 161103_NS500222_0240_HVCC5BGXY, SLX11926

There are 29 sample files and a lostreads file:

- Pool 1:
  - t0 (x3: 1_A375_pool1_rep1_thawedT0, 2_A375_pool1_rep2_thawedT0 and 3_A375_pool1_rep3_thawedT0)
  - t15 (x9):
    - DMSO (x3: 13_A375_pool1_DMSO_rep1, 14_A375_pool1_DMSO_rep2 and 15_A375_pool1_DMSO_rep3)
    - PDS (x3: 16_A375_pool1_PDS_rep1, 17_A375_pool1_PDS_rep2 and 18_A375_pool1_PDS_rep3)
    - PhenDC3 (x3: 19_A375_pool1_PhenDC3_rep1, 20_A375_pool1_PhenDC3_rep2 and 21_A375_pool1_PhenDC3_rep3)
- Pool 2:
  - t0 (x3: 4_A375_pool2_rep1_thawedT0, 5_A375_pool2_rep2_thawedT0 and 6_A375_pool2_rep3_thawedT0)
  - t15 (x9):
    - DMSO (x3: 22_A375_pool2_DMSO_rep1, 23_A375_pool2_DMSO_rep2 and 24_A375_pool2_DMSO_rep3)
    - PDS (x3: 25_A375_pool2_PDS_rep1, 26_A375_pool2_PDS_rep2 and 27_A375_pool2_PDS_rep3)
    - PhenDC3 (x3: 28_A375_pool2_PhenDC3_rep1, 29_A375_pool2_PhenDC3_rep2 and 30_A375_pool2_PhenDC3_rep3)
- Pool 6:
  - t0 (x2: pool_6_t0_replica_2_8 and pool_6_t0_replica_3_9)
- Pool 10
  - t15 (x2):
    - PDS (x1: pool_10_t15_PDS_replica_3_45)
    - PhenDC3 (x1: pool_10_t15_PhenDC3_replica_1_46)
- Pool 12
  - t15 (x1):
    - DMSO (x1: pool_12_t15_DMSO_replica_2_35)


The fastq files from the sequencing are usually available in `sblab-srv001` in `/media/staging`.

First, create receiving folders in `uk-cri-lcst01`:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data
mkdir 20161108
mkdir 20161108/fastq

```

Second, copy them from `sblab-srv001` to `uk-cri-lcst01`,

```bash
rsync --progress /media/staging/161103_NS500222_0240_HVCC5BGXY/fastq/*.fq.gz martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161108/fastq
rsync --progress /media/staging/161103_NS500222_0240_HVCC5BGXY/fastq/SLX-11926.HVCC5BGXY.s_1.barcodesummary.txt martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161108/fastq
rsync --progress /media/staging/161103_NS500222_0240_HVCC5BGXY/fastq/SLX-11926.HVCC5BGXY.s_1.contents.csv martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161108/fastq

```

Looking at `SLX-11926.HVCC5BGXY.s_1.barcodesummary.txt` and comparing them with the previous barcodesummary files. Here we have:

* Total: 404899418. Lost: 11.55%


## Quality check

In uk-cri-lcst01,

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161108/fastq
mkdir ../fastqc

bsub -J "fastqc" -oo ../fastqc/SLX-11926.HVCC5BGXY.log -R "rusage[mem=4096]" "/home/martin03/sw/fastqc/fastqc --noextract --nogroup -q -o ../fastqc SLX-11926.RPI* &&
rsync ../fastqc/*.html martin03@nm149s012925:/home/cri.camres.org/martin03/Desktop"
```

In nm149s012925, having a look at the reports:

```bash
cd /home/cri.camres.org/martin03/Desktop
firefox *.html &
```

There is drop in quality of base 1 too. Rest of things look fine.


## Trimming

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161108/fastq
mkdir ../fastq_trim

for fq in SLX-11926.RPI*;
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
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161108/fastq_trim
IND="/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3"
mkdir ../sam/

for fq in SLX-11926.RPI*;
do
bsub -J $fq -oo ../sam/${fq%.trim.fq.gz}.log -R "rusage[mem=4096]" "zcat $fq | bowtie2 --norc --no-hd -x $IND -U - -S ../sam/${fq%.trim.fq.gz}.sam"
done

cd ../sam

for log in SLX-11926.RPI*.log;
do
echo $log
grep "overall alignment rate" $log
done

```

81-92% of reads with at least one reported alignment. OK.


### Build table of counts

Generating counts:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161108/sam
mkdir ../counts

for sam in SLX-11926.RPI*.sam;
do
bsub -J $sam -R "rusage[mem=4096]" "cat $sam | cut -f3 | sort | uniq -c > ../counts/${sam%.sam}.counts"
done

```

Use python to parse the `.counts` files and put together a table:

```python
import os
import csv

libraries_dict = {}

with open('/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161108/fastq/SLX-11926.HVCC5BGXY.s_1.contents.csv', 'rb') as csvfile:
  reader = csv.reader(csvfile)
  header = reader.next()
  for row in reader:
    libraries_dict["%s.%s.HVCC5BGXY.s_1.r_1" % (row[0], row[1])] = row[3]


hp_lib = {}
hps = set()

for f in os.listdir("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161108/counts"):
  ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161108/counts/%s" % f, "r")
  ilines = ifile.readlines()
  ifile.close()
  print f.replace(".counts", "")
  for l in ilines[1:]:
    fields = l.split()
    c = fields[0]
    hp = fields[1]
    hp_lib[(hp, f.replace(".counts", ""))] = c
    hps.add(hp)


len(hps) # 97850, number of unique hairpins aligned in all libraries, 100 * (97850.0/ 113002) = 87% of the total found in Hs_WG_pool_data_P003_with_12Nic_3.csv. This was 51% for pools 1 and 2, 44% 5 and 10, 77% 3 and 4, 78% 6 and 7, 47% 8 and 9 and 68% 11 and 12.


ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161108/20161108_counts.txt", "w")

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


# Processing of fastq files to count table (Pools 3 (x4 libraries), 4 (x2), 5 (x12), 10 (x12) and 11 (x1)) (Katie)

Pools 5 and 10 have been resequenced completely. Pool 3 (t0 replicates 2 and 3 and PhenDC3 replicates 1 and 3). Pool 4 (t0 replicates 1 and 2). Pool 11 (DMSO replicate 3).


## Location of fastq files

ids: 161106_NS500222_0241_HLVCCBGXY, SLX11927

There are 31 sample files and a lostreads file:

- Pool 3:
  - t0 (x2: 2_A375_pool3_rep2_T0 and 3_A375_pool3_rep3_T0)
  - t15 (x2):
    - PhenDC3 (x2: pool_3_t15_PhenDC3_replica_1_25 and pool_3_t15_PhenDC3_replica_3_27)
- Pool 4:
  - t0 (x2: 4_A375_pool4_rep1_T0 and 5_A375_pool4_rep2_T0)
- Pool 5:
  - t0 (x3: P5_T0_rep1, P5_T0_rep2 and P5_T0_rep3)
  - t15 (x9):
    - DMSO (x3: P5_DMSO_rep1, P5_DMSO_rep2 and P5_DMSO_rep3)
    - PDS (x3: P5_PDS_rep1, P5_PDS_rep2 and P5_PDS_rep3)
    - PhenDC3 (x3: P5_PhenDC3_rep1, P5_PhenDC3_rep2 and P5_PhenDC3_rep3)
- Pool 10:
  - t0 (x3: P10_T0_rep1, P10_T0_rep2 and P10_T0_rep3)
  - t15 (x9):
    - DMSO (x3: P10_DMSO_rep1, P10_DMSO_rep2 and P10_DMSO_rep3)
    - PDS (x3: P10_PDS_rep1, P10_PDS_rep2 and P10_PDS_rep3)
    - PhenDC3 (x3: P10_PhenDC3_rep1, P10_PhenDC3_rep2 and P10_PhenDC3_rep3)
- Pool 11:
  - t15 (x1):
    - DMSO (x1: pool_11_t15_DMSO_replica_3_1)


The fastq files from the sequencing are usually available in `sblab-srv001` in `/media/staging`.

First, create receiving folders in `uk-cri-lcst01`:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data
mkdir 20161109
mkdir 20161109/fastq

```

Second, copy them from `sblab-srv001` to `uk-cri-lcst01`,

```bash
rsync --progress /media/staging/161106_NS500222_0241_HLVCCBGXY/fastq/*.fq.gz martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161109/fastq
rsync --progress /media/staging/161106_NS500222_0241_HLVCCBGXY/fastq/SLX-11927.HLVCCBGXY.s_1.barcodesummary.txt martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161109/fastq
rsync --progress /media/staging/161106_NS500222_0241_HLVCCBGXY/fastq/SLX-11927.HLVCCBGXY.s_1.contents.csv martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161109/fastq

```

Looking at `SLX-11927.HLVCCBGXY.s_1.barcodesummary.txt` and comparing them with the previous barcodesummary files. Here we have:

* Total: 395832364. Lost: 12.06%


## Quality check

In uk-cri-lcst01,

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161109/fastq
mkdir ../fastqc

bsub -J "fastqc" -oo ../fastqc/SLX-11927.HLVCCBGXY.log -R "rusage[mem=4096]" "/home/martin03/sw/fastqc/fastqc --noextract --nogroup -q -o ../fastqc SLX-11927.RPI* &&
rsync ../fastqc/*.html martin03@nm149s012925:/home/cri.camres.org/martin03/Desktop"
```

In nm149s012925, having a look at the reports:

```bash
cd /home/cri.camres.org/martin03/Desktop
firefox *.html &
```

There is drop in quality of base 1 too. Rest of things look fine.


## Trimming

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161109/fastq
mkdir ../fastq_trim

for fq in SLX-11927.RPI*;
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
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161109/fastq_trim
IND="/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3"
mkdir ../sam/

for fq in SLX-11927.RPI*;
do
bsub -J $fq -oo ../sam/${fq%.trim.fq.gz}.log -R "rusage[mem=4096]" "zcat $fq | bowtie2 --norc --no-hd -x $IND -U - -S ../sam/${fq%.trim.fq.gz}.sam"
done

cd ../sam

for log in SLX-11927.RPI*.log;
do
echo $log
grep "overall alignment rate" $log
done

```

80-90% of reads with at least one reported alignment. OK.


### Build table of counts

Generating counts:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161109/sam
mkdir ../counts

for sam in SLX-11927.RPI*.sam;
do
bsub -J $sam -R "rusage[mem=4096]" "cat $sam | cut -f3 | sort | uniq -c > ../counts/${sam%.sam}.counts"
done

```

Use python to parse the `.counts` files and put together a table:

```python
import os
import csv

libraries_dict = {}

with open('/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161109/fastq/SLX-11927.HLVCCBGXY.s_1.contents.csv', 'rb') as csvfile:
  reader = csv.reader(csvfile)
  header = reader.next()
  for row in reader:
    libraries_dict["%s.%s.HLVCCBGXY.s_1.r_1" % (row[0], row[1])] = row[3]


hp_lib = {}
hps = set()

for f in os.listdir("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161109/counts"):
  ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161109/counts/%s" % f, "r")
  ilines = ifile.readlines()
  ifile.close()
  print f.replace(".counts", "")
  for l in ilines[1:]:
    fields = l.split()
    c = fields[0]
    hp = fields[1]
    hp_lib[(hp, f.replace(".counts", ""))] = c
    hps.add(hp)


len(hps) # 98067, number of unique hairpins aligned in all libraries, 100 * (98067.0/ 113002) = 87% of the total found in Hs_WG_pool_data_P003_with_12Nic_3.csv. This was 51% for pools 1 and 2, 44% 5 and 10, 77% 3 and 4, 78% 6 and 7, 47% 8 and 9, 68% 11 and 12 and also 87% for Darcie's rerun above.


ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161109/20161109_counts.txt", "w")

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


# Moving `.fastq` files to archive and remove unnecessary files

Moving `.fastq` files to `/archive/Groups/SBLab/fs05/`. In `nas-srv001`.

```bash
cp /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161108/fastq/* /archive/Groups/SBLab/fs05/martin03/repository/fastq/
cp /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161109/fastq/* /archive/Groups/SBLab/fs05/martin03/repository/fastq/
```

Remove unnecessary files. In `uk-cri-lcst01`.

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161108
rm -rf fastq fastqc fastq_trim sam counts
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161109
rm -rf fastq fastqc fastq_trim sam counts

```

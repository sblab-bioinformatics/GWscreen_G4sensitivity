
# Processing of fastq files to count table


## Location of fastq files

Depending on Genomics, the fastq files might be generated in the following directories:

```bash
/media/solexa3/ProcessedRuns/
/media/staging/
```

As an example we copied the 21 fastq files corresponding to the A375 library here:

```bash
/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/fastq
```

Each file corresponds to both a timepoint and drug treatment:

- t0: 3 (untreated)
- t7: 9 (3 DMSO, 3 PDS and 3 PhenDC3)
- t14: 9 (3 DMSO, 3 PDS and 3 PhenDC3)


## Quality check

In uk-cri-lcst01,

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/fastq
mkdir ../fastqc

bsub -J "fastqc" -oo ../fastqc/SLX-11622.HYLFFBGXX.log -R "rusage[mem=4096]" "/home/martin03/sw/fastqc/fastqc --noextract -q -o ../fastqc SLX-11622.RPI* &&
rsync ../fastqc/*.html martin03@nm149s012925:/home/cri.camres.org/martin03/Desktop
"
```

In nm149s012925, having a look at the reports:

```bash
cd /home/cri.camres.org/martin03/Desktop
firefox *.html &
```


## Trimming

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/fastq
mkdir ../fastq_trim

for fq in SLX-11622.RPI*;
do
  bsub -R "rusage[mem=4096]" "zcat $fq | fastx_trimmer -Q33 -l22 | pigz > ../fastq_trim/${fq%.fq.gz}.trim.fq.gz"
done
```


## Alignment using `bowtie`

### Generate fasta file containing hairpins

```python
with open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3.csv", "rb") as f:
	data = f.readlines()

len(data) # 113003, 113002 hairpins

ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3.fa", "w")

for hp in data[1:]:
  fields=hp.split(",")
  i=fields[4]
  antisense=fields[6][59:81]
  pool="".join([fields[8].split()[0], fields[8].split()[1][1:]])
  ofile.write(">%s\n%s\n" % ("_".join([i, pool]), antisense))

ofile.close()

```


### Generate index for hairpins fasta file using bowtie-build

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420
bowtie-build Hs_WG_pool_data_P003_with_12Nic_3.fa Hs_WG_pool_data_P003_with_12Nic_3
# bowtie -c Hs_WG_pool_data_P003_with_12Nic_3 TGTAATTCTAATGTTTCTGTGT

```


### Align

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/fastq_trim
IND="/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3"
mkdir ../sam

for fq in SLX-11622.RPI*;
do
bsub -J $fq -oo ../sam/${fq%.trim.fq.gz}.log -R "rusage[mem=8192]" "zcat $fq | bowtie -S --sam-nosq -v 3 --norc $IND - ../sam/${fq%.trim.fq.gz}.sam"
done

cd ../sam

for log in SLX-11622.RPI*.log;
do
grep "one reported alignment" $log
done

```

96-98% of reads with at least one reported alignment (allowing 3 mismatches max.)


### Build table of counts

Generating counts:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/sam
mkdir ../counts

for sam in SLX-11622.RPI*.sam;
do
bsub -J $sam -R "rusage[mem=4096]" "tail --lines=+3 $sam | cut -f3 | sort | uniq -c > ../counts/${sam%.sam}.counts"
done

```

Now in an interactive session, use python to parse the `.counts` files and put together a table:

```bash
bsub -R "rusage[mem=8192]" -q interactive -Is bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/counts
```

```python
import os
import csv

libraries_dict = {}

with open('/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/fastq/SLX-11622.HYLFFBGXX.s_1.contents.csv', 'rb') as csvfile:
  reader = csv.reader(csvfile)
  header = reader.next()
  for row in reader:
    libraries_dict["%s.%s.HYLFFBGXX.s_1.r_1" % (row[0], row[1])] = row[3]


hp_lib = {}
hps = set()

for f in os.listdir("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/counts"):
  ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/counts/%s" % f, "r")
  ilines = ifile.readlines()
  ifile.close()
  print f.replace(".counts", "")
  for l in ilines[1:]:
    fields = l.split()
    c = fields[0]
    hp = fields[1]
    hp_lib[(hp, f.replace(".counts", ""))] = c
    hps.add(hp)


len(hps) # 39379, number of unique hairpins aligned in all libraries, 100 * (39379 / 113002) = 35% of the total found in Hs_WG_pool_data_P003_with_12Nic_3.csv (see above)


ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/20160620_counts.txt", "w")

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


## Alignment using `bowtie2`

### Generate fasta file containing hairpins

I used `Hs_WG_pool_data_P003_with_12Nic_3.fa` generated above for `bowtie`.


### Generate index for hairpins fasta file using bowtie2-build

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420
bowtie2-build Hs_WG_pool_data_P003_with_12Nic_3.fa Hs_WG_pool_data_P003_with_12Nic_3
# bowtie2 -c --norc --no-sq -x Hs_WG_pool_data_P003_with_12Nic_3 TGTAATTCTAATGTTTCTGTGT

```


### Align

```bash

cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/fastq_trim
IND="/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3"
rm ../sam/*

for fq in SLX-11622.RPI*;
do
bsub -J $fq -oo ../sam/${fq%.trim.fq.gz}.log -R "rusage[mem=4096]" "zcat $fq | bowtie2 --norc --no-hd -x $IND -U - -S ../sam/${fq%.trim.fq.gz}.sam"
done

cd ../sam

for log in SLX-11622.RPI*.log;
do
grep "overall alignment rate" $log
done

```

92-95% of reads with at least one reported alignment


### Build table of counts

Generating counts:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/sam
rm ../counts/*

for sam in SLX-11622.RPI*.sam;
do
bsub -J $sam -R "rusage[mem=4096]" "cat $sam | cut -f3 | sort | uniq -c > ../counts/${sam%.sam}.counts"
done

```

Now in an interactive session, use python to parse the `.counts` files and put together a table:

```bash
bsub -R "rusage[mem=8192]" -q interactive -Is bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/counts
```

```python
import os
import csv

libraries_dict = {}

with open('/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/fastq/SLX-11622.HYLFFBGXX.s_1.contents.csv', 'rb') as csvfile:
  reader = csv.reader(csvfile)
  header = reader.next()
  for row in reader:
    libraries_dict["%s.%s.HYLFFBGXX.s_1.r_1" % (row[0], row[1])] = row[3]


hp_lib = {}
hps = set()

for f in os.listdir("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/counts"):
  ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/counts/%s" % f, "r")
  ilines = ifile.readlines()
  ifile.close()
  print f.replace(".counts", "")
  for l in ilines[1:]:
    fields = l.split()
    c = fields[0]
    hp = fields[1]
    hp_lib[(hp, f.replace(".counts", ""))] = c
    hps.add(hp)


len(hps) # 38817, number of unique hairpins aligned in all libraries, 100 * (38817 / 113002) = 34% of the total found in Hs_WG_pool_data_P003_with_12Nic_3.csv (see 113002 above)


ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/20160622_counts.txt", "w")

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


## Moving `.fastq` files to archive and remove unnecessary files

Moving `.fastq` files to archive. In `nas-srv001`.

A375:

```bash
cp /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/fastq/* /archive/Groups/SBLab/fs04/martin03/repository/fastq/
```

HT1080:

```bash
cp /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/fastq/*.fq.gz /archive/Groups/SBLab/fs04/martin03/repository/fastq/
```

Save also samplesheets for HT1080. In `sblab-srv001`:

```bash
cd /media/solexa3/ProcessedRuns/160418_NS500222_0156_H2F33BGXY/fastq
rsync SLX-11624.H2F33BGXY.s_1.barcodesummary.txt martin03@nas-srv001:/archive/Groups/SBLab/fs04/martin03/repository/fastq/
rsync SLX-11624.H2F33BGXY.s_1.contents.csv martin03@nas-srv001:/archive/Groups/SBLab/fs04/martin03/repository/fastq/
```

Remove unnecessary files. In `uk-cri-lcst01`.

A375:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503
rm -rf fastq fastqc fastq_trim sam counts
```

HT1080:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420
rm -rf fastq
```

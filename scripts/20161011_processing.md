
I wrote this script on the basis of what was already developed before in [20160627_processing.md](20160627_processing.md) and [20160902_processing.md](20160902_processing.md).

# Processing of fastq files to count table (Pools 8 and 9) (Darcie)

## Location of fastq files

ids: 161008_NS500222_0227_HLWNGBGXY, SLX8948

There are 24 sample files and a lostreads file corresponding to pools 8 and 9

- t0: 6 (3 pool8 and 3 pool9)
- t15: 18 (9 pool8 (3 DMSO, 3 PDS and 3 PhenDC3) and 9 pool9 (3 DMSO, 3 PDS and 3 PhenDC3))

The fastq files from the sequencing are usually available in `sblab-srv001` in `/media/staging`.

First, create receiving folders in `uk-cri-lcst01`:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data
mkdir 20161011
mkdir 20161011/fastq

```

Second, copy them from `sblab-srv001` to `uk-cri-lcst01`,

```bash
rsync --progress /media/staging/161008_NS500222_0227_HLWNGBGXY/fastq/*.fq.gz martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161011/fastq
rsync --progress /media/staging/161008_NS500222_0227_HLWNGBGXY/fastq/SLX-8948.HLWNGBGXY.s_1.barcodesummary.txt martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161011/fastq
rsync --progress /media/staging/161008_NS500222_0227_HLWNGBGXY/fastq/SLX-8948.HLWNGBGXY.s_1.contents.csv martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161011/fastq

```

Looking at `SLX-8948.HLWNGBGXY.s_1.barcodesummary.txt` and comparing them with the previous barcodesummary files. Here we have:

* Pools 8 and 9 (Darcie). Total: 281951161. Lost: 7.89%


## Quality check

In uk-cri-lcst01,

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161011/fastq
mkdir ../fastqc

bsub -J "fastqc" -oo ../fastqc/SLX-8948.HLWNGBGXY.log -R "rusage[mem=4096]" "/home/martin03/sw/fastqc/fastqc --noextract -q -o ../fastqc SLX-8948.RPI* &&
rsync ../fastqc/*.html martin03@nm149s012925:/home/cri.camres.org/martin03/Desktop
"
```

In nm149s012925, having a look at the reports:

```bash
cd /home/cri.camres.org/martin03/Desktop
firefox *.html &
```

There is drop in quality of base 1 too.


## Trimming

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161011/fastq
mkdir ../fastq_trim

for fq in SLX-8948.RPI*;
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
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161011/fastq_trim
IND="/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3"
mkdir ../sam/

for fq in SLX-8948.RPI*;
do
bsub -J $fq -oo ../sam/${fq%.trim.fq.gz}.log -R "rusage[mem=4096]" "zcat $fq | bowtie2 --norc --no-hd -x $IND -U - -S ../sam/${fq%.trim.fq.gz}.sam"
done

cd ../sam

for log in SLX-8948.RPI*.log;
do
echo $log
grep "overall alignment rate" $log
done

```

90-95% of reads with at least one reported alignment. OK.



### Build table of counts

Generating counts:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161011/sam
mkdir ../counts

for sam in SLX-8948.RPI*.sam;
do
bsub -J $sam -R "rusage[mem=4096]" "cat $sam | cut -f3 | sort | uniq -c > ../counts/${sam%.sam}.counts"
done

```

Now in an interactive session, use python to parse the `.counts` files and put together a table:

```bash
bsub -R "rusage[mem=8192]" -q interactive -Is bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161011/counts
```

```python
import os
import csv

libraries_dict = {}

with open('/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161011/fastq/SLX-8948.HLWNGBGXY.s_1.contents.csv', 'rb') as csvfile:
  reader = csv.reader(csvfile)
  header = reader.next()
  for row in reader:
    libraries_dict["%s.%s.HLWNGBGXY.s_1.r_1" % (row[0], row[1])] = row[3]


hp_lib = {}
hps = set()

for f in os.listdir("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161011/counts"):
  ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161011/counts/%s" % f, "r")
  ilines = ifile.readlines()
  ifile.close()
  print f.replace(".counts", "")
  for l in ilines[1:]:
    fields = l.split()
    c = fields[0]
    hp = fields[1]
    hp_lib[(hp, f.replace(".counts", ""))] = c
    hps.add(hp)


len(hps) # 53643, number of unique hairpins aligned in all libraries, 100 * (53643.0/ 113002) = 47% of the total found in Hs_WG_pool_data_P003_with_12Nic_3.csv. This was 51% for pools 1 and 2, 44% 5 and 10, 77% 3 and 4 and 78% 6 and 7.


ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161011/20161011_counts.txt", "w")

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




# Processing of fastq files to count table (Pools 11 and 12) (Katie)

## Location of fastq files

ids: 161017_NS500222_0230_HV5CWBGXY, SLX8944

There are 24 sample files and a lostreads file corresponding to pools 11 and 12

- t0: 6 (3 pool11 and 3 pool12)
- t15: 18 (9 pool11 (3 DMSO, 3 PDS and 3 PhenDC3) and 9 pool12 (3 DMSO, 3 PDS and 3 PhenDC3))

The fastq files from the sequencing are usually available in `sblab-srv001` in `/media/staging`.

First, create receiving folders in `uk-cri-lcst01`:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data
mkdir 20161018
mkdir 20161018/fastq

```

Second, copy them from `sblab-srv001` to `uk-cri-lcst01`,

```bash
rsync --progress /media/staging/161017_NS500222_0230_HV5CWBGXY/fastq/*.fq.gz martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161018/fastq
rsync --progress /media/staging/161017_NS500222_0230_HV5CWBGXY/fastq/SLX-8944.HV5CWBGXY.s_1.barcodesummary.txt martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161018/fastq
rsync --progress /media/staging/161017_NS500222_0230_HV5CWBGXY/fastq/SLX-8944.HV5CWBGXY.s_1.contents.csv martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161018/fastq

```

Looking at `SLX-8948.HLWNGBGXY.s_1.barcodesummary.txt` and comparing them with the previous barcodesummary files. Here we have:

* Pools 11 and 12 (Katie). Total: 323971979. Lost: 10.44%


## Quality check

In uk-cri-lcst01,

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161018/fastq
mkdir ../fastqc

bsub -J "fastqc" -oo ../fastqc/SLX-8944.HV5CWBGXY.log -R "rusage[mem=4096]" "/home/martin03/sw/fastqc/fastqc --noextract -q -o ../fastqc SLX-8944.RPI* &&
rsync ../fastqc/*.html martin03@nm149s012925:/home/cri.camres.org/martin03/Desktop
"
```

In nm149s012925, having a look at the reports:

```bash
cd /home/cri.camres.org/martin03/Desktop
firefox *.html &
```

There is drop in quality of base 1 too.


## Trimming

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161018/fastq
mkdir ../fastq_trim

for fq in SLX-8944.RPI*;
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
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161018/fastq_trim
IND="/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3"
mkdir ../sam/

for fq in SLX-8944.RPI*;
do
bsub -J $fq -oo ../sam/${fq%.trim.fq.gz}.log -R "rusage[mem=4096]" "zcat $fq | bowtie2 --norc --no-hd -x $IND -U - -S ../sam/${fq%.trim.fq.gz}.sam"
done

cd ../sam

for log in SLX-8944.RPI*.log;
do
echo $log
grep "overall alignment rate" $log
done

```

85-93% of reads with at least one reported alignment. OK.


### Build table of counts

Generating counts:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161018/sam
mkdir ../counts

for sam in SLX-8944.RPI*.sam;
do
bsub -J $sam -R "rusage[mem=4096]" "cat $sam | cut -f3 | sort | uniq -c > ../counts/${sam%.sam}.counts"
done

```

Now in an interactive session, use python to parse the `.counts` files and put together a table:

```bash
bsub -R "rusage[mem=8192]" -q interactive -Is bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161018/counts
```

```python
import os
import csv

libraries_dict = {}

with open('/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161018/fastq/SLX-8944.HV5CWBGXY.s_1.contents.csv', 'rb') as csvfile:
  reader = csv.reader(csvfile)
  header = reader.next()
  for row in reader:
    libraries_dict["%s.%s.HV5CWBGXY.s_1.r_1" % (row[0], row[1])] = row[3]


hp_lib = {}
hps = set()

for f in os.listdir("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161018/counts"):
  ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161018/counts/%s" % f, "r")
  ilines = ifile.readlines()
  ifile.close()
  print f.replace(".counts", "")
  for l in ilines[1:]:
    fields = l.split()
    c = fields[0]
    hp = fields[1]
    hp_lib[(hp, f.replace(".counts", ""))] = c
    hps.add(hp)


len(hps) # 77314, number of unique hairpins aligned in all libraries, 100 * (77314.0/ 113002) = 68% of the total found in Hs_WG_pool_data_P003_with_12Nic_3.csv. This was 51% for pools 1 and 2, 44% 5 and 10, 77% 3 and 4, 78% 6 and 7 and 47% 8 and 9.


ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161018/20161018_counts.txt", "w")

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


## Moving `.fastq` files to archive and remove unnecessary files (Pools 8, 9, 11 and 12)

Moving `.fastq` files to `/archive/Groups/SBLab/fs05/`. In `nas-srv001`.

```bash
cp /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161011/fastq/* /archive/Groups/SBLab/fs05/martin03/repository/fastq/
cp /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161018/fastq/* /archive/Groups/SBLab/fs05/martin03/repository/fastq/
```

Remove unnecessary files. In `uk-cri-lcst01`.

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161011
rm -rf fastq fastqc fastq_trim sam counts
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161018
rm -rf fastq fastqc fastq_trim sam counts

```

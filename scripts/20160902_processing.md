
I wrote this script on the basis of what was already developed before in [20160627_processing.md](20160627_processing.md)

# Processing of fastq files to count table (Pools 3 and 4) (Darcie)

## Location of fastq files

ids: 160831_NS500222_0214_HFJVYBGXY, SLX11621

There are 4 lane files for each of the 24 samples corresponding to pools 3 and 4:

- t0: 6 (3 pool3 and 3 pool4)
- t15: 18 (9 pool3 (3 DMSO, 3 PDS and 3 PhenDC3) and 9 pool4 (3 DMSO, 3 PDS and 3 PhenDC3))


The fastq files from the sequencing are usually available in `sblab-srv001` in `/media/staging`, however they taking ages to publish them there.

20160906: after contacting genomics, the sample was not sent through the LIMS system, therefore it would not be published. It is now underway.

I am going to try to get them using [BaseMount](https://basemount.basespace.illumina.com/) and then copy them to uk-cri-lcst01.

First, create receiving folders in uk-cri-lcst01:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data
mkdir 20160902
mkdir 20160902/fastq

```

In nm149s012925 (local computer):

```bash
cd /Users/martin03/Desktop
mkdir BaseSpace
basemount BaseSpace/
cd BaseSpace/Projects/20160831\ pool\ 3\ and\ 4\ T0TF\ SLX11621/Samples/

for folder in `ls`
do
  echo $folder
  cat $folder/Files/*.fastq.gz > ~/Desktop/SLX-11621_$folder.fastq.gz # concatenate reads from different lanes
  rsync --remove-source-files ~/Desktop/SLX-11621_$folder.fastq.gz martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/fastq
done

cd /Users/martin03/Desktop
basemount --unmount BaseSpace/ # unmount the mount point
rm -r BaseSpace # delete folder
```


## Quality check

In uk-cri-lcst01,

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/fastq
mkdir ../fastqc

bsub -J "fastqc" -oo ../fastqc/SLX-11621.HFJVYBGXY.log -R "rusage[mem=4096]" "/home/martin03/sw/fastqc/fastqc --noextract -q -o ../fastqc *.fastq.gz &&
rsync ../fastqc/*.html martin03@nm149s012925:/home/cri.camres.org/martin03/Desktop
"
```

In nm149s012925, having a look at the reports:

```bash
cd /home/cri.camres.org/martin03/Desktop
firefox *.html &
```

Some libraries, e.g. SLX-11621_3C_A375_pool3_rep3_T0.fastq.gz, drop in quality significantly in the 2nd nucleotide


## Trimming

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/fastq
mkdir ../fastq_trim

for fq in SLX-11621*;
do
  bsub -R "rusage[mem=4096]" "zcat $fq | fastx_trimmer -Q33 -l22 | pigz > ../fastq_trim/${fq%.fastq.gz}.trim.fq.gz"
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
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/fastq_trim
IND="/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3"
mkdir ../sam/

for fq in SLX-11621*;
do
bsub -J $fq -oo ../sam/${fq%.trim.fq.gz}.log -R "rusage[mem=4096]" "zcat $fq | bowtie2 --norc --no-hd -x $IND -U - -S ../sam/${fq%.trim.fq.gz}.sam"
done

cd ../sam

for log in SLX-11621*.log;
do
echo $log
grep "overall alignment rate" $log
done

```

90-96% of reads with at least one reported alignment. Good.


### Build table of counts

Generating counts:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/sam
mkdir ../counts

for sam in SLX-11621*.sam;
do
bsub -J $sam -R "rusage[mem=4096]" "cat $sam | cut -f3 | sort | uniq -c > ../counts/${sam%.sam}.counts"
done

```

Now in an interactive session, use python to parse the `.counts` files and put together a table:

```bash
bsub -R "rusage[mem=8192]" -q interactive -Is bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/counts

```

Unlike for the previous pools - see [20160627_processing.md](20160627_processing.md) - here we don't have a `.csv` file to parse to get ids. The ids are already in the name of the `.counts` files.

```python
import os

hp_lib = {}
libs_dict = {}
hps = set()
libs = set()

for f in os.listdir("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/counts"):
  ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/counts/%s" % f, "r")
  ilines = ifile.readlines()
  ifile.close()
  lib_id = "_".join(f.replace(".counts", "").split("_")[1:])
  lib_n = int(f.replace(".counts", "").split("_")[1:][0].replace("C", ""))
  print lib_id
  libs.add(lib_id)
  libs_dict[lib_n] = lib_id
  for l in ilines[1:]:
    fields = l.split()
    c = fields[0]
    hp = fields[1]
    hp_lib[(hp, lib_id)] = c
    hps.add(hp)


len(hps) # 87470, number of unique hairpins aligned in all libraries, 100 * (87470.0/ 113002) = 77% of the total found in Hs_WG_pool_data_P003_with_12Nic_3.csv. This was 51% for pools 1 and 2, and 44% for pools 5 and 10.


ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/20160902_counts.txt", "w")

ofile.write("key\t%s\n" % ("\t".join([libs_dict[lib] for lib in sorted(libs_dict)])))

for hp in sorted(hps):
  ofile.write("%s" % hp)
  for lib_n in sorted(libs_dict):
    if (hp,libs_dict[lib_n]) in hp_lib:
      ofile.write("\t%s" % hp_lib[(hp, libs_dict[lib_n])])
    else:
      ofile.write("\t0")
  ofile.write("\n")

ofile.close()

```





# 20160916 - files for pools 3 and 4 already available through the LIMS!

The fastq files were finally published in `martin03@sblab-srv001:/media/staging`. We are also going to copy them to `/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/fastq`.

The main difference is that the genomics pipeline has preprocessed the reads creating a file of lostreads. I have to rerun the processing to create a new count table and analysis again to have everything processed the same way as in 20160627/, 20160629/ and 20160920/. See below.

In sblab-srv001,

```bash
rsync --progress /media/staging/160831_NS500222_0214_HFJVYBGXY/fastq/*.fq.gz martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/fastq
rsync --progress /media/staging/160831_NS500222_0214_HFJVYBGXY/fastq/SLX-11621.HFJVYBGXY.s_1.barcodesummary.txt martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/fastq
rsync --progress /media/staging/160831_NS500222_0214_HFJVYBGXY/fastq/SLX-11621.HFJVYBGXY.s_1.contents.csv martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/fastq
```

After comparing `SLX-11621.HFJVYBGXY.s_1.barcodesummary.txt` with `SLX-11618.H557MBGXY.s_1.barcodesummary.txt` and `SLX-12004.H3WYFBGXY.s_1.barcodesummary.txt` (in `/archive/Groups/SBLab/fs05/martin03/repository/fastq`), there are less number of reads overall (158489747 compared to 196240664 (1, 2), 249078388 (5, 10)) but better proportion of lost reads (5.00% compared to 9.31% (1, 2) to 11.69% (5, 10)).

* Pools 1 and 2 (Darcie). Total: 196240664. Lost: 9.31%
* Pools 5 and 10 (Katie). Total: 249078388. Lost: 11.69%
* Pools 3 and 4 (Darcie). Total: 158489747. Lost: 5.00%


## Quality check

In uk-cri-lcst01,

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/fastq

bsub -J "fastqc" -oo ../fastqc/SLX-11621.HH7YWBGXY.lims.log -R "rusage[mem=4096]" "/home/martin03/sw/fastqc/fastqc --noextract -q -o ../fastqc SLX-11621.RPI* &&
rsync ../fastqc/SLX-11621.RPI*.html martin03@nm149s012925:/home/cri.camres.org/martin03/Desktop
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
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/fastq

for fq in SLX-11621.RPI*;
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
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/fastq_trim
IND="/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3"

for fq in SLX-11621.RPI*;
do
bsub -J $fq -oo ../sam/${fq%.trim.fq.gz}.log -R "rusage[mem=4096]" "zcat $fq | bowtie2 --norc --no-hd -x $IND -U - -S ../sam/${fq%.trim.fq.gz}.sam"
done

cd ../sam

for log in SLX-11621.RPI*.log;
do
echo $log
grep "overall alignment rate" $log
done

```

91-96% of reads with at least one reported alignment. OK.


### Build table of counts

Generating counts:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/sam

for sam in SLX-11621.RPI*.sam;
do
bsub -J $sam -R "rusage[mem=4096]" "cat $sam | cut -f3 | sort | uniq -c > ../counts/${sam%.sam}.counts"
done

```

Now in an interactive session, use python to parse the `.counts` files and put together a table:

```bash
bsub -R "rusage[mem=8192]" -q interactive -Is bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/counts
```

```python
import os
import csv

libraries_dict = {}

with open('/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/fastq/SLX-11621.HFJVYBGXY.s_1.contents.csv', 'rb') as csvfile:
  reader = csv.reader(csvfile)
  header = reader.next()
  for row in reader:
    libraries_dict["%s.%s.HFJVYBGXY.s_1.r_1" % (row[0], row[1])] = row[3]


hp_lib = {}
hps = set()

for f in os.listdir("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/counts"):
  if "RPI" in f:
    ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/counts/%s" % f, "r")
    ilines = ifile.readlines()
    ifile.close()
    print f.replace(".counts", "")
    for l in ilines[1:]:
      fields = l.split()
      c = fields[0]
      hp = fields[1]
      hp_lib[(hp, f.replace(".counts", ""))] = c
      hps.add(hp)


len(hps) # 87342, number of unique hairpins aligned in all libraries, 100 * (87342.0/ 113002) = 77% of the total found in Hs_WG_pool_data_P003_with_12Nic_3.csv. This was 51% for pools 1 and 2, 44% for pools 5 and 10, 77 for pools 3 and 4 (no LIMS) and 78 for pools 6 and 7.


ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/20160902_counts_lims.txt", "w")

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







# Processing of fastq files to count table (Pools 6 and 7) (Katie)

## Location of fastq files

ids: 160915_NS500222_0217_HH7YWBGXY, SLX12005

There are 24 sample files and a lostreads file corresponding to pools 6 and 7

- t0: 6 (3 pool6 and 3 pool7)
- t15: 18 (9 pool6 (3 DMSO, 3 PDS and 3 PhenDC3) and 9 pool7 (3 DMSO, 3 PDS and 3 PhenDC3))


The fastq files from the sequencing are usually available in `sblab-srv001` in `/media/staging`.

First, create receiving folders in `uk-cri-lcst01`:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data
mkdir 20160920
mkdir 20160920/fastq

```

Second, copy them from `sblab-srv001` to `uk-cri-lcst01`,

```bash
rsync --progress /media/staging/160915_NS500222_0217_HH7YWBGXY/fastq/*.fq.gz martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160920/fastq
rsync --progress /media/staging/160915_NS500222_0217_HH7YWBGXY/fastq/SLX-12005.HH7YWBGXY.s_1.barcodesummary.txt martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160920/fastq
rsync --progress /media/staging/160915_NS500222_0217_HH7YWBGXY/fastq/SLX-12005.HH7YWBGXY.s_1.contents.csv martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160920/fastq

```

Looking at `SLX-12005.HH7YWBGXY.s_1.barcodesummary.txt` and comparing them with the previous barcodesummary files. Here we have:

* Pools 6 and 7 (Katie). Total: 296102200. Lost: 12.71%


## Quality check

In uk-cri-lcst01,

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160920/fastq
mkdir ../fastqc

bsub -J "fastqc" -oo ../fastqc/SLX-12005.HH7YWBGXY.log -R "rusage[mem=4096]" "/home/martin03/sw/fastqc/fastqc --noextract -q -o ../fastqc SLX-12005.RPI* &&
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
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160920/fastq
mkdir ../fastq_trim

for fq in SLX-12005.RPI*;
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
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160920/fastq_trim
IND="/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3"
mkdir ../sam/

for fq in SLX-12005.RPI*;
do
bsub -J $fq -oo ../sam/${fq%.trim.fq.gz}.log -R "rusage[mem=4096]" "zcat $fq | bowtie2 --norc --no-hd -x $IND -U - -S ../sam/${fq%.trim.fq.gz}.sam"
done

cd ../sam

for log in SLX-12005.RPI*.log;
do
echo $log
grep "overall alignment rate" $log
done

```

85-93% of reads with at least one reported alignment. OK.


### Build table of counts

Generating counts:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160920/sam
mkdir ../counts

for sam in SLX-12005.RPI*.sam;
do
bsub -J $sam -R "rusage[mem=4096]" "cat $sam | cut -f3 | sort | uniq -c > ../counts/${sam%.sam}.counts"
done

```

Now in an interactive session, use python to parse the `.counts` files and put together a table:

```bash
bsub -R "rusage[mem=8192]" -q interactive -Is bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160920/counts
```

```python
import os
import csv

libraries_dict = {}

with open('/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160920/fastq/SLX-12005.HH7YWBGXY.s_1.contents.csv', 'rb') as csvfile:
  reader = csv.reader(csvfile)
  header = reader.next()
  for row in reader:
    libraries_dict["%s.%s.HH7YWBGXY.s_1.r_1" % (row[0], row[1])] = row[3]


hp_lib = {}
hps = set()

for f in os.listdir("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160920/counts"):
  ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160920/counts/%s" % f, "r")
  ilines = ifile.readlines()
  ifile.close()
  print f.replace(".counts", "")
  for l in ilines[1:]:
    fields = l.split()
    c = fields[0]
    hp = fields[1]
    hp_lib[(hp, f.replace(".counts", ""))] = c
    hps.add(hp)


len(hps) # 87802, number of unique hairpins aligned in all libraries, 100 * (87802.0/ 113002) = 78% of the total found in Hs_WG_pool_data_P003_with_12Nic_3.csv. This was 51% for pools 1 and 2, 44% for pools 5 and 10, 77 for pools 3 and 4 (no LIMS).


ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160920/20160920_counts.txt", "w")

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




## Moving `.fastq` files to archive and remove unnecessary files (Pools 3, 4, 6 and 7)

Moving `.fastq` files to `/archive/Groups/SBLab/fs05/`. In `nas-srv001`.

```bash
cp /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902/fastq/* /archive/Groups/SBLab/fs05/martin03/repository/fastq/
cp /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160920/fastq/* /archive/Groups/SBLab/fs05/martin03/repository/fastq/
```

Remove unnecessary files. In `uk-cri-lcst01`.

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160902
rm -rf fastq fastqc fastq_trim sam counts
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160920
rm -rf fastq fastqc fastq_trim sam counts

```


This script is writen on the basis of what has already been developed before in [20160627_processing.md](20160627_processing.md), [20160902_processing.md](20160902_processing.md), [20161011_processing.md](20161011_processing.md) and [20161108_processing.md](20161108_processing.md).

There will be one sequencing run performed by Darcie. These data are reruns so this is it complement (in some cases substitute) data that was already generated before.

Pool 3 (PhenDC3 replicates 1 and 3). Pool 11 (DMSO replicate 3). Pool 12 (DMSO replicate 2).


## Location of fastq files

ids: 170104_NS500222_0256_HCMYGBGX2, SLX11929

There are 16 files. Unlike all previous sequencing runs, each of the 4 sample was split in 4 lanes here.

- Pool 3:
  - t15:
    - PhenDC3 (x8: t15_PhenDC3_pool3_rep1 and t15_PhenDC3_pool3_rep3)
- Pool 11:
  - t15:
    - DMSO (x4: t15_DMSO_pool11_rep3)
- Pool 12:
  - t15:
    - DMSO (x4: t15_DMSO_pool12_rep2)


The fastq files from the sequencing are usually available in `sblab-srv001` in `/media/staging`.

First, create receiving folders in `uk-cri-lcst01`:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data
mkdir 20170109
mkdir 20170109/fastq

```

Second, copy them from `sblab-srv001` to `uk-cri-lcst01`,

```bash
rsync --progress /media/staging/170104_NS500222_0256_HCMYGBGX2/fastq/*.fq.gz martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170109/fastq
rsync --progress /media/staging/170104_NS500222_0256_HCMYGBGX2/fastq/SLX-11929.HCMYGBGX2.s_1.barcodesummary.txt martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170109/fastq
rsync --progress /media/staging/170104_NS500222_0256_HCMYGBGX2/fastq/SLX-11929.HCMYGBGX2.s_1.contents.csv martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170109/fastq

```

Looking at `SLX-11929.HCMYGBGX2.s_1.barcodesummary.txt` and comparing them with the previous barcodesummary files. Here we have:

* 310093182 reads. Lost: 7.88%


## Quality check

In uk-cri-lcst01,

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170109/fastq
mkdir ../fastqc

bsub -J "fastqc" -oo ../fastqc/SLX-11929.HCMYGBGX2.log -R "rusage[mem=4096]" "/home/martin03/sw/fastqc/fastqc --noextract --nogroup -q -o ../fastqc SLX-11929.RPI* &&
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
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170109/fastq
mkdir ../fastq_trim

for fq in SLX-11929.RPI*;
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
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170109/fastq_trim
IND="/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3"
mkdir ../sam/

for fq in SLX-11929.RPI*;
do
bsub -J $fq -oo ../sam/${fq%.trim.fq.gz}.log -R "rusage[mem=4096]" "zcat $fq | bowtie2 --norc --no-hd -x $IND -U - -S ../sam/${fq%.trim.fq.gz}.sam"
done

cd ../sam

for log in SLX-11929.RPI*.log;
do
echo $log
grep "overall alignment rate" $log
done

```

91-95% of reads with at least one reported alignment. OK.


### Build table of counts

Generating counts:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170109/sam
mkdir ../counts

for sam in SLX-11929.RPI*.sam;
do
bsub -J $sam -R "rusage[mem=4096]" "cat $sam | cut -f3 | sort | uniq -c > ../counts/${sam%.sam}.counts"
done

```

Use python to parse the `.counts` files and put together a table:

```python
import os
import csv

libraries_dict = {}

with open('/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170109/fastq/SLX-11929.HCMYGBGX2.s_1.contents.csv', 'rb') as csvfile:
  reader = csv.reader(csvfile)
  header = reader.next()
  for row in reader:
    libraries_dict["%s.%s.HCMYGBGX2.s_1.r_1" % (row[0], row[1])] = row[3]


hp_lib = {}
hps = set()

for f in os.listdir("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170109/counts"):
  ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170109/counts/%s" % f, "r")
  ilines = ifile.readlines()
  ifile.close()
  print f.replace(".counts", "")
  for l in ilines[1:]:
    fields = l.split()
    c = fields[0]
    hp = fields[1]
    hp_lib[(hp, f.replace(".counts", ""))] = c
    hps.add(hp)


len(hps) # 96730, number of unique hairpins aligned in all libraries, 100 * (96730.0/ 113002) = 85.6% of the total found in Hs_WG_pool_data_P003_with_12Nic_3.csv.


ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170109/20170109_counts.txt", "w")

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

Moving `.fastq` files to `martin03@nas-srv001:/archive/Groups/SBLab/fs05/`. In `uk-cri-lcst01`,

```bash
rsync -arvuP /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170109/fastq/* martin03@nas-srv001:/archive/Groups/SBLab/fs05/martin03/repository/fastq/
```

Remove unnecessary files. In `uk-cri-lcst01`.

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170109
rm -rf fastq fastqc fastq_trim sam counts

```

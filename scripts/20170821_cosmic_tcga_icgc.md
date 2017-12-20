
From Katie and Darcie:

Here are the databases we were interested in:

- COSMIC http://cancer.sanger.ac.uk/cosmic
- TCGA https://cancergenome.nih.gov/
- ICGC: https://dcc.icgc.org/

Ideally, from our GW 758 G4 sensitiser gene list we would like to know:

- Which genes  have loss-of-function mutations?
- How many loss-of-function mutations per gene are there?
- Which cancer types have these  loss-of-function mutations and how common are they (I.e. Does 50% of patients have these mutations?)
- Are there any cancers that have multiple genes that have loss-of-function mutations that appear in our 758 gene list? If so which mutations/proteins?

From COSMIC:

- Can you say if there is a way to access loss-of-function annotations in COSMIC (or perhaps any other cancer databases)?
- At the moment, there is no direct way to obtain this using Cosmic. We currently have an ongoing project, 'The Cancer Mutation Census', that annotates these and other types of mutations. But you could try the following ways for now:
1. You could select frameshift, nonsense mutations and deletions, which are in most cases loss-of-function (but there are cases where truncating mutation results in gain-of-function so it is not ideal).
2. Select tumour suppressor genes (TSGs). TSG means that inactivating mutations in this gene occur in at least some cancers, but this information is available only for some genes for now as this is under development. Check these TSGs for truncating mutations and maybe to use some algorithm annotate the missense mutations.

In other words, COSMIC describes the type of change caused by the mutation but does not seem to annotate direct functional effects as it stands. We would have to investigate whether TCGA/ICGC annotate loss-of-function mutations directly.

Alternatively, we could either (a) take frameshift, nonsense mutations and deletions as suggested by COSMIC or (b) use mutations in general, not loss-of-function neccesarily. But let's look a bit more into TCGA/ICGC or other databases first.

Regarding the GW 758 G4 sensitiser gene list, I counted 763 (with FC threshold, logFC-1) and 843 (without FC threshold) in the list for literature search that Darcie sent on 27th July (attached). Are we working on the same lists of genes?

In regards to the COSMIC list: could you maybe first pull out all the somatic and germ line mutations known for the 758 genes so we can get an idea of how many there are? Is this a difficult thing to do? As long as there is also information attached in regards to type of mutation: I.e. "frameshift, nonsense mutations and deletions etc “, the type of cancer, frequency of mutation etc, I could start shifting through the data this week with you.

David had a talk with Chris today and they definitely both agree that there is a potential for obtaining some intellectual property from the G4 ligand shRNA screen- so I’d like to start looking at the possibilities of which genes could be part of this sooner rather than later.

Have you written/heard back from TCGA/ICGC about how they annotate their gene mutations?



## Data

### GW 758 G4 sensitiser gene list

See attachment in 20170821 email from Katie.

The project is all so far in `uk-cri-lcst01` however the R package TCGAbiolinks could only be installed in `clust1-headnode` (not `uk-cri-lcst01`) so I created a new branch of the project in `clust1-headnode`:

```bash
cd /scratchb/sblab/martin03/repository
mkdir -p 20160418_shRNAscreen_katie_darcie/data/20170822/tables
cd 20160418_shRNAscreen_katie_darcie/data/20170822/tables
rsync -arvuP martin03@nm149s012925:/Users/martin03/Desktop/A375_758_Gene_List_NEW* .
```


### COSMIC

release v82, 3rd August 2017 downloads [page](http://cancer.sanger.ac.uk/cosmic/download)

In `nm149s012925`,

```bash
cd /Users/martin03/Desktop
sftp "sergio.martinezcuesta@cruk.cam.ac.uk"@sftp-cancer.sanger.ac.uk
#password:ptenstop
get /files/grch38/cosmic/v82/CosmicMutantExport.tsv.gz
```

Copy to `clust1-headnode`,

```bash
cd /scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170822
mkdir cosmic
cd cosmic
rsync -arvuP martin03@nm149s012925:/Users/martin03/Desktop/CosmicMutantExport.tsv.gz .
```


### TCGA

- https://gdc.cancer.gov/access-data
- https://www.biostars.org/p/179077
- TCGA2STAT: Simple TCGA Data Access for Integrated Statistical Analysis in R http://www.liuzlab.org/TCGA2STAT/
- TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data. https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html
- maftools https://bioconductor.org/packages/release/bioc/html/maftools.html
- GenomicDataCommons https://bioconductor.org/packages/release/bioc/html/GenomicDataCommons.html
- RTCGA https://github.com/RTCGA/RTCGA
- https://www.r-bloggers.com/facilitating-exploratory-data-visualization-application-to-tcga-genomic-data/


#### TCGAbiolinks

```r
#source("https://bioconductor.org/biocLite.R")
#biocLite("TCGAbiolinks")

library(TCGAbiolinks)
library(dplyr)

#install.packages("DT")
library(DT)

query.maf.hg19 <- GDCquery(project = "TCGA-CHOL", data.category = "Simple nucleotide variation", data.type = "Simple somatic mutation", access = "open", legacy = TRUE)
datatable(select(getResults(query.maf.hg19),-contains("cases")), filter = 'top', options = list(scrollX = TRUE, keys = TRUE, pageLength = 10), rownames = FALSE)

# useful to download .maf files
query.maf.hg19 <- GDCquery(project = "TCGA-CHOL", data.category = "Simple nucleotide variation", data.type = "Simple somatic mutation", access = "open", file.type = "bcgsc.ca_CHOL.IlluminaHiSeq_DNASeq.1.somatic.maf", legacy = TRUE)
GDCdownload(query.maf.hg19)
maf <- GDCprepare(query.maf.hg19)
datatable(maf[1:20,], filter = 'top', options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), rownames = FALSE)

```


#### TCGA2STAT

```r
#source("https://bioconductor.org/biocLite.R")
#biocLite("CNTools")
#install.packages("/Users/martin03/tmp/TCGA2STAT_1.2.tar.gz", repos = NULL, type = "source")
library(TCGA2STAT)
?getTCGA
allmut.ov <- getTCGA(disease="OV", data.type="Mutation", type="all")
str(allmut.ov)
dim(allmut.ov$dat) # 10060   316
allmut.ov$dat["TET1",]

```


### GDC API

- https://gdc.cancer.gov/developers/gdc-application-programming-interface-api
- https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/

```python
import requests
import json

#status
url = 'https://api.gdc.cancer.gov/status'
response = requests.get(url)
print json.dumps(response.json(), indent=2)
response.json()["tag"] # u'1.10.0'
response.json()["data_release"] # u'Data Release 8.0 - August 22, 2017'

#projects
url = 'https://api.gdc.cancer.gov/projects'
response = requests.get(url)
print json.dumps(response.json(), indent=2)

#cases
url = 'https://api.gdc.cancer.gov/cases'
response = requests.get(url)
print json.dumps(response.json(), indent=2)

#files
url = 'https://api.gdc.cancer.gov/files'
response = requests.get(url)
print json.dumps(response.json(), indent=2)

#annotations
url = 'https://api.gdc.cancer.gov/annotations'
response = requests.get(url)
print json.dumps(response.json(), indent=2)

```

- https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/#request-parameters
- https://portal.gdc.cancer.gov/exploration?facetTab=mutations&searchTableTab=cases


### ICGC





## Analysis

### COSMIC

```python
# Loading screen lists of genes
ifile1 = open("/scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170822/tables/A375_758_Gene_List_NEW.csv", "r")
ilines1 = ifile1.readlines()
ifile1.close()

genename_hgncid = {}

for line in ilines1[1:]:
  fields = line.split("\t")
  genename = fields[0]
  hgncid = fields[6].replace("\n", "")
  genename_hgncid[genename] = hgncid


len(genename_hgncid) # 758
len(set(genename_hgncid.values())) # 758


# Loading cosmic mutation data
import gzip

ifile2 = gzip.open('/scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170822/cosmic/CosmicMutantExport.tsv.gz', 'rb')
ilines2 = ifile2.readlines()
ifile2.close()

hgncid_muts = {}
hgncid_genelength = {}

for line in ilines2[1:]:
  fields = line.split("\t")
  genename = fields[0]
  genelength = fields[2]
  hgncid = fields[3]
  mutationstatus = fields[27]
  fathmmscore = fields[28]
  if hgncid in hgncid_muts:
    hgncid_muts[hgncid].append((mutationstatus, fathmmscore))
    hgncid_genelength[hgncid].add(genelength)
  else:
    hgncid_muts[hgncid] = [(mutationstatus, fathmmscore)]
    hgncid_genelength[hgncid] = set([genelength])


len(hgncid_muts) # 17608
len(hgncid_genelength) # 17608

len(set(hgncid_muts.keys()).intersection(set(genename_hgncid.values()))) # 702, there are 56 genes for which there is no mutation data
set(hgncid_muts.keys()).intersection(set(genename_hgncid.values()))
set(genename_hgncid.values()).difference(set(hgncid_muts.keys()))

hgncid_muts[genename_hgncid["BRCA1"]]
hgncid_genelength[genename_hgncid["BRCA1"]]



# Write summary table of genes
import numpy as np

hgncid_summary = {}

for hgncid in hgncid_muts:
  n_total = len(hgncid_muts[hgncid])
  n_pathogenic = len([mut[1] for mut in hgncid_muts[hgncid] if mut[0] == "PATHOGENIC"])
  n_neutral = len([mut[1] for mut in hgncid_muts[hgncid] if mut[0] == "NEUTRAL"])
  mean_pathogenic = np.mean([float(i) for i in [mut[1] for mut in hgncid_muts[hgncid] if mut[0] == "PATHOGENIC"]])
  mean_neutral = np.mean([float(i) for i in [mut[1] for mut in hgncid_muts[hgncid] if mut[0] == "NEUTRAL"]])
  n_pathogenic_norm = float(n_pathogenic)/int(list(hgncid_genelength[hgncid])[0])
  hgncid_summary[hgncid] = [n_total, n_pathogenic, n_neutral, mean_pathogenic, mean_neutral, n_pathogenic_norm]


n_total = len(hgncid_muts[genename_hgncid["BRCA1"]])
n_pathogenic = len([mut[1] for mut in hgncid_muts[genename_hgncid["BRCA1"]] if mut[0] == "PATHOGENIC"])
n_neutral = len([mut[1] for mut in hgncid_muts[genename_hgncid["BRCA1"]] if mut[0] == "NEUTRAL"])
mean_pathogenic = np.mean([float(i) for i in [mut[1] for mut in hgncid_muts[genename_hgncid["BRCA1"]] if mut[0] == "PATHOGENIC"]])
mean_neutral = np.mean([float(i) for i in [mut[1] for mut in hgncid_muts[genename_hgncid["BRCA1"]] if mut[0] == "NEUTRAL"]])
n_pathogenic_norm = float(n_pathogenic)/int(list(hgncid_genelength[genename_hgncid["BRCA1"]])[0])

hgncid_summary['23380'] # [81, 0, 80, nan, 0.014498374999999999, 0.0]
hgncid_summary['1100'] # [702, 163, 393, 0.88182006134969326, 0.11249539440203563, 0.029148783977110158]
hgncid_summary['1101'] # [1100, 269, 593, 0.90511895910780671, 0.12296883642495783, 0.026225992005459685]


ofile = open("/scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170822/tables/A375_758_Gene_List_NEW_cosmic_summary.txt", "w")

ofile.write("Genes (758)\tHGNC ID\tn_total\tn_pathogenic\tn_neutral\tmean_pathogenic\tmean_neutral\tn_pathogenic_norm\n")

for line in ilines1[1:]:
  fields = line.split("\t")
  genename = fields[0]
  hgncid = fields[6].replace("\n", "")
  features = ["nan", "nan", "nan", "nan", "nan", "nan"]
  if hgncid in hgncid_summary:
    features = [str(i) for i in hgncid_summary[hgncid]]
  ofile.write("%s\t%s\t%s\n" % (genename, hgncid, "\t".join(features)))


ofile.close()



# Write table of mutations
ofile2 = open("/scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170822/tables/A375_758_Gene_List_NEW_cosmic_mutations.txt", "w")

ofile2.write(ilines2[0])

for line in ilines2[1:]:
  fields = line.split("\t")
  hgncid = fields[3]
  if hgncid in genename_hgncid.values():
    ofile2.write(line)

ofile2.close()

```

Compressing output tables:

```bash
cd /scratchb/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170822/tables/
pigz A375_758_Gene_List_NEW_cosmic_summary.txt
pigz A375_758_Gene_List_NEW_cosmic_mutations.txt

zcat A375_758_Gene_List_NEW_cosmic_summary.txt.gz | wc -l # 759, fine
zcat A375_758_Gene_List_NEW_cosmic_mutations.txt.gz | wc -l # 146786
zcat A375_758_Gene_List_NEW_cosmic_mutations.txt.gz | awk '$4 == 1100' # BRCA1
zcat A375_758_Gene_List_NEW_cosmic_mutations.txt.gz | awk '$4 == 14410' # DHX36

```




## TODO

- Explore Intermine tool for list of relevant sensitiser genes

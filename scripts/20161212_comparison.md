This script compares the gene hits obtained in [20161206_integration.md](20161206_integration.md) with datasets provided by Jochen, Barbara and Mike Gormally.

Borrows code from [20161130_comparison.md](20161130_comparison.md)

# Jochen's pulldown data

## Files location

### Jochen's data

In uk-cri-lcst01,

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161130/jochen
ls -lh
#total 4.0K
#-rw-r----- 1 martin03 sblab 1.7K Nov 30 16:44 BG4_RIME_hits.csv

```

### Darcie and Katie's data

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data
mkdir -p 20161212/darcie_katie
cd 20161212/darcie_katie
rsync martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/tables/genes_pass_50pct_or_3_threshold/*.txt .
```


## Comparison

```python

# Convert Jochen's Uniprot accessions to Gene names

import urllib,urllib2
import os

data = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161130/jochen/BG4_RIME_hits.csv", "r")
data_all = data.read().split("\n")
data.close()

accs = [line.split(",")[3] for line in data_all[1:-1]]
len(accs) == len(set(accs)) # True, 74 unique accessionsS

url = 'http://www.uniprot.org/uploadlists/'

params = {
'from':'ACC',
'to':'GENENAME',
'format':'tab',
'query':" ".join(accs)
}

data = urllib.urlencode(params)
request = urllib2.Request(url, data)
contact = "sermarcue@gmail.com" # Please set your email address here to help us debug in case of problems.
request.add_header('User-Agent', 'Python %s' % contact)
response = urllib2.urlopen(request)
mapping = response.read()
print(mapping)

len(mapping.split("\n")[1:-1]) # 77, RL9_HUMAN maps to RPL9, RPL9P7, RPL9P8 and RPL9P9


# Load Darcie and Katie's data

PDSandPhenDC3only_resistant = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161212/darcie_katie/20161207_PDSandPhenDC3only_resistant_genes.txt", "r")
PDSandPhenDC3only_sensitiser = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161212/darcie_katie/20161207_PDSandPhenDC3only_sensitiser_genes.txt", "r")
PDSonly_resistant = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161212/darcie_katie/20161207_PDSonly_resistant_genes.txt", "r")
PDSonly_sensitiser = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161212/darcie_katie/20161207_PDSonly_sensitiser_genes.txt", "r")
PhenDC3only_resistant = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161212/darcie_katie/20161207_PhenDC3only_resistant_genes.txt", "r")
PhenDC3only_sensitiser = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161212/darcie_katie/20161207_PhenDC3only_sensitiser_genes.txt", "r")
DMSO_resistant = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161212/darcie_katie/20161207_DMSO_resistant_genes.txt", "r")
DMSO_sensitiser = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161212/darcie_katie/20161207_DMSO_sensitiser_genes.txt", "r")

PDSandPhenDC3only_resistant_genes = [line.split("\t")[0] for line in PDSandPhenDC3only_resistant.read().split("\n")[1:-1]]
PDSandPhenDC3only_sensitiser_genes = [line.split("\t")[0] for line in PDSandPhenDC3only_sensitiser.read().split("\n")[1:-1]]
PDSonly_resistant_genes = [line.split("\t")[0] for line in PDSonly_resistant.read().split("\n")[1:-1]]
PDSonly_sensitiser_genes = [line.split("\t")[0] for line in PDSonly_sensitiser.read().split("\n")[1:-1]]
PhenDC3only_resistant_genes = [line.split("\t")[0] for line in PhenDC3only_resistant.read().split("\n")[1:-1]]
PhenDC3only_sensitiser_genes = [line.split("\t")[0] for line in PhenDC3only_sensitiser.read().split("\n")[1:-1]]
DMSO_resistant_genes = [line.split("\t")[0] for line in DMSO_resistant.read().split("\n")[1:-1]]
DMSO_sensitiser_genes = [line.split("\t")[0] for line in DMSO_sensitiser.read().split("\n")[1:-1]]

PDSandPhenDC3only_resistant.close()
PDSandPhenDC3only_sensitiser.close()
PDSonly_resistant.close()
PDSonly_sensitiser.close()
PhenDC3only_resistant.close()
PhenDC3only_sensitiser.close()
DMSO_resistant.close()
DMSO_sensitiser.close()

len(PDSandPhenDC3only_resistant_genes) # 71
len(PDSandPhenDC3only_sensitiser_genes) # 107
len(PDSonly_resistant_genes) # 237
len(PDSonly_sensitiser_genes) # 228
len(PhenDC3only_resistant_genes) # 415
len(PhenDC3only_sensitiser_genes) # 505
len(DMSO_resistant_genes) # 277
len(DMSO_sensitiser_genes) # 407

darcie_katie = {('PDS and PhenDC3', 'resistant'):PDSandPhenDC3only_resistant_genes, ('PDS and PhenDC3', 'sensitiser'):PDSandPhenDC3only_sensitiser_genes, ('PDS', 'resistant'):PDSonly_resistant_genes, ('PDS', 'sensitiser'):PDSonly_sensitiser_genes, ('PhenDC3', 'resistant'):PhenDC3only_resistant_genes, ('PhenDC3', 'sensitiser'):PhenDC3only_sensitiser_genes, ('DMSO', 'resistant'):DMSO_resistant_genes, ('DMSO', 'sensitiser'):DMSO_sensitiser_genes}


# Putting it all together

output = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161212/darcie_katie/20161212_BG4_RIME_hits_vs_shRNA_screen.txt", "w")
output.write("uniprot_acc\tgene_name\ttreatment\ttype\n")

for pair in mapping.split("\n")[1:-1]:
  fields=pair.split("\t")
  acc=fields[0]
  gene=fields[1]
  print(gene)
  treatment="NA"
  typ="NA"
  for cond in darcie_katie:
    if gene in darcie_katie[cond]:
      treatment=cond[0]
      typ=cond[1]
  output.write("%s\t%s\t%s\t%s\n" % (acc, gene, treatment, typ))

output.close()

os.system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161212/darcie_katie/20161212_BG4_RIME_hits_vs_shRNA_screen.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161212")

```


## Comparison: based on uniprot accs, rather than gene names

```python

import urllib,urllib2
import os

# Load Jochen's Uniprot accessions

data = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161130/jochen/BG4_RIME_hits.csv", "r")
data_all = data.read().split("\n")
data.close()

accs = [line.split(",")[3] for line in data_all[1:-1]]
len(accs) == len(set(accs)) # True, 74 unique accessionsS


# Load Darcie and Katie's data

PDSandPhenDC3only_resistant = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161212/darcie_katie/20161207_PDSandPhenDC3only_resistant_genes.txt", "r")
PDSandPhenDC3only_sensitiser = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161212/darcie_katie/20161207_PDSandPhenDC3only_sensitiser_genes.txt", "r")
PDSonly_resistant = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161212/darcie_katie/20161207_PDSonly_resistant_genes.txt", "r")
PDSonly_sensitiser = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161212/darcie_katie/20161207_PDSonly_sensitiser_genes.txt", "r")
PhenDC3only_resistant = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161212/darcie_katie/20161207_PhenDC3only_resistant_genes.txt", "r")
PhenDC3only_sensitiser = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161212/darcie_katie/20161207_PhenDC3only_sensitiser_genes.txt", "r")
DMSO_resistant = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161212/darcie_katie/20161207_DMSO_resistant_genes.txt", "r")
DMSO_sensitiser = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161212/darcie_katie/20161207_DMSO_sensitiser_genes.txt", "r")

PDSandPhenDC3only_resistant_genes = [line.split("\t")[0] for line in PDSandPhenDC3only_resistant.read().split("\n")[1:-1]]
PDSandPhenDC3only_sensitiser_genes = [line.split("\t")[0] for line in PDSandPhenDC3only_sensitiser.read().split("\n")[1:-1]]
PDSonly_resistant_genes = [line.split("\t")[0] for line in PDSonly_resistant.read().split("\n")[1:-1]]
PDSonly_sensitiser_genes = [line.split("\t")[0] for line in PDSonly_sensitiser.read().split("\n")[1:-1]]
PhenDC3only_resistant_genes = [line.split("\t")[0] for line in PhenDC3only_resistant.read().split("\n")[1:-1]]
PhenDC3only_sensitiser_genes = [line.split("\t")[0] for line in PhenDC3only_sensitiser.read().split("\n")[1:-1]]
DMSO_resistant_genes = [line.split("\t")[0] for line in DMSO_resistant.read().split("\n")[1:-1]]
DMSO_sensitiser_genes = [line.split("\t")[0] for line in DMSO_sensitiser.read().split("\n")[1:-1]]

PDSandPhenDC3only_resistant.close()
PDSandPhenDC3only_sensitiser.close()
PDSonly_resistant.close()
PDSonly_sensitiser.close()
PhenDC3only_resistant.close()
PhenDC3only_sensitiser.close()
DMSO_resistant.close()
DMSO_sensitiser.close()

len(PDSandPhenDC3only_resistant_genes) # 71
len(PDSandPhenDC3only_sensitiser_genes) # 107
len(PDSonly_resistant_genes) # 237
len(PDSonly_sensitiser_genes) # 228
len(PhenDC3only_resistant_genes) # 415
len(PhenDC3only_sensitiser_genes) # 505
len(DMSO_resistant_genes) # 277
len(DMSO_sensitiser_genes) # 407

darcie_katie = {('PDS and PhenDC3', 'resistant'):PDSandPhenDC3only_resistant_genes, ('PDS and PhenDC3', 'sensitiser'):PDSandPhenDC3only_sensitiser_genes, ('PDS', 'resistant'):PDSonly_resistant_genes, ('PDS', 'sensitiser'):PDSonly_sensitiser_genes, ('PhenDC3', 'resistant'):PhenDC3only_resistant_genes, ('PhenDC3', 'sensitiser'):PhenDC3only_sensitiser_genes, ('DMSO', 'resistant'):DMSO_resistant_genes, ('DMSO', 'sensitiser'):DMSO_sensitiser_genes}


# Get uniprot accessions for Darcie and Katie's data

darcie_katie_accs = {}

url = 'http://www.uniprot.org/uploadlists/'

for (tr, ty) in darcie_katie.keys():
  print tr, ty
  params = {'from':'GENENAME', 'to':'ACC', 'format':'tab', 'query':" ".join(darcie_katie[(tr, ty)])}
  data = urllib.urlencode(params)
  request = urllib2.Request(url, data)
  contact = "sermarcue@gmail.com"
  request.add_header('User-Agent', 'Python %s' % contact)
  response = urllib2.urlopen(request)
  mapping = response.read()
  accessions = []
  for line in mapping.split("\n")[1:-1]:
    fields = line.split("\t")
    genename = fields[0]
    acc = fields[3]
    rev = fields[4]
    org = fields[7]
    if rev == "reviewed" and org.startswith("Homo sapiens"):
      accessions.append(acc)
  darcie_katie_accs[(tr, ty)] = accessions


#test
params = {'from':'GENENAME', 'to':'ACC', 'format':'tab', 'query':'OR2C3 PABPN1', 'reviewed':'yes', 'organism':'9606'}
data = urllib.urlencode(params)
request = urllib2.Request(url, data)
contact = "sermarcue@gmail.com"
request.add_header('User-Agent', 'Python %s' % contact)
response = urllib2.urlopen(request)
mapping = response.read()




```

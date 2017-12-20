
This script aims to perform an analysis of the functions of sensitiser genes that came up significant.

I will be using several tools. Ideas:
- edgeR tools https://bioconductor.org/packages/release/bioc/html/edgeR.html
- goseq http://bioconductor.org/packages/release/bioc/html/goseq.html
- GOstats https://www.bioconductor.org/packages/release/bioc/html/GOstats.html
- gage https://bioconductor.org/packages/release/bioc/html/gage.html
- pathview http://bioconductor.org/packages/release/bioc/html/pathview.html
- GSEA http://software.broadinstitute.org/gsea/msigdb/index.jsp
- DAVID https://david.ncifcrf.gov/

Follow the literature at the bottom of [20170120_integration_finaldataset.md](20170120_integration_finaldataset.md)



## Copy data

In `uk-cri-lcst01`,

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data
mkdir 20170307
rsync -arvuP martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/genes_pass_50pct_or_3_threshold/20170130_*sensitiser*.txt 20170307/
```



## edgeR tools

Parse `Hs_WG_pool_data_P003_with_12Nic_3.csv` to create a mapping between gene name and gene id:

```python
#read file
with open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3.csv", "rb") as f:
    data = f.readlines()

len(data) # 113003, 113002 hairpins

#define a dictionary containing gene names and gene ids
genename_geneid = {}

for hp in data[1:]:
  fields=hp.split(",")
  gene_name=fields[1].split()[0]
  gene_id=fields[2]
  if gene_name in genename_geneid:
    genename_geneid[gene_name].add(gene_id)
  else:
    genename_geneid[gene_name]=set([gene_id])

len(genename_geneid.keys()) # 18937

#check gene names linked to many gene ids

genename_counter = [(genename, len(genename_geneid[genename])) for genename in genename_geneid]

from operator import itemgetter
sorted(genename_counter, key=itemgetter(1), reverse=True)[:20]

ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3.gene_name.gene_id.txt", "w")

for genename in genename_geneid:
  for geneid in list(genename_geneid[genename]):
    ofile.write("%s\t%s\n" % (genename, geneid))

ofile.close()

```

Using the standard `goana()`, `kegga()`, `roast()`, `mroast()`, `fry()`, `camera()`, `romer()` in edgeR - see Chen2015, Ritchie2014 and [20160531_figure_editing.md](https://github.com/sblab-bioinformatics/projects/blob/master/20150501_methylation_brain/20160531_figure_editing.md#analysis-of-expression-of-top-10-5hmc-differences-between-tumour-and-margin)


```R
library(edgeR)
library(data.table)
library(ggplot2)
library(cowplot)

# Adjust width
options(width=250)

# Load lists of sensitiser genes
pds_only <- fread("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170307/20170130_PDSonly_sensitiser_genes.txt")
phendc3_only <- fread("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170307/20170130_PhenDC3only_sensitiser_genes.txt")
pds_phendc3_only <- fread("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170307/20170130_PDSandPhenDC3only_sensitiser_genes.txt")
dmso <- fread("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170307/20170130_DMSO_sensitiser_genes.txt")

nrow(pds_only) # 224
nrow(phendc3_only) # 512
nrow(pds_phendc3_only) # 107
nrow(dmso) # 408

# Load mapping between gene names and gene ids
genename_geneid <- fread("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3.gene_name.gene_id.txt")
setnames(genename_geneid, c("gene_name","gene_id"))

# Merge
pds_only <- merge(pds_only, genename_geneid, by.x = 'gene', by.y='gene_name', all.x=TRUE)
sum(is.na(pds_only$gene_id)) # 2
pds_only[is.na(pds_only$gene_id)] # MARCH3 MARCH6
pds_only[gene=="MARCH3", gene_id := 115123] # http://www.uniprot.org/uniprot/Q86UD3
pds_only[gene=="MARCH6", gene_id := 10299] # http://www.uniprot.org/uniprot/O60337
sum(is.na(pds_only$gene_id)) # 0

phendc3_only <- merge(phendc3_only, genename_geneid, by.x='gene', by.y='gene_name', all.x=TRUE)
sum(is.na(phendc3_only$gene_id)) # 0

pds_phendc3_only <- merge(pds_phendc3_only, genename_geneid, by.x='gene', by.y='gene_name', all.x=TRUE)
sum(is.na(pds_phendc3_only$gene_id)) # 0

dmso <- merge(dmso, genename_geneid, by.x='gene', by.y='gene_name', all.x=TRUE)
sum(is.na(dmso$gene_id)) # 1
dmso[is.na(dmso$gene_id)] # MARCH9
dmso[gene=="MARCH9", gene_id := 92979] # http://www.uniprot.org/uniprot/Q86YJ5
sum(is.na(dmso$gene_id)) # 0

########################################################################################################
# An alternative to using the annotations provided in the sequence file by transomic is to use biomaRt:

#library(biomaRt)

# List biomart databases and versions related to ensembl
#listMarts(host="www.ensembl.org")

# List datasets within database "ENSEMBL_MART_ENSEMBL" as in the December 2016 version
#mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "dec2016.archive.ensembl.org") # http://www.ensembl.org/info/website/archives/index.html
#listDatasets(mart)

# Select dataset hsapiens_gene_ensembl, which is Human genes (GRCh38.p7), and explore
#mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset= 'hsapiens_gene_ensembl', host = "dec2016.archive.ensembl.org")
#head(listAttributes(mart))
#head(listFilters(mart))

# Collect all gene names and entrezgenes
#all_geneid <- data.table(getBM(attributes = c("external_gene_name", "entrezgene"), mart = mart))
#all_geneid <- all_geneid[!is.na(all_geneid$entrezgene),]
#all_geneid <- all_geneid[, list(entrezgene = paste(entrezgene, collapse= "|")), by= external_gene_name]
#stopifnot(length(all_geneid$external_gene_name) == length(unique(all_geneid$external_gene_name)))

# Merge
#pds_only <- merge(pds_only, all_geneid, by.x= 'gene', by.y= 'external_gene_name', all.x= TRUE)
#phendc3_only <- merge(phendc3_only, all_geneid, by.x= 'gene', by.y= 'external_gene_name', all.x= TRUE)
#pds_phendc3_only <- merge(pds_phendc3_only, all_geneid, by.x= 'gene', by.y= 'external_gene_name', all.x= TRUE)
#dmso <- merge(dmso, all_geneid, by.x= 'gene', by.y= 'external_gene_name', all.x= TRUE)

#sum(is.na(pds_only$entrezgene)) # 23
#sum(is.na(phendc3_only$entrezgene)) # 29
#sum(is.na(pds_phendc3_only$entrezgene)) # 9
#sum(is.na(dmso$entrezgene)) # 40

# This is due to synonyms
#head(pds_only, n=20)
# BOD1L is a synonym of BOD1L1, which is the real gene name
#all_geneid[external_gene_name == "BOD1L"]
#Empty data.table (0 rows) of 2 cols: external_gene_name,entrezgene
#all_geneid[external_gene_name == "BOD1L1"]
#   external_gene_name entrezgene
#1:             BOD1L1     259282

# This is not ideal so we might have to use uniprot api to get these geneid queries right
######################################################################################################################

# GO analysis

## Define sensitiser group
sensitiser_geneids <- unique(c(pds_only$gene_id, phendc3_only$gene_id, pds_phendc3_only$gene_id))
length(sensitiser_geneids) # 843
length(pds_only$gene_id) + length(phendc3_only$gene_id) + length(pds_phendc3_only$gene_id) # 843

## Run goana for different groups and correct p-values
### pds_only
go_pds_only <- goana(de = pds_only$gene_id, species = "Hs")
# universe = genename_geneid$gene_id (transomic) OR universe = all_geneid$entrezgene (ensembl) is also an option. The difference is not big though. Leaving it empty uses org.Hs.eg.db I think (see below).
go_pds_only <- data.table(go_pds_only, keep.rownames=T)
go_pds_only[, FDR := p.adjust(P.DE, "fdr")]
#### MF
go_pds_only[Ont=="MF"][order(FDR)][1:20]
nrow(go_pds_only[Ont=="MF"]) # 4263
#### BP
go_pds_only[Ont=="BP"][order(FDR)][1:20]
nrow(go_pds_only[Ont=="BP"]) # 14291

### phendc3_only
go_phendc3_only <- goana(de = phendc3_only$gene_id, species = "Hs")
go_phendc3_only <- data.table(go_phendc3_only, keep.rownames=T)
go_phendc3_only[, FDR := p.adjust(P.DE, "fdr")]
#### MF
go_phendc3_only[Ont=="MF"][order(FDR)][1:20]
nrow(go_phendc3_only[Ont=="MF"]) # 4263
#### BP
go_phendc3_only[Ont=="BP"][order(FDR)][1:20]
nrow(go_phendc3_only[Ont=="BP"]) # 14291

### pds_phendc3_only
go_pds_phendc3_only <- goana(de = pds_phendc3_only$gene_id, species = "Hs")
go_pds_phendc3_only <- data.table(go_pds_phendc3_only, keep.rownames=T)
go_pds_phendc3_only[, FDR := p.adjust(P.DE, "fdr")]
#### MF
go_pds_phendc3_only[Ont=="MF"][order(FDR)][1:20]
nrow(go_pds_phendc3_only[Ont=="MF"]) # 4263
#### BP
go_pds_phendc3_only[Ont=="BP"][order(FDR)][1:20]
nrow(go_pds_phendc3_only[Ont=="BP"]) # 14291

### sensitiser
go_sensitiser <- goana(de = sensitiser_geneids, species = "Hs")
go_sensitiser <- data.table(go_sensitiser, keep.rownames=T)
go_sensitiser[, FDR := p.adjust(P.DE, "fdr")]
#### MF
go_sensitiser[Ont=="MF"][order(FDR)][1:20]
nrow(go_sensitiser[Ont=="MF"]) # 4263
#### BP
go_pds_phendc3_only[Ont=="BP"][order(FDR)][1:20]
nrow(go_pds_phendc3_only[Ont=="BP"]) # 14291

### dmso
go_dmso <- goana(de = dmso$gene_id, species = "Hs")
go_dmso <- data.table(go_dmso, keep.rownames=T)
go_dmso[, FDR := p.adjust(P.DE, "fdr")]
#### MF
go_dmso[Ont=="MF"][order(FDR)][1:20]
nrow(go_dmso[Ont=="MF"]) # 4263
#### BP
go_dmso[Ont=="BP"][order(FDR)][1:20]
nrow(go_dmso[Ont=="BP"]) # 14291


## Write tables
### pds_only
#### MF
write.table(go_pds_only[Ont=="MF"][order(FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170310_GO_pds_only_MF.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170310_GO_pds_only_MF.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170310/tables/")
#### BP
write.table(go_pds_only[Ont=="BP"][order(FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170310_GO_pds_only_BP.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170310_GO_pds_only_BP.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170310/tables/")

### phendc3_only
#### MF
write.table(go_phendc3_only[Ont=="MF"][order(FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170310_GO_phendc3_only_MF.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170310_GO_phendc3_only_MF.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170310/tables/")
#### BP
write.table(go_phendc3_only[Ont=="BP"][order(FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170310_GO_phendc3_only_BP.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170310_GO_phendc3_only_BP.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170310/tables/")

### pds_phendc3_only
#### MF
write.table(go_pds_phendc3_only[Ont=="MF"][order(FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170310_GO_pds_phendc3_only_MF.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170310_GO_pds_phendc3_only_MF.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170310/tables/")
#### BP
write.table(go_pds_phendc3_only[Ont=="BP"][order(FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170310_GO_pds_phendc3_only_BP.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170310_GO_pds_phendc3_only_BP.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170310/tables/")

### sensitiser
#### MF
write.table(go_sensitiser[Ont=="MF"][order(FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170310_GO_sensitiser_MF.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170310_GO_sensitiser_MF.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170310/tables/")
#### BP
write.table(go_sensitiser[Ont=="BP"][order(FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170310_GO_sensitiser_BP.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170310_GO_sensitiser_BP.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170310/tables/")

### dmso
#### MF
write.table(go_dmso[Ont=="MF"][order(FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170310_GO_dmso_MF.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170310_GO_dmso_MF.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170310/tables/")
#### BP
write.table(go_dmso[Ont=="BP"][order(FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170310_GO_dmso_BP.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170310_GO_dmso_BP.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170310/tables/")

#Visually these are the first terms that make most difference between sensitiser and DMSO
#GO:0030554	adenyl nucleotide binding
#GO:0032559	adenyl ribonucleotide binding


## Barplots
### pds_only
#### MF
gg <- ggplot(go_pds_only[Ont=="MF"][order(FDR)][1:20], aes(x = factor(Term, levels = rev(go_pds_only[Ont=="MF"][order(FDR)][1:20]$Term)), y = -log10(FDR))) +
geom_bar(stat="identity") +
xlab("GO term") +
ylab(expression("-log"[10]*"(FDR)")) +
theme_classic() +
coord_flip()

ggg <- ggdraw(switch_axis_position(gg, axis = 'x')) # https://yutannihilation.github.io/allYourFigureAreBelongToUs/cowplot/switch_axis_position/

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170310_GO_pds_only_MF.pdf", width = 20/2.54, height = 16/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170310_GO_pds_only_MF.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170310/figures/")

#### BP
gg <- ggplot(go_pds_only[Ont=="BP"][order(FDR)][1:20], aes(x = factor(Term, levels = rev(go_pds_only[Ont=="BP"][order(FDR)][1:20]$Term)), y = -log10(FDR))) +
geom_bar(stat="identity") +
xlab("GO term") +
ylab(expression("-log"[10]*"(FDR)")) +
theme_classic() +
coord_flip()

ggg <- ggdraw(switch_axis_position(gg, axis = 'x')) # https://yutannihilation.github.io/allYourFigureAreBelongToUs/cowplot/switch_axis_position/

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170310_GO_pds_only_BP.pdf", width = 18/2.54, height = 16/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170310_GO_pds_only_BP.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170310/figures/")


### phendc3_only
#### MF
gg <- ggplot(go_phendc3_only[Ont=="MF"][order(FDR)][1:20], aes(x = factor(Term, levels = rev(go_phendc3_only[Ont=="MF"][order(FDR)][1:20]$Term)), y = -log10(FDR))) +
geom_bar(stat="identity") +
xlab("GO term") +
ylab(expression("-log"[10]*"(FDR)")) +
theme_classic() +
coord_flip()

ggg <- ggdraw(switch_axis_position(gg, axis = 'x')) # https://yutannihilation.github.io/allYourFigureAreBelongToUs/cowplot/switch_axis_position/

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170310_GO_phendc3_only_MF.pdf", width = 16/2.54, height = 16/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170310_GO_phendc3_only_MF.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170310/figures/")

#### BP
gg <- ggplot(go_phendc3_only[Ont=="BP"][order(FDR)][1:20], aes(x = factor(Term, levels = rev(go_phendc3_only[Ont=="BP"][order(FDR)][1:20]$Term)), y = -log10(FDR))) +
geom_bar(stat="identity") +
xlab("GO term") +
ylab(expression("-log"[10]*"(FDR)")) +
theme_classic() +
coord_flip()

ggg <- ggdraw(switch_axis_position(gg, axis = 'x')) # https://yutannihilation.github.io/allYourFigureAreBelongToUs/cowplot/switch_axis_position/

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170310_GO_phendc3_only_BP.pdf", width = 22/2.54, height = 16/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170310_GO_phendc3_only_BP.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170310/figures/")


### pds_phendc3_only
#### MF
gg <- ggplot(go_pds_phendc3_only[Ont=="MF"][order(FDR)][1:20], aes(x = factor(Term, levels = rev(go_pds_phendc3_only[Ont=="MF"][order(FDR)][1:20]$Term)), y = -log10(FDR))) +
geom_bar(stat="identity") +
xlab("GO term") +
ylab(expression("-log"[10]*"(FDR)")) +
theme_classic() +
coord_flip()

ggg <- ggdraw(switch_axis_position(gg, axis = 'x')) # https://yutannihilation.github.io/allYourFigureAreBelongToUs/cowplot/switch_axis_position/

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170310_GO_pds_phendc3_only_MF.pdf", width = 16/2.54, height = 16/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170310_GO_pds_phendc3_only_MF.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170310/figures/")

#### BP
gg <- ggplot(go_pds_phendc3_only[Ont=="BP"][order(FDR)][1:20], aes(x = factor(Term, levels = rev(go_pds_phendc3_only[Ont=="BP"][order(FDR)][1:20]$Term)), y = -log10(FDR))) +
geom_bar(stat="identity") +
xlab("GO term") +
ylab(expression("-log"[10]*"(FDR)")) +
theme_classic() +
coord_flip()

ggg <- ggdraw(switch_axis_position(gg, axis = 'x')) # https://yutannihilation.github.io/allYourFigureAreBelongToUs/cowplot/switch_axis_position/

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170310_GO_pds_phendc3_only_BP.pdf", width = 18/2.54, height = 16/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170310_GO_pds_phendc3_only_BP.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170310/figures/")


### sensitiser
#### MF
gg <- ggplot(go_sensitiser[Ont=="MF"][order(FDR)][1:20], aes(x = factor(Term, levels = rev(go_sensitiser[Ont=="MF"][order(FDR)][1:20]$Term)), y = -log10(FDR))) +
geom_bar(stat="identity") +
xlab("GO term") +
ylab(expression("-log"[10]*"(FDR)")) +
theme_classic() +
coord_flip()

ggg <- ggdraw(switch_axis_position(gg, axis = 'x')) # https://yutannihilation.github.io/allYourFigureAreBelongToUs/cowplot/switch_axis_position/

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170310_GO_sensitiser_MF.pdf", width = 16/2.54, height = 16/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170310_GO_sensitiser_MF.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170310/figures/")

#### BP
gg <- ggplot(go_sensitiser[Ont=="BP"][order(FDR)][1:20], aes(x = factor(Term, levels = rev(go_sensitiser[Ont=="BP"][order(FDR)][1:20]$Term)), y = -log10(FDR))) +
geom_bar(stat="identity") +
xlab("GO term") +
ylab(expression("-log"[10]*"(FDR)")) +
theme_classic() +
coord_flip()

ggg <- ggdraw(switch_axis_position(gg, axis = 'x')) # https://yutannihilation.github.io/allYourFigureAreBelongToUs/cowplot/switch_axis_position/

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170310_GO_sensitiser_BP.pdf", width = 18/2.54, height = 16/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170310_GO_sensitiser_BP.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170310/figures/")


### dmso
#### MF
gg <- ggplot(go_dmso[Ont=="MF"][order(FDR)][1:20], aes(x = factor(Term, levels = rev(go_dmso[Ont=="MF"][order(FDR)][1:20]$Term)), y = -log10(FDR))) +
geom_bar(stat="identity") +
xlab("GO term") +
ylab(expression("-log"[10]*"(FDR)")) +
theme_classic() +
coord_flip()

ggg <- ggdraw(switch_axis_position(gg, axis = 'x')) # https://yutannihilation.github.io/allYourFigureAreBelongToUs/cowplot/switch_axis_position/

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170310_GO_dmso_MF.pdf", width = 19/2.54, height = 16/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170310_GO_dmso_MF.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170310/figures/")

#### BP
gg <- ggplot(go_dmso[Ont=="BP"][order(FDR)][1:20], aes(x = factor(Term, levels = rev(go_dmso[Ont=="BP"][order(FDR)][1:20]$Term)), y = -log10(FDR))) +
geom_bar(stat="identity") +
xlab("GO term") +
ylab(expression("-log"[10]*"(FDR)")) +
theme_classic() +
coord_flip()

ggg <- ggdraw(switch_axis_position(gg, axis = 'x')) # https://yutannihilation.github.io/allYourFigureAreBelongToUs/cowplot/switch_axis_position/

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170310_GO_dmso_BP.pdf", width = 18/2.54, height = 16/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170310_GO_dmso_BP.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170310/figures/")




## Obtain genes related for a given each GO term, e.g. GO:0015106
library(org.Hs.eg.db)

org.Hs.eg()
#DB schema: HUMAN_DB
#DB schema version: 2.1
#Organism: Homo sapiens
#Date for NCBI data: 2015-Sep27
#Date for GO data: 20150919
#Date for KEGG data: 2011-Mar15
#Date for Golden Path data: 2010-Mar22
#Date for Ensembl data: 2015-Jul16

org.Hs.eg_dbInfo()

x <- org.Hs.egGO2ALLEGS
Rkeys(x) <- "GO:0003676"
EG <- mappedLkeys(x)
EG

all_geneid[all_geneid$entrezgene %in% EG,]



# KEGG analysis
kegg <- kegga(de = pds_only$entrezgene, species = "Hs")
topKEGG(kegg, number = 20)


```



## TODO

- An alternative would be to improve mapping between gene names and gene ids using Uniprot, however we managed to do a good job with the transOMIC .csv file a manual curation

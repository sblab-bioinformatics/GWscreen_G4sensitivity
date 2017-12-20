
This script continues with the analysis of the shRNA screen following discussions taking place in Darcie's group meeting on 20161031 and bringing together previous scripts (follows section (G) in [20160914_integration.md](20160914_integration.md) and analyses before).

The aim is to extract PDS only, PhenDC3 only and PDS+PhenDC3 only genes that change significantly from the table of individual shRNAs.

# General questions

- How does edgeR calculate logFC if the three replicates have zero counts (e.g. "NUP54__1__1_Pool12" has zero counts in the PhenDC3 treatment - see pool12)?
It adds a prior.count of 0.125, see [here](https://support.bioconductor.org/p/63901/).


# Analysis

```R
library(edgeR)
library(data.table)
library(ggplot2)
library(VennDiagram)
library(reshape)
library(Biostrings)


# Enlarge the view width when printing tables
options(width = 250)


# Load table of counts into edgeR like in (B)
data1 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160627_Pool1_counts.txt", header = T)
data2 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160627_Pool2_counts.txt", header = T)
data3 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160908_Pool3_counts_lims.txt", header = T)
data4 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160908_Pool4_counts_lims.txt", header = T)
data5 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160629_Pool5_counts.txt", header = T)
data6 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160920_Pool6_counts.txt", header = T)
data7 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160920_Pool7_counts.txt", header = T)
data8 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20161011_Pool8_counts.txt", header = T)
data9 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20161011_Pool9_counts.txt", header = T)
data10 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160629_Pool10_counts.txt", header = T)
data11 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20161018_Pool11_counts.txt", header = T)
data12 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20161018_Pool12_counts.txt", header = T)


# Combine table of counts into a single data structure
colnames(data1) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data2) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data3) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data4) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data5) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data6) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data7) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data8) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data9) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data10) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data11) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data12) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")

data <- rbind(data1, data2, data3, data4, data5, data6, data7, data8, data9, data10, data11, data12)
dim(data) # 110576     12

# 100*(110576/113002) = 97.9% shRNAs detected


# Change to a DGE object
groups <- c("t0_no", "t0_no", "t0_no", "t15_DMSO", "t15_DMSO", "t15_DMSO", "t15_PDS", "t15_PDS", "t15_PDS", "t15_PhenDC3", "t15_PhenDC3", "t15_PhenDC3")
genes_dataframe <- data.frame(shRNA = rownames(data), gene = as.character(sapply(rownames(data), function(x) unlist(strsplit(x, "__"))[1])))
z <- DGEList(counts = data, group = groups, genes = genes_dataframe)
z$genes$pool <- sapply(rownames(data), function(x) gsub("Pool", "", tail(unlist(strsplit(x, "_")), n=1)))


# Vector containing all pools
pools <- unique(z$genes$pool)


# Define empty data.tables to contain shRNAs for all pools
top_DMSOvst0 <- data.table()
top_PDSvst0 <- data.table()
top_PhenDC3vst0 <- data.table()


# Create directory structure for pools
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/pools")


# Loop iterating through all pools
for (p in pools){
  # Print pool number
  print(sprintf("Pool%s", p))

  # Select pool
  z_pool <- z[z$genes$pool == p,]
  print(sprintf("shRNAs: %s", nrow(z_pool$counts)))

  # Redefine lib.size
  z_pool$samples$lib.size <- as.numeric(colSums(z_pool$counts))

  # Filtering
  keep <- rowSums(cpm(z_pool[,1:3])>0.5) >= 3 # A library of 10M will need to have at least 5 counts in each of three t0 replicates individually in order to pass this filter
  z_pool <- z_pool[keep, , keep.lib.sizes=FALSE]
  print(sprintf("shRNAs (after filtering): %s", nrow(z_pool$counts)))

  # Normalisation
  z_pool <- calcNormFactors(z_pool)

  # Differential representation analysis.
  # Set up design matrix for GLM.
  des <- model.matrix(~ 0 + group, data = z_pool$samples)
  colnames(des) <- levels(factor(z_pool$samples$group))

  # Estimate common dispersion
  z_pool_glm <- estimateDisp(z_pool, des)
  print(sprintf("Common dispersion: %s", sqrt(z_pool_glm$common.disp)))

  # Fit negative binomial GLM
  z_pool_fit <- glmFit(z_pool_glm, des)

  # Define matrix of contrasts
  my.contrasts <- makeContrasts(
  DMSOvst0 = t15_DMSO - t0_no,
  PDSvst0 = t15_PDS - t0_no,
  PhenDC3vst0 = t15_PhenDC3 - t0_no,
  #PDSvsDMSO = t15_PDS - t15_DMSO,
  #PhenDC3vsDMSO = t15_PhenDC3 - t15_DMSO,
  #PhenDC3vsPDS = t15_PhenDC3 - t15_PDS,
  #PDSandPhenDC3vst0 = 0.5*(t15_PhenDC3 + t15_PDS) - t0_no,
  #PDSandPhenDC3vsDMSO = 0.5*(t15_PhenDC3 + t15_PDS) - t15_DMSO,
  #PDSandPhenDC3andDMSOvst0 = 1/3*(t15_PhenDC3 + t15_PDS + t15_DMSO) - t0_no,
  levels=des)

  # Comparisons t15 vs t0
  # Carry out Likelihood ratio tests
  lrt_pool_DMSOvst0 <- glmLRT(z_pool_fit, contrast=my.contrasts[,"DMSOvst0"])
  lrt_pool_PDSvst0 <- glmLRT(z_pool_fit, contrast=my.contrasts[,"PDSvst0"])
  lrt_pool_PhenDC3vst0 <- glmLRT(z_pool_fit, contrast=my.contrasts[,"PhenDC3vst0"])
  #lrt_pool_PDSvsDMSO <- glmLRT(z_pool_fit, contrast=my.contrasts[,"PDSvsDMSO"])
  #lrt_pool_PhenDC3vsDMSO <- glmLRT(z_pool_fit, contrast=my.contrasts[,"PhenDC3vsDMSO"])
  #lrt_pool_PhenDC3vsPDS <- glmLRT(z_pool_fit, contrast=my.contrasts[,"PhenDC3vsPDS"])
  #lrt_pool_PDSandPhenDC3vst0 <- glmLRT(z_pool_fit, contrast=my.contrasts[,"PDSandPhenDC3vst0"])
  #lrt_pool_PDSandPhenDC3vsDMSO <- glmLRT(z_pool_fit, contrast=my.contrasts[,"PDSandPhenDC3vsDMSO"])
  #lrt_pool_PDSandPhenDC3andDMSOvst0 <- glmLRT(z_pool_fit, contrast=my.contrasts[,"PDSandPhenDC3andDMSOvst0"])

  # Create folder to produce the output
  system(sprintf("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/pools/%s", p))

  # Get ranked shRNAS, combine shRNAs into a table for all pools and write FDR and LogFC(t15/t0) tables for DMSOvst0, PDSvst0 and PhenDC3vst0
  # DMSOvst0
  top_pool_DMSOvst0 <- topTags(lrt_pool_DMSOvst0, n=Inf)
  top_DMSOvst0 <- rbind(top_DMSOvst0, data.table(top_pool_DMSOvst0$table)[order(FDR),.(shRNA, gene, pool, FDR, logFC)])
  write.table(data.table(top_pool_DMSOvst0$table)[order(FDR),.(shRNA, FDR, logFC)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_DMSOvst0_FDR_logFC.txt", quote = FALSE, sep = "\t", row.names=FALSE)
  system(sprintf("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_DMSOvst0_FDR_logFC.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/pools/%s", p))
  # PDSvst0
  top_pool_PDSvst0 <- topTags(lrt_pool_PDSvst0, n=Inf)
  top_PDSvst0 <- rbind(top_PDSvst0, data.table(top_pool_PDSvst0$table)[order(FDR),.(shRNA, gene, pool, FDR, logFC)])
  write.table(data.table(top_pool_PDSvst0$table)[order(FDR),.(shRNA, FDR, logFC)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_PDSvst0_FDR_logFC.txt", quote = FALSE, sep = "\t", row.names=FALSE)
  system(sprintf("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_PDSvst0_FDR_logFC.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/pools/%s", p))
  # PhenDC3vst0
  top_pool_PhenDC3vst0 <- topTags(lrt_pool_PhenDC3vst0, n=Inf)
  top_PhenDC3vst0 <- rbind(top_PhenDC3vst0, data.table(top_pool_PhenDC3vst0$table)[order(FDR),.(shRNA, gene, pool, FDR, logFC)])
  write.table(data.table(top_pool_PhenDC3vst0$table)[order(FDR),.(shRNA, FDR, logFC)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_PhenDC3vst0_FDR_logFC.txt", quote = FALSE, sep = "\t", row.names=FALSE)
  system(sprintf("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_PhenDC3vst0_FDR_logFC.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/pools/%s", p))
}


# Create directory structure for figures and tables below
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/venn_diagrams")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/venn_diagrams/all")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/venn_diagrams/pass_50_threshold")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/tables")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/tables/all")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/tables/only")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/violin_plot")


# Apply FDR < 1e-3 & (logFC < -1 | logFC > 1) and set vectors containing significant shRNAs for all pools for DMSOvst0, PDSvst0 and PhenDC3vst0
# DMSOvst0
topids_DMSOvst0 <- as.character(top_DMSOvst0[FDR < 1e-3 & (logFC < -1 | logFC > 1)]$shRNA)
# PDSvst0
topids_PDSvst0 <- as.character(top_PDSvst0[FDR < 1e-3 & (logFC < -1 | logFC > 1)]$shRNA)
# PhenDC3vst0
topids_PhenDC3vst0 <- as.character(top_PhenDC3vst0[FDR < 1e-3 & (logFC < -1 | logFC > 1)]$shRNA)


# Venn diagram combining significant shRNAs for all pools
venn.plot <- draw.triple.venn(
  area1 = length(topids_PDSvst0),
  area2 = length(topids_PhenDC3vst0),
  area3 = length(topids_DMSOvst0),
  n12 = length(intersect(topids_PDSvst0, topids_PhenDC3vst0)),
  n23 = length(intersect(topids_PhenDC3vst0, topids_DMSOvst0)),
  n13 = length(intersect(topids_PDSvst0, topids_DMSOvst0)),
  n123 = length(intersect(intersect(topids_PDSvst0, topids_PhenDC3vst0), topids_DMSOvst0)),
  category = c(sprintf("PDS (%i)", length(topids_PDSvst0)), sprintf("PhenDC3 (%i)", length(topids_PhenDC3vst0)), sprintf("DMSO (%i)", length(topids_DMSOvst0))),
  fill = c("red3", "forestgreen", "deepskyblue3"),
  cex = 1.5,
  fontfamily = "sans",
  cat.dist = c(0.08, 0.08, 0.04),
  cat.cex = 1.5,
  cat.col = c("red3", "forestgreen", "deepskyblue3"),
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  print.mode = c("raw", "percent"),
  sigdigs = 2,
  margin = 0.075)

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161102_t15vst0_shRNA_allpools_venn_FDR10minus3_logFC1minus1.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161102_t15vst0_shRNA_allpools_venn_FDR10minus3_logFC1minus1.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/venn_diagrams/all")


# Apply FDR < 1e-3 & (logFC < -1 | logFC > 1) and from the shRNA tables above, set vectors containing significant genes for all pools for DMSOvst0, PDSvst0 and PhenDC3vst0
# DMSOvst0
topids_genes_DMSOvst0 <- unique(as.character(top_DMSOvst0[FDR < 1e-3 & (logFC < -1 | logFC > 1)]$gene))
# PDSvst0
topids_genes_PDSvst0 <- unique(as.character(top_PDSvst0[FDR < 1e-3 & (logFC < -1 | logFC > 1)]$gene))
# PhenDC3vst0
topids_genes_PhenDC3vst0 <- unique(as.character(top_PhenDC3vst0[FDR < 1e-3 & (logFC < -1 | logFC > 1)]$gene))


# Venn diagram combining significant genes for all pools
venn.plot <- draw.triple.venn(
  area1 = length(topids_genes_PDSvst0),
  area2 = length(topids_genes_PhenDC3vst0),
  area3 = length(topids_genes_DMSOvst0),
  n12 = length(intersect(topids_genes_PDSvst0, topids_genes_PhenDC3vst0)),
  n23 = length(intersect(topids_genes_PhenDC3vst0, topids_genes_DMSOvst0)),
  n13 = length(intersect(topids_genes_PDSvst0, topids_genes_DMSOvst0)),
  n123 = length(intersect(intersect(topids_genes_PDSvst0, topids_genes_PhenDC3vst0), topids_genes_DMSOvst0)),
  category = c(sprintf("PDS (%i)", length(topids_genes_PDSvst0)), sprintf("PhenDC3 (%i)", length(topids_genes_PhenDC3vst0)), sprintf("DMSO (%i)", length(topids_genes_DMSOvst0))),
  fill = c("red3", "forestgreen", "deepskyblue3"),
  cex = 1.5,
  fontfamily = "sans",
  cat.dist = c(0.08, 0.08, 0.04),
  cat.cex = 1.5,
  cat.col = c("red3", "forestgreen", "deepskyblue3"),
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  print.mode = c("raw", "percent"),
  sigdigs = 2,
  margin = 0.075)

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161102_t15vst0_genes_allpools_venn_FDR10minus3_logFC1minus1.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161102_t15vst0_genes_allpools_venn_FDR10minus3_logFC1minus1.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/venn_diagrams/all")



# Tables of genes by DMSOvst0, PDSvst0 and PhenDC3vst0 and for resistant (logFC>1) or sensitiser (logFC>-1)
# I checked but there is no logFC=0

# All genes
f <- readDNAStringSet("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3.fa")
all_shRNAs <- names(f)
all_genes <- data.table(table(as.character(sapply(all_shRNAs, function(x) unlist(strsplit(x, "__"))[1])))) # all genes, 18937

# DMSOvst0
detected_genes_DMSOvst0 <- top_DMSOvst0[,.(n_detected_shRNA = uniqueN(shRNA)), by = gene] # detected genes DMSOvst0, 18883, 54 genes with no shRNAs detected, shRNAs filtered are not included in this table either

all_detected_genes_DMSOvst0 <- merge(all_genes, detected_genes_DMSOvst0, by.x="V1", by.y="gene", all=T)
setnames(all_detected_genes_DMSOvst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs"))

significant_resistant_genes_DMSOvst0 <- top_DMSOvst0[FDR < 1e-3 & logFC > 1,.(n_significant_resistant_shRNAs = uniqueN(shRNA), logFC_significant_resistant_shRNAs = mean(logFC)), by = gene] # significant and resistant genes DMSOvst0, 969
significant_sensitiser_genes_DMSOvst0 <- top_DMSOvst0[FDR < 1e-3 & logFC < -1,.(n_significant_sensitiser_shRNAs = uniqueN(shRNA), logFC_significant_sensitiser_shRNAs = mean(logFC)), by = gene] # significant and sensitiser genes DMSOvst0, 996
length(intersect(significant_resistant_genes_DMSOvst0$gene, significant_sensitiser_genes_DMSOvst0$gene)) # 25, e.g. KDM3B

all_detected_significant_resistant_genes_DMSOvst0 <- merge(all_detected_genes_DMSOvst0, significant_resistant_genes_DMSOvst0, by.x="gene", by.y="gene", all=T)
all_detected_significant_resistant_sensitiser_genes_DMSOvst0 <- merge(all_detected_significant_resistant_genes_DMSOvst0, significant_sensitiser_genes_DMSOvst0, by.x="gene", by.y="gene", all=T)

for (j in c("n_detected_shRNAs", "n_significant_resistant_shRNAs", "n_significant_sensitiser_shRNAs")){
  set(all_detected_significant_resistant_sensitiser_genes_DMSOvst0, which(is.na(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[[j]])),j,0)
}

all_detected_significant_resistant_sensitiser_genes_DMSOvst0[, pct_significant_resistant_shRNAs := ifelse(n_detected_shRNAs > 0, 100*n_significant_resistant_shRNAs/n_detected_shRNAs, 0)]
all_detected_significant_resistant_sensitiser_genes_DMSOvst0[, pass_50_threshold_significant_resistant_shRNAs := ifelse(pct_significant_resistant_shRNAs >= 50, "yes", "no")]
all_detected_significant_resistant_sensitiser_genes_DMSOvst0[, pct_significant_sensitiser_shRNAs := ifelse(n_detected_shRNAs > 0, 100*n_significant_sensitiser_shRNAs/n_detected_shRNAs, 0)]
all_detected_significant_resistant_sensitiser_genes_DMSOvst0[, pass_50_threshold_significant_sensitiser_shRNAs := ifelse(pct_significant_sensitiser_shRNAs >= 50, "yes", "no")]

setcolorder(all_detected_significant_resistant_sensitiser_genes_DMSOvst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs", "n_significant_resistant_shRNAs", "logFC_significant_resistant_shRNAs", "pct_significant_resistant_shRNAs", "pass_50_threshold_significant_resistant_shRNAs", "n_significant_sensitiser_shRNAs", "logFC_significant_sensitiser_shRNAs", "pct_significant_sensitiser_shRNAs", "pass_50_threshold_significant_sensitiser_shRNAs"))

all_detected_significant_resistant_sensitiser_genes_DMSOvst0["KDM3B"]
#gene n_total_shRNAs n_detected_shRNAs n_significant_resistant_shRNAs logFC_significant_resistant_shRNAs pct_significant_resistant_shRNAs pass_50_threshold_significant_resistant_shRNAs n_significant_sensitiser_shRNAs
#1: KDM3B             18                17                              1                           1.797318                         5.882353                                             no                               1
#logFC_significant_sensitiser_shRNAs pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs
#1:                           -3.774646                          5.882353                                              no

nrow(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_resistant_shRNAs == "yes"]) # 42 genes are resistant in DMSOvst0, have more than 50% of shRNAS detected with FDR < 1e-3 & logFC > 1
nrow(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"]) # 47 genes are sensitiser in DMSOvst0, have more than 50% of shRNAS detected with FDR < 1e-3 & logFC < -1
nrow(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_resistant_shRNAs == "yes" & pass_50_threshold_significant_sensitiser_shRNAs == "yes"]) # 0 genes are both resistant and sensitiser in DMSOvst0

write.table(all_detected_significant_resistant_sensitiser_genes_DMSOvst0, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_DMSOvst0_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_DMSOvst0_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/tables/all")


# PDSvst0
detected_genes_PDSvst0 <- top_PDSvst0[,.(n_detected_shRNA = uniqueN(shRNA)), by = gene] # detected genes PDSvst0, 18883, 54 genes with no shRNAs detected, shRNAs filtered are not included in this table either

all_detected_genes_PDSvst0 <- merge(all_genes, detected_genes_PDSvst0, by.x="V1", by.y="gene", all=T)
setnames(all_detected_genes_PDSvst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs"))

significant_resistant_genes_PDSvst0 <- top_PDSvst0[FDR < 1e-3 & logFC > 1,.(n_significant_resistant_shRNAs = uniqueN(shRNA), logFC_significant_resistant_shRNAs = mean(logFC)), by = gene] # significant and resistant genes PDSvst0, 1378
significant_sensitiser_genes_PDSvst0 <- top_PDSvst0[FDR < 1e-3 & logFC < -1,.(n_significant_sensitiser_shRNAs = uniqueN(shRNA), logFC_significant_sensitiser_shRNAs = mean(logFC)), by = gene] # significant and sensitiser genes PDSvst0, 1388
length(intersect(significant_resistant_genes_PDSvst0$gene, significant_sensitiser_genes_PDSvst0$gene)) # 71, e.g. KDM3B

all_detected_significant_resistant_genes_PDSvst0 <- merge(all_detected_genes_PDSvst0, significant_resistant_genes_PDSvst0, by.x="gene", by.y="gene", all=T)
all_detected_significant_resistant_sensitiser_genes_PDSvst0 <- merge(all_detected_significant_resistant_genes_PDSvst0, significant_sensitiser_genes_PDSvst0, by.x="gene", by.y="gene", all=T)

for (j in c("n_detected_shRNAs", "n_significant_resistant_shRNAs", "n_significant_sensitiser_shRNAs")){
  set(all_detected_significant_resistant_sensitiser_genes_PDSvst0, which(is.na(all_detected_significant_resistant_sensitiser_genes_PDSvst0[[j]])),j,0)
}

all_detected_significant_resistant_sensitiser_genes_PDSvst0[, pct_significant_resistant_shRNAs := ifelse(n_detected_shRNAs > 0, 100*n_significant_resistant_shRNAs/n_detected_shRNAs, 0)]
all_detected_significant_resistant_sensitiser_genes_PDSvst0[, pass_50_threshold_significant_resistant_shRNAs := ifelse(pct_significant_resistant_shRNAs >= 50, "yes", "no")]
all_detected_significant_resistant_sensitiser_genes_PDSvst0[, pct_significant_sensitiser_shRNAs := ifelse(n_detected_shRNAs > 0, 100*n_significant_sensitiser_shRNAs/n_detected_shRNAs, 0)]
all_detected_significant_resistant_sensitiser_genes_PDSvst0[, pass_50_threshold_significant_sensitiser_shRNAs := ifelse(pct_significant_sensitiser_shRNAs >= 50, "yes", "no")]

setcolorder(all_detected_significant_resistant_sensitiser_genes_PDSvst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs", "n_significant_resistant_shRNAs", "logFC_significant_resistant_shRNAs", "pct_significant_resistant_shRNAs", "pass_50_threshold_significant_resistant_shRNAs", "n_significant_sensitiser_shRNAs", "logFC_significant_sensitiser_shRNAs", "pct_significant_sensitiser_shRNAs", "pass_50_threshold_significant_sensitiser_shRNAs"))

all_detected_significant_resistant_sensitiser_genes_PDSvst0["KDM3B"]
#gene n_total_shRNAs n_detected_shRNAs n_significant_resistant_shRNAs logFC_significant_resistant_shRNAs pct_significant_resistant_shRNAs pass_50_threshold_significant_resistant_shRNAs n_significant_sensitiser_shRNAs
#1: KDM3B             18                17                              2                           1.384025                         11.76471                                             no                               1
#logFC_significant_sensitiser_shRNAs pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs
#1:                           -3.350787                          5.882353                                              no

all_detected_significant_resistant_sensitiser_genes_PDSvst0["DHX36"]
#    gene n_total_shRNAs n_detected_shRNAs n_significant_resistant_shRNAs logFC_significant_resistant_shRNAs pct_significant_resistant_shRNAs pass_50_threshold_significant_resistant_shRNAs n_significant_sensitiser_shRNAs
#1: DHX36              8                 5                              0                                 NA                                0                                             no                               0
#   logFC_significant_sensitiser_shRNAs pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs
#1:                                  NA                                 0                                              no

nrow(all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes"]) # 65 genes are resistant in PDSvst0, have more than 50% of shRNAS detected with FDR < 1e-3 & logFC > 1
nrow(all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"]) # 59 genes are sensitiser in PDSvst0, have more than 50% of shRNAS detected with FDR < 1e-3 & logFC < -1
nrow(all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes" & pass_50_threshold_significant_sensitiser_shRNAs == "yes"]) # 0 genes are both resistant and sensitiser in PDSvst0

write.table(all_detected_significant_resistant_sensitiser_genes_PDSvst0, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_PDSvst0_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_PDSvst0_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/tables/all")


# PhenDC3vst0
detected_genes_PhenDC3vst0 <- top_PhenDC3vst0[,.(n_detected_shRNA = uniqueN(shRNA)), by = gene] # detected genes PhenDC3vst0, 18883, 54 genes with no shRNAs detected, shRNAs filtered are not included in this table either

all_detected_genes_PhenDC3vst0 <- merge(all_genes, detected_genes_PhenDC3vst0, by.x="V1", by.y="gene", all=T)
setnames(all_detected_genes_PhenDC3vst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs"))

significant_resistant_genes_PhenDC3vst0 <- top_PhenDC3vst0[FDR < 1e-3 & logFC > 1,.(n_significant_resistant_shRNAs = uniqueN(shRNA), logFC_significant_resistant_shRNAs = mean(logFC)), by = gene] # significant and resistant genes PhenDC3vst0, 1859
significant_sensitiser_genes_PhenDC3vst0 <- top_PhenDC3vst0[FDR < 1e-3 & logFC < -1,.(n_significant_sensitiser_shRNAs = uniqueN(shRNA), logFC_significant_sensitiser_shRNAs = mean(logFC)), by = gene] # significant and sensitiser genes PhenDC3vst0, 2474
length(intersect(significant_resistant_genes_PhenDC3vst0$gene, significant_sensitiser_genes_PhenDC3vst0$gene)) # 177, e.g. KDM3B

all_detected_significant_resistant_genes_PhenDC3vst0 <- merge(all_detected_genes_PhenDC3vst0, significant_resistant_genes_PhenDC3vst0, by.x="gene", by.y="gene", all=T)
all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0 <- merge(all_detected_significant_resistant_genes_PhenDC3vst0, significant_sensitiser_genes_PhenDC3vst0, by.x="gene", by.y="gene", all=T)

for (j in c("n_detected_shRNAs", "n_significant_resistant_shRNAs", "n_significant_sensitiser_shRNAs")){
  set(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0, which(is.na(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[[j]])),j,0)
}

all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[, pct_significant_resistant_shRNAs := ifelse(n_detected_shRNAs > 0, 100*n_significant_resistant_shRNAs/n_detected_shRNAs, 0)]
all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[, pass_50_threshold_significant_resistant_shRNAs := ifelse(pct_significant_resistant_shRNAs >= 50, "yes", "no")]
all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[, pct_significant_sensitiser_shRNAs := ifelse(n_detected_shRNAs > 0, 100*n_significant_sensitiser_shRNAs/n_detected_shRNAs, 0)]
all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[, pass_50_threshold_significant_sensitiser_shRNAs := ifelse(pct_significant_sensitiser_shRNAs >= 50, "yes", "no")]

setcolorder(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs", "n_significant_resistant_shRNAs", "logFC_significant_resistant_shRNAs", "pct_significant_resistant_shRNAs", "pass_50_threshold_significant_resistant_shRNAs", "n_significant_sensitiser_shRNAs", "logFC_significant_sensitiser_shRNAs", "pct_significant_sensitiser_shRNAs", "pass_50_threshold_significant_sensitiser_shRNAs"))

all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0["KDM3B"]
#gene n_total_shRNAs n_detected_shRNAs n_significant_resistant_shRNAs logFC_significant_resistant_shRNAs pct_significant_resistant_shRNAs pass_50_threshold_significant_resistant_shRNAs n_significant_sensitiser_shRNAs
#1: KDM3B             18                17                              2                           1.645628                         11.76471                                             no                               1
#logFC_significant_sensitiser_shRNAs pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs
#1:                           -4.401803                          5.882353                                              no

all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0["DHX36"]
#gene n_total_shRNAs n_detected_shRNAs n_significant_resistant_shRNAs logFC_significant_resistant_shRNAs pct_significant_resistant_shRNAs pass_50_threshold_significant_resistant_shRNAs n_significant_sensitiser_shRNAs
#1: DHX36              8                 5                              0                                 NA                                0                                             no                               4
#logFC_significant_sensitiser_shRNAs pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs
#1:                            -4.73586                                80                                             yes

nrow(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes"]) # 85 genes are resistant in PhenDC3vst0, have more than 50% of shRNAS detected with FDR < 1e-3 & logFC > 1
nrow(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"]) # 170 genes are sensitiser in PhenDC3vst0, have more than 50% of shRNAS detected with FDR < 1e-3 & logFC < -1
nrow(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes" & pass_50_threshold_significant_sensitiser_shRNAs == "yes"]) # 1 gene (OR5D14, actually an olfatory receptor) is both resistant and sensitiser in PhenDC3vst0 at the same time

write.table(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_PhenDC3vst0_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_PhenDC3vst0_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/tables/all")




# Venn diagram resistant

# Apply 50% of shRNAS cutoff and set vectors containing significant genes for all pools for DMSOvst0, PDSvst0 and PhenDC3vst0
# DMSOvst0
topids_genes_pass_50_threshold_significant_resistant_shRNAs_DMSOvst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_resistant_shRNAs == "yes"]$gene))
# PDSvst0
topids_genes_pass_50_threshold_significant_resistant_shRNAs_PDSvst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes"]$gene))
# PhenDC3vst0
topids_genes_pass_50_threshold_significant_resistant_shRNAs_PhenDC3vst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes"]$gene))

venn.plot <- draw.triple.venn(
  area1 = length(topids_genes_pass_50_threshold_significant_resistant_shRNAs_PDSvst0),
  area2 = length(topids_genes_pass_50_threshold_significant_resistant_shRNAs_PhenDC3vst0),
  area3 = length(topids_genes_pass_50_threshold_significant_resistant_shRNAs_DMSOvst0),
  n12 = length(intersect(topids_genes_pass_50_threshold_significant_resistant_shRNAs_PDSvst0, topids_genes_pass_50_threshold_significant_resistant_shRNAs_PhenDC3vst0)),
  n23 = length(intersect(topids_genes_pass_50_threshold_significant_resistant_shRNAs_PhenDC3vst0, topids_genes_pass_50_threshold_significant_resistant_shRNAs_DMSOvst0)),
  n13 = length(intersect(topids_genes_pass_50_threshold_significant_resistant_shRNAs_PDSvst0, topids_genes_pass_50_threshold_significant_resistant_shRNAs_DMSOvst0)),
  n123 = length(intersect(intersect(topids_genes_pass_50_threshold_significant_resistant_shRNAs_PDSvst0, topids_genes_pass_50_threshold_significant_resistant_shRNAs_PhenDC3vst0), topids_genes_pass_50_threshold_significant_resistant_shRNAs_DMSOvst0)),
  category = c(sprintf("PDS (%i)", length(topids_genes_pass_50_threshold_significant_resistant_shRNAs_PDSvst0)), sprintf("PhenDC3 (%i)", length(topids_genes_pass_50_threshold_significant_resistant_shRNAs_PhenDC3vst0)), sprintf("DMSO (%i)", length(topids_genes_pass_50_threshold_significant_resistant_shRNAs_DMSOvst0))),
  fill = c("red3", "forestgreen", "deepskyblue3"),
  cex = 1.5,
  fontfamily = "sans",
  cat.dist = c(0.08, 0.08, 0.04),
  cat.cex = 1.5,
  cat.col = c("red3", "forestgreen", "deepskyblue3"),
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  print.mode = c("raw", "percent"),
  sigdigs = 2,
  margin = 0.075)

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161102_t15vst0_genes_allpools_venn_FDR10minus3_logFC1minus1_pass_50_threshold_significant_resistant_shRNAs.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161102_t15vst0_genes_allpools_venn_FDR10minus3_logFC1minus1_pass_50_threshold_significant_resistant_shRNAs.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/venn_diagrams/pass_50_threshold")




# Venn diagram sensitiser

# Apply 50% of shRNAS cutoff and set vectors containing significant genes for all pools for DMSOvst0, PDSvst0 and PhenDC3vst0
# DMSOvst0
topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_DMSOvst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"]$gene))
# PDSvst0
topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PDSvst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"]$gene))
# PhenDC3vst0
topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PhenDC3vst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"]$gene))

venn.plot <- draw.triple.venn(
  area1 = length(topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PDSvst0),
  area2 = length(topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PhenDC3vst0),
  area3 = length(topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_DMSOvst0),
  n12 = length(intersect(topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PDSvst0, topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PhenDC3vst0)),
  n23 = length(intersect(topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PhenDC3vst0, topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_DMSOvst0)),
  n13 = length(intersect(topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PDSvst0, topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_DMSOvst0)),
  n123 = length(intersect(intersect(topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PDSvst0, topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PhenDC3vst0), topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_DMSOvst0)),
  category = c(sprintf("PDS (%i)", length(topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PDSvst0)), sprintf("PhenDC3 (%i)", length(topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PhenDC3vst0)), sprintf("DMSO (%i)", length(topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_DMSOvst0))),
  fill = c("red3", "forestgreen", "deepskyblue3"),
  cex = 1.5,
  fontfamily = "sans",
  cat.dist = c(0.08, 0.08, 0.04),
  cat.cex = 1.5,
  cat.col = c("red3", "forestgreen", "deepskyblue3"),
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  print.mode = c("raw", "percent"),
  sigdigs = 2,
  margin = 0.075)

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161102_t15vst0_genes_allpools_venn_FDR10minus3_logFC1minus1_pass_50_threshold_significant_sensitiser_shRNAs.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161102_t15vst0_genes_allpools_venn_FDR10minus3_logFC1minus1_pass_50_threshold_significant_sensitiser_shRNAs.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/venn_diagrams/pass_50_threshold")




# Tables of genes for PDS only, PhenDC3 only and PDS+PhenDC3 only with FDR < 1e-3 & logFC < -1 (sensitiser) or logFC > 1 (resistant) & passing the 50% threshold of significant shRNAS

# resistant
# PDS only
top_PDSonly_resistant <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes"][!(gene %in% c(topids_genes_pass_50_threshold_significant_resistant_shRNAs_DMSOvst0, topids_genes_pass_50_threshold_significant_resistant_shRNAs_PhenDC3vst0))][order(-logFC_significant_resistant_shRNAs),.(gene, logFC_significant_resistant_shRNAs, pct_significant_resistant_shRNAs)]
write.table(top_PDSonly_resistant, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_PDSonly_resistant_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_PDSonly_resistant_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/tables/only")

# PhenDC3 only
top_PhenDC3only_resistant <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes"][!(gene %in% c(topids_genes_pass_50_threshold_significant_resistant_shRNAs_DMSOvst0, topids_genes_pass_50_threshold_significant_resistant_shRNAs_PDSvst0))][order(-logFC_significant_resistant_shRNAs),.(gene, logFC_significant_resistant_shRNAs, pct_significant_resistant_shRNAs)]
write.table(top_PhenDC3only_resistant, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_PhenDC3only_resistant_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_PhenDC3only_resistant_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/tables/only")

# PDS+PhenDC3 only
top_PDSandPhenDC3only_PhenDC3_resistant <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes"][gene %in% topids_genes_pass_50_threshold_significant_resistant_shRNAs_PDSvst0][!(gene %in% topids_genes_pass_50_threshold_significant_resistant_shRNAs_DMSOvst0)][order(-logFC_significant_resistant_shRNAs),.(gene, logFC_significant_resistant_shRNAs, pct_significant_resistant_shRNAs)]
top_PDSandPhenDC3only_PDS_resistant <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes"][gene %in% topids_genes_pass_50_threshold_significant_resistant_shRNAs_PhenDC3vst0][!(gene %in% topids_genes_pass_50_threshold_significant_resistant_shRNAs_DMSOvst0)][order(-logFC_significant_resistant_shRNAs),.(gene, logFC_significant_resistant_shRNAs, pct_significant_resistant_shRNAs)]
top_PDSandPhenDC3only_resistant <- merge(top_PDSandPhenDC3only_PhenDC3_resistant, top_PDSandPhenDC3only_PDS_resistant, by = "gene", suffixes = c(".PhenDC3", ".PDS"))[order(-logFC_significant_resistant_shRNAs.PhenDC3)]
write.table(top_PDSandPhenDC3only_resistant, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_PDSandPhenDC3only_resistant_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_PDSandPhenDC3only_resistant_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/tables/only")


# sensitiser
# PDS only
top_PDSonly_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"][!(gene %in% c(topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_DMSOvst0, topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PhenDC3vst0))][order(logFC_significant_sensitiser_shRNAs),.(gene, logFC_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs)]
write.table(top_PDSonly_sensitiser, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_PDSonly_sensitiser_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_PDSonly_sensitiser_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/tables/only")

# PhenDC3 only
top_PhenDC3only_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"][!(gene %in% c(topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_DMSOvst0, topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PDSvst0))][order(logFC_significant_sensitiser_shRNAs),.(gene, logFC_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs)]
write.table(top_PhenDC3only_sensitiser, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_PhenDC3only_sensitiser_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_PhenDC3only_sensitiser_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/tables/only")

# PDS+PhenDC3 only
top_PDSandPhenDC3only_PhenDC3_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"][gene %in% topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PDSvst0][!(gene %in% topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_DMSOvst0)][order(logFC_significant_sensitiser_shRNAs),.(gene, logFC_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs)]
top_PDSandPhenDC3only_PDS_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"][gene %in% topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PhenDC3vst0][!(gene %in% topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_DMSOvst0)][order(logFC_significant_sensitiser_shRNAs),.(gene, logFC_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs)]
top_PDSandPhenDC3only_sensitiser <- merge(top_PDSandPhenDC3only_PhenDC3_sensitiser, top_PDSandPhenDC3only_PDS_sensitiser, by = "gene", suffixes = c(".PhenDC3", ".PDS"))[order(logFC_significant_sensitiser_shRNAs.PhenDC3)]
write.table(top_PDSandPhenDC3only_sensitiser, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_PDSandPhenDC3only_sensitiser_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161102_PDSandPhenDC3only_sensitiser_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/tables/only")




# Plot distributions to determine optimal threshold for each DMSO, PhenDC3 and PDS, and resistant and sensitiser
# resistant
DMSO_resistant <- all_detected_significant_resistant_sensitiser_genes_DMSOvst0[n_significant_resistant_shRNAs > 0]$pct_significant_resistant_shRNAs
DMSO_sensitiser <- all_detected_significant_resistant_sensitiser_genes_DMSOvst0[n_significant_sensitiser_shRNAs > 0]$pct_significant_sensitiser_shRNAs
PDS_resistant <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[n_significant_resistant_shRNAs > 0]$pct_significant_resistant_shRNAs
PDS_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[n_significant_sensitiser_shRNAs > 0]$pct_significant_sensitiser_shRNAs
PhenDC3_resistant <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[n_significant_resistant_shRNAs > 0]$pct_significant_resistant_shRNAs
PhenDC3_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[n_significant_sensitiser_shRNAs > 0]$pct_significant_sensitiser_shRNAs

resistant <- data.table(stack(list(DMSO_resistant = DMSO_resistant, DMSO_sensitiser = DMSO_sensitiser, PDS_resistant = PDS_resistant, PDS_sensitiser = PDS_sensitiser, PhenDC3_resistant = PhenDC3_resistant, PhenDC3_sensitiser = PhenDC3_sensitiser)))

gg <- ggplot(resistant, aes(ind, values)) +
geom_violin() +
xlab("") +
ylab("% significant resistant shRNAs") +
theme_classic()

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161102_violin.pdf", width = 22/2.54, height = 10/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161102_violin.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161102/violin_plot")


```


# TODO

* Use glmTreat (see section 2.12 in edgeR tutorial and Mccarthy2009) to test genes above lfc=1 in absolute terms. This would be development in next rounds of the screen. Also, the IntensityFilter.R approach developed by Dario could be a good approach.


* individual analyses of pools 3, 4 and 7


* Learning:
edge R manual - see final experiment on CRISPR by pools and section 2.14 Gene testing http://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

shRNA tutorial http://bioinf.wehi.edu.au/shRNAseq/pooledScreenAnalysis.pdf


* Tools: roast, mroast, camera, fry, romer


* Suggestions
David: use plotly to visualise the results of the screen dinamically and interactively (scatterplots and venn.diagrams). Perhaps docker might help here too.


* Publications: Sheridan2015, Yu2013


* Relevant people from Sanger (George Vassiliou, Hannes Ponstingl, William Skarnes, Vivek Iyer)


* Other ideas

Setup a contrast to get PDSonly, PhenDC3only and PDS+PhenDC3only

Perhaps setup models for each comparison of conditions individually

We might need to do thresholding based on logFC higher/lower than 1/-1 of the individual shRNAs, e.g. if more than half of the shRNAs have logFC > 1 or < -1, then select it

Another idea would be to just use the genes found in the toptags of the shRNA analysis directly

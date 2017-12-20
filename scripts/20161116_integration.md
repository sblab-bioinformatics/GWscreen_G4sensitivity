
This script is identical to [20161115_integration.md](20161115_integration.md) but adding the following reads from the resequencing run and then repeating two separate analyses:

* From pool 3: t15_PhenDC3_Pool3_1 and t15_PhenDC3_Pool3_3 (first chunk of analysis)
* From pools 11 and 12: t15_DMSO_Pool11_3 and t15_DMSO_Pool12_2 (second chunk of analysis)


# Analysis adding t15_PhenDC3_Pool3_1 and t15_PhenDC3_Pool3_3 to 20161115

```R
library(edgeR)
library(data.table)
library(ggplot2)
library(VennDiagram)
library(reshape)
library(Biostrings)


# Enlarge the view width when printing tables
options(width = 250)


# Load table of counts into edgeR from the original sequencing runs and add the new sequencing runs to them

# data1
data1 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160627_Pool1_counts.txt", header = T)
data1_rerun <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161110/20161108_Pool1_counts.txt", header = T)

nrow(data1) # 9475
nrow(data1_rerun) # 9478
length(intersect(rownames(data1), rownames(data1_rerun))) # 9473

data1_bind <- rbindlist(list(data.table(data1, keep.rownames=TRUE), data.table(data1_rerun, keep.rownames=TRUE)))
data1_sum <- data1_bind[, lapply(.SD, sum), by=rn]
nrow(data1_sum) # Correct, 9480 = 9473 + (9475-9473) + (9478-9473)

data1_combined <- data.frame(data1_sum)[,2:13]
rownames(data1_combined) <- data1_sum$rn
dim(data1_combined) # 9480   12


# data2
data2 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160627_Pool2_counts.txt", header = T)
data2_rerun <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161110/20161108_Pool2_counts.txt", header = T)

nrow(data2) # 9467
nrow(data2_rerun) # 9468
length(intersect(rownames(data2), rownames(data2_rerun))) # 9467

data2_bind <- rbindlist(list(data.table(data2, keep.rownames=TRUE), data.table(data2_rerun, keep.rownames=TRUE)))
data2_sum <- data2_bind[, lapply(.SD, sum), by=rn]
nrow(data2_sum) # Correct, 9468 = 9467 + (9467-9467) + (9468-9467)

data2_combined <- data.frame(data2_sum)[,2:13]
rownames(data2_combined) <- data2_sum$rn
dim(data2_combined) # 9468   12


# data5
data5 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160629_Pool5_counts.txt", header = T)
data5_rerun <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161110/20161109_Pool5_counts.txt", header = T)

nrow(data5) # 8846
nrow(data5_rerun) # 8859
length(intersect(rownames(data5), rownames(data5_rerun))) # 8843

data5_bind <- rbindlist(list(data.table(data5, keep.rownames=TRUE), data.table(data5_rerun, keep.rownames=TRUE)))
data5_sum <- data5_bind[, lapply(.SD, sum), by=rn]
nrow(data5_sum) # Correct, 8862 = 8843 + (8846-8843) + (8859-8843)

data5_combined <- data.frame(data5_sum)[,2:13]
rownames(data5_combined) <- data5_sum$rn
dim(data5_combined) # 8862   12


# data10
data10 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160629_Pool10_counts.txt", header = T)
data10_rerun1 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161110/20161109_Pool10_counts.txt", header = T)

nrow(data10) # 9486
nrow(data10_rerun1) # 9494
length(intersect(rownames(data10), rownames(data10_rerun1))) # 9484

data10_bind1 <- rbindlist(list(data.table(data10, keep.rownames=TRUE), data.table(data10_rerun1, keep.rownames=TRUE)))
data10_sum1 <- data10_bind1[, lapply(.SD, sum), by=rn]
nrow(data10_sum1) # Correct, 9496 = 9484 + (9486-9484) + (9494-9484)

data10_combined1 <- data.frame(data10_sum1)[,2:13]
rownames(data10_combined1) <- data10_sum1$rn
dim(data10_combined1) # 9496   12

data10_rerun2 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161110/20161108_Pool10_counts.txt", header = T)

nrow(data10_combined1) # 9496
nrow(data10_rerun2) # 9448
length(intersect(rownames(data10_combined1), rownames(data10_rerun2))) # 9447

data10_bind2 <- rbindlist(list(data.table(data10_combined1, keep.rownames=TRUE), data.table(data10_rerun2, keep.rownames=TRUE)), use.names = TRUE, fill = TRUE)
for (j in colnames(data10_bind2)){
  set(data10_bind2, which(is.na(data10_bind2[[j]])),j,0) # change all NAs introduced by rbindlist to 0
}
data10_sum2 <- data10_bind2[, lapply(.SD, sum), by=rn]
nrow(data10_sum2) # Correct, 9497 = 9447 + (9496-9447) + (9448-9447)

data10_combined2 <- data.frame(data10_sum2)[,2:13]
rownames(data10_combined2) <- data10_sum2$rn
dim(data10_combined2) # 9497   12


# data6
data6 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160920_Pool6_counts.txt", header = T)
data6_rerun <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161110/20161108_Pool6_counts.txt", header = T)

nrow(data6) # 9544
nrow(data6_rerun) # 9543
length(intersect(rownames(data6), rownames(data6_rerun))) # 9537

data6_bind <- rbindlist(list(data.table(data6, keep.rownames=TRUE), data.table(data6_rerun, keep.rownames=TRUE)), use.names = TRUE, fill = TRUE)
for (j in colnames(data6_bind)){
  set(data6_bind, which(is.na(data6_bind[[j]])),j,0) # change all NAs introduced by rbindlist to 0
}
data6_sum <- data6_bind[, lapply(.SD, sum), by=rn]
nrow(data6_sum) # Correct, 9550 = 9537 + (9544-9537) + (9543-9537)

data6_combined <- data.frame(data6_sum)[,2:13]
rownames(data6_combined) <- data6_sum$rn
dim(data6_combined) # 9550   12


# data3
data3 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160908_Pool3_counts_lims.txt", header = T)
data3_rerun <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161110/20161109_Pool3_counts.txt", header = T) # here is a difference with 20161115_integration.md, we are adding all four columns here

nrow(data3) # 9518
nrow(data3_rerun) # 9476
length(intersect(rownames(data3), rownames(data3_rerun))) # 9470

data3_bind <- rbindlist(list(data.table(data3, keep.rownames=TRUE), data.table(data3_rerun, keep.rownames=TRUE)), use.names = TRUE, fill = TRUE)
for (j in colnames(data3_bind)){
  set(data3_bind, which(is.na(data3_bind[[j]])),j,0) # change all NAs introduced by rbindlist to 0
}
data3_sum <- data3_bind[, lapply(.SD, sum), by=rn]
nrow(data3_sum) # Correct, 9524 = 9470 + (9518-9470) + (9476-9470)

data3_combined <- data.frame(data3_sum)[,2:13]
rownames(data3_combined) <- data3_sum$rn
dim(data3_combined) # 9524   12


# data4
data4 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160908_Pool4_counts_lims.txt", header = T)
data4_rerun <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161110/20161109_Pool4_counts.txt", header = T)

nrow(data4) # 8384
nrow(data4_rerun) # 7995
length(intersect(rownames(data4), rownames(data4_rerun))) # 7950

data4_bind <- rbindlist(list(data.table(data4, keep.rownames=TRUE), data.table(data4_rerun, keep.rownames=TRUE)), use.names = TRUE, fill = TRUE)
for (j in colnames(data4_bind)){
  set(data4_bind, which(is.na(data4_bind[[j]])),j,0) # change all NAs introduced by rbindlist to 0
}
data4_sum <- data4_bind[, lapply(.SD, sum), by=rn]
nrow(data4_sum) # Correct, 8429 = 7950 + (8384-7950) + (7995-7950)

data4_combined <- data.frame(data4_sum)[,2:13]
rownames(data4_combined) <- data4_sum$rn
dim(data4_combined) # 8429   12


# data7, data8, data9, data11, data12
data7 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160920_Pool7_counts.txt", header = T)
data8 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20161011_Pool8_counts.txt", header = T)
data9 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20161011_Pool9_counts.txt", header = T)
data11 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20161018_Pool11_counts.txt", header = T)
data12 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20161018_Pool12_counts.txt", header = T)


# Combine table of counts into a single data structure
colnames(data1_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data2_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data3_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data4_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data5_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data6_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data7) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data8) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data9) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data10_combined2) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data11) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data12) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")

data <- rbind(data1_combined, data2_combined, data3_combined, data4_combined, data5_combined, data6_combined, data7, data8, data9, data10_combined2, data11, data12)
dim(data) # 110666     12

# 100*(110666/113002) = 97.9% shRNAs detected, without adding the new resequencing we had 110576


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
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/pools")


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
  system(sprintf("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/pools/%s", p))

  # Get ranked shRNAS, combine shRNAs into a table for all pools and write FDR and LogFC(t15/t0) tables for DMSOvst0, PDSvst0 and PhenDC3vst0
  # DMSOvst0
  top_pool_DMSOvst0 <- topTags(lrt_pool_DMSOvst0, n=Inf)
  top_DMSOvst0 <- rbind(top_DMSOvst0, data.table(top_pool_DMSOvst0$table)[order(FDR),.(shRNA, gene, pool, FDR, logFC)])
  write.table(data.table(top_pool_DMSOvst0$table)[order(FDR),.(shRNA, FDR, logFC)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_DMSOvst0_FDR_logFC.txt", quote = FALSE, sep = "\t", row.names=FALSE)
  system(sprintf("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_DMSOvst0_FDR_logFC.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/pools/%s", p))
  # PDSvst0
  top_pool_PDSvst0 <- topTags(lrt_pool_PDSvst0, n=Inf)
  top_PDSvst0 <- rbind(top_PDSvst0, data.table(top_pool_PDSvst0$table)[order(FDR),.(shRNA, gene, pool, FDR, logFC)])
  write.table(data.table(top_pool_PDSvst0$table)[order(FDR),.(shRNA, FDR, logFC)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_PDSvst0_FDR_logFC.txt", quote = FALSE, sep = "\t", row.names=FALSE)
  system(sprintf("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_PDSvst0_FDR_logFC.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/pools/%s", p))
  # PhenDC3vst0
  top_pool_PhenDC3vst0 <- topTags(lrt_pool_PhenDC3vst0, n=Inf)
  top_PhenDC3vst0 <- rbind(top_PhenDC3vst0, data.table(top_pool_PhenDC3vst0$table)[order(FDR),.(shRNA, gene, pool, FDR, logFC)])
  write.table(data.table(top_pool_PhenDC3vst0$table)[order(FDR),.(shRNA, FDR, logFC)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_PhenDC3vst0_FDR_logFC.txt", quote = FALSE, sep = "\t", row.names=FALSE)
  system(sprintf("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_PhenDC3vst0_FDR_logFC.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/pools/%s", p))
}

#[1] "Pool1"
#[1] "shRNAs: 9480"
#[1] "shRNAs (after filtering): 9281"
#[1] "Common dispersion: 0.583118325474814"
#[1] "Pool2"
#[1] "shRNAs: 9468"
#[1] "shRNAs (after filtering): 9285"
#[1] "Common dispersion: 0.573277576626521"
#[1] "Pool3"
#[1] "shRNAs: 9524"
#[1] "shRNAs (after filtering): 9419"
#[1] "Common dispersion: 0.409505095179379"
#[1] "Pool4"
#[1] "shRNAs: 8429"
#[1] "shRNAs (after filtering): 7127"
#[1] "Common dispersion: 0.483172144840967"
#[1] "Pool5"
#[1] "shRNAs: 8862"
#[1] "shRNAs (after filtering): 8454"
#[1] "Common dispersion: 0.629206233924443"
#[1] "Pool6"
#[1] "shRNAs: 9550"
#[1] "shRNAs (after filtering): 9464"
#[1] "Common dispersion: 0.358011336192241"
#[1] "Pool7"
#[1] "shRNAs: 9553"
#[1] "shRNAs (after filtering): 9455"
#[1] "Common dispersion: 0.363310395024059"
#[1] "Pool8"
#[1] "shRNAs: 9581"
#[1] "shRNAs (after filtering): 9512"
#[1] "Common dispersion: 0.319725438980772"
#[1] "Pool9"
#[1] "shRNAs: 9158"
#[1] "shRNAs (after filtering): 8748"
#[1] "Common dispersion: 0.371306793344731"
#[1] "Pool10"
#[1] "shRNAs: 9497"
#[1] "shRNAs (after filtering): 9271"
#[1] "Common dispersion: 0.548651321681364"
#[1] "Pool11"
#[1] "shRNAs: 7979"
#[1] "shRNAs (after filtering): 7726"
#[1] "Common dispersion: 0.360802866789428"
#[1] "Pool12"
#[1] "shRNAs: 9585"
#[1] "shRNAs (after filtering): 9376"
#[1] "Common dispersion: 0.394298205548235"


# Create directory structure for figures and tables below
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/venn_diagrams")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/venn_diagrams/all")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/venn_diagrams/pass_50_threshold")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/tables")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/tables/all")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/tables/only")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/violin_plot")


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

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161116_1_t15vst0_shRNA_allpools_venn_FDR10minus3_logFC1minus1.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161116_1_t15vst0_shRNA_allpools_venn_FDR10minus3_logFC1minus1.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/venn_diagrams/all")


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

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161116_1_t15vst0_genes_allpools_venn_FDR10minus3_logFC1minus1.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161116_1_t15vst0_genes_allpools_venn_FDR10minus3_logFC1minus1.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/venn_diagrams/all")



# Tables of genes by DMSOvst0, PDSvst0 and PhenDC3vst0 and for resistant (logFC>1) or sensitiser (logFC>-1)

# All genes
f <- readDNAStringSet("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3.fa")
all_shRNAs <- names(f)
all_genes <- data.table(table(as.character(sapply(all_shRNAs, function(x) unlist(strsplit(x, "__"))[1])))) # all genes, 18937

# DMSOvst0
detected_genes_DMSOvst0 <- top_DMSOvst0[,.(n_detected_shRNA = uniqueN(shRNA)), by = gene] # detected genes DMSOvst0, 18908, 29 genes with no shRNAs detected, shRNAs filtered are not included in this table either

all_detected_genes_DMSOvst0 <- merge(all_genes, detected_genes_DMSOvst0, by.x="V1", by.y="gene", all=T)
setnames(all_detected_genes_DMSOvst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs"))

significant_resistant_genes_DMSOvst0 <- top_DMSOvst0[FDR < 1e-3 & logFC > 1,.(n_significant_resistant_shRNAs = uniqueN(shRNA), logFC_significant_resistant_shRNAs = mean(logFC)), by = gene] # significant and resistant genes DMSOvst0, 1064
significant_sensitiser_genes_DMSOvst0 <- top_DMSOvst0[FDR < 1e-3 & logFC < -1,.(n_significant_sensitiser_shRNAs = uniqueN(shRNA), logFC_significant_sensitiser_shRNAs = mean(logFC)), by = gene] # significant and sensitiser genes DMSOvst0, 1018
length(intersect(significant_resistant_genes_DMSOvst0$gene, significant_sensitiser_genes_DMSOvst0$gene)) # 28, e.g. KDM3B

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
#1: KDM3B             18                17                              1                           2.505529                         5.882353                                             no                               1
#logFC_significant_sensitiser_shRNAs pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs
#1:                           -3.695722                          5.882353                                              no

nrow(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_resistant_shRNAs == "yes"]) # 30 genes are resistant in DMSOvst0, have more than 50% of shRNAS detected with FDR < 1e-3 & logFC > 1
nrow(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"]) # 30 genes are sensitiser in DMSOvst0, have more than 50% of shRNAS detected with FDR < 1e-3 & logFC < -1
nrow(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_resistant_shRNAs == "yes" & pass_50_threshold_significant_sensitiser_shRNAs == "yes"]) # 0 genes are both resistant and sensitiser in DMSOvst0

write.table(all_detected_significant_resistant_sensitiser_genes_DMSOvst0, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_DMSOvst0_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_DMSOvst0_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/tables/all")


# PDSvst0
detected_genes_PDSvst0 <- top_PDSvst0[,.(n_detected_shRNA = uniqueN(shRNA)), by = gene] # detected genes PDSvst0, 18908, 29 genes with no shRNAs detected, shRNAs filtered are not included in this table either

all_detected_genes_PDSvst0 <- merge(all_genes, detected_genes_PDSvst0, by.x="V1", by.y="gene", all=T)
setnames(all_detected_genes_PDSvst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs"))

significant_resistant_genes_PDSvst0 <- top_PDSvst0[FDR < 1e-3 & logFC > 1,.(n_significant_resistant_shRNAs = uniqueN(shRNA), logFC_significant_resistant_shRNAs = mean(logFC)), by = gene] # significant and resistant genes PDSvst0, 1456
significant_sensitiser_genes_PDSvst0 <- top_PDSvst0[FDR < 1e-3 & logFC < -1,.(n_significant_sensitiser_shRNAs = uniqueN(shRNA), logFC_significant_sensitiser_shRNAs = mean(logFC)), by = gene] # significant and sensitiser genes PDSvst0, 1431
length(intersect(significant_resistant_genes_PDSvst0$gene, significant_sensitiser_genes_PDSvst0$gene)) # 78, e.g. KDM3B

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
#1: KDM3B             18                17                              2                           1.742746                         11.76471                                             no                               1
#logFC_significant_sensitiser_shRNAs pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs
#1:                             -3.2964                          5.882353                                              no

all_detected_significant_resistant_sensitiser_genes_PDSvst0["DHX36"]
#gene n_total_shRNAs n_detected_shRNAs n_significant_resistant_shRNAs logFC_significant_resistant_shRNAs pct_significant_resistant_shRNAs pass_50_threshold_significant_resistant_shRNAs n_significant_sensitiser_shRNAs
#1: DHX36              8                 5                              0                                 NA                                0                                             no                               0
#logFC_significant_sensitiser_shRNAs pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs
#1:                                  NA                                 0                                              no

nrow(all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes"]) # 47 genes are resistant in PDSvst0, have more than 50% of shRNAS detected with FDR < 1e-3 & logFC > 1
nrow(all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"]) # 45 genes are sensitiser in PDSvst0, have more than 50% of shRNAS detected with FDR < 1e-3 & logFC < -1
nrow(all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes" & pass_50_threshold_significant_sensitiser_shRNAs == "yes"]) # 0 genes are both resistant and sensitiser in PDSvst0

write.table(all_detected_significant_resistant_sensitiser_genes_PDSvst0, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_PDSvst0_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_PDSvst0_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/tables/all")


# PhenDC3vst0
detected_genes_PhenDC3vst0 <- top_PhenDC3vst0[,.(n_detected_shRNA = uniqueN(shRNA)), by = gene] # detected genes PhenDC3vst0, 18908, 29 genes with no shRNAs detected, shRNAs filtered are not included in this table either

all_detected_genes_PhenDC3vst0 <- merge(all_genes, detected_genes_PhenDC3vst0, by.x="V1", by.y="gene", all=T)
setnames(all_detected_genes_PhenDC3vst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs"))

significant_resistant_genes_PhenDC3vst0 <- top_PhenDC3vst0[FDR < 1e-3 & logFC > 1,.(n_significant_resistant_shRNAs = uniqueN(shRNA), logFC_significant_resistant_shRNAs = mean(logFC)), by = gene] # significant and resistant genes PhenDC3vst0, 1881
significant_sensitiser_genes_PhenDC3vst0 <- top_PhenDC3vst0[FDR < 1e-3 & logFC < -1,.(n_significant_sensitiser_shRNAs = uniqueN(shRNA), logFC_significant_sensitiser_shRNAs = mean(logFC)), by = gene] # significant and sensitiser genes PhenDC3vst0, 2545
length(intersect(significant_resistant_genes_PhenDC3vst0$gene, significant_sensitiser_genes_PhenDC3vst0$gene)) # 195, e.g. KDM3B

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
#1: KDM3B             18                17                              2                           1.681956                         11.76471                                             no                               1
#logFC_significant_sensitiser_shRNAs pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs
#1:                           -4.795345                          5.882353                                              no

all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0["DHX36"]
#gene n_total_shRNAs n_detected_shRNAs n_significant_resistant_shRNAs logFC_significant_resistant_shRNAs pct_significant_resistant_shRNAs pass_50_threshold_significant_resistant_shRNAs n_significant_sensitiser_shRNAs
#1: DHX36              8                 5                              0                                 NA                                0                                             no                               4
#logFC_significant_sensitiser_shRNAs pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs
#1:                           -4.950305                                80                                             yes

nrow(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes"]) # 63 genes are resistant in PhenDC3vst0, have more than 50% of shRNAS detected with FDR < 1e-3 & logFC > 1
nrow(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"]) # 122 genes are sensitiser in PhenDC3vst0, have more than 50% of shRNAS detected with FDR < 1e-3 & logFC < -1
nrow(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes" & pass_50_threshold_significant_sensitiser_shRNAs == "yes"]) # 0 genes is both resistant and sensitiser in PhenDC3vst0 at the same time

write.table(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_PhenDC3vst0_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_PhenDC3vst0_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/tables/all")



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

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161116_1_t15vst0_genes_allpools_venn_FDR10minus3_logFC1minus1_pass_50_threshold_significant_resistant_shRNAs.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161116_1_t15vst0_genes_allpools_venn_FDR10minus3_logFC1minus1_pass_50_threshold_significant_resistant_shRNAs.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/venn_diagrams/pass_50_threshold")



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

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161116_1_t15vst0_genes_allpools_venn_FDR10minus3_logFC1minus1_pass_50_threshold_significant_sensitiser_shRNAs.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161116_1_t15vst0_genes_allpools_venn_FDR10minus3_logFC1minus1_pass_50_threshold_significant_sensitiser_shRNAs.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/venn_diagrams/pass_50_threshold")




# Tables of genes for PDS only, PhenDC3 only and PDS+PhenDC3 only with FDR < 1e-3 & logFC < -1 (sensitiser) or logFC > 1 (resistant) & passing the 50% threshold of significant shRNAS

# resistant
# PDS only
top_PDSonly_resistant <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes"][!(gene %in% c(topids_genes_pass_50_threshold_significant_resistant_shRNAs_DMSOvst0, topids_genes_pass_50_threshold_significant_resistant_shRNAs_PhenDC3vst0))][order(-logFC_significant_resistant_shRNAs),.(gene, logFC_significant_resistant_shRNAs, pct_significant_resistant_shRNAs)]
write.table(top_PDSonly_resistant, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_PDSonly_resistant_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_PDSonly_resistant_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/tables/only")

# PhenDC3 only
top_PhenDC3only_resistant <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes"][!(gene %in% c(topids_genes_pass_50_threshold_significant_resistant_shRNAs_DMSOvst0, topids_genes_pass_50_threshold_significant_resistant_shRNAs_PDSvst0))][order(-logFC_significant_resistant_shRNAs),.(gene, logFC_significant_resistant_shRNAs, pct_significant_resistant_shRNAs)]
write.table(top_PhenDC3only_resistant, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_PhenDC3only_resistant_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_PhenDC3only_resistant_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/tables/only")

# PDS+PhenDC3 only
top_PDSandPhenDC3only_PhenDC3_resistant <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes"][gene %in% topids_genes_pass_50_threshold_significant_resistant_shRNAs_PDSvst0][!(gene %in% topids_genes_pass_50_threshold_significant_resistant_shRNAs_DMSOvst0)][order(-logFC_significant_resistant_shRNAs),.(gene, logFC_significant_resistant_shRNAs, pct_significant_resistant_shRNAs)]
top_PDSandPhenDC3only_PDS_resistant <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes"][gene %in% topids_genes_pass_50_threshold_significant_resistant_shRNAs_PhenDC3vst0][!(gene %in% topids_genes_pass_50_threshold_significant_resistant_shRNAs_DMSOvst0)][order(-logFC_significant_resistant_shRNAs),.(gene, logFC_significant_resistant_shRNAs, pct_significant_resistant_shRNAs)]
top_PDSandPhenDC3only_resistant <- merge(top_PDSandPhenDC3only_PhenDC3_resistant, top_PDSandPhenDC3only_PDS_resistant, by = "gene", suffixes = c(".PhenDC3", ".PDS"))[order(-logFC_significant_resistant_shRNAs.PhenDC3)]
write.table(top_PDSandPhenDC3only_resistant, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_PDSandPhenDC3only_resistant_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_PDSandPhenDC3only_resistant_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/tables/only")


# sensitiser
# PDS only
top_PDSonly_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"][!(gene %in% c(topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_DMSOvst0, topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PhenDC3vst0))][order(logFC_significant_sensitiser_shRNAs),.(gene, logFC_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs)]
write.table(top_PDSonly_sensitiser, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_PDSonly_sensitiser_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_PDSonly_sensitiser_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/tables/only")

# PhenDC3 only
top_PhenDC3only_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"][!(gene %in% c(topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_DMSOvst0, topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PDSvst0))][order(logFC_significant_sensitiser_shRNAs),.(gene, logFC_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs)]
write.table(top_PhenDC3only_sensitiser, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_PhenDC3only_sensitiser_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_PhenDC3only_sensitiser_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/tables/only")

# PDS+PhenDC3 only
top_PDSandPhenDC3only_PhenDC3_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"][gene %in% topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PDSvst0][!(gene %in% topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_DMSOvst0)][order(logFC_significant_sensitiser_shRNAs),.(gene, logFC_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs)]
top_PDSandPhenDC3only_PDS_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"][gene %in% topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PhenDC3vst0][!(gene %in% topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_DMSOvst0)][order(logFC_significant_sensitiser_shRNAs),.(gene, logFC_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs)]
top_PDSandPhenDC3only_sensitiser <- merge(top_PDSandPhenDC3only_PhenDC3_sensitiser, top_PDSandPhenDC3only_PDS_sensitiser, by = "gene", suffixes = c(".PhenDC3", ".PDS"))[order(logFC_significant_sensitiser_shRNAs.PhenDC3)]
write.table(top_PDSandPhenDC3only_sensitiser, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_PDSandPhenDC3only_sensitiser_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_1_PDSandPhenDC3only_sensitiser_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/tables/only")



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

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161116_1_violin.pdf", width = 22/2.54, height = 10/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161116_1_violin.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_1/violin_plot")

```



# Analysis adding t15_DMSO_Pool11_3 and t15_DMSO_Pool12_2 to above

```R
library(edgeR)
library(data.table)
library(ggplot2)
library(VennDiagram)
library(reshape)
library(Biostrings)


# Enlarge the view width when printing tables
options(width = 250)


# Load table of counts into edgeR from the original sequencing runs and add the new sequencing runs to them

# data1
data1 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160627_Pool1_counts.txt", header = T)
data1_rerun <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161110/20161108_Pool1_counts.txt", header = T)

nrow(data1) # 9475
nrow(data1_rerun) # 9478
length(intersect(rownames(data1), rownames(data1_rerun))) # 9473

data1_bind <- rbindlist(list(data.table(data1, keep.rownames=TRUE), data.table(data1_rerun, keep.rownames=TRUE)))
data1_sum <- data1_bind[, lapply(.SD, sum), by=rn]
nrow(data1_sum) # Correct, 9480 = 9473 + (9475-9473) + (9478-9473)

data1_combined <- data.frame(data1_sum)[,2:13]
rownames(data1_combined) <- data1_sum$rn
dim(data1_combined) # 9480   12


# data2
data2 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160627_Pool2_counts.txt", header = T)
data2_rerun <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161110/20161108_Pool2_counts.txt", header = T)

nrow(data2) # 9467
nrow(data2_rerun) # 9468
length(intersect(rownames(data2), rownames(data2_rerun))) # 9467

data2_bind <- rbindlist(list(data.table(data2, keep.rownames=TRUE), data.table(data2_rerun, keep.rownames=TRUE)))
data2_sum <- data2_bind[, lapply(.SD, sum), by=rn]
nrow(data2_sum) # Correct, 9468 = 9467 + (9467-9467) + (9468-9467)

data2_combined <- data.frame(data2_sum)[,2:13]
rownames(data2_combined) <- data2_sum$rn
dim(data2_combined) # 9468   12


# data5
data5 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160629_Pool5_counts.txt", header = T)
data5_rerun <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161110/20161109_Pool5_counts.txt", header = T)

nrow(data5) # 8846
nrow(data5_rerun) # 8859
length(intersect(rownames(data5), rownames(data5_rerun))) # 8843

data5_bind <- rbindlist(list(data.table(data5, keep.rownames=TRUE), data.table(data5_rerun, keep.rownames=TRUE)))
data5_sum <- data5_bind[, lapply(.SD, sum), by=rn]
nrow(data5_sum) # Correct, 8862 = 8843 + (8846-8843) + (8859-8843)

data5_combined <- data.frame(data5_sum)[,2:13]
rownames(data5_combined) <- data5_sum$rn
dim(data5_combined) # 8862   12


# data10
data10 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160629_Pool10_counts.txt", header = T)
data10_rerun1 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161110/20161109_Pool10_counts.txt", header = T)

nrow(data10) # 9486
nrow(data10_rerun1) # 9494
length(intersect(rownames(data10), rownames(data10_rerun1))) # 9484

data10_bind1 <- rbindlist(list(data.table(data10, keep.rownames=TRUE), data.table(data10_rerun1, keep.rownames=TRUE)))
data10_sum1 <- data10_bind1[, lapply(.SD, sum), by=rn]
nrow(data10_sum1) # Correct, 9496 = 9484 + (9486-9484) + (9494-9484)

data10_combined1 <- data.frame(data10_sum1)[,2:13]
rownames(data10_combined1) <- data10_sum1$rn
dim(data10_combined1) # 9496   12

data10_rerun2 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161110/20161108_Pool10_counts.txt", header = T)

nrow(data10_combined1) # 9496
nrow(data10_rerun2) # 9448
length(intersect(rownames(data10_combined1), rownames(data10_rerun2))) # 9447

data10_bind2 <- rbindlist(list(data.table(data10_combined1, keep.rownames=TRUE), data.table(data10_rerun2, keep.rownames=TRUE)), use.names = TRUE, fill = TRUE)
for (j in colnames(data10_bind2)){
  set(data10_bind2, which(is.na(data10_bind2[[j]])),j,0) # change all NAs introduced by rbindlist to 0
}
data10_sum2 <- data10_bind2[, lapply(.SD, sum), by=rn]
nrow(data10_sum2) # Correct, 9497 = 9447 + (9496-9447) + (9448-9447)

data10_combined2 <- data.frame(data10_sum2)[,2:13]
rownames(data10_combined2) <- data10_sum2$rn
dim(data10_combined2) # 9497   12


# data6
data6 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160920_Pool6_counts.txt", header = T)
data6_rerun <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161110/20161108_Pool6_counts.txt", header = T)

nrow(data6) # 9544
nrow(data6_rerun) # 9543
length(intersect(rownames(data6), rownames(data6_rerun))) # 9537

data6_bind <- rbindlist(list(data.table(data6, keep.rownames=TRUE), data.table(data6_rerun, keep.rownames=TRUE)), use.names = TRUE, fill = TRUE)
for (j in colnames(data6_bind)){
  set(data6_bind, which(is.na(data6_bind[[j]])),j,0) # change all NAs introduced by rbindlist to 0
}
data6_sum <- data6_bind[, lapply(.SD, sum), by=rn]
nrow(data6_sum) # Correct, 9550 = 9537 + (9544-9537) + (9543-9537)

data6_combined <- data.frame(data6_sum)[,2:13]
rownames(data6_combined) <- data6_sum$rn
dim(data6_combined) # 9550   12


# data3
data3 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160908_Pool3_counts_lims.txt", header = T)
data3_rerun <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161110/20161109_Pool3_counts.txt", header = T) # here is a difference with 20161115_integration.md, we are adding all four columns here

nrow(data3) # 9518
nrow(data3_rerun) # 9476
length(intersect(rownames(data3), rownames(data3_rerun))) # 9470

data3_bind <- rbindlist(list(data.table(data3, keep.rownames=TRUE), data.table(data3_rerun, keep.rownames=TRUE)), use.names = TRUE, fill = TRUE)
for (j in colnames(data3_bind)){
  set(data3_bind, which(is.na(data3_bind[[j]])),j,0) # change all NAs introduced by rbindlist to 0
}
data3_sum <- data3_bind[, lapply(.SD, sum), by=rn]
nrow(data3_sum) # Correct, 9524 = 9470 + (9518-9470) + (9476-9470)

data3_combined <- data.frame(data3_sum)[,2:13]
rownames(data3_combined) <- data3_sum$rn
dim(data3_combined) # 9524   12


# data4
data4 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160908_Pool4_counts_lims.txt", header = T)
data4_rerun <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161110/20161109_Pool4_counts.txt", header = T)

nrow(data4) # 8384
nrow(data4_rerun) # 7995
length(intersect(rownames(data4), rownames(data4_rerun))) # 7950

data4_bind <- rbindlist(list(data.table(data4, keep.rownames=TRUE), data.table(data4_rerun, keep.rownames=TRUE)), use.names = TRUE, fill = TRUE)
for (j in colnames(data4_bind)){
  set(data4_bind, which(is.na(data4_bind[[j]])),j,0) # change all NAs introduced by rbindlist to 0
}
data4_sum <- data4_bind[, lapply(.SD, sum), by=rn]
nrow(data4_sum) # Correct, 8429 = 7950 + (8384-7950) + (7995-7950)

data4_combined <- data.frame(data4_sum)[,2:13]
rownames(data4_combined) <- data4_sum$rn
dim(data4_combined) # 8429   12


# data11
data11 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20161018_Pool11_counts.txt", header = T)
data11_rerun <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161110/20161109_Pool11_counts.txt", header = T)

nrow(data11) # 7979
nrow(data11_rerun) # 7981
length(intersect(rownames(data11), rownames(data11_rerun))) # 7964

data11_bind <- rbindlist(list(data.table(data11, keep.rownames=TRUE), data.table(data11_rerun, keep.rownames=TRUE)), use.names = TRUE, fill = TRUE)
for (j in colnames(data11_bind)){
  set(data11_bind, which(is.na(data11_bind[[j]])),j,0) # change all NAs introduced by rbindlist to 0
}
data11_sum <- data11_bind[, lapply(.SD, sum), by=rn]
nrow(data11_sum) # Correct, 7996 = 7964 + (7979-7964) + (7981-7964)

data11_combined <- data.frame(data11_sum)[,2:13]
rownames(data11_combined) <- data11_sum$rn
dim(data11_combined) # 7996   12


# data12
data12 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20161018_Pool12_counts.txt", header = T)
data12_rerun <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161110/20161108_Pool12_counts.txt", header = T)

nrow(data12) # 9585
nrow(data12_rerun) # 9582
length(intersect(rownames(data12), rownames(data12_rerun))) # 9581

data12_bind <- rbindlist(list(data.table(data12, keep.rownames=TRUE), data.table(data12_rerun, keep.rownames=TRUE)), use.names = TRUE, fill = TRUE)
for (j in colnames(data12_bind)){
  set(data12_bind, which(is.na(data12_bind[[j]])),j,0) # change all NAs introduced by rbindlist to 0
}
data12_sum <- data12_bind[, lapply(.SD, sum), by=rn]
nrow(data12_sum) # Correct, 9586 = 9581 + (9585-9581) + (9582-9581)

data12_combined <- data.frame(data12_sum)[,2:13]
rownames(data12_combined) <- data12_sum$rn
dim(data12_combined) # 9586   12


# data7, data8, data9
data7 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160920_Pool7_counts.txt", header = T)
data8 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20161011_Pool8_counts.txt", header = T)
data9 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20161011_Pool9_counts.txt", header = T)


# Combine table of counts into a single data structure
colnames(data1_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data2_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data3_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data4_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data5_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data6_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data7) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data8) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data9) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data10_combined2) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data11_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data12_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")

data <- rbind(data1_combined, data2_combined, data3_combined, data4_combined, data5_combined, data6_combined, data7, data8, data9, data10_combined2, data11_combined, data12_combined)
dim(data) # 110684     12

# 100*(110684/113002) = 97.9% shRNAs detected, without adding the new resequencing we had 110576


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
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/pools")


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
  system(sprintf("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/pools/%s", p))

  # Get ranked shRNAS, combine shRNAs into a table for all pools and write FDR and LogFC(t15/t0) tables for DMSOvst0, PDSvst0 and PhenDC3vst0
  # DMSOvst0
  top_pool_DMSOvst0 <- topTags(lrt_pool_DMSOvst0, n=Inf)
  top_DMSOvst0 <- rbind(top_DMSOvst0, data.table(top_pool_DMSOvst0$table)[order(FDR),.(shRNA, gene, pool, FDR, logFC)])
  write.table(data.table(top_pool_DMSOvst0$table)[order(FDR),.(shRNA, FDR, logFC)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_DMSOvst0_FDR_logFC.txt", quote = FALSE, sep = "\t", row.names=FALSE)
  system(sprintf("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_DMSOvst0_FDR_logFC.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/pools/%s", p))
  # PDSvst0
  top_pool_PDSvst0 <- topTags(lrt_pool_PDSvst0, n=Inf)
  top_PDSvst0 <- rbind(top_PDSvst0, data.table(top_pool_PDSvst0$table)[order(FDR),.(shRNA, gene, pool, FDR, logFC)])
  write.table(data.table(top_pool_PDSvst0$table)[order(FDR),.(shRNA, FDR, logFC)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_PDSvst0_FDR_logFC.txt", quote = FALSE, sep = "\t", row.names=FALSE)
  system(sprintf("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_PDSvst0_FDR_logFC.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/pools/%s", p))
  # PhenDC3vst0
  top_pool_PhenDC3vst0 <- topTags(lrt_pool_PhenDC3vst0, n=Inf)
  top_PhenDC3vst0 <- rbind(top_PhenDC3vst0, data.table(top_pool_PhenDC3vst0$table)[order(FDR),.(shRNA, gene, pool, FDR, logFC)])
  write.table(data.table(top_pool_PhenDC3vst0$table)[order(FDR),.(shRNA, FDR, logFC)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_PhenDC3vst0_FDR_logFC.txt", quote = FALSE, sep = "\t", row.names=FALSE)
  system(sprintf("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_PhenDC3vst0_FDR_logFC.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/pools/%s", p))
}

#[1] "Pool1"
#[1] "shRNAs: 9480"
#[1] "shRNAs (after filtering): 9281"
#[1] "Common dispersion: 0.583118325474814"
#[1] "Pool2"
#[1] "shRNAs: 9468"
#[1] "shRNAs (after filtering): 9285"
#[1] "Common dispersion: 0.573277576626521"
#[1] "Pool3"
#[1] "shRNAs: 9524"
#[1] "shRNAs (after filtering): 9419"
#[1] "Common dispersion: 0.409505095179379"
#[1] "Pool4"
#[1] "shRNAs: 8429"
#[1] "shRNAs (after filtering): 7127"
#[1] "Common dispersion: 0.483172144840967"
#[1] "Pool5"
#[1] "shRNAs: 8862"
#[1] "shRNAs (after filtering): 8454"
#[1] "Common dispersion: 0.629206233924443"
#[1] "Pool6"
#[1] "shRNAs: 9550"
#[1] "shRNAs (after filtering): 9464"
#[1] "Common dispersion: 0.358011336192241"
#[1] "Pool7"
#[1] "shRNAs: 9553"
#[1] "shRNAs (after filtering): 9455"
#[1] "Common dispersion: 0.363310395024059"
#[1] "Pool8"
#[1] "shRNAs: 9581"
#[1] "shRNAs (after filtering): 9512"
#[1] "Common dispersion: 0.319725438980772"
#[1] "Pool9"
#[1] "shRNAs: 9158"
#[1] "shRNAs (after filtering): 8748"
#[1] "Common dispersion: 0.371306793344731"
#[1] "Pool10"
#[1] "shRNAs: 9497"
#[1] "shRNAs (after filtering): 9271"
#[1] "Common dispersion: 0.548651321681364"
#[1] "Pool11"
#[1] "shRNAs: 7996"
#[1] "shRNAs (after filtering): 7726"
#[1] "Common dispersion: 0.385319689976575"
#[1] "Pool12"
#[1] "shRNAs: 9586"
#[1] "shRNAs (after filtering): 9376"
#[1] "Common dispersion: 0.383673380291764"


# Create directory structure for figures and tables below
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/venn_diagrams")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/venn_diagrams/all")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/venn_diagrams/pass_50_threshold")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/tables")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/tables/all")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/tables/only")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/violin_plot")


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

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161116_2_t15vst0_shRNA_allpools_venn_FDR10minus3_logFC1minus1.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161116_2_t15vst0_shRNA_allpools_venn_FDR10minus3_logFC1minus1.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/venn_diagrams/all")


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

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161116_2_t15vst0_genes_allpools_venn_FDR10minus3_logFC1minus1.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161116_2_t15vst0_genes_allpools_venn_FDR10minus3_logFC1minus1.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/venn_diagrams/all")




# Tables of genes by DMSOvst0, PDSvst0 and PhenDC3vst0 and for resistant (logFC>1) or sensitiser (logFC>-1)

# All genes
f <- readDNAStringSet("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3.fa")
all_shRNAs <- names(f)
all_genes <- data.table(table(as.character(sapply(all_shRNAs, function(x) unlist(strsplit(x, "__"))[1])))) # all genes, 18937

# DMSOvst0
detected_genes_DMSOvst0 <- top_DMSOvst0[,.(n_detected_shRNA = uniqueN(shRNA)), by = gene] # detected genes DMSOvst0, 18908, 29 genes with no shRNAs detected, shRNAs filtered are not included in this table either

all_detected_genes_DMSOvst0 <- merge(all_genes, detected_genes_DMSOvst0, by.x="V1", by.y="gene", all=T)
setnames(all_detected_genes_DMSOvst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs"))

significant_resistant_genes_DMSOvst0 <- top_DMSOvst0[FDR < 1e-3 & logFC > 1,.(n_significant_resistant_shRNAs = uniqueN(shRNA), logFC_significant_resistant_shRNAs = mean(logFC)), by = gene] # significant and resistant genes DMSOvst0, 929
significant_sensitiser_genes_DMSOvst0 <- top_DMSOvst0[FDR < 1e-3 & logFC < -1,.(n_significant_sensitiser_shRNAs = uniqueN(shRNA), logFC_significant_sensitiser_shRNAs = mean(logFC)), by = gene] # significant and sensitiser genes DMSOvst0, 801
length(intersect(significant_resistant_genes_DMSOvst0$gene, significant_sensitiser_genes_DMSOvst0$gene)) # 18, e.g. KDM3B

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
#1: KDM3B             18                17                              1                           2.505529                         5.882353                                             no                               1
#logFC_significant_sensitiser_shRNAs pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs
#1:                           -3.695722                          5.882353                                              no

nrow(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_resistant_shRNAs == "yes"]) # 30 genes are resistant in DMSOvst0, have more than 50% of shRNAS detected with FDR < 1e-3 & logFC > 1
nrow(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"]) # 23 genes are sensitiser in DMSOvst0, have more than 50% of shRNAS detected with FDR < 1e-3 & logFC < -1
nrow(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_resistant_shRNAs == "yes" & pass_50_threshold_significant_sensitiser_shRNAs == "yes"]) # 0 genes are both resistant and sensitiser in DMSOvst0

write.table(all_detected_significant_resistant_sensitiser_genes_DMSOvst0, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_DMSOvst0_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_DMSOvst0_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/tables/all")


# PDSvst0
detected_genes_PDSvst0 <- top_PDSvst0[,.(n_detected_shRNA = uniqueN(shRNA)), by = gene] # detected genes PDSvst0, 18908, 29 genes with no shRNAs detected, shRNAs filtered are not included in this table either

all_detected_genes_PDSvst0 <- merge(all_genes, detected_genes_PDSvst0, by.x="V1", by.y="gene", all=T)
setnames(all_detected_genes_PDSvst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs"))

significant_resistant_genes_PDSvst0 <- top_PDSvst0[FDR < 1e-3 & logFC > 1,.(n_significant_resistant_shRNAs = uniqueN(shRNA), logFC_significant_resistant_shRNAs = mean(logFC)), by = gene] # significant and resistant genes PDSvst0, 1355
significant_sensitiser_genes_PDSvst0 <- top_PDSvst0[FDR < 1e-3 & logFC < -1,.(n_significant_sensitiser_shRNAs = uniqueN(shRNA), logFC_significant_sensitiser_shRNAs = mean(logFC)), by = gene] # significant and sensitiser genes PDSvst0, 1379
length(intersect(significant_resistant_genes_PDSvst0$gene, significant_sensitiser_genes_PDSvst0$gene)) # 70, e.g. KDM3B

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
#1: KDM3B             18                17                              1                            2.45804                         5.882353                                             no                               1
#logFC_significant_sensitiser_shRNAs pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs
#1:                             -3.2964                          5.882353                                              no

all_detected_significant_resistant_sensitiser_genes_PDSvst0["DHX36"]
#gene n_total_shRNAs n_detected_shRNAs n_significant_resistant_shRNAs logFC_significant_resistant_shRNAs pct_significant_resistant_shRNAs pass_50_threshold_significant_resistant_shRNAs n_significant_sensitiser_shRNAs
#1: DHX36              8                 5                              0                                 NA                                0                                             no                               0
#logFC_significant_sensitiser_shRNAs pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs
#1:                                  NA                                 0                                              no

nrow(all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes"]) # 46 genes are resistant in PDSvst0, have more than 50% of shRNAS detected with FDR < 1e-3 & logFC > 1
nrow(all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"]) # 41 genes are sensitiser in PDSvst0, have more than 50% of shRNAS detected with FDR < 1e-3 & logFC < -1
nrow(all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes" & pass_50_threshold_significant_sensitiser_shRNAs == "yes"]) # 0 genes are both resistant and sensitiser in PDSvst0

write.table(all_detected_significant_resistant_sensitiser_genes_PDSvst0, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_PDSvst0_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_PDSvst0_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/tables/all")


# PhenDC3vst0
detected_genes_PhenDC3vst0 <- top_PhenDC3vst0[,.(n_detected_shRNA = uniqueN(shRNA)), by = gene] # detected genes PhenDC3vst0, 18908, 29 genes with no shRNAs detected, shRNAs filtered are not included in this table either

all_detected_genes_PhenDC3vst0 <- merge(all_genes, detected_genes_PhenDC3vst0, by.x="V1", by.y="gene", all=T)
setnames(all_detected_genes_PhenDC3vst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs"))

significant_resistant_genes_PhenDC3vst0 <- top_PhenDC3vst0[FDR < 1e-3 & logFC > 1,.(n_significant_resistant_shRNAs = uniqueN(shRNA), logFC_significant_resistant_shRNAs = mean(logFC)), by = gene] # significant and resistant genes PhenDC3vst0, 1826
significant_sensitiser_genes_PhenDC3vst0 <- top_PhenDC3vst0[FDR < 1e-3 & logFC < -1,.(n_significant_sensitiser_shRNAs = uniqueN(shRNA), logFC_significant_sensitiser_shRNAs = mean(logFC)), by = gene] # significant and sensitiser genes PhenDC3vst0, 2496
length(intersect(significant_resistant_genes_PhenDC3vst0$gene, significant_sensitiser_genes_PhenDC3vst0$gene)) # 192, e.g. KDM3B

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
#1: KDM3B             18                17                              2                           1.681785                         11.76471                                             no                               1
#logFC_significant_sensitiser_shRNAs pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs
#1:                           -4.795345                          5.882353                                              no

all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0["DHX36"]
#gene n_total_shRNAs n_detected_shRNAs n_significant_resistant_shRNAs logFC_significant_resistant_shRNAs pct_significant_resistant_shRNAs pass_50_threshold_significant_resistant_shRNAs n_significant_sensitiser_shRNAs
#1: DHX36              8                 5                              0                                 NA                                0                                             no                               4
#logFC_significant_sensitiser_shRNAs pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs
#1:                           -4.950305                                80                                             yes

nrow(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes"]) # 63 genes are resistant in PhenDC3vst0, have more than 50% of shRNAS detected with FDR < 1e-3 & logFC > 1
nrow(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"]) # 120 genes are sensitiser in PhenDC3vst0, have more than 50% of shRNAS detected with FDR < 1e-3 & logFC < -1
nrow(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes" & pass_50_threshold_significant_sensitiser_shRNAs == "yes"]) # 0 genes is both resistant and sensitiser in PhenDC3vst0 at the same time

write.table(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_PhenDC3vst0_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_PhenDC3vst0_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/tables/all")



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

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161116_2_t15vst0_genes_allpools_venn_FDR10minus3_logFC1minus1_pass_50_threshold_significant_resistant_shRNAs.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161116_2_t15vst0_genes_allpools_venn_FDR10minus3_logFC1minus1_pass_50_threshold_significant_resistant_shRNAs.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/venn_diagrams/pass_50_threshold")



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

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161116_2_t15vst0_genes_allpools_venn_FDR10minus3_logFC1minus1_pass_50_threshold_significant_sensitiser_shRNAs.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161116_2_t15vst0_genes_allpools_venn_FDR10minus3_logFC1minus1_pass_50_threshold_significant_sensitiser_shRNAs.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/venn_diagrams/pass_50_threshold")




# Tables of genes for PDS only, PhenDC3 only and PDS+PhenDC3 only with FDR < 1e-3 & logFC < -1 (sensitiser) or logFC > 1 (resistant) & passing the 50% threshold of significant shRNAS

# resistant
# PDS only
top_PDSonly_resistant <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes"][!(gene %in% c(topids_genes_pass_50_threshold_significant_resistant_shRNAs_DMSOvst0, topids_genes_pass_50_threshold_significant_resistant_shRNAs_PhenDC3vst0))][order(-logFC_significant_resistant_shRNAs),.(gene, logFC_significant_resistant_shRNAs, pct_significant_resistant_shRNAs)]
write.table(top_PDSonly_resistant, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_PDSonly_resistant_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_PDSonly_resistant_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/tables/only")

# PhenDC3 only
top_PhenDC3only_resistant <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes"][!(gene %in% c(topids_genes_pass_50_threshold_significant_resistant_shRNAs_DMSOvst0, topids_genes_pass_50_threshold_significant_resistant_shRNAs_PDSvst0))][order(-logFC_significant_resistant_shRNAs),.(gene, logFC_significant_resistant_shRNAs, pct_significant_resistant_shRNAs)]
write.table(top_PhenDC3only_resistant, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_PhenDC3only_resistant_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_PhenDC3only_resistant_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/tables/only")

# PDS+PhenDC3 only
top_PDSandPhenDC3only_PhenDC3_resistant <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes"][gene %in% topids_genes_pass_50_threshold_significant_resistant_shRNAs_PDSvst0][!(gene %in% topids_genes_pass_50_threshold_significant_resistant_shRNAs_DMSOvst0)][order(-logFC_significant_resistant_shRNAs),.(gene, logFC_significant_resistant_shRNAs, pct_significant_resistant_shRNAs)]
top_PDSandPhenDC3only_PDS_resistant <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes"][gene %in% topids_genes_pass_50_threshold_significant_resistant_shRNAs_PhenDC3vst0][!(gene %in% topids_genes_pass_50_threshold_significant_resistant_shRNAs_DMSOvst0)][order(-logFC_significant_resistant_shRNAs),.(gene, logFC_significant_resistant_shRNAs, pct_significant_resistant_shRNAs)]
top_PDSandPhenDC3only_resistant <- merge(top_PDSandPhenDC3only_PhenDC3_resistant, top_PDSandPhenDC3only_PDS_resistant, by = "gene", suffixes = c(".PhenDC3", ".PDS"))[order(-logFC_significant_resistant_shRNAs.PhenDC3)]
write.table(top_PDSandPhenDC3only_resistant, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_PDSandPhenDC3only_resistant_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_PDSandPhenDC3only_resistant_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/tables/only")


# sensitiser
# PDS only
top_PDSonly_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"][!(gene %in% c(topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_DMSOvst0, topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PhenDC3vst0))][order(logFC_significant_sensitiser_shRNAs),.(gene, logFC_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs)]
write.table(top_PDSonly_sensitiser, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_PDSonly_sensitiser_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_PDSonly_sensitiser_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/tables/only")

# PhenDC3 only
top_PhenDC3only_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"][!(gene %in% c(topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_DMSOvst0, topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PDSvst0))][order(logFC_significant_sensitiser_shRNAs),.(gene, logFC_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs)]
write.table(top_PhenDC3only_sensitiser, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_PhenDC3only_sensitiser_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_PhenDC3only_sensitiser_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/tables/only")

# PDS+PhenDC3 only
top_PDSandPhenDC3only_PhenDC3_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"][gene %in% topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PDSvst0][!(gene %in% topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_DMSOvst0)][order(logFC_significant_sensitiser_shRNAs),.(gene, logFC_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs)]
top_PDSandPhenDC3only_PDS_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes"][gene %in% topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_PhenDC3vst0][!(gene %in% topids_genes_pass_50_threshold_significant_sensitiser_shRNAs_DMSOvst0)][order(logFC_significant_sensitiser_shRNAs),.(gene, logFC_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs)]
top_PDSandPhenDC3only_sensitiser <- merge(top_PDSandPhenDC3only_PhenDC3_sensitiser, top_PDSandPhenDC3only_PDS_sensitiser, by = "gene", suffixes = c(".PhenDC3", ".PDS"))[order(logFC_significant_sensitiser_shRNAs.PhenDC3)]
write.table(top_PDSandPhenDC3only_sensitiser, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_PDSandPhenDC3only_sensitiser_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161116_2_PDSandPhenDC3only_sensitiser_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/tables/only")



# Plot distributions to determine optimal threshold for each DMSO, PhenDC3 and PDS, and resistant and sensitiser
DMSO_resistant <- all_detected_significant_resistant_sensitiser_genes_DMSOvst0[n_significant_resistant_shRNAs > 0]$pct_significant_resistant_shRNAs
DMSO_sensitiser <- all_detected_significant_resistant_sensitiser_genes_DMSOvst0[n_significant_sensitiser_shRNAs > 0]$pct_significant_sensitiser_shRNAs
PDS_resistant <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[n_significant_resistant_shRNAs > 0]$pct_significant_resistant_shRNAs
PDS_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[n_significant_sensitiser_shRNAs > 0]$pct_significant_sensitiser_shRNAs
PhenDC3_resistant <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[n_significant_resistant_shRNAs > 0]$pct_significant_resistant_shRNAs
PhenDC3_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[n_significant_sensitiser_shRNAs > 0]$pct_significant_sensitiser_shRNAs

violin <- data.table(stack(list(DMSO_resistant = DMSO_resistant, DMSO_sensitiser = DMSO_sensitiser, PDS_resistant = PDS_resistant, PDS_sensitiser = PDS_sensitiser, PhenDC3_resistant = PhenDC3_resistant, PhenDC3_sensitiser = PhenDC3_sensitiser)))

gg <- ggplot(violin, aes(ind, values)) +
geom_violin() +
xlab("") +
ylab("% significant resistant shRNAs") +
theme_classic()

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161116_2_violin.pdf", width = 22/2.54, height = 10/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161116_2_violin.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161116_2/violin_plot")




```


This script is identical to [20161206_integration.md](20161206_integration.md) but also borrows from the second half of [20161116_integration.md](20161116_integration.md) and adds the following reads from the second resequencing run from 20170105:

* Pool 3: t15_PhenDC3_Pool3_1 and t15_PhenDC3_Pool3_3
* Pools 11 and 12: t15_DMSO_Pool11_3 and t15_DMSO_Pool12_2

The idea is to put together the full dataset (with this, all re-sequencing runs will be complete), analyse it as we did before and compare it to what we did in [20161206_integration.md](20161206_integration.md).

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


# Load table of counts into edgeR from the original sequencing runs and add the first and second sequencing runs to them

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
data3_rerun1 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161110/20161109_Pool3_counts.txt", header = T) # but I only need t0_no_Pool3_2 and t0_no_Pool3_3 from this

nrow(data3) # 9518
nrow(data3_rerun1) # 9476
length(intersect(rownames(data3), rownames(data3_rerun1))) # 9470

data3_bind1 <- rbindlist(list(data.table(data3, keep.rownames=TRUE), data.table(data3_rerun1[,1:2], keep.rownames=TRUE)), use.names = TRUE, fill = TRUE)
for (j in colnames(data3_bind1)){
  set(data3_bind1, which(is.na(data3_bind1[[j]])),j,0) # change all NAs introduced by rbindlist to 0
}
data3_sum1 <- data3_bind1[, lapply(.SD, sum), by=rn]
nrow(data3_sum1) # Correct, 9524 = 9470 + (9518-9470) + (9476-9470)

data3_combined1 <- data.frame(data3_sum1)[,2:13]
rownames(data3_combined1) <- data3_sum1$rn
dim(data3_combined1) # 9524   12

data3_rerun2 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170118/20170105_Pool3_counts.txt", header = T)

nrow(data3_combined1) # 9524
nrow(data3_rerun2) # 9492
length(intersect(rownames(data3_combined1), rownames(data3_rerun2))) # 9491

data3_bind2 <- rbindlist(list(data.table(data3_combined1, keep.rownames=TRUE), data.table(data3_rerun2, keep.rownames=TRUE)), use.names = TRUE, fill = TRUE)
for (j in colnames(data3_bind2)){
  set(data3_bind2, which(is.na(data3_bind2[[j]])),j,0) # change all NAs introduced by rbindlist to 0
}
data3_sum2 <- data3_bind2[, lapply(.SD, sum), by=rn]
nrow(data3_sum2) # Correct, 9525 = 9491 + (9524-9491) + (9492-9491)

data3_combined2 <- data.frame(data3_sum2)[,2:13]
rownames(data3_combined2) <- data3_sum2$rn
dim(data3_combined2) # 9525   12


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
data11_rerun <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170118/20170105_Pool11_counts.txt", header = T)

nrow(data11) # 7979
nrow(data11_rerun) # 7914
length(intersect(rownames(data11), rownames(data11_rerun))) # 7912

data11_bind <- rbindlist(list(data.table(data11, keep.rownames=TRUE), data.table(data11_rerun, keep.rownames=TRUE)), use.names = TRUE, fill = TRUE)
for (j in colnames(data11_bind)){
  set(data11_bind, which(is.na(data11_bind[[j]])),j,0) # change all NAs introduced by rbindlist to 0
}
data11_sum <- data11_bind[, lapply(.SD, sum), by=rn]
nrow(data11_sum) # Correct, 7981 = 7912 + (7979-7912) + (7914-7912)

data11_combined <- data.frame(data11_sum)[,2:13]
rownames(data11_combined) <- data11_sum$rn
dim(data11_combined) # 7981   12


# data12
data12 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20161018_Pool12_counts.txt", header = T)
data12_rerun <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170118/20170105_Pool12_counts.txt", header = T)

nrow(data12) # 9585
nrow(data12_rerun) # 9546
length(intersect(rownames(data12), rownames(data12_rerun))) # 9546

data12_bind <- rbindlist(list(data.table(data12, keep.rownames=TRUE), data.table(data12_rerun, keep.rownames=TRUE)), use.names = TRUE, fill = TRUE)
for (j in colnames(data12_bind)){
  set(data12_bind, which(is.na(data12_bind[[j]])),j,0) # change all NAs introduced by rbindlist to 0
}
data12_sum <- data12_bind[, lapply(.SD, sum), by=rn]
nrow(data12_sum) # Correct, 9585 = 9546 + (9585-9546) + (9546-9546)

data12_combined <- data.frame(data12_sum)[,2:13]
rownames(data12_combined) <- data12_sum$rn
dim(data12_combined) # 9585   12


# data7, data8, data9
data7 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20160920_Pool7_counts.txt", header = T)
data8 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20161011_Pool8_counts.txt", header = T)
data9 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019/20161011_Pool9_counts.txt", header = T)


# Combine table of counts into a single data structure
colnames(data1_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data2_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data3_combined2) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data4_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data5_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data6_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data7) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data8) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data9) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data10_combined2) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data11_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")
colnames(data12_combined) <- c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3")

data <- rbind(data1_combined, data2_combined, data3_combined2, data4_combined, data5_combined, data6_combined, data7, data8, data9, data10_combined2, data11_combined, data12_combined)
dim(data) # 110669     12
# only gain

# 100*(110669/113002) = 97.9% shRNAs detected after the second resequencing, after the first we had pretty much the same amount (110666) but before, without adding any resequencing we had 110576 only.


# We feel confident about this dataset (in fact this will be the table of counts that we will be sharing in papers ...) therefore I create a directory and a table to store this dataset for the future. That will save having to load all the data as shown above again in the future. We would only have to load the table.
system("mkdir /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170125")
write.table(data, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170125/20170125_counts.txt", quote = FALSE, sep = "\t", row.names=TRUE)


# Just to check I will be loading the dataset in a separate objects and compare that the two objects data and data2 are the same
data_2 <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170125/20170125_counts.txt", header = T)
sum(data == data_2) == 110669 * 12 # TRUE
#https://cran.r-project.org/web/packages/compare/index.html
#install.packages("compare")
library(compare)
comparison <- compare(data, data_2, allowAll=TRUE)
comparison # TRUE


# Change to a DGE object
groups <- c("t0_no", "t0_no", "t0_no", "t15_DMSO", "t15_DMSO", "t15_DMSO", "t15_PDS", "t15_PDS", "t15_PDS", "t15_PhenDC3", "t15_PhenDC3", "t15_PhenDC3")
genes_dataframe <- data.frame(shRNA = rownames(data), gene = as.character(sapply(rownames(data), function(x) unlist(strsplit(x, "__"))[1])))
z <- DGEList(counts = data, group = groups, genes = genes_dataframe)
z$genes$pool <- sapply(rownames(data), function(x) gsub("Pool", "", tail(unlist(strsplit(x, "_")), n=1)))


# How many shRNAs are detected per pool?
table(z$genes$pool)
#1   10   11   12    2    3    4    5    6    7    8    9
#9480 9497 7981 9585 9468 9525 8429 8862 9550 9553 9581 9158


# What are the genes with more different shRNAs detected?
head(sort(table(z$genes$gene), decreasing = TRUE), n = 20)
#CARM1  PRKCD    ATR  FKBP5  GTF2B  KAT6B SFMBT2  SP100  KDM1B  KDM2B  PARP1  SETD3 TRIM28  DNMT1   EPC2 FKBP1A HDAC11  HDAC5  KDM3B NAP1L3
#   21     20     19     19     19     19     19     19     18     18     18     18     18     17     17     17     17     17     17     17


# Vector containing all pools
pools <- unique(z$genes$pool)


# Define empty data.tables to contain shRNAs for all pools
top_DMSOvst0 <- data.table()
top_PDSvst0 <- data.table()
top_PhenDC3vst0 <- data.table()


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

  # Get ranked shRNAS, combine shRNAs into a table for all pools and write FDR and LogFC(t15/t0) tables for DMSOvst0, PDSvst0 and PhenDC3vst0
  # DMSOvst0
  top_pool_DMSOvst0 <- topTags(lrt_pool_DMSOvst0, n=Inf)
  top_DMSOvst0 <- rbind(top_DMSOvst0, data.table(top_pool_DMSOvst0$table)[order(FDR),.(gene, shRNA, pool, logFC, FDR)])
  # PDSvst0
  top_pool_PDSvst0 <- topTags(lrt_pool_PDSvst0, n=Inf)
  top_PDSvst0 <- rbind(top_PDSvst0, data.table(top_pool_PDSvst0$table)[order(FDR),.(gene, shRNA, pool, logFC, FDR)])
  # PhenDC3vst0
  top_pool_PhenDC3vst0 <- topTags(lrt_pool_PhenDC3vst0, n=Inf)
  top_PhenDC3vst0 <- rbind(top_PhenDC3vst0, data.table(top_pool_PhenDC3vst0$table)[order(FDR),.(gene, shRNA, pool, logFC, FDR)])
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
#[1] "shRNAs: 9525"
#[1] "shRNAs (after filtering): 9419"
#[1] "Common dispersion: 0.402853816278339"
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
#[1] "shRNAs: 7981"
#[1] "shRNAs (after filtering): 7726"
#[1] "Common dispersion: 0.381118495962652"
#[1] "Pool12"
#[1] "shRNAs: 9585"
#[1] "shRNAs (after filtering): 9376"
#[1] "Common dispersion: 0.402961184359746"


# Create directory structure for figures and tables below
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/genes_all")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/shRNAs_all")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/genes_lost")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/venn_diagrams/pass_50pct_or_3_threshold")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/genes_pass_50pct_or_3_threshold")



# Tables of genes by DMSOvst0, PDSvst0 and PhenDC3vst0 and for resistant (logFC>0) or sensitiser (logFC<0)

# All genes
f <- readDNAStringSet("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3.fa")
all_shRNAs <- names(f)
all_genes <- data.table(table(as.character(sapply(all_shRNAs, function(x) unlist(strsplit(x, "__"))[1])))) # all genes, 18937


# DMSOvst0
detected_genes_DMSOvst0 <- top_DMSOvst0[,.(n_detected_shRNA = uniqueN(shRNA)), by = gene] # detected genes DMSOvst0, 18908, 29 genes with no shRNAs detected, shRNAs filtered are not included in this table either

all_detected_genes_DMSOvst0 <- merge(all_genes, detected_genes_DMSOvst0, by.x="V1", by.y="gene", all=T)
setnames(all_detected_genes_DMSOvst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs"))

significant_resistant_genes_DMSOvst0 <- top_DMSOvst0[FDR <= 0.05 & logFC > 0,.(n_significant_resistant_shRNAs = uniqueN(shRNA), logFC_significant_resistant_shRNAs_average = mean(logFC), FC_significant_resistant_shRNAs_average = 2^(mean(logFC)), logFC_significant_resistant_shRNAs_median = median(logFC), FC_significant_resistant_shRNAs_median = 2^(median(logFC))), by = gene] # significant and resistant genes DMSOvst0, 4809
significant_sensitiser_genes_DMSOvst0 <- top_DMSOvst0[FDR <= 0.05 & logFC < 0,.(n_significant_sensitiser_shRNAs = uniqueN(shRNA), logFC_significant_sensitiser_shRNAs_average = mean(logFC), FC_significant_sensitiser_shRNAs_average = 2^(mean(logFC)), logFC_significant_sensitiser_shRNAs_median = median(logFC), FC_significant_sensitiser_shRNAs_median = 2^(median(logFC))), by = gene] # significant and sensitiser genes DMSOvst0, 4346
length(intersect(significant_resistant_genes_DMSOvst0$gene, significant_sensitiser_genes_DMSOvst0$gene)) # 843

all_detected_significant_resistant_genes_DMSOvst0 <- merge(all_detected_genes_DMSOvst0, significant_resistant_genes_DMSOvst0, by.x="gene", by.y="gene", all=T)
all_detected_significant_resistant_sensitiser_genes_DMSOvst0 <- merge(all_detected_significant_resistant_genes_DMSOvst0, significant_sensitiser_genes_DMSOvst0, by.x="gene", by.y="gene", all=T)

for (j in c("n_detected_shRNAs", "n_significant_resistant_shRNAs", "n_significant_sensitiser_shRNAs")){
  set(all_detected_significant_resistant_sensitiser_genes_DMSOvst0, which(is.na(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[[j]])),j,0)
}

all_detected_significant_resistant_sensitiser_genes_DMSOvst0[, pct_significant_resistant_shRNAs := ifelse(n_detected_shRNAs > 0, 100*n_significant_resistant_shRNAs/n_detected_shRNAs, 0)]
all_detected_significant_resistant_sensitiser_genes_DMSOvst0[, pass_50_threshold_significant_resistant_shRNAs := ifelse(pct_significant_resistant_shRNAs >= 50, "yes", "no")]
all_detected_significant_resistant_sensitiser_genes_DMSOvst0[, pass_3_significant_resistant_shRNAs := ifelse(n_significant_resistant_shRNAs >= 3, "yes", "no")]
all_detected_significant_resistant_sensitiser_genes_DMSOvst0[, pct_significant_sensitiser_shRNAs := ifelse(n_detected_shRNAs > 0, 100*n_significant_sensitiser_shRNAs/n_detected_shRNAs, 0)]
all_detected_significant_resistant_sensitiser_genes_DMSOvst0[, pass_50_threshold_significant_sensitiser_shRNAs := ifelse(pct_significant_sensitiser_shRNAs >= 50, "yes", "no")]
all_detected_significant_resistant_sensitiser_genes_DMSOvst0[, pass_3_significant_sensitiser_shRNAs := ifelse(n_significant_sensitiser_shRNAs >= 3, "yes", "no")]

setcolorder(all_detected_significant_resistant_sensitiser_genes_DMSOvst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs", "n_significant_resistant_shRNAs", "logFC_significant_resistant_shRNAs_average", "FC_significant_resistant_shRNAs_average", "logFC_significant_resistant_shRNAs_median", "FC_significant_resistant_shRNAs_median", "pct_significant_resistant_shRNAs", "pass_50_threshold_significant_resistant_shRNAs", "pass_3_significant_resistant_shRNAs", "n_significant_sensitiser_shRNAs", "logFC_significant_sensitiser_shRNAs_average", "FC_significant_sensitiser_shRNAs_average", "logFC_significant_sensitiser_shRNAs_median", "FC_significant_sensitiser_shRNAs_median", "pct_significant_sensitiser_shRNAs", "pass_50_threshold_significant_sensitiser_shRNAs", "pass_3_significant_sensitiser_shRNAs"))

all_detected_significant_resistant_sensitiser_genes_DMSOvst0["KDM3B"]
#gene n_total_shRNAs n_detected_shRNAs n_significant_resistant_shRNAs logFC_significant_resistant_shRNAs_average FC_significant_resistant_shRNAs_average logFC_significant_resistant_shRNAs_median FC_significant_resistant_shRNAs_median
#1: KDM3B             18                17                              4                                   1.670551                                3.183362                                  1.448301                               2.728865
#pct_significant_resistant_shRNAs pass_50_threshold_significant_resistant_shRNAs pass_3_significant_resistant_shRNAs n_significant_sensitiser_shRNAs logFC_significant_sensitiser_shRNAs_average FC_significant_sensitiser_shRNAs_average
#1:                         23.52941                                             no                                 yes                               2                                   -2.234277                                0.2125278
#logFC_significant_sensitiser_shRNAs_median FC_significant_sensitiser_shRNAs_median pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs pass_3_significant_sensitiser_shRNAs
#1:                                  -2.234277                               0.2125278                          11.76471                                              no                                   no

nrow(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"]) # 284 genes are resistant in DMSOvst0, have more than 50% or at least 3 shRNAS detected with FDR <= 0.05 & logFC > 0
nrow(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"]) # 408 genes are sensitiser in DMSOvst0, have more than 50% or at least 3 shRNAS detected with FDR <= 0.05 & logFC < 0
nrow(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[(pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes") & (pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes")]) # 2 genes are both resistant and sensitiser in DMSOvst0, OR8B8 and ZNF434

write.table(all_detected_significant_resistant_sensitiser_genes_DMSOvst0, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_DMSOvst0_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_DMSOvst0_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/genes_all")


# PDSvst0
detected_genes_PDSvst0 <- top_PDSvst0[,.(n_detected_shRNA = uniqueN(shRNA)), by = gene] # detected genes PDSvst0, 18908, 29 genes with no shRNAs detected, shRNAs filtered are not included in this table either

all_detected_genes_PDSvst0 <- merge(all_genes, detected_genes_PDSvst0, by.x="V1", by.y="gene", all=T)
setnames(all_detected_genes_PDSvst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs"))

significant_resistant_genes_PDSvst0 <- top_PDSvst0[FDR <= 0.05 & logFC > 0,.(n_significant_resistant_shRNAs = uniqueN(shRNA), logFC_significant_resistant_shRNAs_average = mean(logFC), FC_significant_resistant_shRNAs_average = 2^(mean(logFC)), logFC_significant_resistant_shRNAs_median = median(logFC), FC_significant_resistant_shRNAs_median = 2^(median(logFC))), by = gene] # significant and resistant genes PDSvst0, 5743
significant_sensitiser_genes_PDSvst0 <- top_PDSvst0[FDR <= 0.05 & logFC < 0,.(n_significant_sensitiser_shRNAs = uniqueN(shRNA), logFC_significant_sensitiser_shRNAs_average = mean(logFC), FC_significant_sensitiser_shRNAs_average = 2^(mean(logFC)), logFC_significant_sensitiser_shRNAs_median = median(logFC), FC_significant_sensitiser_shRNAs_median = 2^(median(logFC))), by = gene] # significant and sensitiser genes PDSvst0, 5198
length(intersect(significant_resistant_genes_PDSvst0$gene, significant_sensitiser_genes_PDSvst0$gene)) # 1237

all_detected_significant_resistant_genes_PDSvst0 <- merge(all_detected_genes_PDSvst0, significant_resistant_genes_PDSvst0, by.x="gene", by.y="gene", all=T)
all_detected_significant_resistant_sensitiser_genes_PDSvst0 <- merge(all_detected_significant_resistant_genes_PDSvst0, significant_sensitiser_genes_PDSvst0, by.x="gene", by.y="gene", all=T)

for (j in c("n_detected_shRNAs", "n_significant_resistant_shRNAs", "n_significant_sensitiser_shRNAs")){
  set(all_detected_significant_resistant_sensitiser_genes_PDSvst0, which(is.na(all_detected_significant_resistant_sensitiser_genes_PDSvst0[[j]])),j,0)
}

all_detected_significant_resistant_sensitiser_genes_PDSvst0[, pct_significant_resistant_shRNAs := ifelse(n_detected_shRNAs > 0, 100*n_significant_resistant_shRNAs/n_detected_shRNAs, 0)]
all_detected_significant_resistant_sensitiser_genes_PDSvst0[, pass_50_threshold_significant_resistant_shRNAs := ifelse(pct_significant_resistant_shRNAs >= 50, "yes", "no")]
all_detected_significant_resistant_sensitiser_genes_PDSvst0[, pass_3_significant_resistant_shRNAs := ifelse(n_significant_resistant_shRNAs >= 3, "yes", "no")]
all_detected_significant_resistant_sensitiser_genes_PDSvst0[, pct_significant_sensitiser_shRNAs := ifelse(n_detected_shRNAs > 0, 100*n_significant_sensitiser_shRNAs/n_detected_shRNAs, 0)]
all_detected_significant_resistant_sensitiser_genes_PDSvst0[, pass_50_threshold_significant_sensitiser_shRNAs := ifelse(pct_significant_sensitiser_shRNAs >= 50, "yes", "no")]
all_detected_significant_resistant_sensitiser_genes_PDSvst0[, pass_3_significant_sensitiser_shRNAs := ifelse(n_significant_sensitiser_shRNAs >= 3, "yes", "no")]

setcolorder(all_detected_significant_resistant_sensitiser_genes_PDSvst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs", "n_significant_resistant_shRNAs", "logFC_significant_resistant_shRNAs_average", "FC_significant_resistant_shRNAs_average", "logFC_significant_resistant_shRNAs_median", "FC_significant_resistant_shRNAs_median", "pct_significant_resistant_shRNAs", "pass_50_threshold_significant_resistant_shRNAs", "pass_3_significant_resistant_shRNAs", "n_significant_sensitiser_shRNAs", "logFC_significant_sensitiser_shRNAs_average", "FC_significant_sensitiser_shRNAs_average", "logFC_significant_sensitiser_shRNAs_median", "FC_significant_sensitiser_shRNAs_median", "pct_significant_sensitiser_shRNAs", "pass_50_threshold_significant_sensitiser_shRNAs", "pass_3_significant_sensitiser_shRNAs"))

all_detected_significant_resistant_sensitiser_genes_PDSvst0["KDM3B"]
#gene n_total_shRNAs n_detected_shRNAs n_significant_resistant_shRNAs logFC_significant_resistant_shRNAs_average FC_significant_resistant_shRNAs_average logFC_significant_resistant_shRNAs_median FC_significant_resistant_shRNAs_median
#1: KDM3B             18                17                              6                                   1.306184                                2.472866                                  1.113402                               2.163552
#pct_significant_resistant_shRNAs pass_50_threshold_significant_resistant_shRNAs pass_3_significant_resistant_shRNAs n_significant_sensitiser_shRNAs logFC_significant_sensitiser_shRNAs_average FC_significant_sensitiser_shRNAs_average
#1:                         35.29412                                             no                                 yes                               1                                     -3.2964                                0.1017852
#logFC_significant_sensitiser_shRNAs_median FC_significant_sensitiser_shRNAs_median pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs pass_3_significant_sensitiser_shRNAs
#1:                                    -3.2964                               0.1017852                          5.882353                                              no                                   no

all_detected_significant_resistant_sensitiser_genes_PDSvst0["DHX36"]
#gene n_total_shRNAs n_detected_shRNAs n_significant_resistant_shRNAs logFC_significant_resistant_shRNAs_average FC_significant_resistant_shRNAs_average logFC_significant_resistant_shRNAs_median FC_significant_resistant_shRNAs_median
#1: DHX36              8                 5                              0                                         NA                                      NA                                        NA                                     NA
#pct_significant_resistant_shRNAs pass_50_threshold_significant_resistant_shRNAs pass_3_significant_resistant_shRNAs n_significant_sensitiser_shRNAs logFC_significant_sensitiser_shRNAs_average FC_significant_sensitiser_shRNAs_average
#1:                                0                                             no                                  no                               1                                   -1.935227                                0.2614801
#logFC_significant_sensitiser_shRNAs_median FC_significant_sensitiser_shRNAs_median pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs pass_3_significant_sensitiser_shRNAs
#1:                                  -1.935227                               0.2614801                                20                                              no                                   no

nrow(all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"]) # 447 genes are resistant in PDSvst0, have more than 50% or at least 3 shRNAS detected with FDR <= 0.05 & logFC > 0
nrow(all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"]) # 555 genes are sensitiser in PDSvst0, have more than 50% or at least 3 shRNAS detected with FDR <= 0.05 & logFC < 0
nrow(all_detected_significant_resistant_sensitiser_genes_PDSvst0[(pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes") & (pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes")]) # 4 genes are both resistant and sensitiser in PDSvst0, COL8A2, DCUN1D4, KIAA1370, OR6T1

write.table(all_detected_significant_resistant_sensitiser_genes_PDSvst0, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PDSvst0_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PDSvst0_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/genes_all")


# PhenDC3vst0
detected_genes_PhenDC3vst0 <- top_PhenDC3vst0[,.(n_detected_shRNA = uniqueN(shRNA)), by = gene] # detected genes PhenDC3vst0, 18908, 29 genes with no shRNAs detected, shRNAs filtered are not included in this table either

all_detected_genes_PhenDC3vst0 <- merge(all_genes, detected_genes_PhenDC3vst0, by.x="V1", by.y="gene", all=T)
setnames(all_detected_genes_PhenDC3vst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs"))

significant_resistant_genes_PhenDC3vst0 <- top_PhenDC3vst0[FDR <= 0.05 & logFC > 0,.(n_significant_resistant_shRNAs = uniqueN(shRNA), logFC_significant_resistant_shRNAs_average = mean(logFC), FC_significant_resistant_shRNAs_average = 2^(mean(logFC)), logFC_significant_resistant_shRNAs_median = median(logFC), FC_significant_resistant_shRNAs_median = 2^(median(logFC))), by = gene] # significant and resistant genes PhenDC3vst0, 6602
significant_sensitiser_genes_PhenDC3vst0 <- top_PhenDC3vst0[FDR <= 0.05 & logFC < 0,.(n_significant_sensitiser_shRNAs = uniqueN(shRNA), logFC_significant_sensitiser_shRNAs_average = mean(logFC), FC_significant_sensitiser_shRNAs_average = 2^(mean(logFC)), logFC_significant_sensitiser_shRNAs_median = median(logFC), FC_significant_sensitiser_shRNAs_median = 2^(median(logFC))), by = gene] # significant and sensitiser genes PhenDC3vst0, 6537
length(intersect(significant_resistant_genes_PhenDC3vst0$gene, significant_sensitiser_genes_PhenDC3vst0$gene)) # 1985, e.g. KDM3B

all_detected_significant_resistant_genes_PhenDC3vst0 <- merge(all_detected_genes_PhenDC3vst0, significant_resistant_genes_PhenDC3vst0, by.x="gene", by.y="gene", all=T)
all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0 <- merge(all_detected_significant_resistant_genes_PhenDC3vst0, significant_sensitiser_genes_PhenDC3vst0, by.x="gene", by.y="gene", all=T)

for (j in c("n_detected_shRNAs", "n_significant_resistant_shRNAs", "n_significant_sensitiser_shRNAs")){
  set(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0, which(is.na(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[[j]])),j,0)
}

all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[, pct_significant_resistant_shRNAs := ifelse(n_detected_shRNAs > 0, 100*n_significant_resistant_shRNAs/n_detected_shRNAs, 0)]
all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[, pass_50_threshold_significant_resistant_shRNAs := ifelse(pct_significant_resistant_shRNAs >= 50, "yes", "no")]
all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[, pass_3_significant_resistant_shRNAs := ifelse(n_significant_resistant_shRNAs >= 3, "yes", "no")]
all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[, pct_significant_sensitiser_shRNAs := ifelse(n_detected_shRNAs > 0, 100*n_significant_sensitiser_shRNAs/n_detected_shRNAs, 0)]
all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[, pass_50_threshold_significant_sensitiser_shRNAs := ifelse(pct_significant_sensitiser_shRNAs >= 50, "yes", "no")]
all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[, pass_3_significant_sensitiser_shRNAs := ifelse(n_significant_sensitiser_shRNAs >= 3, "yes", "no")]

setcolorder(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs", "n_significant_resistant_shRNAs", "logFC_significant_resistant_shRNAs_average", "FC_significant_resistant_shRNAs_average",  "logFC_significant_resistant_shRNAs_median", "FC_significant_resistant_shRNAs_median", "pct_significant_resistant_shRNAs", "pass_50_threshold_significant_resistant_shRNAs", "pass_3_significant_resistant_shRNAs", "n_significant_sensitiser_shRNAs", "logFC_significant_sensitiser_shRNAs_average", "FC_significant_sensitiser_shRNAs_average", "logFC_significant_sensitiser_shRNAs_median", "FC_significant_sensitiser_shRNAs_median", "pct_significant_sensitiser_shRNAs", "pass_50_threshold_significant_sensitiser_shRNAs", "pass_3_significant_sensitiser_shRNAs"))

all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0["KDM3B"]
#gene n_total_shRNAs n_detected_shRNAs n_significant_resistant_shRNAs logFC_significant_resistant_shRNAs_average FC_significant_resistant_shRNAs_average logFC_significant_resistant_shRNAs_median FC_significant_resistant_shRNAs_median
#1: KDM3B             18                17                              2                                   1.987255                                 3.96482                                  1.987255                                3.96482
#pct_significant_resistant_shRNAs pass_50_threshold_significant_resistant_shRNAs pass_3_significant_resistant_shRNAs n_significant_sensitiser_shRNAs logFC_significant_sensitiser_shRNAs_average FC_significant_sensitiser_shRNAs_average
#1:                         11.76471                                             no                                  no                               1                                   -4.795345                               0.03601284
#logFC_significant_sensitiser_shRNAs_median FC_significant_sensitiser_shRNAs_median pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs pass_3_significant_sensitiser_shRNAs
#1:                                  -4.795345                              0.03601284                          5.882353                                              no                                   no

all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0["DHX36"]
#gene n_total_shRNAs n_detected_shRNAs n_significant_resistant_shRNAs logFC_significant_resistant_shRNAs_average FC_significant_resistant_shRNAs_average logFC_significant_resistant_shRNAs_median FC_significant_resistant_shRNAs_median
#1: DHX36              8                 5                              0                                         NA                                      NA                                        NA                                     NA
#pct_significant_resistant_shRNAs pass_50_threshold_significant_resistant_shRNAs pass_3_significant_resistant_shRNAs n_significant_sensitiser_shRNAs logFC_significant_sensitiser_shRNAs_average FC_significant_sensitiser_shRNAs_average
#1:                                0                                             no                                  no                               4                                   -4.950305                               0.03234518
#logFC_significant_sensitiser_shRNAs_median FC_significant_sensitiser_shRNAs_median pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs pass_3_significant_sensitiser_shRNAs
#1:                                  -4.618994                               0.0406953                                80                                             yes                                  yes

nrow(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"]) # 537 genes are resistant in PhenDC3vst0, have more than 50% or at least 3 shRNAS detected with FDR <= 0.05 & logFC > 0
nrow(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"]) # 777 genes are sensitiser in PhenDC3vst0, have more than 50% or at least 3 shRNAS detected with FDR <= 0.05 & logFC < 0
nrow(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[(pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes") & (pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes")]) # 11 genes are both resistant and sensitiser in PhenDC3vst0 at the same time, C10orf46, HMGB1, LIMCH1, MIER3, OR5I1, OSBPL11, PI3, TSPAN7, TTC14, UTRN, UTS2D

write.table(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PhenDC3vst0_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PhenDC3vst0_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/genes_all")



# Tables of shRNAs by DMSOvst0, PDSvst0 and PhenDC3vst0

scores <- data.table(read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3_scores.txt"))
sum(is.na(as.character(scores$V2))) # 31, these were "Need score in the original table"
scores[V2 == "NULL"] # 2984
scores[V2 == "NULL", V2 := "10"] # see 20161206 14:10 email from Darcie
setnames(scores, c("cshlid","score") )

# DMSOvst0
top_DMSOvst0_scores <- merge(top_DMSOvst0, scores, by.x="shRNA", by.y="cshlid", all.x=T)[order(gene, FDR), .(gene, shRNA, score, pool, logFC, FDR)]
write.table(top_DMSOvst0_scores, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_DMSOvst0_shRNAs.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_DMSOvst0_shRNAs.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/shRNAs_all")

# PDSvst0
top_PDSvst0_scores <- merge(top_PDSvst0, scores, by.x="shRNA", by.y="cshlid", all.x=T)[order(gene, FDR), .(gene, shRNA, score, pool, logFC, FDR)]
write.table(top_PDSvst0_scores, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PDSvst0_shRNAs.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PDSvst0_shRNAs.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/shRNAs_all")

# PhenDC3vst0
top_PhenDC3vst0_scores <- merge(top_PhenDC3vst0, scores, by.x="shRNA", by.y="cshlid", all.x=T)[order(gene, FDR), .(gene, shRNA, score, pool, logFC, FDR)]
write.table(top_PhenDC3vst0_scores, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PhenDC3vst0_shRNAs.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PhenDC3vst0_shRNAs.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/shRNAs_all")



# Loss of genes analysis. This is better way of doing it than in 20161124_lost_genes.md

# Genes found in the raw table of counts (before filtering)
all_genes_counts <- data.table(table(as.character(z$genes$gene))) # 18937 - 18924 = 13 genes not found from Nic's table

# Genes found in the raw table of counts (after filtering)
all_genes_filter <- data.table(table(as.character(top_DMSOvst0$gene))) # 18924 - 18908 = 16 genes removed by t0 filtering

# Genes lost before filtering
genes_lost_before_filtering <- setdiff(all_genes$V1, all_genes_counts$V1)

# Genes lost after filtering
genes_lost_after_filtering <- setdiff(all_genes_counts$V1, all_genes_filter$V1)

# Combine and write output table
l <- list(data.table(gene = genes_lost_before_filtering), data.table(gene = genes_lost_after_filtering))
setattr(l, 'names', c("before_filtering", "after_filtering"))
l_combined <- rbindlist(l, idcol="lost")

write.table(setcolorder(l_combined, c("gene", "lost")), file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_genes_lost_t0_bypool.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_genes_lost_t0_bypool.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/genes_lost")



# Venn diagram resistant

# Apply 50% or at least 3 shRNAS cutoff and set vectors containing significant genes for all pools for DMSOvst0, PDSvst0 and PhenDC3vst0
# DMSOvst0
topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_DMSOvst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"]$gene))
# PDSvst0
topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PDSvst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"]$gene))
# PhenDC3vst0
topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PhenDC3vst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"]$gene))

venn.plot <- draw.triple.venn(
  area1 = length(topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PDSvst0),
  area2 = length(topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PhenDC3vst0),
  area3 = length(topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_DMSOvst0),
  n12 = length(intersect(topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PDSvst0, topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PhenDC3vst0)),
  n23 = length(intersect(topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PhenDC3vst0, topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_DMSOvst0)),
  n13 = length(intersect(topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PDSvst0, topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_DMSOvst0)),
  n123 = length(intersect(intersect(topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PDSvst0, topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PhenDC3vst0), topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_DMSOvst0)),
  category = c(sprintf("PDS (%i)", length(topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PDSvst0)), sprintf("PhenDC3 (%i)", length(topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PhenDC3vst0)), sprintf("DMSO (%i)", length(topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_DMSOvst0))),
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

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170130_t15vst0_genes_allpools_venn_FDR5minus2_pass_50pct_or_3_threshold_significant_resistant_shRNAs.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170130_t15vst0_genes_allpools_venn_FDR5minus2_pass_50pct_or_3_threshold_significant_resistant_shRNAs.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/venn_diagrams/pass_50pct_or_3_threshold")



# Venn diagram sensitiser

# Apply 50% or at least 3 shRNAS cutoff and set vectors containing significant genes for all pools for DMSOvst0, PDSvst0 and PhenDC3vst0
# DMSOvst0
topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_DMSOvst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"]$gene))
# PDSvst0
topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSvst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"]$gene))
# PhenDC3vst0
topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PhenDC3vst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"]$gene))

venn.plot <- draw.triple.venn(
  area1 = length(topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSvst0),
  area2 = length(topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PhenDC3vst0),
  area3 = length(topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_DMSOvst0),
  n12 = length(intersect(topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSvst0, topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PhenDC3vst0)),
  n23 = length(intersect(topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PhenDC3vst0, topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_DMSOvst0)),
  n13 = length(intersect(topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSvst0, topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_DMSOvst0)),
  n123 = length(intersect(intersect(topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSvst0, topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PhenDC3vst0), topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_DMSOvst0)),
  category = c(sprintf("PDS (%i)", length(topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSvst0)), sprintf("PhenDC3 (%i)", length(topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PhenDC3vst0)), sprintf("DMSO (%i)", length(topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_DMSOvst0))),
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

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170130_t15vst0_genes_allpools_venn_FDR5minus2_pass_50pct_or_3_threshold_significant_sensitiser_shRNAs.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170130_t15vst0_genes_allpools_venn_FDR5minus2_pass_50pct_or_3_threshold_significant_sensitiser_shRNAs.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/venn_diagrams/pass_50pct_or_3_threshold")



# Venn diagram sensitiser (same as above and logFC < -1)

# Apply 50% or at least 3 shRNAS cutoff and set vectors containing significant genes for all pools for DMSOvst0, PDSvst0 and PhenDC3vst0
# DMSOvst0
topids_genes_pass_50_3_minus1_threshold_significant_sensitiser_shRNAs_DMSOvst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[(pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes") & logFC_significant_sensitiser_shRNAs_median <= -1]$gene))
# PDSvst0
topids_genes_pass_50_3_minus1_threshold_significant_sensitiser_shRNAs_PDSvst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_PDSvst0[(pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes") & logFC_significant_sensitiser_shRNAs_median <= -1]$gene))
# PhenDC3vst0
topids_genes_pass_50_3_minus1_threshold_significant_sensitiser_shRNAs_PhenDC3vst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[(pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes") & logFC_significant_sensitiser_shRNAs_median <= -1]$gene))

venn.plot <- draw.triple.venn(
  area1 = length(topids_genes_pass_50_3_minus1_threshold_significant_sensitiser_shRNAs_PDSvst0),
  area2 = length(topids_genes_pass_50_3_minus1_threshold_significant_sensitiser_shRNAs_PhenDC3vst0),
  area3 = length(topids_genes_pass_50_3_minus1_threshold_significant_sensitiser_shRNAs_DMSOvst0),
  n12 = length(intersect(topids_genes_pass_50_3_minus1_threshold_significant_sensitiser_shRNAs_PDSvst0, topids_genes_pass_50_3_minus1_threshold_significant_sensitiser_shRNAs_PhenDC3vst0)),
  n23 = length(intersect(topids_genes_pass_50_3_minus1_threshold_significant_sensitiser_shRNAs_PhenDC3vst0, topids_genes_pass_50_3_minus1_threshold_significant_sensitiser_shRNAs_DMSOvst0)),
  n13 = length(intersect(topids_genes_pass_50_3_minus1_threshold_significant_sensitiser_shRNAs_PDSvst0, topids_genes_pass_50_3_minus1_threshold_significant_sensitiser_shRNAs_DMSOvst0)),
  n123 = length(intersect(intersect(topids_genes_pass_50_3_minus1_threshold_significant_sensitiser_shRNAs_PDSvst0, topids_genes_pass_50_3_minus1_threshold_significant_sensitiser_shRNAs_PhenDC3vst0), topids_genes_pass_50_3_minus1_threshold_significant_sensitiser_shRNAs_DMSOvst0)),
  category = c(sprintf("PDS (%i)", length(topids_genes_pass_50_3_minus1_threshold_significant_sensitiser_shRNAs_PDSvst0)), sprintf("PhenDC3 (%i)", length(topids_genes_pass_50_3_minus1_threshold_significant_sensitiser_shRNAs_PhenDC3vst0)), sprintf("DMSO (%i)", length(topids_genes_pass_50_3_minus1_threshold_significant_sensitiser_shRNAs_DMSOvst0))),
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

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170306_t15vst0_genes_allpools_venn_FDR5minus2_pass_50pct_or_3_threshold_significant_sensitiser_shRNAs_logFCminus1.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170306_t15vst0_genes_allpools_venn_FDR5minus2_pass_50pct_or_3_threshold_significant_sensitiser_shRNAs_logFCminus1.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/venn_diagrams/pass_50pct_or_3_threshold")



# Venn diagram sensitiser (shRNAs)

# Apply 50% or at least 3 shRNAS cutoff and set vectors containing significant genes for all pools for DMSOvst0, PDSvst0 and PhenDC3vst0
# DMSOvst0
topids_shRNAs_DMSOvst0 <- top_DMSOvst0[FDR <= 0.05 & logFC < 0]$shRNA
# PDSvst0
topids_shRNAs_PDSvst0 <- top_PDSvst0[FDR <= 0.05 & logFC < 0]$shRNA
# PhenDC3vst0
topids_shRNAs_PhenDC3vst0 <- top_PhenDC3vst0[FDR <= 0.05 & logFC < 0]$shRNA

venn.plot <- draw.triple.venn(
  area1 = length(topids_shRNAs_PDSvst0),
  area2 = length(topids_shRNAs_PhenDC3vst0),
  area3 = length(topids_shRNAs_DMSOvst0),
  n12 = length(intersect(topids_shRNAs_PDSvst0, topids_shRNAs_PhenDC3vst0)),
  n23 = length(intersect(topids_shRNAs_PhenDC3vst0, topids_shRNAs_DMSOvst0)),
  n13 = length(intersect(topids_shRNAs_PDSvst0, topids_shRNAs_DMSOvst0)),
  n123 = length(intersect(intersect(topids_shRNAs_PDSvst0, topids_shRNAs_PhenDC3vst0), topids_shRNAs_DMSOvst0)),
  category = c(sprintf("PDS (%i)", length(topids_shRNAs_PDSvst0)), sprintf("PhenDC3 (%i)", length(topids_shRNAs_PhenDC3vst0)), sprintf("DMSO (%i)", length(topids_shRNAs_DMSOvst0))),
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

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170306_t15vst0_shRNAs_allpools_venn_FDR5minus2.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170306_t15vst0_shRNAs_allpools_venn_FDR5minus2.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/venn_diagrams/pass_50pct_or_3_threshold")





# Tables of genes for PDS only, PhenDC3 only and PDS+PhenDC3 only with FDR <= 0.05 & logFC < 0 (sensitiser) or logFC > 0 (resistant) & passing the 50% or at least 3 threshold of significant shRNAS

# resistant
# PDS only
top_PDSonly_resistant <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"][!(gene %in% c(topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_DMSOvst0, topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PhenDC3vst0))][order(-logFC_significant_resistant_shRNAs_average),.(gene, n_detected_shRNAs, n_significant_resistant_shRNAs, pct_significant_resistant_shRNAs, logFC_significant_resistant_shRNAs_average, FC_significant_resistant_shRNAs_average, logFC_significant_resistant_shRNAs_median, FC_significant_resistant_shRNAs_median)]
write.table(top_PDSonly_resistant, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PDSonly_resistant_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PDSonly_resistant_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/genes_pass_50pct_or_3_threshold")

# PhenDC3 only
top_PhenDC3only_resistant <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"][!(gene %in% c(topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_DMSOvst0, topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PDSvst0))][order(-logFC_significant_resistant_shRNAs_average),.(gene, n_detected_shRNAs, n_significant_resistant_shRNAs, pct_significant_resistant_shRNAs, logFC_significant_resistant_shRNAs_average, FC_significant_resistant_shRNAs_average, logFC_significant_resistant_shRNAs_median, FC_significant_resistant_shRNAs_median)]
write.table(top_PhenDC3only_resistant, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PhenDC3only_resistant_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PhenDC3only_resistant_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/genes_pass_50pct_or_3_threshold")

# PDS+PhenDC3 only
top_PDSandPhenDC3only_PhenDC3_resistant <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"][gene %in% topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PDSvst0][!(gene %in% topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_DMSOvst0)][order(-logFC_significant_resistant_shRNAs_average),.(gene, n_detected_shRNAs, n_significant_resistant_shRNAs, pct_significant_resistant_shRNAs, logFC_significant_resistant_shRNAs_average, FC_significant_resistant_shRNAs_average, logFC_significant_resistant_shRNAs_median, FC_significant_resistant_shRNAs_median)]
top_PDSandPhenDC3only_PDS_resistant <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"][gene %in% topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PhenDC3vst0][!(gene %in% topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_DMSOvst0)][order(-logFC_significant_resistant_shRNAs_average),.(gene, n_detected_shRNAs, n_significant_resistant_shRNAs, pct_significant_resistant_shRNAs, logFC_significant_resistant_shRNAs_average, FC_significant_resistant_shRNAs_average, logFC_significant_resistant_shRNAs_median, FC_significant_resistant_shRNAs_median)]
top_PDSandPhenDC3only_resistant <- merge(top_PDSandPhenDC3only_PhenDC3_resistant, top_PDSandPhenDC3only_PDS_resistant, by = c("gene", "n_detected_shRNAs"), suffixes = c(".PhenDC3", ".PDS"))[order(-logFC_significant_resistant_shRNAs_average.PhenDC3)]
write.table(top_PDSandPhenDC3only_resistant, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PDSandPhenDC3only_resistant_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PDSandPhenDC3only_resistant_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/genes_pass_50pct_or_3_threshold")

# DMSO
top_DMSO_resistant <- all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"][order(-logFC_significant_resistant_shRNAs_average),.(gene, n_detected_shRNAs, n_significant_resistant_shRNAs, pct_significant_resistant_shRNAs, logFC_significant_resistant_shRNAs_average, FC_significant_resistant_shRNAs_average, logFC_significant_resistant_shRNAs_median, FC_significant_resistant_shRNAs_median)]
write.table(top_DMSO_resistant, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_DMSO_resistant_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_DMSO_resistant_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/genes_pass_50pct_or_3_threshold")


# sensitiser
# PDS only
top_PDSonly_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"][!(gene %in% c(topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_DMSOvst0, topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PhenDC3vst0))][order(logFC_significant_sensitiser_shRNAs_average),.(gene, n_detected_shRNAs, n_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs, logFC_significant_sensitiser_shRNAs_average, FC_significant_sensitiser_shRNAs_average, logFC_significant_sensitiser_shRNAs_median, FC_significant_sensitiser_shRNAs_median)]
write.table(top_PDSonly_sensitiser, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PDSonly_sensitiser_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PDSonly_sensitiser_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/genes_pass_50pct_or_3_threshold")

# PhenDC3 only
top_PhenDC3only_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"][!(gene %in% c(topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_DMSOvst0, topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSvst0))][order(logFC_significant_sensitiser_shRNAs_average),.(gene, n_detected_shRNAs, n_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs, logFC_significant_sensitiser_shRNAs_average, FC_significant_sensitiser_shRNAs_average, logFC_significant_sensitiser_shRNAs_median, FC_significant_sensitiser_shRNAs_median)]
write.table(top_PhenDC3only_sensitiser, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PhenDC3only_sensitiser_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PhenDC3only_sensitiser_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/genes_pass_50pct_or_3_threshold")

# PDS+PhenDC3 only
top_PDSandPhenDC3only_PhenDC3_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"][gene %in% topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSvst0][!(gene %in% topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_DMSOvst0)][order(logFC_significant_sensitiser_shRNAs_average),.(gene, n_detected_shRNAs, n_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs, logFC_significant_sensitiser_shRNAs_average, FC_significant_sensitiser_shRNAs_average, logFC_significant_sensitiser_shRNAs_median, FC_significant_sensitiser_shRNAs_median)]
top_PDSandPhenDC3only_PDS_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"][gene %in% topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PhenDC3vst0][!(gene %in% topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_DMSOvst0)][order(logFC_significant_sensitiser_shRNAs_average),.(gene, n_detected_shRNAs, n_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs, logFC_significant_sensitiser_shRNAs_average, FC_significant_sensitiser_shRNAs_average, logFC_significant_sensitiser_shRNAs_median, FC_significant_sensitiser_shRNAs_median)]
top_PDSandPhenDC3only_sensitiser <- merge(top_PDSandPhenDC3only_PhenDC3_sensitiser, top_PDSandPhenDC3only_PDS_sensitiser, by = c("gene", "n_detected_shRNAs"), suffixes = c(".PhenDC3", ".PDS"))[order(logFC_significant_sensitiser_shRNAs_average.PhenDC3)]
write.table(top_PDSandPhenDC3only_sensitiser, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PDSandPhenDC3only_sensitiser_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PDSandPhenDC3only_sensitiser_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/genes_pass_50pct_or_3_threshold")

# DMSO
top_DMSO_sensitiser <- all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"][order(logFC_significant_sensitiser_shRNAs_average),.(gene, n_detected_shRNAs, n_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs, logFC_significant_sensitiser_shRNAs_average, FC_significant_sensitiser_shRNAs_average, logFC_significant_sensitiser_shRNAs_median, FC_significant_sensitiser_shRNAs_median)]
write.table(top_DMSO_sensitiser, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_DMSO_sensitiser_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_DMSO_sensitiser_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/genes_pass_50pct_or_3_threshold")


```









# TODO

- Volcano plots FDR (y-axis) vs. log2FC (x-axis) for genes of interest

- Shiny app to better interface with Katie and Darcie, and also to have it as a resource for the lab to check protein/genes of interest

- Use NLP tool to find protein interactors with G4s / modified bases in the literature

- Move beyond classical statistics and develop an approach for shRNA screen analysis based on bayesian statistics (see recent literature).

- Use tools such as [pathview](http://bioconductor.org/packages/release/bioc/html/pathview.html) and [gage](https://bioconductor.org/packages/release/bioc/html/gage.html) to link relevant genes to functional data.

- Literature:
  - 2017: Winter2017
  - 2016: Jastrzebski2016, Winter2016, Hart2016, Tzelepis2016
  - 2015: Yu2015, Chen2015
  - 2014: Knott2014, Dai2014, Ritchie2014, Chen2014, Li2014
  - 2013: Sheridan2013, Yu2013
  - 2012: McCarthy2012
  - 2011: Zuber2011, Fellman2011, Sims2011
  - 2010: Robinson2010
  - 2006: Vert2006

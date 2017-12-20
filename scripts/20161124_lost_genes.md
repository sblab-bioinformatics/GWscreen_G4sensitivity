
Using the previous data + conservatively adding the sequencing reruns from [20161115_integration.md](20161115_integration.md)). This is the dataset where Darcie and Katie feel more confident. The idea is to generate a list of genes that were lost before the t0 filtering and after the t0 filtering. It will be helpful to know the genes for which we we can say something and the genes for which we can't really say anything.

```R
library(data.table)
library(edgeR)
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
data3_rerun <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161110/20161109_Pool3_counts.txt", header = T) # but I only need t0_no_Pool3_2 and t0_no_Pool3_3 from this

nrow(data3) # 9518
nrow(data3_rerun) # 9476
length(intersect(rownames(data3), rownames(data3_rerun))) # 9470

data3_bind <- rbindlist(list(data.table(data3, keep.rownames=TRUE), data.table(data3_rerun[,1:2], keep.rownames=TRUE)), use.names = TRUE, fill = TRUE)
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


# Filtering
keep <- rowSums(cpm(z[,1:3])>0.5) >= 3 # A library of 10M will need to have at least 5 counts in each of three t0 replicates individually in order to pass this filter
z_filter <- z[keep, , keep.lib.sizes=FALSE]


# All genes
f <- readDNAStringSet("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3.fa")
all_shRNAs <- names(f)
all_genes <- data.table(table(as.character(sapply(all_shRNAs, function(x) unlist(strsplit(x, "__"))[1])))) # all genes, 18937


# Genes found in the raw table of counts (before filtering)
all_genes_counts <- data.table(table(as.character(z$genes$gene))) # 18937 - 18924 = 13 genes not found from Nic's table


# Genes found in the raw table of counts (after filtering)
all_genes_filter <- data.table(table(as.character(z_filter$genes$gene))) # 18924 - 18784 = 140 genes removed by t0 filtering


# Genes lost before filtering
genes_lost_before_filtering <- setdiff(all_genes$V1, all_genes_counts$V1)


# Genes lost after filtering
genes_lost_after_filtering <- setdiff(all_genes_counts$V1, all_genes_filter$V1)


# Combine and write output table
l <- list(data.table(gene = genes_lost_before_filtering), data.table(gene = genes_lost_after_filtering))
setattr(l, 'names', c("before_filtering", "after_filtering"))
l_combined <- rbindlist(l, idcol="lost")

write.table(setcolorder(l_combined, c("gene", "lost")), file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161129_genes_lost_t0.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161129_genes_lost_t0.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161129")

```

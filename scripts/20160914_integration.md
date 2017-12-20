
Ideas to integrate data from different pools and perform analysis of significant hairpins on a gene basis.

* Pools 1 and 2 - Darcie - June
* Pools 5 and 10 - Katie - June
* Pools 3 and 4 - Darcie - September
* Pools 6 and 7 - Katie - September
* Pools 8 and 9 - Darcie - October
* Pools 11 and 12 - Katie - October

Each pool consists of 12 libraries:

* t0 (x3)
* t15_DMSO (x3)
* t15_PDS (x3)
* t15_PhenDC3 (x3)

(A) Combine the 6 pairs of raw tables of counts (20160627_counts.txt, 20160629_counts.txt, 20160902_counts_lims.txt, 20160920_counts.txt, 20161011_counts.txt and 20161018_counts.txt) and add zeros to hairpins not detected in the pair of pools. This will account for the pool8 off-targets (or contamination?) in pools 1 and 2. This will result in a matrix of 144 columns (12 pools x 12 libraries/pool) and 113002 rows (all hairpins). The design can be done using ~ group + pool where group is t0, t15_DMSO, t15_PDS, t15_PhenDC3 and pool is (1, 2, ..., 12) or perhaps simply a ~ group might also do.

(B) Combine the 12 tables of counts for each pool (after R processing) into a matrix of 12 columns (12 libraries/pool) and 113002 rows (all hairpins). The design will still be ~ group.

(C) Same as (B) but with FDR of 10^(-3) / 10^(-2), logFC < -1 and logFC > 1 and 50% of the hairpins coming up significant

(D) Shrink replicates to the one with the highest purity for each condition (t0, t15_DMSO, t15_PDS and t15_PhenDC3) and combine all pools into a matrix of 4 columns and 113002 rows (all hairpins).

(E) Pre-process (sum) the counts for different shRNAs targeting the same gene to create a matrix containing genes as rows and 12 libraries as columns. Then do an analysis similar to what performed for pools individually.

(F) Same as (E) but with FDR of 10^(-3) / 10^(-2), logFC < -1 and logFC > 1 and (50% of the hairpins coming up significant ?)

Strategies (B), (C), (E) and (F) seem more reasonable in the first instance.

(G) Nic's idea: obtain significant genes (FDR < 1e-3, logFC and % significant shRNAs) by analysing each pool individually - say 12 analyses - one per analysis per pool. Then combine significant genes from each pool.




# (A)

Under construction ...



# (B)

## Integrating data

Combine the count tables for each pool (the ones produced so far - see above). On the right, potential transformations that can be done to the data:

```bash
martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20160627/tables/20160627_Pool1_counts.txt
martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20160627/tables/20160627_Pool2_counts.txt
martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20160629/tables/20160629_Pool5_counts.txt
martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20160629/tables/20160629_Pool10_counts.txt # t15_PDS_Pool10_3 -> t15_PDS_Pool10_2, t15_PhenDC3_Pool10_1 -> t15_PhenDC3_Pool10_2
martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20160902/tables/20160908_Pool3_counts_lims.txt # t0_no_Pool3_2 -> t0_no_Pool3_1, t0_no_Pool3_3 -> t0_no_Pool3_1, t15_PhenDC3_Pool3_1 -> t15_PhenDC3_Pool3_2, t15_PhenDC3_Pool3_3 -> t15_PhenDC3_Pool3_2
martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20160902/tables/20160908_Pool4_counts_lims.txt # t0_no_Pool4_1 -> t0_no_Pool4_3, t0_no_Pool4_2 -> t0_no_Pool4_3
martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20160920/tables/20160920_Pool6_counts.txt # t0_no_Pool6_2 -> t0_no_Pool6_1, t0_no_Pool6_3 -> t0_no_Pool6_1 (!)
martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20160920/tables/20160920_Pool7_counts.txt # t0_no_Pool7_2 -> t0_no_Pool7_1, t0_no_Pool7_3 -> t0_no_Pool7_1
martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161011/tables/20161011_Pool8_counts.txt
martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161011/tables/20161011_Pool9_counts.txt
martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161018/tables/20161018_Pool11_counts.txt # t15_DMSO_Pool11_3 -> t15_DMSO_Pool11_2
martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161018/tables/20161018_Pool12_counts.txt

```

Copy tables of counts to a designated folder in uk-cri-lcst01:

```bash
cd /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161019

rsync martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20160627/tables/20160627_Pool1_counts.txt .
rsync martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20160627/tables/20160627_Pool2_counts.txt .
rsync martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20160629/tables/20160629_Pool5_counts.txt .
rsync martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20160629/tables/20160629_Pool10_counts.txt .
rsync martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20160902/tables/20160908_Pool3_counts_lims.txt .
rsync martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20160902/tables/20160908_Pool4_counts_lims.txt .
rsync martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20160920/tables/20160920_Pool6_counts.txt .
rsync martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20160920/tables/20160920_Pool7_counts.txt .
rsync martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161011/tables/20161011_Pool8_counts.txt .
rsync martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161011/tables/20161011_Pool9_counts.txt .
rsync martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161018/tables/20161018_Pool11_counts.txt .
rsync martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161018/tables/20161018_Pool12_counts.txt .

```


## Analysis:

```R
library(edgeR)
library(data.table)
library(ggplot2)
library(VennDiagram)
library(reshape)


# Enlarge the view width when printing tables
options(width = 300)


# Load table of counts into edgeR
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

# 100*(110576/113002) = 97.8 shRNAs detected


# Change to a DGE object
groups <- c("t0_no", "t0_no", "t0_no", "t15_DMSO", "t15_DMSO", "t15_DMSO", "t15_PDS", "t15_PDS", "t15_PDS", "t15_PhenDC3", "t15_PhenDC3", "t15_PhenDC3")
genes_dataframe <- data.frame(shRNA = rownames(data), gene = as.character(sapply(rownames(data), function(x) unlist(strsplit(x, "__"))[1])))
z <- DGEList(counts = data, group = groups, genes = genes_dataframe)


# Which genes have more shRNAs targetting them?
head(sort(table(z$genes$gene), decreasing = T))
# CARM1  PRKCD    ATR  FKBP5  GTF2B  KAT6B SFMBT2  SP100  KDM1B  KDM2B  PARP1  SETD3 TRIM28  DNMT1   EPC2 FKBP1A HDAC11  HDAC5  KDM3B NAP1L3
#    21     20     19     19     19     19     19     19     18     18     18     18     18     17     17     17     17     17     17     17


# Distribution of the number of genes by the number of shRNAs
gg <- ggplot(data.table(table(sort(table(z$genes$gene)))), aes(x = as.numeric(V1), y = N)) +
geom_bar(stat="identity") +
xlab("Number of shRNAs") +
ylab("Number of genes") +
theme_classic()

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161019_numbershRNAs_numbergenes.pdf", width = 10/2.54, height = 10/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161019_numbershRNAs_numbergenes.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/figures/")


# How many counts are of each Pool?
reads_total <- sum(z$counts) # 984949266
for (i in 1:12){
  p <- sprintf("Pool%i", i)
  reads_pool <- sum(z[grepl(paste(p, "$", sep=""), rownames(z$counts)),]$counts)
  reads_pool_pct <- round(100*reads_pool/reads_total, 2)
  hps_pool <- length(rowSums(z[grepl(paste(p, "$", sep=""), rownames(z$counts)),]$counts))
  print(sprintf("%s   %s   %s   %s", p, reads_pool, reads_pool_pct, hps_pool))
}

#[1] "Pool1   46619475   4.73   9475"
#[1] "Pool2   50627013   5.14   9467"
#[1] "Pool3   54337682   5.52   9518"
#[1] "Pool4   62283358   6.32   8384"
#[1] "Pool5   66363885   6.74   8846"
#[1] "Pool6   85638324   8.69   9544"
#[1] "Pool7   112030022   11.37   9553"
#[1] "Pool8   121181066   12.3   9581"
#[1] "Pool9   104799666   10.64   9158"
#[1] "Pool10   53866039   5.47   9486"
#[1] "Pool11   111097070   11.28   7979"
#[1] "Pool12   116105666   11.79   9585"


# How many hairpins were present in the original file provided by Nick
# in bash:
# tail -n+2 /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3.csv | cut -f9 -d "," | sort | uniq -c
#   9600 Pool #1
#   9600 Pool #10
#   8246 Pool #11
#      1 Pool #12
#   9599 Pool #12
#   9600 Pool #2
#   9600 Pool #3
#   9474 Pool #4
#   8906 Pool #5
#   9586 Pool #6
#   9590 Pool #7
#   9600 Pool #8
#   9600 Pool #9


# Pool1 100*(9475/9600) = 98.7%
# Pool2 100*(9467/9600) = 98.6%
# Pool3 100*(9518/9600) = 99.1%
# Pool4 100*(8384/9474) = 88.5%
# Pool5 100*(8846/8906) = 99.3%
# Pool6 100*(9544/9586) = 99.6%
# Pool7 100*(9553/9590) = 99.6%
# Pool8 100*(9581/9600) = 99.8%
# Pool9 100*(9158/9600) = 95.4%
# Pool10 100*(9486/9600) = 98.8%
# Pool11 100*(7979/8246) = 96.8%
# Pool12 100*(9585/9600) = 99.8%


par(mfrow=c(2,1))
barplot(colSums(z$counts), las=2, main="Counts per index", cex.names=0.5, cex.axis=0.8)
barplot(rowSums(z$counts), las=2, main="Counts per hairpin", cex.names=0.5, cex.axis=0.8, axisnames = FALSE)

sort(rowSums(z$counts)[rowSums(z$counts) > 1000000], decreasing=T)
#HEPACAM__sherwood_1u__2__2_Pool2    TAS2R42__sherwood__4__4_Pool2                TAOK1__4__3_Pool6       MKKS__sherwood__3__3_Pool4
#                         3885507                          1428436                          1104848                          1024896


# Pool content table
z_plot <- data.table(z$counts)
setcolorder(z_plot, c("t0_no_1", "t0_no_2", "t0_no_3", "t15_DMSO_1", "t15_DMSO_2", "t15_DMSO_3", "t15_PDS_1", "t15_PDS_2", "t15_PDS_3", "t15_PhenDC3_1", "t15_PhenDC3_2", "t15_PhenDC3_3"))
z_plot[, pool := as.numeric(sapply(rownames(z), function(x) gsub("Pool", "", tail(unlist(strsplit(x, "_")), n=1))))]

z_plot_pcts <- z_plot[, .(t0_no_1 = 100*sum(t0_no_1)/sum(z_plot[, t0_no_1]),
t0_no_2 = 100*sum(t0_no_2)/sum(z_plot[, t0_no_2]),
t0_no_3 = 100*sum(t0_no_3)/sum(z_plot[, t0_no_3]),
t15_DMSO_1 = 100*sum(t15_DMSO_1)/sum(z_plot[, t15_DMSO_1]),
t15_DMSO_2 = 100*sum(t15_DMSO_2)/sum(z_plot[, t15_DMSO_2]),
t15_DMSO_3 = 100*sum(t15_DMSO_3)/sum(z_plot[, t15_DMSO_3]),
t15_PDS_1 = 100*sum(t15_PDS_1)/sum(z_plot[, t15_PDS_1]),
t15_PDS_2 = 100*sum(t15_PDS_2)/sum(z_plot[, t15_PDS_2]),
t15_PDS_3 = 100*sum(t15_PDS_3)/sum(z_plot[, t15_PDS_3]),
t15_PhenDC3_1 = 100*sum(t15_PhenDC3_1)/sum(z_plot[, t15_PhenDC3_1]),
t15_PhenDC3_2 = 100*sum(t15_PhenDC3_2)/sum(z_plot[, t15_PhenDC3_2]),
t15_PhenDC3_3 = 100*sum(t15_PhenDC3_3)/sum(z_plot[, t15_PhenDC3_3])),
keyby = pool]

write.table(z_plot_pcts, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_library_content_table.txt", quote = FALSE, sep = "\t", row.names = FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_library_content_table.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/tables/")


# Pool content plot
z_plot_pcts_melt <- melt(z_plot_pcts, id = "pool", variable.name = "library", value.name = "pct")

gg <- ggplot(z_plot_pcts_melt[order(-pool)], aes(x = library, y = pct, fill = factor(pool))) +
geom_bar(stat="identity") +
xlab("") +
ylab("% Reads") +
theme_classic() +
scale_fill_discrete(name="Pool") +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.y = element_text(size = 16)) +
scale_x_discrete(labels = paste(paste(colnames(z_plot)[1:12], "(n =", as.character(colSums(z_plot)[1:12])), ")", sep=""))

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161019_library_content_plot.pdf", width = 10/2.54, height = 16/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161019_library_content_plot.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/figures/")


# Make MDS plots to visualise relationships between replicate samples

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161019_mds.pdf", width = 16/2.54, height = 16/2.54)
g <- plotMDS(z, labels = c(rep("t0", 3), rep("DMSO", 3), rep("PDS", 3), rep("PhenDC3", 3)), col = c(rep("black", 3), rep("deepskyblue3", 3), rep("red3", 3), rep("forestgreen", 3)), xlab = "Dimension 1", ylab = "Dimension 2",  gene.selection = "common")
legend("topright", legend=c("t0", "DMSO", "PDS", "PhenDC3"), col=c("black", "deepskyblue3", "red3", "forestgreen"), pch=15)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161019_mds.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/figures")


##################
# shRNA analysis #
##################

# Differential representation analysis.
# Set up design matrix for GLM.
des <- model.matrix(~ 0 + group, data = z$samples)
colnames(des) <- levels(factor(z$samples$group))
des


# Estimate dispersions
z_glm <- estimateDisp(z, des)
sqrt(z_glm$common.disp) # 0.6245679


# Plot BCVs versus abundance
plotBCV(z_glm)


# Fit negative binomial GLM
z_fit <- glmFit(z_glm, des)


# Define matrix of contrasts
my.contrasts <- makeContrasts(
t15_DMSOvst0_no = t15_DMSO - t0_no,
t15_PDSvst0_no = t15_PDS - t0_no,
t15_PhenDC3vst0_no = t15_PhenDC3 - t0_no,
t15_PDSvst15_DMSO = t15_PDS - t15_DMSO,
t15_PhenDC3vst15_DMSO = t15_PhenDC3 - t15_DMSO,
t15_PhenDC3vst15_PDS = t15_PhenDC3 - t15_PDS,
levels=des)


# Comparisons t15 vs t0 and within t15
# Carry out Likelihood ratio tests
lrt_t15_DMSOvst0_no <- glmLRT(z_fit, contrast=my.contrasts[,"t15_DMSOvst0_no"])
lrt_t15_PDSvst0_no <- glmLRT(z_fit, contrast=my.contrasts[,"t15_PDSvst0_no"])
lrt_t15_PhenDC3vst0_no <- glmLRT(z_fit, contrast=my.contrasts[,"t15_PhenDC3vst0_no"])
lrt_t15_PDSvst15_DMSO <- glmLRT(z_fit, contrast=my.contrasts[,"t15_PDSvst15_DMSO"])
lrt_t15_PhenDC3vst15_DMSO <- glmLRT(z_fit, contrast=my.contrasts[,"t15_PhenDC3vst15_DMSO"])
lrt_t15_PhenDC3vst15_PDS <- glmLRT(z_fit, contrast=my.contrasts[,"t15_PhenDC3vst15_PDS"])


# Show top ranked hairpins
topTags(lrt_t15_DMSOvst0_no)
topTags(lrt_t15_PDSvst0_no)
topTags(lrt_t15_PhenDC3vst0_no)
topTags(lrt_t15_PDSvst15_DMSO)
topTags(lrt_t15_PhenDC3vst15_DMSO)
topTags(lrt_t15_PhenDC3vst15_PDS)


# Select hairpins with FDR < 0.0001 to highlight on plot
thresh <- 0.0001 # section 4.6.5 Differential representation analysis in Chen2015 (https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)

top_t15_DMSOvst0_no <- topTags(lrt_t15_DMSOvst0_no, n=Inf)
topids_t15_DMSOvst0_no <- as.character(top_t15_DMSOvst0_no$table[top_t15_DMSOvst0_no$table$FDR < thresh,]$shRNA)
length(topids_t15_DMSOvst0_no) # 29

top_t15_PDSvst0_no <- topTags(lrt_t15_PDSvst0_no, n=Inf)
topids_t15_PDSvst0_no <- as.character(top_t15_PDSvst0_no$table[top_t15_PDSvst0_no$table$FDR < thresh,]$shRNA)
length(topids_t15_PDSvst0_no) # 158

top_t15_PhenDC3vst0_no <- topTags(lrt_t15_PhenDC3vst0_no, n=Inf)
topids_t15_PhenDC3vst0_no <- as.character(top_t15_PhenDC3vst0_no$table[top_t15_PhenDC3vst0_no$table$FDR < thresh,]$shRNA)
length(topids_t15_PhenDC3vst0_no) # 617

top_t15_PDSvst15_DMSO <- topTags(lrt_t15_PDSvst15_DMSO, n=Inf)
topids_t15_PDSvst15_DMSO <- as.character(top_t15_PDSvst15_DMSO$table[top_t15_PDSvst15_DMSO$table$FDR < thresh,]$shRNA)
length(topids_t15_PDSvst15_DMSO) # 6

top_t15_PhenDC3vst15_DMSO <- topTags(lrt_t15_PhenDC3vst15_DMSO, n=Inf)
topids_t15_PhenDC3vst15_DMSO <- as.character(top_t15_PhenDC3vst15_DMSO$table[top_t15_PhenDC3vst15_DMSO$table$FDR < thresh,]$shRNA)
length(topids_t15_PhenDC3vst15_DMSO) # 554

top_t15_PhenDC3vst15_PDS <- topTags(lrt_t15_PhenDC3vst15_PDS, n=Inf)
topids_t15_PhenDC3vst15_PDS <- as.character(top_t15_PhenDC3vst15_PDS$table[top_t15_PhenDC3vst15_PDS$table$FDR < thresh,]$shRNA)
length(topids_t15_PhenDC3vst15_PDS) # 498


# Write LogFC and FDR tables
write.table(data.table(as.data.frame(top_t15_DMSOvst0_no))[order(-logFC),.(shRNA, gene, logFC, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_t15_DMSOvst0_logFCdecreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_t15_DMSOvst0_logFCdecreasing_FDR.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/tables/")

write.table(data.table(as.data.frame(top_t15_PDSvst0_no))[order(-logFC),.(shRNA, gene, logFC, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_t15_PDSvst0_logFCdecreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_t15_PDSvst0_logFCdecreasing_FDR.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/tables/")

write.table(data.table(as.data.frame(top_t15_PhenDC3vst0_no))[order(-logFC),.(shRNA, gene, logFC, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_t15_PhenDC3vst0_logFCdecreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_t15_PhenDC3vst0_logFCdecreasing_FDR.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/tables/")

write.table(data.table(as.data.frame(top_t15_PDSvst15_DMSO))[order(-logFC),.(shRNA, gene, logFC, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_t15_PDSvst15_DMSO_logFCdecreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_t15_PDSvst15_DMSO_logFCdecreasing_FDR.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/tables/")

write.table(data.table(as.data.frame(top_t15_PhenDC3vst15_DMSO))[order(-logFC),.(shRNA, gene, logFC, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_t15_PhenDC3vst15_DMSO_logFCdecreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_t15_PhenDC3vst15_DMSO_logFCdecreasing_FDR.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/tables/")

write.table(data.table(as.data.frame(top_t15_PhenDC3vst15_PDS))[order(-logFC),.(shRNA, gene, logFC, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_t15_PhenDC3vst15_PDS_logFCdecreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_t15_PhenDC3vst15_PDS_logFCdecreasing_FDR.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/tables/")


# Plot logFC versus logCPM
plotSmear(lrt_t15_DMSOvst0_no, de.tags=topids_t15_DMSOvst0_no, main = "t15_DMSO vs t0_no")

plotSmear(lrt_t15_PDSvst0_no, de.tags=topids_t15_PDSvst0_no, main = "t15_PDS vs t0_no")
z_t15_PDSvst0_no_cpm <- data.table(cpm(z)[, grepl("t15_PDS|t0_no", colnames(cpm(z)))], keep.rownames=TRUE)
z_t15_PDSvst0_no_cpm[, log2cpmt0 := (log2(t0_no_1) + log2(t0_no_2) + log2(t0_no_3))/3]
z_t15_PDSvst0_no_cpm[, log2FC := lrt_t15_PDSvst0_no$table$logFC]
z_t15_PDSvst0_no_cpm[, hits := ifelse(rn %in% topids_t15_PDSvst0_no, ifelse(rn %in% topids_t15_DMSOvst0_no, "PDS and DMSO", "PDS only"), "not differentially abundant")]
z_t15_PDSvst0_no_cpm[, order := ifelse(hits == "not differentially abundant", 1, 2)]

gg <- ggplot(z_t15_PDSvst0_no_cpm[log2cpmt0 != -Inf][order(order)], aes(x = log2cpmt0, y = log2FC, colour = hits)) +
geom_point(size=0.25) +
xlab(expression("log"[2]*"(t0)")) +
ylab(expression("log"[2]*"(t15/t0)")) +
theme_bw() +
scale_colour_manual(name="",values = c("grey", "deepskyblue3", "red3"))
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161019_log2t0_log2FC_PDS_DMSO.pdf", width = 16/2.54, height = 10/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161019_log2t0_log2FC_PDS_DMSO.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/figures/")

plotSmear(lrt_t15_PhenDC3vst0_no, de.tags=topids_t15_PhenDC3vst0_no, main = "t15_PhenDC3 vs t0_no")
z_t15_PhenDC3vst0_no_cpm <- data.table(cpm(z)[, grepl("t15_PhenDC3|t0_no", colnames(cpm(z)))], keep.rownames=TRUE)
z_t15_PhenDC3vst0_no_cpm[, log2cpmt0 := (log2(t0_no_1) + log2(t0_no_2) + log2(t0_no_3))/3]
z_t15_PhenDC3vst0_no_cpm[, log2FC := lrt_t15_PhenDC3vst0_no$table$logFC]
z_t15_PhenDC3vst0_no_cpm[, hits := ifelse(rn %in% topids_t15_PhenDC3vst0_no, ifelse(rn %in% topids_t15_DMSOvst0_no, "PhenDC3 and DMSO", "PhenDC3 only"), "not differentially abundant")]
z_t15_PhenDC3vst0_no_cpm[, order := ifelse(hits == "not differentially abundant", 1, 2)]

gg <- ggplot(z_t15_PhenDC3vst0_no_cpm[log2cpmt0 != -Inf][order(order)], aes(x = log2cpmt0, y = log2FC, colour = hits)) +
geom_point(size=0.25) +
xlab(expression("log"[2]*"(t0)")) +
ylab(expression("log"[2]*"(t15/t0)")) +
theme_bw() +
scale_colour_manual(name="",values = c("grey", "deepskyblue3", "forestgreen"))
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161019_log2t0_log2FC_PhenDC3_DMSO.pdf", width = 16/2.54, height = 10/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161019_log2t0_log2FC_PhenDC3_DMSO.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/figures/")

plotSmear(lrt_t15_PDSvst15_DMSO, de.tags=topids_t15_PDSvst15_DMSO, main = "t15_PDS vs t15_DMSO")
plotSmear(lrt_t15_PhenDC3vst15_DMSO, de.tags=topids_t15_PhenDC3vst15_DMSO, main = "t15_PhenDC3 vs t15_DMSO")
plotSmear(lrt_t15_PhenDC3vst15_PDS, de.tags=topids_t15_PhenDC3vst15_PDS, main = "t15_PhenDC3 vs t15_PDS")


# Intersect lists of topTags
length(intersect(topids_t15_DMSOvst0_no, topids_t15_PDSvst0_no)) # 18
length(intersect(topids_t15_DMSOvst0_no, topids_t15_PhenDC3vst0_no)) # 13
length(intersect(topids_t15_PDSvst0_no, topids_t15_PhenDC3vst0_no)) # 50

write(topids_t15_DMSOvst0_no, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_topids_t15_DMSOvst0_no.txt")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_topids_t15_DMSOvst0_no.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/tables/")
write(topids_t15_PDSvst0_no, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_topids_t15_PDSvst0_no.txt")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_topids_t15_PDSvst0_no.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/tables/")
write(topids_t15_PhenDC3vst0_no, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_topids_t15_PhenDC3vst0_no.txt")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_topids_t15_PhenDC3vst0_no.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/tables/")

top_t15_PDSvst0_no_PDSonly <- top_t15_PDSvst0_no[1:length(topids_t15_PDSvst0_no),][!(topids_t15_PDSvst0_no %in% c(topids_t15_DMSOvst0_no, topids_t15_PhenDC3vst0_no)),]
write.table(top_t15_PDSvst0_no_PDSonly, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_top_t15_PDSvst0_no_PDSonly.txt", quote = FALSE, sep = "\t", row.names = FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_top_t15_PDSvst0_no_PDSonly.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/tables/")

top_t15_PhenDC3vst0_no_PhenDC3only <- top_t15_PhenDC3vst0_no[1:length(topids_t15_PhenDC3vst0_no),][!(topids_t15_PhenDC3vst0_no %in% c(topids_t15_DMSOvst0_no, topids_t15_PDSvst0_no)),]
write.table(top_t15_PhenDC3vst0_no_PhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_top_t15_PhenDC3vst0_no_PhenDC3only.txt", quote = FALSE, sep = "\t", row.names = FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_top_t15_PhenDC3vst0_no_PhenDC3only.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/tables/")

top_t15_PDSvst0_no_PDSandPhenDC3 <- top_t15_PDSvst0_no[1:length(topids_t15_PDSvst0_no),][(topids_t15_PDSvst0_no %in% topids_t15_PhenDC3vst0_no),]
top_t15_PDSvst0_no_PDSandPhenDC3only <- top_t15_PDSvst0_no_PDSandPhenDC3[!(rownames(top_t15_PDSvst0_no_PDSandPhenDC3) %in% topids_t15_DMSOvst0_no),]
write.table(top_t15_PDSvst0_no_PDSandPhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_top_t15_PDSvst0_no_PDSandPhenDC3only.txt", quote = FALSE, sep = "\t")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_top_t15_PDSvst0_no_PDSandPhenDC3only.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/tables/")

top_t15_DMSOvst0_no_DMSOonly <- top_t15_DMSOvst0_no[1:length(topids_t15_DMSOvst0_no),][!(topids_t15_DMSOvst0_no %in% c(topids_t15_PDSvst0_no, topids_t15_PhenDC3vst0_no)),]


# Venn diagram
venn.plot <- draw.triple.venn(
  area1 = length(topids_t15_PDSvst0_no),
  area2 = length(topids_t15_PhenDC3vst0_no),
  area3 = length(topids_t15_DMSOvst0_no),
  n12 = length(intersect(topids_t15_PDSvst0_no, topids_t15_PhenDC3vst0_no)),
  n23 = length(intersect(topids_t15_PhenDC3vst0_no, topids_t15_DMSOvst0_no)),
  n13 = length(intersect(topids_t15_PDSvst0_no, topids_t15_DMSOvst0_no)),
  n123 = length(intersect(intersect(topids_t15_PDSvst0_no, topids_t15_PhenDC3vst0_no), topids_t15_DMSOvst0_no)),
  category = c(sprintf("PDS (%i)", length(topids_t15_PDSvst0_no)), sprintf("PhenDC3 (%i)", length(topids_t15_PhenDC3vst0_no)), sprintf("DMSO (%i)", length(topids_t15_DMSOvst0_no))),
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


pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161019_t15vst0_venn.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161019_t15vst0_venn.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/figures/")


# Distribution log(counts)
pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161019_t15vst0_distribution_counts.pdf")
hist(log(rowSums(z$counts)), xlab = "log(counts)", main = "", xlim = c(0,12))
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161019_t15vst0_distribution_counts.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/figures/")




#################################
# Gene-set analysis with camera #
#################################

# rows with zero counts need to be removed, this is necessary for camera to work
z_noall0 <- z[!(rowSums(z$counts) == 0), ]


# Estimate dispersions
z_noall0_glm <- estimateDisp(z_noall0, des)


# Fit negative binomial GLM
z_noall0_fit <- glmFit(z_noall0_glm, des)


# Comparisons t15 vs t0 and within t15
# Carry out Likelihood ratio tests
lrt_noall0_t15_DMSOvst0_no <- glmLRT(z_noall0_fit, contrast=my.contrasts[,"t15_DMSOvst0_no"])
lrt_noall0_t15_PDSvst0_no <- glmLRT(z_noall0_fit, contrast=my.contrasts[,"t15_PDSvst0_no"])
lrt_noall0_t15_PhenDC3vst0_no <- glmLRT(z_noall0_fit, contrast=my.contrasts[,"t15_PhenDC3vst0_no"])
lrt_noall0_t15_PDSvst15_DMSO <- glmLRT(z_noall0_fit, contrast=my.contrasts[,"t15_PDSvst15_DMSO"])
lrt_noall0_t15_PhenDC3vst15_DMSO <- glmLRT(z_noall0_fit, contrast=my.contrasts[,"t15_PhenDC3vst15_DMSO"])
lrt_noall0_t15_PhenDC3vst15_PDS <- glmLRT(z_noall0_fit, contrast=my.contrasts[,"t15_PhenDC3vst15_PDS"])


# Top shRNAS
top_noall0_t15_DMSOvst0_no <- topTags(lrt_noall0_t15_DMSOvst0_no, n=Inf)
top_noall0_t15_PDSvst0_no <- topTags(lrt_noall0_t15_PDSvst0_no, n=Inf)
top_noall0_t15_PhenDC3vst0_no <- topTags(lrt_noall0_t15_PhenDC3vst0_no, n=Inf)
top_noall0_t15_PDSvst15_DMSO <- topTags(lrt_noall0_t15_PDSvst15_DMSO, n=Inf)
top_noall0_t15_PhenDC3vst15_DMSO <- topTags(lrt_noall0_t15_PhenDC3vst15_DMSO, n=Inf)
top_noall0_t15_PhenDC3vst15_PDS <- topTags(lrt_noall0_t15_PhenDC3vst15_PDS, n=Inf)


# Link shRNAs and genes
genesymbols <- as.character(z_noall0$genes$gene)
genesymbollist <- list()
unq <- unique(genesymbols)
length(unq) # 18922 genes for 110535 shRNAs

for(i in unq) {
  sel <- genesymbols==i & !is.na(genesymbols)
  genesymbollist[[i]] <- which(sel)
}


# Run camera for all genes - parametric
camera_t15_DMSOvst0_no <- camera(z_noall0_glm, index = genesymbollist, des, contrast=my.contrasts[,"t15_DMSOvst0_no"])
camera_t15_PDSvst0_no <- camera(z_noall0_glm, index = genesymbollist, des, contrast=my.contrasts[,"t15_PDSvst0_no"])
camera_t15_PhenDC3vst0_no <- camera(z_noall0_glm, index = genesymbollist, des, contrast=my.contrasts[,"t15_PhenDC3vst0_no"])
camera_t15_PDSvst15_DMSO <- camera(z_noall0_glm, index = genesymbollist, des, contrast=my.contrasts[,"t15_PDSvst15_DMSO"])
camera_t15_PhenDC3vst15_DMSO <- camera(z_noall0_glm, index = genesymbollist, des, contrast=my.contrasts[,"t15_PhenDC3vst15_DMSO"])
camera_t15_PhenDC3vst15_PDS <- camera(z_noall0_glm, index = genesymbollist, des, contrast=my.contrasts[,"t15_PhenDC3vst15_PDS"])


# Display results for top ranked genes
head(camera_t15_DMSOvst0_no, n = 20)
head(camera_t15_PDSvst0_no, n = 20)
head(camera_t15_PhenDC3vst0_no, n = 20)
head(camera_t15_PDSvst15_DMSO, n = 20)
head(camera_t15_PhenDC3vst15_DMSO, n = 20)
head(camera_t15_PhenDC3vst15_PDS, n = 20)


# Make a barcode plot for an example that ranks highly
# t15_DMSOvst0_no
camera_t15_DMSOvst0_no[rownames(camera_t15_DMSOvst0_no) == "ADSL",]
top_noall0_t15_DMSOvst0_no$table[top_noall0_t15_DMSOvst0_no$table$gene == "ADSL",]
barcodeplot(lrt_noall0_t15_DMSOvst0_no$table$logFC, index=genesymbollist[["ADSL"]], labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)

camera_t15_DMSOvst0_no[rownames(camera_t15_DMSOvst0_no) == "DHX36",]
top_noall0_t15_DMSOvst0_no$table[top_noall0_t15_DMSOvst0_no$table$gene == "DHX36",]
barcodeplot(lrt_noall0_t15_DMSOvst0_no$table$logFC, index=genesymbollist[["DHX36"]], labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)


# t15_PDSvst0_no
camera_t15_PDSvst0_no[rownames(camera_t15_PDSvst0_no) == "ADSL",]
top_noall0_t15_PDSvst0_no$table[top_noall0_t15_PDSvst0_no$table$gene == "ADSL",]
barcodeplot(lrt_noall0_t15_PDSvst0_no$table$logFC, index=genesymbollist[["ADSL"]], labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)

camera_t15_PDSvst0_no[rownames(camera_t15_PDSvst0_no) == "SGOL1",]
top_noall0_t15_PDSvst0_no$table[top_noall0_t15_PDSvst0_no$table$gene == "SGOL1",]
par(mfrow=c(3,1))
barcodeplot(lrt_noall0_t15_PDSvst0_no$table$logFC, index=genesymbollist[["SGOL1"]], main="t15_PDSvst0_no SGOL1", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
barcodeplot(lrt_noall0_t15_PDSvst0_no$table$logFC, index=genesymbollist[["ADSL"]], index2=genesymbollist[["SGOL1"]], main="t15_PDSvst0_no ADSL vs. SGOL1", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
barcodeplot(lrt_noall0_t15_DMSOvst0_no$table$logFC, index=genesymbollist[["SGOL1"]], main="t15_DMSOvst0_no SGOL1", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE) # still close

camera_t15_PDSvst0_no[rownames(camera_t15_PDSvst0_no) == "TOMM40",]
top_noall0_t15_PDSvst0_no$table[top_noall0_t15_PDSvst0_no$table$gene == "TOMM40",]
par(mfrow=c(3,1))
barcodeplot(lrt_noall0_t15_PDSvst0_no$table$logFC, index=genesymbollist[["TOMM40"]], main="t15_PDSvst0_no TOMM40", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
barcodeplot(lrt_noall0_t15_PDSvst0_no$table$logFC, index=genesymbollist[["ADSL"]], index2=genesymbollist[["TOMM40"]], main="t15_PDSvst0_no ADSL vs. TOMM40", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
barcodeplot(lrt_noall0_t15_DMSOvst0_no$table$logFC, index=genesymbollist[["TOMM40"]], main="t15_DMSOvst0_no TOMM40", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE) # still close too

camera_t15_PDSvst0_no[rownames(camera_t15_PDSvst0_no) == "AURKB",]
top_noall0_t15_PDSvst0_no$table[top_noall0_t15_PDSvst0_no$table$gene == "AURKB",]
par(mfrow=c(3,1))
barcodeplot(lrt_noall0_t15_PDSvst0_no$table$logFC, index=genesymbollist[["AURKB"]], main="t15_PDSvst0_no AURKB", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
barcodeplot(lrt_noall0_t15_PDSvst0_no$table$logFC, index=genesymbollist[["ADSL"]], index2=genesymbollist[["AURKB"]], main="t15_PDSvst0_no ADSL vs. AURKB", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
barcodeplot(lrt_noall0_t15_DMSOvst0_no$table$logFC, index=genesymbollist[["AURKB"]], main="t15_DMSOvst0_no AURKB", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE) # fairly convincing

camera_t15_PDSvst0_no[rownames(camera_t15_PDSvst0_no) == "DHX36",]
top_noall0_t15_PDSvst0_no$table[top_noall0_t15_PDSvst0_no$table$gene == "DHX36",]
par(mfrow=c(3,1))
barcodeplot(lrt_noall0_t15_PDSvst0_no$table$logFC, index=genesymbollist[["DHX36"]], main="t15_PDSvst0_no DHX36", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
barcodeplot(lrt_noall0_t15_PDSvst0_no$table$logFC, index=genesymbollist[["ADSL"]], index2=genesymbollist[["DHX36"]], main="t15_PDSvst0_no ADSL vs. DHX36", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
barcodeplot(lrt_noall0_t15_DMSOvst0_no$table$logFC, index=genesymbollist[["DHX36"]], main="t15_DMSOvst0_no DHX36", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)


# t15_PhenDC3vst0_no
camera_t15_PhenDC3vst0_no[rownames(camera_t15_PhenDC3vst0_no) == "ADSL",]
top_noall0_t15_PhenDC3vst0_no$table[top_noall0_t15_PhenDC3vst0_no$table$gene == "ADSL",]
barcodeplot(lrt_noall0_t15_PhenDC3vst0_no$table$logFC, index=genesymbollist[["ADSL"]], labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)

camera_t15_PhenDC3vst0_no[rownames(camera_t15_PhenDC3vst0_no) == "OPA1",]
top_noall0_t15_PhenDC3vst0_no$table[top_noall0_t15_PhenDC3vst0_no$table$gene == "OPA1",]
par(mfrow=c(3,1))
barcodeplot(lrt_noall0_t15_PhenDC3vst0_no$table$logFC, index=genesymbollist[["OPA1"]], main="t15_PhenDC3vst0_no OPA1", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
barcodeplot(lrt_noall0_t15_PhenDC3vst0_no$table$logFC, index=genesymbollist[["ADSL"]], index2=genesymbollist[["OPA1"]], main="t15_PhenDC3vst0_no ADSL vs. OPA1", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
barcodeplot(lrt_noall0_t15_DMSOvst0_no$table$logFC, index=genesymbollist[["OPA1"]], main="t15_DMSOvst0_no OPA1", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE) # This makes me think a parametric camera might not be a good idea

camera_t15_PhenDC3vst0_no[rownames(camera_t15_PhenDC3vst0_no) == "GRSF1",]
top_noall0_t15_PhenDC3vst0_no$table[top_noall0_t15_PhenDC3vst0_no$table$gene == "GRSF1",]
par(mfrow=c(4,1))
barcodeplot(lrt_noall0_t15_PhenDC3vst0_no$table$logFC, index=genesymbollist[["GRSF1"]], main="t15_PhenDC3vst0_no GRSF1", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
barcodeplot(lrt_noall0_t15_PhenDC3vst0_no$table$logFC, index=genesymbollist[["ADSL"]], index2=genesymbollist[["GRSF1"]], main="t15_PhenDC3vst0_no ADSL vs. GRSF1", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
barcodeplot(lrt_noall0_t15_DMSOvst0_no$table$logFC, index=genesymbollist[["GRSF1"]], main="t15_DMSOvst0_no GRSF1", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
barcodeplot(lrt_noall0_t15_PDSvst0_no$table$logFC, index=genesymbollist[["GRSF1"]], main="t15_PDSvst0_no GRSF1", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)

camera_t15_PhenDC3vst0_no[rownames(camera_t15_PhenDC3vst0_no) == "DHX36",]
top_noall0_t15_PhenDC3vst0_no$table[top_noall0_t15_PhenDC3vst0_no$table$gene == "DHX36",]
par(mfrow=c(4,1))
barcodeplot(lrt_noall0_t15_PhenDC3vst0_no$table$logFC, index=genesymbollist[["DHX36"]], main="t15_PhenDC3vst0_no DHX36", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
barcodeplot(lrt_noall0_t15_PhenDC3vst0_no$table$logFC, index=genesymbollist[["ADSL"]], index2=genesymbollist[["DHX36"]], main="t15_PhenDC3vst0_no ADSL vs. DHX36", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
barcodeplot(lrt_noall0_t15_DMSOvst0_no$table$logFC, index=genesymbollist[["DHX36"]], main="t15_DMSOvst0_no DHX36", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
barcodeplot(lrt_noall0_t15_PDSvst0_no$table$logFC, index=genesymbollist[["DHX36"]], main="t15_PDSvst0_no DHX36", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)


# t15_PDSvst15_DMSO
camera_t15_PDSvst15_DMSO[rownames(camera_t15_PDSvst15_DMSO) == "ABCB1",]
top_noall0_t15_PDSvst15_DMSO$table[top_noall0_t15_PDSvst15_DMSO$table$gene == "ABCB1",]
barcodeplot(lrt_noall0_t15_PDSvst15_DMSO$table$logFC, index=genesymbollist[["ABCB1"]], main="t15_PDSvst15_DMSO ABCB1", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)

camera_t15_PDSvst15_DMSO[rownames(camera_t15_PDSvst15_DMSO) == "DHX36",]
top_noall0_t15_PDSvst15_DMSO$table[top_noall0_t15_PDSvst15_DMSO$table$gene == "DHX36",]
barcodeplot(lrt_noall0_t15_PDSvst15_DMSO$table$logFC, index=genesymbollist[["ABCB1"]], index2=genesymbollist[["DHX36"]], main="t15_PDSvst15_DMSO ABCB1 vs. DHX36", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)


# t15_PhenDC3vst15_DMSO
camera_t15_PhenDC3vst15_DMSO[rownames(camera_t15_PhenDC3vst15_DMSO) == "LATS2",]
top_noall0_t15_PhenDC3vst15_DMSO$table[top_noall0_t15_PhenDC3vst15_DMSO$table$gene == "LATS2",]
barcodeplot(lrt_noall0_t15_PhenDC3vst15_DMSO$table$logFC, index=genesymbollist[["LATS2"]], main="t15_PhenDC3vst15_DMSO LATS2", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)

camera_t15_PhenDC3vst15_DMSO[rownames(camera_t15_PhenDC3vst15_DMSO) == "DHX36",]
top_noall0_t15_PhenDC3vst15_DMSO$table[top_noall0_t15_PhenDC3vst15_DMSO$table$gene == "DHX36",]
barcodeplot(lrt_noall0_t15_PhenDC3vst15_DMSO$table$logFC, index=genesymbollist[["LATS2"]], index2=genesymbollist[["DHX36"]], main="t15_PhenDC3vst15_DMSO LATS2 vs. DHX36", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)


# t15_PhenDC3vst15_PDS
camera_t15_PhenDC3vst15_PDS[rownames(camera_t15_PhenDC3vst15_PDS) == "TIAL1",]
top_noall0_t15_PhenDC3vst15_PDS$table[top_noall0_t15_PhenDC3vst15_PDS$table$gene == "TIAL1",]
barcodeplot(lrt_noall0_t15_PhenDC3vst15_PDS$table$logFC, index=genesymbollist[["TIAL1"]], main="t15_PhenDC3vst15_PDS TIAL1", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)

camera_t15_PhenDC3vst15_PDS[rownames(camera_t15_PhenDC3vst15_PDS) == "DHX36",]
top_noall0_t15_PhenDC3vst15_PDS$table[top_noall0_t15_PhenDC3vst15_PDS$table$gene == "DHX36",]
barcodeplot(lrt_noall0_t15_PhenDC3vst15_PDS$table$logFC, index=genesymbollist[["TIAL1"]], index2=genesymbollist[["DHX36"]], main="t15_PhenDC3vst15_PDS TIAL1 vs. DHX36", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)


# Run camera for all genes - ranks - unfortunately it is too slow!
camera_t15_DMSOvst0_no <- camera(z_noall0_glm, index = genesymbollist, des, contrast=my.contrasts[,"t15_DMSOvst0_no"], use.ranks=TRUE) ## too slow!


# Explore positive controls (they are supposed to come up in both PhenDC3 and PDS). This might help to to define the right FDR threshold above
positive <- c("DHX36", "BRCA1", "BRCA2", "HNRNPA1", "WRN", "BLM", "RTEL1", "BRIP1", "FANCA", "NCL", "GRSF1", "NSUN5", "DHX9", "DDX3X", "DDX17", "DDX5", "BRAF")

camera_t15_DMSOvst0_no[rownames(camera_t15_DMSOvst0_no) == "DHX36",]
top_noall0_t15_DMSOvst0_no$table[top_noall0_t15_DMSOvst0_no$table$gene == "DHX36",]
camera_t15_PDSvst0_no[rownames(camera_t15_PDSvst0_no) == "DHX36",]
top_noall0_t15_PDSvst0_no$table[top_noall0_t15_PDSvst0_no$table$gene == "DHX36",]
camera_t15_PhenDC3vst0_no[rownames(camera_t15_PhenDC3vst0_no) == "DHX36",]
top_noall0_t15_PhenDC3vst0_no$table[top_noall0_t15_PhenDC3vst0_no$table$gene == "DHX36",]

# barcodeplot by logFC (left) and FDR (right)
for (gene in positive) {
  print(gene)
  pdf(sprintf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161019_barcodeplot_logFCandFDR_%s.pdf", gene))
  par(mfrow=c(3,2))
  barcodeplot(top_noall0_t15_DMSOvst0_no$table[match(rownames(lrt_noall0_t15_DMSOvst0_no$table), top_noall0_t15_DMSOvst0_no$table[,1]),]$logFC, index=genesymbollist[[gene]], main="DMSO vs. t0 logFC", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
  barcodeplot(top_noall0_t15_DMSOvst0_no$table[match(rownames(lrt_noall0_t15_DMSOvst0_no$table), top_noall0_t15_DMSOvst0_no$table[,1]),]$FDR, index=genesymbollist[[gene]], main="DMSO vs. t0 FDR", labels=c("High FDR", "Low FDR"), quantile=c(0.01,1), worm = FALSE)
  barcodeplot(top_noall0_t15_PDSvst0_no$table[match(rownames(lrt_noall0_t15_PDSvst0_no$table), top_noall0_t15_PDSvst0_no$table[,1]),]$logFC, index=genesymbollist[[gene]], main="PDS vs. t0 logFC", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
  barcodeplot(top_noall0_t15_PDSvst0_no$table[match(rownames(lrt_noall0_t15_PDSvst0_no$table), top_noall0_t15_PDSvst0_no$table[,1]),]$FDR, index=genesymbollist[[gene]], main="PDS vs. t0 FDR", labels=c("High FDR", "Low FDR"), quantile=c(0.01,1), worm = FALSE)
  barcodeplot(top_noall0_t15_PhenDC3vst0_no$table[match(rownames(lrt_noall0_t15_PhenDC3vst0_no$table), top_noall0_t15_PhenDC3vst0_no$table[,1]),]$logFC, index=genesymbollist[[gene]], main="PhenDC3 vs. t0 logFC", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
  barcodeplot(top_noall0_t15_PhenDC3vst0_no$table[match(rownames(lrt_noall0_t15_PhenDC3vst0_no$table), top_noall0_t15_PhenDC3vst0_no$table[,1]),]$FDR, index=genesymbollist[[gene]], main="PhenDC3 vs. t0 FDR", labels=c("High FDR", "Low FDR"), quantile=c(0.01,1), worm = FALSE)
  dev.off()
  system(sprintf("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161019_barcodeplot_logFCandFDR_%s.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/figures/barcodeplots/positive_controls", gene))
}



# Explore negative controls (they are supposed not to change).
negative <- c("OR1D2", "OR3A1", "OR2AG1", "PRSS22")

# barcodeplot by logFC (left) and FDR (right)
for (gene in negative) {
  print(gene)
  pdf(sprintf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161019_barcodeplot_logFCandFDR_%s.pdf", gene))
  par(mfrow=c(3,2))
  barcodeplot(top_noall0_t15_DMSOvst0_no$table[match(rownames(lrt_noall0_t15_DMSOvst0_no$table), top_noall0_t15_DMSOvst0_no$table[,1]),]$logFC, index=genesymbollist[[gene]], main="DMSO vs. t0 logFC", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
  barcodeplot(top_noall0_t15_DMSOvst0_no$table[match(rownames(lrt_noall0_t15_DMSOvst0_no$table), top_noall0_t15_DMSOvst0_no$table[,1]),]$FDR, index=genesymbollist[[gene]], main="DMSO vs. t0 FDR", labels=c("High FDR", "Low FDR"), quantile=c(0.01,1), worm = FALSE)
  barcodeplot(top_noall0_t15_PDSvst0_no$table[match(rownames(lrt_noall0_t15_PDSvst0_no$table), top_noall0_t15_PDSvst0_no$table[,1]),]$logFC, index=genesymbollist[[gene]], main="PDS vs. t0 logFC", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
  barcodeplot(top_noall0_t15_PDSvst0_no$table[match(rownames(lrt_noall0_t15_PDSvst0_no$table), top_noall0_t15_PDSvst0_no$table[,1]),]$FDR, index=genesymbollist[[gene]], main="PDS vs. t0 FDR", labels=c("High FDR", "Low FDR"), quantile=c(0.01,1), worm = FALSE)
  barcodeplot(top_noall0_t15_PhenDC3vst0_no$table[match(rownames(lrt_noall0_t15_PhenDC3vst0_no$table), top_noall0_t15_PhenDC3vst0_no$table[,1]),]$logFC, index=genesymbollist[[gene]], main="PhenDC3 vs. t0 logFC", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
  barcodeplot(top_noall0_t15_PhenDC3vst0_no$table[match(rownames(lrt_noall0_t15_PhenDC3vst0_no$table), top_noall0_t15_PhenDC3vst0_no$table[,1]),]$FDR, index=genesymbollist[[gene]], main="PhenDC3 vs. t0 FDR", labels=c("High FDR", "Low FDR"), quantile=c(0.01,1), worm = FALSE)
  dev.off()
  system(sprintf("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161019_barcodeplot_logFCandFDR_%s.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/figures/barcodeplots/negative_controls", gene))
}


# Explore positive controls (they are supposed to come up in both PhenDC3 and PDS). This might help to to define the right FDR threshold above
interesting <- c("ABCB1", "KDM1A")

# barcodeplot by logFC (left) and FDR (right)
for (gene in interesting) {
  print(gene)
  pdf(sprintf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161019_barcodeplot_logFCandFDR_%s.pdf", gene))
  par(mfrow=c(3,2))
  barcodeplot(top_noall0_t15_DMSOvst0_no$table[match(rownames(lrt_noall0_t15_DMSOvst0_no$table), top_noall0_t15_DMSOvst0_no$table[,1]),]$logFC, index=genesymbollist[[gene]], main="DMSO vs. t0 logFC", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
  barcodeplot(top_noall0_t15_DMSOvst0_no$table[match(rownames(lrt_noall0_t15_DMSOvst0_no$table), top_noall0_t15_DMSOvst0_no$table[,1]),]$FDR, index=genesymbollist[[gene]], main="DMSO vs. t0 FDR", labels=c("High FDR", "Low FDR"), quantile=c(0.01,1), worm = FALSE)
  barcodeplot(top_noall0_t15_PDSvst0_no$table[match(rownames(lrt_noall0_t15_PDSvst0_no$table), top_noall0_t15_PDSvst0_no$table[,1]),]$logFC, index=genesymbollist[[gene]], main="PDS vs. t0 logFC", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
  barcodeplot(top_noall0_t15_PDSvst0_no$table[match(rownames(lrt_noall0_t15_PDSvst0_no$table), top_noall0_t15_PDSvst0_no$table[,1]),]$FDR, index=genesymbollist[[gene]], main="PDS vs. t0 FDR", labels=c("High FDR", "Low FDR"), quantile=c(0.01,1), worm = FALSE)
  barcodeplot(top_noall0_t15_PhenDC3vst0_no$table[match(rownames(lrt_noall0_t15_PhenDC3vst0_no$table), top_noall0_t15_PhenDC3vst0_no$table[,1]),]$logFC, index=genesymbollist[[gene]], main="PhenDC3 vs. t0 logFC", labels=c("Positive logFC", "Negative logFC"), quantile=c(-1,1), worm = FALSE)
  barcodeplot(top_noall0_t15_PhenDC3vst0_no$table[match(rownames(lrt_noall0_t15_PhenDC3vst0_no$table), top_noall0_t15_PhenDC3vst0_no$table[,1]),]$FDR, index=genesymbollist[[gene]], main="PhenDC3 vs. t0 FDR", labels=c("High FDR", "Low FDR"), quantile=c(0.01,1), worm = FALSE)
  dev.off()
  system(sprintf("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161019_barcodeplot_logFCandFDR_%s.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/figures/barcodeplots/interesting_examples", gene))
}


# For individual examples, perhaps dot plots can also be generated for specific examples

```


# (C)

```R
# Take object z from (B)
nrow(z$counts) # 110576


# Any shRNAs with zero counts in all libraries?
nrow(z[rowSums(z$counts) == 0, ]$counts) # 41


# Filtering
keep <- rowSums(cpm(z)>0.5) >= 3
z <- z[keep, , keep.lib.sizes=FALSE]
nrow(z$counts) # 101717


# Normalisation
z <- calcNormFactors(z)


# Differential representation analysis.
# Set up design matrix for GLM.
des <- model.matrix(~ 0 + group, data = z$samples)
colnames(des) <- levels(factor(z$samples$group))
des


# Estimate dispersions
z_glm <- estimateDisp(z, des)
sqrt(z_glm$common.disp) # 0.6061857


# Plot BCVs versus abundance
plotBCV(z_glm)


# Fit negative binomial GLM
z_fit <- glmFit(z_glm, des)


# Define matrix of contrasts
my.contrasts <- makeContrasts(
t15_DMSOvst0_no = t15_DMSO - t0_no,
t15_PDSvst0_no = t15_PDS - t0_no,
t15_PhenDC3vst0_no = t15_PhenDC3 - t0_no,
levels=des)


# Comparisons t15 vs t0 and within t15
# Carry out Likelihood ratio tests
lrt_t15_DMSOvst0_no <- glmLRT(z_fit, contrast=my.contrasts[,"t15_DMSOvst0_no"])
lrt_t15_PDSvst0_no <- glmLRT(z_fit, contrast=my.contrasts[,"t15_PDSvst0_no"])
lrt_t15_PhenDC3vst0_no <- glmLRT(z_fit, contrast=my.contrasts[,"t15_PhenDC3vst0_no"])


# Show top ranked hairpins
topTags(lrt_t15_DMSOvst0_no)
topTags(lrt_t15_PDSvst0_no)
topTags(lrt_t15_PhenDC3vst0_no)


####################################################################################
# Genes with least 50% of the shRNAs with FDR < 1e-3 and (logFC < -1 or logFC > 1) #
####################################################################################

# Select genes
top_t15_DMSOvst0_no <- data.table(topTags(lrt_t15_DMSOvst0_no, n=Inf)$table)[FDR < 1e-3 & (logFC < -1 | logFC > 1), .(n_sig_shRNA = uniqueN(shRNA)), by = gene]
top_t15_DMSOvst0_no[, n_tot_shRNA := sapply(as.character(top_t15_DMSOvst0_no$gene), function(x) sum(z$genes$gene == x))]
topids_t15_DMSOvst0_no <- as.character(top_t15_DMSOvst0_no[n_sig_shRNA/n_tot_shRNA >= 0.5]$gene)
length(topids_t15_DMSOvst0_no) # 3

top_t15_PDSvst0_no <- data.table(topTags(lrt_t15_PDSvst0_no, n=Inf)$table)[FDR < 1e-3 & (logFC < -1 | logFC > 1), .(n_sig_shRNA = uniqueN(shRNA)), by = gene]
top_t15_PDSvst0_no[, n_tot_shRNA := sapply(as.character(top_t15_PDSvst0_no$gene), function(x) sum(z$genes$gene == x))]
topids_t15_PDSvst0_no <- as.character(top_t15_PDSvst0_no[n_sig_shRNA/n_tot_shRNA >= 0.5]$gene)
length(topids_t15_PDSvst0_no) # 9

top_t15_PhenDC3vst0_no <- data.table(topTags(lrt_t15_PhenDC3vst0_no, n=Inf)$table)[FDR < 1e-3 & (logFC < -1 | logFC > 1), .(n_sig_shRNA = uniqueN(shRNA)), by = gene]
top_t15_PhenDC3vst0_no[, n_tot_shRNA := sapply(as.character(top_t15_PhenDC3vst0_no$gene), function(x) sum(z$genes$gene == x))]
topids_t15_PhenDC3vst0_no <- as.character(top_t15_PhenDC3vst0_no[n_sig_shRNA/n_tot_shRNA >= 0.5]$gene)
length(topids_t15_PhenDC3vst0_no) # 41


# Intersect lists of topTags
length(intersect(topids_t15_DMSOvst0_no, topids_t15_PDSvst0_no)) # 1
length(intersect(topids_t15_DMSOvst0_no, topids_t15_PhenDC3vst0_no)) # 0
length(intersect(topids_t15_PDSvst0_no, topids_t15_PhenDC3vst0_no)) # 0

top_t15_PDSvst0_no_PDSonly <- top_t15_PDSvst0_no[n_sig_shRNA/n_tot_shRNA >= 0.5][!(topids_t15_PDSvst0_no %in% c(topids_t15_DMSOvst0_no, topids_t15_PhenDC3vst0_no)),]
write.table(top_t15_PDSvst0_no_PDSonly, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161026_top_t15_PDSvst0_no_PDSonly_FDR10minus3_logFC1minus1_50pctshRNAs.txt", quote = FALSE, sep = "\t", row.names = FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161026_top_t15_PDSvst0_no_PDSonly_FDR10minus3_logFC1minus1_50pctshRNAs.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161026/tables/")

top_t15_PhenDC3vst0_no_PhenDC3only <- top_t15_PhenDC3vst0_no[n_sig_shRNA/n_tot_shRNA >= 0.5][!(topids_t15_PhenDC3vst0_no %in% c(topids_t15_DMSOvst0_no, topids_t15_PDSvst0_no)),]
write.table(top_t15_PhenDC3vst0_no_PhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161026_top_t15_PhenDC3vst0_no_PhenDC3only_FDR10minus3_logFC1minus1_50pctshRNAs.txt", quote = FALSE, sep = "\t", row.names = FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161026_top_t15_PhenDC3vst0_no_PhenDC3only_FDR10minus3_logFC1minus1_50pctshRNAs.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161026/tables/")

top_t15_PDSvst0_no_PDSandPhenDC3 <- top_t15_PDSvst0_no[n_sig_shRNA/n_tot_shRNA >= 0.5][(topids_t15_PDSvst0_no %in% topids_t15_PhenDC3vst0_no),] # empty, no overlap
top_t15_PDSvst0_no_PDSandPhenDC3only <- top_t15_PDSvst0_no_PDSandPhenDC3
write.table(top_t15_PDSvst0_no_PDSandPhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161026_top_t15_PDSvst0_no_PDSandPhenDC3only_FDR10minus3_logFC1minus1_50pctshRNAs.txt", quote = FALSE, sep = "\t")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161026_top_t15_PDSvst0_no_PDSandPhenDC3only_FDR10minus3_logFC1minus1_50pctshRNAs.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161026/tables/")


# Venn diagram
venn.plot <- draw.triple.venn(
  area1 = length(topids_t15_PDSvst0_no),
  area2 = length(topids_t15_PhenDC3vst0_no),
  area3 = length(topids_t15_DMSOvst0_no),
  n12 = length(intersect(topids_t15_PDSvst0_no, topids_t15_PhenDC3vst0_no)),
  n23 = length(intersect(topids_t15_PhenDC3vst0_no, topids_t15_DMSOvst0_no)),
  n13 = length(intersect(topids_t15_PDSvst0_no, topids_t15_DMSOvst0_no)),
  n123 = length(intersect(intersect(topids_t15_PDSvst0_no, topids_t15_PhenDC3vst0_no), topids_t15_DMSOvst0_no)),
  category = c(sprintf("PDS (%i)", length(topids_t15_PDSvst0_no)), sprintf("PhenDC3 (%i)", length(topids_t15_PhenDC3vst0_no)), sprintf("DMSO (%i)", length(topids_t15_DMSOvst0_no))),
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


pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161026_t15vst0_venn_FDR10minus3_logFC1minus1_50pctshRNAs.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161026_t15vst0_venn_FDR10minus3_logFC1minus1_50pctshRNAs.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161026/figures/")



####################################################################################
# Genes with least 50% of the shRNAs with FDR < 1e-2 and (logFC < -1 or logFC > 1) #
####################################################################################

# Select genes
top_t15_DMSOvst0_no <- data.table(topTags(lrt_t15_DMSOvst0_no, n=Inf)$table)[FDR < 1e-2 & (logFC < -1 | logFC > 1), .(n_sig_shRNA = uniqueN(shRNA)), by = gene]
top_t15_DMSOvst0_no[, n_tot_shRNA := sapply(as.character(top_t15_DMSOvst0_no$gene), function(x) sum(z$genes$gene == x))]
topids_t15_DMSOvst0_no <- as.character(top_t15_DMSOvst0_no[n_sig_shRNA/n_tot_shRNA >= 0.5]$gene)
length(topids_t15_DMSOvst0_no) # 13

top_t15_PDSvst0_no <- data.table(topTags(lrt_t15_PDSvst0_no, n=Inf)$table)[FDR < 1e-2 & (logFC < -1 | logFC > 1), .(n_sig_shRNA = uniqueN(shRNA)), by = gene]
top_t15_PDSvst0_no[, n_tot_shRNA := sapply(as.character(top_t15_PDSvst0_no$gene), function(x) sum(z$genes$gene == x))]
topids_t15_PDSvst0_no <- as.character(top_t15_PDSvst0_no[n_sig_shRNA/n_tot_shRNA >= 0.5]$gene)
length(topids_t15_PDSvst0_no) # 37

top_t15_PhenDC3vst0_no <- data.table(topTags(lrt_t15_PhenDC3vst0_no, n=Inf)$table)[FDR < 1e-2 & (logFC < -1 | logFC > 1), .(n_sig_shRNA = uniqueN(shRNA)), by = gene]
top_t15_PhenDC3vst0_no[, n_tot_shRNA := sapply(as.character(top_t15_PhenDC3vst0_no$gene), function(x) sum(z$genes$gene == x))]
topids_t15_PhenDC3vst0_no <- as.character(top_t15_PhenDC3vst0_no[n_sig_shRNA/n_tot_shRNA >= 0.5]$gene)
length(topids_t15_PhenDC3vst0_no) # 140


# Intersect lists of topTags
length(intersect(topids_t15_DMSOvst0_no, topids_t15_PDSvst0_no)) # 8
length(intersect(topids_t15_DMSOvst0_no, topids_t15_PhenDC3vst0_no)) # 5
length(intersect(topids_t15_PDSvst0_no, topids_t15_PhenDC3vst0_no)) # 11

top_t15_PDSvst0_no_PDSonly <- top_t15_PDSvst0_no[n_sig_shRNA/n_tot_shRNA >= 0.5][!(topids_t15_PDSvst0_no %in% c(topids_t15_DMSOvst0_no, topids_t15_PhenDC3vst0_no)),]
write.table(top_t15_PDSvst0_no_PDSonly, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161026_top_t15_PDSvst0_no_PDSonly_FDR10minus2_logFC1minus1_50pctshRNAs.txt", quote = FALSE, sep = "\t", row.names = FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161026_top_t15_PDSvst0_no_PDSonly_FDR10minus2_logFC1minus1_50pctshRNAs.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161026/tables/")

top_t15_PhenDC3vst0_no_PhenDC3only <- top_t15_PhenDC3vst0_no[n_sig_shRNA/n_tot_shRNA >= 0.5][!(topids_t15_PhenDC3vst0_no %in% c(topids_t15_DMSOvst0_no, topids_t15_PDSvst0_no)),]
write.table(top_t15_PhenDC3vst0_no_PhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161026_top_t15_PhenDC3vst0_no_PhenDC3only_FDR10minus2_logFC1minus1_50pctshRNAs.txt", quote = FALSE, sep = "\t", row.names = FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161026_top_t15_PhenDC3vst0_no_PhenDC3only_FDR10minus2_logFC1minus1_50pctshRNAs.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161026/tables/")

top_t15_PDSvst0_no_PDSandPhenDC3 <- top_t15_PDSvst0_no[n_sig_shRNA/n_tot_shRNA >= 0.5][(topids_t15_PDSvst0_no %in% topids_t15_PhenDC3vst0_no),]
top_t15_PDSvst0_no_PDSandPhenDC3only <- top_t15_PDSvst0_no_PDSandPhenDC3[!(top_t15_PDSvst0_no_PDSandPhenDC3$gene %in% topids_t15_DMSOvst0_no),]
write.table(top_t15_PDSvst0_no_PDSandPhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161026_top_t15_PDSvst0_no_PDSandPhenDC3only_FDR10minus2_logFC1minus1_50pctshRNAs.txt", quote = FALSE, sep = "\t")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161026_top_t15_PDSvst0_no_PDSandPhenDC3only_FDR10minus2_logFC1minus1_50pctshRNAs.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161026/tables/")


# Venn diagram
venn.plot <- draw.triple.venn(
  area1 = length(topids_t15_PDSvst0_no),
  area2 = length(topids_t15_PhenDC3vst0_no),
  area3 = length(topids_t15_DMSOvst0_no),
  n12 = length(intersect(topids_t15_PDSvst0_no, topids_t15_PhenDC3vst0_no)),
  n23 = length(intersect(topids_t15_PhenDC3vst0_no, topids_t15_DMSOvst0_no)),
  n13 = length(intersect(topids_t15_PDSvst0_no, topids_t15_DMSOvst0_no)),
  n123 = length(intersect(intersect(topids_t15_PDSvst0_no, topids_t15_PhenDC3vst0_no), topids_t15_DMSOvst0_no)),
  category = c(sprintf("PDS (%i)", length(topids_t15_PDSvst0_no)), sprintf("PhenDC3 (%i)", length(topids_t15_PhenDC3vst0_no)), sprintf("DMSO (%i)", length(topids_t15_DMSOvst0_no))),
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


pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161026_t15vst0_venn_FDR10minus2_logFC1minus1_50pctshRNAs.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161026_t15vst0_venn_FDR10minus2_logFC1minus1_50pctshRNAs.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161026/figures/")

```


# (D)

Under construction



# (E)

## Analysis:

```R

# Sum counts of shRNAs by gene, using z matrix from (B)
data_by_gene <- data.table(z$counts)
data_by_gene[, gene := as.character(sapply(rownames(z$counts), function(x) unlist(strsplit(x, "__"))[1]))]
data_by_gene <- data.frame(data_by_gene[, lapply(.SD, sum), keyby = gene])
rownames(data_by_gene) <- data_by_gene[,1]
data_by_gene <- data_by_gene[, 2:13]
dim(data_by_gene) # 18923    12


# Change to a DGE object
groups <- c("t0_no", "t0_no", "t0_no", "t15_DMSO", "t15_DMSO", "t15_DMSO", "t15_PDS", "t15_PDS", "t15_PDS", "t15_PhenDC3", "t15_PhenDC3", "t15_PhenDC3")
z_gene <- DGEList(counts = data_by_gene, group = groups)


# Which genes have more reads on average?
head(sort(rowSums(z_gene$counts), decreasing = T))
# HEPACAM  TAS2R42    TAOK1     MKKS    CLCA2 KIAA1217
# 3957757  1495175  1309874  1103734   903525   848995
sum(z_gene$samples$lib.size) # 984949266
100*(3957757/984949266) # 0.4018234

par(mfrow=c(2,1))
barplot(colSums(z_gene$counts), las=2, main="Counts per index", cex.names=0.5, cex.axis=0.8)
barplot(rowSums(z_gene$counts), las=2, main="Counts per hairpin", cex.names=0.5, cex.axis=0.8, axisnames = FALSE)


# Any gene with zero counts in all libraries?
z_gene[rowSums(z_gene$counts) == 0, ] # LOC389676


# Filtering
keep <- rowSums(cpm(z_gene)>0.5) >= 6
z_gene <- z_gene[keep, , keep.lib.sizes=FALSE] # 18832


# Normalisation
z_gene <- calcNormFactors(z_gene)


# Make MDS plots to visualise relationships between replicate samples
pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161022_mds.pdf", width = 16/2.54, height = 16/2.54)
g <- plotMDS(z_gene, labels = c(rep("t0", 3), rep("DMSO", 3), rep("PDS", 3), rep("PhenDC3", 3)), col = c(rep("black", 3), rep("deepskyblue3", 3), rep("red3", 3), rep("forestgreen", 3)), xlab = "Dimension 1", ylab = "Dimension 2",  gene.selection = "common")
legend("bottomleft", legend=c("t0", "DMSO", "PDS", "PhenDC3"), col=c("black", "deepskyblue3", "red3", "forestgreen"), pch=15)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161022_mds.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161022/figures")


# Differential representation analysis.
# Set up design matrix for GLM.
des <- model.matrix(~ 0 + group, data = z_gene$samples)
colnames(des) <- levels(factor(z_gene$samples$group))
des


# Estimate dispersions
z_gene_glm <- estimateDisp(z_gene, des)
sqrt(z_gene_glm$common.disp) # 0.3003741


# Plot BCVs versus abundance
plotBCV(z_gene_glm)


# Fit negative binomial GLM
z_gene_fit <- glmFit(z_gene_glm, des)


# Define matrix of contrasts
my.contrasts <- makeContrasts(
t15_DMSOvst0_no = t15_DMSO - t0_no,
t15_PDSvst0_no = t15_PDS - t0_no,
t15_PhenDC3vst0_no = t15_PhenDC3 - t0_no,
t15_PDSvst15_DMSO = t15_PDS - t15_DMSO,
t15_PhenDC3vst15_DMSO = t15_PhenDC3 - t15_DMSO,
t15_PhenDC3vst15_PDS = t15_PhenDC3 - t15_PDS,
levels=des)


# Comparisons t15 vs t0 and within t15
# Carry out Likelihood ratio tests
lrt_t15_DMSOvst0_no <- glmLRT(z_gene_fit, contrast=my.contrasts[,"t15_DMSOvst0_no"])
lrt_t15_PDSvst0_no <- glmLRT(z_gene_fit, contrast=my.contrasts[,"t15_PDSvst0_no"])
lrt_t15_PhenDC3vst0_no <- glmLRT(z_gene_fit, contrast=my.contrasts[,"t15_PhenDC3vst0_no"])
lrt_t15_PDSvst15_DMSO <- glmLRT(z_gene_fit, contrast=my.contrasts[,"t15_PDSvst15_DMSO"])
lrt_t15_PhenDC3vst15_DMSO <- glmLRT(z_gene_fit, contrast=my.contrasts[,"t15_PhenDC3vst15_DMSO"])
lrt_t15_PhenDC3vst15_PDS <- glmLRT(z_gene_fit, contrast=my.contrasts[,"t15_PhenDC3vst15_PDS"])


# Show top ranked hairpins
topTags(lrt_t15_DMSOvst0_no)
topTags(lrt_t15_PDSvst0_no)
topTags(lrt_t15_PhenDC3vst0_no)
topTags(lrt_t15_PDSvst15_DMSO)
topTags(lrt_t15_PhenDC3vst15_DMSO)
topTags(lrt_t15_PhenDC3vst15_PDS)


# Select genes with FDR < 0.0001 to highlight on plot
thresh <- 0.0001 # section 4.6.5 Differential representation analysis in Chen2015 (https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)

top_t15_DMSOvst0_no <- topTags(lrt_t15_DMSOvst0_no, n=Inf)
topids_t15_DMSOvst0_no <- rownames(top_t15_DMSOvst0_no$table[top_t15_DMSOvst0_no$table$FDR < thresh,])
length(topids_t15_DMSOvst0_no) # 65

top_t15_PDSvst0_no <- topTags(lrt_t15_PDSvst0_no, n=Inf)
topids_t15_PDSvst0_no <- rownames(top_t15_PDSvst0_no$table[top_t15_PDSvst0_no$table$FDR < thresh,])
length(topids_t15_PDSvst0_no) # 135

top_t15_PhenDC3vst0_no <- topTags(lrt_t15_PhenDC3vst0_no, n=Inf)
topids_t15_PhenDC3vst0_no <- rownames(top_t15_PhenDC3vst0_no$table[top_t15_PhenDC3vst0_no$table$FDR < thresh,])
length(topids_t15_PhenDC3vst0_no) # 315

top_t15_PDSvst15_DMSO <- topTags(lrt_t15_PDSvst15_DMSO, n=Inf)
topids_t15_PDSvst15_DMSO <- rownames(top_t15_PDSvst15_DMSO$table[top_t15_PDSvst15_DMSO$table$FDR < thresh,])
length(topids_t15_PDSvst15_DMSO) # 21

top_t15_PhenDC3vst15_DMSO <- topTags(lrt_t15_PhenDC3vst15_DMSO, n=Inf)
topids_t15_PhenDC3vst15_DMSO <- rownames(top_t15_PhenDC3vst15_DMSO$table[top_t15_PhenDC3vst15_DMSO$table$FDR < thresh,])
length(topids_t15_PhenDC3vst15_DMSO) # 320

top_t15_PhenDC3vst15_PDS <- topTags(lrt_t15_PhenDC3vst15_PDS, n=Inf)
topids_t15_PhenDC3vst15_PDS <- rownames(top_t15_PhenDC3vst15_PDS$table[top_t15_PhenDC3vst15_PDS$table$FDR < thresh,])
length(topids_t15_PhenDC3vst15_PDS) # 308


# Write LogFC and FDR tables
write.table(data.table(as.data.frame(top_t15_DMSOvst0_no), keep.rownames=TRUE)[order(-logFC),.(rn, logFC, logCPM, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_t15_DMSOvst0_logFCdecreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_t15_DMSOvst0_logFCdecreasing_FDR.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161022/tables/")

write.table(data.table(as.data.frame(top_t15_PDSvst0_no), keep.rownames=TRUE)[order(-logFC),.(rn, logFC, logCPM, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_t15_PDSvst0_logFCdecreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_t15_PDSvst0_logFCdecreasing_FDR.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161022/tables/")

write.table(data.table(as.data.frame(top_t15_PhenDC3vst0_no), keep.rownames=TRUE)[order(-logFC),.(rn, logFC, logCPM, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_t15_PhenDC3vst0_logFCdecreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_t15_PhenDC3vst0_logFCdecreasing_FDR.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161022/tables/")

write.table(data.table(as.data.frame(top_t15_PDSvst15_DMSO), keep.rownames=TRUE)[order(-logFC),.(rn, logFC, logCPM, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_t15_PDSvst15_DMSO_logFCdecreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_t15_PDSvst15_DMSO_logFCdecreasing_FDR.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161022/tables/")

write.table(data.table(as.data.frame(top_t15_PhenDC3vst15_DMSO), keep.rownames=TRUE)[order(-logFC),.(rn, logFC, logCPM, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_t15_PhenDC3vst15_DMSO_logFCdecreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_t15_PhenDC3vst15_DMSO_logFCdecreasing_FDR.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161022/tables/")

write.table(data.table(as.data.frame(top_t15_PhenDC3vst15_PDS), keep.rownames=TRUE)[order(-logFC),.(rn, logFC, logCPM, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_t15_PhenDC3vst15_PDS_logFCdecreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_t15_PhenDC3vst15_PDS_logFCdecreasing_FDR.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161022/tables/")


# Plot logFC versus logCPM
plotSmear(lrt_t15_DMSOvst0_no, de.tags=topids_t15_DMSOvst0_no, main = "t15_DMSO vs t0_no")

plotSmear(lrt_t15_PDSvst0_no, de.tags=topids_t15_PDSvst0_no, main = "t15_PDS vs t0_no")
z_gene_t15_PDSvst0_no_cpm <- data.table(cpm(z_gene)[, grepl("t15_PDS|t0_no", colnames(cpm(z_gene)))], keep.rownames=TRUE)
z_gene_t15_PDSvst0_no_cpm[, log2cpmt0 := (log2(t0_no_1) + log2(t0_no_2) + log2(t0_no_3))/3]
z_gene_t15_PDSvst0_no_cpm[, log2FC := lrt_t15_PDSvst0_no$table$logFC]
z_gene_t15_PDSvst0_no_cpm[, hits := ifelse(rn %in% topids_t15_PDSvst0_no, ifelse(rn %in% topids_t15_DMSOvst0_no, "PDS and DMSO", "PDS only"), "not differentially abundant")]
z_gene_t15_PDSvst0_no_cpm[, order := ifelse(hits == "not differentially abundant", 1, 2)]

gg <- ggplot(z_gene_t15_PDSvst0_no_cpm[log2cpmt0 != -Inf][order(order)], aes(x = log2cpmt0, y = log2FC, colour = hits)) +
geom_point(size=0.5) +
xlab(expression("log"[2]*"(t0)")) +
ylab(expression("log"[2]*"(t15/t0)")) +
theme_bw() +
scale_colour_manual(name="",values = c("grey", "deepskyblue3", "red3"))
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161022_log2t0_log2FC_PDS_DMSO.pdf", width = 16/2.54, height = 10/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161022_log2t0_log2FC_PDS_DMSO.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161022/figures/")

plotSmear(lrt_t15_PhenDC3vst0_no, de.tags=topids_t15_PhenDC3vst0_no, main = "t15_PhenDC3 vs t0_no")
z_gene_t15_PhenDC3vst0_no_cpm <- data.table(cpm(z_gene)[, grepl("t15_PhenDC3|t0_no", colnames(cpm(z_gene)))], keep.rownames=TRUE)
z_gene_t15_PhenDC3vst0_no_cpm[, log2cpmt0 := (log2(t0_no_1) + log2(t0_no_2) + log2(t0_no_3))/3]
z_gene_t15_PhenDC3vst0_no_cpm[, log2FC := lrt_t15_PhenDC3vst0_no$table$logFC]
z_gene_t15_PhenDC3vst0_no_cpm[, hits := ifelse(rn %in% topids_t15_PhenDC3vst0_no, ifelse(rn %in% topids_t15_DMSOvst0_no, "PhenDC3 and DMSO", "PhenDC3 only"), "not differentially abundant")]
z_gene_t15_PhenDC3vst0_no_cpm[, order := ifelse(hits == "not differentially abundant", 1, 2)]

gg <- ggplot(z_gene_t15_PhenDC3vst0_no_cpm[log2cpmt0 != -Inf][order(order)], aes(x = log2cpmt0, y = log2FC, colour = hits)) +
geom_point(size=0.5) +
xlab(expression("log"[2]*"(t0)")) +
ylab(expression("log"[2]*"(t15/t0)")) +
theme_bw() +
scale_colour_manual(name="",values = c("grey", "deepskyblue3", "forestgreen"))
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161022_log2t0_log2FC_PhenDC3_DMSO.pdf", width = 16/2.54, height = 10/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161022_log2t0_log2FC_PhenDC3_DMSO.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161022/figures/")

plotSmear(lrt_t15_PDSvst15_DMSO, de.tags=topids_t15_PDSvst15_DMSO, main = "t15_PDS vs t15_DMSO")
plotSmear(lrt_t15_PhenDC3vst15_DMSO, de.tags=topids_t15_PhenDC3vst15_DMSO, main = "t15_PhenDC3 vs t15_DMSO")
plotSmear(lrt_t15_PhenDC3vst15_PDS, de.tags=topids_t15_PhenDC3vst15_PDS, main = "t15_PhenDC3 vs t15_PDS")


# Intersect lists of topTags
length(intersect(topids_t15_DMSOvst0_no, topids_t15_PDSvst0_no)) # 39
length(intersect(topids_t15_DMSOvst0_no, topids_t15_PhenDC3vst0_no)) # 21
length(intersect(topids_t15_PDSvst0_no, topids_t15_PhenDC3vst0_no)) # 38

write(topids_t15_DMSOvst0_no, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_topids_t15_DMSOvst0_no.txt")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_topids_t15_DMSOvst0_no.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161022/tables/")
write(topids_t15_PDSvst0_no, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_topids_t15_PDSvst0_no.txt")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_topids_t15_PDSvst0_no.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161022/tables/")
write(topids_t15_PhenDC3vst0_no, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_topids_t15_PhenDC3vst0_no.txt")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_topids_t15_PhenDC3vst0_no.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161022/tables/")

top_t15_PDSvst0_no_PDSonly <- top_t15_PDSvst0_no[1:length(topids_t15_PDSvst0_no),][!(topids_t15_PDSvst0_no %in% c(topids_t15_DMSOvst0_no, topids_t15_PhenDC3vst0_no)),]
write.table(top_t15_PDSvst0_no_PDSonly, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_top_t15_PDSvst0_no_PDSonly.txt", quote = FALSE, sep = "\t")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_top_t15_PDSvst0_no_PDSonly.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161022/tables/")

top_t15_PhenDC3vst0_no_PhenDC3only <- top_t15_PhenDC3vst0_no[1:length(topids_t15_PhenDC3vst0_no),][!(topids_t15_PhenDC3vst0_no %in% c(topids_t15_DMSOvst0_no, topids_t15_PDSvst0_no)),]
write.table(top_t15_PhenDC3vst0_no_PhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_top_t15_PhenDC3vst0_no_PhenDC3only.txt", quote = FALSE, sep = "\t")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_top_t15_PhenDC3vst0_no_PhenDC3only.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161022/tables/")

top_t15_PDSvst0_no_PDSandPhenDC3 <- top_t15_PDSvst0_no[1:length(topids_t15_PDSvst0_no),][(topids_t15_PDSvst0_no %in% topids_t15_PhenDC3vst0_no),]
top_t15_PDSvst0_no_PDSandPhenDC3only <- top_t15_PDSvst0_no_PDSandPhenDC3[!(rownames(top_t15_PDSvst0_no_PDSandPhenDC3) %in% topids_t15_DMSOvst0_no),]
write.table(top_t15_PDSvst0_no_PDSandPhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_top_t15_PDSvst0_no_PDSandPhenDC3only.txt", quote = FALSE, sep = "\t")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161022_top_t15_PDSvst0_no_PDSandPhenDC3only.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161022/tables/")

top_t15_DMSOvst0_no_DMSOonly <- top_t15_DMSOvst0_no[1:length(topids_t15_DMSOvst0_no),][!(topids_t15_DMSOvst0_no %in% c(topids_t15_PDSvst0_no, topids_t15_PhenDC3vst0_no)),]


# Venn diagram
venn.plot <- draw.triple.venn(
  area1 = length(topids_t15_PDSvst0_no),
  area2 = length(topids_t15_PhenDC3vst0_no),
  area3 = length(topids_t15_DMSOvst0_no),
  n12 = length(intersect(topids_t15_PDSvst0_no, topids_t15_PhenDC3vst0_no)),
  n23 = length(intersect(topids_t15_PhenDC3vst0_no, topids_t15_DMSOvst0_no)),
  n13 = length(intersect(topids_t15_PDSvst0_no, topids_t15_DMSOvst0_no)),
  n123 = length(intersect(intersect(topids_t15_PDSvst0_no, topids_t15_PhenDC3vst0_no), topids_t15_DMSOvst0_no)),
  category = c(sprintf("PDS (%i)", length(topids_t15_PDSvst0_no)), sprintf("PhenDC3 (%i)", length(topids_t15_PhenDC3vst0_no)), sprintf("DMSO (%i)", length(topids_t15_DMSOvst0_no))),
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


pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161022_t15vst0_venn.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161022_t15vst0_venn.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161022/figures/")

```



# (F)

## Analysis:

```R

# Sum counts of shRNAs by gene, using filtered and normalised z matrix from (C)
data_by_gene <- data.table(z$counts)
data_by_gene[, gene := as.character(sapply(rownames(z$counts), function(x) unlist(strsplit(x, "__"))[1]))]
n <- data_by_gene[, .N, keyby = gene]$N
data_by_gene <- data.frame(data_by_gene[, lapply(.SD, sum), keyby = gene])
rownames(data_by_gene) <- data_by_gene[,1]
data_by_gene <- data_by_gene[, 2:13]
dim(data_by_gene) # 18862   12


# Change to a DGE object
groups <- c("t0_no", "t0_no", "t0_no", "t15_DMSO", "t15_DMSO", "t15_DMSO", "t15_PDS", "t15_PDS", "t15_PDS", "t15_PhenDC3", "t15_PhenDC3", "t15_PhenDC3")
genes_dataframe <- data.frame(gene = rownames(data_by_gene), number_shRNAs = n)
z_gene <- DGEList(counts = data_by_gene, group = groups, genes = genes_dataframe)


# Which genes have more reads on average?
head(sort(rowSums(z_gene$counts), decreasing = T))
# HEPACAM  TAS2R42    TAOK1     MKKS    CLCA2 KIAA1217
# 3957757  1495175  1309874  1103734   903525   848995
sum(z_gene$samples$lib.size) # 983742376

par(mfrow=c(2,1))
barplot(colSums(z_gene$counts), las=2, main="Counts per index", cex.names=0.5, cex.axis=0.8)
barplot(rowSums(z_gene$counts), las=2, main="Counts per hairpin", cex.names=0.5, cex.axis=0.8, axisnames = FALSE)


# Any gene with zero counts in all libraries?
z_gene[rowSums(z_gene$counts) == 0, ] # NO


# Filtering
keep <- rowSums(cpm(z_gene)>0.5) >= 6
z_gene <- z_gene[keep, , keep.lib.sizes=FALSE]
nrow(z_gene$counts) # 18817


# Normalisation
z_gene <- calcNormFactors(z_gene)


# Make MDS plots to visualise relationships between replicate samples
pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161027_mds.pdf", width = 16/2.54, height = 16/2.54)
g <- plotMDS(z_gene, labels = c(rep("t0", 3), rep("DMSO", 3), rep("PDS", 3), rep("PhenDC3", 3)), col = c(rep("black", 3), rep("deepskyblue3", 3), rep("red3", 3), rep("forestgreen", 3)), xlab = "Dimension 1", ylab = "Dimension 2",  gene.selection = "common")
legend("bottomleft", legend=c("t0", "DMSO", "PDS", "PhenDC3"), col=c("black", "deepskyblue3", "red3", "forestgreen"), pch=15)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161027_mds.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161027/figures")


# Differential representation analysis.
# Set up design matrix for GLM.
des <- model.matrix(~ 0 + group, data = z_gene$samples)
colnames(des) <- levels(factor(z_gene$samples$group))
des


# Estimate dispersions
z_gene_glm <- estimateDisp(z_gene, des)
sqrt(z_gene_glm$common.disp) # 0.3024315


# Plot BCVs versus abundance
plotBCV(z_gene_glm)


# Fit negative binomial GLM
z_gene_fit <- glmFit(z_gene_glm, des)


# Define matrix of contrasts
my.contrasts <- makeContrasts(
t15_DMSOvst0_no = t15_DMSO - t0_no,
t15_PDSvst0_no = t15_PDS - t0_no,
t15_PhenDC3vst0_no = t15_PhenDC3 - t0_no,
levels=des)


# Comparisons t15 vs t0 and within t15
# Carry out Likelihood ratio tests
lrt_gene_t15_DMSOvst0_no <- glmLRT(z_gene_fit, contrast=my.contrasts[,"t15_DMSOvst0_no"])
lrt_gene_t15_PDSvst0_no <- glmLRT(z_gene_fit, contrast=my.contrasts[,"t15_PDSvst0_no"])
lrt_gene_t15_PhenDC3vst0_no <- glmLRT(z_gene_fit, contrast=my.contrasts[,"t15_PhenDC3vst0_no"])


# Show top ranked hairpins
topTags(lrt_gene_t15_DMSOvst0_no)
topTags(lrt_gene_t15_PDSvst0_no)
topTags(lrt_gene_t15_PhenDC3vst0_no)


#######################################################
# Genes with FDR < 1e-3 and (logFC < -1 or logFC > 1) #
#######################################################

# Select genes
top_gene_t15_DMSOvst0_no <- data.table(topTags(lrt_gene_t15_DMSOvst0_no, n=Inf)$table)[FDR < 1e-3 & (logFC < -1 | logFC > 1)][,.(gene, number_shRNAs, logFC, logCPM, FDR)]
top_gene_t15_DMSOvst0_no[, type := ifelse(logFC > 0, "resistant", ifelse(logFC < 0, "sensitiser", "none"))]
setcolorder(top_gene_t15_DMSOvst0_no, c("gene", "number_shRNAs", "logFC", "type", "logCPM", "FDR"))
topids_gene_t15_DMSOvst0_no <- as.character(top_gene_t15_DMSOvst0_no$gene)
length(topids_gene_t15_DMSOvst0_no) # 157

top_gene_t15_PDSvst0_no <- data.table(topTags(lrt_gene_t15_PDSvst0_no, n=Inf)$table)[FDR < 1e-3 & (logFC < -1 | logFC > 1)][,.(gene, number_shRNAs, logFC, logCPM, FDR)]
top_gene_t15_PDSvst0_no[, type := ifelse(logFC > 0, "resistant", ifelse(logFC < 0, "sensitiser", "none"))]
setcolorder(top_gene_t15_PDSvst0_no, c("gene", "number_shRNAs", "logFC", "type", "logCPM", "FDR"))
topids_gene_t15_PDSvst0_no <- as.character(top_gene_t15_PDSvst0_no$gene)
length(topids_gene_t15_PDSvst0_no) # 258

top_gene_t15_PhenDC3vst0_no <- data.table(topTags(lrt_gene_t15_PhenDC3vst0_no, n=Inf)$table)[FDR < 1e-3 & (logFC < -1 | logFC > 1)][,.(gene, number_shRNAs, logFC, logCPM, FDR)]
top_gene_t15_PhenDC3vst0_no[, type := ifelse(logFC > 0, "resistant", ifelse(logFC < 0, "sensitiser", "none"))]
setcolorder(top_gene_t15_PhenDC3vst0_no, c("gene", "number_shRNAs", "logFC", "type", "logCPM", "FDR"))
topids_gene_t15_PhenDC3vst0_no <- as.character(top_gene_t15_PhenDC3vst0_no$gene)
length(topids_gene_t15_PhenDC3vst0_no) # 505


# Plot logFC versus logCPM
plotSmear(lrt_gene_t15_PDSvst0_no, de.tags=topids_gene_t15_PDSvst0_no, main = "t15_PDS vs t0_no")
z_gene_t15_PDSvst0_no_cpm <- data.table(cpm(z_gene)[, grepl("t15_PDS|t0_no", colnames(cpm(z_gene)))], keep.rownames=TRUE)
z_gene_t15_PDSvst0_no_cpm[, log2cpmt0 := (log2(t0_no_1) + log2(t0_no_2) + log2(t0_no_3))/3]
z_gene_t15_PDSvst0_no_cpm[, log2FC := lrt_gene_t15_PDSvst0_no$table$logFC]
z_gene_t15_PDSvst0_no_cpm[, hits := ifelse(rn %in% topids_gene_t15_PDSvst0_no, ifelse(rn %in% topids_gene_t15_DMSOvst0_no, "PDS and DMSO", "PDS only"), "not differentially abundant")]
z_gene_t15_PDSvst0_no_cpm[, order := ifelse(hits == "not differentially abundant", 1, 2)]

gg <- ggplot(z_gene_t15_PDSvst0_no_cpm[log2cpmt0 != -Inf][order(order)], aes(x = log2cpmt0, y = log2FC, colour = hits)) +
geom_point(size=0.5) +
xlab(expression("log"[2]*"(t0)")) +
ylab(expression("log"[2]*"(t15/t0)")) +
theme_bw() +
scale_colour_manual(name="",values = c("grey", "deepskyblue3", "red3"))
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161027_log2t0_log2FC_PDS_DMSO_FDR10minus3_logFC1minus1.pdf", width = 16/2.54, height = 10/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161027_log2t0_log2FC_PDS_DMSO_FDR10minus3_logFC1minus1.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161027/figures/")

plotSmear(lrt_gene_t15_PhenDC3vst0_no, de.tags=topids_gene_t15_PhenDC3vst0_no, main = "t15_PhenDC3 vs t0_no")
z_gene_t15_PhenDC3vst0_no_cpm <- data.table(cpm(z_gene)[, grepl("t15_PhenDC3|t0_no", colnames(cpm(z_gene)))], keep.rownames=TRUE)
z_gene_t15_PhenDC3vst0_no_cpm[, log2cpmt0 := (log2(t0_no_1) + log2(t0_no_2) + log2(t0_no_3))/3]
z_gene_t15_PhenDC3vst0_no_cpm[, log2FC := lrt_gene_t15_PhenDC3vst0_no$table$logFC]
z_gene_t15_PhenDC3vst0_no_cpm[, hits := ifelse(rn %in% topids_gene_t15_PhenDC3vst0_no, ifelse(rn %in% topids_gene_t15_DMSOvst0_no, "PhenDC3 and DMSO", "PhenDC3 only"), "not differentially abundant")]
z_gene_t15_PhenDC3vst0_no_cpm[, order := ifelse(hits == "not differentially abundant", 1, 2)]

gg <- ggplot(z_gene_t15_PhenDC3vst0_no_cpm[log2cpmt0 != -Inf][order(order)], aes(x = log2cpmt0, y = log2FC, colour = hits)) +
geom_point(size=0.5) +
xlab(expression("log"[2]*"(t0)")) +
ylab(expression("log"[2]*"(t15/t0)")) +
theme_bw() +
scale_colour_manual(name="",values = c("grey", "deepskyblue3", "forestgreen"))
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161027_log2t0_log2FC_PhenDC3_DMSO_FDR10minus3_logFC1minus1.pdf", width = 16/2.54, height = 10/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161027_log2t0_log2FC_PhenDC3_DMSO_FDR10minus3_logFC1minus1.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161027/figures/")


# Intersect lists of topTags
length(intersect(topids_gene_t15_DMSOvst0_no, topids_gene_t15_PDSvst0_no)) # 80
length(intersect(topids_gene_t15_DMSOvst0_no, topids_gene_t15_PhenDC3vst0_no)) # 55
length(intersect(topids_gene_t15_PDSvst0_no, topids_gene_t15_PhenDC3vst0_no)) # 72

top_gene_t15_PDSvst0_no_PDSonly <- top_gene_t15_PDSvst0_no[!(topids_gene_t15_PDSvst0_no %in% c(topids_gene_t15_DMSOvst0_no, topids_gene_t15_PhenDC3vst0_no)),]
write.table(top_gene_t15_PDSvst0_no_PDSonly, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161027_top_t15_PDSvst0_no_PDSonly_FDR10minus3_logFC1minus1.txt", quote = FALSE, sep = "\t", row.names = FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161027_top_t15_PDSvst0_no_PDSonly_FDR10minus3_logFC1minus1.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161027/tables/")

top_gene_t15_PhenDC3vst0_no_PhenDC3only <- top_gene_t15_PhenDC3vst0_no[!(topids_gene_t15_PhenDC3vst0_no %in% c(topids_gene_t15_DMSOvst0_no, topids_gene_t15_PDSvst0_no)),]
write.table(top_gene_t15_PhenDC3vst0_no_PhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161027_top_t15_PhenDC3vst0_no_PhenDC3only_FDR10minus3_logFC1minus1.txt", quote = FALSE, sep = "\t", row.names = FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161027_top_t15_PhenDC3vst0_no_PhenDC3only_FDR10minus3_logFC1minus1.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161027/tables/")

top_gene_t15_PDSvst0_no_PDSandPhenDC3 <- top_gene_t15_PDSvst0_no[(topids_gene_t15_PDSvst0_no %in% topids_gene_t15_PhenDC3vst0_no),]
top_gene_t15_PDSvst0_no_PDSandPhenDC3only <- top_gene_t15_PDSvst0_no_PDSandPhenDC3[!(gene %in% topids_gene_t15_DMSOvst0_no),]
write.table(top_gene_t15_PDSvst0_no_PDSandPhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161027_top_t15_PDSvst0_no_PDSandPhenDC3only_FDR10minus3_logFC1minus1.txt", quote = FALSE, sep = "\t", row.names = FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161027_top_t15_PDSvst0_no_PDSandPhenDC3only_FDR10minus3_logFC1minus1.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161027/tables/")


# Venn diagram
venn.plot <- draw.triple.venn(
  area1 = length(topids_gene_t15_PDSvst0_no),
  area2 = length(topids_gene_t15_PhenDC3vst0_no),
  area3 = length(topids_gene_t15_DMSOvst0_no),
  n12 = length(intersect(topids_gene_t15_PDSvst0_no, topids_gene_t15_PhenDC3vst0_no)),
  n23 = length(intersect(topids_gene_t15_PhenDC3vst0_no, topids_gene_t15_DMSOvst0_no)),
  n13 = length(intersect(topids_gene_t15_PDSvst0_no, topids_gene_t15_DMSOvst0_no)),
  n123 = length(intersect(intersect(topids_gene_t15_PDSvst0_no, topids_gene_t15_PhenDC3vst0_no), topids_gene_t15_DMSOvst0_no)),
  category = c(sprintf("PDS (%i)", length(topids_gene_t15_PDSvst0_no)), sprintf("PhenDC3 (%i)", length(topids_gene_t15_PhenDC3vst0_no)), sprintf("DMSO (%i)", length(topids_gene_t15_DMSOvst0_no))),
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


pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161027_t15vst0_venn_FDR10minus3_logFC1minus1.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161027_t15vst0_venn_FDR10minus3_logFC1minus1.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161027/figures/")



#######################################################
# Genes with FDR < 1e-2 and (logFC < -1 or logFC > 1) #
#######################################################

# Select genes
top_gene_t15_DMSOvst0_no <- data.table(topTags(lrt_gene_t15_DMSOvst0_no, n=Inf)$table)[FDR < 1e-2 & (logFC < -1 | logFC > 1)][,.(gene, number_shRNAs, logFC, logCPM, FDR)]
top_gene_t15_DMSOvst0_no[, type := ifelse(logFC > 0, "resistant", ifelse(logFC < 0, "sensitiser", "none"))]
setcolorder(top_gene_t15_DMSOvst0_no, c("gene", "number_shRNAs", "logFC", "type", "logCPM", "FDR"))
topids_gene_t15_DMSOvst0_no <- as.character(top_gene_t15_DMSOvst0_no$gene)
length(topids_gene_t15_DMSOvst0_no) # 332

top_gene_t15_PDSvst0_no <- data.table(topTags(lrt_gene_t15_PDSvst0_no, n=Inf)$table)[FDR < 1e-2 & (logFC < -1 | logFC > 1)][,.(gene, number_shRNAs, logFC, logCPM, FDR)]
top_gene_t15_PDSvst0_no[, type := ifelse(logFC > 0, "resistant", ifelse(logFC < 0, "sensitiser", "none"))]
setcolorder(top_gene_t15_PDSvst0_no, c("gene", "number_shRNAs", "logFC", "type", "logCPM", "FDR"))
topids_gene_t15_PDSvst0_no <- as.character(top_gene_t15_PDSvst0_no$gene)
length(topids_gene_t15_PDSvst0_no) # 571

top_gene_t15_PhenDC3vst0_no <- data.table(topTags(lrt_gene_t15_PhenDC3vst0_no, n=Inf)$table)[FDR < 1e-2 & (logFC < -1 | logFC > 1)][,.(gene, number_shRNAs, logFC, logCPM, FDR)]
top_gene_t15_PhenDC3vst0_no[, type := ifelse(logFC > 0, "resistant", ifelse(logFC < 0, "sensitiser", "none"))]
setcolorder(top_gene_t15_PhenDC3vst0_no, c("gene", "number_shRNAs", "logFC", "type", "logCPM", "FDR"))
topids_gene_t15_PhenDC3vst0_no <- as.character(top_gene_t15_PhenDC3vst0_no$gene)
length(topids_gene_t15_PhenDC3vst0_no) # 900


# Plot logFC versus logCPM
plotSmear(lrt_gene_t15_PDSvst0_no, de.tags=topids_gene_t15_PDSvst0_no, main = "t15_PDS vs t0_no")
z_gene_t15_PDSvst0_no_cpm <- data.table(cpm(z_gene)[, grepl("t15_PDS|t0_no", colnames(cpm(z_gene)))], keep.rownames=TRUE)
z_gene_t15_PDSvst0_no_cpm[, log2cpmt0 := (log2(t0_no_1) + log2(t0_no_2) + log2(t0_no_3))/3]
z_gene_t15_PDSvst0_no_cpm[, log2FC := lrt_gene_t15_PDSvst0_no$table$logFC]
z_gene_t15_PDSvst0_no_cpm[, hits := ifelse(rn %in% topids_gene_t15_PDSvst0_no, ifelse(rn %in% topids_gene_t15_DMSOvst0_no, "PDS and DMSO", "PDS only"), "not differentially abundant")]
z_gene_t15_PDSvst0_no_cpm[, order := ifelse(hits == "not differentially abundant", 1, 2)]

gg <- ggplot(z_gene_t15_PDSvst0_no_cpm[log2cpmt0 != -Inf][order(order)], aes(x = log2cpmt0, y = log2FC, colour = hits)) +
geom_point(size=0.5) +
xlab(expression("log"[2]*"(t0)")) +
ylab(expression("log"[2]*"(t15/t0)")) +
theme_bw() +
scale_colour_manual(name="",values = c("grey", "deepskyblue3", "red3"))
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161027_log2t0_log2FC_PDS_DMSO_FDR10minus2_logFC1minus1.pdf", width = 16/2.54, height = 10/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161027_log2t0_log2FC_PDS_DMSO_FDR10minus2_logFC1minus1.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161027/figures/")

plotSmear(lrt_gene_t15_PhenDC3vst0_no, de.tags=topids_gene_t15_PhenDC3vst0_no, main = "t15_PhenDC3 vs t0_no")
z_gene_t15_PhenDC3vst0_no_cpm <- data.table(cpm(z_gene)[, grepl("t15_PhenDC3|t0_no", colnames(cpm(z_gene)))], keep.rownames=TRUE)
z_gene_t15_PhenDC3vst0_no_cpm[, log2cpmt0 := (log2(t0_no_1) + log2(t0_no_2) + log2(t0_no_3))/3]
z_gene_t15_PhenDC3vst0_no_cpm[, log2FC := lrt_gene_t15_PhenDC3vst0_no$table$logFC]
z_gene_t15_PhenDC3vst0_no_cpm[, hits := ifelse(rn %in% topids_gene_t15_PhenDC3vst0_no, ifelse(rn %in% topids_gene_t15_DMSOvst0_no, "PhenDC3 and DMSO", "PhenDC3 only"), "not differentially abundant")]
z_gene_t15_PhenDC3vst0_no_cpm[, order := ifelse(hits == "not differentially abundant", 1, 2)]

gg <- ggplot(z_gene_t15_PhenDC3vst0_no_cpm[log2cpmt0 != -Inf][order(order)], aes(x = log2cpmt0, y = log2FC, colour = hits)) +
geom_point(size=0.5) +
xlab(expression("log"[2]*"(t0)")) +
ylab(expression("log"[2]*"(t15/t0)")) +
theme_bw() +
scale_colour_manual(name="",values = c("grey", "deepskyblue3", "forestgreen"))
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161027_log2t0_log2FC_PhenDC3_DMSO_FDR10minus2_logFC1minus1.pdf", width = 16/2.54, height = 10/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161027_log2t0_log2FC_PhenDC3_DMSO_FDR10minus2_logFC1minus1.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161027/figures/")


# Intersect lists of topTags
length(intersect(topids_gene_t15_DMSOvst0_no, topids_gene_t15_PDSvst0_no)) # 178
length(intersect(topids_gene_t15_DMSOvst0_no, topids_gene_t15_PhenDC3vst0_no)) # 110
length(intersect(topids_gene_t15_PDSvst0_no, topids_gene_t15_PhenDC3vst0_no)) # 152

top_gene_t15_PDSvst0_no_PDSonly <- top_gene_t15_PDSvst0_no[!(topids_gene_t15_PDSvst0_no %in% c(topids_gene_t15_DMSOvst0_no, topids_gene_t15_PhenDC3vst0_no)),]
write.table(top_gene_t15_PDSvst0_no_PDSonly, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161027_top_t15_PDSvst0_no_PDSonly_FDR10minus2_logFC1minus1.txt", quote = FALSE, sep = "\t", row.names = FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161027_top_t15_PDSvst0_no_PDSonly_FDR10minus2_logFC1minus1.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161027/tables/")

top_gene_t15_PhenDC3vst0_no_PhenDC3only <- top_gene_t15_PhenDC3vst0_no[!(topids_gene_t15_PhenDC3vst0_no %in% c(topids_gene_t15_DMSOvst0_no, topids_gene_t15_PDSvst0_no)),]
write.table(top_gene_t15_PhenDC3vst0_no_PhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161027_top_t15_PhenDC3vst0_no_PhenDC3only_FDR10minus2_logFC1minus1.txt", quote = FALSE, sep = "\t", row.names = FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161027_top_t15_PhenDC3vst0_no_PhenDC3only_FDR10minus2_logFC1minus1.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161027/tables/")

top_gene_t15_PDSvst0_no_PDSandPhenDC3 <- top_gene_t15_PDSvst0_no[(topids_gene_t15_PDSvst0_no %in% topids_gene_t15_PhenDC3vst0_no),]
top_gene_t15_PDSvst0_no_PDSandPhenDC3only <- top_gene_t15_PDSvst0_no_PDSandPhenDC3[!(gene %in% topids_gene_t15_DMSOvst0_no),]
write.table(top_gene_t15_PDSvst0_no_PDSandPhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161027_top_t15_PDSvst0_no_PDSandPhenDC3only_FDR10minus2_logFC1minus1.txt", quote = FALSE, sep = "\t", row.names = FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161027_top_t15_PDSvst0_no_PDSandPhenDC3only_FDR10minus2_logFC1minus1.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161027/tables/")


# Venn diagram
venn.plot <- draw.triple.venn(
  area1 = length(topids_gene_t15_PDSvst0_no),
  area2 = length(topids_gene_t15_PhenDC3vst0_no),
  area3 = length(topids_gene_t15_DMSOvst0_no),
  n12 = length(intersect(topids_gene_t15_PDSvst0_no, topids_gene_t15_PhenDC3vst0_no)),
  n23 = length(intersect(topids_gene_t15_PhenDC3vst0_no, topids_gene_t15_DMSOvst0_no)),
  n13 = length(intersect(topids_gene_t15_PDSvst0_no, topids_gene_t15_DMSOvst0_no)),
  n123 = length(intersect(intersect(topids_gene_t15_PDSvst0_no, topids_gene_t15_PhenDC3vst0_no), topids_gene_t15_DMSOvst0_no)),
  category = c(sprintf("PDS (%i)", length(topids_gene_t15_PDSvst0_no)), sprintf("PhenDC3 (%i)", length(topids_gene_t15_PhenDC3vst0_no)), sprintf("DMSO (%i)", length(topids_gene_t15_DMSOvst0_no))),
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


pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161027_t15vst0_venn_FDR10minus2_logFC1minus1.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161027_t15vst0_venn_FDR10minus2_logFC1minus1.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161027/figures/")


```



# (G)

```R
library(edgeR)
library(data.table)
library(ggplot2)
library(VennDiagram)
library(reshape)


# Enlarge the view width when printing tables
options(width = 300)


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

# 100*(110576/113002) = 97.8 shRNAs detected


# Change to a DGE object
groups <- c("t0_no", "t0_no", "t0_no", "t15_DMSO", "t15_DMSO", "t15_DMSO", "t15_PDS", "t15_PDS", "t15_PDS", "t15_PhenDC3", "t15_PhenDC3", "t15_PhenDC3")
genes_dataframe <- data.frame(shRNA = rownames(data), gene = as.character(sapply(rownames(data), function(x) unlist(strsplit(x, "__"))[1])))
z <- DGEList(counts = data, group = groups, genes = genes_dataframe)
z$genes$pool <- sapply(rownames(data), function(x) gsub("Pool", "", tail(unlist(strsplit(x, "_")), n=1)))


# Vector containing all pools
pools <- unique(z$genes$pool)


# e.g. Pool1. This will turn into a for loop on the twelve pools
p <- "1"
sprintf("Pool%s", p)


# Select pool
z_pool <- z[z$genes$pool == p,]
sprintf("shRNAs: %s", nrow(z_pool$counts))


# Redefine lib.size
z_pool$samples$lib.size <- as.numeric(colSums(z_pool$counts))


# Filtering
keep <- rowSums(cpm(z_pool)>0.5) >= 3
z_pool <- z_pool[keep, , keep.lib.sizes=FALSE]
sprintf("shRNAs (after filtering): %s", nrow(z_pool$counts))


# Normalisation
z_pool <- calcNormFactors(z_pool)


# Differential representation analysis.
# Set up design matrix for GLM.
des <- model.matrix(~ 0 + group, data = z_pool$samples)
colnames(des) <- levels(factor(z_pool$samples$group))


# Estimate dispersions
z_pool_glm <- estimateDisp(z_pool, des)
sprintf("Common dispersion: %s", sqrt(z_pool_glm$common.disp))


# Fit negative binomial GLM
z_pool_fit <- glmFit(z_pool_glm, des)


# Define matrix of contrasts
my.contrasts <- makeContrasts(
t15_DMSOvst0_no = t15_DMSO - t0_no,
t15_PDSvst0_no = t15_PDS - t0_no,
t15_PhenDC3vst0_no = t15_PhenDC3 - t0_no,
levels=des)


# Comparisons t15 vs t0 and within t15
# Carry out Likelihood ratio tests
lrt_pool_t15_DMSOvst0_no <- glmLRT(z_pool_fit, contrast=my.contrasts[,"t15_DMSOvst0_no"])
lrt_pool_t15_PDSvst0_no <- glmLRT(z_pool_fit, contrast=my.contrasts[,"t15_PDSvst0_no"])
lrt_pool_t15_PhenDC3vst0_no <- glmLRT(z_pool_fit, contrast=my.contrasts[,"t15_PhenDC3vst0_no"])


# e.g. FDR < 1e-3. This will turn into a for loop on 1e-4, 1e-3, 1e-2
t <- 1e-3

##########
# shRNAs #
##########

# Select
top_t15_DMSOvst0_no <- data.table(topTags(lrt_pool_t15_DMSOvst0_no, n=Inf)$table)[FDR < t & (logFC < -1 | logFC > 1)]
top_t15_DMSOvst0_no[, n_tot_shRNA := sapply(as.character(top_t15_DMSOvst0_no$gene), function(x) sum(z$genes$gene == x))]
topids_t15_DMSOvst0_no <- as.character(top_t15_DMSOvst0_no[n_sig_shRNA/n_tot_shRNA >= 0.5]$gene)
length(topids_t15_DMSOvst0_no) # 3

top_t15_PDSvst0_no <- data.table(topTags(lrt_t15_PDSvst0_no, n=Inf)$table)[FDR < 1e-3 & (logFC < -1 | logFC > 1), .(n_sig_shRNA = uniqueN(shRNA)), by = gene]
top_t15_PDSvst0_no[, n_tot_shRNA := sapply(as.character(top_t15_PDSvst0_no$gene), function(x) sum(z$genes$gene == x))]
topids_t15_PDSvst0_no <- as.character(top_t15_PDSvst0_no[n_sig_shRNA/n_tot_shRNA >= 0.5]$gene)
length(topids_t15_PDSvst0_no) # 9

top_t15_PhenDC3vst0_no <- data.table(topTags(lrt_t15_PhenDC3vst0_no, n=Inf)$table)[FDR < 1e-3 & (logFC < -1 | logFC > 1), .(n_sig_shRNA = uniqueN(shRNA)), by = gene]
top_t15_PhenDC3vst0_no[, n_tot_shRNA := sapply(as.character(top_t15_PhenDC3vst0_no$gene), function(x) sum(z$genes$gene == x))]
topids_t15_PhenDC3vst0_no <- as.character(top_t15_PhenDC3vst0_no[n_sig_shRNA/n_tot_shRNA >= 0.5]$gene)
length(topids_t15_PhenDC3vst0_no) # 41


# Write
top_t15_PDSvst0_no_PDSonly <- top_t15_PDSvst0_no[1:length(topids_t15_PDSvst0_no),][!(topids_t15_PDSvst0_no %in% c(topids_t15_DMSOvst0_no, topids_t15_PhenDC3vst0_no)),]
write.table(top_t15_PDSvst0_no_PDSonly, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_top_t15_PDSvst0_no_PDSonly.txt", quote = FALSE, sep = "\t", row.names = FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_top_t15_PDSvst0_no_PDSonly.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/tables/")

top_t15_PhenDC3vst0_no_PhenDC3only <- top_t15_PhenDC3vst0_no[1:length(topids_t15_PhenDC3vst0_no),][!(topids_t15_PhenDC3vst0_no %in% c(topids_t15_DMSOvst0_no, topids_t15_PDSvst0_no)),]
write.table(top_t15_PhenDC3vst0_no_PhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_top_t15_PhenDC3vst0_no_PhenDC3only.txt", quote = FALSE, sep = "\t", row.names = FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_top_t15_PhenDC3vst0_no_PhenDC3only.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/tables/")

top_t15_PDSvst0_no_PDSandPhenDC3 <- top_t15_PDSvst0_no[1:length(topids_t15_PDSvst0_no),][(topids_t15_PDSvst0_no %in% topids_t15_PhenDC3vst0_no),]
top_t15_PDSvst0_no_PDSandPhenDC3only <- top_t15_PDSvst0_no_PDSandPhenDC3[!(rownames(top_t15_PDSvst0_no_PDSandPhenDC3) %in% topids_t15_DMSOvst0_no),]
write.table(top_t15_PDSvst0_no_PDSandPhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_top_t15_PDSvst0_no_PDSandPhenDC3only.txt", quote = FALSE, sep = "\t")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161019_top_t15_PDSvst0_no_PDSandPhenDC3only.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161019/tables/")


# Venn


############################################
# Select genes from shRNAs with thresholds #
############################################

# Select genes
top_genes_shRNAs_thresholds_t15_DMSOvst0_no <- data.table(topTags(lrt_pool_t15_DMSOvst0_no, n=Inf)$table)[FDR < 1e-3 & (logFC < -1 | logFC > 1), .(n_sig_shRNA = uniqueN(shRNA)), by = gene]
top_genes_shRNAs_thresholds_t15_DMSOvst0_no[, n_tot_shRNA := sapply(as.character(top_genes_shRNAs_thresholds_t15_DMSOvst0_no$gene), function(x) sum(z_pool$genes$gene == x))]
topids_genes_shRNAs_thresholds_t15_DMSOvst0_no <- as.character(top_genes_shRNAs_thresholds_t15_DMSOvst0_no[n_sig_shRNA/n_tot_shRNA >= 0.5]$gene)

top_genes_shRNAs_thresholds_t15_PDSvst0_no <- data.table(topTags(lrt_pool_t15_PDSvst0_no, n=Inf)$table)[FDR < 1e-3 & (logFC < -1 | logFC > 1), .(n_sig_shRNA = uniqueN(shRNA)), by = gene]
top_genes_shRNAs_thresholds_t15_PDSvst0_no[, n_tot_shRNA := sapply(as.character(top_genes_shRNAs_thresholds_t15_PDSvst0_no$gene), function(x) sum(z_pool$genes$gene == x))]
topids_genes_shRNAs_thresholds_t15_PDSvst0_no <- as.character(top_genes_shRNAs_thresholds_t15_PDSvst0_no[n_sig_shRNA/n_tot_shRNA >= 0.5]$gene)
length(topids_t15_PDSvst0_no) # 9

top_t15_PhenDC3vst0_no <- data.table(topTags(lrt_t15_PhenDC3vst0_no, n=Inf)$table)[FDR < 1e-3 & (logFC < -1 | logFC > 1), .(n_sig_shRNA = uniqueN(shRNA)), by = gene]
top_t15_PhenDC3vst0_no[, n_tot_shRNA := sapply(as.character(top_t15_PhenDC3vst0_no$gene), function(x) sum(z$genes$gene == x))]
topids_t15_PhenDC3vst0_no <- as.character(top_t15_PhenDC3vst0_no[n_sig_shRNA/n_tot_shRNA >= 0.5]$gene)
length(topids_t15_PhenDC3vst0_no) # 41


# Write tables PDS only, PhenDC3 only and PDS+PhenDC3 only


# Venn diagram


########################################
# Select genes from shRNAs with camera #
########################################


#######################################
# Select genes from the sum of shRNAs #
#######################################



```


This integrative analysis continues on 20161031_integration.md

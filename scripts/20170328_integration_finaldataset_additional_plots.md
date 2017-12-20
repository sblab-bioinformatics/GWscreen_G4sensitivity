
This script is based on and takes parts from [20170120_integration_finaldataset.md](20170120_integration_finaldataset.md) but also adds the necessary functionality to plot median logFC (by gene) versus logCPM(t0) for PDS only, PhenDC3 only and PDS+PhenDC3 only.

It also adds some functionality to make Volcano plots related to the scatterplots plots. However whether it is easy to calculate medians for LogFC or Logt0 of shRNAs linked to the same gene, how do we combine FDRs?

https://www.biostars.org/p/84242/

https://github.com/brentp/combined-pvalues

https://academic.oup.com/bioinformatics/article/28/22/2986/240603/Comb-p-software-for-combining-analyzing-grouping

http://biorxiv.org/content/biorxiv/early/2016/10/21/082321.full.pdf


# Analysis

```R
library(edgeR)
library(data.table)
library(ggplot2)
library(reshape)
library(Biostrings)
library(ggrepel)


# Enlarge the view width when printing tables
options(width = 250)


# Load data
data <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20170125/20170125_counts.txt", header = T)
dim(data)
# 110669     12


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
  levels=des)

  # Comparisons t15 vs t0
  # Carry out Likelihood ratio tests
  lrt_pool_DMSOvst0 <- glmLRT(z_pool_fit, contrast=my.contrasts[,"DMSOvst0"])
  lrt_pool_PDSvst0 <- glmLRT(z_pool_fit, contrast=my.contrasts[,"PDSvst0"])
  lrt_pool_PhenDC3vst0 <- glmLRT(z_pool_fit, contrast=my.contrasts[,"PhenDC3vst0"])

  # Calculate log2avgcpm t0
  pool_log2avgcpmt0 <- data.frame(log2avgcpmt0 = aveLogCPM(cpm(z_pool)[, grepl("t0_no", colnames(cpm(z_pool)))]))
  rownames(pool_log2avgcpmt0) <- rownames(cpm(z_pool))
  pool_log2avgcpmt0 <- data.table(pool_log2avgcpmt0, keep.rownames=TRUE)

  # Get ranked shRNAS, merge with log2avgcpmt0, combine shRNAs into a table for all pools and write FDR and LogFC(t15/t0) tables for DMSOvst0, PDSvst0 and PhenDC3vst0
  # DMSOvst0
  top_pool_DMSOvst0 <- topTags(lrt_pool_DMSOvst0, n=Inf)
  top_pool_DMSOvst0_log2avgcpmt0 <- merge(data.table(top_pool_DMSOvst0$table), pool_log2avgcpmt0, by.x="shRNA", by.y="rn", all=T, sort = FALSE)
  top_DMSOvst0 <- rbind(top_DMSOvst0, top_pool_DMSOvst0_log2avgcpmt0[order(FDR),.(gene, shRNA, pool, logFC, log2avgcpmt0, FDR)])
  # PDSvst0
  top_pool_PDSvst0 <- topTags(lrt_pool_PDSvst0, n=Inf)
  top_pool_PDSvst0_log2avgcpmt0 <- merge(data.table(top_pool_PDSvst0$table), pool_log2avgcpmt0, by.x="shRNA", by.y="rn", all=T, sort = FALSE)
  top_PDSvst0 <- rbind(top_PDSvst0, top_pool_PDSvst0_log2avgcpmt0[order(FDR),.(gene, shRNA, pool, logFC, log2avgcpmt0, FDR)])
  # PhenDC3vst0
  top_pool_PhenDC3vst0 <- topTags(lrt_pool_PhenDC3vst0, n=Inf)
  top_pool_PhenDC3vst0_log2avgcpmt0 <- merge(data.table(top_pool_PhenDC3vst0$table), pool_log2avgcpmt0, by.x="shRNA", by.y="rn", all=T, sort = FALSE)
  top_PhenDC3vst0 <- rbind(top_PhenDC3vst0, top_pool_PhenDC3vst0_log2avgcpmt0[order(FDR),.(gene, shRNA, pool, logFC, log2avgcpmt0, FDR)])
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

#write.table(all_detected_significant_resistant_sensitiser_genes_DMSOvst0, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_DMSOvst0_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
#system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_DMSOvst0_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/genes_all")

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

#write.table(all_detected_significant_resistant_sensitiser_genes_PDSvst0, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PDSvst0_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
#system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PDSvst0_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/genes_all")


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

#write.table(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PhenDC3vst0_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
#system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20170130_PhenDC3vst0_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170130_integration_all/tables/genes_all")




# Plot median logFC (by gene) versus median logCPM(t0) for PDS only, PhenDC3 only and PDS+PhenDC3 only

#######################
# PDS only sensitiser #
#######################

# Apply 50% or at least 3 shRNAS cutoff and set vectors containing significant genes for all pools for DMSOvst0, PDSvst0 and PhenDC3vst0
# DMSOvst0
topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_DMSOvst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"]$gene))
# PDSvst0
topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSvst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"]$gene))
# PhenDC3vst0
topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PhenDC3vst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"]$gene))

topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSonly <- setdiff(topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSvst0, c(topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_DMSOvst0, topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PhenDC3vst0))

logFC_log2avgcpmt0_genes_PDSvst0 <- top_PDSvst0[, .(logFC_shRNAs_median = median(logFC), log2avgcpmt0_shRNAs_median = median(log2avgcpmt0)), by = gene]
logFC_log2avgcpmt0_genes_PDSvst0[, PDSonly := ifelse(gene %in% topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSonly, "yes", "no")]
logFC_log2avgcpmt0_genes_PDSvst0[, order := ifelse(PDSonly == "no", 1, 2)]

gg <- ggplot(logFC_log2avgcpmt0_genes_PDSvst0[order(order)], aes(x = log2avgcpmt0_shRNAs_median, y = logFC_shRNAs_median, colour = PDSonly)) +
geom_point(size=0.5) +
xlab(expression("median log"[2]*"t0")) +
ylab(expression("median log"[2]*"FC")) +
theme_classic() +
scale_colour_manual(name="",values = c("lightgray", "red3")) +
theme(legend.position="none") +
geom_segment(aes(x = 1, y = 0, xend = 11.75, yend = 0), colour= 'darkgray', linetype= 'solid') +
geom_segment(aes(x = 1, y = -1, xend = 11.75, yend = -1), colour= 'red', linetype= 'dashed') +
coord_cartesian(ylim = c(-7.5, 5))

#geom_text_repel(data = logFC_log2avgcpmt0_genes_PDSvst0[gene %in% c("BRCA1", "BRCA2", "DHX36")], aes(label = gene), color = "black", size = 4, force = 1, fontface = "bold", box.padding = unit(1, "lines"), point.padding = unit(0.25, "lines"), segment.color = "black", nudge_x = 0.5, nudge_y = -0.5)
#https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170404_t15vst0_log2t0_log2FC_significant_sensitiser_genes_PDSonly.pdf", width = 14, height = 14, units = "cm")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170404_t15vst0_log2t0_log2FC_significant_sensitiser_genes_PDSonly.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170404/figures/")

######################
# PDS only resistant #
######################

# Apply 50% or at least 3 shRNAS cutoff and set vectors containing significant genes for all pools for DMSOvst0, PDSvst0 and PhenDC3vst0
# DMSOvst0
topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_DMSOvst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"]$gene))
# PDSvst0
topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PDSvst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"]$gene))
# PhenDC3vst0
topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PhenDC3vst0 <- unique(as.character(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"]$gene))

topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PDSonly <- setdiff(topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PDSvst0, c(topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_DMSOvst0, topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PhenDC3vst0))

logFC_log2avgcpmt0_genes_PDSvst0 <- top_PDSvst0[, .(logFC_shRNAs_median = median(logFC), log2avgcpmt0_shRNAs_median = median(log2avgcpmt0)), by = gene]
logFC_log2avgcpmt0_genes_PDSvst0[, PDSonly := ifelse(gene %in% topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PDSonly, "yes", "no")]
logFC_log2avgcpmt0_genes_PDSvst0[, order := ifelse(PDSonly == "no", 1, 2)]

gg <- ggplot(logFC_log2avgcpmt0_genes_PDSvst0[order(order)], aes(x = log2avgcpmt0_shRNAs_median, y = logFC_shRNAs_median, colour = PDSonly)) +
geom_point(size=0.5) +
xlab(expression("median log"[2]*"t0")) +
ylab(expression("median log"[2]*"FC")) +
theme_classic() +
scale_colour_manual(name="",values = c("lightgray", "red3")) +
theme(legend.position="none") +
geom_segment(aes(x = 1, y = 0, xend = 11.75, yend = 0), colour= 'darkgray', linetype= 'solid') +
geom_segment(aes(x = 1, y = 1, xend = 11.75, yend = 1), colour= 'red', linetype= 'dashed') +
coord_cartesian(ylim = c(-7.5, 5))

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170404_t15vst0_log2t0_log2FC_significant_resistant_genes_PDSonly.pdf", width = 14, height = 14, units = "cm")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170404_t15vst0_log2t0_log2FC_significant_resistant_genes_PDSonly.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170404/figures/")


###########################
# PhenDC3 only sensitiser #
###########################

topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PhenDC3only <- setdiff(topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PhenDC3vst0, c(topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_DMSOvst0, topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSvst0))

logFC_log2avgcpmt0_genes_PhenDC3vst0 <- top_PhenDC3vst0[, .(logFC_shRNAs_median = median(logFC), log2avgcpmt0_shRNAs_median = median(log2avgcpmt0)), by = gene]
logFC_log2avgcpmt0_genes_PhenDC3vst0[, PhenDC3only := ifelse(gene %in% topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PhenDC3only, "yes", "no")]
logFC_log2avgcpmt0_genes_PhenDC3vst0[, order := ifelse(PhenDC3only == "no", 1, 2)]

gg <- ggplot(logFC_log2avgcpmt0_genes_PhenDC3vst0[order(order)], aes(x = log2avgcpmt0_shRNAs_median, y = logFC_shRNAs_median, colour = PhenDC3only)) +
geom_point(size=0.5) +
xlab(expression("median log"[2]*"t0")) +
ylab(expression("median log"[2]*"FC")) +
theme_classic() +
scale_colour_manual(name="",values = c("lightgray", "forestgreen")) +
theme(legend.position="none") +
geom_segment(aes(x = 1, y = 0, xend = 11.75, yend = 0), colour= 'darkgray', linetype= 'solid') +
geom_segment(aes(x = 1, y = -1, xend = 11.75, yend = -1), colour= 'palegreen3', linetype= 'dashed') +
coord_cartesian(ylim = c(-7.5, 5))

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170404_t15vst0_log2t0_log2FC_significant_sensitiser_genes_PhenDC3only.pdf", width = 14, height = 14, units = "cm")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170404_t15vst0_log2t0_log2FC_significant_sensitiser_genes_PhenDC3only.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170404/figures/")

##########################
# PhenDC3 only resistant #
##########################

topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PhenDC3only <- setdiff(topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PhenDC3vst0, c(topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_DMSOvst0, topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PDSvst0))

logFC_log2avgcpmt0_genes_PhenDC3vst0 <- top_PhenDC3vst0[, .(logFC_shRNAs_median = median(logFC), log2avgcpmt0_shRNAs_median = median(log2avgcpmt0)), by = gene]
logFC_log2avgcpmt0_genes_PhenDC3vst0[, PhenDC3only := ifelse(gene %in% topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PhenDC3only, "yes", "no")]
logFC_log2avgcpmt0_genes_PhenDC3vst0[, order := ifelse(PhenDC3only == "no", 1, 2)]

gg <- ggplot(logFC_log2avgcpmt0_genes_PhenDC3vst0[order(order)], aes(x = log2avgcpmt0_shRNAs_median, y = logFC_shRNAs_median, colour = PhenDC3only)) +
geom_point(size=0.5) +
xlab(expression("median log"[2]*"t0")) +
ylab(expression("median log"[2]*"FC")) +
theme_classic() +
scale_colour_manual(name="",values = c("lightgray", "forestgreen")) +
theme(legend.position="none") +
geom_segment(aes(x = 1, y = 0, xend = 11.75, yend = 0), colour= 'darkgray', linetype= 'solid') +
geom_segment(aes(x = 1, y = 1, xend = 11.75, yend = 1), colour= 'palegreen3', linetype= 'dashed') +
coord_cartesian(ylim = c(-7.5, 5))

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170404_t15vst0_log2t0_log2FC_significant_resistant_genes_PhenDC3only.pdf", width = 14, height = 14, units = "cm")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170404_t15vst0_log2t0_log2FC_significant_resistant_genes_PhenDC3only.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170404/figures/")


#################################
# PDSandPhenDC3 only sensitiser #
#################################

topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSandPhenDC3only <- setdiff(intersect(topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PhenDC3vst0, topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSvst0), topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_DMSOvst0)

#PhenDC3coords
logFC_log2avgcpmt0_genes_PhenDC3vst0 <- top_PhenDC3vst0[, .(logFC_shRNAs_median = median(logFC), log2avgcpmt0_shRNAs_median = median(log2avgcpmt0)), by = gene]
logFC_log2avgcpmt0_genes_PhenDC3vst0[, PDSandPhenDC3only := ifelse(gene %in% topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSandPhenDC3only, "yes", "no")]
logFC_log2avgcpmt0_genes_PhenDC3vst0[, order := ifelse(PDSandPhenDC3only == "no", 1, 2)]

gg <- ggplot(logFC_log2avgcpmt0_genes_PhenDC3vst0[order(order)], aes(x = log2avgcpmt0_shRNAs_median, y = logFC_shRNAs_median, colour = PDSandPhenDC3only)) +
geom_point(size=0.5) +
xlab(expression("median log"[2]*"t0")) +
ylab(expression("median log"[2]*"FC")) +
theme_classic() +
scale_colour_manual(name="",values = c("lightgray", "khaki4")) +
theme(legend.position="none") +
geom_segment(aes(x = 1, y = 0, xend = 11.75, yend = 0), colour= 'darkgray', linetype= 'solid') +
geom_segment(aes(x = 1, y = -1, xend = 11.75, yend = -1), colour= 'khaki3', linetype= 'dashed') +
coord_cartesian(ylim = c(-7.5, 5))

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170404_t15vst0_log2t0_log2FC_significant_sensitiser_genes_PDSandPhenDC3only_PhenDC3coords.pdf", width = 14, height = 14, units = "cm")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170404_t15vst0_log2t0_log2FC_significant_sensitiser_genes_PDSandPhenDC3only_PhenDC3coords.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170404/figures/")

#PDScoords
logFC_log2avgcpmt0_genes_PDSvst0 <- top_PDSvst0[, .(logFC_shRNAs_median = median(logFC), log2avgcpmt0_shRNAs_median = median(log2avgcpmt0)), by = gene]
logFC_log2avgcpmt0_genes_PDSvst0[, PDSandPhenDC3only := ifelse(gene %in% topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSandPhenDC3only, "yes", "no")]
logFC_log2avgcpmt0_genes_PDSvst0[, order := ifelse(PDSandPhenDC3only == "no", 1, 2)]

gg <- ggplot(logFC_log2avgcpmt0_genes_PDSvst0[order(order)], aes(x = log2avgcpmt0_shRNAs_median, y = logFC_shRNAs_median, colour = PDSandPhenDC3only)) +
geom_point(size=0.5) +
xlab(expression("median log"[2]*"t0")) +
ylab(expression("median log"[2]*"FC")) +
theme_classic() +
scale_colour_manual(name="",values = c("lightgray", "khaki4")) +
theme(legend.position="none") +
geom_segment(aes(x = 1, y = 0, xend = 11.75, yend = 0), colour= 'darkgray', linetype= 'solid') +
geom_segment(aes(x = 1, y = -1, xend = 11.75, yend = -1), colour= 'khaki3', linetype= 'dashed') +
coord_cartesian(ylim = c(-7.5, 5))

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170404_t15vst0_log2t0_log2FC_significant_sensitiser_genes_PDSandPhenDC3only_PDScoords.pdf", width = 14, height = 14, units = "cm")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170404_t15vst0_log2t0_log2FC_significant_sensitiser_genes_PDSandPhenDC3only_PDScoords.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170404/figures/")


#################################
# PDSandPhenDC3 only sensitiser # with BRCA1 and BRCA2
#################################

topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSandPhenDC3only <- setdiff(intersect(topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PhenDC3vst0, topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSvst0), topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_DMSOvst0)

#PhenDC3coords
logFC_log2avgcpmt0_genes_PhenDC3vst0 <- top_PhenDC3vst0[, .(logFC_shRNAs_median = median(logFC), log2avgcpmt0_shRNAs_median = median(log2avgcpmt0)), by = gene]
logFC_log2avgcpmt0_genes_PhenDC3vst0[, PDSandPhenDC3only := ifelse(gene %in% topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSandPhenDC3only, "yes", "no")]
logFC_log2avgcpmt0_genes_PhenDC3vst0[, order := ifelse(PDSandPhenDC3only == "no", 1, 2)]

gg <- ggplot(logFC_log2avgcpmt0_genes_PhenDC3vst0[order(order)], aes(x = log2avgcpmt0_shRNAs_median, y = logFC_shRNAs_median, colour = PDSandPhenDC3only)) +
geom_point(size=0.5) +
xlab(expression("median log"[2]*"t0")) +
ylab(expression("median log"[2]*"FC")) +
theme_classic() +
scale_colour_manual(name="",values = c("lightgray", "khaki4")) +
theme(legend.position="none") +
geom_segment(aes(x = 1, y = 0, xend = 11.75, yend = 0), colour= 'darkgray', linetype= 'solid') +
geom_segment(aes(x = 1, y = -1, xend = 11.75, yend = -1), colour= 'khaki3', linetype= 'dashed') +
coord_cartesian(ylim = c(-7.5, 5)) +
geom_text_repel(data = logFC_log2avgcpmt0_genes_PhenDC3vst0[gene %in% c("BRCA1", "BRCA2")], aes(label = gene), color = "black", size = 4, force = 1, fontface = "bold", box.padding = unit(1, "lines"), point.padding = unit(0.25, "lines"), segment.color = "black", nudge_x = 0.5, nudge_y = -0.5)

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170404_t15vst0_log2t0_log2FC_significant_sensitiser_genes_PDSandPhenDC3only_PhenDC3coords_BRCA1andBRCA2.pdf", width = 14, height = 14, units = "cm")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170404_t15vst0_log2t0_log2FC_significant_sensitiser_genes_PDSandPhenDC3only_PhenDC3coords_BRCA1andBRCA2.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170404/figures/")

#PDScoords
logFC_log2avgcpmt0_genes_PDSvst0 <- top_PDSvst0[, .(logFC_shRNAs_median = median(logFC), log2avgcpmt0_shRNAs_median = median(log2avgcpmt0)), by = gene]
logFC_log2avgcpmt0_genes_PDSvst0[, PDSandPhenDC3only := ifelse(gene %in% topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSandPhenDC3only, "yes", "no")]
logFC_log2avgcpmt0_genes_PDSvst0[, order := ifelse(PDSandPhenDC3only == "no", 1, 2)]

gg <- ggplot(logFC_log2avgcpmt0_genes_PDSvst0[order(order)], aes(x = log2avgcpmt0_shRNAs_median, y = logFC_shRNAs_median, colour = PDSandPhenDC3only)) +
geom_point(size=0.5) +
xlab(expression("median log"[2]*"t0")) +
ylab(expression("median log"[2]*"FC")) +
theme_classic() +
scale_colour_manual(name="",values = c("lightgray", "khaki4")) +
theme(legend.position="none") +
geom_segment(aes(x = 1, y = 0, xend = 11.75, yend = 0), colour= 'darkgray', linetype= 'solid') +
geom_segment(aes(x = 1, y = -1, xend = 11.75, yend = -1), colour= 'khaki3', linetype= 'dashed') +
coord_cartesian(ylim = c(-7.5, 5)) +
geom_text_repel(data = logFC_log2avgcpmt0_genes_PDSvst0[gene %in% c("BRCA1", "BRCA2")], aes(label = gene), color = "black", size = 4, force = 1, fontface = "bold", box.padding = unit(1, "lines"), point.padding = unit(0.25, "lines"), segment.color = "black", nudge_x = 0.5, nudge_y = -0.5)

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170404_t15vst0_log2t0_log2FC_significant_sensitiser_genes_PDSandPhenDC3only_PDScoords_BRCA1andBRCA2.pdf", width = 14, height = 14, units = "cm")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170404_t15vst0_log2t0_log2FC_significant_sensitiser_genes_PDSandPhenDC3only_PDScoords_BRCA1andBRCA2.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170404/figures/")


################################
# PDSandPhenDC3 only resistant #
################################

topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PDSandPhenDC3only <- setdiff(intersect(topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PhenDC3vst0, topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PDSvst0), topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_DMSOvst0)

#PhenDC3coords
logFC_log2avgcpmt0_genes_PhenDC3vst0 <- top_PhenDC3vst0[, .(logFC_shRNAs_median = median(logFC), log2avgcpmt0_shRNAs_median = median(log2avgcpmt0)), by = gene]
logFC_log2avgcpmt0_genes_PhenDC3vst0[, PDSandPhenDC3only := ifelse(gene %in% topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PDSandPhenDC3only, "yes", "no")]
logFC_log2avgcpmt0_genes_PhenDC3vst0[, order := ifelse(PDSandPhenDC3only == "no", 1, 2)]

gg <- ggplot(logFC_log2avgcpmt0_genes_PhenDC3vst0[order(order)], aes(x = log2avgcpmt0_shRNAs_median, y = logFC_shRNAs_median, colour = PDSandPhenDC3only)) +
geom_point(size=0.5) +
xlab(expression("median log"[2]*"t0")) +
ylab(expression("median log"[2]*"FC")) +
theme_classic() +
scale_colour_manual(name="",values = c("lightgray", "khaki4")) +
theme(legend.position="none") +
geom_segment(aes(x = 1, y = 0, xend = 11.75, yend = 0), colour= 'darkgray', linetype= 'solid') +
geom_segment(aes(x = 1, y = 1, xend = 11.75, yend = 1), colour= 'khaki3', linetype= 'dashed') +
coord_cartesian(ylim = c(-7.5, 5))

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170404_t15vst0_log2t0_log2FC_significant_resistant_genes_PDSandPhenDC3only_PhenDC3coords.pdf", width = 14, height = 14, units = "cm")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170404_t15vst0_log2t0_log2FC_significant_resistant_genes_PDSandPhenDC3only_PhenDC3coords.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170404/figures/")

#PDScoords
logFC_log2avgcpmt0_genes_PDSvst0 <- top_PDSvst0[, .(logFC_shRNAs_median = median(logFC), log2avgcpmt0_shRNAs_median = median(log2avgcpmt0)), by = gene]
logFC_log2avgcpmt0_genes_PDSvst0[, PDSandPhenDC3only := ifelse(gene %in% topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PDSandPhenDC3only, "yes", "no")]
logFC_log2avgcpmt0_genes_PDSvst0[, order := ifelse(PDSandPhenDC3only == "no", 1, 2)]

gg <- ggplot(logFC_log2avgcpmt0_genes_PDSvst0[order(order)], aes(x = log2avgcpmt0_shRNAs_median, y = logFC_shRNAs_median, colour = PDSandPhenDC3only)) +
geom_point(size=0.5) +
xlab(expression("median log"[2]*"t0")) +
ylab(expression("median log"[2]*"FC")) +
theme_classic() +
scale_colour_manual(name="",values = c("lightgray", "khaki4")) +
theme(legend.position="none") +
geom_segment(aes(x = 1, y = 0, xend = 11.75, yend = 0), colour= 'darkgray', linetype= 'solid') +
geom_segment(aes(x = 1, y = 1, xend = 11.75, yend = 1), colour= 'khaki3', linetype= 'dashed') +
coord_cartesian(ylim = c(-7.5, 5))

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170404_t15vst0_log2t0_log2FC_significant_resistant_genes_PDSandPhenDC3only_PDScoords.pdf", width = 14, height = 14, units = "cm")
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20170404_t15vst0_log2t0_log2FC_significant_resistant_genes_PDSandPhenDC3only_PDScoords.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20170404/figures/")

```



# TODO

- Volcano plots in the same way as scatterplots above

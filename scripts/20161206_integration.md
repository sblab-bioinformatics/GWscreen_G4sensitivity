
This script improves the analysis presented in [20161129_integration.md](20161129_integration.md), using the same datasets but adding new features as requested by Darcie and Katie:

* For hairpins, ignore logFC threshold and do it based on an FDR <=0.05 only. For genes, select those with at least 50% or 3 significant shRNAs
* Include median FC and median logFC in final tables
* Add the score from the shERWOOD algorithm to hairpin table
* Add a new table with DMSO genes, in addition to PDS only, PhenDC3 only and PDS+PhenDC3 only.
* Produce barplot of whole genomic screen (shRNA vs LogFC) as sketched in the second A4
* Produce manhattan-like plots for genes and shRNAs like the ones sketched in the first A4. I did not have time to make it in [20161129_integration.md](20161129_integration.md).


# Generate file containing the shRNAid and the shERWOOD algorithm score

```python
with open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3.csv", "rb") as f:
    data = f.readlines()

len(data) # 113003, 113002 hairpins

ofile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3_scores.txt", "w")

for hp in data[1:]:
  fields=hp.split(",")
  i=fields[4]
  score=fields[5]
  if score == "Need score": # for some reason some of the scores in the file are just this "Need score"
    score = "NA"
  pool="".join([fields[8].split()[0], fields[8].split()[1][1:]])
  ofile.write("%s\t%s\n" % ("_".join([i, pool]), score))

ofile.close()

```



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


# Load table of counts as generated in 20161129_integration.md
data <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20161129/20161129_counts.txt", header = T)


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
#[1] "shRNAs: 9524"
#[1] "shRNAs (after filtering): 9419"
#[1] "Common dispersion: 0.394046026915902"
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
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/tables/genes_all")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/tables/shRNAs_all")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/tables/genes_lost")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/venn_diagrams/pass_50pct_or_3_threshold")
system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/tables/genes_pass_50pct_or_3_threshold")




# Tables of genes by DMSOvst0, PDSvst0 and PhenDC3vst0 and for resistant (logFC>0) or sensitiser (logFC<0)

# All genes
f <- readDNAStringSet("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3.fa")
all_shRNAs <- names(f)
all_genes <- data.table(table(as.character(sapply(all_shRNAs, function(x) unlist(strsplit(x, "__"))[1])))) # all genes, 18937


# DMSOvst0
detected_genes_DMSOvst0 <- top_DMSOvst0[,.(n_detected_shRNA = uniqueN(shRNA)), by = gene] # detected genes DMSOvst0, 18908, 29 genes with no shRNAs detected, shRNAs filtered are not included in this table either

all_detected_genes_DMSOvst0 <- merge(all_genes, detected_genes_DMSOvst0, by.x="V1", by.y="gene", all=T)
setnames(all_detected_genes_DMSOvst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs"))

significant_resistant_genes_DMSOvst0 <- top_DMSOvst0[FDR <= 0.05 & logFC > 0,.(n_significant_resistant_shRNAs = uniqueN(shRNA), logFC_significant_resistant_shRNAs_average = mean(logFC), FC_significant_resistant_shRNAs_average = 2^(mean(logFC)), logFC_significant_resistant_shRNAs_median = median(logFC), FC_significant_resistant_shRNAs_median = 2^(median(logFC))), by = gene] # significant and resistant genes DMSOvst0, 4727
significant_sensitiser_genes_DMSOvst0 <- top_DMSOvst0[FDR <= 0.05 & logFC < 0,.(n_significant_sensitiser_shRNAs = uniqueN(shRNA), logFC_significant_sensitiser_shRNAs_average = mean(logFC), FC_significant_sensitiser_shRNAs_average = 2^(mean(logFC)), logFC_significant_sensitiser_shRNAs_median = median(logFC), FC_significant_sensitiser_shRNAs_median = 2^(median(logFC))), by = gene] # significant and sensitiser genes DMSOvst0, 4345
length(intersect(significant_resistant_genes_DMSOvst0$gene, significant_sensitiser_genes_DMSOvst0$gene)) # 817

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
#1: KDM3B             18                17                              4                                   1.670458                                3.183157                                  1.448301                               2.728865
#pct_significant_resistant_shRNAs pass_50_threshold_significant_resistant_shRNAs pass_3_significant_resistant_shRNAs n_significant_sensitiser_shRNAs logFC_significant_sensitiser_shRNAs_average FC_significant_sensitiser_shRNAs_average
#1:                         23.52941                                             no                                 yes                               2                                   -2.240795                                0.2115697
#logFC_significant_sensitiser_shRNAs_median FC_significant_sensitiser_shRNAs_median pct_significant_sensitiser_shRNAs pass_50_threshold_significant_sensitiser_shRNAs pass_3_significant_sensitiser_shRNAs
#1:                                  -2.240795                               0.2115697                          11.76471                                              no                                   no

nrow(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"]) # 277 genes are resistant in DMSOvst0, have more than 50% or at least 3 shRNAS detected with FDR <= 0.05 & logFC > 0
nrow(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"]) # 407 genes are sensitiser in DMSOvst0, have more than 50% or at least 3 shRNAS detected with FDR <= 0.05 & logFC < 0
nrow(all_detected_significant_resistant_sensitiser_genes_DMSOvst0[(pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes") & (pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes")]) # 2 genes are both resistant and sensitiser in DMSOvst0, OR8B8 and ZNF434

write.table(all_detected_significant_resistant_sensitiser_genes_DMSOvst0, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161206_DMSOvst0_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161206_DMSOvst0_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/tables/genes_all")


# PDSvst0
detected_genes_PDSvst0 <- top_PDSvst0[,.(n_detected_shRNA = uniqueN(shRNA)), by = gene] # detected genes PDSvst0, 18908, 29 genes with no shRNAs detected, shRNAs filtered are not included in this table either

all_detected_genes_PDSvst0 <- merge(all_genes, detected_genes_PDSvst0, by.x="V1", by.y="gene", all=T)
setnames(all_detected_genes_PDSvst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs"))

significant_resistant_genes_PDSvst0 <- top_PDSvst0[FDR <= 0.05 & logFC > 0,.(n_significant_resistant_shRNAs = uniqueN(shRNA), logFC_significant_resistant_shRNAs_average = mean(logFC), FC_significant_resistant_shRNAs_average = 2^(mean(logFC)), logFC_significant_resistant_shRNAs_median = median(logFC), FC_significant_resistant_shRNAs_median = 2^(median(logFC))), by = gene] # significant and resistant genes PDSvst0, 5790
significant_sensitiser_genes_PDSvst0 <- top_PDSvst0[FDR <= 0.05 & logFC < 0,.(n_significant_sensitiser_shRNAs = uniqueN(shRNA), logFC_significant_sensitiser_shRNAs_average = mean(logFC), FC_significant_sensitiser_shRNAs_average = 2^(mean(logFC)), logFC_significant_sensitiser_shRNAs_median = median(logFC), FC_significant_sensitiser_shRNAs_median = 2^(median(logFC))), by = gene] # significant and sensitiser genes PDSvst0, 5230
length(intersect(significant_resistant_genes_PDSvst0$gene, significant_sensitiser_genes_PDSvst0$gene)) # 1256

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
#1: KDM3B             18                17                              6                                   1.306119                                2.472755                                  1.113397                               2.163545
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

nrow(all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"]) # 453 genes are resistant in PDSvst0, have more than 50% or at least 3 shRNAS detected with FDR <= 0.05 & logFC > 0
nrow(all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"]) # 561 genes are sensitiser in PDSvst0, have more than 50% or at least 3 shRNAS detected with FDR <= 0.05 & logFC < 0
nrow(all_detected_significant_resistant_sensitiser_genes_PDSvst0[(pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes") & (pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes")]) # 5 genes are both resistant and sensitiser in PDSvst0, C1orf88, COL8A2, DCUN1D4, KIAA1370, OR6T1

write.table(all_detected_significant_resistant_sensitiser_genes_PDSvst0, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161206_PDSvst0_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161206_PDSvst0_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/tables/genes_all")


# PhenDC3vst0
detected_genes_PhenDC3vst0 <- top_PhenDC3vst0[,.(n_detected_shRNA = uniqueN(shRNA)), by = gene] # detected genes PhenDC3vst0, 18908, 29 genes with no shRNAs detected, shRNAs filtered are not included in this table either

all_detected_genes_PhenDC3vst0 <- merge(all_genes, detected_genes_PhenDC3vst0, by.x="V1", by.y="gene", all=T)
setnames(all_detected_genes_PhenDC3vst0, c("gene", "n_total_shRNAs", "n_detected_shRNAs"))

significant_resistant_genes_PhenDC3vst0 <- top_PhenDC3vst0[FDR <= 0.05 & logFC > 0,.(n_significant_resistant_shRNAs = uniqueN(shRNA), logFC_significant_resistant_shRNAs_average = mean(logFC), FC_significant_resistant_shRNAs_average = 2^(mean(logFC)), logFC_significant_resistant_shRNAs_median = median(logFC), FC_significant_resistant_shRNAs_median = 2^(median(logFC))), by = gene] # significant and resistant genes PhenDC3vst0, 6625
significant_sensitiser_genes_PhenDC3vst0 <- top_PhenDC3vst0[FDR <= 0.05 & logFC < 0,.(n_significant_sensitiser_shRNAs = uniqueN(shRNA), logFC_significant_sensitiser_shRNAs_average = mean(logFC), FC_significant_sensitiser_shRNAs_average = 2^(mean(logFC)), logFC_significant_sensitiser_shRNAs_median = median(logFC), FC_significant_sensitiser_shRNAs_median = 2^(median(logFC))), by = gene] # significant and sensitiser genes PhenDC3vst0, 6552
length(intersect(significant_resistant_genes_PhenDC3vst0$gene, significant_sensitiser_genes_PhenDC3vst0$gene)) # 2004, e.g. KDM3B

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
#1: KDM3B             18                17                              2                                   2.002774                                  4.0077                                  2.002774                                 4.0077
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

nrow(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"]) # 545 genes are resistant in PhenDC3vst0, have more than 50% or at least 3 shRNAS detected with FDR <= 0.05 & logFC > 0
nrow(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"]) # 769 genes are sensitiser in PhenDC3vst0, have more than 50% or at least 3 shRNAS detected with FDR <= 0.05 & logFC < 0
nrow(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[(pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes") & (pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes")]) # 11 genes are both resistant and sensitiser in PhenDC3vst0 at the same time, C10orf46, HMGB1, LIMCH1, MIER3, OR5I1, OSBPL11, PI3, TSPAN7, TTC14, UTRN, UTS2D

write.table(all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161206_PhenDC3vst0_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161206_PhenDC3vst0_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/tables/genes_all")



# Tables of shRNAs by DMSOvst0, PDSvst0 and PhenDC3vst0

scores <- data.table(read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160420/Hs_WG_pool_data_P003_with_12Nic_3_scores.txt"))
sum(is.na(as.character(scores$V2))) # 31, these were "Need score in the original table"
scores[V2 == "NULL"] # 2984
scores[V2 == "NULL", V2 := "10"] # see 20161206 14:10 email from Darcie
setnames(scores, c("cshlid","score") )

# DMSOvst0
top_DMSOvst0_scores <- merge(top_DMSOvst0, scores, by.x="shRNA", by.y="cshlid", all.x=T)[order(gene, FDR), .(gene, shRNA, score, pool, logFC, FDR)]
write.table(top_DMSOvst0_scores, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_DMSOvst0_shRNAs.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_DMSOvst0_shRNAs.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/tables/shRNAs_all")

# PDSvst0
top_PDSvst0_scores <- merge(top_PDSvst0, scores, by.x="shRNA", by.y="cshlid", all.x=T)[order(gene, FDR), .(gene, shRNA, score, pool, logFC, FDR)]
write.table(top_PDSvst0_scores, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_PDSvst0_shRNAs.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_PDSvst0_shRNAs.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/tables/shRNAs_all")

# PhenDC3vst0
top_PhenDC3vst0_scores <- merge(top_PhenDC3vst0, scores, by.x="shRNA", by.y="cshlid", all.x=T)[order(gene, FDR), .(gene, shRNA, score, pool, logFC, FDR)]
write.table(top_PhenDC3vst0_scores, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_PhenDC3vst0_shRNAs.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_PhenDC3vst0_shRNAs.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/tables/shRNAs_all")



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

write.table(setcolorder(l_combined, c("gene", "lost")), file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_genes_lost_t0_bypool.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_genes_lost_t0_bypool.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/tables/genes_lost")




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

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161207_t15vst0_genes_allpools_venn_FDR5minus2_pass_50pct_or_3_threshold_significant_resistant_shRNAs.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161207_t15vst0_genes_allpools_venn_FDR5minus2_pass_50pct_or_3_threshold_significant_resistant_shRNAs.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/venn_diagrams/pass_50pct_or_3_threshold")



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

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161207_t15vst0_genes_allpools_venn_FDR5minus2_pass_50pct_or_3_threshold_significant_sensitiser_shRNAs.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()

system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161207_t15vst0_genes_allpools_venn_FDR5minus2_pass_50pct_or_3_threshold_significant_sensitiser_shRNAs.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/venn_diagrams/pass_50pct_or_3_threshold")




# Tables of genes for PDS only, PhenDC3 only and PDS+PhenDC3 only with FDR <= 0.05 & logFC < 0 (sensitiser) or logFC > 0 (resistant) & passing the 50% or at least 3 threshold of significant shRNAS

# resistant
# PDS only
top_PDSonly_resistant <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"][!(gene %in% c(topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_DMSOvst0, topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PhenDC3vst0))][order(-logFC_significant_resistant_shRNAs_average),.(gene, n_detected_shRNAs, n_significant_resistant_shRNAs, pct_significant_resistant_shRNAs, logFC_significant_resistant_shRNAs_average, FC_significant_resistant_shRNAs_average, logFC_significant_resistant_shRNAs_median, FC_significant_resistant_shRNAs_median)]
write.table(top_PDSonly_resistant, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_PDSonly_resistant_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_PDSonly_resistant_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/tables/genes_pass_50pct_or_3_threshold")

# PhenDC3 only
top_PhenDC3only_resistant <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"][!(gene %in% c(topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_DMSOvst0, topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PDSvst0))][order(-logFC_significant_resistant_shRNAs_average),.(gene, n_detected_shRNAs, n_significant_resistant_shRNAs, pct_significant_resistant_shRNAs, logFC_significant_resistant_shRNAs_average, FC_significant_resistant_shRNAs_average, logFC_significant_resistant_shRNAs_median, FC_significant_resistant_shRNAs_median)]
write.table(top_PhenDC3only_resistant, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_PhenDC3only_resistant_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_PhenDC3only_resistant_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/tables/genes_pass_50pct_or_3_threshold")

# PDS+PhenDC3 only
top_PDSandPhenDC3only_PhenDC3_resistant <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"][gene %in% topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PDSvst0][!(gene %in% topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_DMSOvst0)][order(-logFC_significant_resistant_shRNAs_average),.(gene, n_detected_shRNAs, n_significant_resistant_shRNAs, pct_significant_resistant_shRNAs, logFC_significant_resistant_shRNAs_average, FC_significant_resistant_shRNAs_average, logFC_significant_resistant_shRNAs_median, FC_significant_resistant_shRNAs_median)]
top_PDSandPhenDC3only_PDS_resistant <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"][gene %in% topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_PhenDC3vst0][!(gene %in% topids_genes_pass_50_3_threshold_significant_resistant_shRNAs_DMSOvst0)][order(-logFC_significant_resistant_shRNAs_average),.(gene, n_detected_shRNAs, n_significant_resistant_shRNAs, pct_significant_resistant_shRNAs, logFC_significant_resistant_shRNAs_average, FC_significant_resistant_shRNAs_average, logFC_significant_resistant_shRNAs_median, FC_significant_resistant_shRNAs_median)]
top_PDSandPhenDC3only_resistant <- merge(top_PDSandPhenDC3only_PhenDC3_resistant, top_PDSandPhenDC3only_PDS_resistant, by = c("gene", "n_detected_shRNAs"), suffixes = c(".PhenDC3", ".PDS"))[order(-logFC_significant_resistant_shRNAs_average.PhenDC3)]
write.table(top_PDSandPhenDC3only_resistant, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_PDSandPhenDC3only_resistant_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_PDSandPhenDC3only_resistant_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/tables/genes_pass_50pct_or_3_threshold")

# DMSO
top_DMSO_resistant <- all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_resistant_shRNAs == "yes" | pass_3_significant_resistant_shRNAs == "yes"][order(-logFC_significant_resistant_shRNAs_average),.(gene, n_detected_shRNAs, n_significant_resistant_shRNAs, pct_significant_resistant_shRNAs, logFC_significant_resistant_shRNAs_average, FC_significant_resistant_shRNAs_average, logFC_significant_resistant_shRNAs_median, FC_significant_resistant_shRNAs_median)]
write.table(top_DMSO_resistant, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_DMSO_resistant_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_DMSO_resistant_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/tables/genes_pass_50pct_or_3_threshold")


# sensitiser
# PDS only
top_PDSonly_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"][!(gene %in% c(topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_DMSOvst0, topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PhenDC3vst0))][order(logFC_significant_sensitiser_shRNAs_average),.(gene, n_detected_shRNAs, n_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs, logFC_significant_sensitiser_shRNAs_average, FC_significant_sensitiser_shRNAs_average, logFC_significant_sensitiser_shRNAs_median, FC_significant_sensitiser_shRNAs_median)]
write.table(top_PDSonly_sensitiser, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_PDSonly_sensitiser_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_PDSonly_sensitiser_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/tables/genes_pass_50pct_or_3_threshold")

# PhenDC3 only
top_PhenDC3only_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"][!(gene %in% c(topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_DMSOvst0, topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSvst0))][order(logFC_significant_sensitiser_shRNAs_average),.(gene, n_detected_shRNAs, n_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs, logFC_significant_sensitiser_shRNAs_average, FC_significant_sensitiser_shRNAs_average, logFC_significant_sensitiser_shRNAs_median, FC_significant_sensitiser_shRNAs_median)]
write.table(top_PhenDC3only_sensitiser, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_PhenDC3only_sensitiser_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_PhenDC3only_sensitiser_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/tables/genes_pass_50pct_or_3_threshold")

# PDS+PhenDC3 only
top_PDSandPhenDC3only_PhenDC3_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PhenDC3vst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"][gene %in% topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PDSvst0][!(gene %in% topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_DMSOvst0)][order(logFC_significant_sensitiser_shRNAs_average),.(gene, n_detected_shRNAs, n_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs, logFC_significant_sensitiser_shRNAs_average, FC_significant_sensitiser_shRNAs_average, logFC_significant_sensitiser_shRNAs_median, FC_significant_sensitiser_shRNAs_median)]
top_PDSandPhenDC3only_PDS_sensitiser <- all_detected_significant_resistant_sensitiser_genes_PDSvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"][gene %in% topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_PhenDC3vst0][!(gene %in% topids_genes_pass_50_3_threshold_significant_sensitiser_shRNAs_DMSOvst0)][order(logFC_significant_sensitiser_shRNAs_average),.(gene, n_detected_shRNAs, n_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs, logFC_significant_sensitiser_shRNAs_average, FC_significant_sensitiser_shRNAs_average, logFC_significant_sensitiser_shRNAs_median, FC_significant_sensitiser_shRNAs_median)]
top_PDSandPhenDC3only_sensitiser <- merge(top_PDSandPhenDC3only_PhenDC3_sensitiser, top_PDSandPhenDC3only_PDS_sensitiser, by = c("gene", "n_detected_shRNAs"), suffixes = c(".PhenDC3", ".PDS"))[order(logFC_significant_sensitiser_shRNAs_average.PhenDC3)]
write.table(top_PDSandPhenDC3only_sensitiser, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_PDSandPhenDC3only_sensitiser_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_PDSandPhenDC3only_sensitiser_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/tables/genes_pass_50pct_or_3_threshold")

# DMSO
top_DMSO_sensitiser <- all_detected_significant_resistant_sensitiser_genes_DMSOvst0[pass_50_threshold_significant_sensitiser_shRNAs == "yes" | pass_3_significant_sensitiser_shRNAs == "yes"][order(logFC_significant_sensitiser_shRNAs_average),.(gene, n_detected_shRNAs, n_significant_sensitiser_shRNAs, pct_significant_sensitiser_shRNAs, logFC_significant_sensitiser_shRNAs_average, FC_significant_sensitiser_shRNAs_average, logFC_significant_sensitiser_shRNAs_median, FC_significant_sensitiser_shRNAs_median)]
write.table(top_DMSO_sensitiser, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_DMSO_sensitiser_genes.txt", quote = FALSE, sep = "\t", row.names=FALSE)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20161207_DMSO_sensitiser_genes.txt martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/tables/genes_pass_50pct_or_3_threshold")




# Distribution of scores from TRANSOMIC

system("ssh martin03@nm149s012925 mkdir -p /media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/shERWOOD_algorithm_score_distribution")

gg <- ggplot(scores, aes(x = as.numeric(as.vector(SCORE)))) +
stat_bin() +
xlab("shERWOOD algorithm score") +
ylab("Frequency") +
theme_classic() +
coord_cartesian(xlim = c(0, 10))

ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161206_shERWOOD_score.pdf", width = 14/2.54, height = 10/2.54)
system("rsync --remove-source-files /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20161206_shERWOOD_score.pdf martin03@nm149s012925:/media/itshare/research/sblab/group/public_folders/martin03/20160418_shRNAscreen_katie_darcie/20161206/shERWOOD_algorithm_score_distribution")





```

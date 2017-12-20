


This is an analysis of the new data batch obtained by Katie and Darcie.

This script is about the A375 cell line (Pool 8).

Where are the fastq files? Copy them to lustre:


In sblab-srv001:

```bash
cd /media/staging/160426_NS500222_0160_HYLFFBGXX/fastq
rsync --progress *.fq.gz martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/fastq
rsync --progress SLX-11622.HYLFFBGXX.s_1.barcodesummary.txt martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/fastq
rsync --progress SLX-11622.HYLFFBGXX.s_1.contents.csv martin03@uk-cri-lcst01:/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/fastq
```


Nic processed the fastq files, trimmed them, fastqc them and aligned them using bowtie. He produced a table of counts, which I copied here:

/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/all_counts.txt


We started from the table of counts:


```R
library(edgeR)


# Load table of counts into edgeR
data <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/all_counts.txt", header = T)
dim(data) # 44763    22


# Add first column to rownames
rownames(data) <- data[,1]
data <- data[2:22]
dim(data)


# Change to a DGE object
z <- DGEList(counts=data)
sum(grepl("BRCA2", rownames(z$counts))) # 3, this is the number of hairpins targeting BRCA2
sum(grepl("Pool-8", rownames(z$counts))) # 9583, this is the number of pool8 hairpins


# Change column labels z$counts, change row labels z$samples, add group to z$samples, add Replicate to z$samples
treatments <- rep(c(rep("DMSO", 2), rep("PDS", 2), rep("PhenDC3", 2), "no"),3)
timepoints <- rep(c(rep(rep(c("t14", "t7")),3), "t0"),3)
replicates <- c(rep("1", 7), rep("2", 7), rep("3", 7))

colnames(z$counts) <- paste(timepoints, treatments, replicates, sep="_")
rownames(z$samples) <- paste(timepoints, treatments, replicates, sep="_")
z$samples$group <- paste(timepoints, treatments, sep="_")
z$samples$Replicate <- replicates


# Select pool8
z8 <- z[grepl("Pool-8", rownames(z$counts)),]
sum(grepl("BRCA2", rownames(z8$counts))) # 2 - OK


# Write raw counts to table (only pool8)
write.table(z8$counts, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160526_A375_pool8_counts.txt", quote = FALSE, sep = "\t")


# Plot number of hairpins that could be matched per sample
# and total for each hairpin across all samples
par(mfrow=c(2,1))
barplot(colSums(z8$counts), las=2, main="Counts per index", cex.names=0.5, cex.axis=0.8)
barplot(rowSums(z8$counts), las=2, main="Counts per hairpin", cex.names=0.5, cex.axis=0.8, axisnames = FALSE)

sort(rowSums(z8$counts)[rowSums(z8$counts) > 300000], decreasing=T)
#CYP2C9__4__3-Pool-8   MMP8__2__2-Pool-8
#             328860              327617


# Make MDS plots to visualise relationships between replicate samples
plotMDS(z8, labels = z8$samples$group, col = rep(c(1,1,2,2,3,3,4), 3)) # coloured by treatment
legend("bottomright", legend=c("DMSO", "PDS", "PhenDC3", "t0"), col=1:4, pch=15)

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160510_A375_pool8_mds.pdf", width = 9, height = 9)
plotMDS(z8, labels = rownames(z8$samples), col = rep(c(1,1,2,2,3,3,4), 3)) # coloured by treatment
legend("bottomright", legend=c("DMSO", "PDS", "PhenDC3", "t0"), col=1:4, pch=15)
dev.off()

z8_t14vst0 <- z8[,grepl("t14|t0", colnames(z8$counts))]
pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160517_A375_pool8_t14vst0_mds.pdf", width = 9, height = 9)
plotMDS(z8_t14vst0, labels = rownames(z8_t14vst0$samples), col = rep(c(1,2,3,4), 3)) # coloured by treatment
legend("bottomleft", legend=c("DMSO", "PDS", "PhenDC3", "t0"), col=1:4, pch=15)
dev.off()

z8_t7vst0 <- z8[,grepl("t7|t0", colnames(z8$counts))]
pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160517_A375_pool8_t7vst0_mds.pdf", width = 9, height = 9)
plotMDS(z8_t7vst0, labels = rownames(z8_t7vst0$samples), col = rep(c(1,2,3,4), 3)) # coloured by treatment
legend("bottomleft", legend=c("DMSO", "PDS", "PhenDC3", "t0"), col=1:4, pch=15)
dev.off()

plotMDS(z8, labels = z8$samples$group, col = rep(c(1,2,1,2,1,2,3), 3)) # coloured by timepoint
legend("bottomright", legend=c("t14", "t7", "t0"), col=1:3, pch=15)

z8_t0_no_t14_PhenDC3 <- z8[, grepl("t14_PhenDC3|t0_no", colnames(z8$counts))]
plotMDS(z8_t0_no_t14_PhenDC3, labels = rownames(z8_t0_no_t14_PhenDC3$samples), col = rep(c(1,2), 3))

z8_t0_no_t7_PDS <- z8[, grepl("t7_PDS|t0_no", colnames(z8$counts))]
plotMDS(z8_t0_no_t7_PDS, labels = rownames(z8_t0_no_t7_PDS$samples), col = rep(c(1,2), 3))


# Begin differential representation analysis.
# We will use GLMs in edgeR in this case since there are more than 2 groups.
# Set up design matrix for GLM.
des <- model.matrix(~ 0 + group, data = z8$samples) # the 0 + in the model formula is an instruction not to include an intercept column and instead to include a column for each group.
colnames(des) <- levels(factor(z8$samples$group))
des


# Estimate dispersions
z8glm <- estimateDisp(z8, des)
sqrt(z8glm$common.disp) # 0.4124458


# Plot BCVs versus abundance
plotBCV(z8glm)


# Fit negative bionomial GLM
z8fit <- glmFit(z8glm, des)


# Define matrix of contrasts
# edgeR user's guide: https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
my.contrasts <- makeContrasts(
t7_DMSOvst0_no = t7_DMSO - t0_no,
t7_PDSvst0_no = t7_PDS - t0_no,
t7_PhenDC3vst0_no = t7_PhenDC3 - t0_no,
t7_PDSvst7_DMSO = t7_PDS - t7_DMSO,
t7_PhenDC3vst7_DMSO = t7_PhenDC3 - t7_DMSO,
t7_PhenDC3vst7_PDS = t7_PhenDC3 - t7_PDS,
t14_DMSOvst0_no = t14_DMSO - t0_no,
t14_PDSvst0_no = t14_PDS - t0_no,
t14_PhenDC3vst0_no = t14_PhenDC3 - t0_no,
t14_PDSvst14_DMSO = t14_PDS - t14_DMSO,
t14_PhenDC3vst14_DMSO = t14_PhenDC3 - t14_DMSO,
t14_PhenDC3vst14_PDS = t14_PhenDC3 - t14_PDS,
t14_DMSOvst7_DMSO = t14_DMSO - t7_DMSO,
t14_PDSvst7_PDS = t14_PDS - t7_PDS,
t14_PhenDC3vst7_PhenDC3 = t14_PhenDC3 - t7_PhenDC3,
levels=des)


# Comparisons t14 vs t0 and within t14
# Carry out Likelihood ratio tests
lrt_t14_DMSOvst0_no <- glmLRT(z8fit, contrast=my.contrasts[,"t14_DMSOvst0_no"])
lrt_t14_PDSvst0_no <- glmLRT(z8fit, contrast=my.contrasts[,"t14_PDSvst0_no"])
lrt_t14_PhenDC3vst0_no <- glmLRT(z8fit, contrast=my.contrasts[,"t14_PhenDC3vst0_no"])
lrt_t14_PDSvst14_DMSO <- glmLRT(z8fit, contrast=my.contrasts[,"t14_PDSvst14_DMSO"])
lrt_t14_PhenDC3vst14_DMSO <- glmLRT(z8fit, contrast=my.contrasts[,"t14_PhenDC3vst14_DMSO"])
lrt_t14_PhenDC3vst14_PDS <- glmLRT(z8fit, contrast=my.contrasts[,"t14_PhenDC3vst14_PDS"])


# Show top ranked hairpins
topTags(lrt_t14_DMSOvst0_no)
topTags(lrt_t14_PDSvst0_no)
topTags(lrt_t14_PhenDC3vst0_no)
topTags(lrt_t14_PDSvst14_DMSO)
topTags(lrt_t14_PhenDC3vst14_DMSO)
topTags(lrt_t14_PhenDC3vst14_PDS)


# Select hairpins with FDR < 0.05 to highlight on plot
thresh <- 0.05

top_t14_DMSOvst0_no <- topTags(lrt_t14_DMSOvst0_no, n=Inf)
topids_t14_DMSOvst0_no <- rownames(top_t14_DMSOvst0_no$table[top_t14_DMSOvst0_no$table$FDR < thresh,])
length(topids_t14_DMSOvst0_no) # 922

top_t14_PDSvst0_no <- topTags(lrt_t14_PDSvst0_no, n=Inf)
topids_t14_PDSvst0_no <- rownames(top_t14_PDSvst0_no$table[top_t14_PDSvst0_no$table$FDR < thresh,])
length(topids_t14_PDSvst0_no) # 1091

top_t14_PhenDC3vst0_no <- topTags(lrt_t14_PhenDC3vst0_no, n=Inf)
topids_t14_PhenDC3vst0_no <- rownames(top_t14_PhenDC3vst0_no$table[top_t14_PhenDC3vst0_no$table$FDR < thresh,])
length(topids_t14_PhenDC3vst0_no) # 350

top_t14_PDSvst14_DMSO <- topTags(lrt_t14_PDSvst14_DMSO, n=Inf)
topids_t14_PDSvst14_DMSO <- rownames(top_t14_PDSvst14_DMSO$table[top_t14_PDSvst14_DMSO$table$FDR < thresh,])
length(topids_t14_PDSvst14_DMSO) # 696

top_t14_PhenDC3vst14_DMSO <- topTags(lrt_t14_PhenDC3vst14_DMSO, n=Inf)
topids_t14_PhenDC3vst14_DMSO <- rownames(top_t14_PhenDC3vst14_DMSO$table[top_t14_PhenDC3vst14_DMSO$table$FDR < thresh,])
length(topids_t14_PhenDC3vst14_DMSO) # 26

top_t14_PhenDC3vst14_PDS <- topTags(lrt_t14_PhenDC3vst14_PDS, n=Inf)
topids_t14_PhenDC3vst14_PDS <- rownames(top_t14_PhenDC3vst14_PDS$table[top_t14_PhenDC3vst14_PDS$table$FDR < thresh,])
length(topids_t14_PhenDC3vst14_PDS) # 622


# Write LogFC and FDR tables (only pool8)
write.table(data.table(as.data.frame(top_t14_DMSOvst0_no), keep.rownames=TRUE)[order(logFC),.(rn, logFC, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160527_A375_pool8_t14_DMSOvst0_logFCincreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
write.table(data.table(as.data.frame(top_t14_DMSOvst0_no), keep.rownames=TRUE)[order(-logFC),.(rn, logFC, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160527_A375_pool8_t14_DMSOvst0_logFCdecreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
write.table(data.table(as.data.frame(top_t14_PDSvst0_no), keep.rownames=TRUE)[order(logFC),.(rn, logFC, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160527_A375_pool8_t14_PDSvst0_logFCincreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
write.table(data.table(as.data.frame(top_t14_PDSvst0_no), keep.rownames=TRUE)[order(-logFC),.(rn, logFC, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160527_A375_pool8_t14_PDSvst0_logFCdecreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
write.table(data.table(as.data.frame(top_t14_PhenDC3vst0_no), keep.rownames=TRUE)[order(logFC),.(rn, logFC, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160527_A375_pool8_t14_PhenDC3vst0_logFCincreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
write.table(data.table(as.data.frame(top_t14_PhenDC3vst0_no), keep.rownames=TRUE)[order(-logFC),.(rn, logFC, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160527_A375_pool8_t14_PhenDC3vst0_logFCdecreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)


# Plot logFC versus logCPM
plotSmear(lrt_t14_DMSOvst0_no, de.tags=topids_t14_DMSOvst0_no, main = "t14_DMSO vs t0_no")
plotSmear(lrt_t14_PDSvst0_no, de.tags=topids_t14_PDSvst0_no, main = "t14_PDS vs t0_no")
plotSmear(lrt_t14_PhenDC3vst0_no, de.tags=topids_t14_PhenDC3vst0_no, main = "t14_PhenDC3 vs t0_no")
plotSmear(lrt_t14_PDSvst14_DMSO, de.tags=topids_t14_PDSvst14_DMSO, main = "t14_PDS vs t14_DMSO")
plotSmear(lrt_t14_PhenDC3vst14_DMSO, de.tags=topids_t14_PhenDC3vst14_DMSO, main = "t14_PhenDC3 vs t14_DMSO")
plotSmear(lrt_t14_PhenDC3vst14_PDS, de.tags=topids_t14_PhenDC3vst14_PDS, main = "t14_PhenDC3 vs t14_PDS")


# Intersect lists of topTags
# t14_DMSOvst0_no, t14_PDSvst0_no, t14_PhenDC3vst0_no
length(intersect(topids_t14_DMSOvst0_no, topids_t14_PDSvst0_no)) # 365
length(intersect(topids_t14_DMSOvst0_no, topids_t14_PhenDC3vst0_no)) # 256
length(intersect(topids_t14_PDSvst0_no, topids_t14_PhenDC3vst0_no)) # 168


# using online tool venny for venn diagrams
# http://bioinfogp.cnb.csic.es/tools/venny/
# writing t14 topids to file
write(topids_t14_DMSOvst0_no, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160510_A375_pool8_topids_t14_DMSOvst0_no.txt")
write(topids_t14_PDSvst0_no, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160510_A375_pool8_topids_t14_PDSvst0_no.txt")
write(topids_t14_PhenDC3vst0_no, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160510_A375_pool8_topids_t14_PhenDC3vst0_no.txt")
# output here: /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160510_A375_pool8_t14vst0_venn_vennyonline.png


# Intersect lists of topTags from the venn diagrams and print to files
# topids_t14_PDSvst14_DMSO, topids_t14_DMSOvst0_no
# topids_t14_PhenDC3vst14_DMSO, topids_t14_DMSOvst0_no
top_t14_PDSvst0_no_PDSonly <- top_t14_PDSvst0_no[1:length(topids_t14_PDSvst0_no),][!(topids_t14_PDSvst0_no %in% c(topids_t14_DMSOvst0_no, topids_t14_PhenDC3vst0_no)),]
write.table(top_t14_PDSvst0_no_PDSonly, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160512_A375_pool8_top_t14_PDSvst0_no_PDSonly.txt", quote = FALSE, sep = "\t")
head(as.data.frame(top_t14_PDSvst0_no_PDSonly))
z8$counts["AP3M1__sensor__1-Pool-8", c("t0_no_1", "t0_no_2", "t0_no_3", "t14_PDS_1", "t14_PDS_2", "t14_PDS_3")]
cpm(z8)["AP3M1__sensor__1-Pool-8", c("t0_no_1", "t0_no_2", "t0_no_3", "t14_PDS_1", "t14_PDS_2", "t14_PDS_3")]
z8$counts["SERPINB10__4__3-Pool-8", c("t0_no_1", "t0_no_2", "t0_no_3", "t14_PDS_1", "t14_PDS_2", "t14_PDS_3")]
cpm(z8)["SERPINB10__4__3-Pool-8", c("t0_no_1", "t0_no_2", "t0_no_3", "t14_PDS_1", "t14_PDS_2", "t14_PDS_3")]

top_t14_PhenDC3vst0_no_PhenDC3only <- top_t14_PhenDC3vst0_no[1:length(topids_t14_PhenDC3vst0_no),][!(topids_t14_PhenDC3vst0_no %in% c(topids_t14_DMSOvst0_no, topids_t14_PDSvst0_no)),]
write.table(top_t14_PhenDC3vst0_no_PhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160512_A375_pool8_top_t14_PhenDC3vst0_no_PhenDC3only.txt", quote = FALSE, sep = "\t")
head(as.data.frame(top_t14_PhenDC3vst0_no_PhenDC3only)) # HNRNPH1__3__2-Pool-8 at the top
z8$counts["HNRNPH1__3__2-Pool-8", c("t0_no_1", "t0_no_2", "t0_no_3", "t14_PhenDC3_1", "t14_PhenDC3_2", "t14_PhenDC3_3")]
cpm(z8)["HNRNPH1__3__2-Pool-8", c("t0_no_1", "t0_no_2", "t0_no_3", "t14_PhenDC3_1", "t14_PhenDC3_2", "t14_PhenDC3_3")]
tail(as.data.frame(top_t14_PhenDC3vst0_no_PhenDC3only))
z8$counts["HIRA__1__1-Pool-8", c("t0_no_1", "t0_no_2", "t0_no_3", "t14_PhenDC3_1", "t14_PhenDC3_2", "t14_PhenDC3_3")]
cpm(z8)["HIRA__1__1-Pool-8", c("t0_no_1", "t0_no_2", "t0_no_3", "t14_PhenDC3_1", "t14_PhenDC3_2", "t14_PhenDC3_3")]

top_t14_PDSvst0_no_PDSandPhenDC3 <- top_t14_PDSvst0_no[1:length(topids_t14_PDSvst0_no),][(topids_t14_PDSvst0_no %in% topids_t14_PhenDC3vst0_no),]
top_t14_PDSvst0_no_PDSandPhenDC3only <- top_t14_PDSvst0_no_PDSandPhenDC3[!(rownames(top_t14_PDSvst0_no_PDSandPhenDC3) %in% topids_t14_DMSOvst0_no),]
write.table(top_t14_PDSvst0_no_PDSandPhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160512_A375_pool8_top_t14_PDSvst0_no_PDSandPhenDC3only.txt", quote = FALSE, sep = "\t")
top_t14_PDSvst0_no_PDSandPhenDC3only # HNRNPA2B1__6__4-Pool-8 at 7th position

top_t14_DMSOvst0_no_DMSOonly <- top_t14_DMSOvst0_no[1:length(topids_t14_DMSOvst0_no),][!(topids_t14_DMSOvst0_no %in% c(topids_t14_PDSvst0_no, topids_t14_PhenDC3vst0_no)),]



# A slightly different way of doing this
length(intersect(topids_t14_PDSvst14_DMSO, topids_t14_DMSOvst0_no)) # 191
length(setdiff(topids_t14_PDSvst14_DMSO, topids_t14_DMSOvst0_no)) # 505
top_t14_PDSvst14_DMSO[1:length(topids_t14_PDSvst14_DMSO),][!(topids_t14_PDSvst14_DMSO %in% topids_t14_DMSOvst0_no),]

length(intersect(topids_t14_PhenDC3vst14_DMSO, topids_t14_DMSOvst0_no)) # 18
length(setdiff(topids_t14_PhenDC3vst14_DMSO, topids_t14_DMSOvst0_no)) # 8





# Comparisons t7 vs t0 and within t7
# Carry out Likelihood ratio tests
lrt_t7_DMSOvst0_no <- glmLRT(z8fit, contrast=my.contrasts[,"t7_DMSOvst0_no"])
lrt_t7_PDSvst0_no <- glmLRT(z8fit, contrast=my.contrasts[,"t7_PDSvst0_no"])
lrt_t7_PhenDC3vst0_no <- glmLRT(z8fit, contrast=my.contrasts[,"t7_PhenDC3vst0_no"])
lrt_t7_PDSvst7_DMSO <- glmLRT(z8fit, contrast=my.contrasts[,"t7_PDSvst7_DMSO"])
lrt_t7_PhenDC3vst7_DMSO <- glmLRT(z8fit, contrast=my.contrasts[,"t7_PhenDC3vst7_DMSO"])
lrt_t7_PhenDC3vst7_PDS <- glmLRT(z8fit, contrast=my.contrasts[,"t7_PhenDC3vst7_PDS"])


# Show top ranked hairpins
topTags(lrt_t7_DMSOvst0_no)
topTags(lrt_t7_PDSvst0_no)
topTags(lrt_t7_PhenDC3vst0_no)
topTags(lrt_t7_PDSvst7_DMSO)
topTags(lrt_t7_PhenDC3vst7_DMSO)
topTags(lrt_t7_PhenDC3vst7_PDS)


# Select hairpins with FDR < 0.05 to highlight on plot
thresh <- 0.05

top_t7_DMSOvst0_no <- topTags(lrt_t7_DMSOvst0_no, n=Inf)
topids_t7_DMSOvst0_no <- rownames(top_t7_DMSOvst0_no$table[top_t7_DMSOvst0_no$table$FDR < thresh,])
length(topids_t7_DMSOvst0_no) # 138

top_t7_PDSvst0_no <- topTags(lrt_t7_PDSvst0_no, n=Inf)
topids_t7_PDSvst0_no <- rownames(top_t7_PDSvst0_no$table[top_t7_PDSvst0_no$table$FDR < thresh,])
length(topids_t7_PDSvst0_no) # 20

top_t7_PhenDC3vst0_no <- topTags(lrt_t7_PhenDC3vst0_no, n=Inf)
topids_t7_PhenDC3vst0_no <- rownames(top_t7_PhenDC3vst0_no$table[top_t7_PhenDC3vst0_no$table$FDR < thresh,])
length(topids_t7_PhenDC3vst0_no) # 14

top_t7_PDSvst7_DMSO <- topTags(lrt_t7_PDSvst7_DMSO, n=Inf)
topids_t7_PDSvst7_DMSO <- rownames(top_t7_PDSvst7_DMSO$table[top_t7_PDSvst7_DMSO$table$FDR < thresh,])
length(topids_t7_PDSvst7_DMSO) # 65

top_t7_PhenDC3vst7_DMSO <- topTags(lrt_t7_PhenDC3vst7_DMSO, n=Inf)
topids_t7_PhenDC3vst7_DMSO <- rownames(top_t7_PhenDC3vst7_DMSO$table[top_t7_PhenDC3vst7_DMSO$table$FDR < thresh,])
length(topids_t7_PhenDC3vst7_DMSO) # 80

top_t7_PhenDC3vst7_PDS <- topTags(lrt_t7_PhenDC3vst7_PDS, n=Inf)
topids_t7_PhenDC3vst7_PDS <- rownames(top_t7_PhenDC3vst7_PDS$table[top_t7_PhenDC3vst7_PDS$table$FDR < thresh,])
length(topids_t7_PhenDC3vst7_PDS) # 14


# Write LogFC and FDR tables (only pool8)
write.table(data.table(as.data.frame(top_t7_DMSOvst0_no), keep.rownames=TRUE)[order(logFC),.(rn, logFC, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160527_A375_pool8_t7_DMSOvst0_logFCincreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
write.table(data.table(as.data.frame(top_t7_DMSOvst0_no), keep.rownames=TRUE)[order(-logFC),.(rn, logFC, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160527_A375_pool8_t7_DMSOvst0_logFCdecreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
write.table(data.table(as.data.frame(top_t7_PDSvst0_no), keep.rownames=TRUE)[order(logFC),.(rn, logFC, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160527_A375_pool8_t7_PDSvst0_logFCincreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
write.table(data.table(as.data.frame(top_t7_PDSvst0_no), keep.rownames=TRUE)[order(-logFC),.(rn, logFC, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160527_A375_pool8_t7_PDSvst0_logFCdecreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
write.table(data.table(as.data.frame(top_t7_PhenDC3vst0_no), keep.rownames=TRUE)[order(logFC),.(rn, logFC, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160527_A375_pool8_t7_PhenDC3vst0_logFCincreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)
write.table(data.table(as.data.frame(top_t7_PhenDC3vst0_no), keep.rownames=TRUE)[order(-logFC),.(rn, logFC, FDR)], file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160527_A375_pool8_t7_PhenDC3vst0_logFCdecreasing_FDR.txt", quote = FALSE, sep = "\t", row.names=FALSE)


# Plot logFC versus logCPM
plotSmear(lrt_t7_DMSOvst0_no, de.tags=topids_t7_DMSOvst0_no, main = "t7_DMSO vs t0_no")
plotSmear(lrt_t7_PDSvst0_no, de.tags=topids_t7_PDSvst0_no, main = "t7_PDS vs t0_no")
plotSmear(lrt_t7_PhenDC3vst0_no, de.tags=topids_t7_PhenDC3vst0_no, main = "t7_PhenDC3 vs t0_no")
plotSmear(lrt_t7_PDSvst7_DMSO, de.tags=topids_t7_PDSvst7_DMSO, main = "t7_PDS vs t7_DMSO")
plotSmear(lrt_t7_PhenDC3vst7_DMSO, de.tags=topids_t7_PhenDC3vst7_DMSO, main = "t7_PhenDC3 vs t7_DMSO")
plotSmear(lrt_t7_PhenDC3vst7_PDS, de.tags=topids_t7_PhenDC3vst7_PDS, main = "t7_PhenDC3 vs t7_PDS")


# Intersect lists of topTags
# t7_DMSOvst0_no, t7_PDSvst0_no, t7_PhenDC3vst0_no
length(intersect(topids_t7_DMSOvst0_no, topids_t7_PDSvst0_no)) # 10
length(intersect(topids_t7_DMSOvst0_no, topids_t7_PhenDC3vst0_no)) # 5
length(intersect(topids_t7_PDSvst0_no, topids_t7_PhenDC3vst0_no)) # 6


# using online tool venny for venn diagrams
# http://bioinfogp.cnb.csic.es/tools/venny/
# writing t7 topids to file
write(topids_t7_DMSOvst0_no, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160517_A375_pool8_topids_t7_DMSOvst0_no.txt")
write(topids_t7_PDSvst0_no, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160517_A375_pool8_topids_t7_PDSvst0_no.txt")
write(topids_t7_PhenDC3vst0_no, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160517_A375_pool8_topids_t7_PhenDC3vst0_no.txt")
# Style: Colors
# Show %
# Line 3
# Font 28
# Family: Monospace
# output here: /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160517_A375_pool8_t7vst0_venn_vennyonline.png


# Intersect lists of topTags from the venn diagrams and print to files
# topids_t7_PDSvst7_DMSO, topids_t7_DMSOvst0_no
# topids_t7_PhenDC3vst7_DMSO, topids_t7_DMSOvst0_no
top_t7_PDSvst0_no_PDSonly <- top_t7_PDSvst0_no[1:length(topids_t7_PDSvst0_no),][!(topids_t7_PDSvst0_no %in% c(topids_t7_DMSOvst0_no, topids_t7_PhenDC3vst0_no)),]
write.table(top_t7_PDSvst0_no_PDSonly, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160517_A375_pool8_top_t7_PDSvst0_no_PDSonly.txt", quote = FALSE, sep = "\t")
head(as.data.frame(top_t7_PDSvst0_no_PDSonly))
z8$counts["ACADSB__4__4-Pool-8", c("t0_no_1", "t0_no_2", "t0_no_3", "t7_PDS_1", "t7_PDS_2", "t7_PDS_3")]
cpm(z8)["ACADSB__4__4-Pool-8", c("t0_no_1", "t0_no_2", "t0_no_3", "t7_PDS_1", "t7_PDS_2", "t7_PDS_3")]
z8$counts["MRPL33__2__2-Pool-8", c("t0_no_1", "t0_no_2", "t0_no_3", "t7_PDS_1", "t7_PDS_2", "t7_PDS_3")]
cpm(z8)["MRPL33__2__2-Pool-8", c("t0_no_1", "t0_no_2", "t0_no_3", "t7_PDS_1", "t7_PDS_2", "t7_PDS_3")]

top_t7_PhenDC3vst0_no_PhenDC3only <- top_t7_PhenDC3vst0_no[1:length(topids_t7_PhenDC3vst0_no),][!(topids_t7_PhenDC3vst0_no %in% c(topids_t7_DMSOvst0_no, topids_t7_PDSvst0_no)),]
write.table(top_t7_PhenDC3vst0_no_PhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160517_A375_pool8_top_t7_PhenDC3vst0_no_PhenDC3only.txt", quote = FALSE, sep = "\t")
head(as.data.frame(top_t7_PhenDC3vst0_no_PhenDC3only)) # HNRNPL__4__4-Pool-8 is among the top
z8$counts["GNS__3__3-Pool-8", c("t0_no_1", "t0_no_2", "t0_no_3", "t7_PhenDC3_1", "t7_PhenDC3_2", "t7_PhenDC3_3")]
cpm(z8)["GNS__3__3-Pool-8", c("t0_no_1", "t0_no_2", "t0_no_3", "t7_PhenDC3_1", "t7_PhenDC3_2", "t7_PhenDC3_3")]
z8$counts["HNRNPL__4__4-Pool-8", c("t0_no_1", "t0_no_2", "t0_no_3", "t7_PhenDC3_1", "t7_PhenDC3_2", "t7_PhenDC3_3")]
cpm(z8)["HNRNPL__4__4-Pool-8", c("t0_no_1", "t0_no_2", "t0_no_3", "t7_PhenDC3_1", "t7_PhenDC3_2", "t7_PhenDC3_3")]

top_t7_PDSvst0_no_PDSandPhenDC3 <- top_t7_PDSvst0_no[1:length(topids_t7_PDSvst0_no),][(topids_t7_PDSvst0_no %in% topids_t7_PhenDC3vst0_no),]
top_t7_PDSvst0_no_PDSandPhenDC3only <- top_t7_PDSvst0_no_PDSandPhenDC3[!(rownames(top_t7_PDSvst0_no_PDSandPhenDC3) %in% topids_t7_DMSOvst0_no),]
write.table(top_t7_PDSvst0_no_PDSandPhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160517_A375_pool8_top_t7_PDSvst0_no_PDSandPhenDC3only.txt", quote = FALSE, sep = "\t")
top_t7_PDSvst0_no_PDSandPhenDC3only


# How many of the PDSonly, PhenDC3only and PDSandPhenDC3only hairpins in t7 are also present in PDSonly, PhenDC3only and PDSandPhenDC3only hairpins in t14?
intersect(rownames(top_t7_PDSvst0_no_PDSonly), rownames(top_t14_PDSvst0_no_PDSonly))
intersect(rownames(top_t7_PhenDC3vst0_no_PhenDC3only), rownames(top_t14_PhenDC3vst0_no_PhenDC3only))
intersect(rownames(top_t7_PDSvst0_no_PDSandPhenDC3only), rownames(top_t14_PDSvst0_no_PDSandPhenDC3only))

```






The plot "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160510_A375_pool8_mds.pdf" suggests that we should try removing the replicates:
- t14_PhenDC3_3 for the t14 vs. t0 comparison
- t7_PDS_3 for the t7 vs. t0 comparison

In the following R script code chunk we repeat the same as above but removing these two replicates


```R
library(edgeR)


# Load table of counts into edgeR
data <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/all_counts.txt", header = T)
dim(data) # 44763    22


# Add first column to rownames
rownames(data) <- data[,1]
data <- data[2:22]
dim(data)


# Change to a DGE object
z <- DGEList(counts=data)


# Change column labels z$counts, change row labels z$samples, add group to z$samples, add Replicate to z$samples
treatments <- rep(c(rep("DMSO", 2), rep("PDS", 2), rep("PhenDC3", 2), "no"),3)
timepoints <- rep(c(rep(rep(c("t14", "t7")),3), "t0"),3)
replicates <- c(rep("1", 7), rep("2", 7), rep("3", 7))

colnames(z$counts) <- paste(timepoints, treatments, replicates, sep="_")
rownames(z$samples) <- paste(timepoints, treatments, replicates, sep="_")
z$samples$group <- paste(timepoints, treatments, sep="_")
z$samples$Replicate <- replicates


# Select pool8 and remove t14_PhenDC3_3 and t7_PDS_3 replicates
z8 <- z[grepl("Pool-8", rownames(z$counts)), !grepl("t14_PhenDC3_3|t7_PDS_3", colnames(z$counts))]


# Plot number of hairpins that could be matched per sample
# and total for each hairpin across all samples
par(mfrow=c(2,1))
barplot(colSums(z8$counts), las=2, main="Counts per index", cex.names=0.5, cex.axis=0.8)
barplot(rowSums(z8$counts), las=2, main="Counts per hairpin", cex.names=0.5, cex.axis=0.8, axisnames = FALSE)


# Make MDS plots to visualise relationships between replicate samples
plotMDS(z8, labels = z8$samples$group, col = c(1,1,2,2,3,3,4,1,1,2,2,3,3,4,1,1,2,3,4)) # coloured by treatment
legend("bottomright", legend=c("DMSO", "PDS", "PhenDC3", "t0"), col=1:4, pch=15)

pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160519_A375_pool8_without_t14_PhenDC3_3_t7_PDS_3_mds.pdf", width = 9, height = 9)
plotMDS(z8, labels = rownames(z8$samples), col = c(1,1,2,2,3,3,4,1,1,2,2,3,3,4,1,1,2,3,4)) # coloured by treatment
legend("bottomright", legend=c("DMSO", "PDS", "PhenDC3", "t0"), col=1:4, pch=15)
dev.off()

z8_t14vst0 <- z8[,grepl("t14|t0", colnames(z8$counts))]
pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160519_A375_pool8_without_t14_PhenDC3_3_t14vst0_mds.pdf", width = 9, height = 9)
plotMDS(z8_t14vst0, labels = rownames(z8_t14vst0$samples), col = c(1,2,3,4,1,2,3,4,1,2,4)) # coloured by treatment
legend("bottomright", legend=c("DMSO", "PDS", "PhenDC3", "t0"), col=1:4, pch=15)
dev.off()

z8_t7vst0 <- z8[,grepl("t7|t0", colnames(z8$counts))]
pdf("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160519_A375_pool8_without_t7_PDS_3_t7vst0_mds.pdf", width = 9, height = 9)
plotMDS(z8_t7vst0, labels = rownames(z8_t7vst0$samples), col = c(1,2,3,4,1,2,3,4,1,3,4)) # coloured by treatment
legend("bottomright", legend=c("DMSO", "PDS", "PhenDC3", "t0"), col=1:4, pch=15)
dev.off()


# Begin differential representation analysis.
# We will use GLMs in edgeR in this case since there are more than 2 groups.
# Set up design matrix for GLM.
des <- model.matrix(~ 0 + group, data = z8$samples)
colnames(des) <- levels(factor(z8$samples$group))
des


# Estimate dispersions
z8glm <- estimateDisp(z8, des)
sqrt(z8glm$common.disp) # 0.4092689


# Plot BCVs versus abundance
plotBCV(z8glm)


# Fit negative bionomial GLM
z8fit <- glmFit(z8glm, des)


# Define matrix of contrasts
my.contrasts <- makeContrasts(
t7_DMSOvst0_no = t7_DMSO - t0_no,
t7_PDSvst0_no = t7_PDS - t0_no,
t7_PhenDC3vst0_no = t7_PhenDC3 - t0_no,
t7_PDSvst7_DMSO = t7_PDS - t7_DMSO,
t7_PhenDC3vst7_DMSO = t7_PhenDC3 - t7_DMSO,
t7_PhenDC3vst7_PDS = t7_PhenDC3 - t7_PDS,
t14_DMSOvst0_no = t14_DMSO - t0_no,
t14_PDSvst0_no = t14_PDS - t0_no,
t14_PhenDC3vst0_no = t14_PhenDC3 - t0_no,
t14_PDSvst14_DMSO = t14_PDS - t14_DMSO,
t14_PhenDC3vst14_DMSO = t14_PhenDC3 - t14_DMSO,
t14_PhenDC3vst14_PDS = t14_PhenDC3 - t14_PDS,
t14_DMSOvst7_DMSO = t14_DMSO - t7_DMSO,
t14_PDSvst7_PDS = t14_PDS - t7_PDS,
t14_PhenDC3vst7_PhenDC3 = t14_PhenDC3 - t7_PhenDC3,
levels=des)





# Comparisons t14 vs t0 and within t14
# Carry out Likelihood ratio tests
lrt_t14_DMSOvst0_no <- glmLRT(z8fit, contrast=my.contrasts[,"t14_DMSOvst0_no"])
lrt_t14_PDSvst0_no <- glmLRT(z8fit, contrast=my.contrasts[,"t14_PDSvst0_no"])
lrt_t14_PhenDC3vst0_no <- glmLRT(z8fit, contrast=my.contrasts[,"t14_PhenDC3vst0_no"])
lrt_t14_PDSvst14_DMSO <- glmLRT(z8fit, contrast=my.contrasts[,"t14_PDSvst14_DMSO"])
lrt_t14_PhenDC3vst14_DMSO <- glmLRT(z8fit, contrast=my.contrasts[,"t14_PhenDC3vst14_DMSO"])
lrt_t14_PhenDC3vst14_PDS <- glmLRT(z8fit, contrast=my.contrasts[,"t14_PhenDC3vst14_PDS"])


# Show top ranked hairpins
topTags(lrt_t14_DMSOvst0_no)
topTags(lrt_t14_PDSvst0_no)
topTags(lrt_t14_PhenDC3vst0_no)
topTags(lrt_t14_PDSvst14_DMSO)
topTags(lrt_t14_PhenDC3vst14_DMSO)
topTags(lrt_t14_PhenDC3vst14_PDS)


# Select hairpins with FDR < 0.05 to highlight on plot
thresh <- 0.05

top_t14_DMSOvst0_no <- topTags(lrt_t14_DMSOvst0_no, n=Inf)
topids_t14_DMSOvst0_no <- rownames(top_t14_DMSOvst0_no$table[top_t14_DMSOvst0_no$table$FDR < thresh,])
length(topids_t14_DMSOvst0_no) # 976

top_t14_PDSvst0_no <- topTags(lrt_t14_PDSvst0_no, n=Inf)
topids_t14_PDSvst0_no <- rownames(top_t14_PDSvst0_no$table[top_t14_PDSvst0_no$table$FDR < thresh,])
length(topids_t14_PDSvst0_no) # 1085

top_t14_PhenDC3vst0_no <- topTags(lrt_t14_PhenDC3vst0_no, n=Inf)
topids_t14_PhenDC3vst0_no <- rownames(top_t14_PhenDC3vst0_no$table[top_t14_PhenDC3vst0_no$table$FDR < thresh,])
length(topids_t14_PhenDC3vst0_no) # 385

top_t14_PDSvst14_DMSO <- topTags(lrt_t14_PDSvst14_DMSO, n=Inf)
topids_t14_PDSvst14_DMSO <- rownames(top_t14_PDSvst14_DMSO$table[top_t14_PDSvst14_DMSO$table$FDR < thresh,])
length(topids_t14_PDSvst14_DMSO) # 739

top_t14_PhenDC3vst14_DMSO <- topTags(lrt_t14_PhenDC3vst14_DMSO, n=Inf)
topids_t14_PhenDC3vst14_DMSO <- rownames(top_t14_PhenDC3vst14_DMSO$table[top_t14_PhenDC3vst14_DMSO$table$FDR < thresh,])
length(topids_t14_PhenDC3vst14_DMSO) # 25

top_t14_PhenDC3vst14_PDS <- topTags(lrt_t14_PhenDC3vst14_PDS, n=Inf)
topids_t14_PhenDC3vst14_PDS <- rownames(top_t14_PhenDC3vst14_PDS$table[top_t14_PhenDC3vst14_PDS$table$FDR < thresh,])
length(topids_t14_PhenDC3vst14_PDS) # 441


# Plot logFC versus logCPM
plotSmear(lrt_t14_DMSOvst0_no, de.tags=topids_t14_DMSOvst0_no, main = "t14_DMSO vs t0_no")
plotSmear(lrt_t14_PDSvst0_no, de.tags=topids_t14_PDSvst0_no, main = "t14_PDS vs t0_no")
plotSmear(lrt_t14_PhenDC3vst0_no, de.tags=topids_t14_PhenDC3vst0_no, main = "t14_PhenDC3 vs t0_no")
plotSmear(lrt_t14_PDSvst14_DMSO, de.tags=topids_t14_PDSvst14_DMSO, main = "t14_PDS vs t14_DMSO")
plotSmear(lrt_t14_PhenDC3vst14_DMSO, de.tags=topids_t14_PhenDC3vst14_DMSO, main = "t14_PhenDC3 vs t14_DMSO")
plotSmear(lrt_t14_PhenDC3vst14_PDS, de.tags=topids_t14_PhenDC3vst14_PDS, main = "t14_PhenDC3 vs t14_PDS")


# Intersect lists of topTags
# t14_DMSOvst0_no, t14_PDSvst0_no, t14_PhenDC3vst0_no
length(intersect(topids_t14_DMSOvst0_no, topids_t14_PDSvst0_no)) # 386
length(intersect(topids_t14_DMSOvst0_no, topids_t14_PhenDC3vst0_no)) # 284
length(intersect(topids_t14_PDSvst0_no, topids_t14_PhenDC3vst0_no)) # 188


# using online tool venny for venn diagrams
# writing t14 topids to file
write(topids_t14_DMSOvst0_no, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t14_PhenDC3_3_topids_t14_DMSOvst0_no.txt")
write(topids_t14_PDSvst0_no, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t14_PhenDC3_3_topids_t14_PDSvst0_no.txt")
write(topids_t14_PhenDC3vst0_no, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t14_PhenDC3_3_topids_t14_PhenDC3vst0_no.txt")
# Style: Colors
# Show %
# Line 3
# Font 28
# Family: Monospace
# output here: /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160519_A375_pool8_without_t14_PhenDC3_3_t14vst0_venn_vennyonline.png


# Intersect lists of topTags from the venn diagrams and print to files
# topids_t14_PDSvst14_DMSO, topids_t14_DMSOvst0_no
# topids_t14_PhenDC3vst14_DMSO, topids_t14_DMSOvst0_no
top_t14_PDSvst0_no_PDSonly <- top_t14_PDSvst0_no[1:length(topids_t14_PDSvst0_no),][!(topids_t14_PDSvst0_no %in% c(topids_t14_DMSOvst0_no, topids_t14_PhenDC3vst0_no)),]
write.table(top_t14_PDSvst0_no_PDSonly, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t14_PhenDC3_3_top_t14_PDSvst0_no_PDSonly.txt", quote = FALSE, sep = "\t")
head(as.data.frame(top_t14_PDSvst0_no_PDSonly), 20)

top_t14_PhenDC3vst0_no_PhenDC3only <- top_t14_PhenDC3vst0_no[1:length(topids_t14_PhenDC3vst0_no),][!(topids_t14_PhenDC3vst0_no %in% c(topids_t14_DMSOvst0_no, topids_t14_PDSvst0_no)),]
write.table(top_t14_PhenDC3vst0_no_PhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t14_PhenDC3_3_top_t14_PhenDC3vst0_no_PhenDC3only.txt", quote = FALSE, sep = "\t")
head(as.data.frame(top_t14_PhenDC3vst0_no_PhenDC3only), 20)

top_t14_PDSvst0_no_PDSandPhenDC3 <- top_t14_PDSvst0_no[1:length(topids_t14_PDSvst0_no),][(topids_t14_PDSvst0_no %in% topids_t14_PhenDC3vst0_no),]
top_t14_PDSvst0_no_PDSandPhenDC3only <- top_t14_PDSvst0_no_PDSandPhenDC3[!(rownames(top_t14_PDSvst0_no_PDSandPhenDC3) %in% topids_t14_DMSOvst0_no),]
write.table(top_t14_PDSvst0_no_PDSandPhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t14_PhenDC3_3_top_t14_PDSvst0_no_PDSandPhenDC3only.txt", quote = FALSE, sep = "\t")
head(as.data.frame(top_t14_PDSvst0_no_PDSandPhenDC3only), 20)





# Comparisons t7 vs t0 and within t7
# Carry out Likelihood ratio tests
lrt_t7_DMSOvst0_no <- glmLRT(z8fit, contrast=my.contrasts[,"t7_DMSOvst0_no"])
lrt_t7_PDSvst0_no <- glmLRT(z8fit, contrast=my.contrasts[,"t7_PDSvst0_no"])
lrt_t7_PhenDC3vst0_no <- glmLRT(z8fit, contrast=my.contrasts[,"t7_PhenDC3vst0_no"])
lrt_t7_PDSvst7_DMSO <- glmLRT(z8fit, contrast=my.contrasts[,"t7_PDSvst7_DMSO"])
lrt_t7_PhenDC3vst7_DMSO <- glmLRT(z8fit, contrast=my.contrasts[,"t7_PhenDC3vst7_DMSO"])
lrt_t7_PhenDC3vst7_PDS <- glmLRT(z8fit, contrast=my.contrasts[,"t7_PhenDC3vst7_PDS"])


# Show top ranked hairpins
topTags(lrt_t7_DMSOvst0_no)
topTags(lrt_t7_PDSvst0_no)
topTags(lrt_t7_PhenDC3vst0_no)
topTags(lrt_t7_PDSvst7_DMSO)
topTags(lrt_t7_PhenDC3vst7_DMSO)
topTags(lrt_t7_PhenDC3vst7_PDS)


# Select hairpins with FDR < 0.05 to highlight on plot
thresh <- 0.05

top_t7_DMSOvst0_no <- topTags(lrt_t7_DMSOvst0_no, n=Inf)
topids_t7_DMSOvst0_no <- rownames(top_t7_DMSOvst0_no$table[top_t7_DMSOvst0_no$table$FDR < thresh,])
length(topids_t7_DMSOvst0_no) # 169

top_t7_PDSvst0_no <- topTags(lrt_t7_PDSvst0_no, n=Inf)
topids_t7_PDSvst0_no <- rownames(top_t7_PDSvst0_no$table[top_t7_PDSvst0_no$table$FDR < thresh,])
length(topids_t7_PDSvst0_no) # 34

top_t7_PhenDC3vst0_no <- topTags(lrt_t7_PhenDC3vst0_no, n=Inf)
topids_t7_PhenDC3vst0_no <- rownames(top_t7_PhenDC3vst0_no$table[top_t7_PhenDC3vst0_no$table$FDR < thresh,])
length(topids_t7_PhenDC3vst0_no) # 21

top_t7_PDSvst7_DMSO <- topTags(lrt_t7_PDSvst7_DMSO, n=Inf)
topids_t7_PDSvst7_DMSO <- rownames(top_t7_PDSvst7_DMSO$table[top_t7_PDSvst7_DMSO$table$FDR < thresh,])
length(topids_t7_PDSvst7_DMSO) # 11

top_t7_PhenDC3vst7_DMSO <- topTags(lrt_t7_PhenDC3vst7_DMSO, n=Inf)
topids_t7_PhenDC3vst7_DMSO <- rownames(top_t7_PhenDC3vst7_DMSO$table[top_t7_PhenDC3vst7_DMSO$table$FDR < thresh,])
length(topids_t7_PhenDC3vst7_DMSO) # 86

top_t7_PhenDC3vst7_PDS <- topTags(lrt_t7_PhenDC3vst7_PDS, n=Inf)
topids_t7_PhenDC3vst7_PDS <- rownames(top_t7_PhenDC3vst7_PDS$table[top_t7_PhenDC3vst7_PDS$table$FDR < thresh,])
length(topids_t7_PhenDC3vst7_PDS) # 7


# Plot logFC versus logCPM
plotSmear(lrt_t7_DMSOvst0_no, de.tags=topids_t7_DMSOvst0_no, main = "t7_DMSO vs t0_no")
plotSmear(lrt_t7_PDSvst0_no, de.tags=topids_t7_PDSvst0_no, main = "t7_PDS vs t0_no")
plotSmear(lrt_t7_PhenDC3vst0_no, de.tags=topids_t7_PhenDC3vst0_no, main = "t7_PhenDC3 vs t0_no")
plotSmear(lrt_t7_PDSvst7_DMSO, de.tags=topids_t7_PDSvst7_DMSO, main = "t7_PDS vs t7_DMSO")
plotSmear(lrt_t7_PhenDC3vst7_DMSO, de.tags=topids_t7_PhenDC3vst7_DMSO, main = "t7_PhenDC3 vs t7_DMSO")
plotSmear(lrt_t7_PhenDC3vst7_PDS, de.tags=topids_t7_PhenDC3vst7_PDS, main = "t7_PhenDC3 vs t7_PDS")


# Intersect lists of topTags
# t7_DMSOvst0_no, t7_PDSvst0_no, t7_PhenDC3vst0_no
length(intersect(topids_t7_DMSOvst0_no, topids_t7_PDSvst0_no)) # 13
length(intersect(topids_t7_DMSOvst0_no, topids_t7_PhenDC3vst0_no)) # 11
length(intersect(topids_t7_PDSvst0_no, topids_t7_PhenDC3vst0_no)) # 6


# using online tool venny for venn diagrams
# writing t7 topids to file
write(topids_t7_DMSOvst0_no, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t7_PDS_3_topids_t7_DMSOvst0_no.txt")
write(topids_t7_PDSvst0_no, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t7_PDS_3_topids_t7_PDSvst0_no.txt")
write(topids_t7_PhenDC3vst0_no, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t7_PDS_3_topids_t7_PhenDC3vst0_no.txt")
# Style: Colors
# Show %
# Line 3
# Font 28
# Family: Monospace
# output here: /lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160519_A375_pool8_without_t7_PDS_3_t7vst0_venn_vennyonline.png


# Intersect lists of topTags from the venn diagrams and print to files
# topids_t7_PDSvst7_DMSO, topids_t7_DMSOvst0_no
# topids_t7_PhenDC3vst7_DMSO, topids_t7_DMSOvst0_no
top_t7_PDSvst0_no_PDSonly <- top_t7_PDSvst0_no[1:length(topids_t7_PDSvst0_no),][!(topids_t7_PDSvst0_no %in% c(topids_t7_DMSOvst0_no, topids_t7_PhenDC3vst0_no)),]
write.table(top_t7_PDSvst0_no_PDSonly, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t7_PDS_3_top_t7_PDSvst0_no_PDSonly.txt", quote = FALSE, sep = "\t")
top_t7_PDSvst0_no_PDSonly

top_t7_PhenDC3vst0_no_PhenDC3only <- top_t7_PhenDC3vst0_no[1:length(topids_t7_PhenDC3vst0_no),][!(topids_t7_PhenDC3vst0_no %in% c(topids_t7_DMSOvst0_no, topids_t7_PDSvst0_no)),]
write.table(top_t7_PhenDC3vst0_no_PhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t7_PDS_3_top_t7_PhenDC3vst0_no_PhenDC3only.txt", quote = FALSE, sep = "\t")
top_t7_PhenDC3vst0_no_PhenDC3only

top_t7_PDSvst0_no_PDSandPhenDC3 <- top_t7_PDSvst0_no[1:length(topids_t7_PDSvst0_no),][(topids_t7_PDSvst0_no %in% topids_t7_PhenDC3vst0_no),]
top_t7_PDSvst0_no_PDSandPhenDC3only <- top_t7_PDSvst0_no_PDSandPhenDC3[!(rownames(top_t7_PDSvst0_no_PDSandPhenDC3) %in% topids_t7_DMSOvst0_no),]
write.table(top_t7_PDSvst0_no_PDSandPhenDC3only, file = "/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t7_PDS_3_top_t7_PDSvst0_no_PDSandPhenDC3only.txt", quote = FALSE, sep = "\t")
top_t7_PDSvst0_no_PDSandPhenDC3only


# How many of the PDSonly, PhenDC3only and PDSandPhenDC3only hairpins in t7 are also present in PDSonly, PhenDC3only and PDSandPhenDC3only hairpins in t14?
intersect(rownames(top_t7_PDSvst0_no_PDSonly), rownames(top_t14_PDSvst0_no_PDSonly))
intersect(rownames(top_t7_PhenDC3vst0_no_PhenDC3only), rownames(top_t14_PhenDC3vst0_no_PhenDC3only))
intersect(rownames(top_t7_PDSvst0_no_PDSandPhenDC3only), rownames(top_t14_PDSvst0_no_PDSandPhenDC3only))

```




How many of the PDSonly, PhenDC3only and PDSandPhenDC3only hairpins at t7 and t14 (after removing replicates) are found when we keep the replicates?

I can do it directly in python with the lists:

/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160517_A375_pool8_top_t7_PDSvst0_no_PDSonly.txt
/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160517_A375_pool8_top_t7_PhenDC3vst0_no_PhenDC3only.txt
/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160517_A375_pool8_top_t7_PDSvst0_no_PDSandPhenDC3only.txt
/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160512_A375_pool8_top_t14_PDSvst0_no_PDSonly.txt
/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160512_A375_pool8_top_t14_PhenDC3vst0_no_PhenDC3only.txt
/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160512_A375_pool8_top_t14_PDSvst0_no_PDSandPhenDC3only.txt

/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t7_PDS_3_top_t7_PDSvst0_no_PDSonly.txt
/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t7_PDS_3_top_t7_PhenDC3vst0_no_PhenDC3only.txt
/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t7_PDS_3_top_t7_PDSvst0_no_PDSandPhenDC3only.txt
/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t14_PhenDC3_3_top_t14_PDSvst0_no_PDSonly.txt
/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t14_PhenDC3_3_top_t14_PhenDC3vst0_no_PhenDC3only.txt
/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t14_PhenDC3_3_top_t14_PDSvst0_no_PDSandPhenDC3only.txt


```python

# t7 PDS only

t7_pdsonly_ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160517_A375_pool8_top_t7_PDSvst0_no_PDSonly.txt", "r")
t7_pdsonly_lines = t7_pdsonly_ifile.readlines()
t7_pdsonly_ifile.close()

t7_pdsonly_hairpins = []

for line in t7_pdsonly_lines[1:]:
	fields = line.split()
	hairpin = fields[0]
	t7_pdsonly_hairpins.append(hairpin)


t7_pdsonly_not7pds3_ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t7_PDS_3_top_t7_PDSvst0_no_PDSonly.txt", "r")
t7_pdsonly_not7pds3_lines = t7_pdsonly_not7pds3_ifile.readlines()
t7_pdsonly_not7pds3_ifile.close()

t7_pdsonly_not7pds3_hairpins = []

for line in t7_pdsonly_not7pds3_lines[1:]:
	fields = line.split()
	hairpin = fields[0]
	t7_pdsonly_not7pds3_hairpins.append(hairpin)


len(t7_pdsonly_not7pds3_hairpins) # 20
len(t7_pdsonly_hairpins) # 8
len(set(t7_pdsonly_not7pds3_hairpins).intersection(set(t7_pdsonly_hairpins))) # 3



# t7 PhenDC3 only

t7_phendc3only_ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160517_A375_pool8_top_t7_PhenDC3vst0_no_PhenDC3only.txt", "r")
t7_phendc3only_lines = t7_phendc3only_ifile.readlines()
t7_phendc3only_ifile.close()

t7_phendc3only_hairpins = []

for line in t7_phendc3only_lines[1:]:
	fields = line.split()
	hairpin = fields[0]
	t7_phendc3only_hairpins.append(hairpin)


t7_phendc3only_not7pds3_ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t7_PDS_3_top_t7_PhenDC3vst0_no_PhenDC3only.txt", "r")
t7_phendc3only_not7pds3_lines = t7_phendc3only_not7pds3_ifile.readlines()
t7_phendc3only_not7pds3_ifile.close()

t7_phendc3only_not7pds3_hairpins = []

for line in t7_phendc3only_not7pds3_lines[1:]:
	fields = line.split()
	hairpin = fields[0]
	t7_phendc3only_not7pds3_hairpins.append(hairpin)


len(t7_phendc3only_not7pds3_hairpins) # 9
len(t7_phendc3only_hairpins) # 7
len(set(t7_phendc3only_not7pds3_hairpins).intersection(set(t7_phendc3only_hairpins))) # 6



# t7 PDS+PhenDC3 only

t7_pdsandphendc3only_ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160517_A375_pool8_top_t7_PDSvst0_no_PDSandPhenDC3only.txt", "r")
t7_pdsandphendc3only_lines = t7_pdsandphendc3only_ifile.readlines()
t7_pdsandphendc3only_ifile.close()

t7_pdsandphendc3only_hairpins = []

for line in t7_pdsandphendc3only_lines[1:]:
	fields = line.split()
	hairpin = fields[0]
	t7_pdsandphendc3only_hairpins.append(hairpin)


t7_pdsandphendc3only_not7pds3_ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t7_PDS_3_top_t7_PDSvst0_no_PDSandPhenDC3only.txt", "r")
t7_pdsandphendc3only_not7pds3_lines = t7_pdsandphendc3only_not7pds3_ifile.readlines()
t7_pdsandphendc3only_not7pds3_ifile.close()

t7_pdsandphendc3only_not7pds3_hairpins = []

for line in t7_pdsandphendc3only_not7pds3_lines[1:]:
	fields = line.split()
	hairpin = fields[0]
	t7_pdsandphendc3only_not7pds3_hairpins.append(hairpin)


len(t7_pdsandphendc3only_not7pds3_hairpins) # 1
len(t7_pdsandphendc3only_hairpins) # 2
len(set(t7_pdsandphendc3only_not7pds3_hairpins).intersection(set(t7_pdsandphendc3only_hairpins))) # 0



# t14 PDS only

t14_pdsonly_ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160512_A375_pool8_top_t14_PDSvst0_no_PDSonly.txt", "r")
t14_pdsonly_lines = t14_pdsonly_ifile.readlines()
t14_pdsonly_ifile.close()

t14_pdsonly_hairpins = []

for line in t14_pdsonly_lines[1:]:
	fields = line.split()
	hairpin = fields[0]
	t14_pdsonly_hairpins.append(hairpin)


t14_pdsonly_not14phendc3_ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t14_PhenDC3_3_top_t14_PDSvst0_no_PDSonly.txt", "r")
t14_pdsonly_not14phendc3_lines = t14_pdsonly_not14phendc3_ifile.readlines()
t14_pdsonly_not14phendc3_ifile.close()

t14_pdsonly_not14phendc3_hairpins = []

for line in t14_pdsonly_not14phendc3_lines[1:]:
	fields = line.split()
	hairpin = fields[0]
	t14_pdsonly_not14phendc3_hairpins.append(hairpin)


len(t14_pdsonly_not14phendc3_hairpins) # 667
len(t14_pdsonly_hairpins) # 695
len(set(t14_pdsonly_not14phendc3_hairpins).intersection(set(t14_pdsonly_hairpins))) # 617



# t14 PhenDC3 only

t14_phendc3only_ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160512_A375_pool8_top_t14_PhenDC3vst0_no_PhenDC3only.txt", "r")
t14_phendc3only_lines = t14_phendc3only_ifile.readlines()
t14_phendc3only_ifile.close()

t14_phendc3only_hairpins = []

for line in t14_phendc3only_lines[1:]:
	fields = line.split()
	hairpin = fields[0]
	t14_phendc3only_hairpins.append(hairpin)


t14_phendc3only_not14phendc3_ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t14_PhenDC3_3_top_t14_PhenDC3vst0_no_PhenDC3only.txt", "r")
t14_phendc3only_not14phendc3_lines = t14_phendc3only_not14phendc3_ifile.readlines()
t14_phendc3only_not14phendc3_ifile.close()

t14_phendc3only_not14phendc3_hairpins = []

for line in t14_phendc3only_not14phendc3_lines[1:]:
	fields = line.split()
	hairpin = fields[0]
	t14_phendc3only_not14phendc3_hairpins.append(hairpin)


len(t14_phendc3only_not14phendc3_hairpins) # 69
len(t14_phendc3only_hairpins) # 63
len(set(t14_phendc3only_not14phendc3_hairpins).intersection(set(t14_phendc3only_hairpins))) # 29



# t14 PDS+PhenDC3 only

t14_pdsandphendc3only_ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160512_A375_pool8_top_t14_PDSvst0_no_PDSandPhenDC3only.txt", "r")
t14_pdsandphendc3only_lines = t14_pdsandphendc3only_ifile.readlines()
t14_pdsandphendc3only_ifile.close()

t14_pdsandphendc3only_hairpins = []

for line in t14_pdsandphendc3only_lines[1:]:
	fields = line.split()
	hairpin = fields[0]
	t14_pdsandphendc3only_hairpins.append(hairpin)


t14_pdsandphendc3only_not14phendc3_ifile = open("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/tables/20160519_A375_pool8_without_t14_PhenDC3_3_top_t14_PDSvst0_no_PDSandPhenDC3only.txt", "r")
t14_pdsandphendc3only_not14phendc3_lines = t14_pdsandphendc3only_not14phendc3_ifile.readlines()
t14_pdsandphendc3only_not14phendc3_ifile.close()

t14_pdsandphendc3only_not14phendc3_hairpins = []

for line in t14_pdsandphendc3only_not14phendc3_lines[1:]:
	fields = line.split()
	hairpin = fields[0]
	t14_pdsandphendc3only_not14phendc3_hairpins.append(hairpin)


len(t14_pdsandphendc3only_not14phendc3_hairpins) # 32
len(t14_pdsandphendc3only_hairpins) # 31
len(set(t14_pdsandphendc3only_not14phendc3_hairpins).intersection(set(t14_pdsandphendc3only_hairpins))) # 16

```



Produce xyplots of cpms as discussed with Darcie



```R
library(edgeR)
library(data.table)
library(ggplot2)

# Load table of counts into edgeR
data <- read.table("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/data/20160503/all_counts.txt", header = T)

# Add first column to rownames
rownames(data) <- data[,1]
data <- data[2:22]

# Change to a DGE object
z <- DGEList(counts=data)

# Change column labels z$counts, change row labels z$samples, add group to z$samples, add Replicate to z$samples
treatments <- rep(c(rep("DMSO", 2), rep("PDS", 2), rep("PhenDC3", 2), "no"),3)
timepoints <- rep(c(rep(rep(c("t14", "t7")),3), "t0"),3)
replicates <- c(rep("1", 7), rep("2", 7), rep("3", 7))

colnames(z$counts) <- paste(timepoints, treatments, replicates, sep="_")
rownames(z$samples) <- paste(timepoints, treatments, replicates, sep="_")
z$samples$group <- paste(timepoints, treatments, sep="_")
z$samples$Replicate <- replicates

# Select pool8
z8 <- z[grepl("Pool-8", rownames(z$counts)),]
z8$samples$lib.size <- as.vector(colSums(z8$counts))

# t14 vs t0
t14_DMSOvst0 <- data.table(cpm(z8)[, grepl("t14_DMSO|t0", colnames(z8$counts))])
t14_DMSOvst0[, c("t14_DMSO", "t0_no") := list(log2((t14_DMSO_1 + t14_DMSO_2 + t14_DMSO_3)/3), log2((t0_no_1 + t0_no_2 + t0_no_3)/3))]
ggplot(t14_DMSOvst0, aes(x = t0_no, y = t14_DMSO)) + geom_point(alpha = 0.05, size=0.5) + xlab("t0") + ylab("t14 DMSO") + theme_classic() + ggtitle("log2(counts per million)") + geom_abline(intercept=0, slope=1, colour = "red")
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160526_A375_pool8_t14_DMSOvst0_xyplot.pdf", width = 10/2.54, height = 10/2.54)

t14_PDSvst0 <- data.table(cpm(z8)[, grepl("t14_PDS|t0", colnames(z8$counts))])
t14_PDSvst0[, c("t14_PDS", "t0_no") := list(log2((t14_PDS_1 + t14_PDS_2 + t14_PDS_3)/3), log2((t0_no_1 + t0_no_2 + t0_no_3)/3))]
ggplot(t14_PDSvst0, aes(x = t0_no, y = t14_PDS)) + geom_point(alpha = 0.05, size=0.5) + xlab("t0") + ylab("t14 PDS") + theme_classic() + ggtitle("log2(counts per million)") + geom_abline(intercept=0, slope=1, colour = "red")
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160526_A375_pool8_t14_PDSvst0_xyplot.pdf", width = 10/2.54, height = 10/2.54)

t14_PhenDC3vst0 <- data.table(cpm(z8)[, grepl("t14_PhenDC3|t0", colnames(z8$counts))])
t14_PhenDC3vst0[, c("t14_PhenDC3", "t0_no") := list(log2((t14_PhenDC3_1 + t14_PhenDC3_2 + t14_PhenDC3_3)/3), log2((t0_no_1 + t0_no_2 + t0_no_3)/3))]
ggplot(t14_PhenDC3vst0, aes(x = t0_no, y = t14_PhenDC3)) + geom_point(alpha = 0.05, size=0.5) + xlab("t0") + ylab("t14 PhenDC3") + theme_classic() + ggtitle("log2(counts per million)") + geom_abline(intercept=0, slope=1, colour = "red")
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160526_A375_pool8_t14_PhenDC3vst0_xyplot.pdf", width = 10/2.54, height = 10/2.54)

# t7 vs t0
t7_DMSOvst0 <- data.table(cpm(z8)[, grepl("t7_DMSO|t0", colnames(z8$counts))])
t7_DMSOvst0[, c("t7_DMSO", "t0_no") := list(log2((t7_DMSO_1 + t7_DMSO_2 + t7_DMSO_3)/3), log2((t0_no_1 + t0_no_2 + t0_no_3)/3))]
ggplot(t7_DMSOvst0, aes(x = t0_no, y = t7_DMSO)) + geom_point(alpha = 0.05, size=0.5) + xlab("t0") + ylab("t7 DMSO") + theme_classic() + ggtitle("log2(counts per million)") + geom_abline(intercept=0, slope=1, colour = "red")
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160526_A375_pool8_t7_DMSOvst0_xyplot.pdf", width = 10/2.54, height = 10/2.54)

t7_PDSvst0 <- data.table(cpm(z8)[, grepl("t7_PDS|t0", colnames(z8$counts))])
t7_PDSvst0[, c("t7_PDS", "t0_no") := list(log2((t7_PDS_1 + t7_PDS_2 + t7_PDS_3)/3), log2((t0_no_1 + t0_no_2 + t0_no_3)/3))]
ggplot(t7_PDSvst0, aes(x = t0_no, y = t7_PDS)) + geom_point(alpha = 0.05, size=0.5) + xlab("t0") + ylab("t7 PDS") + theme_classic() + ggtitle("log2(counts per million)") + geom_abline(intercept=0, slope=1, colour = "red")
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160526_A375_pool8_t7_PDSvst0_xyplot.pdf", width = 10/2.54, height = 10/2.54)

t7_PhenDC3vst0 <- data.table(cpm(z8)[, grepl("t7_PhenDC3|t0", colnames(z8$counts))])
t7_PhenDC3vst0[, c("t7_PhenDC3", "t0_no") := list(log2((t7_PhenDC3_1 + t7_PhenDC3_2 + t7_PhenDC3_3)/3), log2((t0_no_1 + t0_no_2 + t0_no_3)/3))]
ggplot(t7_PhenDC3vst0, aes(x = t0_no, y = t7_PhenDC3)) + geom_point(alpha = 0.05, size=0.5) + xlab("t0") + ylab("t7 PhenDC3") + theme_classic() + ggtitle("log2(counts per million)") + geom_abline(intercept=0, slope=1, colour = "red")
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160526_A375_pool8_t7_PhenDC3vst0_xyplot.pdf", width = 10/2.54, height = 10/2.54)


```




Produce wave plots - as in Figure 4e in Zuber2011


```R
# t14 vs t0
t14_DMSOvst0_logFCincreasing <- data.table(as.data.frame(top_t14_DMSOvst0_no), keep.rownames=TRUE)[order(logFC)]
t14_DMSOvst0_logFCincreasing$rn2 <- reorder(t14_DMSOvst0_logFCincreasing$rn, t14_DMSOvst0_logFCincreasing$logFC)
ggplot(t14_DMSOvst0_logFCincreasing, aes(rn2, logFC, group = 1)) + geom_line() + xlab("9583 shRNAs") + ylab("log2FC") + theme_classic() + ggtitle("t14 DMSO vs. t0") + geom_abline(intercept=0, slope=0) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank()) + coord_cartesian(ylim = c(-5, 5))
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160527_A375_pool8_t14_DMSOvst0_logFC.pdf", width = 10/2.54, height = 8/2.54)

t14_PDSvst0_logFCincreasing <- data.table(as.data.frame(top_t14_PDSvst0_no), keep.rownames=TRUE)[order(logFC)]
t14_PDSvst0_logFCincreasing$rn2 <- reorder(t14_PDSvst0_logFCincreasing$rn, t14_PDSvst0_logFCincreasing$logFC)
ggplot(t14_PDSvst0_logFCincreasing, aes(rn2, logFC, group = 1)) + geom_line() + xlab("9583 shRNAs") + ylab("log2FC") + theme_classic() + ggtitle("t14 PDS vs. t0") + geom_abline(intercept=0, slope=0) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank()) + coord_cartesian(ylim = c(-5, 5))
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160527_A375_pool8_t14_PDSvst0_logFC.pdf", width = 10/2.54, height = 8/2.54)

t14_PhenDC3vst0_logFCincreasing <- data.table(as.data.frame(top_t14_PhenDC3vst0_no), keep.rownames=TRUE)[order(logFC)]
t14_PhenDC3vst0_logFCincreasing$rn2 <- reorder(t14_PhenDC3vst0_logFCincreasing$rn, t14_PhenDC3vst0_logFCincreasing$logFC)
ggplot(t14_PhenDC3vst0_logFCincreasing, aes(rn2, logFC, group = 1)) + geom_line() + xlab("9583 shRNAs") + ylab("log2FC") + theme_classic() + ggtitle("t14 PhenDC3 vs. t0") + geom_abline(intercept=0, slope=0) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank()) + coord_cartesian(ylim = c(-5, 5))
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160527_A375_pool8_t14_PhenDC3vst0_logFC.pdf", width = 10/2.54, height = 8/2.54)

# t7 vs t0
t7_DMSOvst0_logFCincreasing <- data.table(as.data.frame(top_t7_DMSOvst0_no), keep.rownames=TRUE)[order(logFC)]
t7_DMSOvst0_logFCincreasing$rn2 <- reorder(t7_DMSOvst0_logFCincreasing$rn, t7_DMSOvst0_logFCincreasing$logFC)
ggplot(t7_DMSOvst0_logFCincreasing, aes(rn2, logFC, group = 1)) + geom_line() + xlab("9583 shRNAs") + ylab("log2FC") + theme_classic() + ggtitle("t7 DMSO vs. t0") + geom_abline(intercept=0, slope=0) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank()) + coord_cartesian(ylim = c(-5, 5))
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160527_A375_pool8_t7_DMSOvst0_logFC.pdf", width = 10/2.54, height = 8/2.54)

t7_PDSvst0_logFCincreasing <- data.table(as.data.frame(top_t7_PDSvst0_no), keep.rownames=TRUE)[order(logFC)]
t7_PDSvst0_logFCincreasing$rn2 <- reorder(t7_PDSvst0_logFCincreasing$rn, t7_PDSvst0_logFCincreasing$logFC)
ggplot(t7_PDSvst0_logFCincreasing, aes(rn2, logFC, group = 1)) + geom_line() + xlab("9583 shRNAs") + ylab("log2FC") + theme_classic() + ggtitle("t7 PDS vs. t0") + geom_abline(intercept=0, slope=0) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank()) + coord_cartesian(ylim = c(-5, 5))
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160527_A375_pool8_t7_PDSvst0_logFC.pdf", width = 10/2.54, height = 8/2.54)

t7_PhenDC3vst0_logFCincreasing <- data.table(as.data.frame(top_t7_PhenDC3vst0_no), keep.rownames=TRUE)[order(logFC)]
t7_PhenDC3vst0_logFCincreasing$rn2 <- reorder(t7_PhenDC3vst0_logFCincreasing$rn, t7_PhenDC3vst0_logFCincreasing$logFC)
ggplot(t7_PhenDC3vst0_logFCincreasing, aes(rn2, logFC, group = 1)) + geom_line() + xlab("9583 shRNAs") + ylab("log2FC") + theme_classic() + ggtitle("t7 PhenDC3 vs. t0") + geom_abline(intercept=0, slope=0) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank()) + coord_cartesian(ylim = c(-5, 5))
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160527_A375_pool8_t7_PhenDC3vst0_logFC.pdf", width = 10/2.54, height = 8/2.54)


# All wave plots into one
# t14 vs t0
t14_DMSOvst0_logFCincreasing$rn3 <- 1:9583
t14_DMSOvst0_logFCincreasing$rn4 <- reorder(t14_DMSOvst0_logFCincreasing$rn3, t14_DMSOvst0_logFCincreasing$logFC)
t14_PDSvst0_logFCincreasing$rn3 <- 1:9583
t14_PDSvst0_logFCincreasing$rn4 <- reorder(t14_PDSvst0_logFCincreasing$rn3, t14_PDSvst0_logFCincreasing$logFC)
t14_PhenDC3vst0_logFCincreasing$rn3 <- 1:9583
t14_PhenDC3vst0_logFCincreasing$rn4 <- reorder(t14_PhenDC3vst0_logFCincreasing$rn3, t14_PhenDC3vst0_logFCincreasing$logFC)
l = list(t14_DMSOvst0_logFCincreasing, t14_PDSvst0_logFCincreasing, t14_PhenDC3vst0_logFCincreasing)
setattr(l, 'names', c("DMSO", "PDS", "PhenDC3"))
t14 <- rbindlist(l, use.names=TRUE, idcol="condition")
ggplot(t14, aes(rn4, logFC, group = condition)) + geom_line(aes(colour = condition)) + xlab("9583 shRNAs") + ylab("log2FC") + theme_classic() + ggtitle("t14 vs. t0") + geom_abline(intercept=0, slope=0) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank(), legend.title=element_blank()) + coord_cartesian(ylim = c(-5, 5)) + scale_colour_manual(values=c(rgb(127,255,127, max = 255), "yellow3", rgb(127,127,255, max = 255))) + scale_y_continuous(breaks=c(-5,-3,-1,1,3,5))
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160527_A375_pool8_t14vst0_logFC.pdf", width = 12/2.54, height = 8/2.54)

# t7 vs t0
t7_DMSOvst0_logFCincreasing$rn3 <- 1:9583
t7_DMSOvst0_logFCincreasing$rn4 <- reorder(t7_DMSOvst0_logFCincreasing$rn3, t7_DMSOvst0_logFCincreasing$logFC)
t7_PDSvst0_logFCincreasing$rn3 <- 1:9583
t7_PDSvst0_logFCincreasing$rn4 <- reorder(t7_PDSvst0_logFCincreasing$rn3, t7_PDSvst0_logFCincreasing$logFC)
t7_PhenDC3vst0_logFCincreasing$rn3 <- 1:9583
t7_PhenDC3vst0_logFCincreasing$rn4 <- reorder(t7_PhenDC3vst0_logFCincreasing$rn3, t7_PhenDC3vst0_logFCincreasing$logFC)
l = list(t7_DMSOvst0_logFCincreasing, t7_PDSvst0_logFCincreasing, t7_PhenDC3vst0_logFCincreasing)
setattr(l, 'names', c("DMSO", "PDS", "PhenDC3"))
t7 <- rbindlist(l, use.names=TRUE, idcol="condition")
ggplot(t7, aes(rn4, logFC, group = condition)) + geom_line(aes(colour = condition)) + xlab("9583 shRNAs") + ylab("log2FC") + theme_classic() + ggtitle("t7 vs. t0") + geom_abline(intercept=0, slope=0) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank(), legend.title=element_blank()) + coord_cartesian(ylim = c(-5, 5)) + scale_colour_manual(values=c(rgb(127,255,127, max = 255), "yellow3", rgb(127,127,255, max = 255))) + scale_y_continuous(breaks=c(-5,-3,-1,1,3,5))
ggsave("/lustre/sblab/martin03/repository/20160418_shRNAscreen_katie_darcie/figures/20160527_A375_pool8_t7vst0_logFC.pdf", width = 12/2.54, height = 8/2.54)



```



TODO:

- Find an FDR that minimises the number of DMSO only significant hairpins but still maximises the number of PDS only hairpins for A375 at t14. I have to play with more stringent FDRs (0.01, 0.0001, 0.000001) and see how the overlap looks like.

- edgeR tutorial

- shRNA tutorial

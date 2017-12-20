rm(list=ls())
setwd("~/Downloads/A375_pool8_logFC_and_FDR_tables/")

DMSO_t14 <- read.table("20160527_A375_pool8_t14_DMSOvst0_logFCdecreasing_FDR.txt")
DMSO_t7 <- read.table("20160527_A375_pool8_t7_DMSOvst0_logFCdecreasing_FDR.txt")

PDS_t14 <- read.table("20160527_A375_pool8_t14_PDSvst0_logFCdecreasing_FDR.txt")
PDS_t7 <- read.table("20160527_A375_pool8_t7_PDSvst0_logFCdecreasing_FDR.txt")
PhenDC3_t14 <- read.table("20160527_A375_pool8_t14_PhenDC3vst0_logFCdecreasing_FDR.txt")
PhenDC3_t7 <- read.table("20160527_A375_pool8_t7_PhenDC3vst0_logFCdecreasing_FDR.txt")


t0_counts <- read.csv("20160526_A375_pool8_t0_counts_PDS.txt",stringsAsFactors = F)

t0_counts <- t0_counts[t0_counts$Gene.Identity!="",c(1,12)]
rownames(t0_counts) <- t0_counts$Gene.Identity
t0_counts$average.T0 <- as.numeric(t0_counts$average.T0)

hits_DMSO_t14 <- DMSO_t14[DMSO_t14$FDR < 0.05,]
hits_DMSO_t7 <- DMSO_t7[DMSO_t7$FDR < 0.05,]


PDS_data_t7 <- merge(t0_counts,PDS_t7,by=0)
PDS_data_t14 <- merge(t0_counts,PDS_t14,by=0)

PDS_data_t7$is_DMSO_hit <- PDS_data_t7$Row.names %in% rownames(hits_DMSO_t7)
PDS_data_t14$is_DMSO_hit <- PDS_data_t14$Row.names %in% rownames(hits_DMSO_t14)

PDS_data_t7$hits <- ifelse(PDS_data_t7$FDR<0.05,
                            ifelse(PDS_data_t7$is_DMSO_hit,"PDS and DMSO hit","PDS hit"),"Not a hit")

PDS_data_t14$hits <- ifelse(PDS_data_t14$FDR<0.05,
                            ifelse(PDS_data_t14$is_DMSO_hit,"PDS and DMSO hit","PDS hit"),"Not a hit")
  

PDS_data_t7$timepoint <- "T7"
PDS_data_t14$timepoint <- "T14"
PDS_data <-rbind(PDS_data_t7,PDS_data_t14)
PDS_data$timepoint <- as.factor(PDS_data$timepoint)
PDS_data$timepoint <- factor(PDS_data$timepoint ,levels(PDS_data$timepoint )[c(2,1)])
PDS_data$order <- ifelse(PDS_data$hits=="Not a hit",1,2)

PDS_data <- PDS_data[order(PDS_data$order),]
ggplot(PDS_data,aes(x=log(average.T0),y=logFC,color=hits))+geom_point()+
        scale_colour_manual(name="",values = c("grey","deepskyblue3", "red"))+theme_bw() + ylim(-12,4) +
        facet_wrap(~timepoint)+xlab("log2(T0 counts)")+ylab("log2(TF counts/T0 counts)")+
        ggtitle("PDS") 

ggsave("PDS.png",width=12.7,height=8.51)
?ggsave
PhenDC3_data_t7 <- merge(t0_counts,PhenDC3_t7,by=0)
PhenDC3_data_t14 <- merge(t0_counts,PhenDC3_t14,by=0)

PhenDC3_data_t7$is_DMSO_hit <- PhenDC3_data_t7$Row.names %in% rownames(hits_DMSO_t7)
PhenDC3_data_t14$is_DMSO_hit <- PhenDC3_data_t14$Row.names %in% rownames(hits_DMSO_t14)

PhenDC3_data_t7$hits <- ifelse(PhenDC3_data_t7$FDR<0.05,
                           ifelse(PhenDC3_data_t7$is_DMSO_hit,"PhenDC3 and DMSO hit","PhenDC3 hit"),"Not a hit")

PhenDC3_data_t14$hits <- ifelse(PhenDC3_data_t14$FDR<0.05,
                            ifelse(PhenDC3_data_t14$is_DMSO_hit,"PhenDC3 and DMSO hit","PhenDC3 hit"),"Not a hit")


PhenDC3_data_t7$timepoint <- "T7"
PhenDC3_data_t14$timepoint <- "T14"
PhenDC3_data <-rbind(PhenDC3_data_t7,PhenDC3_data_t14)
PhenDC3_data$timepoint <- as.factor(PhenDC3_data$timepoint)
PhenDC3_data$timepoint <- factor(PhenDC3_data$timepoint ,levels(PhenDC3_data$timepoint )[c(2,1)])
PhenDC3_data$order <- ifelse(PhenDC3_data$hits=="Not a hit",1,2)

PhenDC3_data <- PhenDC3_data[order(PhenDC3_data$order),]
ggplot(PhenDC3_data,aes(x=log(average.T0),y=logFC,color=hits))+geom_point()+
        scale_colour_manual(name="",values = c("grey","deepskyblue3", "red"))+theme_bw() + ylim(-12,4) +
        facet_wrap(~timepoint)+xlab("log2(T0 counts)")+ylab("log2(TF counts/T0 counts)")+
        ggtitle("PhenDC3") 

ggsave("PhenDC3.png",width=12.7,height=8.51)

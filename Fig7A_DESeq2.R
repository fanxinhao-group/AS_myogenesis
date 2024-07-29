setwd("D:\\result\\wangwei\\mouse\\siFXR1\\")

library(tidyverse)
library(DESeq2)
library(dplyr)

mycounts<-read.csv("transcript_count_matrix.csv",header = T)
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]

Run <- c("si_NC_1","si_NC_2","si_NC_3","si_FXR1_2_1","si_FXR1_2_2","si_FXR1_2_3")
condition<- c("control","control","control","treat","treat","treat")
colData <- data.frame(Run,condition)

condition <- factor(c(rep("control",3),rep("treat",3)), levels = c("control","treat"))
condition
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)

res= results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)

write.csv(res,file="DESeq2_transcript_results_NC_vs_siFXR1_2.csv")
table(res$padj<0.05)

diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
diff_gene_deseq2 <- data.frame(diff_gene_deseq2)
diff_gene_deseq2 <- diff_gene_deseq2 %>% mutate(type= if_else(.$log2FoldChange > 0, 'Up', 'Down'))
dim(diff_gene_deseq2)

write.csv(diff_gene_deseq2,file= "DET_NC_vs_siFXR1_2.csv",quote = F,row.names = T)

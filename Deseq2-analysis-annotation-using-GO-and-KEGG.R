#This project involves differential gene analysis of Granulosa cells (PCOS vs Control),Visualisation using Volcano Plot and annotation using GO and KEGG

gc_counts <-read.csv("GSE155489_gc_pcos_counts.csv",row.names=1,check.names=FALSE)
dim(gc_counts)
colnames(gc_counts)
head(gc_counts)

#create a metadata
sample <- c("GC_B7","GC_B8","GC_B15","GC_B16","GC_B13","GC_B14","GC_B2","GC_B30")
condition <- c("PCOS","PCOS","PCOS","PCOS","control","control","control","control")
gc_meta <- data.frame(sample,condition,stringsAsFactors = FALSE)

#check if the meta data matches whats in the gc_column
all(gc_meta$sample %in% colnames(gc_counts))

#Calling DeSeq2
library(DESeq2)

#create an object for Deseq2
ds_gc <-DESeqDataSetFromMatrix(
  countData= gc_counts,
  colData= gc_meta,
  design= ~ condition
)

#accessing DeSeq2
ds_gc                 
head(counts(ds_gc))         #will give the data in countData
colData(ds_gc)              #will show the ColData

#running deseq - Deseq performs Normalisation based on depth, Dispersion estimation of every gene, and Negative Binomial GLM (stores coefficient based on the differences in the two groups)
deseq <- DESeq(ds_gc)

#extracting the results in a table - wald test on the stored coefficients to get Log2foldchange and Padj of the differentially expressed genes
res_gc <- results(deseq,contrast = c("condition","PCOS","control"))
head(res_gc)
res_df <- as.data.frame(res_gc)

#order by significance
res_gc <- res_gc[order(res_gc$padj), ]
head(res_gc)
summary(res_gc)

#subset the genes with padj<0.05 and absolute log2foldchange>1 and extract into a csv file
sig_res <- subset(res_df,padj<0.05 & abs(log2FoldChange)>1)
write.csv(sig_res,"significantPCOSvscontrol.csv")

#volcano plot with marked gene with highest expression in PCOS

library(ggplot2)
library(ggrepel)
res_df$gene <- rownames(res_df)
top_gene <- res_df[which.max(-log10(res_df$padj)), ]
res_df$threshold <- ifelse(!is.na(res_df$padj) & !is.na(res_df$log2FoldChange) & res_df$padj<0.05 & abs(res_df$log2FoldChange) > 1,"Significant","not significant")
plot=ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj),color=res_df$threshold))+
  geom_point(alpha=0.6)+
  geom_text_repel(data=top_gene, aes(label = gene), color="blue")+
  scale_color_manual(values=c("gray","red"))+
  theme_minimal()+
  labs(title="Volcano Plot: PCOS vs Control")+
  theme(plot.title=element_text(hjust=0.5))
plot


#annotation of significant differentially expressed genes

library(clusterProfiler) 
library(org.Hs.eg.db)
sig_genes <- rownames(sig_res)
entrez <- mapIds(org.Hs.eg.db, keys=sig_genes,column="ENTREZID",keytype="SYMBOL",multiVals="first")
entrez <- na.omit(entrez)


#GO
ego <- enrichGO(gene=entrez,OrgDb=org.Hs.eg.db,keyType="ENTREZID",ont="BP",pAdjustMethod = "BH",pvalueCutoff = 0.05, qvalueCutoff = 0.05)
head(ego)

#KEGG
ekegg <- enrichKEGG(gene=entrez,organism="hsa",pvalueCutoff = 0.05)
ekegg_df <- as.data.frame(ekegg)


#GO and KEGG for upregulated genes in PCOS samples
sig_up <- sig_res[sig_res$log2FoldChange>1, ]
sig_up$genes <- rownames(sig_up)
up_genes <- sig_up$genes
entrez_up <- mapIds(org.Hs.eg.db, keys=up_genes, column="ENTREZID", keytype="SYMBOL", multivals="first")
entrez_Up <- na.omit(entrez_up)

ego_up <- enrichGO(gene=entrez_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", pAdjustMethod = "BH",pvalueCutoff = 0.05, qvalueCutoff = 0.05)
head(ego_up)
ego_up <- as.data.frame(ego_up)


kegg_up <- enrichKEGG(gene=entrez_up, organism="hsa", pvalueCutoff = 0.05)
kegg_up <- as.data.frame(kegg_up)


#GO and KEGG for downregulated genes in PCOS samples
sig_down <- sig_res[sig_res$log2FoldChange<1,]
sig_down$gene <- rownames(sig_down)
down_gene <- sig_down$gene
entrez_down <- mapIds(org.Hs.eg.db, keys=down_gene, column="ENTREZID", keytype="SYMBOL", MULTIVALS="first")
entrez_down <- na.omit(entrez_down)

ego_down <- enrichGO( gene =entrez_down, OrgDb =org.Hs.eg.db,keyType ="ENTREZID", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
ego_down <- as.data.frame(ego_down)

kegg_down <- enrichKEGG(gene=entrez_down, organism="hsa", pvalueCutoff = 0.05)
kegg_down <- as.data.frame(kegg_down)

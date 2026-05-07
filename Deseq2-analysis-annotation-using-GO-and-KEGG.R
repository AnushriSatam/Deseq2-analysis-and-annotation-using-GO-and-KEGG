#This project involves differential gene analysis of Granulosa cells (PCOS vs Control),Visualisation using Volcano Plot and Gene Enrichment analysis using GO and KEGG

gc_counts <-read.csv("GSE155489_gc_pcos_counts.csv",row.names=1,check.names=FALSE)
dim(gc_counts)
colnames(gc_counts)
head(gc_counts)

#create a metadata
samples <- c("GC_B7","GC_B8","GC_B15","GC_B16","GC_B13","GC_B14","GC_B2","GC_B30")
gc_meta <- data.frame(
  row.names = samples,
  condition = factor(c("PCOS","PCOS","PCOS","PCOS","control","control","control","control"),
  levels=c("control","PCOS"))
)


#check if the meta data matches whats in the gc_column
all(rownames(gc_meta) %in% colnames(gc_counts))

#Calling libraries
library(DESeq2)

#create an object for Deseq2
ds_gc <-DESeqDataSetFromMatrix(
  countData= gc_counts,
  colData= gc_meta,
  design= ~ condition
)

#accessing DeSeq2              
head(counts(ds_gc))         #will give the data in countData
colData(ds_gc)              #will show the ColData


#filtering low count genes (the genes must have atleast 10 counts, for atleast 1 condition.)
min_group_size <- min(table(gc_meta$condition))
keep <- rowSums(counts(ds_gc) >= 10) >= min_group_size
ds_gc <- ds_gc[keep,]

#running deseq - Deseq performs Median of ratio method for normalisation
deseq <- DESeq(ds_gc)

#extracting the results in a table - wald test on the stored coefficients to get Log2foldchange and Padj of the differentially expressed genes
res_gc <- results(deseq,contrast = c("condition","PCOS","control"))
head(res_gc)
res_df <- as.data.frame(res_gc)
write.csv(res_df,"DeSeq2_result.csv")


#subset the genes with padj<0.05 and extract into csv file.
sig_res <- subset(res_df, !is.na(res_df$padj) & padj<0.05)
sig_res <- sig_res[order(sig_res$padj), ] 
write.csv(sig_res,"significantPCOSvscontrol.csv")

#annotating the significant genes
gene_ids <- rownames(sig_res)
library(biomaRt)
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "useast"
)

annotations <- getBM(
  attributes = c(
    "external_gene_name", 
    "description",
    "chromosome_name",
    "gene_biotype"
  ),
  filters = "external_gene_name",
  values = gene_ids,
  mart = mart
)
sig_res$GeneID <- rownames(sig_res)

res_annotated <- merge(
  sig_res,
  annotations,
  by.x = "GeneID",
  by.y = "external_gene_name",
  all.x = TRUE
)

res_annotated <- res_annotated[order(res_annotated$padj), ]
write.csv(res_annotated,"Significant_genes_annotated.csv")

#volcano plot with marked gene with highest expression in PCOS

library(ggplot2)
library(ggrepel)
res_df$gene <- rownames(res_df)
res_df$gene_type <- "Not Significant"
res_df$gene_type[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange >= 1] <- "upregulated"
res_df$gene_type[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "downregulated"
top_gene <- res_df[which.max(-log10(res_df$padj)), ]

plot=ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj),color=gene_type))+
  geom_point(alpha=0.6)+
  scale_color_manual(values=c(
    "upregulated" = "red",
    "downregulated" = "blue",
    "Not Significant" = "gray"))+
 geom_vline(xintercept = c(-1,1), linetype="dashed",color="black",linewidth=0.7) +
 geom_hline(yintercept = -log10(0.05), linetype="dashed",color="black",linewidth=0.7)+ 
 geom_text_repel(data=top_gene, aes(label = gene), color="blue")+
  labs(title="Volcano Plot: PCOS vs Control",x="Log2 Fold Change",y="-Log10 Adjusted P-value")+
  theme(plot.title=element_text(hjust=0.5))
ggsave("Volcano_Plot_PCOS_vs_Control.png",plot=plot,width=8,height=5,dpi=100)


#Gene Enrichment Analysis (GO and KEGG)
library(clusterProfiler) 
library(org.Hs.eg.db)
sig_genes <- rownames(sig_res)
upregulated <- res_df[res_df$gene_type=="upregulated",]
downregulated <- res_df[res_df$gene_type=="downregulated",]
up_genes <- rownames(upregulated)
down_genes <- rownames(downregulated)

#GO (to access Biological Process, Molecular Function, and Cellular Component)
up_go_BP <- enrichGO(gene=up_genes,OrgDb=org.Hs.eg.db,keyType="SYMBOL",ont="BP",pAdjustMethod = "BH",pvalueCutoff = 0.05, qvalueCutoff = 0.05)
up_go_MF <- enrichGO(gene=up_genes,OrgDb=org.Hs.eg.db,keyType="SYMBOL",ont="MF",pAdjustMethod="BH",pvalueCutoff=0.05, qvalueCutoff = 0.05 )
up_go_CC <- enrichGO(gene=up_genes,OrgDb=org.Hs.eg.db,keyType="SYMBOL",ont="CC",pAdjustMethod="BH",pvalueCutoff=0.05, qvalueCutoff = 0.05 )

down_go_BP <- enrichGO(gene=down_genes,OrgDb=org.Hs.eg.db,keyType="SYMBOL",ont="BP",pAdjustMethod ="BH",pvalueCutoff = 0.05, qvalueCutoff = 0.05)
down_go_MF <- enrichGO(gene=down_genes,OrgDb=org.Hs.eg.db,keyType="SYMBOL",ont="MF",pAdjustMethod="BH",pvalueCutoff=0.05, qvalueCutoff = 0.05 )
down_go_CC <- enrichGO(gene=down_genes,OrgDb=org.Hs.eg.db,keyType="SYMBOL",ont="CC",pAdjustMethod="BH",pvalueCutoff=0.05, qvalueCutoff = 0.05 )

#KEGG (to access enriched pathways)
entrez_up <- mapIds(org.Hs.eg.db, keys=up_genes, column="ENTREZID", keytype="SYMBOL", multiVals="first")
entrez_down <- mapIds(org.Hs.eg.db, keys=down_genes, column="ENTREZID", keytype="SYMBOL", multiVals="first")
entrez_up <- entrez_up[!is.na(entrez_up)]
entrez_down <- entrez_down[!is.na(entrez_down)]
kegg_up <- enrichKEGG(gene=entrez_up,organism="hsa",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
kegg_down <- enrichKEGG(gene=entrez_down, organism="hsa", pvalueCutoff = 0.05,qvalueCutoff= 0.05)


#function to make dotplot
save_dotplot <- function(enrich_result,filename,title) {
  p <- dotplot(enrich_result,showCategory=20)+
  ggtitle(title)+
  theme(plot.title=element_text(hjust=0.5))
  ggsave(filename, plot=p, width=8, height=6, dpi=100)
    message("Saved: ", filename)
  } 

# GO plots
save_dotplot(up_go_BP,   "up_GO_BP.png",   "Upregulated - Biological Process")
save_dotplot(up_go_MF,   "up_GO_MF.png",   "Upregulated - Molecular Function")
save_dotplot(up_go_CC,   "up_GO_CC.png",   "Upregulated - Cellular Component")
save_dotplot(down_go_BP, "down_GO_BP.png", "Downregulated - Biological Process")
save_dotplot(down_go_MF, "down_GO_MF.png", "Downregulated - Molecular Function")
save_dotplot(down_go_CC, "down_GO_CC.png", "Downregulated - Cellular Component")

# KEGG plots
save_dotplot(kegg_up,   "kegg_up.png",   "KEGG - Upregulated")
save_dotplot(kegg_down, "kegg_down.png", "KEGG - Downregulated")



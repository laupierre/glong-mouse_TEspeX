library (openxlsx)
library (DESeq2)
library (ggplot2)


## See Github RNA-Seq_mouse/gene_annotation.R
#system ("cp /projects/ncrrbt_share_la/dev_pipe/gencode.vM32.annotation.txt .")

anno <- read.delim ("gencode.vM32.annotation.txt")

anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)



## normal STAR results from the RNA-Seq IIT pipeline

a <- read.xlsx ("star_gene_raw_counts.xlsx")
a <- a[ ,grep ("Geneid|SND", colnames (a))]

a <- merge (a, anno, by.x="Geneid", by.y="gene_id", all.x=TRUE) 

# remove sample GL9
a <- a[ ,grep ("GL9", colnames (a), invert=TRUE)]


a <- a[grep ("miRNA|Mt_tRNA|Mt_rRNA|rRNA|snRNA|snoRNA|scRNA|sRNA|misc_RNA|scaRNA|ribozyme|IG_|TR_", a$gene_type, invert=TRUE), ]
counts <- annot <- a
row.names (counts) <- counts$Geneid
counts <- counts[ ,-1]
counts <- counts[ ,grep ("SND", colnames (counts))]
head (counts)

annot <- annot[ ,c("Geneid", "gene_name", "gene_type", "mgi_id", "external_gene_name", "description")]


samples <- data.frame (matrix (nrow=dim (counts)[2], ncol=2))
colnames (samples) <- c("sample", "condition")
samples$sample <- colnames (counts)
samples$condition <- as.factor (c("GL15", "GL15", "GL15", rep ("GL3", 6), rep ("GL15", 2), rep ("WT", 6)))
sampleTable <- samples

idx <- match (sampleTable$sample, colnames (counts))
sampleTable <- sampleTable[idx, ]
stopifnot (sampleTable$sample == colnames (counts))


## TEspeX counts from the outfile.txt

tesp <- read.delim ("outfile.txt", row.names=1)
colnames (tesp) <- gsub ("_SND.*", "_SND", colnames (tesp))
colnames (tesp) <- gsub ("IIT_TRM_", "", colnames (tesp))

tesp <- tesp[ ,colnames (tesp) %in% colnames (counts)]
stopifnot (colnames (tesp) == colnames (counts))

counts <- rbind (counts, tesp)


## DESeq2 

dds <- DESeqDataSetFromMatrix(countData = round (counts), colData = sampleTable, design = ~ condition)

keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep,]
dds


## GL15 vs WT
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "GL15", "WT"))

res <- merge (data.frame (res), counts (dds), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid", all.x=TRUE)

res$gene_name [is.na (res$gene_name)] <- res$Row.names [is.na (res$gene_name)]
res$external_gene_name [is.na (res$external_gene_name)] <- res$Row.names [is.na (res$external_gene_name)]
res$gene_type[is.na (res$gene_type)] <- paste ("transposon", gsub (".*#", "", res$Row.names [is.na (res$gene_type)]), sep=":")

res <- res[order (res$padj), ]
colnames (res)[1] <- "Geneid"

# Sanity check
res[res$gene_name == "Snca", ] 
write.xlsx (res, "GLONG_15months_vs_WT_with_transposons_2023.xlsx", rowNames=F)











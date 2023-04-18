library (openxlsx)
library (DESeq2)
library (ggplot2)


## See Github RNA-Seq_mouse/gene_annotation.R
system ("cp /projects/ncrrbt_share_la/dev_pipe/gencode.vM32.annotation.txt .")

anno <- read.delim ("gencode.vM32.annotation.txt")

anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)


## normal STAR results from the RNA-Seq IIT pipeline

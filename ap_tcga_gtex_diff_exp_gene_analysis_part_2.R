rm(list=ls())

library(ggplot2)
library(ggplot2)
library(gridExtra)
library(data.table)
library(ggthemes)

library(DESeq2)
#####
setwd("<mention your directory path>")

tcga_colon_file = "tcga_colon_data.txt"
tcga_colon_data = read.table(tcga_colon_file, sep = "\t",header = T,row.names = 1)
dim(tcga_colon_data)

gtex_colon_file = "gtex_colon_data.txt"
gtex_colon_data = read.table(gtex_colon_file, sep = "\t",header = T,row.names = 1)
dim(gtex_colon_data)

### Merge and normalize the data
tcga_gtex_col_data = merge(tcga_colon_data, gtex_colon_data, by = 'row.names', all = TRUE)
rownames(tcga_gtex_col_data) = tcga_gtex_col_data$Row.names
tcga_gtex_col_data$Row.names = NULL
dim(tcga_gtex_col_data)

# Map ensembl ids to gene
map_file = "probeMap-gencode.v23.annotation.gene.probemap"
map_data = read.table(map_file,sep = "\t",header = T, row.names = 1)
dim(map_data)

mp = merge(tcga_gtex_col_data, map_data, by = 'row.names', all = TRUE)
dim(mp)
fn = cbind.data.frame(mp$gene,mp[2:594])
colnames(fn) = c("Gene",colnames(fn)[2:594])
rownames(fn) = make.names(fn$Gene, unique=TRUE)
new_fn = fn[,-1]
dim(new_fn)

# Removed the genes which are less than zero
dds <- new_fn[rowSums(new_fn[],na.rm = T)>0,]
dim(dds)

data_log <- log(dds+1,base = 2)           
dim(data_log)
#data_log = dds

# Assign the tcga and gtex samples to their associated labels as TCGA and GTEX
cond = c(rep("TCGA",287),rep("GTEX",306))
final_dd = cbind.data.frame(cond,t(data_log))
dim(final_dd)


# Selecting the gene of interest from final_dd data - for gene name.
# Change the gene name at line number - 57, 61 and 66 for plotting.

gene_data = cbind.data.frame(final_dd$cond,final_dd$TP53)
dim(gene_data)
colnames(gene_data) = c("Sample_Type","Gene")

pdf(paste0("TP53_plot_tcga_gtex.pdf"),width = 10,height = 10)

### Boxplot for gene of interest from two group like tcga and gtex
ggplot(gene_data,aes(x=Sample_Type,y=Gene,fill=Sample_Type)) +
  geom_boxplot()+
  ggtitle("Gene - TP53") + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
  theme(legend.position="none") +
  theme(axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold")) +
  theme(axis.text.x = element_text(angle = 45,size = 20,hjust = 1, face = "bold")) +
  theme(axis.text.y=element_text(size = 20, face = "bold")) +
  theme(plot.title = element_text(size = 20, face = "bold")) + 
  labs(y="Expression_Values", face = "bold") +
  labs(x = "Sample_Type", face = "bold") 

dev.off()

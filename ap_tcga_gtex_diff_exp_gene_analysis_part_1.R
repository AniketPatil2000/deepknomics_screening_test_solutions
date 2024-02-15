rm(list=ls())


#####
setwd("<mention your directory path>")

# Reading tcga count data
tcga_data_file = "tcga_gene_expected_count"
tcga = read.table(tcga_data_file, sep = "\t",header = TRUE)
dim(tcga)

# Reading gtex count data
gtex_data_file = "gtex_gene_expected_count"
gtex = read.table(gtex_data_file, sep = "\t",header = TRUE)
dim(gtex)

# Reading TCGA and GTEX phenotypic information for mapping TCGA and GTEX label.
pheno_file = "TcgaTargetGTEX_phenotype.txt"
phenotype = read.table(pheno_file,sep = "\t",header = TRUE)
dim(phenotype)
phenotype$sample = gsub("-",".",phenotype$sample)
dim(phenotype)

# Fetching phenoptyic information
tcga_pri_tum_samples = phenotype[grep("Primary Tumor",phenotype$X_sample_type),]
dim(tcga_pri_tum_samples)
tcga_colon_samples = tcga_pri_tum_samples[grep("Colon", tcga_pri_tum_samples$primary.disease.or.tissue),]
dim(tcga_colon_samples)

gtex_pri_tum_samples = phenotype[grep("GTEX",phenotype$X_study),]
dim(gtex_pri_tum_samples)
gtex_colon_samples = gtex_pri_tum_samples[grep("Colon", gtex_pri_tum_samples$X_primary_site),]
dim(gtex_colon_samples)

# Mapping the TCGA and GTEX labels for TCGA and GTEX count samples for COlon tissue.
######
m1 = match(colnames(tcga),tcga_colon_samples$sample)
w1 = which(!is.na(m1))
tcga_colon_data = t(t(tcga)[w1,])
dim(tcga_colon_data)
tcga_col_data = cbind.data.frame(tcga$sample,tcga_colon_data[,-1:-2])
colnames(tcga_col_data) = c("Genes",colnames(tcga_col_data)[-1])
dim(tcga_col_data)

m2 = match(colnames(gtex),gtex_colon_samples$sample)
w2 = which(!is.na(m2))
gtex_colon_data = t(t(gtex)[w2,])
gtex_col_data = cbind.data.frame(gtex$sample,gtex_colon_data[,-1:-2])
colnames(gtex_col_data) = c("Genes",colnames(gtex_col_data)[-1])
dim(gtex_col_data)

# Writing the tcga and gtex count data for colon into the file.
write.table(tcga_col_data, "tcga_colon_data.txt",sep = "\t",quote = F,row.names = F)
write.table(gtex_col_data, "gtex_colon_data.txt",sep = "\t",quote = F,row.names = F)

################################################################
# Identifying stable and variable sites on the EPIC array
# Olivia Grant
###############################################################


################################################################################
# load packages
################################################################################

if(!require("VennDiagram", character.only = TRUE)){
  install.packages("VennDiagram")
}
library(VennDiagram)


if(!require("seqLogo", character.only = TRUE)){
  install.packages("seqLogo")
}
library(seqLogo)

if(!require("clusterProfiler", character.only = TRUE)){
  install.packages("clusterProfiler")
}
library(clusterProfiler)

if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")

library(ggvenn)

if(!require("rtracklayer", character.only = TRUE)){
  install.packages("rtracklayer")
}
library(rtracklayer)

if(!require("GenomicRanges", character.only = TRUE)){
  install.packages("GenomicRanges")
}
library(GenomicRanges)

#### 450k vs EPIC array

stableAt <- subset(annotationTable,Name %in% replicatedResultsStable)
variableAt <- subset(annotationTable,Name %in% replicatedResultsVariable)

stable450 <- length(which(stableAt$Methyl450_Loci==TRUE))
stableEPIC <- length(which(stableAt$Methyl450_Loci != TRUE))
variable450 <- length(which(variableAt$Methyl450_Loci==TRUE))
variableEPIC <- length(which(variableAt$Methyl450_Loci != TRUE))

vsEpicMatrix <- matrix(0, 2, 2)
colnames(vsEpicMatrix) <- c("450K","EPIC")
rownames(vsEpicMatrix) <- c("Highly stable probes","Highly variable probes")

vsEpicMatrix[1,1] <-  length(which(stableAt$Methyl450_Loci==TRUE))
vsEpicMatrix[1,2] <-  length(which(stableAt$Methyl450_Loci != TRUE))
vsEpicMatrix[2,1] <-  length(which(variableAt$Methyl450_Loci==TRUE))
vsEpicMatrix[2,2] <-  length(which(variableAt$Methyl450_Loci != TRUE))

vsEpicMatrix <- t(vsEpicMatrix)
vsEpicMatrix_aves <- vsEpicMatrix

vsEpicMatrix_aves[,1] <- vsEpicMatrix[,1]/nrow(stableAt)
vsEpicMatrix_aves[,2] <- vsEpicMatrix[,2]/nrow(variableAt)

pdf("/home/og16379/VaryingDNAm/plots/stackedBarplotEPIC450.pdf")
par(mar=c(4,4,4,10.5))
par(xpd = NA)
barplot(vsEpicMatrix_aves, col=cbPalette, yaxt = "none",cex.names=0.8)
axis(2, seq(0, 1, 0.1), las=2,labels = paste0(seq(0, 100, 10), "%"))
legend(2.5,0.75,rownames(vsEpicMatrix_aves),fill=cbPalette, bty = "n")
dev.off()



####  now lets look at housekeeping vs developmental genes 
####  read in data and get hk vs dv  
gene_expression_bpkm <- read.csv("/home/og16379/VaryingDNAm/data/57epigenomes.RPKM.csv",header=TRUE)
# set rownames
rownames(gene_expression_bpkm)<-gene_expression_bpkm$gene_id
#remove column
geneExpBpkm<-gene_expression_bpkm[,-1]
#make matrix and set row and column names
gene_expression_bpkm_percentile40 <- matrix(FALSE, ncol=ncol(geneExpBpkm), nrow=nrow(geneExpBpkm))
rownames(gene_expression_bpkm_percentile40) <- rownames(geneExpBpkm)
colnames(gene_expression_bpkm_percentile40) <- colnames(geneExpBpkm)
for(j in 1:ncol(geneExpBpkm)){
threshold <- quantile(geneExpBpkm[,j], prob=c(0.4))
gene_expression_bpkm_percentile40[geneExpBpkm[,j] >= threshold, j] <- TRUE
                }

gene_expression_bpkm_percentile40_TRUE <- apply(gene_expression_bpkm_percentile40,1,sum)

# plot histogram of gene expression
pdf("/home/og16379/VaryingDNAm/plots/histoGeneExp.pdf")
hist(gene_expression_bpkm_percentile40_TRUE)
dev.off()

#now get just HK and DV genes
geneExpPerc<-as.data.frame(gene_expression_bpkm_percentile40_TRUE)
hk_Genes<- rownames(subset(geneExpPerc,gene_expression_bpkm_percentile40_TRUE == "57"))

#negate not in
`%notin%` <- Negate(`%in%`)
dv_Genes<-as.data.frame(rownames(geneExpBpkm))
dv_Genes<-subset( dv_Genes, rownames(geneExpBpkm) %notin% hk_Genes)
dv_Genes<-dv_Genes[,1]


#now have a list of housekeeping and developmental genes
#so now I want to get the genomic coordinates of the housekeeping and developmental genes
hg19Ref <- import("/home/og16379/VaryingDNAm/data/Homo_sapiens.GRCh37.87.chr.gtf")
hk_Genes_Gr<-subset(hg19Ref, hg19Ref@elementMetadata@listData$gene_id %in% hk_Genes)
seqlevelsStyle(hk_Genes_Gr)<-"UCSC"
dv_Genes_Gr<-subset(hg19Ref, hg19Ref@elementMetadata@listData$gene_id %in% dv_Genes)
seqlevelsStyle(dv_Genes_Gr)<-"UCSC"

replicatedResultsVariableMapped <- import("/home/og16379/VaryingDNAm/objects/replicatedResultsVariableMapped.bed")
replicatedResultsStableMapped <- import("/home/og16379/VaryingDNAm/objects/replicatedResultsStableMapped.bed")


#check how many overlap with hk and dv genes in stable vs variable
stableInHK<-subsetByOverlaps(replicatedResultsStableMapped,hk_Genes_Gr)
variableInHK<-subsetByOverlaps(replicatedResultsVariableMapped,hk_Genes_Gr)
stableInDV<-subsetByOverlaps(replicatedResultsStableMapped,dv_Genes_Gr)
variableInDV<-subsetByOverlaps(replicatedResultsVariableMapped,dv_Genes_Gr)

hkVsDv <- matrix(0, 2, 2)
colnames(hkVsDv) <- c("Housekeeping Genes","Developmental Genes")
rownames(hkVsDv) <- c("Highly Variable Probes","Highly Stable Probes")

hkVsDv[1,1] <- length(variableInHK)
hkVsDv[1,2] <- length(variableInDV)
hkVsDv[2,1] <- length(stableInHK)
hkVsDv[2,2] <- length(stableInDV)

results_matrix3 <- t(hkVsDv)
results_matrix_aves3 <- results_matrix3

results_matrix_aves3[,1] <- results_matrix3[,1]/length(replicatedResultsVariable)
results_matrix_aves3[,2] <- results_matrix3[,2]/length(replicatedResultsStable)

pdf("/home/og16379/VaryingDNAm/plots/stackedBarplot3.pdf")
par(mar=c(4,4,12,10.5))
par(xpd = NA)
barplot(results_matrix_aves3, col=cbPalette, yaxt = "none",cex.names=.7)
axis(2, seq(0, 1, 0.1), las=2,labels = paste0(seq(0, 100, 10), "%"))
legend(2.5,0.55,rownames(results_matrix_aves3),fill=cbPalette, bty = "n")
dev.off()


################################################################
# Variable and stable probes
# Olivia Grant
# now working on improvign the function to annotate genes

################################################################
library(dplyr)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(DESeq2)
library(biomaRt)
# load in data
load("/home/og16379/diff_cpg_fm/data/firstDatasetNormalised.RData")
load("/home/og16379/VaryingDNAm/objects/variableAt.RData")
load("/home/og16379/VaryingDNAm/objects/stableAt.RData")

# annotations
data("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
annotationTable <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
annotationTable<-as.data.frame(annotationTable)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# function to add gene names as a column
countMatrix<-read.delim("/home/og16379/diff_cpg_fm/data/GSE120312_Counts_Matrix.txt")
countData <- data.frame(row.names=(countMatrix$Ensembl.Gene.ID),countMatrix[,11:30])
sex<-c("Female","Male","Male","Male","Female","Female","Male","Male","Female","Male","Female","Female","Female","Female",
"Male","Male","Female","Male","Female","Male")
metaData <- data.frame("id"=colnames(countData),"dex"=sex)

dds <- DESeqDataSetFromMatrix(countData=countData,
                              colData=metaData,
                              design=~dex)

addGeneNames<- function(countData){
  countData$ensembl_gene_id <- rownames(countData)
  geneNames <- getBM(filters="ensembl_gene_id",attributes=c("ensembl_gene_id","external_gene_name"),
  values=rownames(countData),mart=mart)
  geneNames  <- left_join(countData,y=geneNames,by="ensembl_gene_id")
  countData$geneName <- geneNames$external_gene_name
  return(countData)
}


# apply function to countdata
countData<- addGeneNames(countData)
# now we want to subset the count data by groups
# e.g. we want all the count data for the genes annotated to variable probes with low methylation values
# and all the genes annotated to the variable probes with high methylation values
# normalise counts and create subsets for plotting

getCvAndMeans <- function(dds,countData){
  normalisedCounts <- vst(dds)
  normalisedCounts <- assay(normalisedCounts)
  normalisedCounts <- as.data.frame(normalisedCounts)
  mean <- rowMeans(normalisedCounts)
  sd<-apply(normalisedCounts,MARGIN=1,FUN=sd)
  cv <- (sd / mean) *100
  mean2 <- rowMeans(countData[,1:20])
  sd2<-apply(countData[,1:20],MARGIN=1,FUN=sd)
  cv2 <- (sd2 / mean2) *100
  normalisedCounts<- as.data.frame(cbind(normalisedCounts,cv,mean,cv2,mean2))
  normalisedCounts$geneName <- countData$geneName
  normalisedCounts$geneName[normalisedCounts$geneName==""] <-  "unknown"
  normalisedCounts$geneName[is.na(normalisedCounts$geneName)] <- "unknown"
  return(normalisedCounts)
}
# write function now to get the cv and means for all of the probes and genes n the arrays
normalisedCounts <- getCvAndMeans(dds,countData)

# mean levels of probe methylation so can seperate into groups
# also assign a new column with the first gene annotation so we can subset the CD nicely
# using str split for this
# also now want to add the features to each of the grouped annotations
assignAnnotations<- function(annotationTable,dasen_autosome){
means <-rowMeans(dasen_autosome)
sd <-apply(dasen_autosome,MARGIN=1,FUN=sd)
meansSd <- data.frame("Name"=names(means),"averageVal"=means,"methylationSd"=sd)
intermediateMethylation <-filter(meansSd,between(averageVal,0.40,0.60))
highMethylation <- filter(meansSd,between(averageVal,0.80,1))
lowMethylation <- filter(meansSd,between(averageVal,0,0.20))
#annotationTable$featureAnnotation <- annotationTable$featureAnnotation[]
annotationTable$metStatus <- "NA"
annotationTable$methStatus = with(annotationTable, ifelse(Name %in% rownames(intermediateMethylation) , methStatus <- "intermediateMethylation",
                  ifelse(Name %in% rownames(lowMethylation), methStatus <- "lowMethylation",
                  ifelse(Name %in% rownames(highMethylation), methStatus <- "highMethylation", methStatus <- "NA"))))
annotationTable$firstGene <- vapply(strsplit(annotationTable$UCSC_RefGene_Name,";"), `[`, 1, FUN.VALUE=character(1))
cpgAnnotations <- get(load("/storage/st05d/Exeter/ASD_YW_OG/cgAnnotationsDataFrame.RData"))
annotationTable$featureAnnotation<- cpgAnnotations$Feature[match(annotationTable$Name, cpgAnnotations$CG)]
annotationTable <- merge(annotationTable,meansSd,by="Name")
return(annotationTable)
}



# apply function to variable and stable annotation tables
variableAt <- assignAnnotations(variableAt,dasen_autosome)
stableAt <- assignAnnotations(stableAt,dasen_autosome)


addGeneValues <- function(annotationTable,normalisedCounts){
annotationTable$firstGene[is.na(annotationTable$firstGene)] = "noAnno"
ids <- match(annotationTable$firstGene, normalisedCounts$geneName)
ids_notNA <- !is.na(ids)
annotationTable$geneCv <- NA
annotationTable$geneCv[ids_notNA] <- normalisedCounts$cv[ids[ids_notNA]]
annotationTable$averageExp <- NA
annotationTable$averageExp[ids_notNA] <- normalisedCounts$mean[ids[ids_notNA]]
annotationTable$geneCvN <- NA
annotationTable$geneCvN[ids_notNA] <- normalisedCounts$cv2[ids[ids_notNA]]
annotationTable$averageExpN <- NA
annotationTable$averageExpN[ids_notNA] <- normalisedCounts$mean2[ids[ids_notNA]]
return(annotationTable)
}


variableAt_ <- addGeneValues(variableAt,normalisedCounts)
stableAt <- addGeneValues(stableAt,normalisedCounts)
epiAllelesImVAtCleaned <- addGeneValues(epiAllelesImVAtCleaned,normalisedCounts)

variableAt_ <- split(variableAt,variableAt$methStatus)
stableAt_ <- split(stableAt,stableAt$methStatus)
epiAllelesImVAtSplit <- split(epiAllelesImVAtCleaned,epiAllelesImVAtCleaned$featureAnnotation)

#############
# get all cor values for feature annotations

addCorValues <- function(data.frame){
corTest <- cor.test(data.frame$methylationSd,data.frame$geneCv,na.action=na.omit)
corValues <- data.frame("pValue"=corTest$p.value[1],"corValue"=corTest$estimate[1])
return(corValues)
}

epiAllelesImVAtCor <- lapply(epiAllelesImVAtSplit,addCorValues)

corMatrix <- matrix(0,8,2)
colnames(corMatrix) <- c("corValue","pvalue")
rownames(corMatrix) <- names(epiAllelesImVAtCor)

for(i in 1:length(epiAllelesImVAtCor)){
corMatrix[i,1] <- as.numeric(epiAllelesImVAtCor[[i]][2])
corMatrix[i,2] <- as.numeric(epiAllelesImVAtCor[[i]][1])
}

# repeat for island annotations 
# first need to add all shores and shelves together

epiAllelesImVAtCleaned$islandAnnotation = with(epiAllelesImVAtCleaned, ifelse(epiAllelesImVAtCleaned$Relation_to_Island == "S_Shore", islandAnnotation <- "Shore",
                  ifelse(epiAllelesImVAtCleaned$Relation_to_Island == "N_Shore" , islandAnnotation <- "Shore",
                  ifelse(epiAllelesImVAtCleaned$Relation_to_Island == "N_Shelf", islandAnnotation <- "Shelf",
                  ifelse(epiAllelesImVAtCleaned$Relation_to_Island == "S_Shelf", islandAnnotation <- "Shelf",
                  ifelse(epiAllelesImVAtCleaned$Relation_to_Island == "OpenSea" , islandAnnotation <- "OpenSea",
                  ifelse(epiAllelesImVAtCleaned$Relation_to_Island == "Island" , islandAnnotation <- "Island", islandAnnotation <- "NA")))))))


addCorValues <- function(data.frame){
corTest <- cor.test(data.frame$methylationSd,data.frame$geneCv,na.action=na.omit)
corValues <- data.frame("pValue"=corTest$p.value[1],"corValue"=corTest$estimate[1])
return(corValues)
}

epiAllelesImVAtIsland <- split(epiAllelesImVAtCleaned,epiAllelesImVAtCleaned$islandAnnotation)

epiAllelesImVAtIslandCor <- lapply(epiAllelesImVAtIsland,addCorValues)

corMatrixIsland<- matrix(0,4,2)
colnames(corMatrixIsland) <- c("corValue","pvalue")
rownames(corMatrixIsland) <- names(epiAllelesImVAtIslandCor)

for(i in 1:length(epiAllelesImVAtIslandCor)){
corMatrixIsland[i,1] <- as.numeric(epiAllelesImVAtIslandCor[[i]][2])
corMatrixIsland[i,2] <- as.numeric(epiAllelesImVAtIslandCor[[i]][1])
}

allCorValues <- as.data.frame(rbind(corMatrix,corMatrixIsland))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pdf("/home/og16379/VaryingDNAm/plots/divergingPlotEpi.pdf")
color <- ifelse(allCorValues$corValue < 0,cbPalette[2], cbPalette[3])
ggplot(allCorValues, aes(x = reorder(rownames(allCorValues),corValue), y = corValue)) +
  geom_bar(stat = "identity",
           show.legend = FALSE,fill=color) +
  ylim(-1,1)+
theme_minimal() # Remove the legend
dev.off()




if(file.exists("/home/og16379/VaryingDNAm/data/matrixAnnotated.Rdata")){
  load("/home/og16379/VaryingDNAm/data/matrixAnnotated.Rdata")
} else{
addTestValues <- function(betaMatrix){
betaMatrix <- as.data.frame(betaMatrix)
dipTestValues <- apply(betaMatrix,1,dip.test)
fun1 <- function(list, n){
         sapply(list, `[[`, n)
   }
betaMatrix$statistic <- as.numeric(fun1(dipTestValues, 1))
betaMatrix$pValue <- as.numeric(fun1(dipTestValues, 2))
quartileValues <- apply(betaMatrix,1,quantile)
betaMatrix$firstQ <- as.numeric(quartileValues[2,])
betaMatrix$thirdQ <- as.numeric(quartileValues[4,])
betaMatrix$logFC <- betaMatrix$thirdQ - betaMatrix$firstQ 
return(betaMatrix)
}
matrixAnnotated <- addTestValues(dasen_autosome)
  save(matrixAnnotated, file="/home/og16379/VaryingDNAm/data/matrixAnnotated.Rdata")
}





###############################################################
####  Generate Figure 1B  ####
###############################################################

# volcano plot for epialleles
variableImMatrix <- subset(matrixAnnotated,rownames(matrixAnnotated) %in% variableAt_$intermediateMethylation$Name)
epiAllelesImV <- subset(variableImMatrix,pValue < 0.05)
epiAllelesImVAt <- subset(variableAt,Name %in% rownames(epiAllelesImV))

epiAllelesImVAtCleaned <- subset(epiAllelesImVAt,Name %!in% rownames(gapSignals))

variableImMatrix <- variableImMatrix[order(variableImMatrix$pValue),]
variableImMatrix$rank[rownames(variableImMatrix) %in% epiAllelesImVAtCleaned$Name] <- "1"
variableImMatrix$rank[785:15766] <- 2:14983
#variableImMatrix$pValue[786:15766] <- 2:15677
variableImMatrix$rank <- as.numeric(variableImMatrix$rank)
variableImMatrix$epiAllele[rownames(variableImMatrix) %in% epiAllelesImVAtCleaned$Name] <- "Yes"
variableImMatrix$epiAllele[rownames(variableImMatrix) %!in% epiAllelesImVAtCleaned$Name] <- "No"
save(variableImMatrix,file="/home/og16379/VaryingDNAm/variableImMatrix.RData")

pdf("/home/og16379/VaryingDNAm/plots/rankZoomNew.pdf")
ggplot(variableImMatrix,aes(x=rank,y=-log10(pValue))) + geom_point(stat='identity', 
aes(col=epiAllele), size=2) + xlim(14893,-2500) +
scale_colour_manual(values =cbPalette) +theme_classic() + 
coord_cartesian(ylim=c(0,10))
dev.off()

load("/storage/st05d/Exeter/ASD_YW_OG/cgAnnotationsDataFrame.RData")
load("/home/og16379/VaryingDNAm/objects/replicatedResultsVariable.RData")
load("/home/og16379/VaryingDNAm/objects/replicatedResultsStable.RData")

`%!in%` <- Negate(`%in%`)

background.Anno <- subset(df,CG %!in% rownames(epiAllelesImVAt))
background.Anno <- as.data.frame(table(background.Anno$Feature))
epiAllelesImVAt.Anno <- subset(df,CG %in% epiAllelesImVAt$Name)
epiAllelesImVAt.Anno <- as.data.frame(table(epiAllelesImVAt.Anno$Feature))

results_matrix <- matrix(0, 2, 8)
colnames(results_matrix) <- c("3' UTR","5' UTR", "Enhancers","Exon","Introns","other","Introns","TE's") 
rownames(results_matrix) <- c("All Autosomes", "Epialleles")

# for loops to fill in numbers of annotations 
for (i in 1:ncol(results_matrix)) {
  results_matrix[1,i] <- background.Anno$Freq[i]
}

for (i in 1:ncol(results_matrix)) {
  results_matrix[2,i] <- epiAllelesImVAt.Anno$Freq[i]
}



# convert for heat map
hmap <- results_matrix
results_matrix <- t(results_matrix)
results_matrix_aves <- results_matrix

results_matrix_aves[,1] <- results_matrix[,1]/sum(results_matrix[,1])
results_matrix_aves[,2] <- results_matrix[,2]/sum(results_matrix[,2])


#set up colour palettes
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
variableColour<-cbbPalette[3]

pdf("/home/og16379/VaryingDNAm/plots/epiAllele.pdf")
par(mar=c(4,4,4,10.5))
par(xpd = NA)
barplot(results_matrix_aves, col=cbbPalette, yaxt = "none",cex.names=.7)
axis(2, seq(0, 1, 0.1), las=2,labels = paste0(seq(0, 100, 10), "%"))
legend(2.75,0.75,rownames(results_matrix_aves),fill=cbbPalette, bty = "n")
dev.off()



hmap_aves <- hmap

hmap_aves[1,] <- hmap[1,]/sum(hmap[1,])
hmap_aves[2,] <- hmap[2,]/sum(hmap[2,])

hmap_comp <- hmap_aves[2,]
hmap_comp <- log2(hmap_comp/hmap_aves[1,])
hmap_comp[hmap_comp > 2] <-2
hmap_comp[hmap_comp < -2] <- -2

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
custom_at <- seq(-2.5, 2.5, by = (5/30))
cols_contrast <- rev(c(colorRampPalette(c(cbbPalette[7], "white", cbbPalette[6]))(30)))


pdf("/home/og16379/VaryingDNAm/plots/heatmapEpiallele.pdf")
  levelplot(t(hmap_comp),
    at = custom_at,
    xlab = "log2 (obs/exp)",ylab="",
    scales=list(x=list(rot=90)),
    col.regions = cols_contrast)
dev.off()

###############################################################
####  Generate Figure 1F ####
###############################################################

data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
annotationTable <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
annotationTable<-as.data.frame(annotationTable)

epiAlleleFeature <-subset(annotationTable,Name %in% epiAllelesImVAt$Name)
epiAlleleFeature<-as.data.frame(table(epiAlleleFeature$Relation_to_Island))

backgroundFeature <- subset(annotationTable,Name %!in% epiAllelesImVAt$Name)
backgroundFeature<-as.data.frame(table(backgroundFeature$Relation_to_Island))

results_matrix2 <- matrix(0, 2, 6)
colnames(results_matrix2) <- c("Island","N_Shelf","N_Shore","OpenSea","S_Shelf","S_Shore")
rownames(results_matrix2) <- c("All Autosomes", "Epialleles")


# for loops to fill in numbers of annotations 
for (i in 1:ncol(results_matrix2)) {
  results_matrix2[1,i] <- backgroundFeature$Freq[i]
}

for (i in 1:ncol(results_matrix2)) {
  results_matrix2[2,i] <- epiAlleleFeature$Freq[i]
}


hmap2 <- results_matrix2
results_matrix2 <- t(results_matrix2)
results_matrix_aves2 <- results_matrix2

results_matrix_aves2[,1] <- results_matrix2[,1]/sum(results_matrix2[,1])
results_matrix_aves2[,2] <- results_matrix2[,2]/sum(results_matrix2[,2])


pdf("/home/og16379/VaryingDNAm/plots/stackedBarplotEpialleles.pdf")
par(mar=c(4,4,4,10.5))
par(xpd = NA)
barplot(results_matrix_aves2, col=cbbPalette, yaxt = "none",cex.names=.7)
axis(2, seq(0, 1, 0.1), las=2,labels = paste0(seq(0, 100, 10), "%"))
legend(4,0.75,rownames(results_matrix_aves2),fill=cbbPalette, bty = "n")
dev.off()


hmap_aves2 <- hmap2

hmap_aves2[1,] <- hmap2[1,]/sum(hmap2[1,])
hmap_aves2[2,] <- hmap2[2,]/sum(hmap2[2,])


hmap_comp2 <- hmap_aves2[2,]
hmap_comp2[1,] <- log2(hmap_comp2[1,]/hmap_aves2[1,])
hmap_comp[hmap_comp > 2] <- 2
hmap_comp[hmap_comp < -2] <- 2

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
custom_at <- seq(-2.5, 2.5, by = (5/30))
cols_contrast <- rev(c(colorRampPalette(c(cbbPalette[7], "white", cbbPalette[6]))(30)))


pdf("/home/og16379/VaryingDNAm/plots/heatmap2EpiAlleles.pdf")
  levelplot(t(hmap_comp2),
    at = custom_at,
    xlab = "log2 (obs/exp)",ylab="",
    scales=list(x=list(rot=90)),
    col.regions = cols_contrast)
dev.off()


# correlation plots and correlation values

pdf("/home/og16379/VaryingDNAm/plots/Enhancers2.pdf")
ggplot(variableAt_$intermediateMethylation[variableAt_$intermediateMethylation$featureAnnotation=="Enhancers",], aes(x=methylationSd,y=geneCvN))+
    geom_point(na.rm=TRUE,alpha=.5)+
geom_smooth(method=lm) + theme_minimal() +xlab("Standard deviation of beta values") +ylab("CV of target gene ")
ggplot(stableAt_$highMethylation[stableAt_$highMethylation$featureAnnotation=="Enhancers",], aes(x=methylationSd,y=geneCvN))+
    geom_point(na.rm=TRUE,alpha=.5)+
geom_smooth(method=lm) + theme_minimal() +xlab("Standard deviation of beta values") +ylab("CV of target gene ")
ggplot(stableAt_$lowMethylation[stableAt_$lowMethylation$featureAnnotation=="Enhancers",], aes(x=methylationSd,y=geneCvN))+
    geom_point(na.rm=TRUE,alpha=.5)+
geom_smooth(method=lm) + theme_minimal() +xlab("Standard deviation of beta values") +ylab("CV of target gene ")
dev.off()



cor(variableAt_$intermediateMethylation$methylationSd[variableAt_$intermediateMethylation$featureAnnotation=="Enhancers"],
variableAt_$intermediateMethylation$geneCv[variableAt_$intermediateMethylation$featureAnnotation=="Enhancers"],
use="complete.obs")

cor(stableAt_$highMethylation$methylationSd[stableAt_$highMethylation$featureAnnotation=="Enhancers"],
stableAt_$highMethylation$geneCv[stableAt_$highMethylation$featureAnnotation=="Enhancers"],
use="complete.obs")

cor(stableAt_$lowMethylation$methylationSd[stableAt_$lowMethylation$featureAnnotation=="Enhancers"],
stableAt_$lowMethylation$geneCv[stableAt_$lowMethylation$featureAnnotation=="Enhancers"],
use="complete.obs")


cor(epiAllelesImVAt$averageVal[epiAllelesImVAt$featureAnnotation=="3' UTR"],
epiAllelesImVAt$averageExp[epiAllelesImVAt$featureAnnotation=="3' UTR"],
use="complete.obs")



pdf("/home/og16379/VaryingDNAm/plots/shelfEpiAlleles.pdf")
ggplot(epiAllelesImVAt[epiAllelesImVAt$Relation_to_Island== c("N_Shelf","S_Shelf"),], aes(x=methylationSd,y=geneCvN))+
    geom_point(na.rm=TRUE,alpha=.5)+
geom_smooth(method=lm) + theme_minimal() +xlab("Standard deviation of beta values") +ylab("CV of target gene ")
dev.off()



epicVsFour50Pie <- function(annotationTable){
four50 <-length(which(annotationTable$Methyl450_Loci=="TRUE"))
epic <-length(which(annotationTable$Methyl450_Loci==""))
counts <- c(four50,epic)
pie(counts,col=cbPalette,labels=c("450k","EPIC array specific",border=cbPalette))
}
pdf("/home/og16379/VaryingDNAm/plots/epicVs450KEpiAllelesNew.pdf")
epicVsFour50Pie(epiAllelesImVAtCleaned)
dev.off()


###############################################################
####  Generate Figure 1F ####
###############################################################

data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
annotationTable <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
annotationTable<-as.data.frame(annotationTable)
backgroundFeature <- subset(annotationTable,Name %!in% epiAllelesImVAt$Name)
backgroundFeature<-as.data.frame(table(backgroundFeature$Relation_to_Island))

epiAlleleFeature <-subset(annotationTable,Name %in% epiAllelesImVAt$Name)
epiAlleleFeature<-as.data.frame(table(epiAlleleFeature$Relation_to_Island))
epiAlleleFeature$mqtlRelationship <- "NA"
epiAlleleFeature$mqtlRelationship = with(epiAlleleFeature, ifelse(Name %in% nonMqtl$Name , mqtlRelationship <- "non mQTL",
                  ifelse(Name %in% transEpi$gene, mqtlRelationship <- "trans mQTL",
                  ifelse(Name %in% cisEpi$gene, mqtlRelationship <- "cis mQTL", mqtlRelationship <- "NA"))))

nonMqtlAnno <- as.data.frame(table(epiAlleleFeature$Relation_to_Island[epiAlleleFeature$mqtlRelationship=="non mQTL"]))
transEpiAnno <- as.data.frame(table(epiAlleleFeature$Relation_to_Island[epiAlleleFeature$mqtlRelationship=="trans mQTL"]))
cisEpiAnno <- as.data.frame(table(epiAlleleFeature$Relation_to_Island[epiAlleleFeature$mqtlRelationship=="cis mQTL"]))


resultsMatrix <- matrix(0, 4, 6)
colnames(resultsMatrix) <- c("Island","N_Shelf","N_Shore","OpenSea","S_Shelf","S_Shore")
rownames(resultsMatrix) <- c("All Autosomes", "non mQTL","trans mQTLs","cis mQTLs")


# for loops to fill in numbers of annotations 
for (i in 1:ncol(resultsMatrix)) {
  resultsMatrix[1,i] <- backgroundFeature$Freq[i]
  resultsMatrix[2,i] <- nonMqtlAnno$Freq[i]
  resultsMatrix[3,i] <- transEpiAnno$Freq[i]
  resultsMatrix[4,i] <- cisEpiAnno$Freq[i]
}


pdf("/home/og16379/VaryingDNAm/plots/stackedBarplotEpialleles.pdf")
par(mar=c(4,4,4,10.5))
par(xpd = NA)
barplot(resultsMatrixAves, col=cbbPalette, yaxt = "none",cex.names=.7)
axis(2, seq(0, 1, 0.1), las=2,labels = paste0(seq(0, 100, 10), "%"))
legend(5,0.75,rownames(resultsMatrixAves),fill=cbbPalette, bty = "n")
dev.off()


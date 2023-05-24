
################################################################
# identifying stable and variable sites on the EPIC array
# Olivia Grant
###############################################################

###############################################################
#### load libraries
###############################################################

if(!require("bigmelon", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

  BiocManager::install("bigmelon")
}
library(bigmelon)

if(!require("rtracklayer", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

  BiocManager::install("rtracklayer")
}
library(rtracklayer)

if(!require("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
}
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)


if(!require("lattice", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

  BiocManager::install("lattice")
}
library(lattice)

if(!require("dplyr", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

  BiocManager::install("dplyr")
}
library(dplyr)

if(!require("missMethyl", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

  BiocManager::install("missMethyl")
}
library(missMethyl)



if(!require("IlluminaHumanMethylationEPICanno.ilm10b2.hg19", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
}
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)


################################################################
# functions
###############################################################

##get variable cpgs
getVariableCpGs <- function(data, proportion, top=TRUE){
  #no longer need to calculate mean because not using CV
  #mean <- rowMeans(data)
  sd <- apply(data,MARGIN=1,FUN=sd)
  #replacing cv with sd because cv biased towards unmethylated probes
  ##CV <- (sd/mean)*100
  x <- cbind(data,sd)
  x <- x[order(x$sd, decreasing=TRUE),]
  numberToSelect <- round(proportion*nrow(x))
  if(top){
    toSelect <- 1:numberToSelect
  } else{
    toSelect <- (nrow(x)-numberToSelect + 1):nrow(x)
  }
  return(rownames(x[toSelect,]))
}

#create function to downsample the data
#beta matrix and e.g. 0.8 as downsample proportion, assign seperately
#downsample sample number using sample
# random 235 samples
downsampleData <- function(data, downsampleProportion){
  numberToSelect <- round((proportion)*ncol(data))
  return(data[,sample(1:ncol(data), numberToSelect, replace = FALSE)])
}



getRobustVariableCpGs <- function(data, proportion, top = TRUE, downsampleProportion, downsampleTimes = 10, foundIn = 1.0){
  variableCpGs <- matrix("", nrow=round(proportion*nrow(data)), ncol=downsampleTimes)
  for(i in 1:downsampleTimes){
    downsampledData <- downsampleData(data, downsampleProportion)
    variableCpGs[,i] <- getVariableCpGs(downsampledData, proportion, top)
}
variableCpGsCount <- table(as.vector(variableCpGs))
robustVariableCpGs <- names(variableCpGsCount[variableCpGsCount >= floor(foundIn*downsampleTimes)])
return(robustVariableCpGs)
}

proportion <- 0.10
downsampleProportion <- 0.90
downsampleTimes <- 10
foundIn <- 1.0


################################################################
# calculating sites with discovery data set
###############################################################


# loading and cleaning data
# open gfile to remove outliers
gfile <- openfn.gds('/storage/st05d/Exeter/UnderSocMeth/USM1/USM_WF.gds',
                      readonly=T)

# load interpolatedXY normalised data
# https://www.biorxiv.org/content/biorxiv/early/2021/10/01/2021.09.30.462546.full.pdf
load("/home/og16379/diff_cpg_fm/data/firstDatasetNormalised.RData")

# get sex annotation
info <- pData(gfile)[, c('barcode', 'nsex')]
# change class
info$nsex <- gsub('2', 'F', info$nsex)
info$nsex <- gsub('1', 'M', info$nsex)

# remove samples identified as outliers using the DNAme based sex classifer (Wang et al 2020)
# https://www.biorxiv.org/content/10.1101/2020.10.19.345090v1
# assign vector of outlier barcodes

outlier <- c('200611820013_R08C01', '200603220075_R08C01', '200864580018_R02C01', '200611820020_R07C01')

# remove misclassified barcodes from data
info <- info[!(info$barcode %in% outlier), ]
dasen_autosome <- dasen_autosome[, info$barcode]
dasen_autosomeDf <- as.data.frame(dasen_autosome)

################################################################
# NOW, working with the validation data set
###############################################################

# loading and cleaning data
# open gfile to remove outliers
gfile2 <- openfn.gds("/storage/st05d/Exeter/UnderSocMeth/ISER/USM2_WF.gds")

# load interpolatedXY normalised data
# https://www.biorxiv.org/content/biorxiv/early/2021/10/01/2021.09.30.462546.full.pdf
load("/home/og16379/diff_cpg_fm/data/dasen_autosome2round.Rdata")

## get sex annotation
info2 <- pData(gfile2)[, c('barcode','Sex')]
info2$Sex <- gsub('2', 'F', info2$Sex)
info2$Sex <- gsub('1', 'M', info2$Sex)

# remove samples identified as outliers using the DNAme based sex classifer (Wang et al 2020)
# https://www.biorxiv.org/content/10.1101/2020.10.19.345090v1
outlier2 <- c('203991410050_R06C01', '203991410050_R07C01', '204026590012_R07C01', '203998240051_R07C01',
    '204026590078_R04C01', '203994670097_R08C01', '203960330159_R08C01', '204022160070_R03C01',
    '203994670097_R02C01', '203991470090_R05C01', '204026590012_R06C01')

info2 <- info2[!(info2$barcode %in% outlier2), ]
dasen_autosome2 <- dasen_autosome2[, info2$barcode]
dasen_autosome2Df <- as.data.frame(dasen_autosome2)

################################################################
# now cleaning the data
###############################################################

data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
annotationTable <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
annotationTable<-as.data.frame(annotationTable)
annotationTableOrd <- annotationTable[order(annotationTable$chr,annotationTable$pos),]

#removing sex probes
#get probes on x and y chromosome
sexProbes <- annotationTable$Name[annotationTable$chr %in% c("chrX","chrY")]
sexProbes <- as.character(sexProbes)
#remove probes on x and y chromosomes
CpGs_autosomes <- dasen_autosomeDf[!rownames(dasen_autosomeDf) %in% sexProbes,]
SNPAnnotations<-read.delim("/home/og16379/VaryingDNAm/data/EPIC.hg19.manifest.tsv")
SNPAssociatedCpG <- SNPAnnotations %>% filter(MASK_general==TRUE)
SNPAssociatedCpG <- SNPAssociatedCpG[SNPAssociatedCpG$probeID %in% rownames(CpGs_autosomes),]
CpGs_cleaned <- CpGs_autosomes[!rownames(CpGs_autosomes) %in% SNPAssociatedCpG$probeID,]
save(CpGs_cleaned, file = "/home/og16379/VaryingDNAm/objects/CpGs_cleanedDiscovery.RData")


################################################################
# calculate and save the results for the discovery data
###############################################################

highlyVariableRobustCpGsCombined<- getRobustVariableCpGs(CpGs_cleaned, proportion, top = TRUE, downsampleProportion, downsampleTimes, foundIn)
highlyStableRobustCpGsCombined <- getRobustVariableCpGs(CpGs_cleaned, proportion, top = FALSE, downsampleProportion, downsampleTimes, foundIn)
save(highlyVariableRobustCpGsCombined,file="/home/og16379/VaryingDNAm/objects/highlyVariableRobustCpGsDiscovery.Rdata")
save(highlyStableRobustCpGsCombined,file="/home/og16379/VaryingDNAm/objects/highlyStableRobustCpGsDiscovery.Rdata")



################################################################
# calculating sites with validation data set
###############################################################

################################################################
# now cleaning the data
###############################################################

#removing sex probes
CpGs_autosomes2 <- dasen_autosome2Df[!rownames(dasen_autosome2Df) %in% sexProbes,]
SNPAssociatedCpG <- SNPAssociatedCpG[SNPAssociatedCpG$probeID %in% rownames(CpGs_autosomes2),]
CpGs_cleaned2 <- CpGs_autosomes2[!rownames(CpGs_autosomes2) %in% SNPAssociatedCpG$probeID,]
save(CpGs_cleaned2, file = "/home/og16379/VaryingDNAm/objects/CpGs_cleanedValidation.RData")


################################################################
# calculate and save the results for the validation data
###############################################################

highlyVariableRobustCpGsCombined2<- getRobustVariableCpGs(CpGs_cleaned2, proportion, top = TRUE, downsampleProportion, downsampleTimes, foundIn)
highlyStableRobustCpGsCombined2 <- getRobustVariableCpGs(CpGs_cleaned2, proportion, top = FALSE, downsampleProportion, downsampleTimes, foundIn)
save(highlyVariableRobustCpGsCombined2,file="/home/og16379/VaryingDNAm/objects/highlyVariableRobustCpGsValidation.Rdata")
save(highlyStableRobustCpGsCombined2,file="/home/og16379/VaryingDNAm/objects/highlyStable RobustCpGsValidation.Rdata")


################################################################
# see which sites replicate
###############################################################

replicatedResultsStable <- which(highlyStableRobustCpGsCombined2 %in% highlyStableRobustCpGsCombined)
replicatedResultsStable<-highlyStableRobustCpGsCombined2[replicatedResultsStable]

# gap hunter to remove probes = snp 
gfile <- openfn.gds("/storage/st05d/Exeter/UnderSocMeth/USM1/USM_WF.gds")
mSet<-gds2mset(gfile,anno="epic")
mSet<-mapToGenome(mSet,mergeManifest=TRUE)
# get granges
mSet<-granges(mSet)
gHunter <- gaphunter(mSet)
gapSignals <- gHunter$sampleresults

# none of the gap probes are in the stable data set
replicatedResultsVariable <- replicatedResultsVariable[-which(replicatedResultsVariable %in% rownames(gapSignals))]

# save results
write.csv(replicatedResultsStable, file="/home/og16379/VaryingDNAm/csv/stableAdditionalFile.csv")
save(replicatedResultsStable,file="/home/og16379/VaryingDNAm/objects/replicatedResultsStable.RData")

# save results
write.csv(replicatedResultsVariable, file="/home/og16379/VaryingDNAm/csv/variableAdditionalFile.csv")
save(replicatedResultsVariable,file="/home/og16379/VaryingDNAm/objects/replicatedResultsVariable.RData")

#pathway analysis using gometh. accounts for cpg bias. make sure prior prob = t
#get all cpg names
allCpGProbeNames<-annotationTable$Name
#go meth with GO
highlyVariableRobustCpGsPathway_GO<-gometh(replicatedResultsVariable,allCpGProbeNames,collection="GO",array.type="EPIC",prior.prob=T)
#any with FDR below 0.05
highlyVariableRobustCpGsPathway_GO <- subset(highlyVariableRobustCpGsPathway_GO,FDR < 0.05)
#go meth with GO
highlyStableRobustCpGsPathway_GO<-gometh(replicatedResultsStable,allCpGProbeNames,collection="GO",array.type="EPIC",prior.prob=T)
#any with FDR below 0.05
highlyStableRobustCpGsPathway_GO <- subset(highlyStableRobustCpGsPathway_GO,FDR < 0.05)

#go meth with KEGG
highlyVariableRobustCpGsPathway_KEGG<-gometh(replicatedResultsVariable,allCpGProbeNames,collection="KEGG",array.type="EPIC",prior.prob=T)
#any with FDR below 0.05
highlyVariableRobustCpGsPathway_KEGG <- subset(highlyVariableRobustCpGsPathway_KEGG,FDR < 0.05)
highlyVariableRobustCpGsPathway_KEGG <- highlyVariableRobustCpGsPathway_KEGG[order(highlyVariableRobustCpGsPathway_KEGG$FDR),]
#go meth with KEGG
highlyStableRobustCpGsPathway_KEGG<-gometh(replicatedResultsStable,allCpGProbeNames,collection="KEGG",array.type="EPIC",prior.prob=T)
#any with FDR below 0.05
highlyStableRobustCpGsPathway_KEGG <- subset(highlyStableRobustCpGsPathway_KEGG,FDR < 0.05)
highlyStableRobustCpGsPathway_KEGG <- highlyStableRobustCpGsPathway_KEGG[order(highlyStableRobustCpGsPathway_KEGG$FDR),]

# for loop to write csv
enrichmentFiles <- c("highlyStableRobustCpGsPathway_GO","highlyVariableRobustCpGsPathway_GO","highlyStableRobustCpGsPathway_KEGG","highlyVariableRobustCpGsPathway_KEGG")
for(i in 1:length(enrichmentFiles)) {
  write.csv2(get(enrichmentFiles[i]),
             paste0("/home/og16379/VaryingDNAm/csv/",
                    enrichmentFiles[i],
                    ".csv"),
             row.names = FALSE)
}

################################################################
# create data frames containing the sd values in case need them later
###############################################################

# get sd of the replicated probes
  getSD <- function(x,betaMatrix){
    x<- betaMatrix[rownames(betaMatrix) %in% x,]
    sd<-apply(x,MARGIN=1,FUN=sd)
    x<- as.data.frame(cbind(x,sd))
    x<- x[order(x$sd,decreasing=TRUE),]
    x <- x[startsWith(rownames(x),"cg"),]
    return(x)}


stableSd <- getSD(replicatedResultsStable,dasen_autosome)
variableSd <- getSD(replicatedResultsVariable,dasen_autosome)

sd <- apply(dasen_autosome,MARGIN=1,FUN=sd)
dasenAutosomeSd <- as.data.frame(cbind(dasen_autosome,sd))

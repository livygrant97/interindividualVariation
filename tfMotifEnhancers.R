################################################################
# TF motif enrichment for VMPs and SMPs at enhancers
###############################################################

#below is the code to perform the tf motif enrichment

################################################################################
# libraries
################################################################################

if(!require("MotifDb", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("MotifDb")
}
library(MotifDb)


if(!require("PWMEnrich", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("PWMEnrich")
}
library(PWMEnrich)



if(!require("PWMEnrich.Hsapiens.background", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("PWMEnrich.Hsapiens.background")
}
library(PWMEnrich.Hsapiens.background)


if(!require("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
}
library(BSgenome.Hsapiens.UCSC.hg19)

if(!require("bigmelon", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("bigmelon")
}
library(bigmelon)


if(!require("dplyr", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("dplyr")
}
library(dplyr)


################################################################################
# get data
################################################################################

# load in data


################################################################################
# get data
################################################################################
#convert to methylset function
gds2mset <- function(gds, i, j, anno = NULL){
     x <- gds
     if(!is.null(anno)){
         if(!anno %in% c("27k", "450k", "epic")){
         stop("anno needs to be: \'27k\', \'450k\', \'epic\'")
         }
     }
     M <- x[i = i, j = j,   "methylated", name = TRUE, drop = FALSE]
     U <- x[i = i, j = j, "unmethylated", name = TRUE, drop = FALSE]
     pd <- pData(x)[j, , drop = FALSE]
     rownames(pd) <- colnames(x)[j]
     #    pd <- annotatedDataFrameFrom(object = as.matrix(pd), byrow = TRUE)
     if(!is.null(anno)){
         if(anno == "27k"){
             anno <- c("IlluminaHumanMethylation27k", "ilmn12.hg19")
         } else if(anno == "450k"){
             anno <- c("IlluminaHumanMethylation450k", "ilmn12.hg19")
         } else if(anno == "epic"){
             anno <- c("IlluminaHumanMethylationEPIC", "ilm10b4.hg19")
         } else if(anno == "unknown"){
             anno <- c("Unknown", "Unknown")
         }
     }
     # Guess Array Type - will not get correct array if performed on subset.
     if(is.null(anno)){
         nr <- nrow(fData(x))
         # Will guess array type based on number of rows, will fail on subsets!
         if(nr > 50000 & nr < 500000){
             anno <- c("IlluminaHumanMethylation450k", "ilmn12.hg19")
         } else if(nr >= 500000){
             anno <- c("IlluminaHumanMethylationEPIC", "ilm10b4.hg19")
         } else if(nr <=50000){
             anno <- c("IlluminaHumanMethylation27k", "ilmn12.hg19")
         }
     }
     names(anno) <- c("array", "annotation")
     out <- MethylSet(Meth = M, Unmeth = U, colData = pd, annotation = anno)
     out@preprocessMethod <- c(
         rg.norm = "Converted from gdsfmt to MethylSet (bigmelon)",
         minfi = as.character(packageVersion("minfi")),
         manifest = NA #packageVersion(getManifest(anno))
         )
     out
 }

#this function will take a gds file and convert it into a function so that we can map it to the genome
# this will be usefl then so that later we can annotate each probe t the


# convert gds file to methylumi set file
gfile <- openfn.gds("/storage/st05d/Exeter/UnderSocMeth/USM1/USM_WF.gds")
mSet<-gds2mset(gfile,anno="epic")
mSet<-mapToGenome(mSet,mergeManifest=TRUE)
# get granges
mSet<-granges(mSet)

load("/home/og16379/VaryingDNAm/objects/variableAtFeature.RData")
load("/home/og16379/VaryingDNAm/objects/stableAtFeature.RData")


data(PWMLogn.hg19.MotifDb.Hsap)
# register cores
registerCoresPWMEnrich(30)

variableEnhancersMapped<-subset(mSet,mSet@ranges@NAMES %in% variableAtFeature$Enhancers$Name)
stableEnhancersMapped<-subset(mSet,mSet@ranges@NAMES %in% stableAtFeature$Enhancers$Name)



# extend the GRanges by 100bp
variableEnhancersGRanges10<-resize(variableEnhancersMapped,fix="start",width=width(variableEnhancersMapped)+50)
stableEnhancersGRanges10<-resize(stableEnhancersMapped,fix="start",width=width(stableEnhancersMapped)+50)

# make dna string sets
variableEnhancerDMPseq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, variableEnhancersGRanges10)
stableEnhancersDMPseq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, stableEnhancersGRanges10)

# get PWMS
hg19_PWMs <- PWMLogn.hg19.MotifDb.Hsap$pwms


################################################################################
# build backgrounds
################################################################################
if(file.exists("/home/og16379/VaryingDNAm/data/stableEnhancers_background10.RData")){
  load("/home/og16379/VaryingDNAm/data/stableEnhancers_background10.RData")
} else{
  stableEnhancers_seq_bg <- makeBackground(motifs = PWMLogn.hg19.MotifDb.Hsap$pwms, bg.seq=stableEnhancersDMPseq,
                                              type="logn", algorithm="human",  verbose=TRUE)

  save(stableEnhancers_seq_bg, file="/home/og16379/VaryingDNAm/data/stableEnhancers_background10.RData")
}

if(file.exists("/home/og16379/VaryingDNAm/data/variableEnhancers_background10.RData")){
  load("/home/og16379/VaryingDNAm/data/variableEnhancers_background10.RData")
} else{
  variableEnhancers_seq_bg <- makeBackground(motifs = PWMLogn.hg19.MotifDb.Hsap$pwms, bg.seq=variableEnhancerDMPseq,
                                              type="logn", algorithm="human",verbose=TRUE)

  save(variableEnhancers_seq_bg, file="/home/og16379/VaryingDNAm/data/variableEnhancers_background10.RData")
}

################################################################################
# comparison
################################################################################

# motifs differentially enriched in the maintained sequence (with lognormal background correction)
if(file.exists("/home/og16379/VaryingDNAm/data/variableEnhancers_stableEnhancers_enrichment10.RData")){
  load("/home/og16379/VaryingDNAm/data/variableEnhancers_stableEnhancers_enrichment10.RData")
} else{
variableEnhancers_stableEnhancers_enrichment <- motifEnrichment(variableEnhancerDMPseq, stableEnhancers_seq_bg)
  save(variableEnhancers_stableEnhancers_enrichment, file="/home/og16379/VaryingDNAm/data/variableEnhancers_stableEnhancers_enrichment10.RData")
}

if(file.exists("/home/og16379/VaryingDNAm/data/stableEnhancers_variableEnhancers_enrichment10.RData")){
  load("/home/og16379/VaryingDNAm/data/stableEnhancers_variableEnhancers_enrichment10.RData")
} else{
  stableEnhancers_variableEnhancers_enrichment <- motifEnrichment(stableEnhancersDMPseq, variableEnhancers_seq_bg)
  save(stableEnhancers_variableEnhancers_enrichment, file="/home/og16379/VaryingDNAm/data/stableEnhancers_variableEnhancers_enrichment10.RData")
}




# group report
getEnrMotifs <- function(enrichmentData){
  enrichmentDataReport <- groupReport(enrichmentData)
  enrichmentDataReport <- enrichmentDataReport[enrichmentDataReport$p.value < 0.001]
  enrichmentDataReport <- enrichmentDataReport[-grep("UW.Motif.",enrichmentDataReport$target)]
  enrichedMotifs <- unique(enrichmentDataReport$target)
  enrichmentDataReport <-as.data.frame(enrichmentDataReport)
  return(enrichedMotifs)
}

variableEnhancersSpecificMotifs <- getEnrMotifs(variableEnhancers_stableEnhancers_enrichment)
stableEnhancersSpecificMotifs <- getEnrMotifs(stableEnhancers_variableEnhancers_enrichment)


# get list of motifs enriched at sites hypermethylated in variableEnhancerss only
variableEnhancersSpecific<-variableEnhancersSpecificMotifs[-which(variableEnhancersSpecificMotifs %in% stableEnhancersSpecificMotifs)]
# get list of motifs enriched at sites hypermethylated in stableEnhancerss only
stableEnhancersSpecific<-stableEnhancersSpecificMotifs[-which(stableEnhancersSpecificMotifs %in% variableEnhancersSpecificMotifs)]
# get list of motifs enriched at sites hypermethylated in variableEnhancerss and stableEnhancerss
overlap<-variableEnhancersSpecificMotifs[which(variableEnhancersSpecificMotifs %in% stableEnhancersSpecificMotifs)]

write.csv(variableEnhancersSpecific,file="/home/og16379/VaryingDNAm/csv/variableEnhancersSpecific.csv")
write.csv(stableEnhancersSpecific,file="/home/og16379/VaryingDNAm/csv/stableEnhancersSpecific.csv")


getMotifStats <- function(enrichmentData){
  enrichmentDataReport <- groupReport(enrichmentData)
  enrichmentDataReport <- enrichmentDataReport[enrichmentDataReport$p.value < 0.05]
  enrichmentDataReport <- enrichmentDataReport[-grep("UW.Motif.",enrichmentDataReport$target)]
  enrichmentDataReport <-as.data.frame(enrichmentDataReport)
  enrichmentDataReport <- distinct(enrichmentDataReport,target,.keep_all=TRUE)
  return(enrichmentDataReport)
}

variableEnhancersSpecificMotifStats<- getMotifStats(variableEnhancers_stableEnhancers_enrichment)
stableEnhancersSpecificMotifStats <- getMotifStats(stableEnhancers_variableEnhancers_enrichment)
write.csv(variableEnhancersSpecificMotifStats,file="/home/og16379/VaryingDNAm/csv/variableEnhancersSpecificMotifStats.csv")
write.csv(stableEnhancersSpecificMotifStats,file="/home/og16379/VaryingDNAm/csv/stableEnhancersSpecificMotifStats.csv")


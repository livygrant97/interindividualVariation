################################################################
# tf motif enrichment for VMPs and SMPs at promoters
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

variablePromotersMapped<-subset(mSet,mSet@ranges@NAMES %in% variableAtFeature$Promoters$Name)
stablePromotersMapped<-subset(mSet,mSet@ranges@NAMES %in% stableAtFeature$Promoters$Name)



# extend the GRanges by 100bp
variablePromotersGRanges10<-resize(variablePromotersMapped,fix="start",width=width(variablePromotersMapped)+50)
stablePromotersGRanges10<-resize(stablePromotersMapped,fix="start",width=width(stablePromotersMapped)+50)

# make dna string sets
variableEnhancerDMPseq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, variablePromotersGRanges10)
stablePromotersDMPseq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, stablePromotersGRanges10)

# get PWMS
hg19_PWMs <- PWMLogn.hg19.MotifDb.Hsap$pwms


################################################################################
# build backgrounds
################################################################################
if(file.exists("/home/og16379/VaryingDNAm/data/stablePromoters_background10.RData")){
  load("/home/og16379/VaryingDNAm/data/stablePromoters_background10.RData")
} else{
  stablePromoters_seq_bg <- makeBackground(motifs = PWMLogn.hg19.MotifDb.Hsap$pwms, bg.seq=stablePromotersDMPseq,
                                              type="logn", algorithm="human",  verbose=TRUE)

  save(stablePromoters_seq_bg, file="/home/og16379/VaryingDNAm/data/stablePromoters_background10.RData")
}

if(file.exists("/home/og16379/VaryingDNAm/data/variablePromoters_background10.RData")){
  load("/home/og16379/VaryingDNAm/data/variablePromoters_background10.RData")
} else{
  variablePromoters_seq_bg <- makeBackground(motifs = PWMLogn.hg19.MotifDb.Hsap$pwms, bg.seq=variableEnhancerDMPseq,
                                              type="logn", algorithm="human",verbose=TRUE)

  save(variablePromoters_seq_bg, file="/home/og16379/VaryingDNAm/data/variablePromoters_background10.RData")
}

################################################################################
# comparison
################################################################################

# motifs differentially enriched in the maintained sequence (with lognormal background correction)
if(file.exists("/home/og16379/VaryingDNAm/data/variablePromoters_stablePromoters_enrichment10.RData")){
  load("/home/og16379/VaryingDNAm/data/variablePromoters_stablePromoters_enrichment10.RData")
} else{
variablePromoters_stablePromoters_enrichment <- motifEnrichment(variableEnhancerDMPseq, stablePromoters_seq_bg)
  save(variablePromoters_stablePromoters_enrichment, file="/home/og16379/VaryingDNAm/data/variablePromoters_stablePromoters_enrichment10.RData")
}

if(file.exists("/home/og16379/VaryingDNAm/data/stablePromoters_variablePromoters_enrichment10.RData")){
  load("/home/og16379/VaryingDNAm/data/stablePromoters_variablePromoters_enrichment10.RData")
} else{
  stablePromoters_variablePromoters_enrichment <- motifEnrichment(stablePromotersDMPseq, variablePromoters_seq_bg)
  save(stablePromoters_variablePromoters_enrichment, file="/home/og16379/VaryingDNAm/data/stablePromoters_variablePromoters_enrichment10.RData")
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

variablePromotersSpecificMotifs <- getEnrMotifs(variablePromoters_stablePromoters_enrichment)
stablePromotersSpecificMotifs <- getEnrMotifs(stablePromoters_variablePromoters_enrichment)


# get list of motifs enriched at sites hypermethylated in variablePromoterss only
variablePromotersSpecific<-variablePromotersSpecificMotifs[-which(variablePromotersSpecificMotifs %in% stablePromotersSpecificMotifs)]
# get list of motifs enriched at sites hypermethylated in stablePromoterss only
stablePromotersSpecific<-stablePromotersSpecificMotifs[-which(stablePromotersSpecificMotifs %in% variablePromotersSpecificMotifs)]
# get list of motifs enriched at sites hypermethylated in variablePromoterss and stablePromoterss
overlap<-variablePromotersSpecificMotifs[which(variablePromotersSpecificMotifs %in% stablePromotersSpecificMotifs)]

write.csv(variablePromotersSpecific,file="/home/og16379/VaryingDNAm/csv/variablePromotersSpecific.csv")
write.csv(stablePromotersSpecific,file="/home/og16379/VaryingDNAm/csv/stablePromotersSpecific.csv")


getMotifStats <- function(enrichmentData){
  enrichmentDataReport <- groupReport(enrichmentData)
  enrichmentDataReport <- enrichmentDataReport[enrichmentDataReport$p.value < 0.05]
  enrichmentDataReport <- enrichmentDataReport[-grep("UW.Motif.",enrichmentDataReport$target)]
  enrichmentDataReport <-as.data.frame(enrichmentDataReport)
  enrichmentDataReport <- distinct(enrichmentDataReport,target,.keep_all=TRUE)
  return(enrichmentDataReport)
}

variablePromotersSpecificMotifStats<- getMotifStats(variablePromoters_stablePromoters_enrichment)
stablePromotersSpecificMotifStats <- getMotifStats(stablePromoters_variablePromoters_enrichment)
write.csv(variablePromotersSpecificMotifStats,file="/home/og16379/VaryingDNAm/csv/variablePromotersSpecificMotifStats.csv")
write.csv(stablePromotersSpecificMotifStats,file="/home/og16379/VaryingDNAm/csv/stablePromotersSpecificMotifStats.csv")


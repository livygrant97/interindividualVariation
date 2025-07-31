################################################################
# TF motif enrichment for VMPs and SMPs at promoters 
# Olivia Grant
################################################################

### 1. Libraries
libs <- c("MotifDb", "PWMEnrich", "PWMEnrich.Hsapiens.background",
          "BSgenome.Hsapiens.UCSC.hg19", "bigmelon", "dplyr")

for (pkg in libs) {
  if (!require(pkg, character.only = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}

### 2. File Paths
paths <- list(
  gds = "/storage/st05d/Exeter/UnderSocMeth/USM1/USM_WF.gds",
  variable_feature = "/home/og16379/VaryingDNAm/objects/variableAtFeature.RData",
  stable_feature = "/home/og16379/VaryingDNAm/objects/stableAtFeature.RData",
  bg_variable = "/home/og16379/VaryingDNAm/data/variablePromoters_background10.RData",
  bg_stable = "/home/og16379/VaryingDNAm/data/stablePromoters_background10.RData",
  enr_variable = "/home/og16379/VaryingDNAm/data/variablePromoters_stablePromoters_enrichment10.RData",
  enr_stable = "/home/og16379/VaryingDNAm/data/stablePromoters_variablePromoters_enrichment10.RData",
  csv_variable_specific = "/home/og16379/VaryingDNAm/csv/variablePromotersSpecific.csv",
  csv_stable_specific = "/home/og16379/VaryingDNAm/csv/stablePromotersSpecific.csv",
  csv_variable_stats = "/home/og16379/VaryingDNAm/csv/variablePromotersSpecificMotifStats.csv",
  csv_stable_stats = "/home/og16379/VaryingDNAm/csv/stablePromotersSpecificMotifStats.csv"
)

### 3. Load and map data
gfile <- openfn.gds(paths$gds)
mSet <- gds2mset(gfile, anno="epic")
mSet <- mapToGenome(mSet, mergeManifest=TRUE)
mSet <- granges(mSet)

load(paths$variable_feature)
load(paths$stable_feature)

data(PWMLogn.hg19.MotifDb.Hsap)
registerCoresPWMEnrich(30)

variablePromotersMapped <- subset(mSet, mSet@ranges@NAMES %in% variableAtFeature$Promoters$Name)
stablePromotersMapped <- subset(mSet, mSet@ranges@NAMES %in% stableAtFeature$Promoters$Name)

variablePromotersGRanges10 <- resize(variablePromotersMapped, fix="start", width=width(variablePromotersMapped)+50)
stablePromotersGRanges10 <- resize(stablePromotersMapped, fix="start", width=width(stablePromotersMapped)+50)

variableEnhancerDMPseq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, variablePromotersGRanges10)
stablePromotersDMPseq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, stablePromotersGRanges10)

hg19_PWMs <- PWMLogn.hg19.MotifDb.Hsap$pwms

### 4. Backgrounds
if (file.exists(paths$bg_stable)) {
  load(paths$bg_stable)
} else {
  stablePromoters_seq_bg <- makeBackground(motifs = hg19_PWMs, bg.seq=stablePromotersDMPseq,
                                           type="logn", algorithm="human", verbose=TRUE)
  save(stablePromoters_seq_bg, file=paths$bg_stable)
}

if (file.exists(paths$bg_variable)) {
  load(paths$bg_variable)
} else {
  variablePromoters_seq_bg <- makeBackground(motifs = hg19_PWMs, bg.seq=variableEnhancerDMPseq,
                                             type="logn", algorithm="human", verbose=TRUE)
  save(variablePromoters_seq_bg, file=paths$bg_variable)
}

### 5. Motif enrichment
if (file.exists(paths$enr_variable)) {
  load(paths$enr_variable)
} else {
  variablePromoters_stablePromoters_enrichment <- motifEnrichment(variableEnhancerDMPseq, stablePromoters_seq_bg)
  save(variablePromoters_stablePromoters_enrichment, file=paths$enr_variable)
}

if (file.exists(paths$enr_stable)) {
  load(paths$enr_stable)
} else {
  stablePromoters_variablePromoters_enrichment <- motifEnrichment(stablePromotersDMPseq, variablePromoters_seq_bg)
  save(stablePromoters_variablePromoters_enrichment, file=paths$enr_stable)
}

### 6. Enriched motif lists
getEnrMotifs <- function(enrichmentData){
  enrichmentDataReport <- groupReport(enrichmentData)
  enrichmentDataReport <- enrichmentDataReport[enrichmentDataReport$p.value < 0.001]
  enrichmentDataReport <- enrichmentDataReport[-grep("UW.Motif.", enrichmentDataReport$target)]
  enrichedMotifs <- unique(enrichmentDataReport$target)
  enrichmentDataReport <- as.data.frame(enrichmentDataReport)
  return(enrichedMotifs)
}

variablePromotersSpecificMotifs <- getEnrMotifs(variablePromoters_stablePromoters_enrichment)
stablePromotersSpecificMotifs <- getEnrMotifs(stablePromoters_variablePromoters_enrichment)

variablePromotersSpecific <- setdiff(variablePromotersSpecificMotifs, stablePromotersSpecificMotifs)
stablePromotersSpecific <- setdiff(stablePromotersSpecificMotifs, variablePromotersSpecificMotifs)
overlap <- intersect(variablePromotersSpecificMotifs, stablePromotersSpecificMotifs)

write.csv(variablePromotersSpecific, file=paths$csv_variable_specific)
write.csv(stablePromotersSpecific, file=paths$csv_stable_specific)

### 7. Motif statistics
getMotifStats <- function(enrichmentData){
  enrichmentDataReport <- groupReport(enrichmentData)
  enrichmentDataReport <- enrichmentDataReport[enrichmentDataReport$p.value < 0.05]
  enrichmentDataReport <- enrichmentDataReport[-grep("UW.Motif.", enrichmentDataReport$target)]
  enrichmentDataReport <- as.data.frame(enrichmentDataReport)
  enrichmentDataReport <- distinct(enrichmentDataReport, target, .keep_all=TRUE)
  return(enrichmentDataReport)
}

variablePromotersSpecificMotifStats <- getMotifStats(variablePromoters_stablePromoters_enrichment)
stablePromotersSpecificMotifStats <- getMotifStats(stablePromoters_variablePromoters_enrichment)

write.csv(variablePromotersSpecificMotifStats, file=paths$csv_variable_stats)
write.csv(stablePromotersSpecificMotifStats, file=paths$csv_stable_stats)

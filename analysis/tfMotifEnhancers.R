################################################################
# TF motif enrichment for VMPs and SMPs at enhancers (refactored)
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
  bg_variable = "/home/og16379/VaryingDNAm/data/variableEnhancers_background10.RData",
  bg_stable = "/home/og16379/VaryingDNAm/data/stableEnhancers_background10.RData",
  enr_variable = "/home/og16379/VaryingDNAm/data/variableEnhancers_stableEnhancers_enrichment10.RData",
  enr_stable = "/home/og16379/VaryingDNAm/data/stableEnhancers_variableEnhancers_enrichment10.RData",
  csv_variable_specific = "/home/og16379/VaryingDNAm/csv/variableEnhancersSpecific.csv",
  csv_stable_specific = "/home/og16379/VaryingDNAm/csv/stableEnhancersSpecific.csv",
  csv_variable_stats = "/home/og16379/VaryingDNAm/csv/variableEnhancersSpecificMotifStats.csv",
  csv_stable_stats = "/home/og16379/VaryingDNAm/csv/stableEnhancersSpecificMotifStats.csv"
)

### 3. Load data and create GRanges

gfile <- openfn.gds(paths$gds)
mSet <- gds2mset(gfile, anno="epic")
mSet <- mapToGenome(mSet, mergeManifest=TRUE)
mSet <- granges(mSet)

load(paths$variable_feature)
load(paths$stable_feature)

data(PWMLogn.hg19.MotifDb.Hsap)
registerCoresPWMEnrich(30)

variableEnhancersMapped <- subset(mSet, mSet@ranges@NAMES %in% variableAtFeature$Enhancers$Name)
stableEnhancersMapped <- subset(mSet, mSet@ranges@NAMES %in% stableAtFeature$Enhancers$Name)

variableEnhancersGRanges10 <- resize(variableEnhancersMapped, fix="start", width=width(variableEnhancersMapped)+50)
stableEnhancersGRanges10 <- resize(stableEnhancersMapped, fix="start", width=width(stableEnhancersMapped)+50)

variableEnhancerDMPseq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, variableEnhancersGRanges10)
stableEnhancersDMPseq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, stableEnhancersGRanges10)

### 4. Background creation

if (file.exists(paths$bg_stable)) {
  load(paths$bg_stable)
} else {
  stableEnhancers_seq_bg <- makeBackground(motifs = PWMLogn.hg19.MotifDb.Hsap$pwms, bg.seq=stableEnhancersDMPseq,
                                           type="logn", algorithm="human", verbose=TRUE)
  save(stableEnhancers_seq_bg, file=paths$bg_stable)
}

if (file.exists(paths$bg_variable)) {
  load(paths$bg_variable)
} else {
  variableEnhancers_seq_bg <- makeBackground(motifs = PWMLogn.hg19.MotifDb.Hsap$pwms, bg.seq=variableEnhancerDMPseq,
                                             type="logn", algorithm="human", verbose=TRUE)
  save(variableEnhancers_seq_bg, file=paths$bg_variable)
}

### 5. Motif enrichment

if (file.exists(paths$enr_variable)) {
  load(paths$enr_variable)
} else {
  variableEnhancers_stableEnhancers_enrichment <- motifEnrichment(variableEnhancerDMPseq, stableEnhancers_seq_bg)
  save(variableEnhancers_stableEnhancers_enrichment, file=paths$enr_variable)
}

if (file.exists(paths$enr_stable)) {
  load(paths$enr_stable)
} else {
  stableEnhancers_variableEnhancers_enrichment <- motifEnrichment(stableEnhancersDMPseq, variableEnhancers_seq_bg)
  save(stableEnhancers_variableEnhancers_enrichment, file=paths$enr_stable)
}

### 6. Motif filtering

getEnrMotifs <- function(enrichmentData){
  enrichmentDataReport <- groupReport(enrichmentData)
  enrichmentDataReport <- enrichmentDataReport[enrichmentDataReport$p.value < 0.001]
  enrichmentDataReport <- enrichmentDataReport[-grep("UW.Motif.", enrichmentDataReport$target)]
  enrichedMotifs <- unique(enrichmentDataReport$target)
  enrichmentDataReport <- as.data.frame(enrichmentDataReport)
  return(enrichedMotifs)
}

variableEnhancersSpecificMotifs <- getEnrMotifs(variableEnhancers_stableEnhancers_enrichment)
stableEnhancersSpecificMotifs <- getEnrMotifs(stableEnhancers_variableEnhancers_enrichment)

variableEnhancersSpecific <- setdiff(variableEnhancersSpecificMotifs, stableEnhancersSpecificMotifs)
stableEnhancersSpecific <- setdiff(stableEnhancersSpecificMotifs, variableEnhancersSpecificMotifs)
overlap <- intersect(variableEnhancersSpecificMotifs, stableEnhancersSpecificMotifs)

write.csv(variableEnhancersSpecific, file=paths$csv_variable_specific)
write.csv(stableEnhancersSpecific, file=paths$csv_stable_specific)

### 7. Motif stats

getMotifStats <- function(enrichmentData){
  enrichmentDataReport <- groupReport(enrichmentData)
  enrichmentDataReport <- enrichmentDataReport[enrichmentDataReport$p.value < 0.05]
  enrichmentDataReport <- enrichmentDataReport[-grep("UW.Motif.", enrichmentDataReport$target)]
  enrichmentDataReport <- as.data.frame(enrichmentDataReport)
  enrichmentDataReport <- distinct(enrichmentDataReport, target, .keep_all=TRUE)
  return(enrichmentDataReport)
}

variableEnhancersSpecificMotifStats <- getMotifStats(variableEnhancers_stableEnhancers_enrichment)
stableEnhancersSpecificMotifStats <- getMotifStats(stableEnhancers_variableEnhancers_enrichment)

write.csv(variableEnhancersSpecificMotifStats, file=paths$csv_variable_stats)
write.csv(stableEnhancersSpecificMotifStats, file=paths$csv_stable_stats)

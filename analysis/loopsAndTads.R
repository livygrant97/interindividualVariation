###############################################################
# mQTL-TAD-loop overlap analysis
# Olivia Grant
###############################################################

#### SETUP: Paths, Data, and Libraries

# Define file paths
tad_file <- "/home/og16379/VaryingDNAm/data/hicExplorer2/livHiC_min30000_max60000_step40000_thres0.05_delta0.01_fdr_domains.bed"
mqtl_file <- "/home/og16379/VaryingDNAm/data/AllmQTL_Unrelated_with10PCs_1e-10.rdata"
rep_var_file <- "/home/og16379/VaryingDNAm/objects/replicatedResultsVariable.RData"
rep_stable_file <- "/home/og16379/VaryingDNAm/objects/replicatedResultsStable.RData"
plot_dir <- "/home/og16379/VaryingDNAm/plots"

# Load TADs and convert to GRanges
tads <- read.table(tad_file)
colnames(tads) <- c("chr", "start", "end", "V4", "tadScore", "strand", "start1", "end1", "V9")
tads <- makeGRangesFromDataFrame(tads, keep.extra.columns = TRUE)
seqlevelsStyle(tads) <- "UCSC"
tads$TADID <- seq_along(tads)

# Load mQTL and probe sets
load(mqtl_file)
load(rep_var_file)
load(rep_stable_file)

#### FUNCTIONS

# Create GInteractions for SNP–CpG pairs
getAllSnpCpgPairs <- function(mqtlAnno, probesOfInterest) {
  paired <- subset(mqtlAnno, gene %in% probesOfInterest) %>%
    distinct(gene, .keep_all = TRUE) %>%
    as.data.frame()

  rowGR <- GRanges(paired[,1], IRanges(paired[,2], paired[,2]))
  colGR <- GRanges(paired[,9], IRanges(paired[,10], paired[,10]))
  
  gi <- GInteractions(seq_along(rowGR), length(rowGR) + seq_along(colGR),
                      c(rowGR, colGR), mcols = list(paired[,5:22]), mode = "reverse")
  
  anchors <- anchors(gi)
  gi$distClass <- ifelse(distance(anchors$first, anchors$second) <= 500, "cisMqtl", "transMqtl")
  seqlevelsStyle(gi) <- "UCSC"
  return(gi)
}

# Check if anchors fall outside any TADs
anchorsInNoTads <- function(mqtlAnno, probesOfInterest, tads) {
  anchors <- subset(mqtlAnno, gene %in% probesOfInterest) %>% unique(by = "gene")
  a1 <- GRanges(anchors$chrom, IRanges(anchors$start, anchors$start))
  a2 <- GRanges(anchors$chrom2, IRanges(anchors$start2, anchors$start2))
  seqlevelsStyle(a1) <- "UCSC"
  seqlevelsStyle(a2) <- "UCSC"

  olap1 <- findOverlaps(a1, tads)
  olap2 <- findOverlaps(a2, tads)

  inTAD1 <- integer(length(a1))
  inTAD2 <- integer(length(a2))
  inTAD1[queryHits(olap1)] <- subjectHits(olap1)
  inTAD2[queryHits(olap2)] <- subjectHits(olap2)

  return(inTAD1 == 0 | inTAD2 == 0)
}

# Check if both anchors fall in the same TAD
anchorsInSameTads <- function(mqtlAnno, probesOfInterest, tads) {
  anchors <- subset(mqtlAnno, gene %in% probesOfInterest) %>% unique(by = "gene")
  a1 <- GRanges(anchors$chrom, IRanges(anchors$start, anchors$start))
  a2 <- GRanges(anchors$chrom2, IRanges(anchors$start2, anchors$start2))
  seqlevelsStyle(a1) <- "UCSC"
  seqlevelsStyle(a2) <- "UCSC"

  olap1 <- findOverlaps(a1, tads)
  olap2 <- findOverlaps(a2, tads)

  inTAD1 <- integer(length(a1))
  inTAD2 <- integer(length(a2))
  inTAD1[queryHits(olap1)] <- subjectHits(olap1)
  inTAD2[queryHits(olap2)] <- subjectHits(olap2)

  return(inTAD1 != 0 & inTAD1 == inTAD2)
}

#### APPLY: mQTL Sets
backgroundProbes <- annotationTable$Name
chanceMqtl   <- getAllSnpCpgPairs(mqtl, backgroundProbes)
variableMqtl <- getAllSnpCpgPairs(mqtl, replicatedResultsVariable)
stableMqtl   <- getAllSnpCpgPairs(mqtl, replicatedResultsStable)
epiAlleleMqtl <- getAllSnpCpgPairs(mqtl, epiAllelesImVAtCleaned$Name)

# Annotate overlaps
variableMqtl$inSameTad <- anchorsInSameTads(mqtl, replicatedResultsVariable, tads)
stableMqtl$inSameTad   <- anchorsInSameTads(mqtl, replicatedResultsStable, tads)
chanceMqtl$inSameTad   <- anchorsInSameTads(mqtl, backgroundProbes, tads)

variableMqtl$inNoTad <- anchorsInNoTads(mqtl, replicatedResultsVariable, tads)
stableMqtl$inNoTad   <- anchorsInNoTads(mqtl, replicatedResultsStable, tads)
chanceMqtl$inNoTad   <- anchorsInNoTads(mqtl, backgroundProbes, tads)

#### LOOP CONNECTIVITY (assumes `loopsGi` already exists)

# Check if both anchors fall outside loops
anchorsNotConnectedByLoops <- function(mqtlAnno, probesOfInterest, loops) {
  anchors <- subset(mqtlAnno, gene %in% probesOfInterest) %>% distinct(gene, .keep_all = TRUE)
  a1 <- GRanges(anchors[,1], IRanges(anchors[,2], anchors[,2]))
  a2 <- GRanges(anchors[,9], IRanges(anchors[,10], anchors[,10]))
  seqlevelsStyle(a1) <- "UCSC"
  seqlevelsStyle(a2) <- "UCSC"

  olap1 <- findOverlaps(a1, loops)
  olap2 <- findOverlaps(a2, loops)

  l1 <- integer(length(a1))
  l2 <- integer(length(a2))
  l1[queryHits(olap1)] <- subjectHits(olap1)
  l2[queryHits(olap2)] <- subjectHits(olap2)

  return(l1 == 0 & l2 == 0)
}

# Check if both anchors fall in the same loop
anchorsConnectedByLoops <- function(mqtlAnno, probesOfInterest, loops) {
  anchors <- subset(mqtlAnno, gene %in% probesOfInterest) %>% distinct(gene, .keep_all = TRUE)
  a1 <- GRanges(anchors[,1], IRanges(anchors[,2], anchors[,2]))
  a2 <- GRanges(anchors[,9], IRanges(anchors[,10], anchors[,10]))
  seqlevelsStyle(a1) <- "UCSC"
  seqlevelsStyle(a2) <- "UCSC"

  olap1 <- findOverlaps(a1, loops)
  olap2 <- findOverlaps(a2, loops)

  l1 <- integer(length(a1))
  l2 <- integer(length(a2))
  l1[queryHits(olap1)] <- subjectHits(olap1)
  l2[queryHits(olap2)] <- subjectHits(olap2)

  return(l1 != 0 & l1 == l2)
}

# Apply to sets
stableMqtl$cbLoop   <- anchorsConnectedByLoops(mqtl, replicatedResultsStable, loopsGi)
variableMqtl$cbLoop <- anchorsConnectedByLoops(mqtl, replicatedResultsVariable, loopsGi)
chanceMqtl$cbLoop   <- anchorsConnectedByLoops(mqtl, backgroundProbes, loopsGi)

stableMqtl$noLoop   <- anchorsNotConnectedByLoops(mqtl, replicatedResultsStable, loopsGi)
variableMqtl$noLoop <- anchorsNotConnectedByLoops(mqtl, replicatedResultsVariable, loopsGi)
chanceMqtl$noLoop   <- anchorsNotConnectedByLoops(mqtl, backgroundProbes, loopsGi)

#### ANALYSIS AND PLOTTING

# Count matrix based on overlap type
calculateCounts <- function(df, groupCol1, groupCol2, classCol) {
  m <- matrix(0, 3, 3)
  m[1,1] <- sum(df[[groupCol1]] & df[[classCol]] == "transMqtl")
  m[1,2] <- sum(df[[groupCol1]] & df[[classCol]] == "cisMqtl")
  m[1,3] <- sum(df[[groupCol1]])

  m[2,1] <- sum(!df[[groupCol1]] & !df[[groupCol2]] & df[[classCol]] == "transMqtl")
  m[2,2] <- sum(!df[[groupCol1]] & !df[[groupCol2]] & df[[classCol]] == "cisMqtl")
  m[2,3] <- sum(!df[[groupCol1]] & !df[[groupCol2]])

  m[3,1] <- sum(!df[[groupCol1]] & df[[groupCol2]] & df[[classCol]] == "transMqtl")
  m[3,2] <- sum(!df[[groupCol1]] & df[[groupCol2]] & df[[classCol]] == "cisMqtl")
  m[3,3] <- sum(!df[[groupCol1]] & df[[groupCol2]])
  m
}

# Convert counts to percentages
calculatePercentage <- function(counts) {
  apply(counts, 2, function(x) 100 * x / sum(x))
}

# Save barplot to file
createBarplot <- function(results, title, legendNames) {
  pdf(file.path(plot_dir, paste0("groupedBar", title, ".pdf")))
  par(mar = c(4, 4, 4, 10.5), xpd = NA)
  barplot(results, col = cbbPalette[5:7], border = "white", cex.names = 0.4, ylim = c(0, 100), beside = TRUE)
  legend("topright", legend = legendNames, fill = cbbPalette[5:6], bty = "n")
  dev.off()
}

# Apply barplot pipeline
makeAllPlots <- function(df, suffix) {
  counts <- calculateCounts(df, "cbLoop", "noLoop", "distClass")
  perc <- calculatePercentage(counts)
  createBarplot(perc, paste0("Loops", suffix), rownames(perc))
}

makeAllPlots(variableMqtl, "All")
makeAllPlots(variableMqtl, "Trans")
makeAllPlots(variableMqtl, "cis")

# Plot for epialleles
epiResults <- matrix(0, 1, 3)
epiResults[1, 1] <- sum(epiAllelesImVAtCleaned$Name %!in% mqtl$gene)
epiResults[1, 2] <- sum(epiAlleleMqtl@elementMetadata@listData$distClass == "cisMqtl")
epiResults[1, 3] <- sum(epiAlleleMqtl@elementMetadata@listData$distClass == "transMqtl")
epiPerc <- calculatePercentage(epiResults)
createBarplot(epiPerc, "epiMqtlRel", colnames(epiPerc))

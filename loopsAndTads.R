
# read in tad data 
tadsBloodHicExplorerTwo <- read.table("/home/og16379/VaryingDNAm/data/hicExplorer2/livHiC_min30000_max60000_step40000_thres0.05_delta0.01_fdr_domains.bed")
colnames(tadsBloodHicExplorerTwo) <- c("chr","start","end","V4","tadSeperationScore","strand","start1","end1","V9")
tadsBloodHicExplorerTwo <- makeGRangesFromDataFrame(tadsBloodHicExplorerTwo,keep.extra.columns=TRUE)
seqlevelsStyle(tadsBloodHicExplorerTwo) <- "UCSC"
tadsBloodHicExplorerTwo$TADID <- seq(1,length(tadsBloodHicExplorerTwo))
# read in mqtl and probe data 
mqtl <- get(load("/home/og16379/VaryingDNAm/data/AllmQTL_Unrelated_with10PCs_1e-10.rdata"))
load("/home/og16379/VaryingDNAm/objects/replicatedResultsVariable.RData")
load("/home/og16379/VaryingDNAm/objects/replicatedResultsStable.RData")

# get granges for snp and cpg pairs
getAllSnpCpgPairs <- function(mqtlAnno,probesOfInterest){
pairedData <- subset(mqtlAnno,gene %in% probesOfInterest)
pairedData <- as.data.frame(pairedData %>% distinct(gene,.keep_all=TRUE))
convertToGIIndeces <- function(df, mode="reverse"){
row.regions <- GRanges(df[,1], IRanges(df[,2],df[,2]))# interaction start
col.regions <- GRanges(df[,9], IRanges(df[,10],df[,10]))# interaction end
gi <- GInteractions(1:length(row.regions),(length(row.regions)+1):(length(row.regions)+length(col.regions)), c(row.regions,col.regions), mcols=list(df[,5:22]),mode=mode)
giAnchors <- anchors(gi)
gi$distClass = ifelse(distance(giAnchors$first, giAnchors$second) <=500, distClass  <- "cisMqtl",distClass <- "transMqtl")
return(gi)
}
pairedData <- convertToGIIndeces(pairedData)
seqlevelsStyle(pairedData) <-"UCSC"
return(pairedData)
}


# apply function to both sets
backgroundProbes <- annotationTable$Name
chanceMqtl <- getAllSnpCpgPairs(mqtl,backgroundProbes)
variableMqtl <- getAllSnpCpgPairs(mqtl,replicatedResultsVariable)
stableMqtl <- getAllSnpCpgPairs(mqtl,replicatedResultsStable)
epiAlleleMqtl <- getAllSnpCpgPairs(mqtl,epiAllelesImVAtCleaned$Name)



anchorsInNoTads <- function(mqtlAnno, probesOfInterest, tads) {
  anchorBoth <- mqtlAnno[mqtlAnno$gene %in% probesOfInterest, ]
  anchorBoth <- unique(anchorBoth, by = "gene")
  anchorOne <- GRanges(anchorBoth$chrom, IRanges(anchorBoth$start, anchorBoth$start))
  anchorTwo <- GRanges(anchorBoth$chrom2, IRanges(anchorBoth$start2, anchorBoth$start2))
  seqlevelsStyle(anchorOne) <- "UCSC"
  seqlevelsStyle(anchorTwo) <- "UCSC"
  olapA1 <- findOverlaps(anchorOne, tads)
  olapA2 <- findOverlaps(anchorTwo, tads)
  anchorOneTAD <- rep(0, length(anchorOne))
  anchorTwoTAD <- rep(0, length(anchorTwo))
  anchorOneTAD[queryHits(olapA1)] <- subjectHits(olapA1)
  anchorTwoTAD[queryHits(olapA2)] <- subjectHits(olapA2)
  ids_noTAD <- anchorOneTAD == 0 | anchorTwoTAD == 0
  mqtlsInNoTads <- anchorBoth[ids_noTAD, ]
  return(ids_noTAD)
}


anchorsInSameTads <- function(mqtlAnno, probesOfInterest, tads) {
  anchorBoth <- mqtlAnno[mqtlAnno$gene %in% probesOfInterest, ]
  anchorBoth <- unique(anchorBoth, by = "gene")
  anchorOne <- GRanges(anchorBoth$chrom, IRanges(anchorBoth$start, anchorBoth$start))
  anchorTwo <- GRanges(anchorBoth$chrom2, IRanges(anchorBoth$start2, anchorBoth$start2))
  seqlevelsStyle(anchorOne) <- "UCSC"
  seqlevelsStyle(anchorTwo) <- "UCSC"
  olapA1 <- findOverlaps(anchorOne, tads)
  olapA2 <- findOverlaps(anchorTwo, tads)
  anchorOneTAD <- rep(0, length(anchorOne))
  anchorTwoTAD <- rep(0, length(anchorTwo))
  anchorOneTAD[queryHits(olapA1)] <- subjectHits(olapA1)
  anchorTwoTAD[queryHits(olapA2)] <- subjectHits(olapA2)
  ids_sameTAD <- anchorOneTAD != 0 & anchorOneTAD == anchorTwoTAD
  mqtlsInSameTads <- anchorBoth[ids_sameTAD, ]
  return(ids_sameTAD)
}

variableMqtl$inSameTad = anchorsInSameTads(mqtl,replicatedResultsVariable,tadsBloodHicExplorerTwo)
stableMqtl$inSameTad = anchorsInSameTads(mqtl,replicatedResultsStable,tadsBloodHicExplorerTwo)
chanceMqtl$inSameTad = anchorsInSameTads(mqtl,backgroundProbes,tadsBloodHicExplorerTwo)

variableMqtl$inNoTad = anchorsInNoTads(mqtl,replicatedResultsVariable,tadsBloodHicExplorerTwo)
stableMqtl$inNoTad = anchorsInNoTads(mqtl,replicatedResultsStable,tadsBloodHicExplorerTwo)
chanceMqtl$inNoTad = anchorsInNoTads(mqtl,backgroundProbes,tadsBloodHicExplorerTwo)

# Function to calculate the counts for each category
calculateCounts <- function(df, inSameTadCol, inNoTadCol, distClassCol) {
  counts <- matrix(0, 3, 3)
  counts[1, 1] <- length(which(df[[inSameTadCol]] & df[[distClassCol]] == "transMqtl"))
  counts[1, 2] <- length(which(df[[inSameTadCol]] & df[[distClassCol]] == "cisMqtl"))
  counts[1, 3] <- length(which(df[[inSameTadCol]]))
  counts[2, 1] <- length(which(!df[[inSameTadCol]] & !df[[inNoTadCol]] & df[[distClassCol]] == "transMqtl"))
  counts[2, 2] <- length(which(!df[[inSameTadCol]] & !df[[inNoTadCol]] & df[[distClassCol]] == "cisMqtl"))
  counts[2, 3] <- length(which(!df[[inSameTadCol]] & !df[[inNoTadCol]]))
  counts[3, 1] <- length(which(!df[[inSameTadCol]] & df[[inNoTadCol]] & df[[distClassCol]] == "transMqtl"))
  counts[3, 2] <- length(which(!df[[inSameTadCol]] & df[[inNoTadCol]] & df[[distClassCol]] == "cisMqtl"))
  counts[3, 3] <- length(which(!df[[inSameTadCol]] & df[[inNoTadCol]]))
  counts
}

# Function to calculate the percentage
calculatePercentage <- function(counts) {
  apply(counts, 2, function(x) x * 100 / sum(x, na.rm = TRUE))
}

# Function to create the barplot
createBarplot <- function(results, title, legendNames) {
  pdf(paste0("/home/og16379/VaryingDNAm/plots/groupedBar", title, ".pdf"))
  par(mar = c(4, 4, 4, 10.5))
  par(xpd = NA)
  barplot(results, col = cbbPalette[5:7], border = "white", cex.names = 0.4, ylim = c(0, 100), beside = TRUE)
  legend(3.7, 70, legendNames, fill = cbbPalette[5:6], bty = "n")
  dev.off()
}

# Calculate counts and percentages for All data
resultsMqtlAll <- calculateCounts(variableMqtl, "inSameTad", "inNoTad", "distClass")
resultsMqtlPercentageAll <- calculatePercentage(resultsMqtlAll)
createBarplot(resultsMqtlPercentageAll, "All", rownames(resultsMqtlPercentageAll))

# Calculate counts and percentages for trans data
resultsMqtlTrans <- calculateCounts(variableMqtl, "inSameTad", "inNoTad", "distClass")
resultsMqtlPercentageTrans <- calculatePercentage(resultsMqtlTrans)
createBarplot(resultsMqtlPercentageTrans, "Trans", rownames(resultsMqtlPercentageTrans))

# Calculate counts and percentages for cis data
resultsMqtlcis <- calculateCounts(variableMqtl, "inSameTad", "inNoTad", "distClass")
resultsMqtlPercentagecis <- calculatePercentage(resultsMqtlcis)
createBarplot(resultsMqtlPercentagecis, "cis", rownames(resultsMqtlPercentagecis))


# now check if connected by loops
anchorsNotConnectedByLoops <- function(mqtlAnno,probesOfInterest,loopsGi){
anchorBoth <- subset(mqtl,gene %in% probesOfInterest)
anchorBoth <- as.data.frame(anchorBoth %>% distinct(gene,.keep_all=TRUE))
anchorOne <- GRanges(anchorBoth[,1], IRanges(anchorBoth[,2],anchorBoth[,2]))# interaction start
anchorTwo <- GRanges(anchorBoth[,9], IRanges(anchorBoth[,10],anchorBoth[,10]))# interaction start
seqlevelsStyle(anchorOne) <- "UCSC"
seqlevelsStyle(anchorTwo) <- "UCSC"
olapA1 <- findOverlaps(anchorOne, loopsGi)
olapA2 <- findOverlaps(anchorTwo, loopsGi)
anchorOneLoop <- rep(0, length(anchorOne))
anchorTwoLoop <- rep(0, length(anchorTwo))
anchorOneLoop[queryHits(olapA1)] <- subjectHits(olapA1)
anchorTwoLoop[queryHits(olapA2)] <- subjectHits(olapA2)
ids_noLoop <- anchorOneLoop==0 & anchorTwoLoop == 0
return(ids_noLoop)
}


stableMqtl$noLoop <- anchorsNotConnectedByLoops(mqtl,replicatedResultsStable,loopsGi)
variableMqtl$noLoop <- anchorsNotConnectedByLoops(mqtl,replicatedResultsVariable,loopsGi)
chanceMqtl$noLoop <- anchorsNotConnectedByLoops(mqtl,backgroundProbes,loopsGi)


# now check if connected by loops
anchorsConnectedByLoops <- function(mqtlAnno,probesOfInterest,loopsGi){
anchorBoth <- subset(mqtl,gene %in% probesOfInterest)
anchorBoth <- as.data.frame(anchorBoth %>% distinct(gene,.keep_all=TRUE))
anchorOne <- GRanges(anchorBoth[,1], IRanges(anchorBoth[,2],anchorBoth[,2]))# interaction start
anchorTwo <- GRanges(anchorBoth[,9], IRanges(anchorBoth[,10],anchorBoth[,10]))# interaction start
seqlevelsStyle(anchorOne) <- "UCSC"
seqlevelsStyle(anchorTwo) <- "UCSC"
olapA1 <- findOverlaps(anchorOne, loopsGi)
olapA2 <- findOverlaps(anchorTwo, loopsGi)
anchorOneLoop <- rep(0, length(anchorOne))
anchorTwoLoop <- rep(0, length(anchorTwo))
anchorOneLoop[queryHits(olapA1)] <- subjectHits(olapA1)
anchorTwoLoop[queryHits(olapA2)] <- subjectHits(olapA2)
ids_noLoop <- anchorOneLoop!=0 & anchorOneLoop == anchorTwoLoop
return(ids_noLoop)
}


stableMqtl$cbLoop <- anchorsConnectedByLoops(mqtl,replicatedResultsStable,loopsGi)
variableMqtl$cbLoop <- anchorsConnectedByLoops(mqtl,replicatedResultsVariable,loopsGi)
chanceMqtl$cbLoop <- anchorsConnectedByLoops(mqtl,backgroundProbes,loopsGi)


# Function to calculate counts for different datasets
calculateCounts <- function(df, cbLoopCol, noLoopCol, distClassCol) {
  counts <- matrix(0, 3, 3)
  counts[1, 1] <- sum(df[[cbLoopCol]] & df[[distClassCol]] == "transMqtl")
  counts[1, 2] <- sum(df[[cbLoopCol]] & df[[distClassCol]] == "cisMqtl")
  counts[1, 3] <- sum(df[[cbLoopCol]])
  counts[2, 1] <- sum(!df[[cbLoopCol]] & !df[[noLoopCol]] & df[[distClassCol]] == "transMqtl")
  counts[2, 2] <- sum(!df[[cbLoopCol]] & !df[[noLoopCol]] & df[[distClassCol]] == "cisMqtl")
  counts[2, 3] <- sum(!df[[cbLoopCol]] & !df[[noLoopCol]])
  counts[3, 1] <- sum(!df[[cbLoopCol]] & df[[noLoopCol]] & df[[distClassCol]] == "transMqtl")
  counts[3, 2] <- sum(!df[[cbLoopCol]] & df[[noLoopCol]] & df[[distClassCol]] == "cisMqtl")
  counts[3, 3] <- sum(!df[[cbLoopCol]] & df[[noLoopCol]])
  counts
}

# Function to calculate percentages
calculatePercentage <- function(counts) {
  apply(counts, 2, function(x) x * 100 / sum(x, na.rm = TRUE))
}

# Function to create barplots
createBarplot <- function(results, title, legendNames) {
  pdf(paste0("/home/og16379/VaryingDNAm/plots/groupedBar", title, ".pdf"))
  par(mar = c(4, 4, 4, 10.5))
  par(xpd = NA)
  barplot(results, col = cbbPalette[5:7], border = "white", cex.names = 0.4, ylim = c(0, 100), beside = TRUE)
  legend(3.7, 70, legendNames, fill = cbbPalette[5:6], bty = "n")
  dev.off()
}

# Calculate counts and percentages for All data
resultsMqtlAllLoops <- calculateCounts(variableMqtl, "cbLoop", "noLoop", "distClass")
resultsMqtlAllLoopsPercentage <- calculatePercentage(resultsMqtlAllLoops)
createBarplot(resultsMqtlAllLoopsPercentage, "LoopsAll", rownames(resultsMqtlAllLoopsPercentage))

# Calculate counts and percentages for trans data
resultsMqtlTransLoops <- calculateCounts(variableMqtl, "cbLoop", "noLoop", "distClass")
resultsMqtlTransLoopsPercentage <- calculatePercentage(resultsMqtlTransLoops)
createBarplot(resultsMqtlTransLoopsPercentage, "LoopsTrans", rownames(resultsMqtlTransLoopsPercentage))

# Calculate counts and percentages for cis data
resultsMqtlcisLoops <- calculateCounts(variableMqtl, "cbLoop", "noLoop", "distClass")
resultsMqtlcisLoopsPercentage <- calculatePercentage(resultsMqtlcisLoops)
createBarplot(resultsMqtlcisLoopsPercentage, "Loopscis", rownames(resultsMqtlcisLoopsPercentage))

# Calculate counts and percentages for epiallele mQTL relationship
resultsMqtlEpi <- matrix(0, 1, 3)
resultsMqtlEpi[1, 1] <- nrow(subset(epiAllelesImVAtCleaned, Name %!in% mqtl$gene))
resultsMqtlEpi[1, 2] <- sum(epiAlleleMqtl@elementMetadata@listData$distClass == "cisMqtl")
resultsMqtlEpi[1, 3] <- sum(epiAlleleMqtl@elementMetadata@listData$distClass == "transMqtl")
resultsMqtlEpiPercentage <- calculatePercentage(resultsMqtlEpi)
createBarplot(resultsMqtlEpiPercentage, "epiMqtlRel", colnames(resultsMqtlEpiPercentage))



# load packages 

library(BSgenome)
library(GenomicFeatures)
library(GenomeInfoDbData)
library(BSgenome.Hsapiens.UCSC.hg19)

# get chrom info for hg19
hg19 <- getChromInfoFromUCSC("hg19")
hg19Coordinates <- GRanges(seqnames=hg19[,1],ranges=IRanges(end=hg19[,2],width=hg19[,2]))
ranges(hg19Coordinates)
hg19CoordinatesDt <- cbind(hg19[,1],as.data.frame(ranges(hg19Coordinates)[,c(1,2)]))
dim(hg19CoordinatesDt)
# write to file
write.table(file="/home/og16379/VaryingDNAm/data/hg19Coordinates.bed",hg19CoordinatesDt,sep="\t",quote=F,row.names=F, col.names = c("chr","start","end","width"))
# load in data needed
variableProbes <- import("/home/og16379/VaryingDNAm/objects/replicatedResultsVariableMapped.bed")
stableProbes <- import("/home/og16379/VaryingDNAm/objects/replicatedResultsStableMapped.bed")

# set window size
windowSize <- 500000
seqs = seq(hg19CoordinatesDt[1,2],hg19CoordinatesDt[1,3] - windowSize, windowSize)
ranges <- GRanges(hg19CoordinatesDt[1,1],IRanges(seqs,seqs+windowSize-1))

getRangesForChr <- function(input,windowSize,chrNumber){
seqs = seq(hg19CoordinatesDt[as.numeric(paste(chrNumber)),2],hg19CoordinatesDt[as.numeric(paste(chrNumber)),3] - windowSize, windowSize)
ranges <- GRanges(hg19CoordinatesDt[as.numeric(paste(chrNumber)),1],IRanges(seqs,seqs+windowSize-1))
return(ranges)
}


chr1Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,1)
chr2Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,2)
chr3Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,3)
chr4Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,4)
chr5Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,5)
chr6Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,6)
chr7Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,7)
chr8Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,8)
chr9Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,9)
chr10Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,10)
chr11Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,11)
chr12Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,12)
chr13Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,13)
chr14Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,14)
chr15Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,15)
chr16Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,16)
chr17Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,17)
chr18Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,18)
chr19Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,19)
chr20Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,20)
chr21Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,21)
chr22Ranges <- getRangesForChr(hg19CoordinatesDt,windowSize,22)

gL <- GRangesList(chr1Ranges,chr2Ranges,chr3Ranges,chr4Ranges,chr5Ranges,chr6Ranges,chr7Ranges,
chr8Ranges,chr9Ranges,chr10Ranges,chr11Ranges,chr12Ranges,chr13Ranges,chr14Ranges,chr15Ranges,
chr16Ranges,chr17Ranges,chr18Ranges,chr19Ranges,chr20Ranges,chr21Ranges,chr22Ranges)

glistVariable <- gL
glistStable <- gL

gLVariable <- lapply(gL,FUN=countOverlaps,variableProbes)
gLStable <- lapply(gL,FUN=countOverlaps,stableProbes)


for(i in 1:length(gLVariable)){
glistVariable[[i]]$value <- "NA" <- gLVariable[[i]]
}

for(i in 1:length(gLStable)){
glistStable[[i]]$value <- "NA" <- gLStable[[i]]
}


getAllChr <- function(glist){
glist <- do.call("c",glist)
glist <- as.data.frame(glist)
glist <- glist[,-4]
glist <- glist[,-4]
colnames(glist) <- c("chr","start","end","value")
return(glist)
}

glistStable <- getAllChr(glistStable)
glistVariable <- getAllChr(glistVariable)
head(glistStable)
head(glistVariable)


#extract the sequence for each bins
library(BSgenome.Hsapiens.UCSC.hg19)
binsSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg19,gL)
seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19)
binsSeq <- do.call("c",binsSeq)
binsCpgCount <- t(vcountPDict(DNAStringSet(c("CG")),binsSeq))
cpgDensity <- as.data.frame(do.call("c",gL))
cpgDensity$value <- binsCpgCount
cpgDensity <-data.frame("chr"=cpgDensity$seqnames,"start"=cpgDensity$start,"end"=cpgDensity$end,"value"=cpgDensity$value)


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
variableColour<-cbbPalette[3]
stableColour<-cbbPalette[4]


pdf("/home/og16379/VaryingDNAm/plots/circosFig1_2.pdf")
set.seed(123)
circos.initializeWithIdeogram(chromosome.index = paste0("chr", c(1:22)))
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.rect(xlim[1], 0, xlim[2], 1, col = cbPalette[8])
    circos.text(mean(xlim), mean(ylim), chr, cex = 0.4, col = "black",
        facing = "inside", niceFacing = TRUE)
}, track.height = 0.05, bg.border = NA)
circos.par("track.height" = 0.2)
circos.genomicTrack(glistStable,
    panel.fun = function(region, value,...) {
        circos.genomicLines(region,value,type="h",col =stableColour, ...)
})
circos.genomicTrack(glistVariable,
    panel.fun = function(region, value,...) {
        circos.genomicLines(region,value,type="h",col =variableColour, ...)
})
circos.genomicTrack(cpgDensity, 
    panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value,area=FALSE,col=cbPalette[2])
})
dev.off()



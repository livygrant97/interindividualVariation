################################################################
# Characterising Autosomal sex differences between males and females
# Olivia Grant
###############################################################

# below is all the code to produce the figures shown in figure 3A-3D
# previous code must have been run in order to produce these figures


################################################################################
# load libraries
################################################################################

if(!require("bigmelon", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("bigmelon")
}
library(bigmelon)

if(!require("Gviz", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("Gviz")
}
library(Gviz)

if(!require("IlluminaHumanMethylationEPICanno.ilm10b2.hg19", character.only = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
}
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

 
################################################################################
# prepare data for plotting
################################################################################

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

# open gds file
gfile <- openfn.gds('/storage/st05d/Exeter/UnderSocMeth/USM1/USM_WF.gds',
                      readonly=T)
# convert gds to mset object
gfile1<-gds2mset(gfile)

# map it to the genome
msetMapped<-mapToGenome(gfile1)

# convert that to a g range object
msetMapped<-granges(msetMapped)

topEpialleles <- epiAllelesImVAtCleaned[epiAllelesImVAtCleaned$Name %in% rownames(variableImMatrix)[1:112],]
topEpialleles <- msetMapped[msetMapped@ranges@NAMES %in% topEpialleles$Name[1:112],]

epialleleGRanges <- msetMapped[msetMapped@ranges@NAMES %in% epialleleGranges$Name,]

# load in extra libraries
library(rtracklayer)
library(RColorBrewer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Homo.sapiens)
library(BSgenome.Hsapiens.UCSC.hg19)


#get txdb object
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#set genome
gen <- "hg19"
pal <- brewer.pal(2,"Dark2")
#set info for plotting
chrom <- as.character(seqnames(topEpialleles[53]))
start <- as.numeric(start(topEpialleles[53]))
end <- as.numeric(end(topEpialleles[53]+10))
#add extra 25% to plot, change as you see fit for each plot
minbase <- start - 5000 * (end - start)
maxbase <- end + 5000 * (end - start)

# genome axis track
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
#gene region track
allGenes <- genes(txdb)
allGenes$symbol = mapIds(Homo.sapiens, keys=allGenes$gene_id, keytype="ENTREZID",
  column="SYMBOL")
grTrack <- GeneRegionTrack(allGenes, showId=TRUE,shape="arrow")


# i set ideoTrack manually using cytoband data from UCSC
  data<-read.table("/home/og16379/diff_cpg_fm/data/cytoBand.txt",header=F,sep="\t")
  colnames(data) <-c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain')
  ideoTrack <- IdeogramTrack(genome="hg19",chromosome = chrom , bands=data, from=start,to=end)

sTrack <- SequenceTrack(Hsapiens,chrom)

cpgIslands <- read.table("/home/og16379/VaryingDNAm/data/cpgIslands.bed")
colnames(cpgIslands) <- c("bin","chrom","start","end","length","cpgNum","gcNum","perCpg","perGc","obsExp","name")
cpgIslands <-makeGRangesFromDataFrame(cpgIslands,keep.extra.columns=FALSE)



getTopValues <- function(x){
  chr <- as.character(seqnames(x[1:100]))
  start <- as.numeric(start(x[1:100]))
  end <- as.numeric(end(x[1:100])+10)
  minbase <- start  - 10000
  maxbase <- end + 9990
  valuesToPlot <- data.frame(chr,start,end,minbase,maxbase)
  return(valuesToPlot)
}

valuesToPlot <- getTopValues(topEpialleles)

for (i in 1:nrow(valuesToPlot)) {
ideoTrack <- IdeogramTrack(genome="hg19",chromosome = valuesToPlot$chr[i] , bands=data, from=valuesToPlot$start[i],to=valuesToPlot$end[i])
sTrack <- SequenceTrack(Hsapiens,as.character(valuesToPlot$chr[i]))
pdf(paste0("/home/og16379/VaryingDNAm/plots/gvizPlots/2023/",valuesToPlot$chr[i],";",valuesToPlot$minbase[i],"-",valuesToPlot$maxbase[i],".pdf"))
Gviz::plotTracks(list(ideoTrack,gTrack,sTrack,grTrack,AnnotationTrack(cpgIslands,start=valuesToPlot$start[i],end=valuesToPlot$end[i],chr=valuesToPlot$chr[i],name="CpG Islands"),
AnnotationTrack(topEpialleles,start=valuesToPlot$start[i],end=valuesToPlot$end[i],chr=valuesToPlot$chr[i],name="EpiAlleles")),from=valuesToPlot$minbase[i],to=valuesToPlot$maxbase[i])
  dev.off()
}


pdf("/home/og16379/VaryingDNAm/plots/histogram20231.pdf",width=8,height=3.5)
hist(as.numeric(variableImMatrix[1166,1:1175]),breaks=50,col=cbPalette[6])
dev.off()

pdf("/home/og16379/VaryingDNAm/plots/histogram20232.pdf",width=8,height=3.5)
hist(as.numeric(variableImMatrix[1806,1:1175]),breaks=50,col=cbPalette[6])
dev.off()

pdf("/home/og16379/VaryingDNAm/plots/histogram20233.pdf",width=8,height=3.5)
hist(as.numeric(variableImMatrix[1480,1:1175]),breaks=50,col=cbPalette[6])
dev.off()

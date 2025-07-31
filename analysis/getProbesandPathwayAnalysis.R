################################################################
# Identify stable and variable CpG sites on the EPIC array
# Olivia Grant
################################################################

###############################################################
#### Setup: Libraries & Paths
###############################################################

# Define all file paths here to avoid hard-coding
gds_discovery_path <- "/storage/st05d/Exeter/UnderSocMeth/USM1/USM_WF.gds"
gds_validation_path <- "/storage/st05d/Exeter/UnderSocMeth/ISER/USM2_WF.gds"
normalised_discovery_data <- "/home/og16379/diff_cpg_fm/data/firstDatasetNormalised.RData"
normalised_validation_data <- "/home/og16379/diff_cpg_fm/data/dasen_autosome2round.Rdata"
epic_manifest_path <- "/home/og16379/VaryingDNAm/data/EPIC.hg19.manifest.tsv"
output_dir <- "/home/og16379/VaryingDNAm"

# Create output subfolders if not exist
dir.create(file.path(output_dir, "objects"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "csv"), showWarnings = FALSE, recursive = TRUE)

# Function to load packages and install if needed
load_or_install <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Load required packages
pkgs <- c("bigmelon", "rtracklayer", 
          "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", 
          "lattice", "dplyr", "missMethyl", 
          "IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
invisible(lapply(pkgs, load_or_install))

###############################################################
#### Utility Functions
###############################################################

# Identify variable CpGs by standard deviation ranking
getVariableCpGs <- function(data, proportion, top = TRUE) {
  sd <- apply(data, 1, sd)
  x <- cbind(data, sd)
  x <- x[order(x$sd, decreasing = TRUE), ]
  numberToSelect <- round(proportion * nrow(x))
  toSelect <- if (top) 1:numberToSelect else (nrow(x) - numberToSelect + 1):nrow(x)
  rownames(x[toSelect, ])
}

# Randomly downsample columns (samples) from a matrix
downsampleData <- function(data, downsampleProportion) {
  numberToSelect <- round(downsampleProportion * ncol(data))
  data[, sample(1:ncol(data), numberToSelect, replace = FALSE)]
}

# Repeat variable CpG selection on downsampled datasets and keep consensus sites
getRobustVariableCpGs <- function(data, proportion, top = TRUE,
                                  downsampleProportion, downsampleTimes = 10,
                                  foundIn = 1.0) {
  variableCpGs <- matrix("", nrow = round(proportion * nrow(data)), ncol = downsampleTimes)
  for (i in 1:downsampleTimes) {
    downsampledData <- downsampleData(data, downsampleProportion)
    variableCpGs[, i] <- getVariableCpGs(downsampledData, proportion, top)
  }
  counts <- table(as.vector(variableCpGs))
  names(counts[counts >= floor(foundIn * downsampleTimes)])
}

# Get standard deviation values for a set of CpGs
getSD <- function(cpgs, betaMatrix) {
  x <- betaMatrix[rownames(betaMatrix) %in% cpgs, ]
  sd <- apply(x, 1, sd)
  x <- as.data.frame(cbind(x, sd))
  x <- x[order(x$sd, decreasing = TRUE), ]
  x[startsWith(rownames(x), "cg"), ]
}

###############################################################
#### Set Parameters
###############################################################

proportion <- 0.10
downsampleProportion <- 0.90
downsampleTimes <- 10
foundIn <- 1.0

###############################################################
#### Load and clean discovery dataset
###############################################################

# Load and subset discovery data
gfile <- openfn.gds(gds_discovery_path, readonly = TRUE)
load(normalised_discovery_data)

info <- pData(gfile)[, c("barcode", "nsex")]
info$nsex <- gsub("2", "F", info$nsex)
info$nsex <- gsub("1", "M", info$nsex)

outliers_discovery <- c("200611820013_R08C01", "200603220075_R08C01",
                        "200864580018_R02C01", "200611820020_R07C01")
info <- info[!(info$barcode %in% outliers_discovery), ]
dasen_autosome <- dasen_autosome[, info$barcode]
dasen_autosomeDf <- as.data.frame(dasen_autosome)

###############################################################
#### Clean EPIC annotation & remove sex/SNP-affected probes
###############################################################

data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
anno_df <- as.data.frame(anno)

# Remove sex chromosome probes
sexProbes <- as.character(anno_df$Name[anno_df$chr %in% c("chrX", "chrY")])
CpGs_autosomes <- dasen_autosomeDf[!rownames(dasen_autosomeDf) %in% sexProbes, ]

# Remove probes flagged as SNP-associated
SNP_anno <- read.delim(epic_manifest_path)
SNP_masked <- SNP_anno %>% filter(MASK_general == TRUE)
SNP_masked <- SNP_masked[SNP_masked$probeID %in% rownames(CpGs_autosomes), ]
CpGs_cleaned <- CpGs_autosomes[!rownames(CpGs_autosomes) %in% SNP_masked$probeID, ]

save(CpGs_cleaned, file = file.path(output_dir, "objects", "CpGs_cleanedDiscovery.RData"))

###############################################################
#### Get stable and variable CpGs - discovery
###############################################################

hvCpGs <- getRobustVariableCpGs(CpGs_cleaned, proportion, TRUE, downsampleProportion, downsampleTimes, foundIn)
hsCpGs <- getRobustVariableCpGs(CpGs_cleaned, proportion, FALSE, downsampleProportion, downsampleTimes, foundIn)

save(hvCpGs, file = file.path(output_dir, "objects", "highlyVariableRobustCpGsDiscovery.RData"))
save(hsCpGs, file = file.path(output_dir, "objects", "highlyStableRobustCpGsDiscovery.RData"))

###############################################################
#### Repeat for validation dataset
###############################################################

gfile2 <- openfn.gds(gds_validation_path)
load(normalised_validation_data)

info2 <- pData(gfile2)[, c("barcode", "Sex")]
info2$Sex <- gsub("2", "F", info2$Sex)
info2$Sex <- gsub("1", "M", info2$Sex)

outliers_validation <- c("203991410050_R06C01", "203991410050_R07C01", "204026590012_R07C01", 
                         "203998240051_R07C01", "204026590078_R04C01", "203994670097_R08C01", 
                         "203960330159_R08C01", "204022160070_R03C01", "203994670097_R02C01", 
                         "203991470090_R05C01", "204026590012_R06C01")
info2 <- info2[!(info2$barcode %in% outliers_validation), ]
dasen_autosome2 <- dasen_autosome2[, info2$barcode]
dasen_autosome2Df <- as.data.frame(dasen_autosome2)

# Repeat cleaning
CpGs_autosomes2 <- dasen_autosome2Df[!rownames(dasen_autosome2Df) %in% sexProbes, ]
SNP_masked2 <- SNP_masked[SNP_masked$probeID %in% rownames(CpGs_autosomes2), ]
CpGs_cleaned2 <- CpGs_autosomes2[!rownames(CpGs_autosomes2) %in% SNP_masked2$probeID, ]

save(CpGs_cleaned2, file = file.path(output_dir, "objects", "CpGs_cleanedValidation.RData"))

# Repeat variable/stable detection
hvCpGs2 <- getRobustVariableCpGs(CpGs_cleaned2, proportion, TRUE, downsampleProportion, downsampleTimes, foundIn)
hsCpGs2 <- getRobustVariableCpGs(CpGs_cleaned2, proportion, FALSE, downsampleProportion, downsampleTimes, foundIn)

save(hvCpGs2, file = file.path(output_dir, "objects", "highlyVariableRobustCpGsValidation.RData"))
save(hsCpGs2, file = file.path(output_dir, "objects", "highlyStableRobustCpGsValidation.RData"))


# variantSitesRdaFiles.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


createSingleChromosomeVariantSitesRdaFile <- function(
  vcfFilename                 = file.path("/data/galton/users/rpearson", "crossesTesting", "release", "7g8xGb4-qcPlusSamples-0.1.vcf.gz"),
  chromosome                  = "MAL1",
  shouldRemoveInvariant       = TRUE,
#  regionsMask                 = varRegions_v2(), # will remove any variants in these regions. Set to NULL if you don't want to mask any variants out in this way
  regionsMask                 = NULL,
#  filtersToRemove             = c("NoAlternative"),
  filtersToRemove             = NULL,
  samplesToRemove             = NULL,
  possibleMissingValues       = c(".", "./.", ".|."),
  outputRdaFilename           = sub("\\.vcf[\\.gz]*$", paste("\\.variantSites\\.", chromosome,"\\.rda", sep=""), vcfFilename),
  overwriteExisting           = !(file.exists(outputRdaFilename)),
  saveAsRobjectFile           = TRUE
) {
  if(vcfFilename == outputRdaFilename) {
    stop("Input and output filename are the same")
  }
  if(overwriteExisting) {
#    rng <- GRanges(seqnames=chromosome, ranges=IRanges(1, 3.3e+6)) # slightly more than largest Pf chromosome - not very elegant but there you go
    rng <- GRanges(seqnames=chromosome, ranges=IRanges(1, 536870912)) # decided there should be no reason to make the above Pf specific, so went for a maximum very large chromosome
    param <- ScanVcfParam(which=rng, geno=c("GT", "AD"))
    if(grepl("\\.gz", vcfFilename) & !file.exists(paste(vcfFilename, "tbi", sep="."))) {
      indexTabix(vcfFilename, "vcf4")
    }
    vcf <- readVcf(vcfFilename, "Pf", param)
    if(!is.null(samplesToRemove)) {
      sampleToKeep <- setdiff(row.names(colData(vcf)), samplesToRemove)
      vcf <- vcf[, sampleToKeep]
    }
    if(shouldRemoveInvariant) {
#      invariantSNPs <- apply(geno(vcf)[["GT"]], 1, function(x) length(which(x=="0")) == 0 | length(which(x=="1")) == 0)
      invariantSNPs <- apply(geno(vcf)[["GT"]], 1, function(x) length(table(x[!(x %in% possibleMissingValues)], useNA="no"))==1)
      vcf <- vcf[!invariantSNPs]
    }
    if(!is.null(filtersToRemove)) {
      sapply(
        filtersToRemove,
        function(filterToRemove) {
          vcf <<- vcf[!grepl(filterToRemove, filt(vcf))]
        }
      )
    }
    if(!is.null(regionsMask)) {
      vcf <- vcf[is.na(GenomicRanges::match(rowData(vcf), regionsMask))]
    }
    if(saveAsRobjectFile) {
      save(vcf, file=outputRdaFilename)
    }
  } else {
    load(outputRdaFilename)
  }
  return(vcf)
}

variantSitesRdaFiles <- function(
  vcfFilename                 = file.path("/data/galton/users/rpearson", "crossesTesting", "release", "7g8xGb4-qcPlusSamples-0.1.vcf.gz"),
  chromosomes                 = c(paste("MAL", 1:14, sep="")),
#  chromosomes                 = c(paste("MAL", 1:14, sep=""), "MITO1", paste("PLAST", 1:2, sep=""))
  filtersToRemove             = NULL,
  samplesToRemove             = NULL,
  overwriteExisting           = NULL
)  {
  if(is.null(overwriteExisting)) {
    lapply(
      chromosomes,
      function(x) {
        createSingleChromosomeVariantSitesRdaFile(vcfFilename, x, filtersToRemove=filtersToRemove, samplesToRemove=samplesToRemove)
      }
    )
  } else {
    lapply(
      chromosomes,
      function(x) {
        createSingleChromosomeVariantSitesRdaFile(vcfFilename, x, filtersToRemove=filtersToRemove, samplesToRemove=samplesToRemove, overwriteExisting=overwriteExisting)
      }
    )
  }
}

#variantSitesRdaFiles <- function(
#  vcfFilenames                = file.path("/data/galton/users/rpearson", "crossesTesting", "release", paste(c("7g8xGb4", "3d7xHb3", "Hb3xDd2"), "-qcPlusSamples-0.1.vcf.gz", sep="")),
#  chromosomes                 = c(paste("MAL", 1:14, sep="")),
##  chromosomes                 = c(paste("MAL", 1:14, sep=""), "MITO1", paste("PLAST", 1:2, sep=""))
#  overwriteExisting           = NULL
#)  {
#  if(is.null(overwriteExisting)) {
##    clusterLapply(
#    lapply(
#      mapply(
#        function(x, y) {list(vcfFilename=x, chromosome=y)},
#        rep(vcfFilenames, each=length(chromosomes)),
#        rep(chromosomes, length(vcfFilenames)),
#        SIMPLIFY = FALSE,
#        USE.NAMES = FALSE
#      ),
#      "createSingleChromosomeVariantSitesRdaFile"
##      maxRAMinGB                  = 16
#    )
#  } else {
#    lapply(
#      mapply(
#        function(x, y) {list(vcfFilename=x, chromosome=y, overwriteExisting=overwriteExisting)},
#        rep(vcfFilenames, each=length(chromosomes)),
#        rep(chromosomes, length(vcfFilenames)),
#        SIMPLIFY = FALSE,
#        USE.NAMES = FALSE
#      ),
#      "createSingleChromosomeVariantSitesRdaFile"
##      maxRAMinGB                  = 16
#    )
#  }
#}
#

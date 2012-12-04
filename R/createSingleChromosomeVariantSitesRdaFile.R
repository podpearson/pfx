# createSingleChromosomeVariantSitesRdaFile.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


createSingleChromosomeVariantSitesRdaFile <- function(
  vcfFilename                 = file.path("/data/galton/users/rpearson", "crossesTesting", "release", "7g8xGb4-qcPlusSamples-0.1.vcf.gz"),
  chromosome                  = "MAL1",
  genoToLoad                  = c("GT", "AD"),
  shouldRemoveInvariant       = TRUE,
#  regionsMask                 = varRegions_v2(), # will remove any variants in these regions. Set to NULL if you don't want to mask any variants out in this way
  regionsMask                 = NULL,
#  shouldSetMultiallelicFilter = TRUE,
  keepPASSvariantsOnly        = FALSE,
#  filtersToRemove             = c("NoAlternative"),
  filtersToRemove             = NULL,
  samplesToRemove             = NULL,
  additionalInfoFilters       = NULL,
#  additionalInfoFilters     = list(
#    "SVTYPE" = list(operator="%in%", value="SNP")
#  ),
  additionalGenotypeFilters   = NULL,
#  additionalGenotypeFilters     = list(
#    "GT_CONF" = list(operator="<=", value=20),
#    "SITE_CONF" = list(operator="<=", value=200)
#  )
  possibleMissingValues       = c(".", "./.", ".|."),
  outputRdaFilename           = sub("\\.vcf[\\.gz]*$", paste("\\.variantSites\\.", chromosome,"\\.rda", sep=""), vcfFilename),
  overwriteExisting           = !(file.exists(outputRdaFilename)),
  saveAsRobjectFile           = TRUE,
  parentalStrains             = NULL
) {
  if(vcfFilename == outputRdaFilename) {
    stop("Input and output filename are the same")
  }
  if(overwriteExisting) {
#    rng <- GRanges(seqnames=chromosome, ranges=IRanges(1, 3.3e+6)) # slightly more than largest Pf chromosome - not very elegant but there you go
    rng <- GRanges(seqnames=chromosome, ranges=IRanges(1, 536870912)) # decided there should be no reason to make the above Pf specific, so went for a maximum very large chromosome
    param <- ScanVcfParam(which=rng, geno=c(genoToLoad))
    if(grepl("\\.vcf$", vcfFilename)) {
      if(!file.exists(paste(vcfFilename, "gz", sep="."))) {
        bgzip(vcfFilename)
      }
      vcfFilename <- paste(vcfFilename, "gz", sep=".")
    }
    if(grepl("\\.gz$", vcfFilename) & !file.exists(paste(vcfFilename, "tbi", sep="."))) {
      indexTabix(vcfFilename, "vcf4")
    }
    vcf <- readVcf(vcfFilename, "Pf", param)
    vcf <- filterVcf(
      vcf,
      shouldRemoveInvariant       = shouldRemoveInvariant,
      regionsMask                 = regionsMask,
#      shouldSetMultiallelicFilter = shouldSetMultiallelicFilter,
      keepPASSvariantsOnly        = keepPASSvariantsOnly,
      filtersToRemove             = filtersToRemove,
      samplesToRemove             = samplesToRemove,
      additionalInfoFilters       = additionalInfoFilters,
      additionalGenotypeFilters   = additionalGenotypeFilters,
      possibleMissingValues       = possibleMissingValues
    )
#    if(!is.null(samplesToRemove)) {
#      sampleToKeep <- setdiff(row.names(colData(vcf)), samplesToRemove)
#      vcf <- vcf[, sampleToKeep]
#    }
#    if(shouldRemoveInvariant) {
##      invariantSNPs <- apply(geno(vcf)[["GT"]], 1, function(x) length(which(x=="0")) == 0 | length(which(x=="1")) == 0)
#      invariantSNPs <- apply(geno(vcf)[["GT"]], 1, function(x) length(table(x[!(x %in% possibleMissingValues)], useNA="no"))==1)
#      vcf <- vcf[!invariantSNPs]
#    }
#    if(keepPASSvariantsOnly) {
#      vcf <- vcf[filt(vcf)=="PASS"]
#    } else {
#      if(!is.null(filtersToRemove)) {
#        sapply(
#          filtersToRemove,
#          function(filterToRemove) {
#            vcf <<- vcf[!grepl(filterToRemove, filt(vcf))]
#          }
#        )
#      }
#    }
#    if(!is.null(regionsMask)) {
#      vcf <- vcf[is.na(GenomicRanges::match(rowData(vcf), regionsMask))]
#    }
    if(saveAsRobjectFile) {
      save(vcf, file=outputRdaFilename)
    }
  } else {
    load(outputRdaFilename)
  }
  if(!is.null(parentalStrains)) {
    newSampleOrdering <- c(parentalStrains, setdiff(dimnames(vcf)[[2]], parentalStrains))
    vcf <- vcf[, newSampleOrdering]
  }
  return(vcf)
}



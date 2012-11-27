# paintingSeries.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


paintingSeries <- function(
  vcfRdaFilename                = file.path("/data/galton/users/rpearson", "crossesTesting", "release", "7g8xGb4-qcPlusSamples-0.1.vcf.MAL4.rda"),
  pdfFilestem                   = sub("\\.rda", "\\.chromosomePaintingSeries", vcfRdaFilename),
  varMask                       = varRegions_v2(),
#  varMask                       = GRanges(
#    seqnames = "MAL4",
#    ranges   = IRanges(
#      start    = c(1,      552884, 939470, 1150294),
#      end      = c(100256, 617416, 987434, 1204112)
#    )
#  ),
#  col                           = c("white", "blue", "red", "lightblue", "pink", "lightgrey", "black"),
#  breaks                        = c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5),
  samplesToUse                  = NULL,
  width                         = 14,
  height                        = 4,
  filtersToRemove               = c("NoAlternative"),
  discordanceThreshold          = 500,
  verbose                       = TRUE
) {
  if(vcfRdaFilename==pdfFilestem) {
    stop("vcfRdaFilename is the same as pdfFilename")
  }
  if(verbose) {
    cat("paintingSeries: loading vcf...")
  }
  load(vcfRdaFilename)
  if(verbose) {
    cat("done\n")
  }
  if(verbose) {
    cat("paintingSeries: filtering vcf...")
  }
  sapply(
    filtersToRemove,
    function(filterToRemove) {
      vcf <- vcf[!grepl(filterToRemove, filt(vcf))]
    }
  )
  if(verbose) {
    cat("done\n")
  }
  if(verbose) {
    cat("paintingSeries: processing vcf...")
  }
  GTsInt <- genotypeCallsFromGTas012(vcf)  
  if(is.null(samplesToUse)) {
    samplesToUse <- uniqueSamples(vcf, discordanceThreshold=discordanceThreshold)
  }
  GTsIntToKeep <- GTsInt[, samplesToUse]
 
  GTsIntWithNAs <- GTsIntToKeep
  GTsIntWithNAs[GTsIntWithNAs==0] <- NA
  monomorphicRows <- apply(GTsIntWithNAs, 1, function(x) all(x==1, na.rm=TRUE)|all(x==2, na.rm=TRUE))
  rowsToUse <- (!GTsIntToKeep[, 1]==0 | !GTsIntToKeep[, 2]==0) & (!monomorphicRows)
  GTsCFparents <- convertGTsIntToParentBasedGTs(GTsIntToKeep[rowsToUse, ])
  if(verbose) {
    cat("done\n")
  }
  if(verbose) {
    cat("paintingSeries: creating recombination plot 1...")
  }
  pdf(paste(pdfFilestem, "raw.pdf", sep="."), height=height, width=width)
  recombinationPlot(GTsCFparents)
  dev.off()
  if(verbose) {
    cat("done\n")
  }
  missingGenotypesPerVariant <- apply(GTsIntToKeep, 1, function(x) length(which(x==0)))
  rowsToUse <- (!GTsIntToKeep[, 1]==0 | !GTsIntToKeep[, 2]==0) & (!monomorphicRows) & missingGenotypesPerVariant <= 1
  GTsCFparents <- convertGTsIntToParentBasedGTs(GTsIntToKeep[rowsToUse, ])
  if(verbose) {
    cat("done\n")
  }
  if(verbose) {
    cat("paintingSeries: creating recombination plot 2...")
  }
  pdf(paste(pdfFilestem, "removeMissing.pdf", sep="."), height=height, width=width)
  recombinationPlot(GTsCFparents)
  dev.off()
  if(verbose) {
    cat("done\n")
  }
  HeterozygousFilter <- grepl("Heterozygous", filt(vcf))
  rowsToUse <- (!GTsIntToKeep[, 1]==0 | !GTsIntToKeep[, 2]==0) & (!monomorphicRows) & missingGenotypesPerVariant <= 1 & !HeterozygousFilter
  GTsCFparents <- convertGTsIntToParentBasedGTs(GTsIntToKeep[rowsToUse, ])
  if(verbose) {
    cat("done\n")
  }
  if(verbose) {
    cat("paintingSeries: creating recombination plot 3...")
  }
  pdf(paste(pdfFilestem, "heterozygousFilter.pdf", sep="."), height=height, width=width)
  recombinationPlot(GTsCFparents)
  dev.off()
  if(verbose) {
    cat("done\n")
  }
  inVarMask <- is.na(GenomicRanges::match(rowData(vcf), varMask))
  rowsToUse <- (!GTsIntToKeep[, 1]==0 | !GTsIntToKeep[, 2]==0) & (!monomorphicRows) & missingGenotypesPerVariant <= 1 & !HeterozygousFilter & inVarMask
  GTsCFparents <- convertGTsIntToParentBasedGTs(GTsIntToKeep[rowsToUse, ])
  if(verbose) {
    cat("done\n")
  }
  if(verbose) {
    cat("paintingSeries: creating recombination plot 4...")
  }
  pdf(paste(pdfFilestem, "varMask.pdf", sep="."), height=height, width=width)
  recombinationPlot(GTsCFparents)
  dev.off()
  if(verbose) {
    cat("done\n")
  }
  MendelianErrorFilter <- grepl("MendelianErrors", filt(vcf))
  rowsToUse <- (!GTsIntToKeep[, 1]==0 | !GTsIntToKeep[, 2]==0) & (!monomorphicRows) & missingGenotypesPerVariant <= 1 & !HeterozygousFilter & inVarMask & !MendelianErrorFilter
  GTsCFparents <- convertGTsIntToParentBasedGTs(GTsIntToKeep[rowsToUse, ])
  if(verbose) {
    cat("done\n")
  }
  if(verbose) {
    cat("paintingSeries: creating recombination plot 5...")
  }
  pdf(paste(pdfFilestem, "mendelianError.pdf", sep="."), height=height, width=width)
  recombinationPlot(GTsCFparents)
  dev.off()
  if(verbose) {
    cat("done\n")
  }
  segregating <- values(info(vcf))[["SEGREGATING"]]
  rowsToUse <- (!GTsIntToKeep[, 1]==0 | !GTsIntToKeep[, 2]==0) & (!monomorphicRows) & missingGenotypesPerVariant <= 1 & !HeterozygousFilter & inVarMask & !MendelianErrorFilter & segregating
  GTsCFparents <- convertGTsIntToParentBasedGTs(GTsIntToKeep[rowsToUse, ])
  if(verbose) {
    cat("done\n")
  }
  if(verbose) {
    cat("paintingSeries: creating recombination plot 6...")
  }
  pdf(paste(pdfFilestem, "segregating.pdf", sep="."), height=height, width=width)
  recombinationPlot(GTsCFparents)
  dev.off()
  if(verbose) {
    cat("done\n")
  }

  return(TRUE)
}

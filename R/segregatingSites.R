# segregatingSites.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


segregatingSites <- function(
  vcf,
  pdfFilestem                   = paste(meta(exptData(vcf)[["header"]])["DataSetName", "Value"], seqlevels(vcf), "chromosomePaintingSeries", sep="."),
  varMask                       = varRegions_v2(),
  samplesToUse                  = row.names(colData(vcf)),
  width                         = 14,
  height                        = 4,
  filtersToRemove               = c("NoAlternative"),
  discordanceThreshold          = 500,
  suspiciousRegions             = c(
    GRanges("MAL6",  IRanges(1,       74373)),
    GRanges("MAL9",  IRanges(1214920, 1241458)),
    GRanges("MAL11", IRanges(1902245, 1935028)),
    GRanges("MAL12", IRanges(763134,  776156)),
    GRanges("MAL13", IRanges(1578814, 1591642)),
    GRanges("MAL14", IRanges(1,       35424)),
    GRanges("MAL14", IRanges(3193783, 3199050))
  ),
  additionalSuspiciousRegion    = GRanges("MAL14", IRanges(237968, 240153)),
  verbose                       = TRUE
) {
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
  if(!is.null(segregating)) {
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
  } else {
    segregating <- TRUE
  }
  SingleSNPHaplotypeFilter <- grepl("SingleSNPHaplotype", filt(vcf))
  rowsToUse <- (!GTsIntToKeep[, 1]==0 | !GTsIntToKeep[, 2]==0) & (!monomorphicRows) & missingGenotypesPerVariant <= 1 & !HeterozygousFilter & inVarMask & !MendelianErrorFilter & segregating & !SingleSNPHaplotypeFilter
  GTsCFparents <- convertGTsIntToParentBasedGTs(GTsIntToKeep[rowsToUse, ])
  if(verbose) {
    cat("done\n")
  }
  if(verbose) {
    cat("paintingSeries: creating recombination plot 7...")
  }
  pdf(paste(pdfFilestem, "singleSNPHaplotype.pdf", sep="."), height=height, width=width)
  recombinationPlot(GTsCFparents)
  dev.off()
  if(verbose) {
    cat("done\n")
  }
  rowsToUse <- (!GTsIntToKeep[, 1]==0 | !GTsIntToKeep[, 2]==0) & (!monomorphicRows) & missingGenotypesPerVariant <= 1 & !HeterozygousFilter & !MendelianErrorFilter & segregating & !SingleSNPHaplotypeFilter
  GTsCFparents <- convertGTsIntToParentBasedGTs(GTsIntToKeep[rowsToUse, ])
  if(verbose) {
    cat("done\n")
  }
  if(verbose) {
    cat("paintingSeries: creating recombination plot 8...")
  }
  pdf(paste(pdfFilestem, "previousRelease.pdf", sep="."), height=height, width=width)
  recombinationPlot(GTsCFparents)
  dev.off()
  if(verbose) {
    cat("done\n")
  }
  inSuspiciousRegionsMask <- is.na(GenomicRanges::match(rowData(vcf), suspiciousRegions))
  rowsToUse <- (!GTsIntToKeep[, 1]==0 | !GTsIntToKeep[, 2]==0) & (!monomorphicRows) & missingGenotypesPerVariant <= 1 & !HeterozygousFilter & inVarMask & !MendelianErrorFilter & segregating & inSuspiciousRegionsMask
  GTsCFparents <- convertGTsIntToParentBasedGTs(GTsIntToKeep[rowsToUse, ])
  if(verbose) {
    cat("done\n")
  }
  if(verbose) {
    cat("paintingSeries: creating recombination plot 9...")
  }
  pdf(paste(pdfFilestem, "multiSampleRecombinations.pdf", sep="."), height=height, width=width)
  recombinationPlot(GTsCFparents)
  dev.off()
  if(verbose) {
    cat("done\n")
  }
  inAdditionalSuspiciousRegionMask <- is.na(GenomicRanges::match(rowData(vcf), additionalSuspiciousRegion))
  rowsToUse <- (!GTsIntToKeep[, 1]==0 | !GTsIntToKeep[, 2]==0) & (!monomorphicRows) & missingGenotypesPerVariant <= 1 & !HeterozygousFilter & inVarMask & !MendelianErrorFilter & segregating & inSuspiciousRegionsMask & inAdditionalSuspiciousRegionMask
  GTsCFparents <- convertGTsIntToParentBasedGTs(GTsIntToKeep[rowsToUse, ])
  if(verbose) {
    cat("done\n")
  }
  if(verbose) {
    cat("paintingSeries: creating recombination plot 10...")
  }
  pdf(paste(pdfFilestem, "final.pdf", sep="."), height=height, width=width)
  recombinationPlot(GTsCFparents)
  dev.off()
  if(verbose) {
    cat("done\n")
  }

  rowsToUse <- (!GTsIntToKeep[, 1]==0 | !GTsIntToKeep[, 2]==0) & (!monomorphicRows) & missingGenotypesPerVariant <= 1 & !HeterozygousFilter & inVarMask & !MendelianErrorFilter & segregating & inSuspiciousRegionsMask & inAdditionalSuspiciousRegionMask
  return(vcf[rowsToUse, samplesToUse])
}


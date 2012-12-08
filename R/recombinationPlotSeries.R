# recombinationPlotSeries.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


recombinationPlotSeries <- function(
  vcf,
  filters                     = setdiff(unique(unlist(strsplit(filt(vcf), ";"))), c("PASS", ".")),
  chromosomes                 = sprintf("Pf3D7_%02d_v3", 1:14),
  sampleIDcolumn              = "ox_code",
  sampleDuplicates            = NULL,
  plotFilestem                = paste(meta(exptData(vcf)[["header"]])["DataSetName", "Value"], seqlevels(vcf), "chromosomePaintingSeries", sep="."),
  width                       = 14,
  height                      = 4,
  verbose                     = TRUE
) {
  returnMatrix <- sapply(
    chromosomes,
    function(chromosome) {
      GTsInt <- genotypeCallsFromGTas012(vcf[seqnames(vcf)==chromosome])  
      GTsCFparents <- convertGTsIntToParentBasedGTs(GTsInt)
      if(verbose) {
        cat("recombinationPlotSeries: creating recombination plot", chromosome, "raw", "\n")
      }
      pdf(paste(plotFilestem, chromosome, "raw", "pdf", sep="."), height=height, width=width)
      recombinationPlot(GTsCFparents)
      dev.off()
      returnValue <- dim(GTsCFparents)[1]
      names(returnValue) <- "raw"
      c(
        returnValue,
        sapply(
          seq(along=filters),
          function(numberOfFilters) {
            filtersToUse <- filters[1:numberOfFilters]
            filtersJoined <- paste(filtersToUse, collapse="_")
            variantsToRemoveMatrix <- sapply(
              filtersToUse,
              function(filter) {
                variantsToRemove <- grepl(filter, filt(vcf[seqnames(vcf)==chromosome]))
              }
            )
            variantsToRemove <- apply(variantsToRemoveMatrix, 1, any)
            GTsCFparents <- convertGTsIntToParentBasedGTs(GTsInt[!variantsToRemove, ])
            sampleIDmappings <- createSampleIDmappings(
              sampleIDs=dimnames(GTsCFparents)[[2]],
              sampleIDcolumn=sampleIDcolumn,
              sampleDuplicates=sampleDuplicates
            )
            if(verbose) {
              cat("recombinationPlotSeries: creating recombination plot", chromosome, filtersJoined, "\n")
            }
            pdf(paste(plotFilestem, chromosome, filtersJoined, "pdf", sep="."), height=height, width=width)
            recombinationPlot(GTsCFparents)
            dev.off()
            returnValue <- length(which(!variantsToRemove))
            names(returnValue) <- filtersJoined
            returnValue
          }
        )
      )
    },
    USE.NAMES=TRUE
  )
  return(returnMatrix)
}

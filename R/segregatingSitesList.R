# segregatingSitesList.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


segregatingSitesList <- function(
  vcfList,
  samplesToUse                  = row.names(colData(vcfList[[1]])),
  pdfFilestem                   = NULL,
  verbose                       = TRUE
) {
  sapply(
    vcfList,
    function(vcf) {
      if(is.null(pdfFilestem)) {
        segregatingSites(vcf, samplesToUse=samplesToUse, verbose=verbose, pdfFilestem=paste(meta(exptData(vcf)[["header"]])["DataSetName", "Value"], seqlevels(vcf), "chromosomePaintingSeries", sep="."))
      } else {
        segregatingSites(vcf, samplesToUse=samplesToUse, verbose=verbose, pdfFilestem=paste(pdfFilestem, seqlevels(vcf), "chromosomePaintingSeries", sep="."))
      }
    },
    simplify=FALSE,
    USE.NAMES=TRUE
  )
}


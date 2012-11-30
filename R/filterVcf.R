# filterVcf.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


filterVcf <- function(
  vcf,
  shouldRemoveInvariant       = TRUE,
#  regionsMask                 = varRegions_v2(), # will remove any variants in these regions. Set to NULL if you don't want to mask any variants out in this way
  regionsMask                 = NULL,
  keepPASSvariantsOnly        = FALSE,
#  filtersToRemove             = c("NoAlternative"),
  filtersToRemove             = NULL,
  samplesToRemove             = NULL,
  additionalInfoFilters         = NULL,
#  additionalInfoFilters     = list(
#    "SVTYPE" = list(operator="%in%", value="SNP")
#  ),
  additionalGenotypeFilters     = NULL,
#  additionalGenotypeFilters     = list(
#    "GT_CONF" = list(operator="<=", value=20),
#    "SITE_CONF" = list(operator="<=", value=200)
#  )
  possibleMissingValues       = c(".", "./.", ".|.")
) {
  if(!is.null(samplesToRemove)) {
    sampleToKeep <- setdiff(row.names(colData(vcf)), samplesToRemove)
    vcf <- vcf[, sampleToKeep]
  }
  if(shouldRemoveInvariant) {
    invariantSNPs <- apply(geno(vcf)[["GT"]], 1, function(x) length(table(x[!(x %in% possibleMissingValues)], useNA="no"))==1)
    vcf <- vcf[!invariantSNPs]
  }
  if(keepPASSvariantsOnly) {
    vcf <- vcf[filt(vcf)=="PASS"]
  } else {
    if(!is.null(filtersToRemove)) {
      sapply(
        filtersToRemove,
        function(filterToRemove) {
          vcf <<- vcf[!grepl(filterToRemove, filt(vcf))]
        }
      )
    }
  }
  if(!is.null(regionsMask)) {
    vcf <- vcf[is.na(GenomicRanges::match(rowData(vcf), regionsMask))]
  }
  if(!is.null(additionalInfoFilters)) {
    sapply(
      names(additionalInfoFilters),
      function(filterName) {
        if(additionalInfoFilters[[filterName]][["operator"]] == "%in%") {
          vcf <<- vcf[values(info(vcf))[[filterName]] %in% additionalInfoFilters[[filterName]][["value"]]]
        }
      }
    )
  }
  browser()
  if(!is.null(additionalGenotypeFilters)) {
    sapply(
      names(additionalGenotypeFilters),
      function(filterName) {
        if(additionalGenotypeFilters[[filterName]][["operator"]] == "<=") {
          geno(vcf)[["GT"]][geno(vcf)[[filterName]] <= additionalGenotypeFilters[[filterName]][["value"]]] <<- NA
        }
      }
    )
  }

  return(vcf)
}

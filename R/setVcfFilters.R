# setVcfFilters.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


setVcfFilters <- function(
  vcf,
  shouldSetInvariantFilter    = FALSE,
#  shouldSetInvariantFilter    = TRUE,
#  regionsMask                 = varRegions_v2(), # will filter any variants in these regions. Set to NULL if you don't want to mask any variants out in this way
#  regionsMask                 = varRegions_v3,
  regionsMask                 = NULL,
  regionsMaskFilterName       = "InVarRegion",
  shouldSetMultiallelicFilter = FALSE,
  markDotsAsPass              = TRUE,
  keepPASSvariantsOnly        = FALSE,
#  additionalInfoFilters       = NULL,
  additionalInfoFilters     = list(
    "LowQD" = list(column="QD", operator="<=", value=36)
  )
) {
  if(markDotsAsPass) {
    filt(vcf)[filt(vcf) == "."] <- "PASS"
  }
  if(shouldSetMultiallelicFilter) {
    multiAllelicVariants <- elementLengths(alt(vcf)) > 1
    filt(vcf)[!(filt(vcf) %in% c("PASS", ".")) & multiAllelicVariants] <- paste(filt(vcf)[!(filt(vcf) %in% c("PASS", ".")) & multiAllelicVariants], "MultiAllelic", sep=";")
    filt(vcf)[filt(vcf) %in% c("PASS", ".") & multiAllelicVariants] <- "MultiAllelic"
  }
  if(shouldSetInvariantFilter) {
    invariantSNPs <- apply(geno(vcf)[["GT"]], 1, function(x) length(table(x[!(x %in% possibleMissingValues)], useNA="no"))<=1)
    filt(vcf)[!(filt(vcf) %in% c("PASS", ".")) & invariantSNPs] <- paste(filt(vcf)[!(filt(vcf) %in% c("PASS", ".")) & invariantSNPs], "Invariant", sep=";")
    filt(vcf)[filt(vcf) %in% c("PASS", ".") & invariantSNPs] <- "Invariant"
  }
  if(!is.null(additionalInfoFilters)) {
    sapply(
      names(additionalInfoFilters),
      function(filterName) {
        if(additionalInfoFilters[[filterName]][["operator"]] == "%in%") {
          filteredVariants <- values(info(vcf))[[additionalInfoFilters[[filterName]][["column"]]]] %in% additionalInfoFilters[[filterName]][["value"]]
        } else if(additionalInfoFilters[[filterName]][["operator"]] == "<=") {
          filteredVariants <- values(info(vcf))[[additionalInfoFilters[[filterName]][["column"]]]] <= additionalInfoFilters[[filterName]][["value"]]
        } else if(additionalInfoFilters[[filterName]][["operator"]] == "<") {
          filteredVariants <- values(info(vcf))[[additionalInfoFilters[[filterName]][["column"]]]] < additionalInfoFilters[[filterName]][["value"]]
        } else if(additionalInfoFilters[[filterName]][["operator"]] == ">=") {
          filteredVariants <- values(info(vcf))[[additionalInfoFilters[[filterName]][["column"]]]] >= additionalInfoFilters[[filterName]][["value"]]
        } else if(additionalInfoFilters[[filterName]][["operator"]] == ">") {
          filteredVariants <- values(info(vcf))[[additionalInfoFilters[[filterName]][["column"]]]] > additionalInfoFilters[[filterName]][["value"]]
        }
        browser()
        filt(vcf)[!(filt(vcf) %in% c("PASS", ".")) & filteredVariants] <<- paste(filt(vcf)[!(filt(vcf) %in% c("PASS", ".")) & filteredVariants], filterName, sep=";")
        filt(vcf)[filt(vcf) %in% c("PASS", ".") & filteredVariants] <<- filterName
      }
    )
  }
  if(!is.null(regionsMask)) {
    maskedVariants <- is.na(GenomicRanges::match(rowData(vcf), regionsMask))
    filt(vcf)[!(filt(vcf) %in% c("PASS", ".")) & maskedVariants] <- paste(filt(vcf)[!(filt(vcf) %in% c("PASS", ".")) & maskedVariants], regionsMaskFilterName, sep=";")
    filt(vcf)[filt(vcf) %in% c("PASS", ".") & maskedVariants] <- regionsMaskFilterName
  }
  return(vcf)
}

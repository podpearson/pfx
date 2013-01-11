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
#  regionsMask                 = varRegions_v3(),
  regionsMask                 = NULL,
  regionsMaskFilterName       = "InVarRegion",
  shouldSetMultiallelicFilter = FALSE,
  shouldSetNonSegregatingFilt = FALSE,
  shouldSetMaxNoCallsFilter   = FALSE,
  setMonomorphicProgenyFilter = FALSE,
#  monomorphicSkipChromosomes  = "Pf3D7_13_v3", # This is necessary in 3d7_hb3 as all progeny clearly inherit from 3D7 for part of this chromosome 
  monomorphicSkipChromosomes  = NULL,
  maxNoCallsAllowed           = 1,
  markDotsAsPass              = TRUE,
  keepPASSvariantsOnly        = FALSE,
#  additionalInfoFilters       = NULL,
  additionalInfoFilters     = list(
    "LowQD" = list(column="QD", operator="<=", value=36)
  ),
  possibleMissingValues       = c(".", "./.", ".|."),
  parentalIDs                 = dimnames(vcf)[[2]][1:2]
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
        if(!is.null(additionalInfoFilters[[filterName]][["filterOutNAs"]])) {
          if(additionalInfoFilters[[filterName]][["filterOutNAs"]]) {
            filteredVariants[is.na(filteredVariants)] <- TRUE
          } else {
            filteredVariants[is.na(filteredVariants)] <- FALSE
          }
        }
#        browser()
        filt(vcf)[!(filt(vcf) %in% c("PASS", ".")) & filteredVariants] <<- paste(filt(vcf)[!(filt(vcf) %in% c("PASS", ".")) & filteredVariants], filterName, sep=";")
        filt(vcf)[filt(vcf) %in% c("PASS", ".") & filteredVariants] <<- filterName
      }
    )
  }
  if(!is.null(regionsMask)) {
    maskedVariants <- rowData(vcf) %in% regionsMask
#    maskedVariants <- is.na(GenomicRanges::match(rowData(vcf), regionsMask))
    filt(vcf)[!(filt(vcf) %in% c("PASS", ".")) & maskedVariants] <- paste(filt(vcf)[!(filt(vcf) %in% c("PASS", ".")) & maskedVariants], regionsMaskFilterName, sep=";")
    filt(vcf)[filt(vcf) %in% c("PASS", ".") & maskedVariants] <- regionsMaskFilterName
  }
  if(shouldSetNonSegregatingFilt) {
    nonSegregatingVariants <- (
      geno(vcf)[["GT"]][, parentalIDs[1]] %in% possibleMissingValues |
      geno(vcf)[["GT"]][, parentalIDs[2]] %in% possibleMissingValues |
      geno(vcf)[["GT"]][, parentalIDs[1]] == geno(vcf)[["GT"]][, parentalIDs[2]]
    )
    filt(vcf)[!(filt(vcf) %in% c("PASS", ".")) & nonSegregatingVariants] <- paste(filt(vcf)[!(filt(vcf) %in% c("PASS", ".")) & nonSegregatingVariants], "NonSegregating", sep=";")
    filt(vcf)[filt(vcf) %in% c("PASS", ".") & nonSegregatingVariants] <- "NonSegregating"
  }
  if(setMonomorphicProgenyFilter) {
    progenyIDs <- setdiff(dimnames(vcf)[[2]], parentalIDs)
    invariantSNPs <- apply(geno(vcf)[["GT"]][, progenyIDs], 1, function(x) length(table(x[!(x %in% possibleMissingValues)], useNA="no"))<=1)
    if(!is.null(monomorphicSkipChromosomes)) {
      invariantSNPs[as.character(seqnames(vcf)) %in% monomorphicSkipChromosomes] <- FALSE
    }
    filt(vcf)[!(filt(vcf) %in% c("PASS", ".")) & invariantSNPs] <- paste(filt(vcf)[!(filt(vcf) %in% c("PASS", ".")) & invariantSNPs], "Invariant", sep=";")
    filt(vcf)[filt(vcf) %in% c("PASS", ".") & invariantSNPs] <- "MonomorphicInProgeny"
  }
  if(shouldSetMaxNoCallsFilter) {
    aboveMaxNoCalls <- apply(geno(vcf)[["GT"]], 1, function(x) length(which(x %in% possibleMissingValues)) > maxNoCallsAllowed)
    filt(vcf)[!(filt(vcf) %in% c("PASS", ".")) & aboveMaxNoCalls] <- paste(filt(vcf)[!(filt(vcf) %in% c("PASS", ".")) & maskedVariants], "ExcessiveNoCalls", sep=";")
    filt(vcf)[filt(vcf) %in% c("PASS", ".") & aboveMaxNoCalls] <- "ExcessiveNoCalls"
  }
  if(keepPASSvariantsOnly) {
    vcf <- vcf[filt(vcf) == "PASS"]
  }
  return(vcf)
}

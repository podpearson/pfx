# filterGenotypes.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


filterGenotypes <- function(
  vcf,
  genotypeFilters = list(
    "LowGQ" = list(column="GQ", operator="<", value=99),
    "LowDP" = list(column="DP", operator="<", value=10)
  ),
  shouldSetFilteredGTtoMissing=TRUE,
  missingGTValue              = "."
) {
  if(!("FT" %in% names(geno(vcf)))) {
    geno(vcf)[["FT"]] <<- matrix("PASS", nrow=nrow(geno(vcf)[["GT"]]), ncol=ncol(geno(vcf)[["GT"]]), dimnames=dimnames(geno(vcf)[["GT"]]))
  }
  if(!is.null(genotypeFilters)) {
    sapply(
      names(genotypeFilters),
      function(filterName) {
        if(genotypeFilters[[filterName]][["operator"]] == "%in%") {
          filteredGenotypes <- geno(vcf)[[genotypeFilters[[filterName]][["column"]]]] %in% genotypeFilters[[filterName]][["value"]]
        } else if(genotypeFilters[[filterName]][["operator"]] == "<=") {
          filteredGenotypes <- geno(vcf)[[genotypeFilters[[filterName]][["column"]]]] <= genotypeFilters[[filterName]][["value"]]
        } else if(genotypeFilters[[filterName]][["operator"]] == "<") {
          filteredGenotypes <- geno(vcf)[[genotypeFilters[[filterName]][["column"]]]] < genotypeFilters[[filterName]][["value"]]
        } else if(genotypeFilters[[filterName]][["operator"]] == ">=") {
          filteredGenotypes <- geno(vcf)[[genotypeFilters[[filterName]][["column"]]]] >= genotypeFilters[[filterName]][["value"]]
        } else if(genotypeFilters[[filterName]][["operator"]] == ">") {
          filteredGenotypes <- geno(vcf)[[genotypeFilters[[filterName]][["column"]]]] > genotypeFilters[[filterName]][["value"]]
        }
        if(!is.null(genotypeFilters[[filterName]][["filterOutNAs"]])) {
          if(genotypeFilters[[filterName]][["filterOutNAs"]]) {
            filteredGenotypes[is.na(filteredGenotypes)] <- TRUE
          } else {
            filteredGenotypes[is.na(filteredGenotypes)] <- FALSE
          }
        }
#        browser()
        geno(vcf)[["FT"]][!(geno(vcf)[["FT"]] %in% c("PASS")) & filteredGenotypes] <<- paste(geno(vcf)[["FT"]][!(geno(vcf)[["FT"]] %in% c("PASS")) & filteredGenotypes], filterName, sep=";")
        geno(vcf)[["FT"]][geno(vcf)[["FT"]] %in% c("PASS") & filteredGenotypes] <<- filterName
#        filt(vcf)[!(filt(vcf) %in% c("PASS", ".")) & filteredVariants] <<- paste(filt(vcf)[!(filt(vcf) %in% c("PASS", ".")) & filteredVariants], filterName, sep=";")
#        filt(vcf)[filt(vcf) %in% c("PASS", ".") & filteredVariants] <<- filterName
      }
    )
  }
  if(shouldSetFilteredGTtoMissing) {
    geno(vcf)[["GT"]][geno(vcf)[["FT"]] != "PASS"] <- missingGTValue
  }
  return(vcf)
}


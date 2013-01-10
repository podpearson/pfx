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
    "LowGQ" = list(column="GQ", operator="<", value=99, filterOutNAs=TRUE),
    "LowDP" = list(column="DP", operator="<", value=10, filterOutNAs=TRUE),
    "HighMAF" = list(column="MAF", operator=">", value=0.1, filterOutNAs=TRUE)
  ),
  shouldSetFilteredGTtoMissing=FALSE,
  missingGTValue              = ".",
  shouldSetINFOcolumn         = TRUE,
  shouldAlsoSetSNPfilters     = FALSE,
  shouldRemoveFilteredSNPs    = FALSE,
  additionalInfoFilters       = list(
    "filteredGenotypes" = list(column="numFilteredGenotypes", operator=">", value=2)
  )
) {
  if(!("FT" %in% names(geno(vcf)))) {
    geno(vcf)[["FT"]] <- matrix("PASS", nrow=nrow(geno(vcf)[["GT"]]), ncol=ncol(geno(vcf)[["GT"]]), dimnames=dimnames(geno(vcf)[["GT"]]))
  }
  if(!is.null(genotypeFilters)) {
    sapply(
      names(genotypeFilters),
      function(filterName) {
        if(genotypeFilters[[filterName]][["column"]] == "DP" & !("DP" %in% names(geno(vcf))) & ("AD" %in% names(geno(vcf)))) {
          RefReads <- matrix(
            sapply(geno(vcf)[["AD"]], function(x) x[1]),
            ncol=dim(geno(vcf)[["AD"]])[2],
            dimnames=dimnames(geno(vcf)[["AD"]])
          )
          FirstAltReads <- matrix(
            sapply(geno(vcf)[["AD"]], function(x) x[2]),
            ncol=dim(geno(vcf)[["AD"]])[2],
            dimnames=dimnames(geno(vcf)[["AD"]])
          )
          AllReads <- RefReads + FirstAltReads
          AllReads[is.na(AllReads)] <- 0
          geno(vcf)[["DP"]] <<- AllReads
        }
        if(genotypeFilters[[filterName]][["column"]] == "MAF" & !("MAF" %in% names(geno(vcf))) & ("AD" %in% names(geno(vcf)))) {
          RefReads <- matrix(
            sapply(geno(vcf)[["AD"]], function(x) x[1]),
            ncol=dim(geno(vcf)[["AD"]])[2],
            dimnames=dimnames(geno(vcf)[["AD"]])
          )
          FirstAltReads <- matrix(
            sapply(geno(vcf)[["AD"]], function(x) x[2]),
            ncol=dim(geno(vcf)[["AD"]])[2],
            dimnames=dimnames(geno(vcf)[["AD"]])
          )
          MAF <- pmin(RefReads/(RefReads+FirstAltReads), FirstAltReads/(RefReads+FirstAltReads))
          geno(vcf)[["MAF"]] <<- MAF
        }
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
  if(shouldSetINFOcolumn) {
    numFilteredGenotypes <- apply(geno(vcf)[["FT"]], 1, function(x) length(which(x != "PASS" | is.na(x))))
    info(vcf) <- cbind(
      values(info(vcf)),
      DataFrame(
        numFilteredGenotypes = numFilteredGenotypes
      )
    )
    exptData(vcf)[["header"]] <- VCFHeader(
      reference=reference(exptData(vcf)[["header"]]),
      samples=samples(exptData(vcf)[["header"]]),
      header=DataFrameList(
        META=meta(exptData(vcf)[["header"]]),
        FILTER=fixed(exptData(vcf)[["header"]])[["FILTER"]],
        FORMAT=geno(exptData(vcf)[["header"]]),
        INFO=rbind(
          info(exptData(vcf)[["header"]]),
          DataFrame(Number="1", Type="Integer", Description="Number of samples that have a filtered genotype at this position", row.names="numFilteredGenotypes")
        )
      )
    )
  }
  if(shouldAlsoSetSNPfilters) {
    vcf <- setVcfFilters(vcf, additionalInfoFilters=additionalInfoFilters)
  }
  if(shouldRemoveFilteredSNPs) {
    vcf <- filterVcf(vcf, filtersToRemove = names(additionalInfoFilters))
  }
  return(vcf)
}


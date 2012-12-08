# recombinationPoints.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


recombinationPoints <- function(
  vcf,
  gffGRL                      = gffGRL,
  GTsToIntMapping             = c("0"=1, "1"=2, "."=0),
  shouldCharacterise          = TRUE
) {
#  GTsCFparents <- convertGTsIntToParentBasedGTs(genotypeCallsFromGTas012(vcf))
  chromosomeLevelResults <- sapply(
    seqlevels(vcf),
    function(chromosome) {
      cat(".")
#      GTsCFparents <- convertGTsIntToParentBasedGTs(genotypeCallsFromGTas012(vcf[seqnames(vcf)==chromosome]))
#      recombinationsGRL <- findRecombinations(GTsCFparents, chromosome)
      recombinationsGRL <- findRecombinations(vcf[seqnames(vcf)==chromosome], chromosome, gffGRL=gffGRL[seqnames(gffGRL)==chromosome], GTsToIntMapping=GTsToIntMapping, shouldCharacterise=shouldCharacterise)
      widthsList <- BiocGenerics::sapply(recombinationsGRL, width)
      list(recombinationsGRL=recombinationsGRL, widthsList=widthsList)
    },
    simplify=FALSE,
    USE.NAMES=TRUE
  )
  recombinationsList <- sapply(chromosomeLevelResults, function(x) x[["recombinationsGRL"]], USE.NAMES=TRUE, simplify=FALSE)
  recombinationsUncertainties <- sapply(
    chromosomeLevelResults,
    function(x) {
      unlist(x[["widthsList"]])
    }
  )
  sampleLevelResults <- sapply(
    names(recombinationsList[[1]]),
    function(sampleID) {
        lapply(
          recombinationsList,
          function(x) {
            x[[sampleID]]
          }
        )
    },
    simplify=FALSE,
    USE.NAMES=TRUE
  )
  recombinationsBySampleGRlist <- GRangesList(
    sapply(
      names(recombinationsList[[1]]),
      function(sampleID) {
        GenomicRanges::unlist(do.call(GRangesList, sampleLevelResults[[sampleID]][sapply(sampleLevelResults[[sampleID]], length) > 0]))
      },
      simplify=FALSE,
      USE.NAMES=TRUE
    )
  )
  haplotypeLengths <- BiocGenerics::sapply(recombinationsBySampleGRlist, function(x) width(gaps(x)))
  return(
    list(
      chromosomeLevelResults       = chromosomeLevelResults,
      recombinationsUncertainties  = recombinationsUncertainties,
      sampleLevelResults           = sampleLevelResults,
      recombinationsBySampleGRlist = recombinationsBySampleGRlist
    )
  )
}

findRecombinations <- function(
  vcf,
#  parentBasedGTs,
  chromosome,
  gffGRL                      = readGffAsGRangesList("/data/galton/mirror/nfs/team112/annotation/plasmodium/falciparum/Pfalciparum_PlasmoDB-7.2.gff"),
  GTsToIntMapping             = c("0"=1, "1"=2, "."=0),
  removeChromsomeFromRownames = FALSE,
  shouldCharacterise          = TRUE
) {
  parentBasedGTs <- convertGTsIntToParentBasedGTs(genotypeCallsFromGTas012(vcf, GTsToIntMapping=GTsToIntMapping), return0asNA=TRUE)
  markerPositions <- start(ranges(rowData(vcf)))
#  if(removeChromsomeFromRownames) {
#    markerPositions <- as.integer(sub("^[^:]*:([0-9]+)$", "\\1", dimnames(parentBasedGTs)[[1]]))
#  } else {
#    markerPositions <- as.integer(dimnames(parentBasedGTs)[[1]])
#  }
  grl <- GRangesList(
    apply(
      parentBasedGTs,
      2,
      function(x) {
#        browser()
        startIndicesOfHaplotypes <- which(c(x, NA)!=c(NA, x))
        endIndicesOfHaplotypes <- startIndicesOfHaplotypes-1
        if(length(startIndicesOfHaplotypes) == 0) {
          gr <- GRanges()
        } else {
          gr <- GRanges(
            seqnames = chromosome,
            ranges = IRanges(
              start = markerPositions[endIndicesOfHaplotypes],
              end = markerPositions[startIndicesOfHaplotypes]
            ),
            seqinfo = seqinfo(gffGRL)[chromosome]
#            uncertinaty = markerPositions[startIndicesOfHaplotypes] - markerPositions[endIndicesOfHaplotypes]
          )
        }
        if(shouldCharacterise && length(gr) > 0) {
          if("CDS" %in% names(gffGRL)) {
            gffGRL[["exons"]] <- gffGRL[["CDS"]]
          }
          values(gr)[["uncertainty"]] <- width(gr)
          values(gr)[["midpoint"]] <- start(gr) + round(width(gr)/2)
          startGR <- GRanges(seqnames=seqnames(gr), ranges=IRanges(start=start(gr), width=1))
          endGR <- GRanges(seqnames=seqnames(gr), ranges=IRanges(start=end(gr), width=1))
          startForwardGeneMatches <- GenomicRanges::match(startGR, gffGRL[["genes"]][strand(gffGRL[["genes"]])=="+"])
          startReverseGeneMatches <- GenomicRanges::match(startGR, gffGRL[["genes"]][strand(gffGRL[["genes"]])=="-"])
          endForwardGeneMatches <- GenomicRanges::match(endGR, gffGRL[["genes"]][strand(gffGRL[["genes"]])=="+"])
          endReverseGeneMatches <- GenomicRanges::match(endGR, gffGRL[["genes"]][strand(gffGRL[["genes"]])=="-"])
          startForwardExonMatches <- GenomicRanges::match(startGR, gffGRL[["exons"]][strand(gffGRL[["exons"]])=="+"])
          startReverseExonMatches <- GenomicRanges::match(startGR, gffGRL[["exons"]][strand(gffGRL[["exons"]])=="-"])
          endForwardExonMatches <- GenomicRanges::match(endGR, gffGRL[["exons"]][strand(gffGRL[["exons"]])=="+"])
          endReverseExonMatches <- GenomicRanges::match(endGR, gffGRL[["exons"]][strand(gffGRL[["exons"]])=="-"])
          values(gr)[["startForwardGene"]] <- sapply(
            startForwardGeneMatches,
            function(x) {
              if(is.na(x)) {
                ""
              } else {
                values(gffGRL[["genes"]][strand(gffGRL[["genes"]])=="+"])[x, "Name"]
              }
            }
          )
          values(gr)[["startReverseGene"]] <- sapply(
            startReverseGeneMatches,
            function(x) {
              if(is.na(x)) {
                ""
              } else {
                values(gffGRL[["genes"]][strand(gffGRL[["genes"]])=="-"])[x, "Name"]
              }
            }
          )
          values(gr)[["endForwardGene"]] <- sapply(
            endForwardGeneMatches,
            function(x) {
              if(is.na(x)) {
                ""
              } else {
                values(gffGRL[["genes"]][strand(gffGRL[["genes"]])=="+"])[x, "Name"]
              }
            }
          )
          values(gr)[["endReverseGene"]] <- sapply(
            endReverseGeneMatches,
            function(x) {
              if(is.na(x)) {
                ""
              } else {
                values(gffGRL[["genes"]][strand(gffGRL[["genes"]])=="-"])[x, "Name"]
              }
            }
          )
          values(gr)[["startForwardExon"]] <- sapply(
            startForwardExonMatches,
            function(x) {
              if(is.na(x)) {
                ""
              } else {
                values(gffGRL[["exons"]][strand(gffGRL[["exons"]])=="+"])[x, "ID"]
              }
            }
          )
          values(gr)[["startReverseExon"]] <- sapply(
            startReverseExonMatches,
            function(x) {
              if(is.na(x)) {
                ""
              } else {
                values(gffGRL[["exons"]][strand(gffGRL[["exons"]])=="-"])[x, "ID"]
              }
            }
          )
          values(gr)[["endForwardExon"]] <- sapply(
            endForwardExonMatches,
            function(x) {
              if(is.na(x)) {
                ""
              } else {
                values(gffGRL[["exons"]][strand(gffGRL[["exons"]])=="+"])[x, "ID"]
              }
            }
          )
          values(gr)[["endReverseExon"]] <- sapply(
            endReverseExonMatches,
            function(x) {
              if(is.na(x)) {
                ""
              } else {
                values(gffGRL[["exons"]][strand(gffGRL[["exons"]])=="-"])[x, "ID"]
              }
            }
          )
          if(length(gr) > 1) {
            values(gr)[["nearestMidpointLeft"]] <- numeric(length(gr))
            values(gr)[["nearestMidpointRight"]] <- numeric(length(gr))
            values(gr)[seq(2, length(gr)), "nearestMidpointLeft"] <- values(gr)[seq(1, length(gr)-1), "midpoint"]
            values(gr)[seq(1, length(gr)-1), "nearestMidpointRight"] <- values(gr)[seq(2, length(gr)), "midpoint"]
          }
          values(gr)[1, "nearestMidpointLeft"] <- as.integer(NA)
          values(gr)[length(gr), "nearestMidpointRight"] <- as.integer(NA)
        }
        return(gr)
      }
    )
  )
#  values(grl)[["uncertainty"]] <- width(grl)
#  browser()
}

#vcfMAL4 <- vcfSegregating[seqnames(vcfSegregating)=="MAL4"]
#temp <- findRecombinations(vcfMAL4, "MAL4", gffGRL[seqnames(gffGRL)=="MAL4"])
#
#vcfMAL1 <- vcfSegregating[seqnames(vcfSegregating)=="MAL1"]
#temp1 <- findRecombinations(vcfMAL1, "MAL1", gffGRL[seqnames(gffGRL)=="MAL1"])
#
#sapply(
#  names(subjectRecombinations[["chromosomeLevelResults"]]),
#  function(chromosome) {
#    sapply(
#      names(subjectRecombinations[["chromosomeLevelResults"]][[chromosome]][["recombinationsGRL"]]),
#      function(ID) {
#        class(
#          values(subjectRecombinations[["chromosomeLevelResults"]][[chromosome]][["recombinationsGRL"]][[ID]])[["nearestMidpointLeft"]]
#        )
#      }
#    )
#  }
#)

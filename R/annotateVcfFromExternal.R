# annotateVcfFromExternal.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


annotateVcfFromExternal <- function(
  vcf,
  chromosome                  = "Pf3D7_01_v3",
  externalFileDetails         = list(
    list(fileFmt="data/genome/sanger/version3/September_2012/homopolymers/%s.homopolymer_proximity.5.txt.gz", columnsInFile="prox", columnsInVcf="homopolymer5Proximity", chromColumn="chr", posColumn="pos")
#    list(fileFmt="data/genome/sanger/version3/September_2012/homopolymers/%s.homopolymer_proximity.10.txt.gz", columnInFile="prox", columnInVcf="homopolymer10Proximity"),
#    list(fileFmt="data/genome/sanger/version3/September_2012/homopolymers/%s.homopolymer_proximity.15.txt.gz", columnInFile="prox", columnInVcf="homopolymer15Proximity"),
  ),
#  regionsBedFilename          = "meta/regions_v3.bed",
#  UQfilename                  = paste("data/genome/sanger/version3/September_2012/", chromosome, ".uniqueness_to500.txt", sep=""),
#  GC100filename               = paste("data/genome/sanger/version3/September_2012/gc/", chromosome, ".gc.100.txt", sep=""),
#  GC500filename               = paste("data/genome/sanger/version3/September_2012/gc/", chromosome, ".gc.500.txt", sep=""),
#  gffFilename                 = paste("data/genome/sanger/version3/September_2012/", chromosome, ".gff.gz", sep=""),
  verbose                     = TRUE
) {
  lapply(
    externalFileDetails,
    function(currentFile) {
      hasHeader <- if(is.integer(currentFile[["chromColumn"]])) FALSE else TRUE
      valuesDF <- read.table(sprintf(currentFile[["fileFmt"]], chromosome), header=hasHeader, as.is=TRUE)
      valuesGR <- GRanges(valuesDF[[currentFile[["chromColumn"]]]], IRanges(currentFile[["posColumn"]], width=1))
      values(valuesGR) <- DataFrame(valuesDF[, currentFile[["columnsInFile"]]])
      names(values(valuesGR)) <- currentFile[["columnsInVcf"]]
      newInfo <- values(valuesGR[rowData(vcf)])
      info(vcf) <- cbind(
        values(info(vcf)),
        newInfo
      )
    }
  )
  browser()
  
  homopolymer5Proximity <- read.delim(paste("data/genome/sanger/version3/September_2012/homopolymers/", chromosome, ".homopolymer_proximity.5.txt.gz", sep=""))
#  homopolymer5Proximity[["homopolymer5Proximity"]] <- pmin(homopolymer5Proximity[["lprox"]], homopolymer5Proximity[["rprox"]], na.rm=TRUE)
  homopolymer5Proximity[["homopolymer5Proximity"]] <- homopolymer5Proximity[["prox"]]
  homopolymer5Proximity[["homopolymer5proxhombase"]] <- homopolymer5Proximity[["proxhombase"]]
  homopolymer5Proximity[["homopolymer5proxhomlength"]] <- homopolymer5Proximity[["proxhomlength"]]
  homopolymer5Proximity[["homopolymer5proxdirection"]] <- ifelse(homopolymer5Proximity[["lprox"]] == homopolymer5Proximity[["rprox"]], "e", ifelse(homopolymer5Proximity[["lprox"]] < homopolymer5Proximity[["rprox"]], "l", ifelse(homopolymer5Proximity[["lprox"]] > homopolymer5Proximity[["rprox"]], "r", "?")))
  if(verbose) {
    cat("Reading homopolymer10Proximity\n")
  }
  homopolymer10Proximity <- read.delim(paste("data/genome/sanger/version3/September_2012/homopolymers/", chromosome, ".homopolymer_proximity.10.txt.gz", sep=""))
#  homopolymer10Proximity[["homopolymer10Proximity"]] <- pmin(homopolymer10Proximity[["lprox"]], homopolymer5Proximity[["rprox"]], na.rm=TRUE)
  homopolymer10Proximity[["homopolymer10Proximity"]] <- homopolymer10Proximity[["prox"]]
  homopolymer10Proximity[["homopolymer10proxhombase"]] <- homopolymer10Proximity[["proxhombase"]]
  homopolymer10Proximity[["homopolymer10proxhomlength"]] <- homopolymer10Proximity[["proxhomlength"]]
  homopolymer10Proximity[["homopolymer10proxdirection"]] <- ifelse(homopolymer10Proximity[["lprox"]] == homopolymer10Proximity[["rprox"]], "e", ifelse(homopolymer10Proximity[["lprox"]] < homopolymer10Proximity[["rprox"]], "l", ifelse(homopolymer10Proximity[["lprox"]] > homopolymer10Proximity[["rprox"]], "r", "?")))
  if(verbose) {
    cat("Reading homopolymer15Proximity\n")
  }
  homopolymer15Proximity <- read.delim(paste("data/genome/sanger/version3/September_2012/homopolymers/", chromosome, ".homopolymer_proximity.15.txt.gz", sep=""))
#  homopolymer15Proximity[["homopolymer15Proximity"]] <- pmin(homopolymer15Proximity[["lprox"]], homopolymer5Proximity[["rprox"]], na.rm=TRUE)
  homopolymer15Proximity[["homopolymer15Proximity"]] <- homopolymer15Proximity[["prox"]]
  homopolymer15Proximity[["homopolymer15proxhombase"]] <- homopolymer15Proximity[["proxhombase"]]
  homopolymer15Proximity[["homopolymer15proxhomlength"]] <- homopolymer15Proximity[["proxhomlength"]]
  homopolymer15Proximity[["homopolymer15proxdirection"]] <- ifelse(homopolymer15Proximity[["lprox"]] == homopolymer15Proximity[["rprox"]], "e", ifelse(homopolymer15Proximity[["lprox"]] < homopolymer15Proximity[["rprox"]], "l", ifelse(homopolymer15Proximity[["lprox"]] > homopolymer15Proximity[["rprox"]], "r", "?")))
  if(verbose) {
    cat("Reading regions bed\n")
  }
  regionsBed <- read.table(regionsBedFilename, header=FALSE, sep="\t", as.is=TRUE)
  names(regionsBed) <- c("chr", "startPos", "endPos", "region")
  regionsBedThisChromosome <- subset(regionsBed, chr==chromosome)
  regionsDF <- do.call(
    rbind,
    sapply(
      seq(along=regionsBedThisChromosome[["region"]]),
      function(rowNumber) {
        data.frame(
          pos = seq(regionsBedThisChromosome[rowNumber, "startPos"] + 1, regionsBedThisChromosome[rowNumber, "endPos"]),
          region = rep(regionsBedThisChromosome[rowNumber, "region"], regionsBedThisChromosome[rowNumber, "endPos"]-regionsBedThisChromosome[rowNumber, "startPos"])
#          stringsAsFactors = FALSE
        )
      },
      simplify=FALSE
    )
  )
  if(verbose) {
    cat("Reading uniqueness\n")
  }
  UQ <- read.table(UQfilename, header=FALSE, sep="\t", as.is=TRUE)
  names(UQ) <- c("chrom", "pos", "UQ")
  if(verbose) {
    cat("Reading GC 100\n")
  }
  GC100 <- read.table(GC100filename, header=FALSE, sep="\t", as.is=TRUE)
  names(GC100) <- c("chrom", "pos", "GC100")
  if(verbose) {
    cat("Reading GC 500\n")
  }
  GC500 <- read.table(GC500filename, header=FALSE, sep="\t", as.is=TRUE)
  names(GC500) <- c("chrom", "pos", "GC500")


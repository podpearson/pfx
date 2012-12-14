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
    list(fileFmt="/data/malariagen2/plasmodium/pf-crosses/data/genome/sanger/version3/September_2012/homopolymers/%s.homopolymer_proximity.5.txt.gz", columnsInFile="prox", columnsInVcf="homopolymer5Proximity", chromColumn="chr", posColumn="pos"),
    list(fileFmt="/data/malariagen2/plasmodium/pf-crosses/data/genome/sanger/version3/September_2012/homopolymers/%s.homopolymer_proximity.10.txt.gz", columnsInFile="prox", columnsInVcf="homopolymer10Proximity", chromColumn="chr", posColumn="pos"),
    list(fileFmt="/data/malariagen2/plasmodium/pf-crosses/data/genome/sanger/version3/September_2012/homopolymers/%s.homopolymer_proximity.15.txt.gz", columnsInFile="prox", columnsInVcf="homopolymer15Proximity", chromColumn="chr", posColumn="pos"),
    list(fileFmt="/data/malariagen2/plasmodium/pf-crosses/data/genome/sanger/version3/September_2012/%s.uniqueness_to500.txt", columnsInFile=3, columnsInVcf="UQ", chromColumn=1, posColumn=2),
    list(fileFmt="/data/malariagen2/plasmodium/pf-crosses/data/genome/sanger/version3/September_2012/gc/%s.gc.500.txt", columnsInFile=3, columnsInVcf="GC500", chromColumn=1, posColumn=2)
  ),
  verbose                     = TRUE
) {
  lapply(
    externalFileDetails,
    function(currentFile) {
      hasHeader <- if(is.integer(currentFile[["chromColumn"]])) FALSE else TRUE
      if(verbose) {
        cat("Reading", sprintf(currentFile[["fileFmt"]], chromosome), "\n")
      }
      valuesDF <- read.table(sprintf(currentFile[["fileFmt"]], chromosome), header=hasHeader, sep="\t", as.is=TRUE)
      valuesGR <- GRanges(valuesDF[[currentFile[["chromColumn"]]]], IRanges(valuesDF[[currentFile[["posColumn"]]]], width=1))
      values(valuesGR) <- DataFrame(valuesDF[, currentFile[["columnsInFile"]]])
      names(values(valuesGR)) <- currentFile[["columnsInVcf"]]
      newInfo <- values(valuesGR[valuesGR %in% rowData(vcf)])
      if(verbose) {
        cat("Merging", sprintf(currentFile[["fileFmt"]], chromosome), "\n")
      }
      info(vcf) <<- cbind(
        values(info(vcf)),
        newInfo
      )
    }
  )
  vcf
}

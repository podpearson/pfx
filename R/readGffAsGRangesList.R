# readGffAsGRangesList.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


readGffAsGRangesList <- function(
  gffFilename                 = "/data/galton/mirror/nfs/team112/annotation/plasmodium/falciparum/Pfalciparum_PlasmoDB-7.2.gff",
  outputRdaFilename           = "Pfalciparum_PlasmoDB-7.2.gff.rda",
  overwriteExisting           = !file.exists(outputRdaFilename),
  chromsomeNames              = paste("MAL", 1:14, sep=""),
  verbose                     = TRUE
) {
  require(VEDA)
  if(overwriteExisting) {
    PlasmoDB <- readGffAsGRanges(gffFilename)
    seqlevels(PlasmoDB)[2:15] <- chromsomeNames
    gffGRL <- GRangesList(
      all         = PlasmoDB,
      chromosomes = PlasmoDB[values(PlasmoDB)[["type"]]=="supercontig"],
      genes       = PlasmoDB[values(PlasmoDB)[["type"]]=="gene"],
      exons       = PlasmoDB[values(PlasmoDB)[["type"]]=="exon"]
    )
    save(gffGRL, file=outputRdaFilename)
  } else {
    load(outputRdaFilename)
  }
  seqlengths(gffGRL)[as.character(seqnames(gffGRL[["chromosomes"]]))] <- end(ranges(gffGRL[["chromosomes"]]))
  return(gffGRL)
}


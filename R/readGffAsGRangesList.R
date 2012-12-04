# readGffAsGRangesList.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


readGffAsGRangesList <- function(
  gffFilename                 = "/data/galton/mirror/nfs/team112/annotation/plasmodium/falciparum/Pfalciparum_PlasmoDB-7.2.gff",
  outputRdaFilename           = paste(gffFilename, "rda", sep="."),
  overwriteExisting           = !file.exists(outputRdaFilename),
  chromsomeNames              = paste("MAL", 1:14, sep=""),
  shouldChangeChromosomeNames = grepl("PlasmoDB", gffFilename),
  verbose                     = TRUE
) {
  require(VEDA)
  if(overwriteExisting) {
    PlasmoDB <- readGffAsGRanges(gffFilename)
    if(shouldChangeChromosomeNames) {
      seqlevels(PlasmoDB)[2:15] <- chromsomeNames
    }
    gffGRL <- GRangesList(
      all           = PlasmoDB,
      chromosomes   = PlasmoDB[values(PlasmoDB)[["type"]]=="supercontig"],
      genes         = PlasmoDB[values(PlasmoDB)[["type"]]=="gene"],
      exons         = PlasmoDB[values(PlasmoDB)[["type"]]=="exon"],
      CDS           = PlasmoDB[values(PlasmoDB)[["type"]]=="CDS"],
      repeat_region = PlasmoDB[values(PlasmoDB)[["type"]]=="repeat_region"]
    )
    seqlengths2 <- sapply(
      names(seqlengths(gffGRL)),
      function(chromosome) {
        max(end(PlasmoDB[seqnames(PlasmoDB) == chromosome]))
      }
    )
    seqlengths(gffGRL) <- seqlengths2
    save(gffGRL, file=outputRdaFilename)
  } else {
    load(outputRdaFilename)
  }
#  seqlengths1 <- sapply(
#    names(seqlengths(gffGRL)),
#    function(chromosome) {
#      max(end(gffGRL[["repeat_region"]][seqnames(gffGRL[["repeat_region"]]) == chromosome]))
#    }
#  )
#  seqlengths(gffGRL)[as.character(seqnames(gffGRL[["chromosomes"]]))] <- end(ranges(gffGRL[["chromosomes"]]))
  return(gffGRL)
}


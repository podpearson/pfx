# varRegionsBedFile.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


varRegionsBedFile <- function(
  gffFilename                 = "/data/galton/mirror/nfs/team112/annotation/plasmodium/falciparum/Pfalciparum_PlasmoDB-7.2.gff",
  chromsomeNames              = paste("MAL", 1:14, sep="")
) {
  require(VEDA)
  PlasmoDB <- readGffAsGRanges(gffFilename)
  seqlevels(PlasmoDB)[2:15] <- chromsomeNames
  PlasmoDB_genes <- PlasmoDB[values(PlasmoDB)[["type"]]=="gene"]
  PlasmoDB_genes[seqnames(PlasmoDB_genes)=="MAL4"][order(start(ranges(PlasmoDB_genes[seqnames(PlasmoDB_genes)=="MAL4"])))][, c("Name", "description", "size", "web_id", "Alias")][1:18]
  allGenesSubset <- function(PlasmoDB_genes, chromosome="MAL1") {
    PlasmoDB_genes[
      seqnames(PlasmoDB_genes)==chromosome
    ][
      order(
        start(
          ranges(
            PlasmoDB_genes[
              seqnames(PlasmoDB_genes)==chromosome
            ]
          )
        )
      )
    ][, c("Name", "description", "size", "web_id", "Alias")]
  }
  varRifStevorSubset <- function(
    PlasmoDB_genes,
    chromosome                = "MAL1",
    removePseudogenes         = TRUE
  ) {
    genesToReturn <- PlasmoDB_genes[
      seqnames(PlasmoDB_genes)==chromosome
    ][
      order(
        start(
          ranges(
            PlasmoDB_genes[
              seqnames(PlasmoDB_genes)==chromosome
            ]
          )
        )
      )
    ][
      grep(
        "PfEMP1|rifin|stevor",
#        "VAR|RIF|stevor",
        values(
          PlasmoDB_genes[
            seqnames(PlasmoDB_genes)==chromosome
          ][
            order(
              start(
                ranges(
                  PlasmoDB_genes[
                    seqnames(PlasmoDB_genes)==chromosome
                  ]
                )
              )
            )
          ]
        )[["description"]]
#        )[["Alias"]]
      )
    ][, c("Name", "description", "size", "web_id", "Alias")]
    if(removePseudogenes) {
      genesToReturn <- genesToReturn[-grep("pseudogene", values(genesToReturn)[["description"]])]
    }
    return(genesToReturn)
  }
  MAL1 <- varRifStevorSubset(PlasmoDB_genes, "MAL1")
  

}

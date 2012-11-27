# varRegions_v2.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################

#  see varRegionsBedFile.R for details of how these were determined
# e.g. varRifStevorSubset(PlasmoDB_genes, "MAL1")
# e.g. varRifStevorSubset(PlasmoDB_genes, "MAL6")
# e.g. allGenesSubset(PlasmoDB_genes, "MAL6")[11:29]
# end points of chromosomes determined with PlasmoDB[3:16][, "Name"]
# Note that we have excluded the following 2 single internal genes:
# MAL4 [165737, 166772]
# MAL12 [ 766647,  774190]

varRegions_v2 <- function() {
  require(GenomicRanges)
  c(
    GRanges(seqnames = "MAL1",  ranges = IRanges(start = 1,       end = 91653  ), region="S"),
    GRanges(seqnames = "MAL1",  ranges = IRanges(start = 565425 , end = 643292 ), region="S"),
    GRanges(seqnames = "MAL2",  ranges = IRanges(start = 1,       end = 67545  ), region="S"),
    GRanges(seqnames = "MAL2",  ranges = IRanges(start = 860464 , end = 947102 ), region="S"),
    GRanges(seqnames = "MAL3",  ranges = IRanges(start = 1,       end = 65508  ), region="S"),
    GRanges(seqnames = "MAL3",  ranges = IRanges(start = 1012465, end = 1060087), region="S"),
    GRanges(seqnames = "MAL4",  ranges = IRanges(start = 1,       end = 100256 ), region="S"),
    GRanges(seqnames = "MAL4",  ranges = IRanges(start = 552884 , end = 617416 ), region="I"),
    GRanges(seqnames = "MAL4",  ranges = IRanges(start = 939470 , end = 987434 ), region="I"),
    GRanges(seqnames = "MAL4",  ranges = IRanges(start = 1150294, end = 1204112), region="S"),
    GRanges(seqnames = "MAL5",  ranges = IRanges(start = 1,       end = 37564  ), region="S"),
    GRanges(seqnames = "MAL5",  ranges = IRanges(start = 1322395, end = 1343552), region="S"),
    GRanges(seqnames = "MAL6",  ranges = IRanges(start = 1,       end = 33244  ), region="S"),
    GRanges(seqnames = "MAL6",  ranges = IRanges(start = 723114 , end = 742403 ), region="I"),
    GRanges(seqnames = "MAL6",  ranges = IRanges(start = 1329774, end = 1418244), region="S"),
    GRanges(seqnames = "MAL7",  ranges = IRanges(start = 1,       end = 132431 ), region="S"),
    GRanges(seqnames = "MAL7",  ranges = IRanges(start = 567328 , end = 658759 ), region="I"),
    GRanges(seqnames = "MAL7",  ranges = IRanges(start = 1441011, end = 1501717), region="S"),
    GRanges(seqnames = "MAL8",  ranges = IRanges(start = 1,       end = 58827  ), region="S"),
    GRanges(seqnames = "MAL8",  ranges = IRanges(start = 432287 , end = 467688 ), region="I"),
    GRanges(seqnames = "MAL8",  ranges = IRanges(start = 1349128, end = 1419563), region="S"),
    GRanges(seqnames = "MAL9",  ranges = IRanges(start = 1,       end = 78874  ), region="S"),
    GRanges(seqnames = "MAL9",  ranges = IRanges(start = 1474419, end = 1541723), region="S"),
    GRanges(seqnames = "MAL10", ranges = IRanges(start = 1,       end = 61366  ), region="S"),
    GRanges(seqnames = "MAL10", ranges = IRanges(start = 1599899, end = 1687655), region="S"),
    GRanges(seqnames = "MAL11", ranges = IRanges(start = 1,       end = 74442  ), region="S"),
    GRanges(seqnames = "MAL11", ranges = IRanges(start = 2006703, end = 2038337), region="S"),
    GRanges(seqnames = "MAL12", ranges = IRanges(start = 1,       end = 56805  ), region="S"),
    GRanges(seqnames = "MAL12", ranges = IRanges(start = 1694139, end = 1743396), region="I"),
    GRanges(seqnames = "MAL12", ranges = IRanges(start = 2190467, end = 2271478), region="S"),
    GRanges(seqnames = "MAL13", ranges = IRanges(start = 1,       end = 63650  ), region="S"),
    GRanges(seqnames = "MAL13", ranges = IRanges(start = 2826695, end = 2895605), region="S"),
    GRanges(seqnames = "MAL14", ranges = IRanges(start = 1,       end = 28429  ), region="S"),
    GRanges(seqnames = "MAL14", ranges = IRanges(start = 3266071, end = 3291871), region="S")
  )
}

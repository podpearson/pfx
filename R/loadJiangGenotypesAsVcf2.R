# loadJiangGenotypesAsVcf2.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


loadJiangGenotypesAsVcf2 <- function(
  jiangGenotypesFilename      = "data/jiang_et_al/gb-2011-12-4-r33-s3.txt",
#  jiangGenotypesFilename      = "/data/ib/users/rpearson/malariagen/PfCrossesPaper/externalData/gb-2011-12-4-r33-s3.txt",
  jiangVcfFilename            = "data/jiang_et_al/gb-2011-12-4-r33-s3.v3.vcf",
  jiangRdaFilename            = "data/jiang_et_al/gb-2011-12-4-r33-s3.v3.vcf.rda",
  jiangV2BedFilename          = "data/jiang_et_al/gb-2011-12-4-r33-s3.bed",
  jiangV2FixedBedFilename     = "data/jiang_et_al/gb-2011-12-4-r33-s3.fixed.bed",
  jiangV3BedFilename          = "data/jiang_et_al/gb-2011-12-4-r33-s3.v3.bed",
  jiangV3UnmappedFilename     = "data/jiang_et_al/gb-2011-12-4-r33-s3.v3.unmapped",
  v2tov3chainFile             = "data/jiang_et_al/2to3.liftOver",
  shouldSaveVcfFile           = TRUE, # Note that due to a bug in VariantAnnotation::.makeVcfGeno it is necessary to add "DUMMY" genotypes data
  shouldSaveRdaFile           = TRUE
) {
  require(rtracklayer)
  jiangGenotypes <- read.table(jiangGenotypesFilename, skip=1, header=TRUE, sep="\t", as.is=TRUE, nrows=3452)
  jiangGenotypes <- subset(jiangGenotypes, Marker!="centromere")
  jiangGenotypes <- subset(jiangGenotypes, substr(Marker, 1, 2)=="Pf")
  jiangGenotypes[["pos"]] <- as.integer(as.numeric(jiangGenotypes[["Position..Kb."]])*1000)
#  jiangGenotypes[["pos"]] <- as.integer(sub("\\.", "", jiangGenotypes[["Position..Kb."]]))
  row.names(jiangGenotypes) <- paste(sprintf("Pf3D7_%02d", jiangGenotypes[["Chr"]]), ":", jiangGenotypes[["pos"]], sep="")
#  row.names(jiangGenotypes) <- paste("MAL", jiangGenotypes[["Chr"]], ":", jiangGenotypes[["pos"]], sep="")
  GT = cbind(
    matrix("7", nrow=dim(jiangGenotypes)[1], ncol=1, dimnames=list(dimnames(jiangGenotypes)[[1]], "7G8")),
    matrix("G", nrow=dim(jiangGenotypes)[1], ncol=1, dimnames=list(dimnames(jiangGenotypes)[[1]], "GB4")),
    as.matrix(jiangGenotypes[, 6:37])
  )
#  GT <- matrix(
#    as.integer(as.factor(GT))-1,
#    ncol=ncol(GT),
#    dimnames=dimnames(GT)
#  )
#  Couldn't get the following to work...
#  chain <- import.chain(v2tov3chainFile)
#  rowDataLiftedOver <- liftOver(rowData, chain)
#  browser()
  jiangV2Bed <- read.table(jiangV2BedFilename, header=FALSE, sep="\t", col.names=c("Chr", "Start", "End", "ID"), as.is=TRUE)
  jiangV2Bed[["Start"]] <- jiangV2Bed[["Start"]]-1
  write.table(jiangV2Bed, file=jiangV2FixedBedFilename, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  system(paste("data/jiang_et_al/liftOver", jiangV2FixedBedFilename, v2tov3chainFile, jiangV3BedFilename, jiangV3UnmappedFilename))  
  jiangV3Bed <- read.table(jiangV3BedFilename, header=FALSE, sep="\t", col.names=c("Chr", "Start", "End", "ID"), as.is=TRUE)

#  couple of sanity checks. Firstly are all IDs still the same (they should be), and secondly have any variants changed chromosomes (only a few moved from 8 to 7)
  identical(jiangV2Bed[["ID"]], jiangV3Bed[["ID"]])
  table(jiangV2Bed[["Chr"]], jiangV3Bed[["Chr"]])
  
  jiangRowData  = GRanges(
    seqnames = jiangV3Bed[["Chr"]],
    ranges   = IRanges(
      start = jiangV3Bed[["End"]],
      width = 1,
      names  = jiangV3Bed[["ID"]]
    )
  )

  jiangVcf <- VCF(
    rowData = jiangRowData,
#    rowData  = GRanges(
#      seqnames = paste("MAL", jiangGenotypes[["Chr"]], sep=""),
#      ranges   = IRanges(
#        start = as.integer(jiangGenotypes[["pos"]]),
#        width = 1,
#        names  = row.names(jiangGenotypes)
#      )
##      paramRangeID = factor(rep(NA, dim(jiangGenotypes)[1]))
#    ),
    colData  = BiocGenerics::rbind(
      DataFrame(
        Samples=1:2,
        row.names=c("_7G8", "GB4")
      ),
      DataFrame(
        Samples=seq(3, dim(jiangGenotypes[, 6:37])[2]+2),
        row.names = dimnames(jiangGenotypes)[[2]][6:37]
      )
    ),
    exptData = SimpleList(
      header = VCFHeader(
        samples = c("_7G8", "GB4", dimnames(jiangGenotypes)[[2]][6:37]),
        header  = DataFrameList(
          META = rbind(
            DataFrame(Value = "VCFv4.0", row.names="fileformat"),
            DataFrame(Value = "JiangEtAl", row.names="ProjectName"),
            DataFrame(Value = "_7G8", row.names="PARENT"),
            DataFrame(Value = "GB4", row.names="PARENT.1")
          ),
          FILTER = DataFrame(Descrption=character()),
#          FILTER = DataFrame(Descrption="PASS", row.names="PASS"),
          FORMAT = rbind(
            DataFrame(Number = "1", Type="String", Description="Genotype (7 means matches 7G8, G means matches GB4)", row.names="GT"),
            DataFrame(Number = "2", Type="String", Description="Genotype (7 means matches 7G8, G means matches GB4)", row.names="AD")
          ),
          INFO = rbind(
            DataFrame(Number = "1", Type = "Float", Description="Marker.Distance..Kb.", row.names="Marker.Distance..Kb."),
            DataFrame(Number = "1", Type = "Integer", Description="Crossovers", row.names="Crossovers"),
            DataFrame(Number = ".", Type = "Float", Description="Flanking.Size..Kb.", row.names="Flanking.Size..Kb."),
            DataFrame(Number = ".", Type = "String", Description="Valedated.Genotype", row.names="Valedated.Genotype"),
            DataFrame(Number = ".", Type = "String", Description="Forward.Primer", row.names="Forward.Primer"),
            DataFrame(Number = ".", Type = "String", Description="Reverse.Primer", row.names="Reverse.Primer")
          )
        )
      )
    ),
    fixed    = DataFrame(
      REF = DNAStringSet(rep("A", dim(jiangGenotypes)[1])),
      ALT = DNAStringSetList("C"),
#      ALT = DNAStringSetList(rep("C", dim(jiangGenotypes)[1])),
      QUAL = rep(0.0, dim(jiangGenotypes)[1]),
      FILTER = rep("PASS", dim(jiangGenotypes)[1])
    ),
    info     = DataFrame(
      jiangGenotypes[, c(4:5, 38:41)]
    ),
    geno     = SimpleList(
      GT = GT,
      AD = matrix(sapply(as.vector(GT), function(x) list(c(x, x))), ncol=ncol(GT), dimnames=dimnames(GT))
#      GT = cbind(
#        matrix("7", nrow=dim(jiangGenotypes)[1], ncol=1, dimnames=list(dimnames(jiangGenotypes)[[1]], "7G8")),
#        matrix("G", nrow=dim(jiangGenotypes)[1], ncol=1, dimnames=list(dimnames(jiangGenotypes)[[1]], "GB4")),
#        as.matrix(jiangGenotypes[, 6:37])
#      ),
#      GT2 = cbind(
#        matrix("7", nrow=dim(jiangGenotypes)[1], ncol=1, dimnames=list(dimnames(jiangGenotypes)[[1]], "7G8")),
#        matrix("G", nrow=dim(jiangGenotypes)[1], ncol=1, dimnames=list(dimnames(jiangGenotypes)[[1]], "GB4")),
#        as.matrix(jiangGenotypes[, 6:37])
#      )
    )
  )
#  genome(jiangVcf) <- "Pf"
  jiangVcf <- jiangVcf[order(rowData(jiangVcf))]
  if(shouldSaveRdaFile) {
    save(jiangVcf, file=jiangRdaFilename)
  }
  if(shouldSaveVcfFile) {
    writeVcf(jiangVcf, jiangVcfFilename, index=TRUE)
#    writeVcf(jiangVcf, jiangVcfFilename)
  }
  return(jiangVcf)
}



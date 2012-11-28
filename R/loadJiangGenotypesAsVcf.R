# loadJiangGenotypesAsVcf.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


loadJiangGenotypesAsVcf <- function(
  jiangGenotypesFilename      = "/data/malariagen2/users/rpearson/pfCrosses/externalData/gb-2011-12-4-r33-s3.txt",
#  jiangGenotypesFilename      = "/data/ib/users/rpearson/malariagen/PfCrossesPaper/externalData/gb-2011-12-4-r33-s3.txt",
  jiangVcfFilename            = "/data/malariagen2/users/rpearson/pfCrosses/externalData/gb-2011-12-4-r33-s3.vcf",
  jiangRdaFilename            = "/data/malariagen2/users/rpearson/pfCrosses/externalData/gb-2011-12-4-r33-s3.vcf.rda",
  shouldSaveVcfFile           = FALSE,
  shouldSaveRdaFile           = FALSE
) {
  jiangGenotypes <- read.table(jiangGenotypesFilename, skip=1, header=TRUE, sep="\t", as.is=TRUE, nrows=3452)
  jiangGenotypes <- subset(jiangGenotypes, Marker!="centromere")
  jiangGenotypes <- subset(jiangGenotypes, substr(Marker, 1, 2)=="Pf")
  jiangGenotypes[["pos"]] <- as.integer(sub("\\.", "", jiangGenotypes[["Position..Kb."]]))
  row.names(jiangGenotypes) <- paste("MAL", jiangGenotypes[["Chr"]], ":", jiangGenotypes[["pos"]], sep="")
#  GT = cbind(
#    matrix("7", nrow=dim(jiangGenotypes)[1], ncol=1, dimnames=list(dimnames(jiangGenotypes)[[1]], "7G8")),
#    matrix("G", nrow=dim(jiangGenotypes)[1], ncol=1, dimnames=list(dimnames(jiangGenotypes)[[1]], "GB4")),
#    as.matrix(jiangGenotypes[, 6:37])
#  )
#  GT <- matrix(
#    as.integer(as.factor(GT))-1,
#    ncol=ncol(GT),
#    dimnames=dimnames(GT)
#  )
  jiangVcf <- VCF(
    rowData  = GRanges(
      seqnames = paste("MAL", jiangGenotypes[["Chr"]], sep=""),
      ranges   = IRanges(
        start = as.integer(jiangGenotypes[["pos"]]),
        width = 1,
        names  = row.names(jiangGenotypes)
      )
#      paramRangeID = factor(rep(NA, dim(jiangGenotypes)[1]))
    ),
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
#          FILTER = DataFrame(Descrption=character()),
          FILTER = DataFrame(Descrption="PASS", row.names="PASS"),
          FORMAT = rbind(
            DataFrame(Number = "1", Type="String", Description="Genotype (7 means matches 7G8, G means matches GB4)", row.names="GT")
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
      ALT = DNAStringSetList(rep("C", dim(jiangGenotypes)[1])),
      QUAL = rep(0.0, dim(jiangGenotypes)[1]),
      FILTER = rep("PASS", dim(jiangGenotypes)[1])
    ),
    info     = DataFrame(
      jiangGenotypes[, c(4:5, 38:41)]
    ),
    geno     = SimpleList(
#      GT = GT
      GT = cbind(
        matrix("7", nrow=dim(jiangGenotypes)[1], ncol=1, dimnames=list(dimnames(jiangGenotypes)[[1]], "7G8")),
        matrix("G", nrow=dim(jiangGenotypes)[1], ncol=1, dimnames=list(dimnames(jiangGenotypes)[[1]], "GB4")),
        as.matrix(jiangGenotypes[, 6:37])
      )
    )
  )
  genome(jiangVcf) <- "Pf"
  if(shouldSaveVcfFile) {
    writeVcf(jiangVcf, jiangVcfFilename, index=TRUE)
  }
  if(shouldSaveRdaFile) {
    save(jiangVcf, file=jiangRdaFilename)
  }
  return(jiangVcf)
}


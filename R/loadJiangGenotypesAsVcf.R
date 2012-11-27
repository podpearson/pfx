# loadJiangGenotypesAsVcf.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


loadJiangGenotypesAsVcf <- function(
  jiangGenotypesFilename      = "/data/malariagen2/users/rpearson/pfCrosses/externalData/gb-2011-12-4-r33-s3.txt"
#  jiangGenotypesFilename      = "/data/ib/users/rpearson/malariagen/PfCrossesPaper/externalData/gb-2011-12-4-r33-s3.txt"
) {
  jiangGenotypes <- read.table(jiangGenotypesFilename, skip=1, header=TRUE, sep="\t", as.is=TRUE, nrows=3452)
  jiangGenotypes <- subset(jiangGenotypes, Marker!="centromere")
  jiangGenotypes <- subset(jiangGenotypes, substr(Marker, 1, 2)=="Pf")
  jiangGenotypes[["pos"]] <- as.integer(sub("\\.", "", jiangGenotypes[["Position..Kb."]]))
  row.names(jiangGenotypes) <- paste("MAL", jiangGenotypes[["Chr"]], ":", jiangGenotypes[["pos"]], sep="")
  jiangVcf <- VCF(
    rowData  = GRanges(
      seqnames = paste("MAL", jiangGenotypes[["Chr"]], sep=""),
      ranges   = IRanges(
        start = as.integer(jiangGenotypes[["pos"]]),
        width = 1,
        names  = row.names(jiangGenotypes)
      )
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
#    exptData = exptData(vcfList[[1]]),
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
      GT = cbind(
        matrix("7", nrow=dim(jiangGenotypes)[1], ncol=1, dimnames=list(dimnames(jiangGenotypes)[[1]], "7G8")),
        matrix("G", nrow=dim(jiangGenotypes)[1], ncol=1, dimnames=list(dimnames(jiangGenotypes)[[1]], "GB4")),
        as.matrix(jiangGenotypes[, 6:37])
      )
    )
  )
  return(jiangVcf)
}


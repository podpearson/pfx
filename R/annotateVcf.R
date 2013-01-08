# annotateVcf.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


annotateVcf <- function(
  vcf,
  parentalIDs                 = dimnames(vcf)[[2]][1:2]
) {
  GTsInt <- genotypeCallsFromGTas012(vcf)
  columnIndexesOfParents <- match(parentalIDs, dimnames(GTsInt)[[2]])
  columnIndexesOfProgeny <- setdiff(seq(along=dimnames(GTsInt)[[2]]), match(parentalIDs, dimnames(GTsInt)[[2]]))
#  segregating <- (
#    (GTsInt[, parentalIDs[1]] == 1 & GTsInt[, parentalIDs[2]] == 2) |
#    (GTsInt[, parentalIDs[1]] == 2 & GTsInt[, parentalIDs[2]] == 1)
#  )
  ADsArray <- array(
    unlist(geno(vcf)[["AD"]]),
    dim=c(2, dim(geno(vcf)[["AD"]])[1], dim(geno(vcf)[["AD"]])[2]),
    dimnames=list(c("Ref", "Nonref"), dimnames(geno(vcf)[["AD"]])[[1]], dimnames(geno(vcf)[["AD"]])[[2]])
#  Was originally doing this as below which is incorrect
#    dim=c(2, dim(geno(vcf)[["AD"]])[2], dim(geno(vcf)[["AD"]])[1]),
#    dimnames=list(c("Ref", "Nonref"), dimnames(geno(vcf)[["AD"]])[[2]], dimnames(geno(vcf)[["AD"]])[[1]])
  )
  RefReads <- matrix(
    sapply(geno(vcf)[["AD"]], function(x) x[1]),
    ncol=dim(geno(vcf)[["AD"]])[2],
    dimnames=dimnames(geno(vcf)[["AD"]])
  )
  FirstAltReads <- matrix(
    sapply(geno(vcf)[["AD"]], function(x) x[2]),
    ncol=dim(geno(vcf)[["AD"]])[2],
    dimnames=dimnames(geno(vcf)[["AD"]])
  )
  MAF <- pmin(RefReads, FirstAltReads)/(RefReads+FirstAltReads)
#  MAF <- pmin(ADsArray[1,,], ADsArray[2,,])/(ADsArray[1,,]+ADsArray[2,,])
  meanVariantMAFs <- rowMeans(MAF, na.rm = TRUE)
  maxVariantMAFs <- apply(MAF, 1, function(x) max(x, na.rm=TRUE))
  parent1MAF <- MAF[, parentalIDs[1]]
  parent2MAF <- MAF[, parentalIDs[2]]
  maxParentMAF <- pmax(parent1MAF, parent2MAF, na.rm=TRUE)
  
  missingness <- apply(GTsInt, 1, function(x) length(which(x==0)))
  
  ADas0123 <- genotypeCallsFromADas0123(vcf)
  missingnessPerVariant <- apply(ADas0123, 1, function(x) length(which(is.na(x) | x==0)))
  heterozygosityPerVariant <- apply(ADas0123, 1, function(x) length(which(x==3)))
  
  AllReads <- RefReads + FirstAltReads
  AllReads[is.na(AllReads)] <- 0
#  depthsPerVariant <- rowSums(AllReads)
#  
#  GCbins <- round(values(info(vcf))[["GC500"]])
#  GCbins[GCbins<10] <- 9
#  GCbins[GCbins>40] <- 41
#  
#  medianDepthPerGCbin <- by(depthsPerVariant, GCbins, median)
#  plotDF <- data.frame(
#    GCbin = GCbins,
#    log10depth = log10(depthsPerVariant)
#  )
#  require(ggplot2)
#  pdf("~/log10depthByGCbin.pdf", height=12, width=8)
#  print(qplot(log10depth, facets=GCbin~., fill=GCbin, data=plotDF) + theme_bw())
#  dev.off()
#  
#  HighQD <- values(info(vcf))[["QD"]] > 36
##  GCbins_HighQD <- round(values(info(vcf))[["GC500"]][HighQD])
#  GCbins_HighQD <- GCbins[HighQD]
#  depthsPerVariant_HighQD <- depthsPerVariant[HighQD]
#  depthsPerVariant_HighQD <- values(info(vcf))[["DP"]][HighQD]
#  plotDF_HighQD <- data.frame(
#    GCbin = GCbins_HighQD,
#    log10depth = log10(depthsPerVariant_HighQD)
#  )
#  require(ggplot2)
#  pdf("~/log10depthByGCbin_HighQD.pdf", height=12, width=8)
#  print(qplot(log10depth, facets=GCbin~., fill=GCbin, data=plotDF_HighQD) + theme_bw())
#  dev.off()
#  medianDepthPerGCbin_HighQD <- by(depthsPerVariant_HighQD, GCbins_HighQD, median)
#  
#  scaledDepths <- t(scale(t(AllReads)))
  scaledDepths <- scale(AllReads)
#  log10depthSDs <- log10(apply(scaledDepths, 1, sd))
#  chr11 <- as.logical(seqnames(vcf)=="Pf3D7_11_v3" & HighQD)
#  pdf("~/chr11_log10depthSDs_heatmap.pdf", height=8, width=12)
#  print(qplot(start(rowData(vcf[chr11])), log10depthSDs[chr11]) + stat_bin2d() + theme_bw())
#  dev.off()
#  pdf("~/chr11_log10depthSDs.pdf", height=8, width=12)
#  print(qplot(start(rowData(vcf[chr11])), log10depthSDs[chr11]) + theme_bw())
#  dev.off()
  
  scaledDepthSD <- apply(scaledDepths, 1, sd)
  
  
  
# debugging stuff - ignore
#  ADsArray[,1,2]
#  geno(vcf)[["AD"]][1,2]
#  ADsArray[,2,1]
#  geno(vcf)[["AD"]][2,1]
#  geno(vcf)[["AD"]][26,1]
#  geno(vcf)[["AD"]][27,1]
#  ADsArray[,26,1]
#  ADsArray[,27,1]
#  RefReads[26,1]
#  RefReads[27,1]
#  FirstAltReads[26,1]
#  FirstAltReads[27,1]
#  v2 <- unlist(geno(vcf)[["AD"]][2,])
#  ADsTemp <- array(
#    unlist(geno(vcf)[["AD"]][2, ]),
#    dim=c(2, 1, dim(geno(vcf)[["AD"]])[1]),
#    dimnames=list(c("Ref", "Nonref"), dimnames(geno(vcf)[["AD"]])[[2]][2], dimnames(geno(vcf)[["AD"]])[[1]])
#  )
#  browser()
  info(vcf) <- cbind(
    values(info(vcf)),
    DataFrame(
      meanMAF = meanVariantMAFs,
      maxMAF = maxVariantMAFs,
      parent1MAF = parent1MAF,
      parent2MAF = parent2MAF,
      maxParentMAF = maxParentMAF,
      missingness = missingness,
      missingness2 = missingnessPerVariant,
      heterozgosity = heterozygosityPerVariant,
      scaledDepthSD = scaledDepthSD,
      ProperPair = values(info(vcf))[["DPProperPair"]] / values(info(vcf))[["DPAll"]],
      MateUnmapped = values(info(vcf))[["DPMateUnmapped"]] / values(info(vcf))[["DPAll"]],
      MateOtherChrom = values(info(vcf))[["DPMateOtherChrom"]] / values(info(vcf))[["DPAll"]],
      MateSameStrand = values(info(vcf))[["DPMateSameStrand"]] / values(info(vcf))[["DPAll"]],
      FaceAway = values(info(vcf))[["DPFaceAway"]] / values(info(vcf))[["DPAll"]],
      SoftClipped = values(info(vcf))[["DPSoftClipped"]] / values(info(vcf))[["DPAll"]],
      RepeatCopies1 = sapply(values(info(vcf))[["RepeatCopies"]], function(x) if(is.na(x)) 0 else x[1]),
      RepeatPeriod1 = sapply(values(info(vcf))[["RepeatPeriod"]], function(x) if(is.na(x)) 0 else x[1]),
      RepeatScore1 = sapply(values(info(vcf))[["RepeatScore"]], function(x) if(is.na(x)) 0 else x[1]),
      RepeatSize1 = sapply(values(info(vcf))[["RepeatSize"]], function(x) if(is.na(x)) 0 else x[1]),
      RepeatEntropy1 = sapply(values(info(vcf))[["RepeatEntropy"]], function(x) if(is.na(x)) 0 else x[1]),
      QUAL = qual(vcf),
      QUALbyDP = qual(vcf)/values(info(vcf))[["DP"]],
      QUALperSample = qual(vcf)/dim(vcf)[2],
      SEGREGATING=(
        (GTsInt[, parentalIDs[1]] == 1 & GTsInt[, parentalIDs[2]] == 2) |
        (GTsInt[, parentalIDs[1]] == 2 & GTsInt[, parentalIDs[2]] == 1)
      ),
      MendelianErrors=apply(
        GTsInt,
        1,
        function(genotypesForVariant) {
          sum(
            genotypesForVariant[parentalIDs[1]] == 1 &
              genotypesForVariant[parentalIDs[2]] == 1 &
              genotypesForVariant[columnIndexesOfProgeny] != 1
          ) +
          sum(
            genotypesForVariant[parentalIDs[1]] == 2 &
              genotypesForVariant[parentalIDs[2]] == 2 &
              genotypesForVariant[columnIndexesOfProgeny] != 2
          )
        }
      )
    )
  )
  exptData(vcf)[["header"]] <- VCFHeader(
    reference=reference(exptData(vcf)[["header"]]),
    samples=samples(exptData(vcf)[["header"]]),
    header=DataFrameList(
      META=meta(exptData(vcf)[["header"]]),
      FILTER=fixed(exptData(vcf)[["header"]])[["FILTER"]],
#      fixed=fixed(exptData(vcf)[["header"]]),
      FORMAT=geno(exptData(vcf)[["header"]]),
      INFO=rbind(
        info(exptData(vcf)[["header"]]),
        DataFrame(Number="1", Type="Float", Description="Mean across samples of proportion of minor allele reads (something like a heterozygosity score)", row.names="meanMAF"),
        DataFrame(Number="1", Type="Float", Description="Maximum across samples of proportion of minor allele reads (something like a heterozygosity score)", row.names="maxMAF"),
        DataFrame(Number="1", Type="Float", Description="Proportion of minor allele reads for first parent (something like a heterozygosity score)", row.names="parent1MAF"),
        DataFrame(Number="1", Type="Float", Description="Proportion of minor allele reads for second parent (something like a heterozygosity score)", row.names="parent2MAF"),
        DataFrame(Number="1", Type="Float", Description="Maximum proportion of minor allele reads in parents (something like a heterozygosity score)", row.names="maxParentMAF"),
        DataFrame(Number="1", Type="Integer", Description="Number of samples with a missing genotype call", row.names="missingness"),
        DataFrame(Number="1", Type="Integer", Description="Number of samples with zero depth", row.names="missingness2"),
        DataFrame(Number="1", Type="Integer", Description="Number of samples with at least 2 ref and 2 alt reads, and at least 5 reads in total", row.names="heterozgosity"),
        DataFrame(Number="1", Type="Integer", Description="Standard deviation of the depth across samples after depths have been normalised (mean 0, sd 1) by sample. Higher values suggest copy number differences between samples.", row.names="scaledDepthSD"),
        DataFrame(Number="1", Type="Float", Description="Proportion of reads than are in a proper pair", row.names="ProperPair"),
        DataFrame(Number="1", Type="Float", Description="Proportion of reads that have an unmapped mate", row.names="MateUnmapped"),
        DataFrame(Number="1", Type="Float", Description="Proportion of reads that have a mate on another chromosome", row.names="MateOtherChrom"),
        DataFrame(Number="1", Type="Float", Description="Proportion of reads that have a mate on the same strand", row.names="MateSameStrand"),
        DataFrame(Number="1", Type="Float", Description="Proportion of reads that have a faceaway mate", row.names="FaceAway"),
        DataFrame(Number="1", Type="Float", Description="Proportion of reads than are soft clipped", row.names="SoftClipped"),
        DataFrame(Number=".", Type="Float", Description="Number of copies aligned with the first consensus pattern", row.names="RepeatCopies1"),
        DataFrame(Number=".", Type="Integer", Description="Period size of the first repeat", row.names="RepeatPeriod1"),
        DataFrame(Number=".", Type="Integer", Description="Alignment score of the first repeat", row.names="RepeatScore1"),
        DataFrame(Number=".", Type="Integer", Description="Size of consensus pattern of the first repeat (may differ slightly from the period size)", row.names="RepeatSize1"),
        DataFrame(Number=".", Type="Float", Description="Entropy measure of first repeat based on percent composition", row.names="RepeatEntropy1"),
        DataFrame(Number="1", Type="Float", Description="Quality (same as value in QUAL column)", row.names="QUAL"),
        DataFrame(Number="0", Type="Flag", Description="Is this a segregating site (i.e. do parents have different genotypes", row.names="SEGREGATING"),
        DataFrame(Number="1", Type="Integer", Description="Number of Mendelian errors (parents have same genotype, progeny has differenet genotype) in progeny", row.names="MendelianErrors")
      )
    )
  )
  return(vcf)
}

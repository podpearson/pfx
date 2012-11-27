# uniqueSamples.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


uniqueSamples <- function(
  vcf,
  discordanceThreshold        = 10000,
  stemThreshold               = 2*discordanceThreshold,
  plotFilestem                = meta(exptData(vcf)[["header"]])["DataSetName", "Value"],
  shouldCreatePlots           = TRUE,
  shouldCalcMissingnessAndHet = "AD" %in% names(geno(vcf)),
  GTsToIntMapping             = c("0"=1, "1"=2, "."=0),
#  GTsToIntMapping             = c("7"=1, "G"=2, "."=0), # values for jiangVcf
  verbose                     = TRUE
) {
  require(ggplot2)
  require(reshape2)
  if(shouldCalcMissingnessAndHet) {
    ADas0123 <- genotypeCallsFromADas0123(vcf)
    missingnessPerSample <- apply(ADas0123, 2, function(x) length(which(x==0)))
    heterozygosityPerSample <- apply(ADas0123, 2, function(x) length(which(x==3)))
    if(shouldCreatePlots) {
      pdf(paste(plotFilestem, "missingnessPerSample.pdf", sep="."), height=5, width=8)
      print(
        qplot(
          x=names(missingnessPerSample),
          y=missingnessPerSample,
          xlab="Sample ID",
          ylab="Number of SNPs with missing genotype calls",
          geom="bar", stat="identity"
        )
        + theme_bw()
        + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
      )
      dev.off()
      pdf(paste(plotFilestem, "heterozygosityPerSample.pdf", sep="."), height=5, width=8)
      print(
        qplot(
          x=names(heterozygosityPerSample),
          y=heterozygosityPerSample,
          xlab="Sample ID",
          ylab="Number of SNPs with heterozygous genotype calls",
          geom="bar", stat="identity"
        )
        + theme_bw()
        + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
      )
      dev.off()
    }
  }
  discordance <- function(x, y) length(which(x!=y))
  vecDiscordance <- Vectorize(discordance)
  GTsInt <- genotypeCallsFromGTas012(vcf, GTsToIntMapping=GTsToIntMapping)  
  GTsIntWithNAs <- GTsInt
  GTsIntWithNAs[GTsIntWithNAs==0] <- NA
  GTsIntList <- split(GTsIntWithNAs, col(GTsIntWithNAs))
  names(GTsIntList) <- dimnames(GTsIntWithNAs)[[2]]
  GTsIntDiscordanceMatrix <- outer(GTsIntList, GTsIntList, vecDiscordance)
  diag(GTsIntDiscordanceMatrix) <- NA
  if(verbose) {
    stem(GTsIntDiscordanceMatrix)
    stem(GTsIntDiscordanceMatrix[GTsIntDiscordanceMatrix < stemThreshold])
  }
  if(shouldCreatePlots) {
    pdf(paste(plotFilestem, "pairwiseConcordanceHistogramAll.pdf", sep="."), height=4, width=6)
    print(
      qplot(
        as.vector(GTsIntDiscordanceMatrix),
        xlab="Number of discordant SNPs between pairwise sample comparisons",
        ylab="Frequency (number of sample pairs)"
      ) +
      geom_vline(xintercept = discordanceThreshold, colour="red") +
      theme_bw()
    )
    dev.off()
    if(any(GTsIntDiscordanceMatrix < discordanceThreshold, na.rm=TRUE)) {
      pdf(paste(plotFilestem, "pairwiseConcordanceHistogramDuplicates.pdf", sep="."), height=4, width=6)
      print(
        qplot(
          as.vector(GTsIntDiscordanceMatrix[GTsIntDiscordanceMatrix < discordanceThreshold]),
          xlab="Number of discordant SNPs between pairwise sample comparisons",
          ylab="Frequency (number of sample pairs)"
        ) +
        theme_bw()
      )
      dev.off()
    }
    pdf(paste(plotFilestem, "concordanceHeatmapAll.pdf", sep="."), height=6, width=8)
    print(
      ggplot(
        melt(GTsIntDiscordanceMatrix, value.name="Discordances"),
        aes(x=Var1, y=Var2, fill=Discordances)
      )
      + geom_tile()
      + scale_fill_gradient2(low="red", high="blue", midpoint=median(GTsIntDiscordanceMatrix, na.rm=TRUE))
#      + scale_fill_gradient(low="yellow", high="red")
      + theme_bw()
      + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
      + theme(axis.title.x = element_blank())
      + theme(axis.title.y = element_blank())
    )
    dev.off()
  }
  which(GTsIntDiscordanceMatrix < discordanceThreshold, arr.ind=TRUE)
  lowDiscordancePairs <- matrix(dimnames(GTsIntDiscordanceMatrix)[[1]][which(GTsIntDiscordanceMatrix < discordanceThreshold, arr.ind=TRUE)], ncol=2)
  duplicateSamplePairs <- matrix(lowDiscordancePairs[lowDiscordancePairs[, 1] != lowDiscordancePairs[, 2]], ncol=2)
  uniqueSamplePairs <- unique(apply(duplicateSamplePairs, 1, function(x) paste(sort(x), collapse="_and_")))
  if(verbose) {
    cat(paste(uniqueSamplePairs, collapse="\n"))
  }
  if(shouldCalcMissingnessAndHet) {
    duplicateSamplePairs <- cbind(duplicateSamplePairs, (missingnessPerSample[duplicateSamplePairs[, 1]] <= missingnessPerSample[duplicateSamplePairs[, 2]])+1)
    duplicateSamplePairs <- cbind(duplicateSamplePairs, (heterozygosityPerSample[duplicateSamplePairs[, 1]] <= heterozygosityPerSample[duplicateSamplePairs[, 2]])+1)
  } else {
    duplicateSamplePairs <- cbind(duplicateSamplePairs, (duplicateSamplePairs[, 1] <= duplicateSamplePairs[, 2])+1)
    duplicateSamplePairs <- cbind(duplicateSamplePairs, (duplicateSamplePairs[, 1] <= duplicateSamplePairs[, 2])+1)
  }
  samplesToRemove <- unique(duplicateSamplePairs[cbind(as.integer(seq(length=dim(duplicateSamplePairs)[1])), as.integer(duplicateSamplePairs[, 4]))])
  samplesToKeep <- setdiff(dimnames(GTsInt)[[2]], samplesToRemove)
  return(samplesToKeep)
}


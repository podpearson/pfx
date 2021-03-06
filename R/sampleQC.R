# sampleQC.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


sampleQC <- function(
  vcf,
  discordanceThreshold        = 10000,
  stemThreshold               = 2*discordanceThreshold,
  plotFilestem                = meta(exptData(vcf)[["header"]])["DataSetName", "Value"],
  shouldCreatePlots           = TRUE,
  shouldCalcMissingnessAndHet = "AD" %in% names(geno(vcf)),
  shouldCalcRecombinations    = TRUE,
  gffFilename                 = "/data/malariagen2/plasmodium/pf-crosses/data/genome/sanger/version3/September_2012/Pf3D7_v3.gff",
  gffGRL                      = readGffAsGRangesList(gffFilename, chromsomeNames=sprintf("Pf3D7_%02d_v3", 1:14)),
  GTsToIntMapping             = c("0"=1, "1"=2, "."=0, "./."=0),
#  GTsToIntMapping             = c("7"=1, "G"=2, "."=0), # values for jiangVcf
  parentalIDs                 = dimnames(vcf)[[2]][1:2],
  keepPASSvariantsOnly        = TRUE,
  shouldRenameSamples         = TRUE,
  sampleIDcolumn              = "ena_run_accession",
  sampleIDmappingsColumn      = sampleIDcolumn,
  thresholdLowDepth           = "3sd",
  thresholdHeterozgosity      = "3sd",
  thresholdMAF                = "3sd",
  thresholdMendelianErrors    = "3sd",
#  thresholdNoCalls            = "3sd",
  thresholdThirdOrFourthAllele= "3sd",
#  thresholdMissingness        = 50,
#  thresholdHeterozgosity      = 1000,
#  thresholdMendelianErrors    = 25,
  thresholdNoCalls            = 300,
#  thresholdThirdOrFourthAllele= 50,
  thresholdRecombinations     = "3sd",
  verbose                     = TRUE
) {
  require(ggplot2)
  require(reshape2)
  if(shouldRenameSamples) {
    sampleIDmappings <- createSampleIDmappings(
      sampleIDs                   = dimnames(vcf)[[2]],
      sampleIDcolumn              = sampleIDcolumn,
      sampleIDmappingsColumn      = sampleIDmappingsColumn
    )
    sampleNames <- names(sampleIDmappings)
    names(sampleNames) <- dimnames(vcf)[[2]]
#    vcf <- renameSamples(vcf)
  } else {
    sampleNames <- dimnames(vcf)[[2]]
    names(sampleNames) <- dimnames(vcf)[[2]]
  }
  if(keepPASSvariantsOnly) {
    vcf <- vcf[filt(vcf) %in% c("PASS", ".")]
  }
  if(shouldCalcMissingnessAndHet) {
    ADas0123 <- genotypeCallsFromADas0123(vcf)
    lowDepthVariantsPerSample <- apply(ADas0123, 2, function(x) length(which(x==0)))
    if(is.character(thresholdLowDepth)) {
      meanValue <- mean(lowDepthVariantsPerSample)
      sdValue <- sd(lowDepthVariantsPerSample)
      numberOfSDs <- as.integer(sub("sd", "", thresholdLowDepth))
      thresholdLowDepth <- meanValue + (numberOfSDs * sdValue)
    }
    heterozygosityPerSample <- apply(ADas0123, 2, function(x) length(which(x==3)))
    if(is.character(thresholdHeterozgosity)) {
      meanValue <- mean(heterozygosityPerSample)
      sdValue <- sd(heterozygosityPerSample)
      numberOfSDs <- as.integer(sub("sd", "", thresholdHeterozgosity))
      thresholdHeterozgosity <- meanValue + (numberOfSDs * sdValue)
    }
#    ADsArray <- array(
#      unlist(geno(vcf)[["AD"]]),
#      dim=c(2, dim(geno(vcf)[["AD"]])[1], dim(geno(vcf)[["AD"]])[2]),
#      dimnames=list(c("Ref", "Nonref"), dimnames(geno(vcf)[["AD"]])[[1]], dimnames(geno(vcf)[["AD"]])[[2]])
#    )
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
    meanMAFperSample <- colMeans(MAF, na.rm = TRUE)
    if(is.character(thresholdMAF)) {
      meanValue <- mean(meanMAFperSample)
      sdValue <- sd(meanMAFperSample)
      numberOfSDs <- as.integer(sub("sd", "", thresholdMAF))
      thresholdMAF <- meanValue + (numberOfSDs * sdValue)
    }

    if(shouldCreatePlots) {
      pdf(paste(plotFilestem, "lowDepthVariantsPerSample.pdf", sep="."), height=8, width=8)
      print(
        qplot(
          x=factor(sampleNames[names(lowDepthVariantsPerSample)], levels=sampleNames[names(lowDepthVariantsPerSample)]),
          y=lowDepthVariantsPerSample,
          xlab="Sample ID",
          ylab="Number of low depth SNPs",
          geom="bar", stat="identity"
        )
        + geom_hline(yintercept = thresholdLowDepth, colour="red")
        + theme_bw()
        + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
      )
      dev.off()
      pdf(paste(plotFilestem, "heterozygosityPerSample.pdf", sep="."), height=8, width=8)
      print(
        qplot(
          x=factor(sampleNames[names(heterozygosityPerSample)], levels=sampleNames[names(heterozygosityPerSample)]),
          y=heterozygosityPerSample,
          xlab="Sample ID",
          ylab="Number of heterozygous SNPs (>=2 ref and alt reads)",
          geom="bar", stat="identity"
        )
        + geom_hline(yintercept = thresholdHeterozgosity, colour="red")
        + theme_bw()
        + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
      )
      dev.off()
      pdf(paste(plotFilestem, "meanMAFperSample.pdf", sep="."), height=8, width=8)
      print(
        qplot(
          x=factor(sampleNames[names(meanMAFperSample)], levels=sampleNames[names(meanMAFperSample)]),
          y=meanMAFperSample,
          xlab="Sample ID",
          ylab="Mean number of minor allele reads per SNP",
          geom="bar", stat="identity"
        )
        + geom_hline(yintercept = thresholdMAF, colour="red")
        + theme_bw()
        + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
      )
      dev.off()
    }
  }
  discordance <- function(x, y) length(which(x!=y))
  vecDiscordance <- Vectorize(discordance)
  GTsInt <- genotypeCallsFromGTas012(vcf, GTsToIntMapping=GTsToIntMapping)
  if(shouldCreatePlots) {
    columnIndexesOfParents <- match(parentalIDs, dimnames(GTsInt)[[2]])
    columnIndexesOfProgeny <- setdiff(seq(along=dimnames(GTsInt)[[2]]), match(parentalIDs, dimnames(GTsInt)[[2]]))
    nocallGenotypesPerSample=apply(
      GTsInt,
      2,
      function(genotypesForVariant) {
        sum(genotypesForVariant==0, na.rm=TRUE)
      }
    )
    if(is.character(thresholdNoCalls)) {
      meanValue <- mean(nocallGenotypesPerSample)
      sdValue <- sd(nocallGenotypesPerSample)
      numberOfSDs <- as.integer(sub("sd", "", thresholdNoCalls))
      thresholdNoCalls <- meanValue + (numberOfSDs * sdValue)
    }
    pdf(paste(plotFilestem, "nocallGenotypesPerSample.pdf", sep="."), height=8, width=8)
    print(
      qplot(
        x=factor(sampleNames[names(nocallGenotypesPerSample)], levels=sampleNames[names(nocallGenotypesPerSample)]),
        y=nocallGenotypesPerSample,
        xlab="Sample ID",
        ylab="Number of no call genotypes",
        geom="bar", stat="identity"
      )
      + geom_hline(yintercept = thresholdNoCalls, colour="red")
      + theme_bw()
      + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
    )
    dev.off()
    thirdOrFourthAllelesPerSample=apply(
      GTsInt,
      2,
      function(genotypesForVariant) {
        sum(is.na(genotypesForVariant))
      }
    )
    if(is.character(thresholdThirdOrFourthAllele)) {
      meanValue <- mean(thirdOrFourthAllelesPerSample)
      sdValue <- sd(thirdOrFourthAllelesPerSample)
      numberOfSDs <- as.integer(sub("sd", "", thresholdThirdOrFourthAllele))
      thresholdThirdOrFourthAllele <- meanValue + (numberOfSDs * sdValue)
    }
    pdf(paste(plotFilestem, "thirdOrFourthAllelesPerSample.pdf", sep="."), height=8, width=8)
    print(
      qplot(
        x=factor(sampleNames[names(thirdOrFourthAllelesPerSample)], levels=sampleNames[names(thirdOrFourthAllelesPerSample)]),
        y=thirdOrFourthAllelesPerSample,
        xlab="Sample ID",
        ylab="Number of third or fourth allele genotypes",
        geom="bar", stat="identity"
      )
      + geom_hline(yintercept = thresholdThirdOrFourthAllele, colour="red")
      + theme_bw()
      + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
    )
    dev.off()
    MendelianErrorsPerSample=apply(
      GTsInt,
      2,
      function(genotypesForVariant) {
        sum(
          GTsInt[, parentalIDs[1]] == 1 &
            GTsInt[, parentalIDs[2]] == 1 &
            genotypesForVariant != 1,
          na.rm=TRUE
        ) +
        sum(
          GTsInt[, parentalIDs[1]] == 2 &
            GTsInt[, parentalIDs[2]] == 2 &
            genotypesForVariant != 2,
          na.rm=TRUE
        )
      }
    )
    if(is.character(thresholdMendelianErrors)) {
      meanValue <- mean(MendelianErrorsPerSample)
      sdValue <- sd(MendelianErrorsPerSample)
      numberOfSDs <- as.integer(sub("sd", "", thresholdMendelianErrors))
      thresholdMendelianErrors <- meanValue + (numberOfSDs * sdValue)
    }
    pdf(paste(plotFilestem, "mendelianErrorsPerSample.pdf", sep="."), height=8, width=8)
    print(
      qplot(
        x=factor(sampleNames[names(MendelianErrorsPerSample)], levels=sampleNames[names(MendelianErrorsPerSample)]),
        y=MendelianErrorsPerSample,
        xlab="Sample ID",
        ylab="Number of SNPs with Mendelian errors",
        geom="bar", stat="identity"
      )
      + geom_hline(yintercept = thresholdMendelianErrors, colour="red")
      + theme_bw()
      + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
    )
    dev.off()
  }
  GTsIntWithNAs <- GTsInt
  GTsIntWithNAs[GTsIntWithNAs==0] <- NA
  GTsIntList <- split(GTsIntWithNAs, col(GTsIntWithNAs))
  names(GTsIntList) <- dimnames(GTsIntWithNAs)[[2]]
#  names(GTsIntList) <- sampleNames[dimnames(GTsIntWithNAs)[[2]]]
  GTsIntDiscordanceMatrix <- outer(GTsIntList, GTsIntList, vecDiscordance)
  diag(GTsIntDiscordanceMatrix) <- NA
  if(verbose) {
    stem(GTsIntDiscordanceMatrix)
    if(length(which(GTsIntDiscordanceMatrix < stemThreshold)) > 0) {
      stem(GTsIntDiscordanceMatrix[GTsIntDiscordanceMatrix < stemThreshold])
    }
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
    discordanceDF <- melt(GTsIntDiscordanceMatrix, value.name="Discordances")
    discordanceDF[["sample1"]] <- sampleNames[as.character(discordanceDF[["Var1"]])]
    discordanceDF[["sample2"]] <- sampleNames[as.character(discordanceDF[["Var2"]])]
    pdf(paste(plotFilestem, "concordanceHeatmapAll.pdf", sep="."), height=10, width=12)
    print(
      ggplot(
        discordanceDF,
#        melt(GTsIntDiscordanceMatrix, value.name="Discordances"),
        aes(x=sample1, y=sample2, fill=Discordances)
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
  lowDiscordancePairsDF <- cbind(data.frame(lowDiscordancePairs), data.frame(discordantSNPs = GTsIntDiscordanceMatrix[lowDiscordancePairs]))
  duplicateSamplePairs <- matrix(lowDiscordancePairs[lowDiscordancePairs[, 1] != lowDiscordancePairs[, 2]], ncol=2)
  uniqueSamplePairs <- unique(apply(duplicateSamplePairs, 1, function(x) paste(sort(x), collapse="_and_")))
  if(verbose) {
    cat(paste(uniqueSamplePairs, collapse="\n"), "\n")
  }
  if(shouldCalcMissingnessAndHet) {
    lowDepthVariantsPerSample[parentalIDs] <- 0 # This is a fix to ensure parental strains are always selected (even if the quality of these is sometimes lower than that of a duplicate)
    heterozygosityPerSample[parentalIDs] <- 0   # This is a fix to ensure parental strains are always selected (even if the quality of these is sometimes lower than that of a duplicate)
    duplicateSamplePairs <- cbind(
      duplicateSamplePairs,
      (
        (lowDepthVariantsPerSample[duplicateSamplePairs[, 1]] < lowDepthVariantsPerSample[duplicateSamplePairs[, 2]]) |
        (lowDepthVariantsPerSample[duplicateSamplePairs[, 1]] == lowDepthVariantsPerSample[duplicateSamplePairs[, 2]] & duplicateSamplePairs[, 1] < duplicateSamplePairs[, 2])
      )+1
    )
    duplicateSamplePairs <- cbind(
      duplicateSamplePairs,
      (
        (heterozygosityPerSample[duplicateSamplePairs[, 1]] < heterozygosityPerSample[duplicateSamplePairs[, 2]]) |
        (heterozygosityPerSample[duplicateSamplePairs[, 1]] == heterozygosityPerSample[duplicateSamplePairs[, 2]] & duplicateSamplePairs[, 1] < duplicateSamplePairs[, 2])
      )+1
    )
  } else {
    duplicateSamplePairs <- cbind(duplicateSamplePairs, (duplicateSamplePairs[, 1] < duplicateSamplePairs[, 2])+1)
    duplicateSamplePairs <- cbind(duplicateSamplePairs, (duplicateSamplePairs[, 1] < duplicateSamplePairs[, 2])+1)
  }
  if(dim(duplicateSamplePairs)[1] > 0) {
    samplesToRemove <- unique(duplicateSamplePairs[cbind(as.integer(seq(length=dim(duplicateSamplePairs)[1])), as.integer(duplicateSamplePairs[, 3]))]) # remove duplicate samples that have non-lowest lowDepthVariants
  } else {
    samplesToRemove <- NULL
  }
#  samplesToRemove <- unique(duplicateSamplePairs[cbind(as.integer(seq(length=dim(duplicateSamplePairs)[1])), as.integer(duplicateSamplePairs[, 4]))])
  uniqueSamples <- setdiff(dimnames(GTsInt)[[2]], samplesToRemove)
  sampleDuplicates <- sapply(
    dimnames(duplicateSamplePairs)[[1]],
    function(sampleID) {
      paste(
        sort(
          unique(
            c(
              sampleID,
              duplicateSamplePairs[duplicateSamplePairs[,1]==sampleID, 2],
              duplicateSamplePairs[duplicateSamplePairs[,2]==sampleID, 1]
            )
          )
        ),
        collapse="_"
      )
    }
  )
      
      
  
  if(shouldCalcRecombinations) {
    mgRecombinations <- recombinationPoints(vcf, gffGRL, shouldCharacterise=FALSE, GTsToIntMapping=GTsToIntMapping)
    recombinationsPerSample <- rev(
      sapply(
        mgRecombinations[["sampleLevelResults"]],
        function(x) {
          sum(sapply(x, length))
        }
      )
    )
    if(is.character(thresholdRecombinations)) {
      meanRecombinations <- mean(recombinationsPerSample)
      sdRecombinations <- sd(recombinationsPerSample)
      numberOfSDs <- as.integer(sub("sd", "", thresholdRecombinations))
      thresholdRecombinations <- meanRecombinations + (numberOfSDs * sdRecombinations)
    }
    
    if(shouldCreatePlots) {
      pdf(paste(plotFilestem, "recombinationsPerSample.pdf", sep="."), height=8, width=8)
      print(
        qplot(
          x=factor(sampleNames[names(recombinationsPerSample)], levels=sampleNames[names(recombinationsPerSample)]),
          y=recombinationsPerSample,
          xlab="Sample ID",
          ylab="Number of apparent recombinations",
          geom="bar", stat="identity"
        )
        + geom_hline(yintercept = thresholdRecombinations, colour="red")
        + theme_bw()
        + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
      )
      dev.off()
    }
    qcFailedSamples <- names(
      which(
#      missingnessPerSample > thresholdMissingness |
#        heterozygosityPerSample > thresholdHeterozgosity |
#        MendelianErrorsPerSample > thresholdMendelianErrors |
        nocallGenotypesPerSample> thresholdNoCalls |
        recombinationsPerSample > thresholdRecombinations
      )
    )
  } else {
    qcFailedSamples <- NULL
    mgRecombinations <- NULL
  }
  
  return(
    list(
      uniqueSamples         = uniqueSamples,
      sampleDuplicates      = sampleDuplicates,
      qcFailedSamples       = qcFailedSamples,
      mgRecombinations      = mgRecombinations,
      discordanceMatrix     = GTsIntDiscordanceMatrix,
      lowDiscordancePairsDF = lowDiscordancePairsDF
    )
  )
}

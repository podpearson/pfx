# compareCalls.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


compareCalls <- function(
  subjectVcf,
  comparisonVcf,
  subjectName                 = "MalariaGEN",
  comparisonName              = "Jiang et al",
  distanceThresholds          = c(0, 22),
  discordanceThreshold        = 100,
  discordanceProportionThreshold = 0.15,
  malariagenDiscordanceThreshold = 200,
  malariagenDiscordanceProportionThreshold = 0.15,
  sampleAnnotationFilename    = "/data/malariagen2/plasmodium/pf-crosses/meta/qcmeta_annotated.tsv",
  comparisonDSthreshold       = NULL, # This is used specifically for Uberchip calls
  plotFilestem                = paste(meta(exptData(subjectVcf)[["header"]])["DataSetName", "Value"], "comparison", sep="."),
  shouldRenameSubjectSamples  = TRUE,
  IDparent1                   = dimnames(subjectVcf)[[2]][1],
  IDparent2                   = dimnames(subjectVcf)[[2]][2],
  IDcomparisonParent1         = dimnames(comparisonName)[[2]][1],
  IDcomparisonParent2         = dimnames(comparisonName)[[2]][2],
  shouldSubsetToBialleleic    = FALSE,
  shouldCompareRefsAndAlts    = FALSE,
  GTsToCompare                = c("parentBased", "asVcf"),
  GTsToIntMapping             = c("7"=1, "G"=2, "."=0),
#  GTsToIntMapping             = c("0"=1, "0/0"=1, "0|0"=1, "1"=2, "1/1"=2, "1|1"=2, "."=0, "./."=0, "./."=0, "2"=0, "3"=0, "0/1"=0, "1/0"=0, "0|1"=0, "1|0"=0) # "./." is needed as sometimes this is output by GATK's UG (presumably a bug). "2", "3", needed for the case of multi-allelic sites
  expectedMatches             = NULL
) {
  require(reshape2)
  if(shouldRenameSubjectSamples) {
    subjectVcf <- renameSamples(subjectVcf)
  }
#  sampleAnnotation <- readSampleAnnotation(sampleAnnotationFilename)
#  sampleAnnotation <- read.delim(sampleAnnotationFilename)
#  sampleAnnotation[["cross"]] <- ifelse(sampleAnnotation[["project_code"]]=="PFproj1", "Hb3xDd2", ifelse(sampleAnnotation[["project_code"]]=="PFproj2", "3d7xHb3", "7g8xGb4"))
#  sampleAnnotation[["qcStatus"]] <- ifelse(sampleAnnotation[["ajm_qc"]]=="fail", "fail", "pass")
#  sampleAnnotation <- sampleAnnotation[!grepl("not a cross progeny sample", sampleAnnotation[["ajm_qc_notes"]]), ]
#  sampleAnnotation <- sampleAnnotation[!duplicated(sampleAnnotation[["ox_code"]]), ]
#  row.names(sampleAnnotation) <- sampleAnnotation[["ox_code"]]
##  row.names(sampleAnnotation) <- gsub("-", "\\.", sampleAnnotation[["ox_code"]])
  
  if(shouldSubsetToBialleleic) {
    subjectVcf <- subjectVcf[elementLengths(alt(subjectVcf)) == 1]
    comparisonVcf <- comparisonVcf[elementLengths(alt(comparisonVcf)) == 1]
  }
  if(!is.null(comparisonDSthreshold)) {
    geno(comparisonVcf)[["GT"]][geno(comparisonVcf)[["DS"]] > comparisonDSthreshold] <- "."
  }
  comparisonNearestSubjectGR <- rowData(subjectVcf)[nearest(rowData(comparisonVcf), rowData(subjectVcf))]
  table(start(ranges(comparisonNearestSubjectGR)) - start(ranges(rowData(comparisonVcf))))
  table(start(ranges(comparisonNearestSubjectGR))-start(ranges(rowData(comparisonVcf))) >=0 & start(ranges(comparisonNearestSubjectGR))-start(ranges(rowData(comparisonVcf)))<=22)
  subjectVsComparisonPositionDifferences <- start(ranges(comparisonNearestSubjectGR))-start(ranges(rowData(comparisonVcf)))
  require(ggplot2)
  pdf(paste(plotFilestem, "positionDifferences.pdf", sep="."), height=4, width=6)
  print(
    qplot(subjectVsComparisonPositionDifferences, binwidth=1, xlim=c(-50, 50), xlab=paste("Distance (bp) from", comparisonName, "SNP to nearest", subjectName, "SNP"), ylab="Frequency (number of SNPs)")
    + theme_bw()
  )
  dev.off()
  
  comparisonGTsas012 <- genotypeCallsFromGTas012(comparisonVcf, GTsToIntMapping = GTsToIntMapping)
  segregatingComparisonVariants <- which(
    comparisonGTsas012[, IDcomparisonParent1] != 0 &
    comparisonGTsas012[, IDcomparisonParent2] != 0 &
    comparisonGTsas012[, IDcomparisonParent1] != comparisonGTsas012[, IDcomparisonParent2]
  )
  
  comparisonRowsWithMatchesInSubject <- which(
    subjectVsComparisonPositionDifferences >= distanceThresholds[1] &
    subjectVsComparisonPositionDifferences <= distanceThresholds[2]
  )
#  positionsWithinThresholdsText <- paste(length(comparisonRowsWithMatchesInSubject), "of", length(subjectVsComparisonPositionDifferences), "SNPs within thresholds")
  positionsWithinThresholdsText <- paste(length(comparisonRowsWithMatchesInSubject), "of", length(subjectVsComparisonPositionDifferences), comparisonName, "SNPs are in", subjectName)
  pdf(paste(plotFilestem, "positionDifferencesWithThresholds.pdf", sep="."), height=4, width=6)
  print(
    qplot(subjectVsComparisonPositionDifferences, binwidth=1, xlim=c(-50, 50), xlab=paste("Distance (bp) from", comparisonName, "SNP to nearest", subjectName, "SNP"), ylab="Frequency (number of SNPs)")
    + geom_vline(xintercept = distanceThresholds[1], colour="red")
    + geom_vline(xintercept = distanceThresholds[2]+1, colour="red")
    + theme_bw()
    + labs(title = positionsWithinThresholdsText)
#    + geom_text(data = data.frame(), aes(-50, 100, label = positionsWithinThresholdsText))
  )
  dev.off()
  
  sensitivityResults <- numeric(0)
  sensitivityResults["numberOfSegregatingSNPsInComparison"] <- length(segregatingComparisonVariants)
  sensitivityResults["numberOfSegregatingSNPsWithPositionMatches"] <- length(intersect(segregatingComparisonVariants, comparisonRowsWithMatchesInSubject))

  subjectRowsWithMatchesInComparison <- nearest(rowData(comparisonVcf), rowData(subjectVcf))[comparisonRowsWithMatchesInSubject]
#  comparisonRowsWithMatchesInSubject <- which(subjectVsComparisonPositionDifferences >=distanceThresholds[1] & subjectVsComparisonPositionDifferences<=distanceThresholds[2])
#  subjectRowsWithMatchesInComparison <- nearest(rowData(comparisonVcf), rowData(subjectVcf))[comparisonRowsWithMatchesInSubject]
  table(as.character(ref(subjectVcf[subjectRowsWithMatchesInComparison])), useNA="ifany")
  table(as.character(ref(comparisonVcf[comparisonRowsWithMatchesInSubject])), useNA="ifany")
  table(
    as.character(ref(subjectVcf[subjectRowsWithMatchesInComparison])),
    as.character(ref(comparisonVcf[comparisonRowsWithMatchesInSubject])),
    useNA="ifany"
  )
  if(shouldCompareRefsAndAlts) {
    comparisonRowsWithMatchesInSubject2 <- comparisonRowsWithMatchesInSubject[which(as.character(ref(comparisonVcf[comparisonRowsWithMatchesInSubject])) == as.character(ref(subjectVcf[subjectRowsWithMatchesInComparison])))]
    subjectRowsWithMatchesInComparison2 <- subjectRowsWithMatchesInComparison[which(as.character(ref(comparisonVcf[comparisonRowsWithMatchesInSubject])) == as.character(ref(subjectVcf[subjectRowsWithMatchesInComparison])))]
    table(
      as.character(ref(subjectVcf[subjectRowsWithMatchesInComparison2])),
      as.character(ref(comparisonVcf[comparisonRowsWithMatchesInSubject2])),
      useNA="ifany"
    )
    table(
      as.character(unlist(alt(subjectVcf[subjectRowsWithMatchesInComparison2]))),
      as.character(unlist(alt(comparisonVcf[comparisonRowsWithMatchesInSubject2]))),
      useNA="ifany"
    )
    comparisonRowsWithMatchesInSubject3 <- comparisonRowsWithMatchesInSubject2[which(as.character(unlist(alt(comparisonVcf[comparisonRowsWithMatchesInSubject2]))) == as.character(unlist(alt(subjectVcf[subjectRowsWithMatchesInComparison2]))))]
    subjectRowsWithMatchesInComparison3 <- subjectRowsWithMatchesInComparison2[which(as.character(unlist(alt(comparisonVcf[comparisonRowsWithMatchesInSubject2]))) == as.character(unlist(alt(subjectVcf[subjectRowsWithMatchesInComparison2]))))]
    table(
      as.character(ref(subjectVcf[subjectRowsWithMatchesInComparison3])),
      as.character(ref(comparisonVcf[comparisonRowsWithMatchesInSubject3])),
      useNA="ifany"
    )
    table(
      as.character(unlist(alt(subjectVcf[subjectRowsWithMatchesInComparison3]))),
      as.character(unlist(alt(comparisonVcf[comparisonRowsWithMatchesInSubject3]))),
      useNA="ifany"
    )
    allComparisonRowsWithMathcesInSubject <- comparisonRowsWithMatchesInSubject
    comparisonRowsWithMatchesInSubject <- comparisonRowsWithMatchesInSubject3
    subjectRowsWithMatchesInComparison <- subjectRowsWithMatchesInComparison3
    sensitivityResults["numberOfSegregatingSNPsWithPositionAndAlleleMatches"] <- length(intersect(segregatingComparisonVariants, comparisonRowsWithMatchesInSubject3))
  }
  
  if(GTsToCompare[1] == "parentBased") {
    subjectGTsCFparents <- convertGTsIntToParentBasedGTs(genotypeCallsFromGTas012(subjectVcf), IDparent1=IDparent1, IDparent2=IDparent2)
    comparisonGTsCFparents <- convertGTsIntToParentBasedGTs(genotypeCallsFromGTas012(comparisonVcf, GTsToIntMapping = GTsToIntMapping))
    comparisonMatchesGTs <- comparisonGTsCFparents[comparisonRowsWithMatchesInSubject,]
    subjectMatchesGTs <- subjectGTsCFparents[subjectRowsWithMatchesInComparison, ]
  #  jiangGTsGWInt <- -(matrix(as.integer(factor(jiangMatchesGTsGW)), nrow=nrow(jiangMatchesGTsGW))-2)
  #  dimnames(jiangGTsGWInt) <- dimnames(jiangMatchesGTsGW)
  } else if(GTsToCompare[1] == "asVcf") {
    subjectGTs <- genotypeCallsFromGTas012(subjectVcf)
    comparisonGTs <- genotypeCallsFromGTas012(comparisonVcf, GTsToIntMapping = GTsToIntMapping)
    comparisonMatchesGTs <- comparisonGTs[comparisonRowsWithMatchesInSubject,]
    subjectMatchesGTs <- subjectGTs[subjectRowsWithMatchesInComparison, ]
    comparisonMatchesGTs[comparisonMatchesGTs==0] <- NA
    subjectMatchesGTs[subjectMatchesGTs==0] <- NA
  } else {
    stop("GTsToCompare must be either \"parentBased\" or \"asVcf\"")
  }
  
  return(
    list(
      comparisonVsSubjectDiscordanceMatrix = callComparisonPlots(
        subjectMatchesGTs,
        comparisonMatchesGTs,
        comparison="discordances",
        threshold=discordanceThreshold,
        plotFilestem=plotFilestem,
        expectedMatches=expectedMatches[["comparisonVsSubject"]],
        subjectName=subjectName,
        comparisonName=comparisonName
      ),
      comparisonVsSubjectDiscordanceProportionMatrix = callComparisonPlots(
        subjectMatchesGTs,
        comparisonMatchesGTs,
        comparison="discordanceProportions",
        threshold=discordanceProportionThreshold,
        plotFilestem=plotFilestem,
        expectedMatches=expectedMatches[["comparisonVsSubject"]],
        subjectName=subjectName,
        comparisonName=comparisonName
      ),
      subjectVsSubjectDiscordanceMatrix = callComparisonPlots(
        subjectMatchesGTs,
        subjectMatchesGTs,
        comparison="discordances",
        threshold=discordanceThreshold,
        plotFilestem=paste(plotFilestem, "subject", sep="."),
        expectedMatches=expectedMatches[["subjectVsSubject"]],
        subjectName=subjectName,
        comparisonName=subjectName,
        discordanceHeatmapWidth=12
      ),
      subjectVsSubjectDiscordanceProportionMatrix = callComparisonPlots(
        subjectMatchesGTs,
        subjectMatchesGTs,
        comparison="discordanceProportions",
        threshold=discordanceProportionThreshold,
        plotFilestem=paste(plotFilestem, "subject", sep="."),
        expectedMatches=expectedMatches[["subjectVsSubject"]],
        subjectName=subjectName,
        comparisonName=subjectName,
        discordanceHeatmapWidth=12
      ),
      comparisonVsComparisonDiscordanceMatrix = callComparisonPlots(
        comparisonMatchesGTs,
        comparisonMatchesGTs,
        comparison="discordances",
        threshold=discordanceThreshold,
        plotFilestem=paste(plotFilestem, "comparison", sep="."),
        expectedMatches=expectedMatches[["comparisonVsComparison"]],
        subjectName=comparisonName,
        comparisonName=comparisonName,
        discordanceHeatmapHeight=8
      ),
      comparisonVsComparisonDiscordanceProportionMatrix = callComparisonPlots(
        comparisonMatchesGTs,
        comparisonMatchesGTs,
        comparison="discordanceProportions",
        threshold=discordanceProportionThreshold,
        plotFilestem=paste(plotFilestem, "comparison", sep="."),
        expectedMatches=expectedMatches[["comparisonVsComparison"]],
        subjectName=comparisonName,
        comparisonName=comparisonName,
        discordanceHeatmapHeight=8
      ),
      sensitivityResults = sensitivityResults
    )
  )
#
#  callComparisonPlots(subjectMatchesGTs, comparisonMatchesGTs, comparison="discordances", threshold=discordanceThreshold, plotFilestem=plotFilestem, expectedMatches=expectedMatches[["comparisonVsSubject"]], subjectName=subjectName, comparisonName=comparisonName)
#  callComparisonPlots(subjectMatchesGTs, comparisonMatchesGTs, comparison="discordanceProportions", threshold=discordanceProportionThreshold, plotFilestem=plotFilestem, expectedMatches=expectedMatches[["comparisonVsSubject"]], subjectName=subjectName, comparisonName=comparisonName)
#  callComparisonPlots(subjectMatchesGTs, subjectMatchesGTs, comparison="discordances", threshold=discordanceThreshold, plotFilestem=paste(plotFilestem, "subject", sep="."), expectedMatches=expectedMatches[["subjectVsSubject"]], subjectName=subjectName, comparisonName=subjectName, discordanceHeatmapWidth=12)
#  callComparisonPlots(subjectMatchesGTs, subjectMatchesGTs, comparison="discordanceProportions", threshold=discordanceProportionThreshold, plotFilestem=paste(plotFilestem, "subject", sep="."), expectedMatches=expectedMatches[["subjectVsSubject"]], subjectName=subjectName, comparisonName=subjectName, discordanceHeatmapWidth=12)
#  callComparisonPlots(comparisonMatchesGTs, comparisonMatchesGTs, comparison="discordances", threshold=discordanceThreshold, plotFilestem=paste(plotFilestem, "comparison", sep="."), expectedMatches=expectedMatches[["comparisonVsComparison"]], subjectName=comparisonName, comparisonName=comparisonName, discordanceHeatmapHeight=8)
#  callComparisonPlots(comparisonMatchesGTs, comparisonMatchesGTs, comparison="discordanceProportions", threshold=discordanceProportionThreshold, plotFilestem=paste(plotFilestem, "comparison", sep="."), expectedMatches=expectedMatches[["comparisonVsComparison"]], subjectName=comparisonName, comparisonName=comparisonName, discordanceHeatmapHeight=8)
#  
#  discordance <- function(x, y) length(which(x!=y))
#  vecDiscordance <- Vectorize(discordance)
#  comparisonVsSubjectDiscordanceMatrix <- outer(data.frame(comparisonMatchesGTs), data.frame(subjectMatchesGTs), vecDiscordance)
#  discordanceProportion <- function(x, y) length(which(x!=y)) / length(which(!is.na(x) & !is.na(y)))
#  vecDiscordanceProportion <- Vectorize(discordanceProportion)
#  comparisonVsSubjectDiscordanceProportionMatrix <- outer(data.frame(comparisonMatchesGTs), data.frame(subjectMatchesGTs), vecDiscordanceProportion)
#
#  discordancesText <- paste("Median pairwise discordance =", median(comparisonVsSubjectDiscordanceMatrix), "\nMedian discordance of putative matches = ", median(comparisonVsSubjectDiscordanceMatrix[comparisonVsSubjectDiscordanceMatrix<discordanceThreshold]))
#  pdf(paste(plotFilestem, "discordances.pdf", sep="."), height=4, width=6)
#  print(
#    qplot(as.vector(comparisonVsSubjectDiscordanceMatrix), binwidth=50, xlab=paste("Number of discordant SNPs between", comparisonName, "and", subjectName, "samples"), ylab="Frequency (number of sample pairs)")
#    + geom_vline(xintercept = discordanceThreshold, colour="red")
#    + theme_bw()
#    + labs(title=discordancesText)
#  )
#  dev.off()
#  pdf(paste(plotFilestem, "discordancesDuplicates.pdf", sep="."), height=4, width=6)
#  print(
#    qplot(as.vector(comparisonVsSubjectDiscordanceMatrix[comparisonVsSubjectDiscordanceMatrix<discordanceThreshold]), binwidth=1, xlab=paste("Number of discordant SNPs between", comparisonName, "and", subjectName, "samples"), ylab="Frequency (number of sample pairs)")
#    + theme_bw()
#  )
#  dev.off()
#  if(!is.null(expectedMatches) && "comparisonVsSubject" %in% names(expectedMatches)) {
#    comparisonVsSubjectDiscordanceMatrixExpectMatches <- comparisonVsSubjectDiscordanceMatrix[expectedMatches[["comparisonVsSubject"]]]
#    discordancesText <- paste("Median discordances in expected matches =", median(comparisonVsSubjectDiscordanceMatrixExpectMatches))
#    pdf(paste(plotFilestem, "discordancesExpectedMatches.pdf", sep="."), height=4, width=6)
#    print(
#      qplot(as.vector(comparisonVsSubjectDiscordanceMatrixExpectMatches), xlab=paste("Number of discordant SNPs between", comparisonName, "and", subjectName, "expected match samples"), ylab="Frequency (number of sample pairs)")
#      + geom_vline(xintercept = discordanceThreshold, colour="red")
#      + theme_bw()
#      + labs(title=discordancesText)
#    )
#    dev.off()
#    discordancesText <- paste("Median discordances in true expected matches = ", median(comparisonVsSubjectDiscordanceMatrix[comparisonVsSubjectDiscordanceMatrix<discordanceThreshold]))
#    pdf(paste(plotFilestem, "discordancesDuplicatesExpectedMatches.pdf", sep="."), height=4, width=6)
#    print(
#      qplot(as.vector(comparisonVsSubjectDiscordanceMatrixExpectMatches[comparisonVsSubjectDiscordanceMatrixExpectMatches<discordanceThreshold]), xlab=paste("Number of discordant SNPs between", comparisonName, "and", subjectName, "expected match samples"), ylab="Frequency (number of sample pairs)")
#      + theme_bw()
#      + labs(title=discordancesText)
#    )
#    dev.off()
#  }
#  
#  discordanceDF <- melt(comparisonVsSubjectDiscordanceMatrix, value.name="Discordances")
#  discordanceDF[["putativeDuplicateSample"]] <- discordanceDF[["Discordances"]] <= discordanceThreshold
#  if(!is.null(expectedMatches) && "comparisonVsSubject" %in% names(expectedMatches)) {
#    comparisonVsSubjectExpectedMatches <- comparisonVsSubjectDiscordanceMatrix
#    comparisonVsSubjectExpectedMatches[,] <- 0
#    comparisonVsSubjectExpectedMatches[expectedMatches[["comparisonVsSubject"]]] <- 1
#    expectedMatchesDF <- melt(comparisonVsSubjectExpectedMatches, value.name="ExpectedMatch")
#    discordanceDF[["expectedDuplicateSample"]] <- as.logical(expectedMatchesDF[["ExpectedMatch"]])
#    discordanceDF[discordanceDF[["expectedDuplicateSample"]] & discordanceDF[["putativeDuplicateSample"]], "SamplePairStatus"] <- "Expected match"
#    discordanceDF[!discordanceDF[["expectedDuplicateSample"]] & discordanceDF[["putativeDuplicateSample"]], "SamplePairStatus"] <- "Unexpected match"
#    discordanceDF[discordanceDF[["expectedDuplicateSample"]] & !discordanceDF[["putativeDuplicateSample"]], "SamplePairStatus"] <- "Unexpected non-match"
#    discordanceDF[!discordanceDF[["expectedDuplicateSample"]] & !discordanceDF[["putativeDuplicateSample"]], "SamplePairStatus"] <- "Expected non-match"
#  }
#  
#  pdf(paste(plotFilestem, "discordanceHeatmap.pdf", sep="."), height=10, width=10)
#  print(
#    ggplot(
#      discordanceDF,
##      melt(comparisonVsSubjectDiscordanceMatrix, value.name="Discordances"),
#      aes(x=Var1, y=Var2, fill=Discordances)
#    )
#    + geom_tile()
#    + scale_fill_gradient2(low="red", high="blue", midpoint=median(comparisonVsSubjectDiscordanceMatrix, na.rm=TRUE))
#    + theme_bw()
#    + xlab(paste(comparisonName, "sample ID"))
#    + ylab(paste(subjectName, "sample ID"))
#    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
##    + theme(axis.title.x = text(labels="Jiang et al. sample ID"))
##    + theme(axis.title.y = element_blank())
#  )
#  dev.off()
#
#  pdf(paste(plotFilestem, "putativeIdenticalSamples.pdf", sep="."), height=10, width=10)
#  print(
#    ggplot(
#      discordanceDF,
#      aes(x=Var1, y=Var2, fill=putativeDuplicateSample)
#    )
#    + geom_tile()
#    + scale_fill_grey()
#    + theme_bw()
#    + xlab(paste(comparisonName, "sample ID"))
#    + ylab(paste(subjectName, "sample ID"))
#    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
#  )
#  dev.off()
#  
#  if(!is.null(expectedMatches) && "comparisonVsSubject" %in% names(expectedMatches)) {
#    pdf(paste(plotFilestem, "samplePairStatus.pdf", sep="."), height=10, width=10)
#    print(
#      ggplot(
#        discordanceDF,
#        aes(x=Var1, y=Var2, fill=SamplePairStatus)
#      )
#      + geom_tile()
#      + scale_fill_manual(values = c("green", "grey90", "orange", "red"))
#      + theme_bw()
#      + xlab(paste(comparisonName, "sample ID"))
#      + ylab(paste(subjectName, "sample ID"))
#      + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
#    )
#    dev.off()
#  }
#  
#  discordancesText <- paste("Median pairwise discordance proportion =", round(median(comparisonVsSubjectDiscordanceProportionMatrix), 3), "\nMedian discordance proportion of putative matches = ", round(median(comparisonVsSubjectDiscordanceProportionMatrix[comparisonVsSubjectDiscordanceProportionMatrix<discordanceProportionThreshold]), 3))
#  pdf(paste(plotFilestem, "discordanceProportions.pdf", sep="."), height=4, width=6)
#  print(
#    qplot(as.vector(comparisonVsSubjectDiscordanceProportionMatrix), xlab=paste("Proportion of discordant SNPs between", comparisonName, "and", subjectName, "samples"), ylab="Frequency (number of sample pairs)")
#    + geom_vline(xintercept = discordanceProportionThreshold, colour="red")
#    + theme_bw()
#    + labs(title=discordancesText)
#  )
#  dev.off()
#  pdf(paste(plotFilestem, "discordanceProportionsDuplicates.pdf", sep="."), height=4, width=6)
#  print(
#    qplot(as.vector(comparisonVsSubjectDiscordanceProportionMatrix[comparisonVsSubjectDiscordanceProportionMatrix<discordanceProportionThreshold]), xlab=paste("Proportion of discordant SNPs between", comparisonName, "and", subjectName, "samples"), ylab="Frequency (number of sample pairs)")
#    + theme_bw()
#  )
#  dev.off()
#    discordancesText <- paste("Median discordances in expected matches =", median(comparisonVsSubjectDiscordanceMatrixExpectMatches))
#
#  discordanceProportionDF <- melt(comparisonVsSubjectDiscordanceProportionMatrix, value.name="DiscordanceProportions")
#  discordanceProportionDF[["putativeDuplicateSample"]] <- discordanceProportionDF[["DiscordanceProportions"]] <= discordanceProportionThreshold
#  if(!is.null(expectedMatches) && "comparisonVsSubject" %in% names(expectedMatches)) {
#    comparisonVsSubjectExpectedMatches <- comparisonVsSubjectDiscordanceMatrix
#    comparisonVsSubjectExpectedMatches[,] <- 0
#    comparisonVsSubjectExpectedMatches[expectedMatches[["comparisonVsSubject"]]] <- 1
#    expectedMatchesDF <- melt(comparisonVsSubjectExpectedMatches, value.name="ExpectedMatch")
#    discordanceProportionDF[["expectedDuplicateSample"]] <- as.logical(expectedMatchesDF[["ExpectedMatch"]])
#    discordanceProportionDF[discordanceProportionDF[["expectedDuplicateSample"]] & discordanceProportionDF[["putativeDuplicateSample"]], "SamplePairStatus"] <- "Expected match"
#    discordanceProportionDF[!discordanceProportionDF[["expectedDuplicateSample"]] & discordanceProportionDF[["putativeDuplicateSample"]], "SamplePairStatus"] <- "Unexpected match"
#    discordanceProportionDF[discordanceProportionDF[["expectedDuplicateSample"]] & !discordanceProportionDF[["putativeDuplicateSample"]], "SamplePairStatus"] <- "Unexpected non-match"
#    discordanceProportionDF[!discordanceProportionDF[["expectedDuplicateSample"]] & !discordanceProportionDF[["putativeDuplicateSample"]], "SamplePairStatus"] <- "Expected non-match"
#  }
#  
#  pdf(paste(plotFilestem, "discordanceProportionHeatmap.pdf", sep="."), height=10, width=10)
#  print(
#    ggplot(
#      discordanceProportionDF,
##      melt(comparisonVsSubjectDiscordanceMatrix, value.name="Discordances"),
#      aes(x=Var1, y=Var2, fill=DiscordanceProportions)
#    )
#    + geom_tile()
#    + scale_fill_gradient2(low="red", high="blue", midpoint=median(comparisonVsSubjectDiscordanceProportionMatrix, na.rm=TRUE))
#    + theme_bw()
#    + xlab(paste(comparisonName, "sample ID"))
#    + ylab(paste(subjectName, "sample ID"))
#    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
##    + theme(axis.title.x = text(labels="Jiang et al. sample ID"))
##    + theme(axis.title.y = element_blank())
#  )
#  dev.off()
#
#  pdf(paste(plotFilestem, "putativeIdenticalSamplesFromProportions.pdf", sep="."), height=10, width=10)
#  print(
#    ggplot(
#      discordanceProportionDF,
#      aes(x=Var1, y=Var2, fill=putativeDuplicateSample)
#    )
#    + geom_tile()
#    + scale_fill_grey()
#    + theme_bw()
#    + xlab(paste(comparisonName, "sample ID"))
#    + ylab(paste(subjectName, "sample ID"))
#    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
#  )
#  dev.off()
#  
#  if(!is.null(expectedMatches) && "comparisonVsSubject" %in% names(expectedMatches)) {
#    pdf(paste(plotFilestem, "samplePairStatusUsingProportions.pdf", sep="."), height=10, width=10)
#    print(
#      ggplot(
#        discordanceProportionDF,
#        aes(x=Var1, y=Var2, fill=SamplePairStatus)
#      )
#      + geom_tile()
#      + scale_fill_manual(values = c("green", "grey90", "orange", "red"))
#      + theme_bw()
#      + xlab(paste(comparisonName, "sample ID"))
#      + ylab(paste(subjectName, "sample ID"))
#      + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
#    )
#    dev.off()
#  }
#  
#  
#  
#  subjectVsSubjectDiscordanceMatrix <- outer(data.frame(subjectMatchesGTs), data.frame(subjectMatchesGTs), vecDiscordance)
#  subjectVsSubjectDiscordanceProportionMatrix <- outer(data.frame(subjectMatchesGTs), data.frame(subjectMatchesGTs), vecDiscordanceProportion)
#
#  pdf(paste(plotFilestem, "malariagenDiscordances.pdf", sep="."), height=4, width=6)
#  print(
#    qplot(as.vector(subjectVsSubjectDiscordanceMatrix), binwidth=50, xlab=paste("Number of discordant SNPs between pairs of", subjectName, "samples"), ylab="Frequency (number of sample pairs)")
#    + geom_vline(xintercept = malariagenDiscordanceThreshold, colour="red")
#    + theme_bw()
#  )
#  dev.off()
#  pdf(paste(plotFilestem, "malariagenDiscordancesDuplicates.pdf", sep="."), height=4, width=6)
#  print(
#    qplot(as.vector(subjectVsSubjectDiscordanceMatrix[subjectVsSubjectDiscordanceMatrix<malariagenDiscordanceThreshold]), binwidth=1, xlab=paste("Number of discordant SNPs between pairs of", subjectName, "samples"), ylab="Frequency (number of sample pairs)")
#    + theme_bw()
#  )
#  dev.off()
#
#  subjectDiscordanceDF <- melt(subjectVsSubjectDiscordanceMatrix, value.name="Discordances")
#  subjectDiscordanceDF[["putativeDuplicateSample"]] <- subjectDiscordanceDF[["Discordances"]] <= malariagenDiscordanceThreshold
#  
#  pdf(paste(plotFilestem, "malariagenDiscordanceHeatmap.pdf", sep="."), height=10, width=12)
#  print(
#    ggplot(
#      subjectDiscordanceDF,
##      melt(comparisonVsSubjectDiscordanceMatrix, value.name="Discordances"),
#      aes(x=Var1, y=Var2, fill=Discordances)
#    )
#    + geom_tile()
#    + scale_fill_gradient2(low="red", high="blue", midpoint=median(subjectVsSubjectDiscordanceMatrix, na.rm=TRUE))
#    + theme_bw()
#    + xlab(paste(subjectName, "sample ID"))
#    + ylab(paste(subjectName, "sample ID"))
#    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
#  )
#  dev.off()
#
#  pdf(paste(plotFilestem, "malariagenPutativeIdenticalSamples.pdf", sep="."), height=10, width=12)
#  print(
#    ggplot(
#      subjectDiscordanceDF,
#      aes(x=Var1, y=Var2, fill=putativeDuplicateSample)
#    )
#    + geom_tile()
#    + scale_fill_grey()
#    + theme_bw()
#    + xlab(paste(subjectName, "sample ID"))
#    + ylab(paste(subjectName, "sample ID"))
#    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
#  )
#  dev.off()
#  
#  pdf(paste(plotFilestem, "malariagenDiscordanceProportions.pdf", sep="."), height=4, width=6)
#  print(
#    qplot(as.vector(subjectVsSubjectDiscordanceProportionMatrix), xlab=paste("Proportion of discordant SNPs between pairs of", subjectName, "samples"), ylab="Frequency (number of sample pairs)")
#    + geom_vline(xintercept = malariagenDiscordanceProportionThreshold, colour="red")
#    + theme_bw()
#  )
#  dev.off()
#  pdf(paste(plotFilestem, "malariagenDiscordanceProportionsDuplicates.pdf", sep="."), height=4, width=6)
#  print(
#    qplot(as.vector(subjectVsSubjectDiscordanceProportionMatrix[subjectVsSubjectDiscordanceProportionMatrix<malariagenDiscordanceProportionThreshold]), binwidth=1, xlab=paste("Proportion of discordant SNPs between pairs of", subjectName, "samples"), ylab="Frequency (number of sample pairs)")
#    + theme_bw()
#  )
#  dev.off()
#
#  subjectDiscordanceProportionDF <- melt(subjectVsSubjectDiscordanceProportionMatrix, value.name="DiscordanceProportions")
#  subjectDiscordanceProportionDF[["putativeDuplicateSample"]] <- subjectDiscordanceProportionDF[["DiscordanceProportions"]] <= malariagenDiscordanceProportionThreshold
#  
#  pdf(paste(plotFilestem, "malariagenDiscordanceProportionHeatmap.pdf", sep="."), height=10, width=12)
#  print(
#    ggplot(
#      subjectDiscordanceProportionDF,
##      melt(comparisonVsSubjectDiscordanceMatrix, value.name="Discordances"),
#      aes(x=Var1, y=Var2, fill=DiscordanceProportions)
#    )
#    + geom_tile()
#    + scale_fill_gradient2(low="red", high="blue", midpoint=median(subjectVsSubjectDiscordanceProportionMatrix, na.rm=TRUE))
#    + theme_bw()
#    + xlab(paste(subjectName, "sample ID"))
#    + ylab(paste(subjectName, "sample ID"))
#    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
#  )
#  dev.off()
#
#  pdf(paste(plotFilestem, "malariagenPutativeIdenticalSamplesFromProportions.pdf", sep="."), height=10, width=12)
#  print(
#    ggplot(
#      subjectDiscordanceProportionDF,
#      aes(x=Var1, y=Var2, fill=putativeDuplicateSample)
#    )
#    + geom_tile()
#    + scale_fill_grey()
#    + theme_bw()
#    + xlab(paste(subjectName, "sample ID"))
#    + ylab(paste(subjectName, "sample ID"))
#    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
#  )
#  dev.off()
#  
#  
#  
#  comparisonVsComparisonDiscordanceMatrix <- outer(data.frame(comparisonMatchesGTs), data.frame(comparisonMatchesGTs), vecDiscordance)
#  comparisonVsComparisonDiscordanceProportionMatrix <- outer(data.frame(comparisonMatchesGTs), data.frame(comparisonMatchesGTs), vecDiscordanceProportion)
#
#  pdf(paste(plotFilestem, "comparisonDiscordances.pdf", sep="."), height=4, width=6)
#  print(
#    qplot(as.vector(comparisonVsComparisonDiscordanceMatrix), binwidth=50, xlab=paste("Number of discordant SNPs between pairs of", comparisonName, "samples"), ylab="Frequency (number of sample pairs)")
#    + geom_vline(xintercept = discordanceThreshold, colour="red")
#    + theme_bw()
#  )
#  dev.off()
#  pdf(paste(plotFilestem, "comparisonDiscordancesDuplicates.pdf", sep="."), height=4, width=6)
#  print(
#    qplot(as.vector(comparisonVsComparisonDiscordanceMatrix[comparisonVsComparisonDiscordanceMatrix<discordanceThreshold]), binwidth=1, xlab=paste("Number of discordant SNPs between pairs of", comparisonName, "samples"), ylab="Frequency (number of sample pairs)")
#    + theme_bw()
#  )
#  dev.off()
#
#  comparisonDiscordanceDF <- melt(comparisonVsComparisonDiscordanceMatrix, value.name="Discordances")
#  comparisonDiscordanceDF[["putativeDuplicateSample"]] <- comparisonDiscordanceDF[["Discordances"]] <= discordanceThreshold
#  
#  pdf(paste(plotFilestem, "comparisonDiscordanceHeatmap.pdf", sep="."), height=8, width=10)
#  print(
#    ggplot(
#      comparisonDiscordanceDF,
##      melt(comparisonVsSubjectDiscordanceMatrix, value.name="Discordances"),
#      aes(x=Var1, y=Var2, fill=Discordances)
#    )
#    + geom_tile()
#    + scale_fill_gradient2(low="red", high="blue", midpoint=median(comparisonVsSubjectDiscordanceMatrix, na.rm=TRUE))
#    + theme_bw()
#    + xlab(paste(comparisonName, "sample ID"))
#    + ylab(paste(comparisonName, "sample ID"))
#    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
##    + theme(axis.title.x = text(labels="Jiang et al. sample ID"))
##    + theme(axis.title.y = element_blank())
#  )
#  dev.off()
#
#  pdf(paste(plotFilestem, "comparisonPutativeIdenticalSamples.pdf", sep="."), height=8, width=10)
#  print(
#    ggplot(
#      comparisonDiscordanceDF,
#      aes(x=Var1, y=Var2, fill=putativeDuplicateSample)
#    )
#    + geom_tile()
#    + scale_fill_grey()
#    + theme_bw()
#    + xlab(paste(comparisonName, "sample ID"))
#    + ylab(paste(comparisonName, "sample ID"))
#    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
#  )
#  dev.off()
#
#  pdf(paste(plotFilestem, "comparisonDiscordanceProportions.pdf", sep="."), height=4, width=6)
#  print(
#    qplot(as.vector(comparisonVsComparisonDiscordanceProportionMatrix), xlab=paste("Proportion of discordant SNPs between pairs of", comparisonName, "samples"), ylab="Frequency (number of sample pairs)")
#    + geom_vline(xintercept = discordanceProportionThreshold, colour="red")
#    + theme_bw()
#  )
#  dev.off()
#  pdf(paste(plotFilestem, "comparisonDiscordanceProportionsDuplicates.pdf", sep="."), height=4, width=6)
#  print(
#    qplot(as.vector(comparisonVsComparisonDiscordanceProportionMatrix[comparisonVsComparisonDiscordanceProportionMatrix<discordanceProportionThreshold]), binwidth=1, xlab=paste("Proportion of discordant SNPs between pairs of", comparisonName, "samples"), ylab="Frequency (number of sample pairs)")
#    + theme_bw()
#  )
#  dev.off()
#
#  comparisonDiscordanceProportionDF <- melt(comparisonVsComparisonDiscordanceProportionMatrix, value.name="DiscordanceProportions")
#  comparisonDiscordanceProportionDF[["putativeDuplicateSample"]] <- comparisonDiscordanceProportionDF[["DiscordanceProportions"]] <= discordanceProportionThreshold
#  
#  pdf(paste(plotFilestem, "comparisonDiscordanceProportionHeatmap.pdf", sep="."), height=8, width=10)
#  print(
#    ggplot(
#      comparisonDiscordanceProportionDF,
##      melt(comparisonVsSubjectDiscordanceMatrix, value.name="Discordances"),
#      aes(x=Var1, y=Var2, fill=DiscordanceProportions)
#    )
#    + geom_tile()
#    + scale_fill_gradient2(low="red", high="blue", midpoint=median(comparisonVsSubjectDiscordanceProportionMatrix, na.rm=TRUE))
#    + theme_bw()
#    + xlab(paste(comparisonName, "sample ID"))
#    + ylab(paste(comparisonName, "sample ID"))
#    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
##    + theme(axis.title.x = text(labels="Jiang et al. sample ID"))
##    + theme(axis.title.y = element_blank())
#  )
#  dev.off()
#
#  pdf(paste(plotFilestem, "comparisonPutativeIdenticalSamplesFromProportions.pdf", sep="."), height=8, width=10)
#  print(
#    ggplot(
#      comparisonDiscordanceProportionDF,
#      aes(x=Var1, y=Var2, fill=putativeDuplicateSample)
#    )
#    + geom_tile()
#    + scale_fill_grey()
#    + theme_bw()
#    + xlab(paste(comparisonName, "sample ID"))
#    + ylab(paste(comparisonName, "sample ID"))
#    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
#  )
#  dev.off()
#
#  
#
#  mean(comparisonVsSubjectDiscordanceMatrix[comparisonVsSubjectDiscordanceMatrix<discordanceThreshold])
#  genotypeErrorRateIdentical <- sum(comparisonVsSubjectDiscordanceMatrix[comparisonVsSubjectDiscordanceMatrix<discordanceThreshold]) / (dim(comparisonMatchesGTs)[1] * length(which(comparisonVsSubjectDiscordanceMatrix<discordanceThreshold)))
#  genotypeConcordance <- 1-genotypeErrorRateIdentical
#  
#  duplicates <- which(comparisonVsSubjectDiscordanceMatrix<discordanceThreshold, arr.ind=TRUE)
#  duplicates <- duplicates[duplicates[, 1]!=duplicates[, 2], ]
#  diffs <- comparisonMatchesGTs[, duplicates[, 1]] != subjectMatchesGTs[, duplicates[, 2]]
#  table(rowSums(diffs))
#  table(subjectMatchesGTs[24, duplicates[, 2]], comparisonMatchesGTs[24, duplicates[, 1]], useNA="ifany")
###  which variants are monomorphic in progeny?
##  malariagenProgeny <- setdiff(dimnames(subjectMatchesGTs)[[2]], c("7G8_NIH/PG0083-C/ERR027099", "GB4_NIH/PG0084-C/ERR027100"))
##  table(apply(subjectMatchesGTs[, malariagenProgeny], 1, function(x) length(unique(x))), useNA="ifany")
##  jiangProgeny <- setdiff(dimnames(comparisonMatchesGTs)[[2]], c("_7G8", "GB4"))
##  table(apply(comparisonMatchesGTs[, jiangProgeny], 1, function(x) length(unique(x))), useNA="ifany")
##  No vairants are monomorphic in progeny - this was a dead end...
#
#  mostDiscordantSNPs <- which(rowSums(diffs) > 3)
#  sapply(
#    mostDiscordantSNPs,
#    function(SNPindex) {
#      names(which(subjectMatchesGTs[SNPindex, duplicates[, 2]]!=comparisonMatchesGTs[SNPindex, duplicates[, 1]]))
#    }
#  )
#  lapply(
#    mostDiscordantSNPs,
#    function(SNPindex) {
#      subjectMatchesGTs[(SNPindex-2):(SNPindex+2), ]
#    }
#  )
#      
##  browser()
#  return(
#    list(
#      comparisonVsSubjectDiscordanceMatrix = comparisonVsSubjectDiscordanceMatrix,
#      comparisonVsSubjectDiscordanceProportionMatrix = comparisonVsSubjectDiscordanceProportionMatrix,
#      subjectVsSubjectDiscordanceMatrix = subjectVsSubjectDiscordanceMatrix,
#      subjectVsSubjectDiscordanceProportionMatrix = subjectVsSubjectDiscordanceProportionMatrix,
#      comparisonVsComparisonDiscordanceMatrix = comparisonVsComparisonDiscordanceMatrix,
#      comparisonVsComparisonDiscordanceProportionMatrix = comparisonVsComparisonDiscordanceProportionMatrix
#    )
#  )
}

# callComparisonPlots.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


callComparisonPlots <- function(
  subjectMatchesGTs,
  comparisonMatchesGTs,
  comparison                  = c("discordances", "discordanceProportions"),
  threshold,
  plotFilestem,
  expectedMatches,
  subjectName,
  comparisonName,
  subjectIsComparison         = identical(subjectMatchesGTs, comparisonMatchesGTs) && subjectName==comparisonName,
  ...,
  discordanceHistogramHeight  = 4,
  discordanceHistogramWidth   = 6,
#  discordanceHistogramBinwidth= 50
  discordanceHeatmapHeight    = 10,
  discordanceHeatmapWidth     = 10,
  heatmapMargin               = 10,
  verbose                     = TRUE
)
{
  require(gplots)

  if(subjectIsComparison) {
    expectedMatches <- rbind(expectedMatches, cbind(expectedMatches[, 2], expectedMatches[, 1]))
  }
  comparison <- comparison[1]
  if(comparison == "discordances") {
    discordance <- function(x, y) length(which(x!=y))
    discordanceText <- "Number"
  } else if(comparison == "discordanceProportions") {
    discordance <- function(x, y) length(which(x!=y)) / length(which(!is.na(x) & !is.na(y)))
    discordanceText <- "Proportion"
  } else {
    stop("discordance must be either \"discordances\" or \"discordanceProportions\"")
  }
  vecDiscordance <- Vectorize(discordance)
  comparisonMatchesGTsDF <- data.frame(comparisonMatchesGTs)
  names(comparisonMatchesGTsDF) <- dimnames(comparisonMatchesGTs)[[2]]
  subjectMatchesGTsDF <- data.frame(subjectMatchesGTs)
  names(subjectMatchesGTsDF) <- dimnames(subjectMatchesGTs)[[2]]
#  comparisonVsSubjectDiscordanceMatrix <- outer(data.frame(comparisonMatchesGTs), data.frame(subjectMatchesGTs), vecDiscordance)
  comparisonVsSubjectDiscordanceMatrix <- outer(comparisonMatchesGTsDF, subjectMatchesGTsDF, vecDiscordance)

  discordancesText <- paste(
    "Median pairwise",
    comparison,
    "=",
    formatC(median(comparisonVsSubjectDiscordanceMatrix), 2, format="fg"),
    "\nMedian",
    comparison,
    "of putative matches = ",
    formatC(median(comparisonVsSubjectDiscordanceMatrix[comparisonVsSubjectDiscordanceMatrix<threshold]), 2, format="fg")
  )
  pdf(
    paste(plotFilestem, comparison, "histogram", "pdf", sep="."),
    height=discordanceHistogramHeight,
    width=discordanceHistogramWidth
  )
  if(subjectIsComparison) {
    xlabText <- paste(discordanceText, "of discordant SNPs between pairs of", subjectName, "samples")
  } else {
    xlabText <- paste(discordanceText, "of discordant SNPs between", comparisonName, "and", subjectName, "samples")
  }
  print(
    qplot(
      as.vector(comparisonVsSubjectDiscordanceMatrix),
#      binwidth=discordanceHistogramBinwidth,
      xlab=xlabText,
      ylab="Frequency (number of sample pairs)"
    )
    + geom_vline(xintercept = threshold, colour="red")
    + theme_bw()
    + labs(title=discordancesText)
  )
  dev.off()
  pdf(
    paste(plotFilestem, comparison, "putativeMatchesHistogram", "pdf", sep="."),
    height=discordanceHistogramHeight,
    width=discordanceHistogramWidth
  )
  print(
    qplot(
      as.vector(comparisonVsSubjectDiscordanceMatrix[comparisonVsSubjectDiscordanceMatrix<threshold]),
#      binwidth=1,
      xlab=xlabText,
      ylab="Frequency (number of sample pairs)"
    )
    + theme_bw()
  )
  dev.off()
#  if(!is.null(expectedMatches) && "comparisonVsSubject" %in% names(expectedMatches)) {
#    comparisonVsSubjectDiscordanceMatrixExpectMatches <- comparisonVsSubjectDiscordanceMatrix[expectedMatches[["comparisonVsSubject"]]]
  if(!is.null(expectedMatches)) {
    comparisonVsSubjectDiscordanceMatrixExpectMatches <- comparisonVsSubjectDiscordanceMatrix[expectedMatches]
    discordancesText <- paste(
      "Median",
      comparison,
      "in expected matches =",
      formatC(median(comparisonVsSubjectDiscordanceMatrixExpectMatches), 2, format="fg")
    )
    pdf(
      paste(plotFilestem, comparison, "expectedMatchesHistogram.pdf", sep="."),
      height=discordanceHistogramHeight,
      width=discordanceHistogramWidth
    )
    print(
      qplot(
          as.vector(comparisonVsSubjectDiscordanceMatrixExpectMatches),
          xlab=xlabText,
#          xlab=paste(discordanceText, "of discordant SNPs between", comparisonName, "and", subjectName, "expected match samples"),
          ylab="Frequency (number of sample pairs)"
      )
      + geom_vline(xintercept = threshold, colour="red")
      + theme_bw()
      + labs(title=discordancesText)
    )
    dev.off()
    discordancesText <- paste(
      "Median",
      comparison,
      "in true expected matches = ",
      formatC(median(comparisonVsSubjectDiscordanceMatrix[comparisonVsSubjectDiscordanceMatrix<threshold]), 2, format="fg")
    )
    pdf(
      paste(plotFilestem, comparison, "DuplicatesExpectedMatches.pdf", sep="."),
      height=discordanceHistogramHeight,
      width=discordanceHistogramWidth
    )
    print(
      qplot(
          as.vector(comparisonVsSubjectDiscordanceMatrixExpectMatches[comparisonVsSubjectDiscordanceMatrixExpectMatches<threshold]),
          xlab=xlabText,
#          xlab=paste(discordanceText, "of discordant SNPs between", comparisonName, "and", subjectName, "expected match samples"),
          ylab="Frequency (number of sample pairs)"
      )
      + theme_bw()
      + labs(title=discordancesText)
    )
    dev.off()
  }
  
  discordanceDF <- melt(comparisonVsSubjectDiscordanceMatrix, value.name="Discordances")
  discordanceDF[["putativeDuplicateSample"]] <- discordanceDF[["Discordances"]] <= threshold
  discordanceDF[["putativeDuplicateSample"]][is.na(discordanceDF[["putativeDuplicateSample"]])] <- FALSE
  if(!is.null(expectedMatches)) {
    comparisonVsSubjectExpectedMatches <- comparisonVsSubjectDiscordanceMatrix
    comparisonVsSubjectExpectedMatches[,] <- 0
    comparisonVsSubjectExpectedMatches[expectedMatches] <- 1
    expectedMatchesDF <- melt(comparisonVsSubjectExpectedMatches, value.name="ExpectedMatch")
    discordanceDF[["expectedDuplicateSample"]] <- as.logical(expectedMatchesDF[["ExpectedMatch"]])
    discordanceDF[discordanceDF[["expectedDuplicateSample"]] & discordanceDF[["putativeDuplicateSample"]], "SamplePairStatus"] <- "Expected match"
    discordanceDF[!discordanceDF[["expectedDuplicateSample"]] & discordanceDF[["putativeDuplicateSample"]], "SamplePairStatus"] <- "Unexpected match"
    discordanceDF[discordanceDF[["expectedDuplicateSample"]] & !discordanceDF[["putativeDuplicateSample"]], "SamplePairStatus"] <- "Unexpected non-match"
    discordanceDF[!discordanceDF[["expectedDuplicateSample"]] & !discordanceDF[["putativeDuplicateSample"]], "SamplePairStatus"] <- "Expected non-match"
    discordanceDF[is.na(discordanceDF[["Discordances"]]), "SamplePairStatus"] <- "No matching SNPs"
    if(subjectIsComparison) {
      discordanceDF[discordanceDF[["Var1"]] == discordanceDF[["Var2"]], "SamplePairStatus"] <- "Identical sample"
    }
#    discordanceDF[["SamplePairStatus"]] <- factor(discordanceDF[["SamplePairStatus"]], levels=c("Expected match", "Expected non-match", "Unexpected match", "Unexpected non-match", "No matching SNPs", "Identical sample"))
#    browser()
  }
  
  pdf(
    paste(plotFilestem, comparison, "Heatmap.pdf", sep="."),
    height=discordanceHeatmapHeight,
    width=discordanceHeatmapWidth
  )
  print(
    ggplot(
      discordanceDF,
      aes(x=Var1, y=Var2, fill=Discordances)
    )
    + geom_tile()
    + scale_fill_gradient2(low="red", high="blue", midpoint=median(comparisonVsSubjectDiscordanceMatrix, na.rm=TRUE))
    + theme_bw()
    + xlab(paste(comparisonName, "sample ID"))
    + ylab(paste(subjectName, "sample ID"))
    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  )
  dev.off()

  pdf(
    paste(plotFilestem, comparison, "Heatmap2.pdf", sep="."),
    height=discordanceHeatmapHeight,
    width=discordanceHeatmapWidth
  )
  heatmap.2(comparisonVsSubjectDiscordanceMatrix, margins=c(heatmapMargin, heatmapMargin))
  dev.off()
  
  print(
    ggplot(
      discordanceDF,
      aes(x=Var1, y=Var2, fill=Discordances)
    )
    + geom_tile()
    + scale_fill_gradient2(low="red", high="blue", midpoint=median(comparisonVsSubjectDiscordanceMatrix, na.rm=TRUE))
    + theme_bw()
    + xlab(paste(comparisonName, "sample ID"))
    + ylab(paste(subjectName, "sample ID"))
    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  )
  dev.off()

#  pdf(
#    paste(plotFilestem, comparison, "putativeIdenticalSamples.pdf", sep="."),
#    height=10,
#    width=10
#  )
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
  
  if(!is.null(expectedMatches)) {
    pdf(
      paste(plotFilestem, comparison, "samplePairStatus.pdf", sep="."),
      height=discordanceHeatmapHeight,
      width=discordanceHeatmapWidth
    )
    print(
      ggplot(
        discordanceDF,
        aes(x=Var1, y=Var2, fill=SamplePairStatus)
      )
      + geom_tile()
      + scale_fill_manual(values = c("Expected match"="green", "Expected non-match"="grey90", "Unexpected match"="orange", "Unexpected non-match"="red", "No matching SNPs"="black", "Identical sample"="white"))
      + theme_bw()
      + xlab(paste(comparisonName, "sample ID"))
      + ylab(paste(subjectName, "sample ID"))
      + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
    )
    dev.off()
  }
  
  return(comparisonVsSubjectDiscordanceMatrix)

}

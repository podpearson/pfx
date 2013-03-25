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
  malariagenDiscordanceThreshold        = 200,
  sampleAnnotationFilename    = "/data/malariagen2/plasmodium/pf-crosses/meta/qcmeta_annotated.tsv",
  plotFilestem                = paste(meta(exptData(subjectVcf)[["header"]])["DataSetName", "Value"], "comparison", sep="."),
  shouldRenameSubjectSamples  = TRUE,
  IDparent1                   = dimnames(subjectVcf)[[2]][1],
  IDparent2                   = dimnames(subjectVcf)[[2]][2]
) {
  require(reshape2)
  if(shouldRenameSubjectSamples) {
    subjectVcf <- renameSamples(subjectVcf)
  }
  sampleAnnotation <- readSampleAnnotation(sampleAnnotationFilename)
#  sampleAnnotation <- read.delim(sampleAnnotationFilename)
#  sampleAnnotation[["cross"]] <- ifelse(sampleAnnotation[["project_code"]]=="PFproj1", "Hb3xDd2", ifelse(sampleAnnotation[["project_code"]]=="PFproj2", "3d7xHb3", "7g8xGb4"))
#  sampleAnnotation[["qcStatus"]] <- ifelse(sampleAnnotation[["ajm_qc"]]=="fail", "fail", "pass")
#  sampleAnnotation <- sampleAnnotation[!grepl("not a cross progeny sample", sampleAnnotation[["ajm_qc_notes"]]), ]
#  sampleAnnotation <- sampleAnnotation[!duplicated(sampleAnnotation[["ox_code"]]), ]
#  row.names(sampleAnnotation) <- sampleAnnotation[["ox_code"]]
##  row.names(sampleAnnotation) <- gsub("-", "\\.", sampleAnnotation[["ox_code"]])
  
  subjectGTsCFparents <- convertGTsIntToParentBasedGTs(genotypeCallsFromGTas012(subjectVcf), IDparent1=IDparent1, IDparent2=IDparent2)
#  browser()
#  dimnames(subjectGTsCFparents)[[2]] <- paste(sampleAnnotation[dimnames(subjectGTsCFparents)[[2]], "source_code"], " (", dimnames(subjectGTsCFparents)[[2]], ")", sep="")
  comparisonGTsCFparents <- convertGTsIntToParentBasedGTs(genotypeCallsFromGTas012(comparisonVcf, GTsToIntMapping = c("7"=1, "G"=2, "."=0)))
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
  pdf(paste(plotFilestem, "positionDifferencesWithThresholds.pdf", sep="."), height=4, width=6)
  print(
    qplot(subjectVsComparisonPositionDifferences, binwidth=1, xlim=c(-50, 50), xlab="Distance (bp) from", comparisonName, "SNP to nearest", subjectName, "SNP", ylab="Frequency (number of SNPs)")
    + geom_vline(xintercept = distanceThresholds[1], colour="red")
    + geom_vline(xintercept = distanceThresholds[2], colour="red")
    + theme_bw()
  )
  dev.off()

  comparisonRowsWithMatchesInSubject <- which(subjectVsComparisonPositionDifferences >=distanceThresholds[1] & subjectVsComparisonPositionDifferences<=distanceThresholds[2])
  subjectRowsWithMatchesInComparison <- nearest(rowData(comparisonVcf), rowData(subjectVcf))[comparisonRowsWithMatchesInSubject]
  comparisonMatchesGTs <- comparisonGTsCFparents[comparisonRowsWithMatchesInSubject,]
  subjectMatchesGTs <- subjectGTsCFparents[subjectRowsWithMatchesInComparison, ]
#  jiangGTsGWInt <- -(matrix(as.integer(factor(jiangMatchesGTsGW)), nrow=nrow(jiangMatchesGTsGW))-2)
#  dimnames(jiangGTsGWInt) <- dimnames(jiangMatchesGTsGW)
  
  discordance <- function(x, y) length(which(x!=y))
  vecDiscordance <- Vectorize(discordance)
  comparisonVsSubjectDiscordanceMatrix <- outer(data.frame(comparisonMatchesGTs), data.frame(subjectMatchesGTs), vecDiscordance)

  pdf(paste(plotFilestem, "discordances.pdf", sep="."), height=4, width=6)
  print(
    qplot(as.vector(comparisonVsSubjectDiscordanceMatrix), binwidth=50, xlab=paste("Number of discordant SNPs between", comparisonName, "and", subjectName, "samples"), ylab="Frequency (number of sample pairs)")
    + geom_vline(xintercept = discordanceThreshold, colour="red")
    + theme_bw()
  )
  dev.off()
  pdf(paste(plotFilestem, "discordancesDuplicates.pdf", sep="."), height=4, width=6)
  print(
    qplot(as.vector(comparisonVsSubjectDiscordanceMatrix[comparisonVsSubjectDiscordanceMatrix<discordanceThreshold]), binwidth=1, xlab=paste("Number of discordant SNPs between", comparisonName, "and", subjectName, "samples"), ylab="Frequency (number of sample pairs)")
    + theme_bw()
  )
  dev.off()

  discordanceDF <- melt(comparisonVsSubjectDiscordanceMatrix, value.name="Discordances")
  discordanceDF[["putativeDuplicateSample"]] <- discordanceDF[["Discordances"]] <= discordanceThreshold
  
  pdf(paste(plotFilestem, "discordanceHeatmap.pdf", sep="."), height=10, width=10)
  print(
    ggplot(
      discordanceDF,
#      melt(comparisonVsSubjectDiscordanceMatrix, value.name="Discordances"),
      aes(x=Var1, y=Var2, fill=Discordances)
    )
    + geom_tile()
    + scale_fill_gradient2(low="red", high="blue", midpoint=median(comparisonVsSubjectDiscordanceMatrix, na.rm=TRUE))
    + theme_bw()
    + xlab(paste(comparisonName, "sample ID"))
    + ylab(paste(subjectName, "sample ID"))
    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
#    + theme(axis.title.x = text(labels="Jiang et al. sample ID"))
#    + theme(axis.title.y = element_blank())
  )
  dev.off()

  pdf(paste(plotFilestem, "putativeIdenticalSamples.pdf", sep="."), height=10, width=10)
  print(
    ggplot(
      discordanceDF,
      aes(x=Var1, y=Var2, fill=putativeDuplicateSample)
    )
    + geom_tile()
    + scale_fill_grey()
    + theme_bw()
    + xlab(paste(comparisonName, "sample ID"))
    + ylab(paste(subjectName, "sample ID"))
    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  )
  dev.off()
  
  
  
  subjectVsSubjectDiscordanceMatrix <- outer(data.frame(subjectMatchesGTs), data.frame(subjectMatchesGTs), vecDiscordance)

  pdf(paste(plotFilestem, "malariagenDiscordances.pdf", sep="."), height=4, width=6)
  print(
    qplot(as.vector(subjectVsSubjectDiscordanceMatrix), binwidth=50, xlab=paste("Number of discordant SNPs between pairs of", subjectName, "samples"), ylab="Frequency (number of sample pairs)")
    + geom_vline(xintercept = malariagenDiscordanceThreshold, colour="red")
    + theme_bw()
  )
  dev.off()
  pdf(paste(plotFilestem, "malariagenDiscordancesDuplicates.pdf", sep="."), height=4, width=6)
  print(
    qplot(as.vector(subjectVsSubjectDiscordanceMatrix[subjectVsSubjectDiscordanceMatrix<malariagenDiscordanceThreshold]), binwidth=1, xlab=paste("Number of discordant SNPs between pairs of", subjectName, "samples"), ylab="Frequency (number of sample pairs)")
    + theme_bw()
  )
  dev.off()

  subjectDiscordanceDF <- melt(subjectVsSubjectDiscordanceMatrix, value.name="Discordances")
  subjectDiscordanceDF[["putativeDuplicateSample"]] <- subjectDiscordanceDF[["Discordances"]] <= malariagenDiscordanceThreshold
  
  pdf(paste(plotFilestem, "malariagenDiscordanceHeatmap.pdf", sep="."), height=10, width=12)
  print(
    ggplot(
      subjectDiscordanceDF,
#      melt(comparisonVsSubjectDiscordanceMatrix, value.name="Discordances"),
      aes(x=Var1, y=Var2, fill=Discordances)
    )
    + geom_tile()
    + scale_fill_gradient2(low="red", high="blue", midpoint=median(subjectVsSubjectDiscordanceMatrix, na.rm=TRUE))
    + theme_bw()
    + xlab(paste(subjectName, "sample ID"))
    + ylab(paste(subjectName, "sample ID"))
    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  )
  dev.off()

  pdf(paste(plotFilestem, "malariagenPutativeIdenticalSamples.pdf", sep="."), height=10, width=12)
  print(
    ggplot(
      subjectDiscordanceDF,
      aes(x=Var1, y=Var2, fill=putativeDuplicateSample)
    )
    + geom_tile()
    + scale_fill_grey()
    + theme_bw()
    + xlab(paste(subjectName, "sample ID"))
    + ylab(paste(subjectName, "sample ID"))
    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  )
  dev.off()
  
  
  
  comparisonVsComparisonDiscordanceMatrix <- outer(data.frame(comparisonMatchesGTs), data.frame(comparisonMatchesGTs), vecDiscordance)

  pdf(paste(plotFilestem, "comparisonDiscordances.pdf", sep="."), height=4, width=6)
  print(
    qplot(as.vector(comparisonVsComparisonDiscordanceMatrix), binwidth=50, xlab=paste("Number of discordant SNPs between pairs of", comparisonName, "samples"), ylab="Frequency (number of sample pairs)")
    + geom_vline(xintercept = discordanceThreshold, colour="red")
    + theme_bw()
  )
  dev.off()
  pdf(paste(plotFilestem, "comparisonDiscordancesDuplicates.pdf", sep="."), height=4, width=6)
  print(
    qplot(as.vector(comparisonVsComparisonDiscordanceMatrix[comparisonVsComparisonDiscordanceMatrix<discordanceThreshold]), binwidth=1, xlab=paste("Number of discordant SNPs between pairs of", comparisonName, "samples"), ylab="Frequency (number of sample pairs)")
    + theme_bw()
  )
  dev.off()

  comparisonDiscordanceDF <- melt(comparisonVsComparisonDiscordanceMatrix, value.name="Discordances")
  comparisonDiscordanceDF[["putativeDuplicateSample"]] <- comparisonDiscordanceDF[["Discordances"]] <= discordanceThreshold
  
  pdf(paste(plotFilestem, "comparisonDiscordanceHeatmap.pdf", sep="."), height=8, width=10)
  print(
    ggplot(
      comparisonDiscordanceDF,
#      melt(comparisonVsSubjectDiscordanceMatrix, value.name="Discordances"),
      aes(x=Var1, y=Var2, fill=Discordances)
    )
    + geom_tile()
    + scale_fill_gradient2(low="red", high="blue", midpoint=median(comparisonVsSubjectDiscordanceMatrix, na.rm=TRUE))
    + theme_bw()
    + xlab(paste(comparisonName, "sample ID"))
    + ylab(paste(comparisonName, "sample ID"))
    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
#    + theme(axis.title.x = text(labels="Jiang et al. sample ID"))
#    + theme(axis.title.y = element_blank())
  )
  dev.off()

  pdf(paste(plotFilestem, "comparisonPutativeIdenticalSamples.pdf", sep="."), height=8, width=10)
  print(
    ggplot(
      comparisonDiscordanceDF,
      aes(x=Var1, y=Var2, fill=putativeDuplicateSample)
    )
    + geom_tile()
    + scale_fill_grey()
    + theme_bw()
    + xlab(paste(comparisonName, "sample ID"))
    + ylab(paste(comparisonName, "sample ID"))
    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  )
  dev.off()

  

  mean(comparisonVsSubjectDiscordanceMatrix[comparisonVsSubjectDiscordanceMatrix<discordanceThreshold])
  genotypeErrorRateIdentical <- sum(comparisonVsSubjectDiscordanceMatrix[comparisonVsSubjectDiscordanceMatrix<discordanceThreshold]) / (dim(comparisonMatchesGTs)[1] * length(which(comparisonVsSubjectDiscordanceMatrix<discordanceThreshold)))
  genotypeConcordance <- 1-genotypeErrorRateIdentical
  
  duplicates <- which(comparisonVsSubjectDiscordanceMatrix<discordanceThreshold, arr.ind=TRUE)
  duplicates <- duplicates[duplicates[, 1]!=duplicates[, 2], ]
  diffs <- comparisonMatchesGTs[, duplicates[, 1]] != subjectMatchesGTs[, duplicates[, 2]]
  table(rowSums(diffs))
  table(subjectMatchesGTs[24, duplicates[, 2]], comparisonMatchesGTs[24, duplicates[, 1]], useNA="ifany")
##  which variants are monomorphic in progeny?
#  malariagenProgeny <- setdiff(dimnames(subjectMatchesGTs)[[2]], c("7G8_NIH/PG0083-C/ERR027099", "GB4_NIH/PG0084-C/ERR027100"))
#  table(apply(subjectMatchesGTs[, malariagenProgeny], 1, function(x) length(unique(x))), useNA="ifany")
#  jiangProgeny <- setdiff(dimnames(comparisonMatchesGTs)[[2]], c("_7G8", "GB4"))
#  table(apply(comparisonMatchesGTs[, jiangProgeny], 1, function(x) length(unique(x))), useNA="ifany")
#  No vairants are monomorphic in progeny - this was a dead end...

  mostDiscordantSNPs <- which(rowSums(diffs) > 3)
  sapply(
    mostDiscordantSNPs,
    function(SNPindex) {
      names(which(subjectMatchesGTs[SNPindex, duplicates[, 2]]!=comparisonMatchesGTs[SNPindex, duplicates[, 1]]))
    }
  )
  lapply(
    mostDiscordantSNPs,
    function(SNPindex) {
      subjectMatchesGTs[(SNPindex-2):(SNPindex+2), ]
    }
  )
      
#  browser()
  return(comparisonVsSubjectDiscordanceMatrix)
}

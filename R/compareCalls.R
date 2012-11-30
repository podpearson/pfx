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
  distanceThresholds          = c(0, 22),
  discordanceThreshold        = 100,
  sampleAnnotationFilename    = "/data/malariagen2/plasmodium/pf-crosses/meta/qcmeta_annotated.tsv",
  plotFilestem                = paste(meta(exptData(subjectVcf)[["header"]])["DataSetName", "Value"], "comparison", sep=".")
) {
  sampleAnnotation <- readSampleAnnotation(sampleAnnotationFilename)
#  sampleAnnotation <- read.delim(sampleAnnotationFilename)
#  sampleAnnotation[["cross"]] <- ifelse(sampleAnnotation[["project_code"]]=="PFproj1", "Hb3xDd2", ifelse(sampleAnnotation[["project_code"]]=="PFproj2", "3d7xHb3", "7g8xGb4"))
#  sampleAnnotation[["qcStatus"]] <- ifelse(sampleAnnotation[["ajm_qc"]]=="fail", "fail", "pass")
#  sampleAnnotation <- sampleAnnotation[!grepl("not a cross progeny sample", sampleAnnotation[["ajm_qc_notes"]]), ]
#  sampleAnnotation <- sampleAnnotation[!duplicated(sampleAnnotation[["ox_code"]]), ]
#  row.names(sampleAnnotation) <- sampleAnnotation[["ox_code"]]
##  row.names(sampleAnnotation) <- gsub("-", "\\.", sampleAnnotation[["ox_code"]])
  
  subjectGTsCFparents <- convertGTsIntToParentBasedGTs(genotypeCallsFromGTas012(subjectVcf))
  dimnames(subjectGTsCFparents)[[2]] <- paste(sampleAnnotation[dimnames(subjectGTsCFparents)[[2]], "source_code"], " (", dimnames(subjectGTsCFparents)[[2]], ")", sep="")
  comparisonGTsCFparents <- convertGTsIntToParentBasedGTs(genotypeCallsFromGTas012(comparisonVcf, GTsToIntMapping = c("7"=1, "G"=2, "."=0)))
  comparisonNearestSubjectGR <- rowData(subjectVcf)[nearest(rowData(comparisonVcf), rowData(subjectVcf))]
  table(start(ranges(comparisonNearestSubjectGR)) - start(ranges(rowData(comparisonVcf))))
  table(start(ranges(comparisonNearestSubjectGR))-start(ranges(rowData(comparisonVcf))) >=0 & start(ranges(comparisonNearestSubjectGR))-start(ranges(rowData(comparisonVcf)))<=22)
  subjectVsComparisonPositionDifferences <- start(ranges(comparisonNearestSubjectGR))-start(ranges(rowData(comparisonVcf)))
  require(ggplot2)
  pdf(paste(plotFilestem, "positionDifferences.pdf", sep="."), height=4, width=6)
  print(
    qplot(subjectVsComparisonPositionDifferences, binwidth=1, xlim=c(-50, 50), xlab="Distance (bp) from Jiang et al SNP to nearest MalariaGEN SNP", ylab="Frequency (number of SNPs)")
    + theme_bw()
  )
  dev.off()
  pdf(paste(plotFilestem, "positionDifferencesWithThresholds.pdf", sep="."), height=4, width=6)
  print(
    qplot(subjectVsComparisonPositionDifferences, binwidth=1, xlim=c(-50, 50), xlab="Distance (bp) from Jiang et al SNP to nearest MalariaGEN SNP", ylab="Frequency (number of SNPs)")
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
    qplot(as.vector(comparisonVsSubjectDiscordanceMatrix), binwidth=50, xlab="Number of discordant SNPs between Jiang et al and malariagen samples", ylab="Frequency (number of sample pairs)")
    + geom_vline(xintercept = discordanceThreshold, colour="red")
    + theme_bw()
  )
  dev.off()
  pdf(paste(plotFilestem, "discordancesDuplicates.pdf", sep="."), height=4, width=6)
  print(
    qplot(as.vector(comparisonVsSubjectDiscordanceMatrix[comparisonVsSubjectDiscordanceMatrix<discordanceThreshold]), binwidth=1, xlab="Number of discordant SNPs between Jiang et al and malariagen samples", ylab="Frequency (number of sample pairs)")
    + theme_bw()
  )
  dev.off()

  pdf(paste(plotFilestem, "discordanceHeatmap.pdf", sep="."), height=6, width=10)
  print(
    ggplot(
      melt(comparisonVsSubjectDiscordanceMatrix, value.name="Discordances"),
      aes(x=Var1, y=Var2, fill=Discordances)
    )
    + geom_tile()
    + scale_fill_gradient2(low="red", high="blue", midpoint=median(comparisonVsSubjectDiscordanceMatrix, na.rm=TRUE))
    + theme_bw()
    + xlab("Jiang et al. sample ID")
    + ylab("Malariagen sample ID")
    + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
#    + theme(axis.title.x = text(labels="Jiang et al. sample ID"))
#    + theme(axis.title.y = element_blank())
  )
  dev.off()

  mean(comparisonVsSubjectDiscordanceMatrix[comparisonVsSubjectDiscordanceMatrix<discordanceThreshold])
  genotypeErrorRateIdentical <- sum(comparisonVsSubjectDiscordanceMatrix[comparisonVsSubjectDiscordanceMatrix<discordanceThreshold]) / (dim(comparisonMatchesGTs)[1] * length(which(comparisonVsSubjectDiscordanceMatrix<discordanceThreshold)))
  genotypeConcordance <- 1-genotypeErrorRateIdentical
  return(genotypeConcordance)
}

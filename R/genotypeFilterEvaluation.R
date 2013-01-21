# genotypeFilterEvaluation.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


#  MAF_genotypeFilterEvaluation_3d7_hb3_snps <- genotypeFilterEvaluation("3d7_hb3", "snps")
#  MAF_genotypeFilterEvaluation_7g8_gb4_snps <- genotypeFilterEvaluation("7g8_gb4", "snps", parentalIDs = c("ERR027099", "ERR027100", "ERR029410", "ERR045626"))
#  MAF_genotypeFilterEvaluation_hb3_dd2_snps <- genotypeFilterEvaluation("hb3_dd2", "snps", parentalIDs = c("ERR012788", "ERR012840", "ERR022939"))
#  MAF_genotypeFilterEvaluation_3d7_hb3_indels <- genotypeFilterEvaluation("3d7_hb3", "indels", minMeanMAFtoConsiderContam=0.02)
#  MAF_genotypeFilterEvaluation_7g8_gb4_indels <- genotypeFilterEvaluation("7g8_gb4", "indels", minMeanMAFtoConsiderContam=0.05, parentalIDs = c("ERR027099", "ERR027100", "ERR029410", "ERR045626"))
#  MAF_genotypeFilterEvaluation_hb3_dd2_indels <- genotypeFilterEvaluation("hb3_dd2", "indels", minMeanMAFtoConsiderContam=0.05, parentalIDs = c("ERR012788", "ERR012840", "ERR022939"))

#  DP_genotypeFilterEvaluation_3d7_hb3_snps <- genotypeFilterEvaluation("3d7_hb3", "snps", MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1, 0))
#  DP_genotypeFilterEvaluation_7g8_gb4_snps <- genotypeFilterEvaluation("7g8_gb4", "snps", parentalIDs = c("ERR027099", "ERR027100", "ERR029410", "ERR045626"), MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1, 0))
#  DP_genotypeFilterEvaluation_hb3_dd2_snps <- genotypeFilterEvaluation("hb3_dd2", "snps", parentalIDs = c("ERR012788", "ERR012840", "ERR022939"), MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1, 0))
#  DP_genotypeFilterEvaluation_3d7_hb3_indels <- genotypeFilterEvaluation("3d7_hb3", "indels", minMeanMAFtoConsiderContam=0.02, MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1, 0))
#  DP_genotypeFilterEvaluation_7g8_gb4_indels <- genotypeFilterEvaluation("7g8_gb4", "indels", minMeanMAFtoConsiderContam=0.05, parentalIDs = c("ERR027099", "ERR027100", "ERR029410", "ERR045626"), MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1, 0))
#  DP_genotypeFilterEvaluation_hb3_dd2_indels <- genotypeFilterEvaluation("hb3_dd2", "indels", minMeanMAFtoConsiderContam=0.05, parentalIDs = c("ERR012788", "ERR012840", "ERR022939"), MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1, 0))

#  GQ_genotypeFilterEvaluation_3d7_hb3_snps <- genotypeFilterEvaluation("3d7_hb3", "snps", MAFthresholds=NULL, GQthresholds=c(99, 0), plotFilestemExtra="GQ")
#  GQ_genotypeFilterEvaluation_7g8_gb4_snps <- genotypeFilterEvaluation("7g8_gb4", "snps", parentalIDs = c("ERR027099", "ERR027100", "ERR029410", "ERR045626"), MAFthresholds=NULL, GQthresholds=c(99, 0), plotFilestemExtra="GQ")
#  GQ_genotypeFilterEvaluation_hb3_dd2_snps <- genotypeFilterEvaluation("hb3_dd2", "snps", parentalIDs = c("ERR012788", "ERR012840", "ERR022939"), MAFthresholds=NULL, GQthresholds=c(99, 0), plotFilestemExtra="GQ")
#  GQ_genotypeFilterEvaluation_3d7_hb3_indels <- genotypeFilterEvaluation("3d7_hb3", "indels", minMeanMAFtoConsiderContam=0.02, MAFthresholds=NULL, GQthresholds=c(99, 0), plotFilestemExtra="GQ")
#  GQ_genotypeFilterEvaluation_7g8_gb4_indels <- genotypeFilterEvaluation("7g8_gb4", "indels", minMeanMAFtoConsiderContam=0.05, parentalIDs = c("ERR027099", "ERR027100", "ERR029410", "ERR045626"), MAFthresholds=NULL, GQthresholds=c(99, 0), plotFilestemExtra="GQ")
#  GQ_genotypeFilterEvaluation_hb3_dd2_indels <- genotypeFilterEvaluation("hb3_dd2", "indels", minMeanMAFtoConsiderContam=0.05, parentalIDs = c("ERR012788", "ERR012840", "ERR022939"), MAFthresholds=NULL, GQthresholds=c(99, 0), plotFilestemExtra="GQ")

#  DP_QUALbyDP5_genotypeFilterEvaluation_3d7_hb3_snps <- genotypeFilterEvaluation("3d7_hb3", "snps", MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1, 0), filters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), plotFilestemExtra="DP_QUALbyDP5")
#  DP_QUALbyDP5_genotypeFilterEvaluation_7g8_gb4_snps <- genotypeFilterEvaluation("7g8_gb4", "snps", parentalIDs = c("ERR027099", "ERR027100", "ERR029410", "ERR045626"), MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1, 0), filters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), plotFilestemExtra="DP_QUALbyDP5")
#  DP_QUALbyDP5_genotypeFilterEvaluation_hb3_dd2_snps <- genotypeFilterEvaluation("hb3_dd2", "snps", parentalIDs = c("ERR012788", "ERR012840", "ERR022939"), MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1, 0), filters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), plotFilestemExtra="DP_QUALbyDP5")
#  DP_QUALbyDP5_genotypeFilterEvaluation_3d7_hb3_indels <- genotypeFilterEvaluation("3d7_hb3", "indels", minMeanMAFtoConsiderContam=0.02, MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1, 0), filters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), plotFilestemExtra="DP_QUALbyDP5")
#  DP_QUALbyDP5_genotypeFilterEvaluation_7g8_gb4_indels <- genotypeFilterEvaluation("7g8_gb4", "indels", minMeanMAFtoConsiderContam=0.05, parentalIDs = c("ERR027099", "ERR027100", "ERR029410", "ERR045626"), MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1, 0), filters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), plotFilestemExtra="DP_QUALbyDP5")
#  DP_QUALbyDP5_genotypeFilterEvaluation_hb3_dd2_indels <- genotypeFilterEvaluation("hb3_dd2", "indels", minMeanMAFtoConsiderContam=0.05, parentalIDs = c("ERR012788", "ERR012840", "ERR022939"), MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1, 0), filters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), plotFilestemExtra="DP_QUALbyDP5")

#  DP_QUALbyDP5_genotypeFilterEvaluation_3d7_hb3_snps_3 <- genotypeFilterEvaluation("3d7_hb3", "snps", MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1), filters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), plotFilestemExtra="DP_QUALbyDP5", monomorphicSkipChromosomes  = "Pf3D7_13_v3")
#  DP_QUALbyDP5_genotypeFilterEvaluation_7g8_gb4_snps_3 <- genotypeFilterEvaluation("7g8_gb4", "snps", parentalIDs = c("ERR027099", "ERR027100", "ERR029410", "ERR045626"), MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1), filters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), plotFilestemExtra="DP_QUALbyDP5")
#  DP_QUALbyDP5_genotypeFilterEvaluation_hb3_dd2_snps_3 <- genotypeFilterEvaluation("hb3_dd2", "snps", parentalIDs = c("ERR012788", "ERR012840", "ERR022939"), MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1), filters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), plotFilestemExtra="DP_QUALbyDP5")
#  DP_QUALbyDP5_genotypeFilterEvaluation_3d7_hb3_indels_3 <- genotypeFilterEvaluation("3d7_hb3", "indels", minMeanMAFtoConsiderContam=0.02, MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1), filters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), plotFilestemExtra="DP_QUALbyDP5", monomorphicSkipChromosomes  = "Pf3D7_13_v3")
#  DP_QUALbyDP5_genotypeFilterEvaluation_7g8_gb4_indels_3 <- genotypeFilterEvaluation("7g8_gb4", "indels", minMeanMAFtoConsiderContam=0.05, parentalIDs = c("ERR027099", "ERR027100", "ERR029410", "ERR045626"), MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1), filters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), plotFilestemExtra="DP_QUALbyDP5")
#  DP_QUALbyDP5_genotypeFilterEvaluation_hb3_dd2_indels_3 <- genotypeFilterEvaluation("hb3_dd2", "indels", minMeanMAFtoConsiderContam=0.05, parentalIDs = c("ERR012788", "ERR012840", "ERR022939"), MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1), filters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), plotFilestemExtra="DP_QUALbyDP5")

#  DP_QUALbyDP5_genotypeFilterEvaluation_3d7_hb3_snps_4 <- genotypeFilterEvaluation("3d7_hb3", "snps", MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1), filters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), plotFilestemExtra="DP_QUALbyDP5", monomorphicSkipChromosomes  = "Pf3D7_13_v3")
#  DP_QUALbyDP5_genotypeFilterEvaluation_7g8_gb4_snps_4 <- genotypeFilterEvaluation("7g8_gb4", "snps", parentalIDs = c("ERR027099", "ERR027100", "ERR029410", "ERR045626"), MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1), filters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), plotFilestemExtra="DP_QUALbyDP5")
#  DP_QUALbyDP5_genotypeFilterEvaluation_hb3_dd2_snps_4 <- genotypeFilterEvaluation("hb3_dd2", "snps", parentalIDs = c("ERR012788", "ERR012840", "ERR022939"), MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1), filters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), plotFilestemExtra="DP_QUALbyDP5")
#  DP_QUALbyDP5_genotypeFilterEvaluation_3d7_hb3_indels_4 <- genotypeFilterEvaluation("3d7_hb3", "indels", minMeanMAFtoConsiderContam=0.02, MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1), filters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), plotFilestemExtra="DP_QUALbyDP5", monomorphicSkipChromosomes  = "Pf3D7_13_v3")
#  DP_QUALbyDP5_genotypeFilterEvaluation_7g8_gb4_indels_4 <- genotypeFilterEvaluation(
#    "7g8_gb4",
#    "indels",
#    minMeanMAFtoConsiderContam=0.05,
#    parentalIDs = c("ERR027099", "ERR027100", "ERR029410", "ERR045626"),
#    MAFthresholds=NULL,
#    DPthresholds=c(5),
#    filters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)),
#    plotFilestemExtra="DP_QUALbyDP5",
#    shouldCreateQCFilteringPlots=TRUE
##    sampleSets=list("FinalSamples" = list(additionalInfoFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), sampleSetName="final"))
#  )
#  DP_QUALbyDP5_genotypeFilterEvaluation_hb3_dd2_indels_4 <- genotypeFilterEvaluation(
#    "hb3_dd2",
#    "indels",
#    minMeanMAFtoConsiderContam=0.05,
#    parentalIDs = c("ERR012788", "ERR012840", "ERR022939"),
#    MAFthresholds=NULL,
#    DPthresholds=c(5),
#    filters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)),
#    plotFilestemExtra="DP_QUALbyDP5",
#    shouldCreateQCFilteringPlots=TRUE
##    sampleSets=list("FinalSamples" = list(additionalInfoFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), sampleSetName="final"))
#  )
#  DP_QUALbyDP5_genotypeFilterEvaluation_hb3_dd2_indels_4 <- genotypeFilterEvaluation("hb3_dd2", "indels", minMeanMAFtoConsiderContam=0.05, parentalIDs = c("ERR012788", "ERR012840", "ERR022939"), MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1), filters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), plotFilestemExtra="DP_QUALbyDP5")

  DP_QUALbyDP5_genotypeFilterEvaluation_3d7_hb3_snps_5 <- genotypeFilterEvaluation("3d7_hb3", "snps", MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1), plotFilestemExtra="DP_QUALbyDP5", monomorphicSkipChromosomes  = "Pf3D7_13_v3")
  DP_QUALbyDP5_genotypeFilterEvaluation_7g8_gb4_snps_5 <- genotypeFilterEvaluation("7g8_gb4", "snps", parentalIDs = c("ERR027099", "ERR027100", "ERR029410", "ERR045626"), MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1), plotFilestemExtra="DP_QUALbyDP5")
  DP_QUALbyDP5_genotypeFilterEvaluation_hb3_dd2_snps_5 <- genotypeFilterEvaluation("hb3_dd2", "snps", parentalIDs = c("ERR012788", "ERR012840", "ERR022939"), MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1), plotFilestemExtra="DP_QUALbyDP5")
  DP_QUALbyDP5_genotypeFilterEvaluation_3d7_hb3_indels_5 <- genotypeFilterEvaluation("3d7_hb3", "indels", minMeanMAFtoConsiderContam=0.02, MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1), plotFilestemExtra="DP_QUALbyDP5", monomorphicSkipChromosomes  = "Pf3D7_13_v3")
  DP_QUALbyDP5_genotypeFilterEvaluation_7g8_gb4_indels_5 <- genotypeFilterEvaluation("7g8_gb4", "indels", minMeanMAFtoConsiderContam=0.05, parentalIDs = c("ERR027099", "ERR027100", "ERR029410", "ERR045626"), MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1), plotFilestemExtra="DP_QUALbyDP5")
  DP_QUALbyDP5_genotypeFilterEvaluation_hb3_dd2_indels_5 <- genotypeFilterEvaluation("hb3_dd2", "indels", minMeanMAFtoConsiderContam=0.05, parentalIDs = c("ERR012788", "ERR012840", "ERR022939"), MAFthresholds=NULL, DPthresholds=c(15, 10, 8, 6, 5, 4, 3, 2, 1), plotFilestemExtra="DP_QUALbyDP5")

#  genotypeFilterEvaluation2_3d7_hb3_snps <- genotypeFilterEvaluation("3d7_hb3", "snps", monomorphicSkipChromosomes  = "Pf3D7_13_v3")
#  genotypeFilterEvaluation2_7g8_gb4_snps <- genotypeFilterEvaluation("7g8_gb4", "snps")
#  genotypeFilterEvaluation2_hb3_dd2_snps <- genotypeFilterEvaluation("hb3_dd2", "snps")
#  genotypeFilterEvaluation2_3d7_hb3_indels <- genotypeFilterEvaluation("3d7_hb3", "indels", monomorphicSkipChromosomes  = "Pf3D7_13_v3", minMeanMAFtoConsiderContam=0.02)
#  genotypeFilterEvaluation2_7g8_gb4_indels <- genotypeFilterEvaluation("7g8_gb4", "indels", minMeanMAFtoConsiderContam=0.02)
#  genotypeFilterEvaluation2_hb3_dd2_indels <- genotypeFilterEvaluation("hb3_dd2", "indels", minMeanMAFtoConsiderContam=0.02)

genotypeFilterEvaluation <- function(
  cross                       = "3d7_hb3",
  variantType                 = "snps",
  analysisDirectory           = "/data/malariagen2/plasmodium/pf-crosses/data/3d7_v3/bwa_n0.01_k4_l32/genotypes_analysis_20120107/gatk",
#  filters                     = NULL,
  filters=list(
    "QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5),
    "QD14" = list(column="QD", operator="<", value=14),
    "BaseQRankSum45" = list(column="BaseQRankSum", operator="<", value=-45, filterOutNAs=FALSE),
    "numFilteredGenotypes10" = list(column="numFilteredGenotypes", operator=">", value=10)
  ),
#  filters=list(
##      "biallelic" = list(column="biallelic", operator="==", value=TRUE),
#    "meanMAF0.2" = list(column="meanMAF", operator=">", value=0.2, filterOutNAs=TRUE),
#    "QD" = list(column="QUALbyDP", operator=">", value=33),
#    "missingness1" = list(column="missingness", operator=">", value=1),
#    "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2),
#    "MQ0Fraction" = list(column="MQ0Fraction", operator=">", value=0.01),
#    "MateOtherChrom" = list(column="MateOtherChrom", operator=">", value=0.1),
#    "SoftClipped" = list(column="MateOtherChrom", operator=">", value=0.2),
#    "ProperPair" = list(column="ProperPair", operator="<", value=0.9),
#    "MQ" = list(column="MQ", operator="<", value=38),
#    "HRun" = list(column="HRun", operator=">", value=4, filterOutNAs=TRUE),
#    "FS" = list(column="FS", operator=">", value=600)
#  ),
#  regionsMask                 = NULL,
  regionsMask                 = varRegions_v3("/data/malariagen2/plasmodium/pf-crosses/meta/regions_v3_rdp.bed"),
  parentalIDs                 = NULL,
  setMonomorphicProgenyFilter = TRUE,
  monomorphicSkipChromosomes  = NULL,
  GQthresholds                = NULL,
  DPthresholds                = NULL,
#  GQthresholds                = c(99, 50, 5),
#  DPthresholds                = c(10, 5, 1),
  MAFthresholds               = c(0.1, 0, 0.02, 0.05, 0.2, 0.35, 0.5),
#  GQthresholds                = c(99, seq(95, 0, -5)),
#  DPthresholds                = seq(20, 1, -1),
#  MAFthresholds               = c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5),
  GQthresholdDefault          = 99,
  DPthresholdDefault          = 10,
  MAFthresholdDefault         = 0.1,
  plotFilestemExtra           = "MAF",
  shouldCreateQCFilteringPlots= TRUE,
  genotypeFiltersList         = c(
    lapply(
      MAFthresholds,
      function(MAFthreshold) {
        list(
          "LowGQ" = list(column="GQ", operator="<", value=GQthresholdDefault, filterOutNAs=TRUE),
          "LowDP" = list(column="DP", operator="<", value=DPthresholdDefault, filterOutNAs=TRUE),
          "HighMAF" = list(column="MAF", operator=">", value=MAFthreshold, filterOutNAs=TRUE)
        )
      }
    ),
    lapply(
      DPthresholds,
      function(DPthreshold) {
        list(
          "LowGQ" = list(column="GQ", operator="<", value=GQthresholdDefault, filterOutNAs=TRUE),
          "LowDP" = list(column="DP", operator="<", value=DPthreshold, filterOutNAs=TRUE),
          "HighMAF" = list(column="MAF", operator=">", value=MAFthresholdDefault, filterOutNAs=TRUE)
        )
      }
    ),
    lapply(
      GQthresholds,
      function(GQthreshold) {
        list(
          "LowGQ" = list(column="GQ", operator="<", value=GQthreshold, filterOutNAs=TRUE),
          "LowDP" = list(column="DP", operator="<", value=DPthresholdDefault, filterOutNAs=TRUE),
          "HighMAF" = list(column="MAF", operator=">", value=MAFthresholdDefault, filterOutNAs=TRUE)
        )
      }
    )
  ),
  maxNumFilteredGenotypes     = 2,
  minMeanMAFtoConsiderContam  = 0.01,
  sampleSets                  = list(
    "FinalSamples"                      = list(additionalInfoFilters=filters, sampleSetName="final"),
    "BestReplicate"                     = list(additionalInfoFilters=filters, sampleSetName="bestReplicate")
  ),
#  sampleSets                  = list(
#    "FinalSamples"                      = list(additionalInfoFilters=filters, sampleSetName="final"),
#    "FinalSamplesGenotypeOnly"          = list(additionalInfoFilters=NULL,    sampleSetName="final"),
#    "BestReplicate"                     = list(additionalInfoFilters=filters, sampleSetName="bestReplicate"),
#    "BestReplicateGenotypeOnly"         = list(additionalInfoFilters=NULL,    sampleSetName="bestReplicate"),
#    "UncontaminatedSamples"             = list(additionalInfoFilters=filters, sampleSetName="uncontaminated"),
#    "UncontaminatedSamplesGenotypeOnly" = list(additionalInfoFilters=NULL,    sampleSetName="uncontaminated")
#  ),
  shouldReturnVcfOnly         = FALSE,
  shouldUseExistingRda        = FALSE
)
{
  initialSampleQCresultsFilename <- file.path(analysisDirectory, cross, variantType, paste(cross, ".initialSampleQCresults.rda", sep=""))
  if(file.exists(initialSampleQCresultsFilename)) {
    load(initialSampleQCresultsFilename)
  } else {
    load(file.path(analysisDirectory, cross, variantType, paste(cross, ".vcfInitialFiltered.rda", sep="")))
    vcfSegregating <- filterVcf(vcfInitialFiltered, keepPASSvariantsOnly=TRUE)
    initialSampleQCresults <- sampleQC(
      vcfSegregating,
      discordanceThreshold=1000,
      shouldCreatePlots=FALSE,
      sampleIDcolumn="ena_run_accession",
      sampleIDmappingsColumn="ena_run_accession"
    )
    save(initialSampleQCresults, file=initialSampleQCresultsFilename)
  }
  vcfAnnotatedFinalSamplesFilename <- file.path(analysisDirectory, cross, variantType, paste(cross, ".vcfAnnotatedFinalSamples.rda", sep=""))
  vcfAnnotatedBestReplicateSamplesFilename <- file.path(analysisDirectory, cross, variantType, paste(cross, ".vcfAnnotatedBestReplicateSamples.rda", sep=""))
  vcfAnnotatedUncontaminatedSamplesFilename <- file.path(analysisDirectory, cross, variantType, paste(cross, ".vcfAnnotatedUncontaminatedSamples.rda", sep=""))
  vcfList <- list()
  sampleSetNames <- sapply(sampleSets, function(x) x[["sampleSetName"]])
  if(shouldUseExistingRda && file.exists(vcfAnnotatedFinalSamplesFilename) && file.exists(vcfAnnotatedBestReplicateSamplesFilename) && file.exists(vcfAnnotatedUncontaminatedSamplesFilename)) {
    load(vcfAnnotatedFinalSamplesFilename)
    load(vcfAnnotatedBestReplicateSamplesFilename)
    load(vcfAnnotatedUncontaminatedSamplesFilename)
    vcfList[["final"]] <- vcfAnnotatedFinalSamples
    vcfList[["bestReplicate"]] <- vcfAnnotatedBestReplicateSamples
    vcfList[["uncontaminated"]] <- vcfAnnotatedUncontaminatedSamples
  } else {
    if(is.null(regionsMask)) {
      load(file.path(analysisDirectory, cross, variantType, paste(cross, ".vcfVariant.rda", sep="")))
    } else {
      load(file.path(analysisDirectory, cross, variantType, paste(cross, ".vcfCoreFinalSamples.rda", sep="")))
      vcfVariant <- vcfCoreFinalSamples
      filt(vcfVariant) <- "."
      rm(vcfCoreFinalSamples)
      gc()
    }
    if("final" %in% sampleSetNames) {
      finalSamples <- setdiff(dimnames(vcfVariant)[[2]], initialSampleQCresults[["qcFailedSamples"]])
      vcfList[["final"]] <-  annotateVcf(vcfVariant[, finalSamples])
#      save(vcfList[["final"]], file=vcfAnnotatedFinalSamplesFilename)
    }
    if("bestReplicate" %in% sampleSetNames) {
  #    bestReplicateSamples <- setdiff(dimnames(vcfVariant)[[2]], initialSampleQCresults[["uniqueSamples"]])
      bestReplicateSamples <- setdiff(initialSampleQCresults[["uniqueSamples"]], initialSampleQCresults[["qcFailedSamples"]])
      vcfList[["bestReplicate"]] <-  annotateVcf(vcfVariant[, bestReplicateSamples])
#      save(vcfList[["bestReplicate"]], file=vcfAnnotatedBestReplicateSamplesFilename)
    }
    if("uncontaminated" %in% sampleSetNames) {
      if(!exists("vcfInitialFiltered")) {
        load(file.path(analysisDirectory, cross, variantType, paste(cross, ".vcfInitialFiltered.rda", sep="")))
      }
      vcfInitialFilteredPASS <- vcfInitialFiltered[filt(vcfInitialFiltered)=="PASS"]
      RefReads <- matrix(
        sapply(geno(vcfInitialFilteredPASS)[["AD"]], function(x) x[1]),
        ncol=dim(geno(vcfInitialFilteredPASS)[["AD"]])[2],
        dimnames=dimnames(geno(vcfInitialFilteredPASS)[["AD"]])
      )
      FirstAltReads <- matrix(
        sapply(geno(vcfInitialFilteredPASS)[["AD"]], function(x) x[2]),
        ncol=dim(geno(vcfInitialFilteredPASS)[["AD"]])[2],
        dimnames=dimnames(geno(vcfInitialFilteredPASS)[["AD"]])
      )
      MAF <- pmin(RefReads, FirstAltReads)/(RefReads+FirstAltReads)
      meanMAFperSample <- colMeans(MAF, na.rm = TRUE)
      uncontaminatedSamples <- intersect(names(which(meanMAFperSample < minMeanMAFtoConsiderContam)), dimnames(vcfVariant)[[2]])
      vcfList[["uncontaminated"]] <-  annotateVcf(vcfVariant[, uncontaminatedSamples])
#      save(vcfList[["uncontaminated"]], file=vcfAnnotatedUncontaminatedSamplesFilename)
    }
  }
  
#  browser()
#  
#  filterResults <- evaluateGenotypeFilters(
#    vcfAnnotatedFinalSamples,
##    vcfCoreFinalSamples,
#    plotFilestem                = file.path(analysisDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", sep=".")),
#    additionalInfoFilters       = filters,
#    regionsMask                 = regionsMask,
#    genotypeFiltersList         = genotypeFiltersList,
#    setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
#    monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
#    maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
#    sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
#    shouldReturnVcfOnly         = shouldReturnVcfOnly
#  )
  filterResultsList <- sapply(
    names(sampleSets),
    function(sampleSet) {
      if(is.null(parentalIDs)) {
        filterResults <- evaluateGenotypeFilters(
          vcfList[[sampleSets[[sampleSet]][["sampleSetName"]]]],
#          vcfAnnotatedFinalSamples,
      #    vcfCoreFinalSamples,
          plotFilestem                = file.path(analysisDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", plotFilestemExtra, sampleSets[[sampleSet]][["sampleSetName"]], sep=".")),
          additionalInfoFilters       = sampleSets[[sampleSet]][["additionalInfoFilters"]],
          regionsMask                 = regionsMask,
          shouldCreateQCFilteringPlots=shouldCreateQCFilteringPlots,
          genotypeFiltersList         = genotypeFiltersList,
          setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
          monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
          maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
          sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
          shouldReturnVcfOnly         = shouldReturnVcfOnly
        )
      } else {
        filterResults <- evaluateGenotypeFilters(
          vcfList[[sampleSets[[sampleSet]][["sampleSetName"]]]],
#          vcfAnnotatedFinalSamples,
      #    vcfCoreFinalSamples,
          plotFilestem                = file.path(analysisDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", plotFilestemExtra, sampleSets[[sampleSet]][["sampleSetName"]], sep=".")),
          additionalInfoFilters       = sampleSets[[sampleSet]][["additionalInfoFilters"]],
          regionsMask                 = regionsMask,
          shouldCreateQCFilteringPlots=shouldCreateQCFilteringPlots,
          genotypeFiltersList         = genotypeFiltersList,
          setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
          monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
          maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
          sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
          parentalIDs                 = parentalIDs,
          shouldReturnVcfOnly         = shouldReturnVcfOnly
        )
      }
      names(filterResults) <- paste(sampleSet, names(filterResults))
      filterResults
    },
    simplify=FALSE,
    USE.NAMES=TRUE
  )
#  if("FinalSamples" %in% sampleSets) {
#    if(is.null(parentalIDs)) {
#      filterResultsList[["FinalSamples"]] <- evaluateGenotypeFilters(
#        vcfAnnotatedFinalSamples,
#    #    vcfCoreFinalSamples,
#        plotFilestem                = file.path(analysisDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", plotFilestemExtra, "final", sep=".")),
#        additionalInfoFilters       = filters,
#        regionsMask                 = regionsMask,
#        genotypeFiltersList         = genotypeFiltersList,
#        setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
#        monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
#        maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
#        sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
#        shouldReturnVcfOnly         = shouldReturnVcfOnly
#      )
#    } else {
#      filterResultsList[["FinalSamples"]] <- evaluateGenotypeFilters(
#        vcfAnnotatedFinalSamples,
#    #    vcfCoreFinalSamples,
#        plotFilestem                = file.path(analysisDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", plotFilestemExtra, "final", sep=".")),
#        additionalInfoFilters       = filters,
#        regionsMask                 = regionsMask,
#        genotypeFiltersList         = genotypeFiltersList,
#        setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
#        monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
#        maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
#        sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
#        parentalIDs                 = parentalIDs,
#        shouldReturnVcfOnly         = shouldReturnVcfOnly
#      )
#    }
#    names(filterResultsList[["FinalSamples"]]) <- paste("FinalSamples", names(filterResultsList[["FinalSamples"]]))
#  }
#  if(is.null(parentalIDs)) {
#    filterResultsFinalSamplesGenotypeOnly <- evaluateGenotypeFilters(
#      vcfAnnotatedFinalSamples,
#  #    vcfCoreFinalSamples,
#      plotFilestem                = file.path(analysisDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", plotFilestemExtra, "final", sep=".")),
#      additionalInfoFilters       = NULL,
#      regionsMask                 = regionsMask,
#      genotypeFiltersList         = genotypeFiltersList,
#      setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
#      monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
#      maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
#      sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
#      shouldReturnVcfOnly         = shouldReturnVcfOnly
#    )
#  } else {
#    filterResultsFinalSamplesGenotypeOnly <- evaluateGenotypeFilters(
#      vcfAnnotatedFinalSamples,
#  #    vcfCoreFinalSamples,
#      plotFilestem                = file.path(analysisDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", plotFilestemExtra, "final", sep=".")),
#      additionalInfoFilters       = NULL,
#      regionsMask                 = regionsMask,
#      genotypeFiltersList         = genotypeFiltersList,
#      setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
#      monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
#      maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
#      sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
#      parentalIDs                 = parentalIDs,
#      shouldReturnVcfOnly         = shouldReturnVcfOnly
#    )
#  }
#  names(filterResultsFinalSamplesGenotypeOnly) <- paste("FinalSamplesGenotypeOnly", names(filterResultsFinalSamplesGenotypeOnly))
#  if(is.null(parentalIDs)) {
#    filterResultsBestReplicateSamples <- evaluateGenotypeFilters(
#      vcfAnnotatedBestReplicateSamples,
#  #    vcfCoreFinalSamples,
#      plotFilestem                = file.path(analysisDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", plotFilestemExtra, "bestReplicate", sep=".")),
#      additionalInfoFilters       = filters,
#      regionsMask                 = regionsMask,
#      genotypeFiltersList         = genotypeFiltersList,
#      setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
#      monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
#      maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
#      sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
#      shouldReturnVcfOnly         = shouldReturnVcfOnly
#    )
#  } else {
#    filterResultsBestReplicateSamples <- evaluateGenotypeFilters(
#      vcfAnnotatedBestReplicateSamples,
#  #    vcfCoreFinalSamples,
#      plotFilestem                = file.path(analysisDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", plotFilestemExtra, "bestReplicate", sep=".")),
#      additionalInfoFilters       = filters,
#      regionsMask                 = regionsMask,
#      genotypeFiltersList         = genotypeFiltersList,
#      setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
#      monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
#      maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
#      sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
#      parentalIDs                 = parentalIDs,
#      shouldReturnVcfOnly         = shouldReturnVcfOnly
#    )
#  }
#  names(filterResultsBestReplicateSamples) <- paste("BestReplicate", names(filterResultsBestReplicateSamples))
#  if(is.null(parentalIDs)) {
#    filterResultsBestReplicateSamplesGenotypeOnly <- evaluateGenotypeFilters(
#      vcfAnnotatedBestReplicateSamples,
#  #    vcfCoreFinalSamples,
#      plotFilestem                = file.path(analysisDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", plotFilestemExtra, "bestReplicate", sep=".")),
#      additionalInfoFilters       = NULL,
#      regionsMask                 = regionsMask,
#      genotypeFiltersList         = genotypeFiltersList,
#      setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
#      monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
#      maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
#      sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
#      shouldReturnVcfOnly         = shouldReturnVcfOnly
#    )
#  } else {
#    filterResultsBestReplicateSamplesGenotypeOnly <- evaluateGenotypeFilters(
#      vcfAnnotatedBestReplicateSamples,
#  #    vcfCoreFinalSamples,
#      plotFilestem                = file.path(analysisDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", plotFilestemExtra, "bestReplicate", sep=".")),
#      additionalInfoFilters       = NULL,
#      regionsMask                 = regionsMask,
#      genotypeFiltersList         = genotypeFiltersList,
#      setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
#      monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
#      maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
#      sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
#      parentalIDs                 = parentalIDs,
#      shouldReturnVcfOnly         = shouldReturnVcfOnly
#    )
#  }
#  names(filterResultsBestReplicateSamplesGenotypeOnly) <- paste("BestReplicateGenotypeOnly", names(filterResultsBestReplicateSamplesGenotypeOnly))
#  if(is.null(parentalIDs)) {
#    filterResultsUncontaminatedSamples <- evaluateGenotypeFilters(
#      vcfAnnotatedUncontaminatedSamples,
#  #    vcfCoreFinalSamples,
#      plotFilestem                = file.path(analysisDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", plotFilestemExtra, "uncontaminated", sep=".")),
#      additionalInfoFilters       = filters,
#      regionsMask                 = regionsMask,
#      genotypeFiltersList         = genotypeFiltersList,
#      setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
#      monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
#      maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
#      sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
#      shouldReturnVcfOnly         = shouldReturnVcfOnly
#    )
#  } else {
#    filterResultsUncontaminatedSamples <- evaluateGenotypeFilters(
#      vcfAnnotatedUncontaminatedSamples,
#  #    vcfCoreFinalSamples,
#      plotFilestem                = file.path(analysisDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", plotFilestemExtra, "uncontaminated", sep=".")),
#      additionalInfoFilters       = filters,
#      regionsMask                 = regionsMask,
#      genotypeFiltersList         = genotypeFiltersList,
#      setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
#      monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
#      maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
#      sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
#      parentalIDs                 = parentalIDs,
#      shouldReturnVcfOnly         = shouldReturnVcfOnly
#    )
#  }
#  names(filterResultsUncontaminatedSamples) <- paste("UncontaminatedSamples", names(filterResultsUncontaminatedSamples))
#  if(is.null(parentalIDs)) {
#    filterResultsUncontaminatedSamplesGenotypeOnly <- evaluateGenotypeFilters(
#      vcfAnnotatedUncontaminatedSamples,
#  #    vcfCoreFinalSamples,
#      plotFilestem                = file.path(analysisDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", plotFilestemExtra, "uncontaminated", sep=".")),
#      additionalInfoFilters       = NULL,
#      regionsMask                 = regionsMask,
#      genotypeFiltersList         = genotypeFiltersList,
#      setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
#      monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
#      maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
#      sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
#      shouldReturnVcfOnly         = shouldReturnVcfOnly
#    )
#  } else {
#    filterResultsUncontaminatedSamplesGenotypeOnly <- evaluateGenotypeFilters(
#      vcfAnnotatedUncontaminatedSamples,
#  #    vcfCoreFinalSamples,
#      plotFilestem                = file.path(analysisDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", plotFilestemExtra, "uncontaminated", sep=".")),
#      additionalInfoFilters       = NULL,
#      regionsMask                 = regionsMask,
#      genotypeFiltersList         = genotypeFiltersList,
#      setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
#      monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
#      maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
#      sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
#      parentalIDs                 = parentalIDs,
#      shouldReturnVcfOnly         = shouldReturnVcfOnly
#    )
#  }
#  names(filterResultsUncontaminatedSamplesGenotypeOnly) <- paste("UncontaminatedSamplesGenotypeOnly", names(filterResultsUncontaminatedSamplesGenotypeOnly))
  filterResults <- do.call(cbind, filterResultsList)
#  filterResults <- cbind(
#    filterResultsFinalSamples, filterResultsBestReplicateSamples, filterResultsUncontaminatedSamples,
#    filterResultsFinalSamplesGenotypeOnly, filterResultsBestReplicateSamplesGenotypeOnly, filterResultsUncontaminatedSamplesGenotypeOnly
#  )
  return(filterResults)
}

# allFilterEvaluationSeries.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################

#  allFilterEvaluationSeriesResults <- allFilterEvaluationSeries()
#  allFilterEvaluationSeriesResults_7g8_gb4_noSD <- allFilterEvaluationSeries(
#    callsets = list(
#      snps_7g8_gb4 = c(cross = "7g8_gb4", variantType = "snps")
#    ),
#    filters = list(
#      "hp5prox1" = list(column="homopolymer5Proximity", operator="%in%", value=1),
#      "hp15prox0to20" = list(column="homopolymer15Proximity", operator="%in%", value=0:20),
#      "MQ0" = list(column="MQ0", operator=">", value=0),
#      "TRF" = list(column="RepeatPeriod1", operator=">", value=0),
#      "ReadPosMinus2" = list(column="ReadPosRankSum", operator="<", value=-2, filterOutNAs=TRUE),
#      "ProperPair0.95" = list(column="ProperPair", operator="<", value=0.95),
#      "SoftClip0.1" = list(column="SoftClipped", operator=">", value=0.1),
#      "UQ35" = list(column="UQ", operator=">", value=35),
#      "QUAL12000" = list(column="QUAL", operator="<", value=12000)
#    )
#  )
#  allFilterEvaluationSeriesResultsQDandSB <- allFilterEvaluationSeries(
#    filters=list(
#      "depthSD0.5" = list(column="scaledDepthSD", operator=">", value=0.5),
#      "QD36" = list(column="QD", operator="<=", value=36),
#      "SBminus6000" = list(column="SB", operator=">=", value=-6000)
#    )
#  )
#  allFilterEvaluationSeriesResultsQDandSB_3d7_hb3_SD1 <- allFilterEvaluationSeries(
#    callsets = list(
#      snps_3d7_hb3 = c(cross = "3d7_hb3", variantType = "snps")
#    ),
#    filters=list(
#      "depthSD1" = list(column="scaledDepthSD", operator=">", value=1),
#      "QD36" = list(column="QD", operator="<=", value=36),
#      "SBminus6000" = list(column="SB", operator=">=", value=-6000)
#    )
#  )
#  filterEvaluations_snps_QD30 <- allFilterEvaluationSeries(
#    callsets = list(
#      snps_3d7_hb3 = c(cross = "3d7_hb3", variantType = "snps"),
#      snps_7g8_gb4 = c(cross = "7g8_gb4", variantType = "snps"),
#      snps_hb3_dd2 = c(cross = "hb3_dd2", variantType = "snps")
#    ),
#    filters=list(
#      "unfiltered" = list(column="QUAL", operator="<", value=0),
#      "QD30" = list(column="QD", operator="<=", value=30)
#    )
#  )
#  filterEvaluations_snps_QD30_UQ999 <- allFilterEvaluationSeries(
#    callsets = list(
#      snps_3d7_hb3 = c(cross = "3d7_hb3", variantType = "snps"),
#      snps_7g8_gb4 = c(cross = "7g8_gb4", variantType = "snps"),
#      snps_hb3_dd2 = c(cross = "hb3_dd2", variantType = "snps")
#    ),
#    filters=list(
#      "unfiltered" = list(column="QUAL", operator="<", value=0),
#      "QD30" = list(column="QD", operator="<=", value=30),
#      "UQ999" = list(column="UQ", operator=">", value=998)
#    )
#  )
#  filterEvaluations_snps_meanMAF0.01_QD30 <- allFilterEvaluationSeries(
#    callsets = list(
#      snps_3d7_hb3 = c(cross = "3d7_hb3", variantType = "snps"),
#      snps_7g8_gb4 = c(cross = "7g8_gb4", variantType = "snps"),
#      snps_hb3_dd2 = c(cross = "hb3_dd2", variantType = "snps")
#    ),
#    filters=list(
#      "unfiltered" = list(column="QUAL", operator="<", value=0),
#      "meanMAF0.01" = list(column="meanMAF", operator=">", value=0.01),
#      "QD30" = list(column="QD", operator="<=", value=30)
#    )
#  )
#  filterEvaluations_snps_meanMAF0.01 <- allFilterEvaluationSeries(
#    callsets = list(
#      snps_3d7_hb3 = c(cross = "3d7_hb3", variantType = "snps"),
#      snps_7g8_gb4 = c(cross = "7g8_gb4", variantType = "snps"),
#      snps_hb3_dd2 = c(cross = "hb3_dd2", variantType = "snps")
#    ),
#    filters=list(
#      "meanMAF0.01" = list(column="meanMAF", operator=">", value=0.01)
#    )
#  )
#  filterEvaluations_snps_QUALperSample300_meanMAF0.01 <- allFilterEvaluationSeries(
#    callsets = list(
#      snps_3d7_hb3 = c(cross = "3d7_hb3", variantType = "snps"),
#      snps_7g8_gb4 = c(cross = "7g8_gb4", variantType = "snps"),
#      snps_hb3_dd2 = c(cross = "hb3_dd2", variantType = "snps")
#    ),
#    filters=list(
#      "QUALperSample300" = list(column="QUALperSample", operator="<", value=300),
#      "meanMAF0.01" = list(column="meanMAF", operator=">", value=0.01)
#    )
#  )
#  filterEvaluations_snps_withoutHeterozygosity <- allFilterEvaluationSeries(
#    callsets = list(
#      snps_3d7_hb3 = c(cross = "3d7_hb3", variantType = "snps"),
#      snps_7g8_gb4 = c(cross = "7g8_gb4", variantType = "snps"),
#      snps_hb3_dd2 = c(cross = "hb3_dd2", variantType = "snps")
#    ),
#    filters = list(
#      "hp5prox1" = list(column="homopolymer5Proximity", operator="%in%", value=1),
#      "hp15prox0to20" = list(column="homopolymer15Proximity", operator="%in%", value=0:20),
#      "MQ0" = list(column="MQ0", operator=">", value=0),
#      "TRF" = list(column="RepeatPeriod1", operator=">", value=0),
#      "ReadPosMinus2" = list(column="ReadPosRankSum", operator="<", value=-2, filterOutNAs=TRUE),
#      "ProperPair0.95" = list(column="ProperPair", operator="<", value=0.95),
#      "SoftClip0.1" = list(column="SoftClipped", operator=">", value=0.1),
#      "UQ35" = list(column="UQ", operator=">", value=35),
#      "depthSD1" = list(column="scaledDepthSD", operator=">", value=1),
#      "QUAL12000" = list(column="QUAL", operator="<", value=12000)
#    )
#  )


allFilterEvaluationSeries <- function(
  callsets = list(
    snps_3d7_hb3 = c(cross = "3d7_hb3", variantType = "snps"),
    snps_7g8_gb4 = c(cross = "7g8_gb4", variantType = "snps"),
    snps_hb3_dd2 = c(cross = "hb3_dd2", variantType = "snps"),
    indels_3d7_hb3 = c(cross = "3d7_hb3", variantType = "indels"),
    indels_7g8_gb4 = c(cross = "7g8_gb4", variantType = "indels"),
    indels_hb3_dd2 = c(cross = "hb3_dd2", variantType = "indels")
  ),
  analysisDirectory           = "/data/malariagen2/plasmodium/pf-crosses/data/3d7_v3/bwa_n0.01_k4_l32/genotypes_analysis_20120107/gatk",
  filters = list(
    "hp5prox1" = list(column="homopolymer5Proximity", operator="%in%", value=1),
    "hp15prox0to20" = list(column="homopolymer15Proximity", operator="%in%", value=0:20),
    "MQ0" = list(column="MQ0", operator=">", value=0),
    "TRF" = list(column="RepeatPeriod1", operator=">", value=0),
    "ReadPosMinus2" = list(column="ReadPosRankSum", operator="<", value=-2, filterOutNAs=TRUE),
    "depthSD1" = list(column="scaledDepthSD", operator=">", value=1),
    "SoftClip0.1" = list(column="SoftClipped", operator=">", value=0.1),
    "UQ35" = list(column="UQ", operator=">", value=35),
    "QUAL12000" = list(column="QUAL", operator="<", value=12000)
  )
) {
  resultsList <- sapply(
    callsets,
    function(callset) {
      filterEvaluationSeries(
        cross                       = callset["cross"],
        variantType                 = callset["variantType"],
        analysisDirectory           = analysisDirectory,
        filters = filters
      )
    },
    simplify=FALSE,
    USE.NAMES=TRUE
  )
  return(resultsList)
}

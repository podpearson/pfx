# allFilterEvaluations.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################

#  filterEvaluations_snps_QUALperSample300_maxMAF0.4_QUALbyDP40_missingness1 <- allFilterEvaluations(
#    callsets = list(
##      snps_3d7_hb3 = c(cross = "3d7_hb3", variantType = "snps"),
#      snps_7g8_gb4 = c(cross = "7g8_gb4", variantType = "snps"),
#      snps_hb3_dd2 = c(cross = "hb3_dd2", variantType = "snps")
#    ),
#    filters=list(
#      "QUALperSample300" = list(column="QUALperSample", operator="<", value=300),
#      "maxMAF0.4" = list(column="maxMAF", operator=">", value=0.4),
#      "QUALbyDP40" = list(column="QUALbyDP", operator=">", value=40),
#      "missingness1" = list(column="missingness", operator=">", value=1)
#    )
#  )
#  filterEvaluations_snps_QUALperSample300_maxMAF0.2_QUALbyDP40_missingness1 <- allFilterEvaluations(
#    callsets = list(
#      snps_3d7_hb3 = c(cross = "3d7_hb3", variantType = "snps"),
#      snps_7g8_gb4 = c(cross = "7g8_gb4", variantType = "snps"),
#      snps_hb3_dd2 = c(cross = "hb3_dd2", variantType = "snps")
#    ),
#    filters=list(
#      "QUALperSample300" = list(column="QUALperSample", operator="<", value=300),
#      "maxMAF0.2" = list(column="maxMAF", operator=">", value=0.2),
#      "QUALbyDP40" = list(column="QUALbyDP", operator=">", value=40),
#      "missingness1" = list(column="missingness", operator=">", value=1)
#    )
#  )
#  filterEvaluations_snps_QUALperSample300_maxMAF0.2_QUALbyDP29_missingness1 <- allFilterEvaluations(
#    callsets = list(
#      snps_3d7_hb3 = c(cross = "3d7_hb3", variantType = "snps"),
#      snps_7g8_gb4 = c(cross = "7g8_gb4", variantType = "snps"),
#      snps_hb3_dd2 = c(cross = "hb3_dd2", variantType = "snps")
#    ),
#    filters=list(
#      "QUALperSample300" = list(column="QUALperSample", operator="<", value=300),
#      "maxMAF0.2" = list(column="maxMAF", operator=">", value=0.2),
#      "QUALbyDP29" = list(column="QUALbyDP", operator=">", value=29),
#      "missingness1" = list(column="missingness", operator=">", value=1)
#    )
#  )
#  filterEvaluations_snps_QUALperSample300_maxMAF0.2_QUALbyDP29_missingness1_numFilteredGenotypes2 <- allFilterEvaluations(
#    callsets = list(
#      snps_3d7_hb3 = c(cross = "3d7_hb3", variantType = "snps"),
#      snps_7g8_gb4 = c(cross = "7g8_gb4", variantType = "snps"),
#      snps_hb3_dd2 = c(cross = "hb3_dd2", variantType = "snps")
#    ),
#    genotypeFilters = list(
#      "LowGQ" = list(column="GQ", operator="<", value=99, filterOutNAs=TRUE),
#      "LowDP" = list(column="DP", operator="<", value=10, filterOutNAs=TRUE),
#      "HighMAF" = list(column="MAF", operator=">", value=0.1, filterOutNAs=TRUE)
#    ),
#    filters=list(
#      "QUALperSample300" = list(column="QUALperSample", operator="<", value=300),
#      "maxMAF0.2" = list(column="maxMAF", operator=">", value=0.2),
#      "QUALbyDP29" = list(column="QUALbyDP", operator=">", value=29),
#      "missingness1" = list(column="missingness", operator=">", value=1),
#      "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2)
#    )
#  )

#  filterEvaluations_snps_QUALperSample300_maxMAF0.2_missingness1 <- allFilterEvaluations(
#    callsets = list(
#      snps_3d7_hb3 = c(cross = "3d7_hb3", variantType = "snps"),
#      snps_7g8_gb4 = c(cross = "7g8_gb4", variantType = "snps"),
#      snps_hb3_dd2 = c(cross = "hb3_dd2", variantType = "snps")
#    ),
#    filters=list(
#      "QUALperSample300" = list(column="QUALperSample", operator="<", value=300),
#      "maxMAF0.2" = list(column="maxMAF", operator=">", value=0.2),
#      "missingness1" = list(column="missingness", operator=">", value=1)
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


allFilterEvaluations <- function(
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
  ),
  shouldReturnVcfOnly         = FALSE
) {
  resultsList <- sapply(
    callsets,
    function(callset) {
      filterEvaluation(
        cross                       = callset["cross"],
        variantType                 = callset["variantType"],
        analysisDirectory           = analysisDirectory,
        filters                     = filters,
        genotypeFilters             = genotypeFilters,
        shouldReturnVcfOnly         = shouldReturnVcfOnly
      )
    },
    simplify=FALSE,
    USE.NAMES=TRUE
  )
  return(resultsList)
}

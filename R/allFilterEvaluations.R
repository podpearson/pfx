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
#  filterEvaluations_snps_QUALperSample300_maxMAF0.2_QUALbyDP29_missingness1_numFilteredGenotypes2_DP5 <- allFilterEvaluations(
#    callsets = list(
#      snps_3d7_hb3 = c(cross = "3d7_hb3", variantType = "snps"),
#      snps_7g8_gb4 = c(cross = "7g8_gb4", variantType = "snps"),
#      snps_hb3_dd2 = c(cross = "hb3_dd2", variantType = "snps")
#    ),
#    genotypeFilters = list(
#      "LowGQ" = list(column="GQ", operator="<", value=99, filterOutNAs=TRUE),
#      "LowDP" = list(column="DP", operator="<", value=5, filterOutNAs=TRUE),
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
#  filterEvaluations_alistair_20130110_2 <- allFilterEvaluations(
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
##      "biallelic" = list(column="biallelic", operator="==", value=TRUE),
#      "meanMAF0.2" = list(column="meanMAF", operator=">", value=0.2),
#      "QD" = list(column="QUALbyDP", operator=">", value=33),
#      "missingness1" = list(column="missingness", operator=">", value=1),
#      "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2),
#      "MQ0Fraction" = list(column="MQ0Fraction", operator=">", value=0.01),
#      "MateOtherChrom" = list(column="MateOtherChrom", operator=">", value=0.1),
#      "SoftClipped" = list(column="MateOtherChrom", operator=">", value=0.2),
#      "ProperPair" = list(column="ProperPair", operator="<", value=0.9),
#      "MQ" = list(column="MQ", operator="<", value=38),
#      "HRun" = list(column="HRun", operator=">", value=4, filterOutNAs=TRUE),
#      "FS" = list(column="FS", operator=">", value=600)
#    ),
#    regionsMask                 = NULL,
#  )
#  filterEvaluations_snps_GT99_10_0.1 <- allFilterEvaluations(
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
#    ),
#    regionsMask                 = NULL,
#  )
#  filterEvaluations_snps_GT99_10_0.1_max3 <- allFilterEvaluations(
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
#    filters=list(),
#    maxNumFilteredGenotypes=3,
#    regionsMask                 = NULL,
#  )
#  filterEvaluations_snps_GT99_10_0.1_max50 <- allFilterEvaluations(
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
#    filters=list(),
#    maxNumFilteredGenotypes=50,
#    regionsMask                 = NULL,
#  )
#  filterEvaluations_snps_GT99_5_0.1_max2 <- allFilterEvaluations(
#    callsets = list(
#      snps_3d7_hb3 = c(cross = "3d7_hb3", variantType = "snps"),
#      snps_7g8_gb4 = c(cross = "7g8_gb4", variantType = "snps"),
#      snps_hb3_dd2 = c(cross = "hb3_dd2", variantType = "snps")
#    ),
#    genotypeFilters = list(
#      "LowGQ" = list(column="GQ", operator="<", value=99, filterOutNAs=TRUE),
#      "LowDP" = list(column="DP", operator="<", value=5, filterOutNAs=TRUE),
#      "HighMAF" = list(column="MAF", operator=">", value=0.1, filterOutNAs=TRUE)
#    ),
#    filters=list(),
#    maxNumFilteredGenotypes=2,
#    regionsMask                 = NULL,
#  )
#  filterEvaluations_snps_GT99_10_0.05_max2 <- allFilterEvaluations(
#    callsets = list(
#      snps_3d7_hb3 = c(cross = "3d7_hb3", variantType = "snps"),
#      snps_7g8_gb4 = c(cross = "7g8_gb4", variantType = "snps"),
#      snps_hb3_dd2 = c(cross = "hb3_dd2", variantType = "snps")
#    ),
#    genotypeFilters = list(
#      "LowGQ" = list(column="GQ", operator="<", value=99, filterOutNAs=TRUE),
#      "LowDP" = list(column="DP", operator="<", value=10, filterOutNAs=TRUE),
#      "HighMAF" = list(column="MAF", operator=">", value=0.05, filterOutNAs=TRUE)
#    ),
#    filters=list(),
#    maxNumFilteredGenotypes=2,
#    regionsMask                 = NULL,
#  )
#  filterEvaluations_snps_GT50_10_0.1_max2 <- allFilterEvaluations(
#    callsets = list(
#      snps_3d7_hb3 = c(cross = "3d7_hb3", variantType = "snps"),
#      snps_7g8_gb4 = c(cross = "7g8_gb4", variantType = "snps"),
#      snps_hb3_dd2 = c(cross = "hb3_dd2", variantType = "snps")
#    ),
#    genotypeFilters = list(
#      "LowGQ" = list(column="GQ", operator="<", value=50, filterOutNAs=TRUE),
#      "LowDP" = list(column="DP", operator="<", value=10, filterOutNAs=TRUE),
#      "HighMAF" = list(column="MAF", operator=">", value=0.1, filterOutNAs=TRUE)
#    ),
#    filters=list(),
#    maxNumFilteredGenotypes=2,
#    regionsMask                 = NULL,
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
  regionsMask                 = varRegions_v3(),
  genotypeFilters             = list(
    "LowGQ" = list(column="GQ", operator="<", value=99, filterOutNAs=TRUE),
    "LowDP" = list(column="DP", operator="<", value=10, filterOutNAs=TRUE),
    "HighMAF" = list(column="MAF", operator=">", value=0.1, filterOutNAs=TRUE)
  ),
  maxNumFilteredGenotypes     = 2,
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
        maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
        shouldReturnVcfOnly         = shouldReturnVcfOnly
      )
    },
    simplify=FALSE,
    USE.NAMES=TRUE
  )
  return(resultsList)
}

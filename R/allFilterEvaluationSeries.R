# allFilterEvaluationSeries.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


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
    "depthSD0.5" = list(column="scaledDepthSD", operator=">", value=0.5),
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

# compare7g8_gb4WithJiang2.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


compare7g8_gb4WithJiang2 <- function(
  malariagenVcfList             = list(
    Intersection = loadCallsSubset(
      vcfFilename                 = "data/release/1.0.combined.RC1/hb3_dd2.combined.vcf.gz",
      subsetVcfFilename           = "analysis/release/1.0.combined.RC1/hb3_dd2.combined.Intersection.vcf",
      subsetRdaFilename           = "analysis/release/1.0.combined.RC1/hb3_dd2.combined.Intersection.vcf.rda",
      subsetGrep                  = "set=Intersection"
    ),
    GATK = loadCallsSubset(
      vcfFilename                 = "data/release/1.0.GATK.RC6/hb3_dd2.gatk.both.final.vcf.gz",
      subsetVcfFilename           = "analysis/release/1.0.GATK.RC6/hb3_dd2.gatk.both.final.PASS.vcf",
      subsetRdaFilename           = "analysis/release/1.0.GATK.RC6/hb3_dd2.gatk.both.final.PASS.vcf.rda",
      subsetGrep                  = "PASS"
    ),
    Cortex = loadCallsSubset(
      vcfFilename                 = "data/release/1.0.cortex.RC1/hb3_dd2.cortex.final.vcf.gz",
      subsetVcfFilename           = "analysis/release/1.0.cortex.RC1/hb3_dd2.cortex.final.PASS.vcf",
      subsetRdaFilename           = "analysis/release/1.0.cortex.RC1/hb3_dd2.cortex.final.PASS.vcf.rda",
      subsetGrep                  = "PASS"
    )
  ),
  jiangVcf                    = loadJiangGenotypesAsVcf2(),
  discordanceThreshold        = 100,
  discordanceProportionThreshold = 0.15,
  plotFilestem                = "analysis/release/1.0.combined.RC1/jiang",
  IDparent1                   = "7G8_NIH/PG0083-C/ERR027099",
  IDparent2                   = "GB4_NIH/PG0084-C/ERR027100",
  GTsToIntMapping             = c("7"=1, "G"=2, "."=0),
  expectedMatches = list(
    comparisonVsSubject = rbind(
      c("XG10", "XG10.PG0109.C.ERR029405"),
      c("XF12", "XF12_18_05_11.PG0102.CW.ERR045635"),
      c("XE7", "XE7.PG0106.C.ERR029407"),
      c("XD8", "XD8_13_05_11.PG0105.CW.ERR045628"),
      c("XB3", "XB3.PG0093.C.ERR029105"), # 11 discordances
      c("WF12", "WF12.PG0097.C.ERR027109"),
      c("WE2", "WE2.PG0085.C.ERR027101"),
      c("WC4", "WC4.PG0082.C.ERR029093"),
      c("TF1", "TF1.PG0080.C.ERR027103"), # 19 discordances
      c("QF5", "QF5.PG0078.C.ERR029092"),
      c("NIC", "NIC_18_05_11.PG0095.CW.ERR045631"),
      c("NF10", "NF10.PG0096.C.ERR027108"), # 16 discordances
      c("LC12", "LC12_5_5_11.PG0110.CW.ERR045641"),
      c("LA10", "LA10_13_05_11.PG0086.CW.ERR045629"),
      c("KH7", "KH7.PG0088.C.ERR027111"),
      c("KC5", "KC5.PG0101.C.ERR029147"),
      c("KC2", "KC2_5_5_11.PG0090.CW.ERR045640"),
      c("KB8", "KB8_5_5_11.PG0104.CW.ERR045642"),
      c("KA6", "KA6.PG0091.C.ERR027117"),
      c("JON", "JON.PG0107.C.ERR029408"),
      c("JF6", "JF6.PG0079.CW.ERR045637"),
      c("JE11", "JE11_18.05_11.PG0100.CW.ERR045630"),
      c("JC9", "JC9_18_05_11.PG0111.CW.ERR045634"),
      c("JC3", "JC3.PG0077.CW.ERR045636"),
      c("JB8", "JB8.PG0087.C.ERR029091"),
      c("JB12", "JB12.PG0099.C.ERR029146"), # 73 discordances
  #    c("GB4", "GB4_NIH.PG0084.C.ERR027100"),
      c("DEV", "DEV_18_05_11.PG0081.CW.ERR045633"),
      c("DAN", "DAN.PG0098.C.ERR027110"),
      c("AUD", "AUD_5_5_11.PG0112.CW.ERR045639"),
      c("AL2", "AL2_13_05_11.PG0103.CW.ERR045627")
  #    c("7G8", "X7G8_NIH.PG0083.C.ERR027099")
    )
#    comparisonVsComparison = rbind(
#      c("X7C126_1", "X7C126_2"),
#      c("X7C126_1", "X7C126_3"),
#      c("X7C126_2", "X7C126_3"),
#      c("X7C424_1", "X7C424_2"),
#      c("CH3_61_1", "CH3_61_2"),
#      c("Dd2_1", "Dd2_2"),
#      c("Dd2_1", "Dd2_3"),
#      c("Dd2_1", "Dd2_4"),
#      c("Dd2_1", "Dd2_5"),
#      c("Dd2_2", "Dd2_3"),
#      c("Dd2_2", "Dd2_4"),
#      c("Dd2_2", "Dd2_5"),
#      c("Dd2_3", "Dd2_4"),
#      c("Dd2_3", "Dd2_5"),
#      c("Dd2_4", "Dd2_5"),
#      c("HB3_1", "HB3_2"),
#      c("HB3_1", "HB3_3"),
#      c("HB3_1", "HB3_4"),
#      c("HB3_2", "HB3_3"),
#      c("HB3_2", "HB3_4"),
#      c("HB3_3", "HB3_4"),
#      c("QC23_1", "QC23_2"),
#      c("SC05_1", "SC05_2")
#    ),
#    subjectVsSubject = rbind(
#      c("CH3_61.PG0033.C.ERR022940", "CH3_61.PG0033.Cx.ERR175544"),
#      c("GC06.PG0028.C.ERR012894", "GC06.PG0028.C.ERR015456"),
#      c("X1BB5.PG0023.C.ERR012893", "X1BB5.PG0023.C.ERR015449"),
#      c("X7C126.PG0047.C.ERR012891", "X7C126.PG0047.C.ERR015452"),
#      c("X7C46.PG0046.C.ERR022939", "X7C46.PG0046.Cx.ERR107476")
#    )
  )
) {
#  malariagenVcf <- malariagenVcf[as.character(unlist(alt(malariagenVcf))) %in% c("A", "C", "T", "G")]
  
  discordanceMatricesList <- sapply(
    names(malariagenVcfList),
    function(callSetName) {
      malariagenVcf <- malariagenVcfList[[callSetName]][geno(malariagenVcfList[[callSetName]])[["GT"]][, IDparent1] != geno(malariagenVcfList[[callSetName]])[["GT"]][, IDparent2]]
      compareCalls(
        malariagenVcf,
        uberchipVcf,
        subjectName                 = paste("MalariaGEN", callSetName, sep="."),
        comparisonName              = "Jiang et al",
        distanceThresholds          = c(0, 22),
        discordanceThreshold        = discordanceThreshold,
        discordanceProportionThreshold = discordanceProportionThreshold,
        comparisonDSthreshold       = comparisonDSthreshold,
        plotFilestem                = paste(plotFilestem, "Vs", callSetName, "_DS<=", comparisonDSthreshold, sep=""),
        IDparent1                   = IDparent1,
        IDparent2                   = IDparent2,
        shouldSubsetToBialleleic    = TRUE,
        shouldCompareRefsAndAlts    = TRUE,
        GTsToCompare                = "parentBased",
        GTsToIntMapping             = GTsToIntMapping,
        expectedMatches             = expectedMatches
      )
    },
    USE.NAMES = TRUE,
    simplify = FALSE
  )
}


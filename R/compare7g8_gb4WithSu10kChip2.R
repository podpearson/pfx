# compare7g8_gb4WithSu10kChip2.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


compare7g8_gb4WithSu10kChip2 <- function(
  malariagenVcfList             = list(
    Intersection = loadCallsSubset(
      vcfFilename                 = "data/release/1.0.combined.RC1/7g8_gb4.combined.vcf.gz",
      subsetVcfFilename           = "analysis/release/1.0.combined.RC1/7g8_gb4.combined.Intersection.vcf",
      subsetRdaFilename           = "analysis/release/1.0.combined.RC1/7g8_gb4.combined.Intersection.vcf.rda",
      subsetGrep                  = "set=Intersection"
    ),
    GATK = loadCallsSubset(
      vcfFilename                 = "data/release/1.0.GATK.RC6/7g8_gb4.gatk.both.final.vcf.gz",
      subsetVcfFilename           = "analysis/release/1.0.GATK.RC6/7g8_gb4.gatk.both.final.PASS.vcf",
      subsetRdaFilename           = "analysis/release/1.0.GATK.RC6/7g8_gb4.gatk.both.final.PASS.vcf.rda",
      subsetGrep                  = "PASS"
    ),
    Cortex = loadCallsSubset(
      vcfFilename                 = "data/release/1.0.cortex.RC1/7g8_gb4.cortex.final.vcf.gz",
      subsetVcfFilename           = "analysis/release/1.0.cortex.RC1/7g8_gb4.cortex.final.PASS.vcf",
      subsetRdaFilename           = "analysis/release/1.0.cortex.RC1/7g8_gb4.cortex.final.PASS.vcf.rda",
      subsetGrep                  = "PASS"
    )
  ),
  su10kChipVcf                = loadSu10kChipAsVcf(),
  discordanceThreshold        = 25,
  discordanceProportionThreshold = 0.15,
  plotFilestem                = "analysis/release/1.0.combined.RC1/su10kChip",
  IDparent1                   = "7G8_NIH/PG0083-C/ERR027099",
  IDparent2                   = "GB4_NIH/PG0084-C/ERR027100",
  GTsToIntMapping             = c("0"=1, "1"=2, "."=0, "./."=0),
  expectedMatches = list(
    comparisonVsSubject = rbind(
      c("X7G8_.a.4095667.82315", "X7G8_NIH.PG0083.C.ERR027099"),
      c("X7G8_.a.4095667.82251", "X7G8_NIH.PG0083.C.ERR027099"),
      c("AL2_.a.4095667.82930", "AL2_13_05_11.PG0103.CW.ERR045627"),
      c("AUD_.a.4095667.83267", "AUD.PG0112.C.ERR029406"),
      c("AUD_.a.4095667.83267", "AUD_5_5_11.PG0112.CW.ERR045639"),
      c("D2_.a.4095667.82289", "D2.PG0094.C.ERR027106"),
      c("D2_.a.4095667.82289", "D2_18_05_11.PG0094.CW.ERR045632"),
      c("DAN_.a.4095667.82179", "DAN.PG0098.C.ERR027110"),
      c("DEV_.a.4095667.82951", "DEV.PG0081.C.ERR027104"),
      c("DEV_.a.4095667.82951", "DEV_18_05_11.PG0081.CW.ERR045633"),
      c("GB4_.a.4095667.82338", "GB4_NIH.PG0084.C.ERR027100"),
      c("GB4_.a.4095667.82103", "GB4_NIH.PG0084.C.ERR027100"),
      c("JB12_.a.4095667.82318", "JB12.PG0099.C.ERR029146"),
      c("JB8_.a.4095667.82047", "JB8.PG0087.C.ERR029091"),
      c("JC3_.a.4095667.82947", "JC3.PG0077.C.ERR027112"),
      c("JC3_.a.4095667.82947", "JC3.PG0077.CW.ERR045636"),
      c("JE11_.a.4095667.82996", "JE11.PG0100.C.ERR029404"),
      c("JE11_.a.4095667.82996", "JE11_18.05_11.PG0100.CW.ERR045630"),
      c("JH6_.a.4095667.82998", "JH6.PG0113.C.ERR029410"),
      c("JH6_.a.4095667.82998", "JH6_12_05_11.PG0113.CW.ERR045626"),
      c("KA6_.a.4095667.82134", "KA6.PG0091.C.ERR027117"),
      c("KB8_.a.4095667.82997", "KB8.PG0104.C.ERR029148"),
      c("KB8_.a.4095667.82997", "KB8_5_5_11.PG0104.CW.ERR045642"),
      c("KC2_.a.4095667.82147", "KC2.PG0090.C.ERR027116"),
      c("KC2_.a.4095667.82147", "KC2_5_5_11.PG0090.CW.ERR045640"),
      c("LA10_.a.4095667.82142", "LA10.PG0086.C.ERR029090"),
      c("LA10_.a.4095667.82142", "LA10_13_05_11.PG0086.CW.ERR045629"),
      c("LC12_.a.4095667.82931", "LC12.PG0110.C.ERR171454"),
      c("LC12_.a.4095667.82931", "LC12_5_5_11.PG0110.CW.ERR045641"),
      c("NF10_.a.4095667.82283", "NF10.PG0096.C.ERR027108"),
      c("QF5_.a.4095667.83109", "QF5.PG0078.C.ERR029092"),
      c("QF5_.a.4095667.83109", "QF5.PG0078.CW.ERR045638"),
      c("WC4_.a.4095667.82184", "WC4.PG0082.C.ERR029093"),
      c("WE2_.a.4095667.82209", "WE2.PG0085.C.ERR027101"),
      c("WF12_.a.4095667.82323", "WF12.PG0097.C.ERR027109"),
      c("XE7_.a.4095667.83293", "XE7.PG0106.C.ERR029407"),
      c("XF12_.a.4095667.82106", "XF12.PG0102.C.ERR029143"),
      c("XF12_.a.4095667.82106", "XF12_18_05_11.PG0102.CW.ERR045635"),
      c("XG10_.a.4095667.82160", "XG10.PG0109.C.ERR029405")
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
        su10kChipVcf,
        subjectName                 = paste("MalariaGEN", callSetName, sep="."),
        comparisonName              = "Su10kChip",
        distanceThresholds          = c(0, 0),
        discordanceThreshold        = discordanceThreshold,
        discordanceProportionThreshold = discordanceProportionThreshold,
        plotFilestem                = paste(plotFilestem, "Vs", callSetName, sep=""),
        IDparent1                   = IDparent1,
        IDparent2                   = IDparent2,
        shouldSubsetToBialleleic    = TRUE,
        shouldCompareRefsAndAlts    = TRUE,
        GTsToCompare                = "asVcf",
        GTsToIntMapping             = GTsToIntMapping,
        expectedMatches             = expectedMatches
      )
    },
    USE.NAMES = TRUE,
    simplify = FALSE
  )
}


# compare7g8_gb4WithSu10kChip.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


compare7g8_gb4WithSu10kChip <- function(
  malariagenVcf               = loadCallsSubset(
    vcfFilename                 = "data/release/1.0.combined.RC1/7g8_gb4.combined.vcf.gz",
    subsetVcfFilename           = "analysis/release/1.0.combined.RC1/7g8_gb4.combined.Intersection.vcf",
    subsetRdaFilename           = "analysis/release/1.0.combined.RC1/7g8_gb4.combined.Intersection.vcf.rda",
    subsetGrep                  = "set=Intersection"
  ),
  su10kChipVcf                = loadSu10kChipAsVcf(),
  discordanceThreshold        = 25,
  plotFilestem                = "analysis/release/1.0.combined.RC1/su10kChipVsIntersection",
  IDparent1                   = "7G8_NIH/PG0083-C/ERR027099",
  IDparent2                   = "GB4_NIH/PG0084-C/ERR027100",
  GTsToIntMapping             = c("0"=1, "1"=2, "."=0, "./."=0)
) {
  malariagenVcf <- malariagenVcf[geno(malariagenVcf)[["GT"]][, IDparent1] != geno(malariagenVcf)[["GT"]][, IDparent2]]
  
  comparisonVsSubjectDiscordanceMatrix <- compareCalls(
    malariagenVcf,
    su10kChipVcf,
    subjectName                 = "MalariaGEN",
    comparisonName              = "Su10kChip",
    distanceThresholds          = c(0, 0),
    discordanceThreshold        = discordanceThreshold,
    plotFilestem                = plotFilestem,
    IDparent1                   = IDparent1,
    IDparent2                   = IDparent2,
    shouldSubsetToBialleleic    = TRUE,
    shouldCompareRefsAndAlts    = TRUE,
    GTsToCompare                = "asVcf",
    GTsToIntMapping             = GTsToIntMapping
  )
  browser()
  
#  #  attempt to track down which sample pairs have suspect concordance
#  which(comparisonVsSubjectDiscordanceMatrix>20 & comparisonVsSubjectDiscordanceMatrix<100, arr.ind=TRUE)
#  comparisonVsSubjectDiscordanceMatrix[28,] #JB12 is not particularly good match
#  comparisonVsSubjectDiscordanceMatrix[29,] #DEV matches better to -CW than -C version (which has low coverage)
#  which(comparisonVsSubjectDiscordanceMatrix>10 & comparisonVsSubjectDiscordanceMatrix<=20, arr.ind=TRUE)
#  comparisonVsSubjectDiscordanceMatrix[17,] #KC5 has better match to KC5 (3 discordant SNPs) than to NH11
#  comparisonVsSubjectDiscordanceMatrix[19,] #KB8 has better match to KB8 and KC5 (3 discordant SNPs) than to NH11
  
  expectedMatches <- rbind(
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
  stem(comparisonVsSubjectDiscordanceMatrix[expectedMatches])
  median(comparisonVsSubjectDiscordanceMatrix[expectedMatches])
  expectedMatches[which(comparisonVsSubjectDiscordanceMatrix[expectedMatches]>100, arr.ind=TRUE), ]
  goodMatches <- expectedMatches[-(which(comparisonVsSubjectDiscordanceMatrix[expectedMatches]>100, arr.ind=TRUE)), ]
  stem(comparisonVsSubjectDiscordanceMatrix[goodMatches])
  median(comparisonVsSubjectDiscordanceMatrix[goodMatches])
  1-median(comparisonVsSubjectDiscordanceMatrix[expectedMatches])/dim(subjectMatchesGTs)[1]
  1-median(comparisonVsSubjectDiscordanceMatrix[expectedMatches])/(dim(subjectMatchesGTs)[1]/2)
  
#  double-checking less good matches don't have a better match
  comparisonVsSubjectDiscordanceMatrix[, "JB12.PG0099.C.ERR029146"]
  comparisonVsSubjectDiscordanceMatrix[, "XB3.PG0093.C.ERR029105"]
  comparisonVsSubjectDiscordanceMatrix[, "TF1.PG0080.C.ERR027103"]
  comparisonVsSubjectDiscordanceMatrix[, "NF10.PG0096.C.ERR027108"]
  
#  double-checking less good matches don't have a better match
  comparisonVsSubjectDiscordanceMatrixList[[1]]["CH3_61_2", ]
  stem(comparisonVsSubjectDiscordanceMatrixList[[1]]["CH3_61_2", ])
  min(comparisonVsSubjectDiscordanceMatrixList[[1]]["CH3_61_2", ])
  which.min(comparisonVsSubjectDiscordanceMatrixList[[1]]["CH3_61_2", ])
  comparisonVsSubjectDiscordanceMatrixList[[1]]["QC23_2", ]
  stem(comparisonVsSubjectDiscordanceMatrixList[[1]]["QC23_2", ])
  min(comparisonVsSubjectDiscordanceMatrixList[[1]]["QC23_2", ])
  which.min(comparisonVsSubjectDiscordanceMatrixList[[1]]["QC23_2", ])

}


# compareHb3_Dd2WithUberchip.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


compareHb3_Dd2WithUberchip <- function(
  malariagenVcf               = loadCallsSubset(
    vcfFilename                 = "data/release/1.0.combined.RC1/hb3_dd2.combined.vcf.gz",
    subsetVcfFilename           = "analysis/release/1.0.combined.RC1/hb3_dd2.combined.Intersection.vcf",
    subsetRdaFilename           = "analysis/release/1.0.combined.RC1/hb3_dd2.combined.Intersection.vcf.rda",
    subsetGrep                  = "set=Intersection",
  ),
  uberchipVcf                 = loadUberchipAsVcf(),
  plotFilestem                = "analysis/release/1.0.combined.RC1/uberchipVsIntersection",
  IDparent1                   = "HB3_Ferdig/PG0004-CW/ERR012788",
  IDparent2                   = "DD2_Ferdig/PG0008-CW/ERR012840",
  GTsToIntMapping             = c("0"=1, "1"=2, "."=0, "./."=0)
) {
  malariagenVcf <- malariagenVcf[geno(malariagenVcf)[["GT"]][, IDparent1] != geno(malariagenVcf)[["GT"]][, IDparent2]]
#  malariagenVcf <- malariagenVcf[as.character(unlist(alt(malariagenVcf))) %in% c("A", "C", "T", "G")]
  comparisonVsSubjectDiscordanceMatrix <- compareCalls(
    malariagenVcf,
    uberchipVcf,
    subjectName                 = "MalariaGEN",
    comparisonName              = "Uberchip",
    distanceThresholds          = c(0, 0),
    plotFilestem                = plotFilestem,
    IDparent1                   = IDparent1,
    IDparent2                   = IDparent2,
    GTsToIntMapping             = GTsToIntMapping
  )
  browser()
  
  #  attempt to track down which sample pairs have suspect concordance
  which(comparisonVsSubjectDiscordanceMatrix>20 & comparisonVsSubjectDiscordanceMatrix<100, arr.ind=TRUE)
  comparisonVsSubjectDiscordanceMatrix[28,] #JB12 is not particularly good match
  comparisonVsSubjectDiscordanceMatrix[29,] #DEV matches better to -CW than -C version (which has low coverage)
  which(comparisonVsSubjectDiscordanceMatrix>10 & comparisonVsSubjectDiscordanceMatrix<=20, arr.ind=TRUE)
  comparisonVsSubjectDiscordanceMatrix[17,] #KC5 has better match to KC5 (3 discordant SNPs) than to NH11
  comparisonVsSubjectDiscordanceMatrix[19,] #KB8 has better match to KB8 and KC5 (3 discordant SNPs) than to NH11
  
  expectedMatches <- rbind(
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
  stem(comparisonVsSubjectDiscordanceMatrix[expectedMatches])
  median(comparisonVsSubjectDiscordanceMatrix[expectedMatches])
  1-median(comparisonVsSubjectDiscordanceMatrix[expectedMatches])/dim(subjectMatchesGTs)[1]
  1-median(comparisonVsSubjectDiscordanceMatrix[expectedMatches])/(dim(subjectMatchesGTs)[1]/2)
  
#  double-checking less good matches don't have a better match
  comparisonVsSubjectDiscordanceMatrix[, "JB12.PG0099.C.ERR029146"]
  comparisonVsSubjectDiscordanceMatrix[, "XB3.PG0093.C.ERR029105"]
  comparisonVsSubjectDiscordanceMatrix[, "TF1.PG0080.C.ERR027103"]
  comparisonVsSubjectDiscordanceMatrix[, "NF10.PG0096.C.ERR027108"]

}


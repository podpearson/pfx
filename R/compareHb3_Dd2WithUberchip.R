# compareHb3_Dd2WithUberchip.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


compareHb3_Dd2WithUberchip <- function(
  malariagenVcfList             = list(
    Intersection = loadCallsSubset(
      vcfFilename                 = "data/release/1.0.combined.RC1/hb3_dd2.combined.vcf.gz",
      subsetVcfFilename           = "analysis/release/1.0.combined.RC1/hb3_dd2.combined.Intersection.vcf",
      subsetRdaFilename           = "analysis/release/1.0.combined.RC1/hb3_dd2.combined.Intersection.vcf.rda",
      subsetGrep                  = "set=Intersection"
    ),
    Union = loadCallsSubset(
      vcfFilename                 = "data/release/1.0.combined.RC1/hb3_dd2.combined.vcf.gz",
      subsetVcfFilename           = "analysis/release/1.0.combined.RC1/hb3_dd2.combined.Union.vcf",
      subsetRdaFilename           = "analysis/release/1.0.combined.RC1/hb3_dd2.combined.Union.vcf.rda",
      subsetGrep                  = ""
    ),
    GATK = loadCallsSubset(
      vcfFilename                 = "data/release/1.0.GATK.RC6/hb3_dd2.gatk.both.final.vcf.gz",
      subsetVcfFilename           = "analysis/release/1.0.GATK.RC6/hb3_dd2.gatk.both.final.PASS.vcf",
      subsetRdaFilename           = "analysis/release/1.0.GATK.RC6/hb3_dd2.gatk.both.final.PASS.vcf.rda",
      subsetGrep                  = "PASS"
    ),
    GATKunfiltered = loadCallsSubset(
      vcfFilename                 = "data/release/1.0.GATK.RC6/hb3_dd2.gatk.both.final.vcf.gz",
      subsetVcfFilename           = "analysis/release/1.0.GATK.RC6/hb3_dd2.gatk.both.final.unfiltered.vcf",
      subsetRdaFilename           = "analysis/release/1.0.GATK.RC6/hb3_dd2.gatk.both.final.unfiltered.vcf.rda",
      subsetGrep                  = "",
      reload                      = TRUE
    ),
    Cortex = loadCallsSubset(
      vcfFilename                 = "data/release/1.0.cortex.RC1/hb3_dd2.cortex.final.vcf.gz",
      subsetVcfFilename           = "analysis/release/1.0.cortex.RC1/hb3_dd2.cortex.final.PASS.vcf",
      subsetRdaFilename           = "analysis/release/1.0.cortex.RC1/hb3_dd2.cortex.final.PASS.vcf.rda",
      subsetGrep                  = "PASS"
    ),
    CortexUnfiltered = loadCallsSubset(
      vcfFilename                 = "data/release/1.0.cortex.RC1/hb3_dd2.cortex.final.vcf.gz",
      subsetVcfFilename           = "analysis/release/1.0.cortex.RC1/hb3_dd2.cortex.final.unfiltered.vcf",
      subsetRdaFilename           = "analysis/release/1.0.cortex.RC1/hb3_dd2.cortex.final.unfiltered.vcf.rda",
      subsetGrep                  = "",
      reload                      = TRUE
    )
  ),
#  malariagenVcf               = loadCallsSubset(
#    vcfFilename                 = "data/release/1.0.combined.RC1/hb3_dd2.combined.vcf.gz",
#    subsetVcfFilename           = "analysis/release/1.0.combined.RC1/hb3_dd2.combined.Intersection.vcf",
#    subsetRdaFilename           = "analysis/release/1.0.combined.RC1/hb3_dd2.combined.Intersection.vcf.rda",
#    subsetGrep                  = "set=Intersection"
#  ),
  uberchipVcf                 = loadUberchipAsVcf(),
  discordanceThreshold        = 200,
  discordanceProportionThreshold = 0.15,
  comparisonDSthresholds      = c(0.5),
#  comparisonDSthresholds      = c(1.0, 0.5, 0.2, 0.1),
  plotFilestem                = "analysis/release/1.0.combined.RC1/uberchip",
  IDparent1                   = "HB3_Ferdig/PG0004-CW/ERR012788",
  IDparent2                   = "DD2_Ferdig/PG0008-CW/ERR012840",
  IDcomparisonParent1         = "HB3_3",
  IDcomparisonParent2         = "Dd2_5",
  GTsToIntMapping             = c("0"=1, "1"=2, "."=0, "./."=0),
  expectedMatches = list(
    comparisonVsSubject = rbind(
      c("CH3_61_1", "CH3_61.PG0033.C.ERR022940"),
      c("CH3_61_1", "CH3_61.PG0033.Cx.ERR175544"),
      c("CH3_61_2", "CH3_61.PG0033.C.ERR022940"),
      c("CH3_61_2", "CH3_61.PG0033.Cx.ERR175544"),
      c("Dd2_1", "DD2_Ferdig.PG0008.CW.ERR012840"),
      c("Dd2_2", "DD2_Ferdig.PG0008.CW.ERR012840"),
      c("Dd2_3", "DD2_Ferdig.PG0008.CW.ERR012840"),
      c("Dd2_4", "DD2_Ferdig.PG0008.CW.ERR012840"),
      c("Dd2_5", "DD2_Ferdig.PG0008.CW.ERR012840"),
      c("HB3_1", "HB3_Ferdig.PG0004.CW.ERR012788"),
      c("HB3_2", "HB3_Ferdig.PG0004.CW.ERR012788"),
      c("HB3_3", "HB3_Ferdig.PG0004.CW.ERR012788"),
      c("HB3_4", "HB3_Ferdig.PG0004.CW.ERR012788"),
      c("QC23_1", "QC23.PG0045.C.ERR012892"),
      c("QC23_2", "QC23.PG0045.C.ERR012892"),
      c("SC05_1", "SC05.PG0019.C.ERR019051"),
      c("SC05_2", "SC05.PG0019.C.ERR019051"),
      c("X7C126_1", "X7C126.PG0047.C.ERR015452"),
      c("X7C126_1", "X7C126.PG0047.C.ERR012891"),
      c("X7C126_2", "X7C126.PG0047.C.ERR015452"),
      c("X7C126_2", "X7C126.PG0047.C.ERR012891"),
      c("X7C126_3", "X7C126.PG0047.C.ERR015452"),
      c("X7C126_3", "X7C126.PG0047.C.ERR012891"),
      c("X7C424_1", "X7C424.PG0044.C.ERR019043"),
      c("X7C424_2", "X7C424.PG0044.C.ERR019043")
    ),
    comparisonVsComparison = rbind(
      c("X7C126_1", "X7C126_2"),
      c("X7C126_1", "X7C126_3"),
      c("X7C126_2", "X7C126_3"),
      c("X7C424_1", "X7C424_2"),
      c("CH3_61_1", "CH3_61_2"),
      c("Dd2_1", "Dd2_2"),
      c("Dd2_1", "Dd2_3"),
      c("Dd2_1", "Dd2_4"),
      c("Dd2_1", "Dd2_5"),
      c("Dd2_2", "Dd2_3"),
      c("Dd2_2", "Dd2_4"),
      c("Dd2_2", "Dd2_5"),
      c("Dd2_3", "Dd2_4"),
      c("Dd2_3", "Dd2_5"),
      c("Dd2_4", "Dd2_5"),
      c("HB3_1", "HB3_2"),
      c("HB3_1", "HB3_3"),
      c("HB3_1", "HB3_4"),
      c("HB3_2", "HB3_3"),
      c("HB3_2", "HB3_4"),
      c("HB3_3", "HB3_4"),
      c("QC23_1", "QC23_2"),
      c("SC05_1", "SC05_2")
    ),
    subjectVsSubject = rbind(
      c("CH3_61.PG0033.C.ERR022940", "CH3_61.PG0033.Cx.ERR175544"),
      c("GC06.PG0028.C.ERR012894", "GC06.PG0028.C.ERR015456"),
      c("X1BB5.PG0023.C.ERR012893", "X1BB5.PG0023.C.ERR015449"),
      c("X7C126.PG0047.C.ERR012891", "X7C126.PG0047.C.ERR015452"),
      c("X7C46.PG0046.C.ERR022939", "X7C46.PG0046.Cx.ERR107476")
    )
  )
) {
#  malariagenVcf <- malariagenVcf[as.character(unlist(alt(malariagenVcf))) %in% c("A", "C", "T", "G")]
  
  discordanceMatricesList <- sapply(
    names(malariagenVcfList),
    function(callSetName) {
      malariagenVcf <- malariagenVcfList[[callSetName]][geno(malariagenVcfList[[callSetName]])[["GT"]][, IDparent1] != geno(malariagenVcfList[[callSetName]])[["GT"]][, IDparent2]]
      comparisonVsSubjectDiscordanceMatrixList <- sapply(
        comparisonDSthresholds,
        function(comparisonDSthreshold) {
          compareCalls(
            malariagenVcf,
            uberchipVcf,
            subjectName                 = paste("MalariaGEN", callSetName, sep="."),
            comparisonName              = "Uberchip",
            distanceThresholds          = c(0, 0),
            discordanceThreshold        = discordanceThreshold,
            discordanceProportionThreshold = discordanceProportionThreshold,
            comparisonDSthreshold       = comparisonDSthreshold,
            plotFilestem                = paste(plotFilestem, "Vs", callSetName, "_DS<=", comparisonDSthreshold, sep=""),
            IDparent1                   = IDparent1,
            IDparent2                   = IDparent2,
            IDcomparisonParent1         = IDcomparisonParent1,
            IDcomparisonParent2         = IDcomparisonParent2,
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
      names(comparisonVsSubjectDiscordanceMatrixList) <- comparisonDSthresholds
      comparisonVsSubjectDiscordanceMatrixList
    },
    USE.NAMES = TRUE,
    simplify = FALSE
  )
  
  return(discordanceMatricesList)
  
#  comparisonVsSubjectDiscordanceMatrixList <- sapply(
#    comparisonDSthresholds,
#    function(comparisonDSthreshold) {
#      compareCalls(
#        malariagenVcf,
#        uberchipVcf,
#        subjectName                 = "MalariaGEN",
#        comparisonName              = "Uberchip",
#        distanceThresholds          = c(0, 0),
#        discordanceThreshold        = discordanceThreshold,
#        discordanceProportionThreshold = discordanceProportionThreshold,
#        comparisonDSthreshold       = comparisonDSthreshold,
#        plotFilestem                = paste(plotFilestem, comparisonDSthreshold, sep=""),
#        IDparent1                   = IDparent1,
#        IDparent2                   = IDparent2,
#        shouldSubsetToBialleleic    = TRUE,
#        shouldCompareRefsAndAlts    = TRUE,
#        GTsToCompare                = "asVcf",
#        GTsToIntMapping             = GTsToIntMapping
#      )
#    },
#    USE.NAMES = TRUE,
#    simplify = FALSE
#  )
  browser()
  
#  #  attempt to track down which sample pairs have suspect concordance
#  which(comparisonVsSubjectDiscordanceMatrix>20 & comparisonVsSubjectDiscordanceMatrix<100, arr.ind=TRUE)
#  comparisonVsSubjectDiscordanceMatrix[28,] #JB12 is not particularly good match
#  comparisonVsSubjectDiscordanceMatrix[29,] #DEV matches better to -CW than -C version (which has low coverage)
#  which(comparisonVsSubjectDiscordanceMatrix>10 & comparisonVsSubjectDiscordanceMatrix<=20, arr.ind=TRUE)
#  comparisonVsSubjectDiscordanceMatrix[17,] #KC5 has better match to KC5 (3 discordant SNPs) than to NH11
#  comparisonVsSubjectDiscordanceMatrix[19,] #KB8 has better match to KB8 and KC5 (3 discordant SNPs) than to NH11
  
  expectedMatches <- rbind(
    c("CH3_61_1", "CH3_61.PG0033.C.ERR022940"),
    c("CH3_61_1", "CH3_61.PG0033.Cx.ERR175544"),
    c("CH3_61_2", "CH3_61.PG0033.C.ERR022940"),
    c("CH3_61_2", "CH3_61.PG0033.Cx.ERR175544"),
    c("Dd2_1", "DD2_Ferdig.PG0008.CW.ERR012840"),
    c("Dd2_2", "DD2_Ferdig.PG0008.CW.ERR012840"),
    c("Dd2_3", "DD2_Ferdig.PG0008.CW.ERR012840"),
    c("Dd2_4", "DD2_Ferdig.PG0008.CW.ERR012840"),
    c("Dd2_5", "DD2_Ferdig.PG0008.CW.ERR012840"),
    c("HB3_1", "HB3_Ferdig.PG0004.CW.ERR012788"),
    c("HB3_2", "HB3_Ferdig.PG0004.CW.ERR012788"),
    c("HB3_3", "HB3_Ferdig.PG0004.CW.ERR012788"),
    c("HB3_4", "HB3_Ferdig.PG0004.CW.ERR012788"),
    c("QC23_1", "QC23.PG0045.C.ERR012892"),
    c("QC23_2", "QC23.PG0045.C.ERR012892"),
    c("SC05_1", "SC05.PG0019.C.ERR019051"),
    c("SC05_2", "SC05.PG0019.C.ERR019051"),
    c("X7C126_1", "X7C126.PG0047.C.ERR015452"),
    c("X7C126_1", "X7C126.PG0047.C.ERR012891"),
    c("X7C126_2", "X7C126.PG0047.C.ERR015452"),
    c("X7C126_2", "X7C126.PG0047.C.ERR012891"),
    c("X7C126_3", "X7C126.PG0047.C.ERR015452"),
    c("X7C126_3", "X7C126.PG0047.C.ERR012891"),
    c("X7C424_1", "X7C424.PG0044.C.ERR019043"),
    c("X7C424_2", "X7C424.PG0044.C.ERR019043")
  )
  
  numberOfComparisons <- list()
  numberOfComparisons[["Intersection"]] <- discordanceMatricesList[["Intersection"]][["0.5"]][["comparisonVsSubjectDiscordanceMatrix"]][expectedMatches] / discordanceMatricesList[["Intersection"]][["0.5"]][["comparisonVsSubjectDiscordanceProportionMatrix"]][expectedMatches]
  numberOfComparisons[["GATK"]] <- discordanceMatricesList[["GATK"]][["0.5"]][["comparisonVsSubjectDiscordanceMatrix"]][expectedMatches] / discordanceMatricesList[["GATK"]][["0.5"]][["comparisonVsSubjectDiscordanceProportionMatrix"]][expectedMatches]
  numberOfComparisons[["Cortex"]] <- discordanceMatricesList[["Cortex"]][["0.5"]][["comparisonVsSubjectDiscordanceMatrix"]][expectedMatches] / discordanceMatricesList[["Cortex"]][["0.5"]][["comparisonVsSubjectDiscordanceProportionMatrix"]][expectedMatches]
  
  stem(discordanceMatricesList[["Intersection"]][["0.5"]][["comparisonVsSubjectDiscordanceProportionMatrix"]][expectedMatches])
  stem(discordanceMatricesList[["GATK"]][["0.5"]][["comparisonVsSubjectDiscordanceProportionMatrix"]][expectedMatches])
  stem(discordanceMatricesList[["Cortex"]][["0.5"]][["comparisonVsSubjectDiscordanceProportionMatrix"]][expectedMatches])
  median(discordanceMatricesList[["Intersection"]][["0.5"]][["comparisonVsSubjectDiscordanceProportionMatrix"]][expectedMatches])
  median(discordanceMatricesList[["GATK"]][["0.5"]][["comparisonVsSubjectDiscordanceProportionMatrix"]][expectedMatches])
  median(discordanceMatricesList[["Cortex"]][["0.5"]][["comparisonVsSubjectDiscordanceProportionMatrix"]][expectedMatches])
  
  median(discordanceMatricesList[["Intersection"]][["0.5"]][["comparisonVsSubjectDiscordanceMatrix"]][expectedMatches])
  median(discordanceMatricesList[["GATK"]][["0.5"]][["comparisonVsSubjectDiscordanceMatrix"]][expectedMatches])
  median(discordanceMatricesList[["Cortex"]][["0.5"]][["comparisonVsSubjectDiscordanceMatrix"]][expectedMatches])
  
  
  stem(comparisonVsSubjectDiscordanceMatrixList[[1]][expectedMatches])
  median(comparisonVsSubjectDiscordanceMatrixList[[1]][expectedMatches])
  expectedMatches[which(comparisonVsSubjectDiscordanceMatrixList[[1]][expectedMatches] > 1000, arr.ind=TRUE), ]
#  1-median(comparisonVsSubjectDiscordanceMatrixList[[1]][expectedMatches])/dim(subjectMatchesGTs)[1]
#  1-median(comparisonVsSubjectDiscordanceMatrixList[[1]][expectedMatches])/(dim(subjectMatchesGTs)[1]/2)
  
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


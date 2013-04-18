# comparisonsWithArrayData.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


comparisonsWithArrayData <- function(
  reload                      = FALSE
) {
  if(!file.exists("analysis/release/comparisonsWithArrayData")) {
    dir.create("analysis/release/comparisonsWithArrayData")
  }
  if(reload || !file.exists("analysis/release/comparisonsWithArrayData/Jiang_results.rda")) {
    Jiang_results <- compare7g8_gb4WithJiang2()
    save(Jiang_results, file="analysis/release/comparisonsWithArrayData/Jiang_results.rda")
  } else {
    load("analysis/release/comparisonsWithArrayData/Jiang_results.rda")
  }
  gc()
  if(reload || !file.exists("analysis/release/comparisonsWithArrayData/Su10k_results.rda")) {
    Su10k_results <- compare7g8_gb4WithSu10kChip2()
    save(Su10k_results, file="analysis/release/comparisonsWithArrayData/Su10k_results.rda")
  } else {
    load("analysis/release/comparisonsWithArrayData/Su10k_results.rda")
  }
  Su10k_sensitivity <- data.frame(t(sapply(Su10k_results, function(x) x[["sensitivityResults"]])))
  print(Su10k_sensitivity)
  Su10k_sensitivity[["Position matches"]] <- paste(
    Su10k_sensitivity[["numberOfSegregatingSNPsWithPositionMatches"]],
    " (",
    round((Su10k_sensitivity[["numberOfSegregatingSNPsWithPositionMatches"]]/Su10k_sensitivity[["numberOfSegregatingSNPsInComparison"]])*100),
    "%)",
    sep=""
  )
  Su10k_sensitivity[["Position and allele matches"]] <- paste(
    Su10k_sensitivity[["numberOfSegregatingSNPsWithPositionAndAlleleMatches"]],
    " (",
    round((Su10k_sensitivity[["numberOfSegregatingSNPsWithPositionAndAlleleMatches"]]/Su10k_sensitivity[["numberOfSegregatingSNPsInComparison"]])*100),
    "%)",
    sep=""
  )
  Su10k_sensitivity[["numberOfSegregatingSNPsWithPositionMatches"]] <- NULL
  Su10k_sensitivity[["numberOfSegregatingSNPsWithPositionAndAlleleMatches"]] <- NULL
  Su10k_sensitivity[["numberOfSegregatingSNPsInComparison"]] <- NULL
  write.table(
    Su10k_sensitivity,
    file="analysis/release/comparisonsWithArrayData/Su10k_sensitivity.txt",
    sep="\t"
  )
  gc()
  
  if(reload || !file.exists("analysis/release/comparisonsWithArrayData/Uberchip_results.rda")) {
    Uberchip_results <- compareHb3_Dd2WithUberchip()
    save(Uberchip_results, file="analysis/release/comparisonsWithArrayData/Uberchip_results.rda")
  } else {
    load("analysis/release/comparisonsWithArrayData/Uberchip_results.rda")
  }
  Uberchip_sensitivity <- data.frame(t(sapply(Uberchip_results, function(x) x[["0.5"]][["sensitivityResults"]])))
  print(Uberchip_sensitivity)
  Uberchip_sensitivity[["Position matches"]] <- paste(
    Uberchip_sensitivity[["numberOfSegregatingSNPsWithPositionMatches"]],
    " (",
    round((Uberchip_sensitivity[["numberOfSegregatingSNPsWithPositionMatches"]]/Uberchip_sensitivity[["numberOfSegregatingSNPsInComparison"]])*100),
    "%)",
    sep=""
  )
  Uberchip_sensitivity[["Position and allele matches"]] <- paste(
    Uberchip_sensitivity[["numberOfSegregatingSNPsWithPositionAndAlleleMatches"]],
    " (",
    round((Uberchip_sensitivity[["numberOfSegregatingSNPsWithPositionAndAlleleMatches"]]/Uberchip_sensitivity[["numberOfSegregatingSNPsInComparison"]])*100),
    "%)",
    sep=""
  )
  Uberchip_sensitivity[["numberOfSegregatingSNPsWithPositionMatches"]] <- NULL
  Uberchip_sensitivity[["numberOfSegregatingSNPsWithPositionAndAlleleMatches"]] <- NULL
  Uberchip_sensitivity[["numberOfSegregatingSNPsInComparison"]] <- NULL
  write.table(
    Uberchip_sensitivity,
    file="analysis/release/comparisonsWithArrayData/Uberchip_sensitivity.txt",
    sep="\t"
  )
  gc()
  
  return(
    list(
      Jiang_results=Jiang_results,
      Su10k_results=Su10k_results,
      Uberchip_results=Uberchip_results
    )
  )
}

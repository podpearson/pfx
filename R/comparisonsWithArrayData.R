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
  print(t(sapply(Su10k_results, function(x) x[["sensitivityResults"]])))
  write.table(
    t(sapply(Su10k_results, function(x) x[["sensitivityResults"]])),
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
  print(t(sapply(Uberchip_results, function(x) x[["0.5"]][["sensitivityResults"]])))
  write.table(
    t(sapply(Su10k_results, function(x) x[["sensitivityResults"]])),
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

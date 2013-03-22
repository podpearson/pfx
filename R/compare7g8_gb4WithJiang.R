# compare7g8_gb4WithJiang.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


compare7g8_gb4WithJiang <- function(
  malariagenVcf               = loadIntersectionCalls(),
  jiangVcf                    = loadJiangGenotypesAsVcf2(),
  plotFilestem                = "analysis/release/1.0.combined.RC1/jiangVsIntersection"
) {
  malariagenVcf <- malariagenVcf[geno(malariagenVcf)[["GT"]][, "7G8_NIH/PG0083-C/ERR027099"] != geno(malariagenVcf)[["GT"]][, "GB4_NIH/PG0084-C/ERR027100"]]
#  malariagenVcf <- malariagenVcf[as.character(unlist(alt(malariagenVcf))) %in% c("A", "C", "T", "G")]
  compareCalls(malariagenVcf, jiangVcf, plotFilestem=plotFilestem, IDparent1="7G8_NIH/PG0083-C/ERR027099", IDparent2="GB4_NIH/PG0084-C/ERR027100")
}

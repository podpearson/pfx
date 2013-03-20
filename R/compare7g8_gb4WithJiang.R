# compare7g8_gb4WithJiang.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


compare7g8_gb4WithJiang <- function(
  malariagenVcf               = loadIntersectionCalls(),
  jiangVcf                    = loadJiangGenotypesAsVcf2()
) {
  compareCalls(malariagenVcf, jiangVcf)
}

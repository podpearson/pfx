# setHaplotypeLengths.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


setHaplotypeLengths <- function(
  vcf,
  GTsToIntMapping             = c("0"=1, "1"=2, "."=0, "./."=0, "2"=0, "3"=0) # "./." is needed as sometimes this is output by GATK's UG (presumably a bug). "2", "3", needed for the case of multi-allelic sites
) {
  parentBasedGTs <- convertGTsIntToParentBasedGTs(
    genotypeCallsFromGTas012(
      vcf,
      GTsToIntMapping             = GTsToIntMapping
    ),
    return0asNA                 = TRUE
  )
  browser()
  apply(
    parentBasedGTs,
    2,
    function(snpsForThisSample) {
      
    }
  )

}

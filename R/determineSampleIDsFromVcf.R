# determineSampleIDsFromVcf.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


determineSampleIDsFromVcf <- function(
  vcfFilename                 = "/data/galton/users/rpearson/crossesTesting/release/7g8xGb4-qcPlusSamples-0.1.vcf.gz"
) {
  vcfHeaderLine <- strsplit(paste(system(paste("zcat", vcfFilename, "| head -1000 | grep ^# | tail -1"), intern=TRUE), collapse=""), "\t")[[1]]
  sampleIDs <- vcfHeaderLine[-(1:9)]
  return(sampleIDs)
}

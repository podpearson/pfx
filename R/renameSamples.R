# renameSamples.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


renameSamples <- function(
  vcf,
  sampleIDmappings            = createSampleIDmappings(
    sampleIDs                   = dimnames(vcf)[[2]],
    sampleIDcolumn              = "ena_run_accession",
    sampleIDmappingsColumn      = sampleIDcolumn
  )
) {
  dimnames(vcf)[[2]] <- names(sampleIDmappings)
  vcf
}

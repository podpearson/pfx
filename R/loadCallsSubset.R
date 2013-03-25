# loadCallsSubset.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


loadCallsSubset <- function(
  vcfFilename                 = "data/release/1.0.combined.RC1/7g8_gb4.combined.vcf.gz",
  subsetVcfFilename           = "analysis/release/1.0.combined.RC1/7g8_gb4.combined.Intersection.vcf",
  subsetRdaFilename           = "analysis/release/1.0.combined.RC1/7g8_gb4.combined.Intersection.vcf.rda",
  subsetGrep                  = "set=Intersection",
  reload                      = FALSE
) {
  if(reload || !file.exists(subsetVcfFilename)) {
    systemCommand <- sprintf("(zcat %s | head -1000 | grep ^#; zgrep '%s' %s) > %s", vcfFilename, subsetGrep, vcfFilename, subsetVcfFilename)
    cat(systemCommand)
    system(systemCommand)
  }
  if(reload || !file.exists(subsetRdaFilename)) {
    vcf <- readVcf(subsetVcfFilename, genome="Pf")
    save(vcf, file=subsetRdaFilename)
  } else {
    load(subsetRdaFilename)
  }
  return(vcf)
}


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
  if(reload || !file.exists(intersectionVcfFilename)) {
    systemCommand <- sprintf("(zcat %s | head -1000 | grep ^#; zgrep '%s' %s) > %s", vcfFilename, subsetGrep, vcfFilename, intersectionVcfFilename)
    cat(systemCommand)
    system(systemCommand)
  }
  if(reload || !file.exists(intersectionRdaFilename)) {
    vcf <- readVcf(intersectionVcfFilename, genome="Pf")
    save(vcf, file=intersectionRdaFilename)
  } else {
    load(intersectionRdaFilename)
  }
  return(vcf)
}


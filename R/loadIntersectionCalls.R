# loadIntersectionCalls.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


loadIntersectionCalls <- function(
  vcfFilename                 = "data/release/1.0.combined.RC1/7g8_gb4.combined.vcf.gz",
  intersectionVcfFilename     = "analysis/release/1.0.combined.RC1/7g8_gb4.combined.Intersection.vcf",
  intersectionRdaFilename     = "analysis/release/1.0.combined.RC1/7g8_gb4.combined.Intersection.vcf.rda",
  reload                      = FALSE
) {
  if(reload || !file.exists(intersectionVcfFilename)) {
    systemCommand <- sprintf("(zcat %s | head -1000 | grep ^#; zgrep 'set=Intersection' %s) > %s", vcfFilename, vcfFilename, intersectionVcfFilename)
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

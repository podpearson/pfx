# combineVcfListIntoVcf.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


combineVcfListIntoVcf <- function(
  vcfList
) {
  VCF(
    rowData  = do.call(c, lapply(vcfList, rowData)),
    colData  = colData(vcfList[[1]]),
    exptData = exptData(vcfList[[1]]),
    fixed    = do.call(BiocGenerics::rbind, lapply(vcfList, function(vcf) values(fixed(vcf)))),
    info     = do.call(BiocGenerics::rbind, lapply(vcfList, function(vcf) values(info(vcf)))),
    geno     = SimpleList(
      sapply(
        names(geno(vcfList[[1]])),
        function(genoField) {
          do.call(BiocGenerics::rbind, lapply(vcfList, function(vcf) geno(vcf)[[genoField]]))
        },
        simplify=FALSE,
        USE.NAMES=TRUE
      )
    )
  )
}

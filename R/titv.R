# titv.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


titv <- function(
  vcf,
  shouldIncludeATtransversions = TRUE
) {
  REFs <- as.character(ref(vcf))
  ALTs <- as.character(unlist(alt(vcf)))[do.call(c, sapply(elementLengths(alt(vcf)), seq, simplify=FALSE))==1]
  transitions <- REFs=="A" & ALTs=="G" | REFs=="G" & ALTs=="A" | REFs=="C" & ALTs=="T" | REFs=="T" & ALTs=="C"
  if(shouldIncludeATtransversions) {
    transversions <- REFs=="A" & ALTs %in% c("T", "C") | REFs=="G" & ALTs %in% c("T", "C") | REFs=="C" & ALTs %in% c("A", "G") | REFs=="T" & ALTs %in% c("A", "G")
  } else {
    transversions <- REFs=="A" & ALTs %in% c("C") | REFs=="G" & ALTs %in% c("T", "C") | REFs=="C" & ALTs %in% c("A", "G") | REFs=="T" & ALTs %in% c("G")
  }
  titvRatio <- length(which(transitions)) / length(which(transversions))
  return(titvRatio)
}

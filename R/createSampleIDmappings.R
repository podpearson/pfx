# createSampleIDmappings.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


createSampleIDmappings <- function(
  vcfFilename                 = "/data/galton/users/rpearson/crossesTesting/release/7g8xGb4-qcPlusSamples-0.1.vcf.gz",
  sampleAnnotationFilename    = "/data/malariagen2/plasmodium/pf-crosses/meta/qcmeta_annotated.tsv"
) {
  sampleIDs <- determineSampleIDsFromVcf(vcfFilename)
  sampleAnnotation <- readSampleAnnotation(sampleAnnotationFilename)
  sampleIDmappings <- sub(
    "^([^_]+)_.*$",
    "\\1",
    sampleAnnotation[sampleIDs, "source_code"]
  )
  names(sampleIDmappings) <- paste(
    sampleAnnotation[sampleIDs, "source_code"],
    " (",
    sampleAnnotation[sampleIDs, "ox_code"],
    ", ",
    sampleAnnotation[sampleIDs, "ena_run_accession"],
    ")",
    sep=""
  )
  return(sampleIDmappings)
}

# createSampleIDmappings.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


createSampleIDmappings <- function(
  vcfFilename                 = "/data/galton/users/rpearson/crossesTesting/release/7g8xGb4-qcPlusSamples-0.1.vcf.gz",
  shouldUseSampleAnnotation   = TRUE,
  sampleAnnotationFilename    = "/data/malariagen2/plasmodium/pf-crosses/meta/qcmeta_annotated.tsv",
  sampleDuplicates        = c(
    "JF6" = "JF6_KC2",
    "KC2" = "JF6_KC2",
    "AUD" = "AUD_LC12",
    "LC12" = "AUD_LC12",
    "KB8" = "KB8_KC5_NH11",
    "KC5" = "KB8_KC5_NH11",
    "NH11" = "KB8_KC5_NH11",
    "D2" = "D2_TF1",
    "TF1" = "D2_TF1",
    "7G8" = "7G8_JH6",
    "JH6" = "7G8_JH6"
  )
) {
  sampleIDs <- determineSampleIDsFromVcf(vcfFilename)
  if(shouldUseSampleAnnotation) {
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
  } else {
    sampleIDmappings <- sampleIDs
    names(sampleIDmappings) <- sampleIDs
  }
  if(!is.null(sampleDuplicates)) {
    matches <- match(sampleIDmappings, names(sampleDuplicates))
    hasMatch <- which(!is.na(matches))
    if(length(hasMatch) > 0) {
      sampleIDmappings[hasMatch] <- sampleDuplicates[matches[!is.na(matches)]]
    }
  }
  return(sampleIDmappings)
}

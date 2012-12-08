# readSampleAnnotation.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


readSampleAnnotation <- function(
  sampleAnnotationFilename    = "/data/malariagen2/plasmodium/pf-crosses/meta/qcmeta_annotated.tsv",
  sampleIDcolumn              = "ox_code"
) {
  sampleAnnotation <- read.delim(sampleAnnotationFilename)
  sampleAnnotation[["cross"]] <- ifelse(sampleAnnotation[["project_code"]]=="PFproj1", "Hb3xDd2", ifelse(sampleAnnotation[["project_code"]]=="PFproj2", "3d7xHb3", "7g8xGb4"))
  sampleAnnotation[["qcStatus"]] <- ifelse(sampleAnnotation[["ajm_qc"]]=="fail", "fail", "pass")
  sampleAnnotation <- sampleAnnotation[!grepl("not a cross progeny sample", sampleAnnotation[["ajm_qc_notes"]]), ]
  sampleAnnotation <- sampleAnnotation[!duplicated(sampleAnnotation[[sampleIDcolumn]]), ]
  row.names(sampleAnnotation) <- sampleAnnotation[[sampleIDcolumn]]
  return(sampleAnnotation)
}

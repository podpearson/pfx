# readSingleChromosomeWithRegionsMask.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


readSingleChromosomeWithRegionsMask <- function(
  vcfFilename                 = "/data/malariagen2/plasmodium/pf-crosses/data/3d7_v3/bwa_n0.01_k4_l32/genotypes_analysis_20130121/per_sample_realigned/gatk/3d7_hb3/snps/",
  chromosome                  = "Pf3D7_01_v3",
  genoToLoad                  = c("GT", "AD", "GQ", "MQ0"),
  nonVarRegions               = nonVarRegions_v3(),
  outputRdaFilename           = sub("\\.vcf[\\.gz]*$", paste("\\.coreAllSamples\\.rda", sep=""), vcfFilename),
  overwriteExisting           = !(file.exists(outputRdaFilename)),
  saveAsRobjectFile           = TRUE,
  parentalStrains             = NULL
) {
  if(vcfFilename == outputRdaFilename) {
    stop("Input and output filename are the same")
  }
  if(overwriteExisting) {
    rng <- nonVarRegions[seqnames(nonVarRegions)==chromosome]
#    rng <- GRanges(seqnames=chromosome, ranges=nonVarRegions[seqnames(nonVarRegions)==chromosome])
    param <- ScanVcfParam(which=rng, geno=c(genoToLoad))
    if(grepl("\\.vcf$", vcfFilename)) {
      if(!file.exists(paste(vcfFilename, "gz", sep="."))) {
        bgzip(vcfFilename)
      }
      vcfFilename <- paste(vcfFilename, "gz", sep=".")
    }
    if(grepl("\\.gz$", vcfFilename) & !file.exists(paste(vcfFilename, "tbi", sep="."))) {
      indexTabix(vcfFilename, "vcf4")
    }
    vcf <- readVcf(vcfFilename, "Pf", param)
    if(saveAsRobjectFile) {
      save(vcf, file=outputRdaFilename)
    }
  } else {
    load(outputRdaFilename)
  }
  if(!is.null(parentalStrains)) {
    newSampleOrdering <- c(parentalStrains, setdiff(dimnames(vcf)[[2]], parentalStrains))
    vcf <- vcf[, newSampleOrdering]
  }
  return(vcf)
}

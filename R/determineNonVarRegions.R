# determineNonVarRegions.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


determineNonVarRegions <- function(
  vcfFilename                 = "/data/malariagen2/plasmodium/pf-crosses/data/3d7_v3/bwa_n0.01_k4_l32/genotypes_analysis_20130121/per_sample_realigned/gatk/3d7_hb3/snps/Pf3D7_01_v3.annotated.vcf.gz",
  nonVarCriterion             = list(column="RegionType", operator="==", value="Conserved")
) {
  param <- ScanVcfParam(geno=NA, info=nonVarCriterion[["column"]])
  vcf <- readVcf(vcfFilename, "Pf", param)
  rng <- rowData(vcf[values(info(vcf))[[nonVarCriterion[["column"]]]] == nonVarCriterion[["value"]]])
  return(rng)
}

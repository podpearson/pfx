# haplotypesPerSample.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


haplotypesPerSample <- function(
  vcf,
  GTsToIntMapping             = c("0"=1, "1"=2, "."=0, "./."=0, "2"=0, "3"=0), # "./." is needed as sometimes this is output by GATK's UG (presumably a bug). "2", "3", needed for the case of multi-allelic sites
  createInfoFieldSummaries    = TRUE,
  maxSNPsInShortHaplotype     = 10
) {
  parentBasedGTs <- convertGTsIntToParentBasedGTs(
    genotypeCallsFromGTas012(
      vcf,
      GTsToIntMapping             = GTsToIntMapping
    ),
    reverseOrderOfSamples       = FALSE,
    return0asNA                 = TRUE,
    returnNonSegregatingAsNA    = TRUE
  )
  parentBasedGTsSplitByChromosome <- split(as.data.frame(parentBasedGTs), as.character(seqnames(vcf)))
  haplotypesPerSampleAndChromosomeMatrix <- sapply(
    parentBasedGTsSplitByChromosome,
    function(parentBasedGTs) {
      apply(
        parentBasedGTs,
        2,
        function(snpsForThisSample) {
          NAGTs <- is.na(snpsForThisSample)
          GTsRle <- Rle(snpsForThisSample[!NAGTs])
          nrun(GTsRle)
        }
      )
    }
  )
}


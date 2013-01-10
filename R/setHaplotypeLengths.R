# setHaplotypeLengths.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


setHaplotypeLengths <- function(
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
  HLsMatrixList <- lapply(
    parentBasedGTsSplitByChromosome,
    function(parentBasedGTs) {
      apply(
        parentBasedGTs,
        2,
        function(snpsForThisSample) {
          NAGTs <- is.na(snpsForThisSample)
          GTsRle <- Rle(snpsForThisSample[!NAGTs])
          haplotypeLengths <- rep(runLength(GTsRle), runLength(GTsRle))
          HLs <- rep(as.integer(NA), length(snpsForThisSample))
          HLs[!NAGTs] <- haplotypeLengths
          HLs
        }
      )
    }
  )
  HLsMatrix <- do.call(rbind, HLsMatrixList)
  if(!identical(dimnames(HLsMatrix)[[2]], dimnames(geno(vcf)[["GT"]])[[2]])) {
    stop("setHaplotypeLengths: somehow the column names got changed")
  }
  dimnames(HLsMatrix) <- dimnames(geno(vcf)[["GT"]])
  geno(vcf)[["HL"]] <- HLsMatrix
  exptData(vcf)[["header"]] <- VCFHeader(
    reference=reference(exptData(vcf)[["header"]]),
    samples=samples(exptData(vcf)[["header"]]),
    header=DataFrameList(
      META=meta(exptData(vcf)[["header"]]),
      FILTER=fixed(exptData(vcf)[["header"]])[["FILTER"]],
      FORMAT=rbind(
        geno(exptData(vcf)[["header"]]),
        DataFrame(Number="1", Type="Integer", Description="Length of haplotype in terms of number of SNPs", row.names="HL")
      ),
      INFO=info(exptData(vcf)[["header"]])
    )
  )
  if(createInfoFieldSummaries) {
    numSingleSNPhaplotypes <- apply(HLsMatrix, 1, function(x) length(which(x==1)))
    singleSNPhaplotype <- numSingleSNPhaplotypes > 0
    numShortHaplotypes <- apply(HLsMatrix, 1, function(x) length(which(x<maxSNPsInShortHaplotype)))
    shortHaplotype <- numShortHaplotypes > 0
    info(vcf) <- cbind(
      values(info(vcf)),
      DataFrame(
        numSingleSNPhaplotypes = numSingleSNPhaplotypes,
        singleSNPhaplotype = singleSNPhaplotype,
        numShortHaplotypes = numShortHaplotypes,
        shortHaplotype = shortHaplotype
      )
    )
    shortHaplotypeDescription <- sprintf("Number of samples that have a short haplotype (max %d SNPs) at this position", maxSNPsInShortHaplotype)
    exptData(vcf)[["header"]] <- VCFHeader(
      reference=reference(exptData(vcf)[["header"]]),
      samples=samples(exptData(vcf)[["header"]]),
      header=DataFrameList(
        META=meta(exptData(vcf)[["header"]]),
        FILTER=fixed(exptData(vcf)[["header"]])[["FILTER"]],
  #      fixed=fixed(exptData(vcf)[["header"]]),
        FORMAT=geno(exptData(vcf)[["header"]]),
        INFO=rbind(
          info(exptData(vcf)[["header"]]),
          DataFrame(Number="1", Type="Integer", Description="Number of samples that have a single SNP haplotype at this position", row.names="numSingleSNPhaplotypes"),
          DataFrame(Number="0", Type="Flag", Description="Do any samples have a single SNP haplotype at this position", row.names="singleSNPhaplotype"),
          DataFrame(Number="1", Type="Integer", Description=shortHaplotypeDescription, row.names="numShortHaplotypes"),
          DataFrame(Number="0", Type="Flag", Description="Do any samples have a short haplotype at this position", row.names="shortHaplotype")
        )
      )
    )
  }
  return(vcf)
}

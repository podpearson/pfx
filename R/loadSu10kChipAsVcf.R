# loadSu10kChipAsVcf.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


loadSu10kChipAsVcf <- function(
  su10kChipFilename           = "data/Su10kChip/10K_Genotypes_HJiang_80%MinAssayCallRate.txt",
  su10kChipVcfFilename        = "data/Su10kChip/10K_Genotypes_HJiang_80%MinAssayCallRate.v3.vcf",
  su10kChipRdaFilename        = "data/Su10kChip/10K_Genotypes_HJiang_80%MinAssayCallRate.v3.vcf.rda",
  su10kChipV2BedFilename      = "data/Su10kChip/10K_Genotypes_HJiang_80%MinAssayCallRate.bed",
  su10kChipV3BedFilename      = "data/Su10kChip/10K_Genotypes_HJiang_80%MinAssayCallRate.v3.bed",
  su10kChipV3UnmappedFilename = "data/Su10kChip/10K_Genotypes_HJiang_80%MinAssayCallRate.v3.unmapped",
  v2tov3chainFile             = "opt/liftover/2to3.liftOver",
  liftoverExecutable          = "opt/liftover/liftOver",
  shouldSaveVcfFile           = TRUE, # Note that due to a bug in VariantAnnotation::.makeVcfGeno it is necessary to add "DUMMY" genotypes data
  shouldSaveRdaFile           = TRUE,
  reload                      = FALSE,
  verbose                     = TRUE
) {
  if(reload || !file.exists(su10kChipRdaFilename)) {
    sampleDetails <- read.table(su10kChipFilename, header=FALSE, sep="\t", as.is=TRUE, nrows=6)
    sampleManifest <- as.data.frame(t(sampleDetails)[13:51, ], stringsAsFactors=FALSE)
    names(sampleManifest) <- sampleDetails[, 12]
    row.names(sampleManifest) <- paste(sampleManifest[["Sample Name "]], sampleManifest[["Experiment Name"]], sep="_")
    
    su10kChipV2calls <- read.table(su10kChipFilename, header=TRUE, sep="\t", as.is=TRUE, skip=6)
    row.names(su10kChipV2calls) <- su10kChipV2calls[["External.Id"]]
    genotypes <- as.matrix(su10kChipV2calls[, 13:51])
    dimnames(genotypes) <- list(su10kChipV2calls[["External.Id"]], row.names(sampleManifest))
    su10kChipV2calls[["RefGenotype"]] <- ifelse(genotypes[, "3D7_(a)4095667-83237"] == genotypes[, "3D7_(a)4095667-82282"] & genotypes[, "3D7_(a)4095667-83237"] %in% c("A", "C", "T", "G"), genotypes[, "3D7_(a)4095667-83237"], NA)
    goodSNPs <- !is.na(su10kChipV2calls[["RefGenotype"]])
    su10kChipV2calls <- su10kChipV2calls[goodSNPs, ]
    sampleIDs <- row.names(sampleManifest[sampleManifest[["Experiment Call Rate %"]] > 0, ])
    genotypes <- genotypes[goodSNPs, sampleIDs]
    
    REFs <- matrix(rep(su10kChipV2calls[["RefGenotype"]], dim(genotypes)[2]), ncol=dim(genotypes)[2], dimnames=dimnames(genotypes))
    ALTs <- matrix(
      rep(
        ifelse(su10kChipV2calls[["Right"]]==su10kChipV2calls[["RefGenotype"]], su10kChipV2calls[["Left"]], su10kChipV2calls[["Right"]]),
        dim(genotypes)[2]
      ),
      ncol=dim(genotypes)[2],
      dimnames=dimnames(genotypes)
    )
    GT <- ifelse(genotypes==REFs, "0", ifelse(genotypes==ALTs, "1", "."))
  
    su10kChipV2Bed <- data.frame(
      Chr = sprintf("Pf3D7_%02d", su10kChipV2calls[["Chrom.Name"]]),
      Start = su10kChipV2calls[["Chrom.Position"]]-1,
      End = su10kChipV2calls[["Chrom.Position"]],
      ID = su10kChipV2calls[["External.Id"]],
      stringsAsFactors = FALSE
    )
    row.names(su10kChipV2Bed) <- su10kChipV2Bed[["ID"]]
    write.table(su10kChipV2Bed, file=su10kChipV2BedFilename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    system(paste(liftoverExecutable, su10kChipV2BedFilename, v2tov3chainFile, su10kChipV3BedFilename, su10kChipV3UnmappedFilename))  
    su10kChipV3Bed <- read.table(su10kChipV3BedFilename, header=FALSE, sep="\t", col.names=c("Chr", "Start", "End", "ID"), as.is=TRUE)
    row.names(su10kChipV3Bed) <- su10kChipV3Bed[["ID"]]
  
  #  couple of sanity checks. Firstly are all IDs still the same (they should be), and secondly have any variants changed chromosomes (only a few moved from 8 to 7)
    identical(su10kChipV2Bed[row.names(su10kChipV3Bed), "ID"], su10kChipV3Bed[["ID"]])
    table(su10kChipV2Bed[row.names(su10kChipV3Bed), "Chr"], su10kChipV3Bed[["Chr"]])
    
    su10kChipRowData  = GRanges(
      seqnames = su10kChipV3Bed[["Chr"]],
      ranges   = IRanges(
        start = su10kChipV3Bed[["End"]],
        width = 1,
        names  = su10kChipV3Bed[["ID"]]
      )
    )
      
    su10kChipVcf <- VCF(
      rowData = su10kChipRowData,
      colData = DataFrame(
        Samples = seq(along=sampleIDs),
        row.names = sampleIDs
      ),
      exptData = SimpleList(
        header = VCFHeader(
          samples = sampleIDs,
          header  = DataFrameList(
            META = rbind(
              DataFrame(Value = "VCFv4.0", row.names="fileformat"),
              DataFrame(Value = "Su10kChip", row.names="ProjectName"),
              DataFrame(Value = "7G8_(a)4095667-82251", row.names="PARENT"),
              DataFrame(Value = "GB4_(a)4095667-82103", row.names="PARENT.1")
            ),
  #          FILTER = DataFrame(Descrption=character()),
            FILTER = DataFrame(Descrption="PASS", row.names="PASS"),
            FORMAT = rbind(
              DataFrame(Number = "1", Type="String", Description="Genotype", row.names="GT")
            ),
            INFO = rbind(
              DataFrame(Number = "1", Type = "Integer", Description="Assay.Index", row.names="Assay.Index"),
              DataFrame(Number = "1", Type = "Integer", Description="Assay.Id", row.names="Assay.Id"),
              DataFrame(Number = "1", Type = "Float", Description="Assay.Call.Rate..", row.names="Assay.Call.Rate..")
            )
          )
        )
      ),
      fixed    = DataFrame(
        REF = DNAStringSet(REFs[row.names(su10kChipV3Bed), 1]),
        ALT = DNAStringSetList(ALTs[row.names(su10kChipV3Bed), 1]),
        QUAL = rep(0.0, dim(su10kChipV3Bed)[1]),
        FILTER = rep("PASS", dim(su10kChipV3Bed)[1])
      ),
      info = DataFrame(
        su10kChipV2calls[row.names(su10kChipV3Bed), c(1, 2, 11), ]
      ),
      geno     = SimpleList(
        GT = GT[row.names(su10kChipV3Bed), ]
      )
    )
  #  genome(su10kChipVcf) <- "Pf"
    su10kChipVcf <- su10kChipVcf[order(rowData(su10kChipVcf))]
    if(shouldSaveRdaFile) {
      save(su10kChipVcf, file=su10kChipRdaFilename)
    }
  } else {
    load(su10kChipRdaFilename)
  }
  if(shouldSaveVcfFile) {
    writeVcf(su10kChipVcf, su10kChipVcfFilename, index=TRUE)
#    writeVcf(su10kChipVcf, su10kChipVcfFilename)
  }
  return(su10kChipVcf)
}




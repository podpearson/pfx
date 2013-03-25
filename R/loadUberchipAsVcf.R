# loadUberchipAsVcf.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


loadUberchipAsVcf <- function(
  uberchipDirectory           = "data/uberchip/SNPchip_Data/SNP Calls",
  uberchipVcfFilename         = "data/uberchip/SNPchip_Data/SNP Calls/uberchip.v3.vcf",
  uberchipRdaFilename         = "data/uberchip/SNPchip_Data/SNP Calls/uberchip.v3.vcf.rda",
  uberchipV2BedFilename       = "data/uberchip/SNPchip_Data/uberchip.v2.bed",
#  uberchipV2FixedBedFilename  = "data/uberchip_et_al/gb-2011-12-4-r33-s3.fixed.bed",
  uberchipV3BedFilename       = "data/uberchip/SNPchip_Data/uberchip.v3.bed",
  uberchipV3UnmappedFilename  = "data/uberchip/SNPchip_Data/uberchip.v3.unmapped",
  liftoverExecutable          = "opt/liftover/liftOver",
  v2tov3chainFile             = "opt/liftover/2to3.liftOver",
  shouldSaveVcfFile           = FALSE,
  shouldSaveRdaFile           = TRUE,
  reload                      = FALSE,
  verbose                     = TRUE
) {
  if(reload || !file.exists(uberchipRdaFilename)) {
    snpCallFilenames <- dir(uberchipDirectory, pattern="txt$")
    uberchipV2calls <- sapply(
      snpCallFilenames,
      function(snpCallFilename) {
        if(verbose) cat(snpCallFilename, "\n")
        read.delim(file.path(uberchipDirectory, snpCallFilename))
      },
      USE.NAMES=TRUE,
      simplify=FALSE
    )
    sampleIDs <- sub("_SNPCalled.txt", "", names(uberchipV2calls))
    Chromosomes <- do.call(
      cbind,
      lapply(uberchipV2calls, function(x) x[["Chromosome"]])
    )
    Positions <- do.call(
      cbind,
      lapply(uberchipV2calls, function(x) x[["Position"]])
    )
    REFs <- do.call(
      cbind,
      lapply(uberchipV2calls, function(x) as.character(x[["X3D7.Allele"]]))
    )
    row.names(REFs) <- sprintf("uberchipSNP_%05d", seq(along=Positions[, 1]))
    ALTs <- do.call(
      cbind,
      lapply(uberchipV2calls, function(x) as.character(x[["Base.Call"]]))
    )
    table(apply(REFs, 1, function(x) length(unique(x[!is.na(x)]))))
    table(apply(ALTs, 1, function(x) length(unique(x[!is.na(x)]))))
    table(apply(ALTs, 1, function(x) length(unique(x[!is.na(x)]))))
    ALTsNoRef <- ALTs
    ALTsNoRef[ALTsNoRef==REFs] <- NA
    ALT <- apply(
      ALTsNoRef,
      1,
      function(x) {
        sortedALTs <- names(sort(table(x, useNA="ifany")))
        sortedALTs <- sortedALTs[!is.na(sortedALTs)]
  #      if(is.na(sortedALTs)) sortedALTs<- character(0)
        DNAStringSetList(sortedALTs)
      }
    )
  #  ALT <- apply(
  #    ALTs,
  #    1,
  #    function(x) {
  #      sortedALTs <- names(sort(table(x)))
  #      if(is.null(sortedALTs)) sortedALTs <- character(0)
  #      DNAStringSet(sortedALTs)
  #    }
  #  )
    ALTSSL <- do.call(c, ALT)
    table(elementLengths(ALTSSL))
    names(ALTSSL) <- sprintf("uberchipSNP_%05d", seq(along=Positions[, 1]))
  #  ALT_SSL <- do.call(c, ALT)
  #  ALT_SSL <- DNAStringSetList(ALT)
    
    GTs <- mapply(
      function(ref, alt, alleles) {
        match(alleles, c(ref, as.character(unlist(alt))))-1
      },
      REFs[, 1],
      ALT,
      split(ALTs, row(ALTs))
    )
    DSP <- do.call(
      cbind,
      lapply(uberchipV2calls, function(x) round(x[["Positive.Strand.Disc..Score"]], 2))
    )
    row.names(DSP) <- sprintf("uberchipSNP_%05d", seq(along=Positions[, 1]))
    DSN <- do.call(
      cbind,
      lapply(uberchipV2calls, function(x) round(x[["Negative.Strand.Disc..Score"]], 2))
    )
    row.names(DSN) <- sprintf("uberchipSNP_%05d", seq(along=Positions[, 1]))
    DS <- (DSP+DSN)/2
  
    uberchipV2Bed <- data.frame(
      Chr = sprintf("Pf3D7_%02d", Chromosomes[, 1]),
      Start = Positions[, 1]-1,
      End = Positions[, 1],
      ID = sprintf("uberchipSNP_%05d", seq(along=Positions[, 1])),
      stringsAsFactors = FALSE
    )
    row.names(uberchipV2Bed) <- uberchipV2Bed[["ID"]]
    write.table(uberchipV2Bed, file=uberchipV2BedFilename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    system(paste(liftoverExecutable, uberchipV2BedFilename, v2tov3chainFile, uberchipV3BedFilename, uberchipV3UnmappedFilename))  
    uberchipV3Bed <- read.table(uberchipV3BedFilename, header=FALSE, sep="\t", col.names=c("Chr", "Start", "End", "ID"), as.is=TRUE)
    row.names(uberchipV3Bed) <- uberchipV3Bed[["ID"]]
  
  #  couple of sanity checks. Firstly are all IDs still the same (they should be), and secondly have any variants changed chromosomes (only a few moved from 8 to 7)
    identical(uberchipV2Bed[row.names(uberchipV3Bed), "ID"], uberchipV3Bed[["ID"]])
    table(uberchipV2Bed[row.names(uberchipV3Bed), "Chr"], uberchipV3Bed[["Chr"]])
    
    uberchipRowData  = GRanges(
      seqnames = uberchipV3Bed[["Chr"]],
      ranges   = IRanges(
        start = uberchipV3Bed[["End"]],
        width = 1,
        names  = uberchipV3Bed[["ID"]]
      )
    )
    
    GTsText <- t(matrix(as.character(GTs), ncol=ncol(GTs), dimnames=dimnames(GTs)))
    GTsText[is.na(GTsText)] <- "."
    row.names(GTsText) <- sprintf("uberchipSNP_%05d", seq(along=Positions[, 1]))
  
    uberchipVcf <- VCF(
      rowData = uberchipRowData,
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
              DataFrame(Value = "Uberchip", row.names="ProjectName"),
              DataFrame(Value = "HB3_1", row.names="PARENT"),
              DataFrame(Value = "Dd2_1", row.names="PARENT.1")
            ),
  #          FILTER = DataFrame(Descrption=character()),
            FILTER = DataFrame(Descrption="PASS", row.names="PASS"),
            FORMAT = rbind(
              DataFrame(Number = "1", Type="String", Description="Genotype", row.names="GT"),
              DataFrame(Number = "1", Type="Float", Description="Dscore (discrimination score) on the positive strand", row.names="DSP"),
              DataFrame(Number = "1", Type="Float", Description="Dscore (discrimination score) on the negative strand", row.names="DSN"),
              DataFrame(Number = "1", Type="Float", Description="Mean Dscore (discrimination score) between strands", row.names="DS")
            ),
            INFO = rbind(
              DataFrame(Number = "1", Type = "String", Description="Dummy variable, please ignire", row.names="Dummy")
            )
          )
        )
      ),
      fixed    = DataFrame(
        REF = DNAStringSet(REFs[row.names(uberchipV3Bed), 1]),
        ALT = ALTSSL[row.names(uberchipV3Bed)],
        QUAL = rep(0.0, dim(uberchipV3Bed)[1]),
        FILTER = rep("PASS", dim(uberchipV3Bed)[1])
      ),
      info     = DataFrame(
        Dummy = rep("dummy", dim(uberchipV3Bed)[1])
      ),
      geno     = SimpleList(
        GT = GTsText[row.names(uberchipV3Bed), ],
        DSP = DSP[row.names(uberchipV3Bed), ],
        DSN = DSN[row.names(uberchipV3Bed), ],
        DS = DS[row.names(uberchipV3Bed), ]
      )
    )
  #  genome(uberchipVcf) <- "Pf"
    uberchipVcf <- uberchipVcf[order(rowData(uberchipVcf))]
    if(shouldSaveRdaFile) {
      save(uberchipVcf, file=uberchipRdaFilename)
    }
  } else {
    load(uberchipRdaFilename)
  }
  if(shouldSaveVcfFile) {
    writeVcf(uberchipVcf, uberchipVcfFilename, index=TRUE)
#    writeVcf(uberchipVcf, uberchipVcfFilename)
  }
  return(uberchipVcf)
}




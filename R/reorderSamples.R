# reorderSamples.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


reorderSamples <- function(
  GTsInt,
  parentalIDs,
  sampleIDmappings
) {
  parentalSampleIndexes <- which(
    sampleIDmappings %in% sampleIDmappings[
      which(dimnames(GTsInt)[[2]] %in% parentalIDs)
    ]
  )
  parentalSampleGTsInt <- GTsInt[
    ,
    parentalSampleIndexes
  ][
    ,
    order(sampleIDmappings[parentalSampleIndexes])
  ]
  dimnames(parentalSampleGTsInt)[[2]] <- names(sampleIDmappings)[
    parentalSampleIndexes
  ][
    order(sampleIDmappings[parentalSampleIndexes])
  ]
  progenySampleGTsInt <- GTsInt[
    ,
    -(parentalSampleIndexes)
  ][
    ,
    order(sampleIDmappings[-(parentalSampleIndexes)])
  ]
  dimnames(progenySampleGTsInt)[[2]] <- names(sampleIDmappings)[
    -(parentalSampleIndexes)
  ][
    order(sampleIDmappings[-(parentalSampleIndexes)])
  ]
  reorderedGTsInt <- cbind(parentalSampleGTsInt, progenySampleGTsInt)
  uniqueSampleIDs <- c(
    sampleIDmappings[parentalSampleIndexes][
      order(sampleIDmappings[parentalSampleIndexes])
    ],
    sampleIDmappings[-(parentalSampleIndexes)][
      order(sampleIDmappings[-(parentalSampleIndexes)])
    ]
  )
  linePositions <- which(
    as.integer(
      factor(uniqueSampleIDs)
    )!=c(
      as.integer(factor(uniqueSampleIDs))[-1],
      NA
    )
  )
  
  return(
    list(
      GTsInt=reorderedGTsInt,
      linePositions=linePositions,
      parentalIDs=names(
        sampleIDmappings[which(dimnames(GTsInt)[[2]] %in% parentalIDs)]
      )
    )
  )
}

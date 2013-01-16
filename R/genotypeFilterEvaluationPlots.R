# genotypeFilterEvaluationPlots.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################

#  lapply(c("final", "bestReplicate", "uncontaminated"), function(x) genotypeFilterEvaluationPlots(sampleSet="x", infoFilterString = paste(x, "0.2.33.1.2.0.01.0.1.0.2.0.9.38.4.600", sep=".")))
#  lapply(c("final", "bestReplicate", "uncontaminated"), function(x) genotypeFilterEvaluationPlots(sampleSet="x", infoFilterString = x))
#  lapply(c("final", "bestReplicate", "uncontaminated"), function(x) genotypeFilterEvaluationPlots(sampleSet="x", variablesToPlot = list("DPforGQ99andMAF0.1" = c("15"=1, "10"=2, "8"=3, "6"=4, "5"=5, "4"=6, "3"=7, "2"=8, "1"=9, "0"=10)), infoFilterString = paste(x, "0.2.33.1.2.0.01.0.1.0.2.0.9.38.4.600", sep=".")))
#  lapply(c("final", "bestReplicate", "uncontaminated"), function(x) genotypeFilterEvaluationPlots(sampleSet="x", variablesToPlot = list("DPforGQ99andMAF0.1" = c("15"=1, "10"=2, "8"=3, "6"=4, "5"=5, "4"=6, "3"=7, "2"=8, "1"=9, "0"=10)), infoFilterString = x))

genotypeFilterEvaluationPlots <- function(
  analysisDirectory           = "/data/malariagen2/plasmodium/pf-crosses/data/3d7_v3/bwa_n0.01_k4_l32/genotypes_analysis_20120107/gatk",
  crosses                     = c("3d7_hb3", "7g8_gb4", "hb3_dd2"),
  variantTypes                = c("snps", "indels"),
#  crosses                     = dir(analysisDirectory),
#  variantTypes                = rev(dir(file.path(analysisDirectory, crosses[1]))), # rev just because I want snps before indels
  sampleSet                   = c("final", "bestReplicate", "uncontaminated"),
#  sampleSet                   = c("FinalSamples", "BestReplicate", "Uncontaminated"),
  infoFilterString            = sampleSet[1],
#  infoFilterString            = paste(sampleSet, "0.2.33.1.2.0.01.0.1.0.2.0.9.38.4.600", sep="."),
  variablesToPlot             = list("MAFforGQ99andDP10" = c("0" = 2, "0.02"=3, "0.05"=4, "0.1"=1, "0.2"=5, "0.35"=6, "0.5"=7)),
  metricsToExclude            = c(
    "haplotypeParent1", "haplotypeParent2", "cross", "variantType", "value",
    "totalHaplotypesInProgeny", "medianHaplotypesPerProgeny", "minHaplotypesPerProgeny", "maxHaplotypesPerProgeny",
    "totalRecombinations", "whichMinRecombinationsPerSample", "whichMaxRecombinationsPerSample",
    "whichMinSingleSNPhaplotypesPerSample", "whichMaxSingleSNPhaplotypesPerSample",
    "totalPutativeCrossovers", "whichMinPutativeCrossoversPerSample", "whichMaxPutativeCrossoversPerSample",
    "numberOfVariants"
  ),
  height                      = 20,
  width                       = height*(16.53/11.69)
) {
  require(ggplot2)
  require(gridExtra)
  require(RColorBrewer)
  crossesVector <- rep(crosses, times=length(variantTypes))
  variantTypesVector <- rep(variantTypes, each=length(crosses))
  crossesByVariantTypes <- paste(
    rep(crosses, times=length(variantTypes)),
    rep(variantTypes, each=length(crosses)),
    sep="."
  )
  fullReturnDFfilenames <- file.path(analysisDirectory, paste(crossesByVariantTypes, "evaluateGenotypeFilters", infoFilterString, "fullReturnDF.rda", sep="."))
  allResultsList <- sapply(
    names(variablesToPlot),
    function(variableToPlot) {
      do.call(
        rbind,
        lapply(
          crossesByVariantTypes,
          function(crossVariantType) {
            cross <- strsplit(crossVariantType, "\\.")[[1]][1]
            variantType <- strsplit(crossVariantType, "\\.")[[1]][2]
            load(file.path(analysisDirectory, cross, variantType, paste(crossVariantType, "evaluateGenotypeFilters", infoFilterString, "fullReturnDF.rda", sep=".")))
            fullReturnDF[["cross"]] <- cross
            fullReturnDF[["variantType"]] <- variantType
            fullReturnDF[["crossVariantType"]] <- crossVariantType
            fullReturnDF[variablesToPlot[[variableToPlot]], "value"] <- factor(names(variablesToPlot[[variableToPlot]]), levels=names(variablesToPlot[[variableToPlot]]))
            fullReturnDF[variablesToPlot[[variableToPlot]], ]
          }
        )
      )
    },
    simplify=FALSE,
    USE.NAMES=TRUE
  )
  sapply(
    names(allResultsList),
    function(variableToPlot) {
      plotsToCreate <- setdiff(names(allResultsList[[variableToPlot]]), metricsToExclude)
      pdf(file.path(analysisDirectory, paste("genotypeFilterEvaluationPlots", variableToPlot, infoFilterString, "pdf", sep=".")), height=height, width=width)
      plots <- lapply(
        plotsToCreate,
        function(plotToCreate) {
          cat(".")
          qplotStatement <- paste(
            "qplot(value,",
            plotToCreate,
            ", data=allResultsList[[variableToPlot]], geom=\"line\", group=crossVariantType, colour=crossVariantType) + theme_bw() + theme(legend.position=\"none\") + scale_colour_brewer(palette=\"Paired\")"
          )
          eval(parse(text=qplotStatement))
#          qplot(
#            Quantile,
#            ProportionOfMendelianErrors,
#    #        fill=Annotation,
#    #        facets=Annotation~.,
#            data=plotDFquantiles[[plotDFquantileIndex]],
#            geom="bar",
#            stat="identity",
#            fill=I(col[plotDFquantileIndex]),
#    #        fill=I(brewer.pal(12, "Set3")[plotDFquantileIndex]),
#            main=names(plotDFquantiles)[plotDFquantileIndex],
#            ylim = c(0, maxHeight),
#            ylab = paste("Proportion of", errorVariable, ">", errorThreshold)
#          ) +
#    #      scale_fill_brewer(palette="Set3")[plotDFquantileIndex] +
#          theme_bw() +
#          theme(legend.position = c(0.5, 1)) + 
#          theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#          theme(axis.title.x = element_blank())
        }
      )
      print(
        do.call(
          grid.arrange,  c(plots, ncol=ceiling(sqrt(length(plots))))
        )
      )
      dev.off()
    }
  )
}

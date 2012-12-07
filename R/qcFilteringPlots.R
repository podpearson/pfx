# qcFilteringPlots.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


qcFilteringPlots <- function(
  vcf,
  plotFilestem                = "7G8xGB4",
  variablesToPlot             = c(
#    "AC"             = "highIsGood",
    "BaseQRankSum"   = "highIsGood",
#    "DS"             = "lowIsGood",
#    "Dels"           = "lowIsGood",
    "FS"             = "lowIsGood",
    "HaplotypeScore" = "lowIsGood",
    "MQ"             = "highIsGood",
    "MQ0"            = "lowIsGood",
    "MQRankSum"      = "highIsGood",
    "QD"             = "highIsGood",
#    "RPA"            = "lowIsGood",
    "ReadPosRankSum" = "highIsGood",
    "SB"             = "lowIsGood",
    "meanMAF"        = "lowIsGood",
#    "missingness"    = "lowIsGood",
    "missingness2"   = "lowIsGood",
    "heterozgosity"  = "lowIsGood"
  ),
  variablesToPlotQuantiles = c(
#    "AC"             = "highIsGood",
    "BaseQRankSum"   = "highIsGood",
#    "DS"             = "lowIsGood",
#    "Dels"           = "lowIsGood",
    "FS"             = "lowIsGood",
    "HaplotypeScore" = "lowIsGood",
    "MQ"             = "highIsGood",
    "MQ0"            = "lowIsGood",
    "MQRankSum"      = "highIsGood",
    "QD"             = "highIsGood",
#    "RPA"            = "lowIsGood",
    "ReadPosRankSum" = "highIsGood",
    "SB"             = "lowIsGood",
    "meanMAF"        = "lowIsGood",
#    "missingness"    = "lowIsGood"
    "missingness2"   = "lowIsGood",
    "heterozgosity"  = "lowIsGood"
  ),
  subsetToBiallelic           = FALSE,
#  subsetToBiallelic           = TRUE,
  regionsToMask               = NULL,
#  regionsToMask               = varRegions_v3(),
  numberOfQuantiles           = 50,
  verbose                     = TRUE
) {
  if(!is.null(regionsToMask)) {
    vcf <- vcf[!(rowData(vcf) %in% regionsToMask)]
  }
  if(subsetToBiallelic) {
    vcf <- vcf[elementLengths(alt(vcf)) == 1]
  }
  plotDF <- do.call(
    rbind,
    lapply(
      names(variablesToPlot),
      function(variableToPlot) {
        if(verbose) {
          cat(variableToPlot, "\n")
        }
        data.frame(
          Annotation     = variableToPlot,
          Log10ErrorRate = log10(
            cumsum(
              values(info(vcf))[["MendelianErrors"]][
                order(
                  values(info(vcf))[[variableToPlot]],
                  decreasing=variablesToPlot[variableToPlot]=="highIsGood"
                )
              ] > 0
            ) / seq(along=values(info(vcf))[["MendelianErrors"]] > 0)
          ),
          NumberOfSegregatingSites = seq(along=values(info(vcf))[["MendelianErrors"]] > 0)
        )
      }
    )
  )
  plotDFquantiles <- # do.call(
#    rbind,
    sapply(
      names(variablesToPlotQuantiles),
      function(variableToPlot) {
        if(verbose) {
          cat(variableToPlot, "\n")
        }
#        quantiles <- cut(values(info(vcf))[[variableToPlot]], quantile(values(info(vcf))[[variableToPlot]], probs=seq(0, 1, 1/100)))
#        quantiles <- ggplot2::cut_number(values(info(vcf))[[variableToPlot]], numberOfQuantiles)
        set.seed(12345)
        quantiles <- ceiling(rank(values(info(vcf))[[variableToPlot]], ties.method="random")/(dim(vcf)[1]/numberOfQuantiles))
        proportions <- by(values(info(vcf))[["MendelianErrors"]], quantiles, function(x) length(which(x>0))/length(x))
        lowerBounds <- by(values(info(vcf))[[variableToPlot]], quantiles, function(x) min(x, na.rm=TRUE))
        upperBounds <- by(values(info(vcf))[[variableToPlot]], quantiles, function(x) max(x, na.rm=TRUE))
        xaxisLabels <- paste(seq(along=lowerBounds), ": [", format(lowerBounds, digits=3), ", ", format(upperBounds, digits=3), "]", sep="")
        data.frame(
          Annotation                  = variableToPlot,
          Quantile                    = factor(xaxisLabels, levels=xaxisLabels),
#          Quantile                    = factor(as.character(seq(along=proportions)), levels=as.character(seq(along=proportions))),
          ProportionOfMendelianErrors = as.vector(proportions)
        )
      },
      simplify=FALSE,
      USE.NAMES=TRUE
    )
#  )
  require(ggplot2)
  require(gridExtra)
  require(RColorBrewer)
  pdf(paste(plotFilestem, "binnedErrorRates.pdf", sep="."), height=20, width=30)
  maxHeight <- max(sapply(plotDFquantiles, function(x) max(x[["ProportionOfMendelianErrors"]])))
  plots <- lapply(
    seq(along = plotDFquantiles),
#    names(plotDFquantiles),
    function(plotDFquantileIndex) {
      qplot(
        Quantile,
        ProportionOfMendelianErrors,
#        fill=Annotation,
#        facets=Annotation~.,
        data=plotDFquantiles[[plotDFquantileIndex]],
        geom="bar",
        stat="identity",
        fill=I(brewer.pal(12, "Set3")[plotDFquantileIndex]),
        main=names(plotDFquantiles)[plotDFquantileIndex],
        ylim = c(0, maxHeight)
      ) +
#      scale_fill_brewer(palette="Set3")[plotDFquantileIndex] +
      theme_bw() +
      theme(legend.position = c(0.5, 1)) + 
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
      theme(axis.title.x = element_blank())
    }
  )
  print(
    do.call(
      grid.arrange,  c(plots, ncol=4)
    )
  )
#  print(
#    qplot(
#      Quantile,
#      ProportionOfMendelianErrors,
#      fill=Annotation,
#      facets=Annotation~.,
#      data=plotDFquantiles,
#      geom="bar"
##      xlab="Quantile",
##      ylab="log10 (Mendelian/SingleSNPhaplotype error rate)",
##      ylim=c(-2,0)
#    )
##    + geom_bar()
#    + scale_fill_brewer(palette="Set3")
#    + theme_bw()
#  )
  dev.off()
  pdf(paste(plotFilestem, "log10ErrorRates.pdf", sep="."), height=6, width=10)
  print(
    qplot(
      NumberOfSegregatingSites,
      Log10ErrorRate,
      colour=Annotation,
      data=plotDF,
      xlab="# segregating sites",
      ylab="log10 (Mendelian/SingleSNPhaplotype error rate)",
      ylim=c(-4,0)
    )
    + scale_colour_brewer(palette="Set3")
    + theme_bw()
  )
  dev.off()
  
#  some debuggined stuff - ignore
#  head(subset(plotDF, Annotation=="HaplotypeScore"))
#  subset(plotDF, Annotation=="HaplotypeScore")[c(1, 10, 30, 50, 100, 200, 300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 30000, 50000, 100000, 200000, 251099), ]
#  plotDF[
#    c(
#      c(1, 10, 30, 50, 100, 200, 300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 30000, 50000, 100000, 200000, 251099),
#      251099+c(1, 10, 30, 50, 100, 200, 300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 30000, 50000, 100000, 200000, 251099)
#    ),
#  ]
  return(list(quantilesDF=plotDFquantiles, filteringDF=plotDF))
}

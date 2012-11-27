# recombinationPlot.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


recombinationPlot = function(
  genotypesMatrix,
  positions                   = as.integer(sub("^[^:]*:(.*)$", "\\1", dimnames(genotypesMatrix)[[1]])),
  sampleIDs                   = dimnames(genotypesMatrix)[[2]],
  col                         = c("white", "blue", "red", "lightblue", "pink", "lightgrey", "black"),
  breaks                      = c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5),
  main                        = NULL,
  thin                        = round(length(positions)/200),
  linePositions               = c(2),
  xaxisLabelQuantiles         = seq(0, 1, 0.1)
)
{
  layout(matrix(c(1,2),byrow = TRUE, ncol = 1),heights = c(4,2))
  num.samples = dim(genotypesMatrix)[2]
  par(mar=c(0,7,1,2))
  par(cex.lab=1.3)
  if(is.null(main)) {
    image(genotypesMatrix,col = col,breaks=breaks,axes=FALSE,ylab="Samples")
  } else {
    image(genotypesMatrix,col = col,breaks=breaks,axes=FALSE,ylab="Samples",main=main )
  }
  lineYcoords <- 1 - (linePositions/(length(sampleIDs)-1)) + (0.5/length(sampleIDs))
  sapply(
    lineYcoords,
    function(lineYcoord) {
      lines(x=c(0,1), y=c(lineYcoord, lineYcoord), ,col ="black")
    }
  )
  num.snps = length(positions)
  axis(2, labels = sampleIDs, at = seq(0.0,num.samples,length.out= num.samples)/(num.samples), cex.axis = 0.4, las = 2)
  box()
  max.val = max(positions,na.rm=TRUE);
  min.val = min(positions,na.rm=TRUE);
  par(mar=c(3,7,0,2))
  plot(1,col="white", type="n", xlab= "SNP Position", axes=FALSE,xaxs="i",yaxs="i",ylim=c(0,1),xlim=c(min.val,max.val),ylab="")
  box();
  for (k in seq(1,num.snps,by = thin))
  {
    lines(x = c(positions[k],min.val+(k/num.snps*(max.val-min.val))),y=c(0,1),col ="black")
  }
  axis(1,at = round(min.val+xaxisLabelQuantiles*(max.val-min.val)), label = round(min.val+xaxisLabelQuantiles*(max.val-min.val)))
}


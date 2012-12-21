# randomColours.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


randomColours <- Vectorize(
  function(inputInteger) {
    set.seed(inputInteger)
    rgb(runif(1), runif(1), runif(1))
  }
)


# pipelineFilterEvaluation.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


# pipelineFilterEvaluation(chromosomes = sprintf("Pf3D7_%02d_v3", 1), DPthresholds=5, shouldCreateQCFilteringPlots=TRUE)
# pipelineFilterEvaluation(DPthresholds=5, shouldUseSavedVersions=TRUE)

# final sets
#pipelineFilterEvaluation("3d7_hb3", "indels", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), sampleSets=list("FinalSamples" = list(additionalInfoFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), sampleSetName="final"), "BestReplicate" = list(additionalInfoFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), sampleSetName="bestReplicate")))


#pipelineFilterEvaluation("3d7_hb3", "snps")
#pipelineFilterEvaluation("7g8_gb4", "snps")
#pipelineFilterEvaluation("hb3_dd2", "snps")
#pipelineFilterEvaluation("3d7_hb3", "indels")
#pipelineFilterEvaluation("7g8_gb4", "indels")
#pipelineFilterEvaluation("hb3_dd2", "indels")

#pipelineFilterEvaluation("3d7_hb3", "snps", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)))
#pipelineFilterEvaluation("7g8_gb4", "snps", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)))
#pipelineFilterEvaluation("hb3_dd2", "snps", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)))
#pipelineFilterEvaluation("3d7_hb3", "indels", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)))
#pipelineFilterEvaluation("7g8_gb4", "indels", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), DPthresholds=10)
#pipelineFilterEvaluation("hb3_dd2", "indels", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)))

#pipelineFilterEvaluation("3d7_hb3", "snps", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2)))
#pipelineFilterEvaluation("7g8_gb4", "snps", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2)))
#pipelineFilterEvaluation("hb3_dd2", "snps", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2)))
#pipelineFilterEvaluation("3d7_hb3", "indels", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5, "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2))))
#pipelineFilterEvaluation("7g8_gb4", "indels", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2)))
#pipelineFilterEvaluation("hb3_dd2", "indels", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2)))
#
#pipelineFilterEvaluation("3d7_hb3", "snps", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41)))
#pipelineFilterEvaluation("7g8_gb4", "snps", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41)))
#pipelineFilterEvaluation("hb3_dd2", "snps", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41)))
#pipelineFilterEvaluation("3d7_hb3", "indels", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41)))
#pipelineFilterEvaluation("7g8_gb4", "indels", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41)))
#pipelineFilterEvaluation("hb3_dd2", "indels", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41)))
#
#pipelineFilterEvaluation("3d7_hb3", "snps", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("7g8_gb4", "snps", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("hb3_dd2", "snps", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("3d7_hb3", "indels", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("7g8_gb4", "indels", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("hb3_dd2", "indels", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#
#pipelineFilterEvaluation("3d7_hb3", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2)))
#pipelineFilterEvaluation("7g8_gb4", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2)))
#pipelineFilterEvaluation("hb3_dd2", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2)))
#pipelineFilterEvaluation("3d7_hb3", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6, "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2))))
#pipelineFilterEvaluation("7g8_gb4", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2)))
#pipelineFilterEvaluation("hb3_dd2", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2)))
#
## The following is the set I emailed Alistair about 2013/01/22
#pipelineFilterEvaluation("3d7_hb3", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("7g8_gb4", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("hb3_dd2", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("3d7_hb3", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("7g8_gb4", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("hb3_dd2", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#
#pipelineFilterEvaluation("3d7_hb3", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), sampleSets=list("FinalSamples" = list(additionalInfoFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), sampleSetName="final")))
#pipelineFilterEvaluation("7g8_gb4", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), sampleSets=list("FinalSamples" = list(additionalInfoFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), sampleSetName="final")))
#pipelineFilterEvaluation("hb3_dd2", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), sampleSets=list("FinalSamples" = list(additionalInfoFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), sampleSetName="final")))
#pipelineFilterEvaluation("3d7_hb3", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), sampleSets=list("FinalSamples" = list(additionalInfoFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), sampleSetName="final")))
#pipelineFilterEvaluation("7g8_gb4", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), sampleSets=list("FinalSamples" = list(additionalInfoFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), sampleSetName="final")))
#pipelineFilterEvaluation("hb3_dd2", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), sampleSets=list("FinalSamples" = list(additionalInfoFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), sampleSetName="final")))
#
#pipelineFilterEvaluation("3d7_hb3", "snps", DPthresholds=3, variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("7g8_gb4", "snps", DPthresholds=3, variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("hb3_dd2", "snps", DPthresholds=3, variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("3d7_hb3", "indels", DPthresholds=3, variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("7g8_gb4", "indels", DPthresholds=3, variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("hb3_dd2", "indels", DPthresholds=3, variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#
#pipelineFilterEvaluation("3d7_hb3", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), setMonomorphicProgenyFilter=FALSE)
#pipelineFilterEvaluation("7g8_gb4", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), setMonomorphicProgenyFilter=FALSE)
#pipelineFilterEvaluation("hb3_dd2", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), setMonomorphicProgenyFilter=FALSE)
#pipelineFilterEvaluation("3d7_hb3", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), setMonomorphicProgenyFilter=FALSE)
#pipelineFilterEvaluation("7g8_gb4", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), setMonomorphicProgenyFilter=FALSE)
#pipelineFilterEvaluation("hb3_dd2", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), setMonomorphicProgenyFilter=FALSE)
#
#
#vcf <- pipelineFilterEvaluation("7g8_gb4", "snps", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), shouldReturnVcfOnly=TRUE)
#vcf7 <- pipelineFilterEvaluation("7g8_gb4", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), shouldReturnVcfOnly=TRUE)[[1]][[1]]
#vcfhd <- pipelineFilterEvaluation("hb3_dd2", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), shouldReturnVcfOnly=TRUE)[[1]][[1]]
#
#vcf3 <- pipelineFilterEvaluation("3d7_hb3", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), shouldReturnVcfOnly=TRUE, regionsMask=NULL)[[1]][[1]]

## Rerun of set I emailed Alistair about 2013/01/22 using Alisatir's definitions of nonVar regions
#pipelineFilterEvaluation("3d7_hb3", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), shouldUseSavedVersions=FALSE)
#pipelineFilterEvaluation("7g8_gb4", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), shouldUseSavedVersions=FALSE)
#pipelineFilterEvaluation("hb3_dd2", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), shouldUseSavedVersions=FALSE)
#pipelineFilterEvaluation("3d7_hb3", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), shouldUseSavedVersions=FALSE)
#pipelineFilterEvaluation("7g8_gb4", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), shouldUseSavedVersions=FALSE)
#pipelineFilterEvaluation("hb3_dd2", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), shouldUseSavedVersions=FALSE)

## Remove MAF restriction on indels
#pipelineFilterEvaluation("3d7_hb3", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), shouldUseSavedVersions=FALSE)
#pipelineFilterEvaluation("7g8_gb4", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), shouldUseSavedVersions=FALSE)
#pipelineFilterEvaluation("hb3_dd2", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), shouldUseSavedVersions=FALSE)
#pipelineFilterEvaluation("3d7_hb3", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), genotypeFiltersList=list(list("LowGQ" = list(column="GQ", operator="<", value=99, filterOutNAs=TRUE), "LowDP" = list(column="DP", operator="<", value=5, filterOutNAs=TRUE), "HighMAF" = list(column="MAF", operator=">", value=Inf, filterOutNAs=FALSE))), overwriteExisting=TRUE, shouldUseSavedVersions=FALSE)
#pipelineFilterEvaluation("7g8_gb4", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), genotypeFiltersList=list(list("LowGQ" = list(column="GQ", operator="<", value=99, filterOutNAs=TRUE), "LowDP" = list(column="DP", operator="<", value=5, filterOutNAs=TRUE), "HighMAF" = list(column="MAF", operator=">", value=Inf, filterOutNAs=FALSE))), overwriteExisting=TRUE, shouldUseSavedVersions=FALSE)
#pipelineFilterEvaluation("hb3_dd2", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), genotypeFiltersList=list(list("LowGQ" = list(column="GQ", operator="<", value=99, filterOutNAs=TRUE), "LowDP" = list(column="DP", operator="<", value=5, filterOutNAs=TRUE), "HighMAF" = list(column="MAF", operator=">", value=Inf, filterOutNAs=FALSE))), overwriteExisting=TRUE, shouldUseSavedVersions=FALSE)

## Rerun of set I emailed Alistair about 2013/01/22 using Alisatir's definitions of nonVar regions and DP rather than AD for depth calculations
#pipelineFilterEvaluation("3d7_hb3", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), overwriteExisting=TRUE, shouldUseSavedVersions=FALSE)
#pipelineFilterEvaluation("7g8_gb4", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), overwriteExisting=TRUE, shouldUseSavedVersions=FALSE)
#pipelineFilterEvaluation("hb3_dd2", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), overwriteExisting=TRUE, shouldUseSavedVersions=FALSE)
#pipelineFilterEvaluation("3d7_hb3", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), overwriteExisting=TRUE, shouldUseSavedVersions=FALSE)
#pipelineFilterEvaluation("7g8_gb4", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), overwriteExisting=TRUE, shouldUseSavedVersions=FALSE)
#pipelineFilterEvaluation("hb3_dd2", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), overwriteExisting=TRUE, shouldUseSavedVersions=FALSE)

## The following is the set I emailed Alistair about 2013/01/22
#pipelineFilterEvaluation("3d7_hb3", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("7g8_gb4", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("hb3_dd2", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("3d7_hb3", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("7g8_gb4", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("hb3_dd2", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))

## The following is the set I emailed Alistair about 2013/01/22, but without filtering out AD=0,0 variants
#pipelineFilterEvaluation("3d7_hb3", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), filterOutMAFNAs=FALSE)
#pipelineFilterEvaluation("7g8_gb4", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), filterOutMAFNAs=FALSE)
#pipelineFilterEvaluation("hb3_dd2", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), filterOutMAFNAs=FALSE)
#pipelineFilterEvaluation("3d7_hb3", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), filterOutMAFNAs=FALSE)
#pipelineFilterEvaluation("7g8_gb4", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), filterOutMAFNAs=FALSE)
#pipelineFilterEvaluation("hb3_dd2", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), filterOutMAFNAs=FALSE)

## The following is the set I emailed Alistair about 2013/01/22
#pipelineFilterEvaluation("3d7_hb3", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("7g8_gb4", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("hb3_dd2", "snps", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("3d7_hb3", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("7g8_gb4", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))
#pipelineFilterEvaluation("hb3_dd2", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)))


#vcf_3d7_hbs_indels <- pipelineFilterEvaluation("3d7_hb3", "indels", variantFilters=list("QUALbyDP6" = list(column="QUALbyDP", operator="<", value=6), "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2), "UQ41" = list(column="UQ", operator=">", value=41), "DPSoftClippedFraction0.25" = list(column="DPSoftClippedFraction", operator=">", value=0.25)), genotypeFiltersList=list(list("LowGQ" = list(column="GQ", operator="<", value=99, filterOutNAs=TRUE), "LowDP" = list(column="DP", operator="<", value=5, filterOutNAs=TRUE), "HighMAF" = list(column="MAF", operator=">", value=Inf, filterOutNAs=FALSE))), shouldReturnVcfOnly=TRUE)

#pipelineFilterEvaluation("7g8_gb4", "indels", variantFilters=list("QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5)), DPthresholds=10)
pipelineFilterEvaluation <- function(
  cross                       = "3d7_hb3",
  variantType                 = "snps",
  genotypesDirectory          = "data/3d7_v3/bwa_n0.01_k4_l32/genotypes/per_sample_realigned/gatk_2.3.9",
  chromosomes                 = sprintf("Pf3D7_%02d_v3", 1:14),
  outputDirectory             = "data/3d7_v3/bwa_n0.01_k4_l32/genotypes_analysis_20130129/per_sample_realigned/gatk_2.3.9",
  genotypesFileFmt            = "%s.annotated.vcf",
  initialSampleQCresultsFile  = file.path("data/3d7_v3/bwa_n0.01_k4_l32/genotypes_analysis_20130121/per_sample_realigned/gatk", cross, variantType, paste(cross, ".initialSampleQCresults.rda", sep="")),
  variantFilters              = NULL,
#  variantFilters=list(
#    "QUALbyDP5" = list(column="QUALbyDP", operator="<", value=5),
#    "QD14" = list(column="QD", operator="<", value=14),
#    "BaseQRankSum45" = list(column="BaseQRankSum", operator="<", value=-45, filterOutNAs=FALSE),
#    "numFilteredGenotypes10" = list(column="numFilteredGenotypes", operator=">", value=10)
#  ),
  sampleSets                  = list(
#    "FinalSamples"              = list(additionalInfoFilters=variantFilters, sampleSetName="final"),
    "BestReplicate"             = list(additionalInfoFilters=variantFilters, sampleSetName="bestReplicate")
  ),
  parentalStrains             = if(cross=="3d7_hb3") {
    c("ERR019061", "ERR019054")
  } else if(cross=="7g8_gb4") {
    c("ERR027099", "ERR027100")
  } else if(cross=="hb3_dd2") {
    c("ERR012788", "ERR012840")
  },
  GQthresholds                = NULL,
#  GQthresholds                = c(99, 50, 5),
#  DPthresholds                = NULL,
  DPthresholds                = c(5),
#  DPthresholds                = c(10, 5, 1),
  MAFthresholds               = NULL,
#  MAFthresholds               = c(0.1, 0, 0.02, 0.05, 0.2, 0.35, 0.5),
  GQthresholdDefault          = 99,
  DPthresholdDefault          = 5,
  MAFthresholdDefault         = 0.1,
  filterOutMAFNAs             = TRUE,
  genotypeFiltersList         = c(
    lapply(
      MAFthresholds,
      function(MAFthreshold) {
        list(
          "LowGQ" = list(column="GQ", operator="<", value=GQthresholdDefault, filterOutNAs=TRUE),
          "LowDP" = list(column="DP", operator="<", value=DPthresholdDefault, filterOutNAs=TRUE),
          "HighMAF" = list(column="MAF", operator=">", value=MAFthreshold, filterOutNAs=filterOutMAFNAs)
        )
      }
    ),
    lapply(
      DPthresholds,
      function(DPthreshold) {
        list(
          "LowGQ" = list(column="GQ", operator="<", value=GQthresholdDefault, filterOutNAs=TRUE),
          "LowDP" = list(column="DP", operator="<", value=DPthreshold, filterOutNAs=TRUE),
          "HighMAF" = list(column="MAF", operator=">", value=MAFthresholdDefault, filterOutNAs=filterOutMAFNAs)
        )
      }
    ),
    lapply(
      GQthresholds,
      function(GQthreshold) {
        list(
          "LowGQ" = list(column="GQ", operator="<", value=GQthreshold, filterOutNAs=TRUE),
          "LowDP" = list(column="DP", operator="<", value=DPthresholdDefault, filterOutNAs=TRUE),
          "HighMAF" = list(column="MAF", operator=">", value=MAFthresholdDefault, filterOutNAs=filterOutMAFNAs)
        )
      }
    )
  ),
#  overwriteExisting           = TRUE,
  overwriteExisting           = NULL,
  shouldUseSavedVersions      = TRUE,
  shouldUseExistingRda        = FALSE,
  minMeanMAFtoConsiderContam  = 0.01,
  plotFilestemExtra           = if(!is.null(MAFthresholds)) {
    "varyingMAF"
  } else if(!is.null(DPthresholds)) {
    "varyingDP"
  } else if(!is.null(GQthresholds)) {
    "varyingGQ"
  },
  setMonomorphicProgenyFilter = TRUE,
  monomorphicSkipChromosomes  = if(cross=="3d7_hb3") "Pf3D7_13_v3" else NULL,
  regionsMask                 = varRegions_v3("/data/malariagen2/plasmodium/pf-crosses/meta/regions_v3_rdp.bed"),
  shouldCreateQCFilteringPlots= TRUE,
  vcfCoreAllSamplesRda        = file.path(outputDirectory, cross, variantType, paste(cross, "vcfCoreFinalSamples.rda", sep=".")),
  shouldReturnVcfOnly         = FALSE
) {
  if(file.exists(vcfCoreAllSamplesRda) & shouldUseSavedVersions) {
    load(vcfCoreAllSamplesRda)
  } else {
    if(is.null(overwriteExisting)) {
      vcfList <- readAllChromosomesWithRegionsMask(
        cross                       = cross,
        variantType                 = variantType,
        genotypesDirectory          = genotypesDirectory,
        chromosomes                 = chromosomes,
        outputDirectory             = outputDirectory,
        genotypesFileFmt            = genotypesFileFmt,
        parentalStrains             = parentalStrains
      )
    } else {
      vcfList <- readAllChromosomesWithRegionsMask(
        cross                       = cross,
        variantType                 = variantType,
        genotypesDirectory          = genotypesDirectory,
        chromosomes                 = chromosomes,
        outputDirectory             = outputDirectory,
        genotypesFileFmt            = genotypesFileFmt,
        overwriteExisting           = overwriteExisting,
        parentalStrains             = parentalStrains
      )
    }
    vcfCoreAllSamples <- combineVcfListIntoVcf(vcfList)
    save(vcfCoreAllSamples, file=vcfCoreAllSamplesRda)
  }
  load(initialSampleQCresultsFile)
  vcfAnnotatedFinalSamplesFilename <- file.path(outputDirectory, cross, variantType, paste(cross, ".vcfAnnotatedFinalSamples.rda", sep=""))
  vcfAnnotatedBestReplicateSamplesFilename <- file.path(outputDirectory, cross, variantType, paste(cross, ".vcfAnnotatedBestReplicateSamples.rda", sep=""))
  vcfAnnotatedUncontaminatedSamplesFilename <- file.path(outputDirectory, cross, variantType, paste(cross, ".vcfAnnotatedUncontaminatedSamples.rda", sep=""))
  vcfList <- list()
  sampleSetNames <- sapply(sampleSets, function(x) x[["sampleSetName"]])
  if(shouldUseExistingRda && file.exists(vcfAnnotatedFinalSamplesFilename) && file.exists(vcfAnnotatedBestReplicateSamplesFilename) && file.exists(vcfAnnotatedUncontaminatedSamplesFilename)) {
    load(vcfAnnotatedFinalSamplesFilename)
    load(vcfAnnotatedBestReplicateSamplesFilename)
    load(vcfAnnotatedUncontaminatedSamplesFilename)
    vcfList[["final"]] <- vcfAnnotatedFinalSamples
    vcfList[["bestReplicate"]] <- vcfAnnotatedBestReplicateSamples
    vcfList[["uncontaminated"]] <- vcfAnnotatedUncontaminatedSamples
    rm(vcfAnnotatedFinalSamplesFilename)
    rm(vcfAnnotatedBestReplicateSamples)
    rm(vcfAnnotatedUncontaminatedSamples)
    gc()
  } else {
    if("final" %in% sampleSetNames) {
      finalSamples <- setdiff(dimnames(vcfCoreAllSamples)[[2]], initialSampleQCresults[["qcFailedSamples"]])
      vcfList[["final"]] <-  annotateVcf(vcfCoreAllSamples[, finalSamples])
    }
    if("bestReplicate" %in% sampleSetNames) {
      bestReplicateSamples <- setdiff(initialSampleQCresults[["uniqueSamples"]], initialSampleQCresults[["qcFailedSamples"]])
      vcfList[["bestReplicate"]] <-  annotateVcf(vcfCoreAllSamples[, bestReplicateSamples])
    }
    if("uncontaminated" %in% sampleSetNames) {
#      if(!exists("vcfInitialFiltered")) {
#        load(file.path(analysisDirectory, cross, variantType, paste(cross, ".vcfInitialFiltered.rda", sep="")))
#      }
#      vcfInitialFilteredPASS <- vcfInitialFiltered[filt(vcfInitialFiltered)=="PASS"]
      RefReads <- matrix(
        sapply(geno(vcfCoreAllSamples)[["AD"]], function(x) x[1]),
        ncol=dim(geno(vcfCoreAllSamples)[["AD"]])[2],
        dimnames=dimnames(geno(vcfCoreAllSamples)[["AD"]])
      )
      FirstAltReads <- matrix(
        sapply(geno(vcfCoreAllSamples)[["AD"]], function(x) x[2]),
        ncol=dim(geno(vcfCoreAllSamples)[["AD"]])[2],
        dimnames=dimnames(geno(vcfCoreAllSamples)[["AD"]])
      )
      MAF <- pmin(RefReads, FirstAltReads)/(RefReads+FirstAltReads)
      meanMAFperSample <- colMeans(MAF, na.rm = TRUE)
      uncontaminatedSamples <- intersect(names(which(meanMAFperSample < minMeanMAFtoConsiderContam)), dimnames(vcfCoreAllSamples)[[2]])
      vcfList[["uncontaminated"]] <-  annotateVcf(vcfCoreAllSamples[, uncontaminatedSamples])
    }
  }
  filterResultsList <- sapply(
    names(sampleSets),
    function(sampleSet) {
      if(is.null(parentalStrains)) {
        filterResults <- evaluateGenotypeFilters(
          vcfList[[sampleSets[[sampleSet]][["sampleSetName"]]]],
          plotFilestem                = file.path(outputDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", plotFilestemExtra, sampleSets[[sampleSet]][["sampleSetName"]], sep=".")),
          additionalInfoFilters       = sampleSets[[sampleSet]][["additionalInfoFilters"]],
          regionsMask                 = regionsMask,
          shouldCreateQCFilteringPlots=shouldCreateQCFilteringPlots,
          genotypeFiltersList         = genotypeFiltersList,
          setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
          monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
          sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
          shouldReturnVcfOnly         = shouldReturnVcfOnly
        )
      } else {
        filterResults <- evaluateGenotypeFilters(
          vcfList[[sampleSets[[sampleSet]][["sampleSetName"]]]],
          plotFilestem                = file.path(outputDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", plotFilestemExtra, sampleSets[[sampleSet]][["sampleSetName"]], sep=".")),
          additionalInfoFilters       = sampleSets[[sampleSet]][["additionalInfoFilters"]],
          regionsMask                 = regionsMask,
          shouldCreateQCFilteringPlots=shouldCreateQCFilteringPlots,
          genotypeFiltersList         = genotypeFiltersList,
          setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
          monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
          sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
          parentalIDs                 = parentalStrains,
          shouldReturnVcfOnly         = shouldReturnVcfOnly
        )
      }
      names(filterResults) <- paste(sampleSet, names(filterResults))
      filterResults
    },
    simplify=FALSE,
    USE.NAMES=TRUE
  )
  if(shouldReturnVcfOnly) {
    return(filterResultsList)
  } else {
    filterResults <- do.call(cbind, filterResultsList)
    return(filterResults)
  }
}


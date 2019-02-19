library(FourCSeq)
library(NOISeq)
library(BSgenome.Mmusculus.UCSC.mm9)
library(tidyverse)

metadata <- list(projectPath = "DataPPARG",  # directory to store results
                 fragmentDir = "re_fragments", # directory to store restriction fragment info
                 # subdir of projectPath
                 referenceGenomeFile = BSgenome.Mmusculus.UCSC.mm9, #path to the reference genome or a BSgenome object.
                 reSequence1 = "GAATTC", # restriction enzyme recognition pattern of the first restriction enzyme
                 reSequence2 = "GTAC", # reSequence2, restriction enzyme recognition pattern of the second restriction enzyme
                 primerFile = "data/PPARG/primerfile_pparg.fa", # path to file with primer sequences for 4C library prep
                 bamFilePath = "data/PPARG" #path to dir where bam files are stored
)


colData <- DataFrame(viewpoint = "testdata", 
                     condition = factor(rep(c("cd_zt6", 
                                              "cd_zt18", 
                                              "hf_zt6", 
                                              "hf_zt18"), 
                                            each = 2), 
                                        levels = c("cd_zt6", 
                                                   "cd_zt18", 
                                                   "hf_zt6", 
                                                   "hf_zt18")
                     ),
                     replicate = rep(c(1,2),4), 
                     bamFile = c("CD6_1_PPARG_mm9.bam",
                                 "CD6_2_PPARG_mm9.bam", 
                                 "CD18_1_PPARG_mm9.bam", 
                                 "CD18_2_PPARG_mm9.bam", 
                                 "HF6_1_PPARG_mm9.bam",
                                 "HF6_2_PPARG_mm9.bam", 
                                 "HF18_1_PPARG_mm9.bam", 
                                 "HF18_2_PPARG_mm9.bam"), 
                     sequencingPrimer = "first")

fc <- FourC(colData = colData, metadata = metadata)


###
fc<- addFragments(fc, minSize = 20, filter = TRUE, save=FALSE)
findViewpointFragments(fc)
fc <- addViewpointFrags(fc)
###
fc <- countFragmentOverlaps(fc, trim = 6, minMapq = 30)
fc <- combineFragEnds(fc)
###
writeTrackFiles(fc, format = "bw")
writeTrackFiles(fc, format='bedGraph')
###
fc <- smoothCounts(fc, assay = "counts", binWidth = 5)

# plotScatter(fc[, c("testdata_hf_zt6_1", 
#                    "testdata_hf_zt6_2")], 
#             xlab= "hfzt6_1", 
#             ylab="hfzt6_2", 
#             asp = 1)

fc@assays$data$counts[, c("testdata_hf_zt6_1",
                          "testdata_hf_zt6_2")
                      ]

# qqplot(x = fc@assays$data$counts[,"testdata_hf_zt6_1"], 
#        y = fc@assays$data$counts[,"testdata_hf_zt6_2"]
#        )

########

fcf <-getZScores(fc, 
                 removeZeros = TRUE, 
                 minCount = 40, 
                 minDist = NULL, 
                 fitFun = "distFitMonotoneSymmetric", 
                 sdFun = mad
                 )

zScore <- assay(fcf, "zScore")
my_samples <- colnames(zScore)
names(my_samples) <- colnames(zScore)

####
# my_histograms <- lapply(X = my_samples, FUN = function(i){
#   p <- qplot(zScore[,i], 
#              geom="histogram", 
#              binwidth = 0.05, 
#              xlab = i) 
# })
# sapply(my_histograms, plot)
####

# lapply(X = my_samples, FUN = function(i){
#   p <- qqnorm(zScore[,i], 
#               main = paste0("Normal Q-Q Plot - ",  i)
#               )
#   
#   abline(a=0, b=1)
#   plot(p)
# })
# for(i in my_samples){
#   qqnorm(zScore[,i], 
#          main = paste0("Normal Q-Q Plot - ",  i)
#   )
#   abline(a=0, b=1)
# }

###
fcf <- addPeaks(fcf, 
                zScoreThresh = 2.0, 
                fdrThresh = 0.05
                )

# for(i in seq_along(my_samples)){
#   plotFits(fcf[,i], main=my_samples[i])
# }
###
# plotZScores(fcf, 
#             cols = NULL, 
#             plotWindows = c(1e+05, 1e+06), 
#             controls = NULL, 
#             textsize = 20, 
#             txdb = TxDb.Mmusculus.UCSC.mm9.knownGene, 
#             plotSingle = FALSE
#             )

####
# my_rows <- as.data.frame(rowRanges(fcf)) %>% 
#   dplyr::select(seqnames, start, end, mid)
# 
# my_trafo <- assay(fcf, "trafo") 
# my_peaks <- assay(fcf, "peaks") 
# 
# my_trafo <- cbind(my_rows, my_trafo)
# my_peaks <- cbind(my_rows, my_peaks)
# 
# writeTrackFiles(object = fcf, 
#                 assay = "trafo", 
#                 folder = "tracks", 
#                 format = "bw", 
#                 removeZeros = FALSE)
# 
# writeTrackFiles(object = fcf, 
#                 assay = "trafo", 
#                 folder = "tracks", 
#                 format = "bedGraph", 
#                 removeZeros = FALSE)
# 
# names(assays(fcf))
# ###
# 
# fcf <- normalizeRPM(object = fcf, 
#                     assay = "counts", 
#                     normalized = "rpm")
# 
# ###
# writeTrackFiles(object = fcf, 
#                 assay = "fit", 
#                 folder = "tracks", 
#                 format = "bedGraph", 
#                 removeZeros = FALSE)
# 
# writeTrackFiles(object = fcf, 
#                 assay = "zScore", 
#                 folder = "tracks", 
#                 format = "bedGraph", 
#                 removeZeros = FALSE)
# 

###


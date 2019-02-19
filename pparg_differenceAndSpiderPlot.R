source("spider_plots.r")
library(FourCSeq)

#run everything in pipeline "pparg"
source("pipeline_pparg.R")

#get differences 
fcf <- getDifferences(fcf,referenceCondition="cd_zt6")

#analysis plots
plotDispEsts(fcf)

plotNormalizationFactors(fcf)

plotMA(results(fcf, 
               contrast=c("condition", "WE_68h", "MESO_68h")
               ),
       alpha=0.01,
       xlab="Mean 4C signal",
       ylab="log2 fold change",
       ylim=c(-3.1,3.1)
       )

#extract results

#here, the conditions should be modified for each experimental question
test_results <- as.data.frame(results(fcf, 
                            contrast =  c("condition", 
                                          "cd_zt6", 
                                          "hf_zt6"),
                            format = "GRanges")
                    )

#filter by log2FoldChange 

significant_results <- 
test_results %>% 
  filter(abs(log2FoldChange) >= 2 )


#extract info of the viewpoint fragment 
my_viewpoint <- fcf@colData[1,c("start", "end")]

#get midpoint of viewpoint
my_midpoint <- mean(my_viewpoint$start, my_viewpoint$end)

#make spider plot

my_significant_spider <- significant_results[,c("start", "end")]

min(my_significant_spider)
max(my_significant_spider)

(max(my_significant_spider) - min(my_significant_spider))*2

makeSpiderGramSingle(dom=my_significant_spider,#cbind(x1,x2), 
                     chrom.len=(max(my_significant_spider) - min(my_significant_spider))*2, #150e6, 
                     vp.loc=my_midpoint,#130e6, 
                     col=rgb(0,0,1)
                     ) 
